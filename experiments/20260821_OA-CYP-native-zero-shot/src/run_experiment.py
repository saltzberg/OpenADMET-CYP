#!/usr/bin/env python3
"""Run the released OpenADMET CYP model with no challenge-label fitting."""
from __future__ import annotations
import hashlib, json, os, platform, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.stats import kendalltau, spearmanr
import sklearn
from sklearn.metrics import mean_absolute_error, r2_score

ROOT=Path(__file__).resolve().parents[3]
EXP=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/cyp-challenge-train-test"
MODEL=ROOT/"cyp_feature_store/models/openadmet_cyp_chemeleon_v1/anvil_training"
OPENADMET=Path("/home/dan/swr/venvs/openadmet-models/bin/openadmet")
ENDPOINTS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
TARGET={e:f"{e}_pIC50_direct_inhibition" for e in ENDPOINTS}
PRED={e:f"OADMET_PRED_chemprop_OPENADMET_LOGAC50_{e}" for e in ENDPOINTS}
MODEL_REVISION="ef24cf941ae21c7d7a64df378a846bd2066eceda"
CREATED="2026-08-21T01:09:25Z"

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()

def st_error(p,lo,hi): return np.maximum(0,np.maximum(lo-p,p-hi))
def st_rae(y,p,lo,hi): return float(st_error(p,lo,hi).sum()/np.abs(y-y.mean()).sum())
def axes(ax):
 ax.spines[['top','right']].set_visible(False); ax.spines[['left','bottom']].set_color('#666'); ax.grid(False)
def save(fig,base):
 fig.savefig(base.with_suffix('.png'),dpi=300,bbox_inches='tight',facecolor='white'); fig.savefig(base.with_suffix('.svg'),bbox_inches='tight',facecolor='white'); plt.close(fig)
def run_native(source,output):
 env=os.environ.copy(); env.update({'TABPFN_TELEMETRY_OPTOUT':'1','OADMET_NO_RICH_LOGGING':'1'})
 command=[str(OPENADMET),'predict','--input-path',str(source),'--input-col','SMILES','--model-dir',str(MODEL),'--output-csv',str(output),'--accelerator','gpu']
 completed=subprocess.run(command,cwd=ROOT,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (EXP/'artifacts'/f'{output.stem}.log').write_text(completed.stdout)
 if completed.returncode: raise RuntimeError(completed.stdout)
 return command

def flowchart():
 fig,ax=plt.subplots(figsize=(10,2.8)); ax.set_xlim(0,10); ax.set_ylim(0,3); ax.axis('off')
 items=[(0.5,'Challenge\nSMILES'),(2.5,'Released ChEMBL 37\nCYP-finetuned\nCheMeleon encoder'),(5.3,'Released native FFN\n2048 → 1024 → 1024\n→ 1024 → 4 CYP outputs'),(8.3,'Native CYP\npIC50 predictions')]
 for x,text in items: ax.text(x,1.7,text,ha='center',va='center',fontsize=10)
 for (x,_),(x2,_) in zip(items[:-1],items[1:]): ax.annotate('',xy=(x2-.75,1.7),xytext=(x+.75,1.7),arrowprops={'arrowstyle':'->','lw':1,'color':'#555'})
 ax.text(4.8,.35,'No encoder, head, calibration, or post-processing is fitted on challenge labels.',ha='center',fontsize=10,color='#7a3e00')
 ax.text(8.3,.35,'Training labels are joined only after prediction for evaluation.',ha='center',fontsize=9,color='#555')
 save(fig,EXP/'figures/01_prediction_flow')
def plots(long,blind):
 colors={'Observed training':'#222','Native train prediction':'#2f5d8a','Native blind prediction':'#a6531a'}
 fig,axs=plt.subplots(2,2,figsize=(8,6.5),sharex=True,sharey=True)
 for ax,e in zip(axs.flat,ENDPOINTS):
  g=long[long.endpoint==e]
  vals={'Observed training':g.observed.to_numpy(float),'Native train prediction':g.prediction.to_numpy(float),'Native blind prediction':blind.loc[blind.endpoint==e,'prediction'].to_numpy(float)}
  for name,x in vals.items():
   x=np.sort(x); ax.plot(x,np.arange(1,len(x)+1)/len(x),color=colors[name],lw=1.3)
  ax.set_title(e,loc='left'); ax.set_xlim(1.5,8); ax.set_ylim(0,1); axes(ax)
 axs[1,0].set_xlabel('pIC50'); axs[1,1].set_xlabel('pIC50'); axs[0,0].set_ylabel('Cumulative fraction'); axs[1,0].set_ylabel('Cumulative fraction')
 fig.legend([Line2D([0],[0],color=c,lw=1.4) for c in colors.values()],list(colors),loc='lower center',ncol=3,frameon=False,fontsize=8.5); fig.subplots_adjust(bottom=.13,hspace=.28,wspace=.18); save(fig,EXP/'figures/02_distributions')
 fig,axs=plt.subplots(2,2,figsize=(7.2,7.2),sharex=True,sharey=True)
 for ax,e in zip(axs.flat,ENDPOINTS):
  g=long[long.endpoint==e]; ax.scatter(g.observed,g.prediction,s=8,alpha=.42,color='#2f5d8a',linewidths=0,rasterized=True); ax.plot([1.5,8],[1.5,8],color='#555',lw=.7); ax.set_xlim(1.5,8); ax.set_ylim(1.5,8); ax.set_aspect('equal'); ax.set_title(e,loc='left'); axes(ax)
 axs[1,0].set_xlabel('Observed pIC50'); axs[1,1].set_xlabel('Observed pIC50'); axs[0,0].set_ylabel('Native predicted pIC50'); axs[1,0].set_ylabel('Native predicted pIC50'); fig.subplots_adjust(hspace=.2,wspace=.18); save(fig,EXP/'figures/03_observed_vs_native')
def render(metrics):
 rows=''.join(f'<tr><td>{r.endpoint}</td><td>{r.st_rae:.4f}</td><td>{r.mae:.4f}</td><td>{r.spearman_rho:.4f}</td><td>{r.n}</td></tr>' for r in metrics.itertuples())
 edited=datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00','Z')
 html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>OpenADMET native zero-shot CYP predictions</title><style>body{{color:#222;font:16px/1.5 Georgia,serif;margin:0}}main{{max-width:62rem;margin:3rem auto;padding:0 1rem}}h1{{font-size:2.5rem}}h2{{margin-top:2.3rem;border-bottom:1px solid #ccc}}.meta{{color:#666;font:13px ui-monospace,monospace}}img{{width:min(100%,900px)}}table{{border-collapse:collapse;width:100%;font:14px system-ui}}th,td{{padding:.4rem;border-bottom:1px solid #ddd;text-align:right}}th:first-child,td:first-child{{text-align:left}}</style></head><body><main><div class="meta">20260821_OA-CYP-native-zero-shot · Created {CREATED} · Last edited {edited}</div><h1>OpenADMET CYP native zero-shot predictions</h1><p><strong>Purpose:</strong> evaluate the released OpenADMET CYP-ChEMBL-finetuned CheMeleon encoder and native prediction head on challenge training and blind SMILES without fitting on challenge labels.</p><h2>Prediction construction</h2><figure><img src="figures/01_prediction_flow.png"><figcaption>Released weights are used unchanged. Challenge labels enter only after prediction for evaluation.</figcaption></figure><h2>Training-set metrics</h2><table><thead><tr><th>Endpoint</th><th>ST-RAE</th><th>MAE</th><th>Spearman ρ</th><th>n</th></tr></thead><tbody>{rows}</tbody></table><h2>Observed and predicted distributions</h2><figure><img src="figures/02_distributions.png"></figure><h2>Observed versus native prediction</h2><figure><img src="figures/03_observed_vs_native.png"></figure><h2>Artifacts</h2><p><a href="artifacts/metrics.csv">metrics</a> · <a href="artifacts/native_blind_predictions.csv">blind predictions</a> · <a href="artifacts/run_manifest.json">manifest</a> · <a href="README.md">canonical record</a></p></main></body></html>'''
 (EXP/'index.html').write_text(html)
def main():
 for d in [EXP/'artifacts',EXP/'figures']: d.mkdir(exist_ok=True)
 train_path=DATA/'cyp-challenge-TRAIN_inhibition.csv'; blind_path=DATA/'cyp-challenge-TEST-BLINDED.csv'
 raw_train=EXP/'artifacts/native_train_raw.csv'; raw_blind=EXP/'artifacts/native_blind_raw.csv'
 commands=[run_native(train_path,raw_train),run_native(blind_path,raw_blind)]
 train=pd.read_csv(raw_train); source=pd.read_csv(train_path); blind=pd.read_csv(raw_blind); blind_source=pd.read_csv(blind_path)
 assert train.Molecule_Name.tolist()==source.Molecule_Name.tolist() and blind.Molecule_Name.tolist()==blind_source.Molecule_Name.tolist()
 parts=[]; bparts=[]; rows=[]
 for e in ENDPOINTS:
  t=TARGET[e]; m=source[t].notna(); y=source.loc[m,t].to_numpy(float); p=train.loc[m,PRED[e]].to_numpy(float); lo=source.loc[m,f'{t}_conf_low'].to_numpy(float); hi=source.loc[m,f'{t}_conf_high'].to_numpy(float)
  parts.append(pd.DataFrame({'compound_id':source.loc[m,'Molecule_Name'],'endpoint':e,'observed':y,'conf_low':lo,'conf_high':hi,'prediction':p}))
  bparts.append(pd.DataFrame({'compound_id':blind_source.Molecule_Name,'SMILES':blind_source.SMILES,'endpoint':e,'prediction':blind[PRED[e]]}))
  rows.append({'endpoint':e,'st_rae':st_rae(y,p,lo,hi),'mae':mean_absolute_error(y,p),'r2':r2_score(y,p),'spearman_rho':spearmanr(y,p).statistic,'kendall_tau':kendalltau(y,p).statistic,'n':len(y),'prediction_sd':np.std(p,ddof=1),'observed_sd':np.std(y,ddof=1)})
 long=pd.concat(parts,ignore_index=True); blong=pd.concat(bparts,ignore_index=True); metrics=pd.DataFrame(rows)
 long.to_parquet(EXP/'artifacts/native_train_scored.parquet',index=False); blong.to_csv(EXP/'artifacts/native_blind_predictions.csv',index=False); metrics.to_csv(EXP/'artifacts/metrics.csv',index=False)
 flowchart(); plots(long,blong); render(metrics)
 outputs={str(p.relative_to(EXP)):sha(p) for p in sorted(EXP.rglob('*')) if p.is_file() and p.name!='run_manifest.json' and '__pycache__' not in p.parts}
 manifest={'status':'completed_zero_shot','created_at_utc':datetime.now(timezone.utc).isoformat(),'experiment_id':EXP.name,'challenge_label_training':False,'prediction_flow':['challenge SMILES','released CYP-finetuned CheMeleon encoder','released native 3-layer FFN','four native predictions'],'model_revision':MODEL_REVISION,'commands':commands,'source_files':{str(train_path.relative_to(ROOT)):sha(train_path),str(blind_path.relative_to(ROOT)):sha(blind_path),str((MODEL/'model.pth').relative_to(ROOT)):sha(MODEL/'model.pth'),str((MODEL/'model.json').relative_to(ROOT)):sha(MODEL/'model.json')},'runtime':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__},'outputs':outputs}
 (EXP/'artifacts/run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); print(metrics.to_string(index=False)); return 0
if __name__=='__main__': raise SystemExit(main())
