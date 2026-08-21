#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
E=Path(__file__).resolve().parents[1]; R=E.parents[1]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 m=json.loads((E/'artifacts/run_manifest.json').read_text()); assert m['status']=='completed_zero_shot' and m['challenge_label_training'] is False
 assert m['prediction_flow']==['challenge SMILES','released CYP-finetuned CheMeleon encoder','released native 3-layer FFN','four native predictions']
 for p,h in m['source_files'].items(): assert sha(R/p)==h,p
 for p,h in m['outputs'].items(): assert sha(E/p)==h,p
 train=pd.read_parquet(E/'artifacts/native_train_scored.parquet'); blind=pd.read_csv(E/'artifacts/native_blind_predictions.csv'); met=pd.read_csv(E/'artifacts/metrics.csv')
 assert len(train)==6525 and len(blind)==3000 and len(met)==4
 assert not train.duplicated(['compound_id','endpoint']).any() and not blind.duplicated(['compound_id','endpoint']).any()
 assert np.isfinite(train.prediction).all() and np.isfinite(blind.prediction).all()
 assert set(blind.endpoint)=={'CYP1A2','CYP2C9','CYP2D6','CYP3A4'} and blind.groupby('endpoint').size().eq(750).all()
 assert len(m['commands'])==2 and all('--model-dir' in c and '--input-col' in c for c in m['commands'])
 html=(E/'index.html').read_text(); assert html.count('<figure>')==3 and 'No encoder, head, calibration, or post-processing is fitted on challenge labels.' not in html
 for stem in ['01_prediction_flow','02_distributions','03_observed_vs_native']:
  for suffix in ['.png','.svg']: assert (E/'figures'/f'{stem}{suffix}').stat().st_size>1000
 print('PASS zero-shot invariant, hashes, 6,525 scored rows, 3,000 blind predictions, figures, HTML')
 return 0
if __name__=='__main__': raise SystemExit(main())
