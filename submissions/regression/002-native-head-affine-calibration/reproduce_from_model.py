#!/usr/bin/env python3
"""Reproduce native predictions, apply fixed affine maps, and write submission CSV."""
from __future__ import annotations
import argparse, hashlib, os, subprocess, tempfile
from pathlib import Path
import pandas as pd

TEST_SHA256="a342f8444a8dcb531ca12f3685293f0bd6c36ae9073f491e44a9bc1cc4b741f9"
ENDPOINTS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
NATIVE={e:f"OADMET_PRED_chemprop_OPENADMET_LOGAC50_{e}" for e in ENDPOINTS}
OUTPUT={e:f"{e}_pIC50_direct_inhibition" for e in ENDPOINTS}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(); p.add_argument('--test-csv',type=Path,required=True); p.add_argument('--model-dir',type=Path,required=True); p.add_argument('--coefficients',type=Path,default=Path(__file__).with_name('affine_coefficients.csv')); p.add_argument('--output',type=Path,required=True); p.add_argument('--accelerator',choices=['cpu','gpu'],default='cpu'); a=p.parse_args()
 if sha(a.test_csv)!=TEST_SHA256: raise SystemExit('blinded test SHA-256 mismatch')
 test=pd.read_csv(a.test_csv); coeff=pd.read_csv(a.coefficients); coeff=coeff[coeff.fit_scope.eq('all_training_labels_for_blind_application')].set_index('endpoint')
 if set(coeff.index)!=set(ENDPOINTS): raise SystemExit('final coefficient rows missing')
 with tempfile.TemporaryDirectory() as tmp:
  raw=Path(tmp)/'native.csv'; env=os.environ.copy(); env.update({'TABPFN_TELEMETRY_OPTOUT':'1','OADMET_NO_RICH_LOGGING':'1'})
  subprocess.run(['openadmet','predict','--input-path',str(a.test_csv),'--input-col','SMILES','--model-dir',str(a.model_dir),'--output-csv',str(raw),'--accelerator',a.accelerator],check=True,env=env)
  native=pd.read_csv(raw)
 result=test[['SMILES','Molecule_Name']].copy()
 for e in ENDPOINTS: result[OUTPUT[e]]=float(coeff.loc[e,'intercept'])+float(coeff.loc[e,'slope'])*native[NATIVE[e]].to_numpy(float)
 a.output.parent.mkdir(parents=True,exist_ok=True); result.to_csv(a.output,index=False,lineterminator='\n'); print(a.output); print(sha(a.output)); return 0
if __name__=='__main__': raise SystemExit(main())
