#!/usr/bin/env python3
"""Build the canonical submission from calibrated long-form predictions."""
from __future__ import annotations
import argparse, hashlib
from pathlib import Path
import pandas as pd

REPOSITORY = Path(__file__).resolve().parents[3]
TEST = REPOSITORY / "data/cyp-challenge-train-test/cyp-challenge-TEST-BLINDED.csv"
LONG = REPOSITORY / "experiments/20260821_OA-CYP-native-affine-calibration/artifacts/calibrated_blind_predictions.csv"
TEST_SHA256 = "a342f8444a8dcb531ca12f3685293f0bd6c36ae9073f491e44a9bc1cc4b741f9"
ENDPOINTS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
COLUMNS = [f"{endpoint}_pIC50_direct_inhibition" for endpoint in ENDPOINTS]

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main() -> int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--test-csv",type=Path,default=TEST)
    parser.add_argument("--calibrated-long",type=Path,default=LONG)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if sha256(args.test_csv)!=TEST_SHA256: raise SystemExit("blinded test SHA-256 mismatch")
    test=pd.read_csv(args.test_csv)
    long=pd.read_csv(args.calibrated_long)
    required={"compound_id","SMILES","endpoint","calibrated_prediction"}
    if not required.issubset(long): raise SystemExit(f"missing long columns: {sorted(required-set(long))}")
    if len(long)!=3000 or long.duplicated(["compound_id","endpoint"]).any(): raise SystemExit("expected one row per 750 compound × 4 endpoints")
    if set(long.endpoint)!=set(ENDPOINTS): raise SystemExit("endpoint mismatch")
    smiles=long.groupby("compound_id").SMILES.nunique()
    if not smiles.eq(1).all(): raise SystemExit("conflicting SMILES in long predictions")
    wide=long.pivot(index="compound_id",columns="endpoint",values="calibrated_prediction").rename(columns=dict(zip(ENDPOINTS,COLUMNS)))
    output=test[["SMILES","Molecule_Name"]].merge(wide,left_on="Molecule_Name",right_index=True,how="left",validate="one_to_one")
    output=output[["SMILES","Molecule_Name",*COLUMNS]]
    if len(output)!=750 or output[COLUMNS].isna().any().any(): raise SystemExit("incomplete canonical output")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    output.to_csv(args.output,index=False,lineterminator="\n")
    print(args.output)
    print(sha256(args.output))
    return 0
if __name__=="__main__": raise SystemExit(main())
