# Constant 4.0 baseline

This model predicts a direct-inhibition pIC50 of `4.0` for every compound and every CYP endpoint: CYP1A2, CYP2C9, CYP2D6, and CYP3A4.

## Model

Just set everything to 4.0

## Purpose

This is to test the submission scripts.  

## Reproduction

```bash
python submissions/regression/001-constant-4-baseline/create_submission.py \
  --timestamp 20260819T191426Z
```
