# Constant 4.0 baseline

This model predicts a direct-inhibition pIC50 of `4.0` for every compound and every CYP endpoint: CYP1A2, CYP2C9, CYP2D6, and CYP3A4.

## Model

There is no training step and the model does not use molecular structure or assay data. It is a constant predictor at the assay's lower activity threshold.

## Purpose

The model is a submission-system smoke test and a deliberately simple reference point. It checks that compound identifiers, endpoint columns, file validation, scoring, and reporting all work end to end.

Because every prediction is identical, the model cannot rank compounds. Its Spearman correlations should therefore carry no useful discrimination signal.

## Reproduction

```bash
python submissions/regression/001-constant-4-baseline/create_submission.py \
  --timestamp 20260819T191426Z
```
