# All-False baseline

This model predicts that every compound is not a time-dependent inhibitor for both CYP2D6 and CYP3A4.

## Model

There is no training step and the model does not use molecular structure or assay data. It is a constant negative classifier.

## Purpose

The model is a submission-system smoke test and a deliberately simple reference point. It checks that compound identifiers, Boolean endpoint columns, file validation, scoring, and reporting all work end to end.

The model has no ability to distinguish compounds and will miss every positive TDI example.

## Reproduction

```bash
python submissions/classification/001-all-false-baseline/create_submission.py \
  --timestamp 20260819T191426Z
```
