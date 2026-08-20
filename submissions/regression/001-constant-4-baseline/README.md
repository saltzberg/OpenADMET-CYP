# Regression submission 1: constant 4.0 baseline

Status: prepared; not submitted.

The script reads the shared blinded test set from `data/cyp-challenge-train-test/cyp-challenge-TEST-BLINDED.csv`, preserves its identifiers and row order, and sets all four direct-inhibition predictions to `4.0`.

- Prepared: `2026-08-19T19:14:26Z`
- Submitted: not submitted
- Leaderboard identity: `dargason`
- Model report: `https://github.com/saltzberg/OpenADMET-CYP/blob/main/submissions/regression/001-constant-4-baseline/README.md`
- Artifact: `openadmet-cyp_regression-sub_dargason_20260819T191426Z.csv`
- SHA-256: `42a27cbb800352d8c3d0506cbdc2a0199e6785d3c260b928a2ee81e2ea1666a2`
- Proprietary data: no
- Open-source declaration: no

The script uses the repository's blinded test CSV when present. Otherwise it downloads the public blinded test set from Hugging Face and requires SHA-256 `a342f8444a8dcb531ca12f3685293f0bd6c36ae9073f491e44a9bc1cc4b741f9`.

Reproduce the exact artifact:

```bash
python submissions/regression/001-constant-4-baseline/create_submission.py \
  --timestamp 20260819T191426Z
```

The output is `openadmet-cyp_regression-sub_dargason_20260819T191426Z.csv` and must match the recorded SHA-256.
