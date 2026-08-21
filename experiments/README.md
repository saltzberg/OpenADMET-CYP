# Experiment intents

These dated records capture hypotheses and planned evaluations for the OpenADMET
CYP challenge. They are not submissions. Each experiment keeps its canonical
editable intent in `README.md` and a standalone human-readable view in
`index.html`. Promotion into `submissions/` is a separate, explicitly approved
workflow.

Experiment directories use `YYYYMMDD_<simple-name>`, where the date is the UTC
intent date.

## Experiments

- [`20260821_OA-CYP-native-zero-shot`](20260821_OA-CYP-native-zero-shot/README.md)
  — run the released CYP-finetuned CheMeleon encoder and native four-output head
  unchanged; challenge labels enter only after prediction for evaluation.
- [`20260821_OA-CYP-native-affine-calibration`](20260821_OA-CYP-native-affine-calibration/README.md)
  — fit endpoint-specific affine maps to native predictions using grouped OOF
  challenge labels while leaving the released encoder and native head unchanged.
- [`20260821_OA-CYP-PLS-probe`](20260821_OA-CYP-PLS-probe/README.md)
  — compare a nested PLS frozen-encoder probe with the matched ridge probe.
- [`20260821_OA-CYP-interval-linear-probe`](20260821_OA-CYP-interval-linear-probe/README.md)
  — preregister an interval-loss linear probe; retain the run as invalid/blocked
  after its mandatory optimizer convergence gate failed.

- [`20260820_OA-CYP-finetuned-chemeleon`](20260820_OA-CYP-finetuned-chemeleon/README.md)
  — test OpenADMET's ChEMBL 37 CYP-finetuned CheMeleon representation for the
  four direct-inhibition pIC50 endpoints.
- [`20260820_public-cyp-tdi-models`](20260820_public-cyp-tdi-models/README.md)
  — test public CYP substrate, inhibitor, and reactant scores against CYP2D6
  and CYP3A4 TDI labels; current scope completed with DeepP450 deferred.
- [`20260820_cyp-series-validation`](20260820_cyp-series-validation/README.md)
  — repeat the matched representation comparison under broader deterministic
  ECFP4/Butina chemical-series groups.
- [`20260820_cyp-nested-affine-calibration`](20260820_cyp-nested-affine-calibration/README.md)
  — test leakage-safe affine range calibration of the CYP-finetuned ridge probe.
- [`20260820_cyp-pca-ridge-readout`](20260820_cyp-pca-ridge-readout/README.md)
  — test a prespecified PCA(64) plus ridge readout on the CYP-finetuned embedding.
- [`20260820_cyp-incremental-over-rdkit2d`](20260820_cyp-incremental-over-rdkit2d/README.md)
  — test whether CYP-finetuned predictions add nested residual signal beyond
  RDKit 2D.
