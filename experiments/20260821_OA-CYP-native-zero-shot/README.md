# OpenADMET CYP native zero-shot predictions

> **Status:** Completed native zero-shot evaluation  
> **Intent date (UTC):** 20260821  
> **Experiment ID:** `20260821_OA-CYP-native-zero-shot`

## Summary

Use the released OpenADMET ChEMBL 37 CYP CheMeleon encoder and its released
four-output neural prediction head unchanged. Challenge SMILES are inputs;
challenge pIC50 labels are joined only after prediction to evaluate the training
subset. No model stage is fitted, selected, calibrated, clipped, or stacked using
challenge labels.

## Hypothesis

The released OpenADMET ChEMBL 37 CYP CheMeleon model and its native four-output neural head transfer to the challenge direct-inhibition endpoints without fitting any stage on challenge labels.

## Endpoints

- `CYP1A2_pIC50_direct_inhibition`
- `CYP2C9_pIC50_direct_inhibition`
- `CYP2D6_pIC50_direct_inhibition`
- `CYP3A4_pIC50_direct_inhibition`

## Representation and inputs

- Challenge training and blinded canonical input SMILES.
- Released model revision `ef24cf941ae21c7d7a64df378a846bd2066eceda`.
- Released CYP-finetuned CheMeleon message-passing encoder.
- Released native three-layer FFN and four CYP outputs.

## Planned evaluation

1. Run the released model on training SMILES without providing labels to the
   prediction process.
2. Run the identical released model on blinded SMILES.
3. Join training predictions to reported point estimates and confidence
   intervals only after prediction.
4. Report ST-RAE, MAE, R², Spearman's rho, Kendall's tau, spread, and coverage.
5. Preserve raw outputs, commands, input/model hashes, and a prediction-flow
   figure.

## Comparisons

- Observed training distributions versus native training predictions.
- Native training versus native blinded prediction distributions.
- No foundation-model CYP comparison is claimed because the untouched
  CheMeleon encoder has no released native four-CYP pIC50 head.

## Risks and confounders

- Audit overlap between representation pretraining/fine-tuning data and challenge records before interpreting validation performance.
- ChEMBL pIC50 and the challenge assays are not interchangeable measurement
  contexts.
- Source-model training used no held-out split; source overlap limits a pure
  external-transfer interpretation.

## Decision rule

Treat this as the vanilla zero-shot anchor. It describes the released model's
unaltered transfer behavior; it is not selected or promoted using leaderboard
feedback. Any challenge-trained head or calibration belongs in a separate
experiment.

## Provenance

- Intent recorded on 20260821 UTC.
- No predictions, scores, or submission artifact are asserted by this record.

## Results

The released model was run unchanged on 4,905 challenge-training SMILES and 750
blinded SMILES. No challenge label, confidence interval, fold, fitted scaler,
probe, calibrator, clipping rule, or stack entered prediction generation.

| Endpoint | ST-RAE | MAE | R² | Spearman's rho | Prediction SD | Observed SD |
|---|---:|---:|---:|---:|---:|---:|
| CYP1A2 | 0.6695 | 0.7692 | -0.0643 | 0.2827 | 0.4361 | 1.0305 |
| CYP2C9 | 0.8862 | 0.8517 | -0.8481 | 0.2841 | 0.3874 | 0.7823 |
| CYP2D6 | 1.0216 | 0.8668 | -0.5362 | 0.2200 | 0.4331 | 0.9161 |
| CYP3A4 | 1.0152 | 1.3109 | -1.1588 | 0.3656 | 0.4368 | 1.0932 |

The native heads retain weak-to-moderate rank signal but are compressed and
miscalibrated for the challenge assays. This is the zero-shot anchor; it is not
a candidate selected using challenge performance.

## Outputs

- `index.html` — report with explicit prediction construction.
- `figures/01_prediction_flow.*` — data/model flow.
- `artifacts/native_train_raw.csv` and `native_blind_raw.csv` — unmodified CLI outputs.
- `artifacts/native_train_scored.parquet` — post-prediction training-label join.
- `artifacts/native_blind_predictions.csv` — 3,000 native blind predictions.
- `artifacts/run_manifest.json` — commands, model/input hashes, and zero-shot invariant.
