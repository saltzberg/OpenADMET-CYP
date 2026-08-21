# Native-head affine calibration

This model applies one endpoint-specific affine map to the released OpenADMET ChEMBL 37 CYP CheMeleon model's native predictions. The released molecular encoder and native four-output neural head are unchanged.

## Model

```text
challenge SMILES
→ released CYP-finetuned CheMeleon encoder [unchanged]
→ released native four-output FFN [unchanged]
→ native endpoint prediction
→ intercept + slope [fitted on challenge training labels]
→ calibrated pIC50 prediction
```

For endpoint `e`, the submitted prediction is

```text
calibrated_e = intercept_e + slope_e × native_e
```

The final coefficients are:

| Endpoint | Intercept | Slope |
|---|---:|---:|
| CYP1A2 | 1.232525445548 | 0.697933314620 |
| CYP2C9 | 1.293264085224 | 0.617372751411 |
| CYP2D6 | 2.112812270247 | 0.489984338819 |
| CYP3A4 | -0.867715737338 | 0.929491692017 |

No clipping, nonlinear model, encoder fine-tuning, native-head fine-tuning, ensemble, or residual correction is applied.

## Inputs and outputs

Input is the challenge-provided `SMILES` for each of the 750 blinded compounds. The released OpenADMET model produces four native CYP pIC50 values. The fixed affine maps produce the four direct-inhibition pIC50 output columns required by the regression track.

## Training and data

The OpenADMET source model was released as `openadmet/cyp1a2-cyp2d6-cyp3a4-cyp2c9-chemeleon-v1`, revision `ef24cf941ae21c7d7a64df378a846bd2066eceda`, and was trained by OpenADMET on ChEMBL 37 CYP pIC50 records.

Only the affine intercepts and slopes were fitted in this work. They were evaluated using five grouped outer folds over the challenge training compounds. For every OOF prediction, the affine map was fitted without that compound's outer fold. Final blinded maps were then fitted from all available challenge training labels for the corresponding endpoint.

OOF macro ST-RAE changed from `0.8981` for the native predictions to `0.5734` after affine calibration. Every endpoint improved in the paired grouped OOF analysis. The affine model does not improve molecular ordering within a fitted map; it adjusts endpoint offset and scale.

## Expected behavior

The source model's native predictions are shifted and compressed relative to the challenge training assays. The affine maps primarily correct systematic endpoint calibration while retaining the native model's ordering within each final endpoint map.

## Limitations

- ChEMBL pIC50 records and the challenge assays are not interchangeable measurement contexts.
- Affine calibration cannot recover ranking information absent from the native model.
- Source-model training overlap limits a completely external transfer interpretation.
- OOF validation used exact Bemis–Murcko groups, which are mostly singletons in this dataset.
- Blinded performance is unknown until challenge evaluation.

## Reproduction

Install the OpenADMET model runtime, download the released model directory, then run:

```bash
python submissions/regression/002-native-head-affine-calibration/reproduce_from_model.py \
  --test-csv data/cyp-challenge-train-test/cyp-challenge-TEST-BLINDED.csv \
  --model-dir /path/to/cyp1a2-cyp2d6-cyp3a4-cyp2c9-chemeleon-v1/anvil_training \
  --coefficients submissions/regression/002-native-head-affine-calibration/affine_coefficients.csv \
  --output reproduced.csv \
  --accelerator cpu
```

To rebuild the canonical CSV directly from the preserved calibrated long-form artifact:

```bash
python submissions/regression/002-native-head-affine-calibration/create_submission.py \
  --output submissions/regression/002-native-head-affine-calibration/openadmet-cyp_regression-sub_dargason_20260821T021138Z.csv
```
