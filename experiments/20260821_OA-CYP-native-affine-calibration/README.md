# OpenADMET native-head affine calibration

> **Status:** Completed experiment — not a submission  
> **Intent date (UTC):** 20260821  
> **Experiment ID:** `20260821_OA-CYP-native-affine-calibration`

## Summary

Experiment level: **adapted native model**. Start from the released model's
already-generated native predictions; fit only one intercept and slope per
endpoint using challenge-training labels. The released encoder and native FFN
remain unchanged.

## Hypothesis

A leakage-safe endpoint-specific affine map can correct offset and scale in native zero-shot OpenADMET predictions without changing molecular ordering.

## Endpoints

- CYP1A2, CYP2C9, CYP2D6, and CYP3A4 direct-inhibition pIC50.

## Representation and inputs

```text
challenge SMILES
→ released CYP-finetuned CheMeleon encoder [released; frozen]
→ released native four-output FFN [released; frozen]
→ native endpoint prediction
→ affine intercept + slope [fitted on outer challenge-training folds]
→ calibrated prediction
```

The parent native predictions are from
`experiments/20260821_OA-CYP-native-zero-shot/`. Challenge labels fit only the
affine map; they do not fit the encoder or native FFN.

## Planned evaluation

1. Retain the parent's native training predictions unchanged.
2. Use the established five grouped outer folds.
3. Within each endpoint and outer fold, fit ordinary least-squares intercept and
   slope on outer-training native predictions and labels only.
4. Apply that map to the untouched outer holdout and assemble OOF predictions.
5. Compare native versus calibrated ST-RAE, rank metrics, spread, and paired
   scaffold-bootstrap differences.
6. Only after OOF evaluation, fit one affine map on all endpoint labels and apply
   it to native blinded predictions.

## Comparisons

- Unchanged native zero-shot predictions on identical endpoint rows.

## Risks and confounders

- Audit overlap between representation pretraining/fine-tuning data and challenge records before interpreting validation performance.
- Affine maps can correct offset and scale but cannot repair within-map molecular
  ordering.
- Five independently fitted outer maps can slightly change assembled global OOF
  ranks even when every fitted slope is positive.

## Decision rule

Support requires a paired macro ST-RAE improvement without a reproducible
endpoint regression. Ranking should remain unchanged apart from ties.

## Provenance

- Intent recorded on 20260821 UTC.
- No predictions, scores, or submission artifact are asserted by this record.

## Results

<!-- EXECUTED_RESULTS_START -->
| Endpoint | Raw ST-RAE | Affine ST-RAE | Δ affine−raw | paired 95% CI low | paired 95% CI high | Raw Spearman ρ | Affine Spearman ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| CYP1A2 | 0.669469545 | 0.619721829 | -0.049747716 | -0.073725568 | -0.025027275 | 0.282740991 | 0.280922070 |
| CYP2C9 | 0.886225859 | 0.501177151 | -0.385048707 | -0.442084277 | -0.329627972 | 0.284065022 | 0.279077565 |
| CYP2D6 | 1.021610249 | 0.644680449 | -0.376929800 | -0.432750141 | -0.321856103 | 0.219967563 | 0.218872120 |
| CYP3A4 | 1.015180659 | 0.527839720 | -0.487340939 | -0.538695280 | -0.435733278 | 0.365595734 | 0.364563912 |
| MA | 0.898121578 | 0.573354787 | -0.324766791 | -0.349625332 | -0.297364519 | 0.288092327 | 0.285858917 |

The unweighted endpoint-macro ST-RAE difference was **-0.324766791** (2,000-replicate paired scaffold-bootstrap 95% CI -0.349625332 to -0.297364519). A reproducible endpoint regression was defined as an endpoint paired interval wholly above zero; none occurred. The preregistered no-endpoint-regression decision rule was **MET**.

### Final all-training affine maps applied to blinded native predictions

| Endpoint | Intercept | Slope | Training rows |
|---|---:|---:|---:|
| CYP1A2 | 1.232525445548 | 0.697933314620 | 1412 |
| CYP2C9 | 1.293264085224 | 0.617372751411 | 1285 |
| CYP2D6 | 2.112812270247 | 0.489984338819 | 1493 |
| CYP3A4 | -0.867715737338 | 0.929491692017 | 2335 |

All fitted slopes were positive. Therefore each individual outer-fold map and each final all-training endpoint map preserved ordering (including ties) within the rows to which that one map was applied. The assembled OOF endpoint ranks changed slightly because five independently fitted maps place different outer holdouts on different affine scales; the exact raw and affine Spearman values are reported above. No clipping was applied. Blinded labels were not available or used.

### Reproduction

```text
/home/dan/swr/miniconda3/envs/cheminf/bin/python src/run_experiment.py
/home/dan/swr/miniconda3/envs/cheminf/bin/python src/test_experiment.py
/home/dan/swr/miniconda3/envs/cheminf/bin/python src/verify_outputs.py
```
<!-- EXECUTED_RESULTS_END -->
