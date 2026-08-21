# 002-cypreact-raw

This submission uses the frozen [CypReact v1.2](https://github.com/Le0nT1/CypReact_old/releases/tag/1.2) endpoint models, maps their public `R` (reactant) and `N` (non-reactant) calls directly to `True` and `False`, and outputs CYP2D6 and CYP3A4 TDI classifications.

[CypReact](https://doi.org/10.1021/acs.jcim.8b00035) predicts whether a molecule is metabolized by specific human cytochrome P450 enzymes. The input is molecular fragment descriptors and MACCS structural keys and the output is a binary reactant or non-reactant call.

## Model

**Information flow**

```text
challenge SMILES
→ ID-labelled SDF
→ released CypReact v1.2 endpoint models [frozen]
→ public R/N reactant calls [frozen interface]
→ R=True, N=False [direct mapping]
→ CYP2D6_is_TDI and CYP3A4_is_TDI

Challenge-fitted or challenge-selected stages: none
```

**Training details**

- **Base model:** CypReact v1.2, using the base CYP2D6 and CYP3A4 models.
- **Prediction interface:** the frozen public binary `R`/`N` interface; it exposes neither probabilities nor confidence scores.
- **Challenge-data use:** no challenge labels were used for any training

The adapter preserves blinded-test row order and identifiers and fails closed if either endpoint is missing or invalid. Exact inputs, runtime, model assets, commands, and checksums are recorded in [`inference_manifest.json`](inference_manifest.json).

## Results and discussion

Performance on OpenADMET training labels was evaluated using 2,000-replicate Bemis–Murcko scaffold-cluster bootstrap intervals. 

| Endpoint | Coverage | Average precision (95% CI) | ROC-AUC (95% CI) | Positive calls |
|---|---:|---:|---:|---:|
| CYP2D6 | 1.000 | 0.219 (0.198–0.241) | 0.508 (0.497–0.518) | 96.3% |
| CYP3A4 | 1.000 | 0.214 (0.199–0.230) | 0.502 (0.500–0.504) | 99.4% |

CypReact predicts almost everything to be reactive.  Thus, as a TDI predictor, this is fairly uninformative.  

Are the negative calls informative, however?

| Endpoint | # CypReact_non-reactive | TDI-positive | TDI-negative | TDI-negative rate | Baseline TDI-negative rate | Absolute % Enrichment |
|---|---:|---:|---:|---:|---:|---:|
| CYP2D6 | 55 | 8 | 47 | 85.5% | 78.4% | +7.1% |
| CYP3A4 | 21 | 2 | 19 | 90.5% | 78.7% | +11.8% |

**Maybe:** CypReact non-reactive calls are above the baseline TDI-negative rates, but there are false-positives.  So a negative call by CypReact is, perhaps, slightly informative for negative-TDI.  

## Limitations

- CypReact predicts metabolic reactant status, whereas the challenge endpoint is time-dependent inhibition.
- The public CypReact interface provides only Yes/No predictions, so we cannot assess confidence or do any sort of ranking within the large `R` group.

## Reproduction

From the repository root, with the local CypReact v1.2 bundle and cheminformatics environment available:

```bash
/home/dan/swr/miniconda3/envs/cheminf/bin/python submissions/classification/002-cypreact-raw/create_submission.py --timestamp 20260821T021607Z
```

The expected output is `submissions/classification/002-cypreact-raw/openadmet-cyp_classification-sub_dargason_20260821T021607Z.csv`. See the [full experiment record](../../../experiments/20260820_public-cyp-tdi-models/README.md) for model evaluation, provenance, and supporting artifacts.
