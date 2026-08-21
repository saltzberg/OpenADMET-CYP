# 002-cypreact-raw

This submission uses the frozen CypReact v1.2 endpoint models, maps their public `R` (reactant) and `N` (non-reactant) calls directly to `True` and `False`, and outputs CYP2D6 and CYP3A4 TDI classifications.

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

- **Base model:** CypReact v1.2, using the released CYP2D6 and CYP3A4 models and bundled support files.
- **Prediction interface:** the frozen public binary `R`/`N` interface; it exposes neither probabilities nor confidence scores.
- **Challenge-data use:** no challenge labels were used for training, selection, calibration, threshold fitting, or post-processing.

The adapter preserves blinded-test row order and identifiers and fails closed if either endpoint is missing or invalid. Exact inputs, runtime, model assets, commands, and checksums are recorded in [`inference_manifest.json`](inference_manifest.json).

## Results and discussion

Association with released OpenADMET training labels was evaluated separately from the blinded submission using 2,000-replicate Bemis–Murcko scaffold-cluster bootstrap intervals. These values describe released-training association, not blind or leaderboard performance.

| Endpoint | Coverage | Average precision (95% CI) | ROC-AUC (95% CI) | Reactant calls |
|---|---:|---:|---:|---:|
| CYP2D6 | 1.000 | 0.219 (0.198–0.241) | 0.508 (0.497–0.518) | 96.3% |
| CYP3A4 | 1.000 | 0.214 (0.199–0.230) | 0.502 (0.500–0.504) | 99.4% |

Average precision was close to the label prevalence (0.216 for CYP2D6 and 0.213 for CYP3A4), while ROC-AUC was near 0.5. The binary output was also strongly saturated toward `R`. Together, these released-training results indicate near-random TDI discrimination and little ability to prioritize compounds, rather than evidence about blind performance. Exact metrics are in the [experiment table](../../../experiments/20260820_public-cyp-tdi-models/artifacts/metrics.csv).

## Limitations

- CypReact predicts metabolic reactant status, whereas the challenge endpoint is time-dependent inhibition. The direct mapping is therefore a transfer baseline, not a claim that reactant status and TDI are biologically equivalent.
- The public interface provides only saturated binary calls, so it cannot express confidence or support useful ranking within the large `R` group.
- The historical Java/Weka runtime and model bundle constrain portability.

## Reproduction

From the repository root, with the local CypReact v1.2 bundle and cheminformatics environment available:

```bash
/home/dan/swr/miniconda3/envs/cheminf/bin/python submissions/classification/002-cypreact-raw/create_submission.py --timestamp 20260821T021607Z
```

The expected output is `submissions/classification/002-cypreact-raw/openadmet-cyp_classification-sub_dargason_20260821T021607Z.csv`. See the [full experiment record](../../../experiments/20260820_public-cyp-tdi-models/README.md) for model evaluation, provenance, and supporting artifacts.
