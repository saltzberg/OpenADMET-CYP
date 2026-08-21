# DeepMetab substrate score versus OpenADMET TDI

## Purpose

Assess whether the public DeepMetab substrate classifier ranks time-dependent
inhibitors for the matching CYP2D6 or CYP3A4 endpoint. This is an external-score
association test, not TDI-model training.

## Result

DeepMetab substrate scores do not provide useful positive TDI enrichment in
this dataset. CYP2D6 is indistinguishable from random ranking. CYP3A4 orders
TDI positives slightly worse than random, including significant depletion in
the top 20% by the scaffold-bootstrap interval.

| Endpoint | Labels | TDI prevalence | Coverage | Average precision (95% CI) | ROC-AUC (95% CI) |
|---|---:|---:|---:|---:|---:|
| CYP2D6 | 1,497 | 0.2164 | 1.000 | 0.2140 (0.1879–0.2457) | 0.4917 (0.4562–0.5260) |
| CYP3A4 | 3,584 | 0.2132 | 1.000 | 0.1974 (0.1794–0.2178) | 0.4635 (0.4391–0.4877) |

Average precision is the primary precision–recall statistic. The no-skill
reference is the endpoint's positive prevalence: 0.2164 for CYP2D6 and 0.2132
for CYP3A4.

## Enrichment

| Endpoint | Fraction | Ranked molecules | TDI positives | Top positive rate | Enrichment (95% CI) |
|---|---:|---:|---:|---:|---:|
| CYP2D6 | 5% | 75 | 17 | 0.2267 | 1.047 (0.650–1.507) |
| CYP2D6 | 10% | 150 | 30 | 0.2000 | 0.924 (0.688–1.229) |
| CYP2D6 | 20% | 300 | 63 | 0.2100 | 0.970 (0.777–1.165) |
| CYP3A4 | 5% | 180 | 37 | 0.2056 | 0.964 (0.719–1.259) |
| CYP3A4 | 10% | 359 | 67 | 0.1866 | 0.875 (0.698–1.057) |
| CYP3A4 | 20% | 717 | 127 | 0.1771 | 0.831 (0.703–0.945) |

## Inputs and inference

- DeepMetab revision: `78c7511327f1a4042b61a64c44d35abb4c4b6b9c`
- Input: 4,822 unique TDI-labelled-union molecules using original supplied SMILES
- Checkpoints: five public multitask substrate classifiers
- Score: unthresholded mean of the five classification scores
- Inference coverage: 4,822/4,822; zero failures
- Runtime: Python 3.8.20, ChemProp 1.5.2, PyTorch 2.0.1 CPU, RDKit 2022.09.5
- Prediction SHA-256: `8552da25bfd5a1b3460fa66fc59d5698f8adefa2a2c09012fcabfdea53dd15fe`

The public checkpoint state dictionaries name their trained output layers
`readout.*`; compatible ChemProp 1.x names those layers `ffn.*`. Unmodified
ChemProp loading silently skips the trained output head. The benchmark adapter
maps those names, verifies exact key and tensor-shape parity, and performs a
strict full-state load.

## Uncertainty

Confidence intervals are 2,000-replicate percentile intervals from a
Bemis–Murcko scaffold-cluster bootstrap over the full labelled endpoint
universe. Acyclic molecules use deterministic molecular-identity singleton
groups. Seed: `20260820`.

## Artifacts

Generated local artifacts:

- `predictions/raw/deepmetab.csv`: all five checkpoint scores, ensemble means,
  and ensemble variances for all nine CYP tasks
- `predictions/raw/deepmetab.manifest.json`: checkpoint, runtime, input, and
  output provenance
- `predictions/raw/deepmetab.failures.csv`: empty failure ledger
- `outputs/deepmetab_v1/metrics.csv`
- `outputs/deepmetab_v1/enrichment.csv`
- `outputs/deepmetab_v1/bootstrap_metrics.csv`
- `outputs/deepmetab_v1/joined_predictions.csv`
- `outputs/deepmetab_v1/run_manifest.json`

## Interpretation boundary

This result says that DeepMetab's substrate ranking does not transfer to the
OpenADMET TDI label. It does not imply that DeepMetab is a poor substrate
classifier, nor that substrate recognition is unrelated to TDI mechanistically.
The endpoints, assay construction, training sources, and class definitions are
different. The score should remain an external comparator rather than a TDI
feature promoted on prior plausibility alone.
