# CypReact public reactant call versus OpenADMET TDI

## Purpose

Assess whether the faithful public CypReact v1.2 reactant classification ranks
OpenADMET CYP2D6 or CYP3A4 time-dependent inhibition (TDI) positives. The
released CLI exports only `R`/`N`, so the benchmark encodes `R=1`, `N=0` and
does not claim a continuous probability.

## Result

CypReact is effectively an almost-all-positive classifier on this compound
universe. It has full CYP2D6 label coverage and misses one CYP3A4-labelled
compound. Its PR-AUC is at prevalence and ROC-AUC is approximately 0.5 for both
endpoints. The public R/N output does not provide useful TDI discrimination.

| Endpoint | Labels | Scored | Coverage | Average precision (95% CI) | ROC-AUC (95% CI) |
|---|---:|---:|---:|---:|---:|
| CYP2D6 | 1,497 | 1,497 | 1.0000 | 0.2191 (0.1977–0.2415) | 0.5077 (0.4971–0.5177) |
| CYP3A4 | 3,584 | 3,583 | 0.9997 | 0.2139 (0.1994–0.2300) | 0.5021 (0.5004–0.5037) |

The no-skill average-precision references are 0.2164 for CYP2D6 and 0.2132 for
CYP3A4.

## Public-call behavior on labelled molecules

| Endpoint | TDI positive called R | TDI positive called N | TDI negative called R | TDI negative called N |
|---|---:|---:|---:|---:|
| CYP2D6 | 316 | 8 | 1,126 | 47 |
| CYP3A4 | 762 | 2 | 2,800 | 19 |

Thus CypReact calls 96.3% of scored CYP2D6-labelled molecules and 99.4% of
scored CYP3A4-labelled molecules reactants. This high call rate explains why
ranking metrics remain near random despite high TDI-positive recall.

## Tie-aware enrichment

The public score has only two values. A top-k cutoff therefore falls inside a
large tied `R` group. Enrichment uses fractional expected inclusion within the
boundary tie rather than arbitrary molecule-ID ordering.

| Endpoint | Top fraction | Expected TDI positives | Enrichment (95% CI) | Boundary tied R calls |
|---|---:|---:|---:|---:|
| CYP2D6 | 5% | 16.44/75 | 1.013 (0.995–1.029) | 1,442 |
| CYP2D6 | 10% | 32.87/150 | 1.013 (0.995–1.029) | 1,442 |
| CYP2D6 | 20% | 65.74/300 | 1.013 (0.995–1.029) | 1,442 |
| CYP3A4 | 5% | 38.51/180 | 1.003 (1.001–1.006) | 3,562 |
| CYP3A4 | 10% | 76.80/359 | 1.003 (1.001–1.006) | 3,562 |
| CYP3A4 | 20% | 153.38/717 | 1.003 (1.001–1.006) | 3,562 |

The numerically narrow CYP3A4 interval reflects the enormous tied R group; the
lift is only 0.3% and is not practically useful for prioritization.

## Inputs and inference

- Public release: CypReact v1.2 historical runnable bundle.
- Release ZIP SHA-256: `2778b65dc786125a9f6af489d599f6de6bbe81952f3d484285d0dc73370da48f`.
- Runtime: OpenJDK 1.8.0_492; ID-preserving RDKit SDF adapter.
- Input: 4,822 unique molecules using original supplied SMILES.
- Returned: 4,821 molecules with both CYP calls.
- Failure: `OCNT-2328840` was omitted by the public tool and was not imputed.
- Prediction SHA-256: `e364e4a38380d52862877153cc94b85d678294372673c493979c3c6b4c99a47d`.

Across the complete 4,821 returned molecules, CypReact emitted 4,602 CYP2D6 R
calls and 4,790 CYP3A4 R calls.

## Uncertainty

Confidence intervals are 2,000-replicate percentile intervals from a
Bemis–Murcko scaffold-cluster bootstrap over the full labelled endpoint
universe. Acyclic molecules use deterministic molecular-identity singleton
groups. Seed: `20260820`.

## Artifacts

- `predictions/raw/cypreact.csv`
- `predictions/raw/cypreact.failures.csv`
- `predictions/raw/cypreact.manifest.json`
- `predictions/raw/cypreact_work/cypreact_input.sdf`
- `predictions/raw/cypreact_work/cypreact_output.csv`
- `predictions/raw/cypreact_work/cypreact.log`
- `outputs/public_models_v3/metrics.csv`
- `outputs/public_models_v3/enrichment.csv`
- `outputs/public_models_v3/bootstrap_metrics.csv`
- `outputs/public_models_v3/joined_predictions.csv`
- `outputs/public_models_v3/run_manifest.json`

## Interpretation boundary

This result evaluates the released binary reactant call, not an unexposed Weka
class distribution or CYP2D6 hard-vote fraction. It does not test whether a
continuous internal CypReact score might contain more ranking information. The
public output itself is too saturated to prioritize TDI compounds.
