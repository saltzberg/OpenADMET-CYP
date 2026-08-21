# CYPMol substrate and inhibitor scores versus OpenADMET TDI

## Purpose

Assess whether the public CYPMol substrate or inhibitor ensemble scores rank
OpenADMET CYP2D6 and CYP3A4 time-dependent inhibition (TDI) positives. The
published models are used without fitting or threshold tuning against TDI.

## Result

None of the four CYPMol score/endpoint comparisons provides useful positive TDI
prioritization. CYP2D6 substrate is approximately random. CYP2D6 inhibitor and
CYP3A4 substrate are inversely associated with TDI. CYP3A4 inhibitor has a
small positive global ROC-AUC but materially depletes TDI positives in the
highest-scoring tail.

| Score | Endpoint | Scored | Coverage | Average precision (95% CI) | ROC-AUC (95% CI) |
|---|---|---:|---:|---:|---:|
| Substrate | CYP2D6 | 1,497 | 1.0000 | 0.2297 (0.2005–0.2668) | 0.5104 (0.4740–0.5455) |
| Inhibitor | CYP2D6 | 1,497 | 1.0000 | 0.1962 (0.1741–0.2251) | 0.4687 (0.4355–0.5025) |
| Substrate | CYP3A4 | 3,582 | 0.9994 | 0.1865 (0.1708–0.2035) | 0.4773 (0.4538–0.5009) |
| Inhibitor | CYP3A4 | 3,582 | 0.9994 | 0.2069 (0.1900–0.2255) | 0.5321 (0.5090–0.5541) |

The no-skill average-precision references are 0.2164 for CYP2D6 and 0.2132 for
CYP3A4.

## Top-fraction enrichment

| Score | Endpoint | Top 5% | Top 10% | Top 20% |
|---|---|---:|---:|---:|
| Substrate | CYP2D6 | 0.862 (0.541–1.380) | 1.016 (0.741–1.327) | 1.032 (0.819–1.208) |
| Inhibitor | CYP2D6 | 0.616 (0.296–0.975) | 0.770 (0.482–1.030) | 0.755 (0.588–0.947) |
| Substrate | CYP3A4 | 0.208 (0.076–0.364) | 0.261 (0.149–0.386) | 0.504 (0.383–0.626) |
| Inhibitor | CYP3A4 | 0.443 (0.236–0.660) | 0.522 (0.366–0.690) | 0.745 (0.607–0.868) |

The CYP3A4 inhibitor score illustrates why global ROC-AUC and top-tail
prioritization must be reported separately: its overall pair ordering is
slightly positive, but its highest scores are strongly depleted in TDI
positives.

## Inputs and inference

- CYPMol Git revision: `0c6657dc882f23293ff9a43977c8f31818f8a6f8`.
- Hugging Face model revision: `3810a93e04d1dbab4cde73792a1e3830b1f17e41`.
- Ten substrate checkpoints and ten inhibitor checkpoints, approximately 56.2
  GB total.
- Model score: mean positive-class score over the ten task-specific checkpoints.
- Input: 4,822 unique molecules using original supplied SMILES.
- Molecule preparation: authors' AddHs, ETKDG-style embedding with random seed
  42 and random coordinates, MMFF optimization, removal of hydrogens, coordinate
  centering, and Uni-Mol tokenization.
- Protein context: published CYP2D6/CYP3A4 sequences, pocket residues, and (for
  substrate) substrate-recognition-site intervals.
- Runtime: Python 3.10.20, PyTorch 2.11.0+cu128, fair-esm 2.0.0, RDKit 2026.03.2,
  NVIDIA RTX 5070 Ti.
- Execution: strict full-state checkpoint loading, one checkpoint resident on
  the GPU at a time.

Two molecules failed the published conformer-generation recipe and were not
rescued or imputed:

- `OCNT-0453746`
- `OCNT-0453782`

Both are CYP3A4-labelled, producing 3,582/3,584 endpoint coverage. Neither is
in the CYP2D6 label set.

## Reproducibility repairs

The public repository is not directly runnable: it contains hard-coded author
paths, requires missing token dictionaries and base-model initialization files,
and has no ensemble driver. Each released checkpoint contains the complete
Uni-Mol, ESM-2, fusion, and classifier state, so the adapter constructs the
published architecture without external initialization weights and requires
exact key and tensor-shape parity before loading.

The original inference repeatedly encodes identical CYP sequences for every
molecule. The adapter computes the endpoint-specific protein context once per
checkpoint, then reuses it for every molecular batch. This is algebraically
identical in evaluation mode for fixed endpoint batches and avoids redundant
ESM-2 work.

The public web API was also retested and remained unavailable: both substrate
and inhibitor endpoints returned HTTP 500 because their backend connection was
refused.

## Uncertainty

Confidence intervals are 2,000-replicate percentile intervals from a
Bemis–Murcko scaffold-cluster bootstrap over the full labelled endpoint
universe. Acyclic molecules use deterministic molecular-identity singleton
groups. Seed: `20260820`.

## Artifacts

- `predictions/raw/cypmol_substrate.csv`
- `predictions/raw/cypmol_substrate.manifest.json`
- `predictions/raw/cypmol_inhibitor.csv`
- `predictions/raw/cypmol_inhibitor.manifest.json`
- `predictions/raw/cypmol.failures.csv`
- `outputs/public_models_v4/metrics.csv`
- `outputs/public_models_v4/enrichment.csv`
- `outputs/public_models_v4/bootstrap_metrics.csv`
- `outputs/public_models_v4/joined_predictions.csv`
- `outputs/public_models_v4/run_manifest.json`

## Interpretation boundary

This result tests transfer of published substrate/inhibitor rankings to a
different TDI assay endpoint. It does not assess CYPMol on its own training
claims. The GitHub source repository has no detected license; the Hugging Face
model repository declares Apache-2.0. Keep this internal experiment separate
from redistribution or commercial integration until source-code terms are
clarified.
