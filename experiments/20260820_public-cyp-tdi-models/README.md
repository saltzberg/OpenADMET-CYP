# Public CYP models versus OpenADMET TDI

> **Status:** Completed current scope — DeepP450 deferred
> **Intent date (UTC):** 20260820
> **Experiment ID:** `20260820_public-cyp-tdi-models`

## Summary

Test whether public CYP substrate, inhibitor, or reactant scores rank OpenADMET
CYP2D6 and CYP3A4 time-dependent inhibition (TDI) positives. These external
scores are compared with TDI labels without retraining or threshold tuning.

The completed models do not provide useful TDI prioritization.
DeepMetab is random for CYP2D6 and slightly inverse for CYP3A4. CypReact's
public R/N output calls almost every molecule a reactant. CYPMol's substrate
and inhibitor scores are random or inverse in the high-scoring TDI tail.

## Procedure

1. Build a long-form label table from the released OpenADMET TDI training data.
2. Preserve the original supplied SMILES for public-model inference. The broader
   feature store retains canonical and enumerated molecular states for later
   sensitivity work.
3. Run each public model without fitting to TDI labels. Preserve raw component
   scores, ensemble aggregation, failures, runtime, revisions, and hashes.
4. Evaluate each endpoint separately using average precision as the primary
   precision–recall statistic, ROC-AUC, top 5%/10%/20% enrichment, and coverage.
5. Estimate 95% percentile intervals with a 2,000-replicate Bemis–Murcko
   scaffold-cluster bootstrap over the full labelled endpoint universe.

## DeepMetab implementation

- Upstream revision: `78c7511327f1a4042b61a64c44d35abb4c4b6b9c`.
- Five public multitask ChemProp substrate checkpoints.
- Score: unthresholded mean of the five checkpoint classification scores.
- Runtime: Python 3.8.20, ChemProp 1.5.2, PyTorch 2.0.1 CPU, RDKit 2022.09.5.
- Inference: 4,822/4,822 molecules; zero failures.

The released checkpoint output layers are named `readout.*`, while compatible
ChemProp 1.x builds `ffn.*`. Unmodified ChemProp loading silently skips the
trained output head. The experiment adapter maps the names, verifies complete
key and tensor-shape parity, then requires a strict state load.

## CypReact implementation

- Public artifact: historical runnable CypReact v1.2 bundle.
- Release ZIP SHA-256:
  `2778b65dc786125a9f6af489d599f6de6bbe81952f3d484285d0dc73370da48f`.
- Input: original supplied SMILES converted to an ID-labelled RDKit SDF.
- Output: faithful public `R`/`N` call encoded `R=1`, `N=0`; no continuous
  probability or CYP2D6 vote fraction is exposed.
- Runtime: OpenJDK 1.8.0_492.
- Returned: 4,821/4,822 molecules. `OCNT-2328840` was omitted by the public
  tool and retained as missing.

## CYPMol implementation

- Git revision: `0c6657dc882f23293ff9a43977c8f31818f8a6f8`.
- Hugging Face revision: `3810a93e04d1dbab4cde73792a1e3830b1f17e41`.
- Ten substrate and ten inhibitor checkpoints; approximately 56.2 GB total.
- Strict full-state loading on an RTX 5070 Ti with PyTorch 2.11.0+cu128.
- Original supplied SMILES with the authors' AddHs, random-seed-42 embedding,
  MMFF, heavy-atom, centered-coordinate Uni-Mol preparation.
- Fixed published CYP2D6/CYP3A4 sequence and pocket contexts; substrate also
  uses published substrate-recognition-site intervals.
- Two conformer failures (`OCNT-0453746`, `OCNT-0453782`) were not rescued or
  imputed. Both affect CYP3A4 coverage only.

## Results

| Model | Endpoint | Labels | Coverage | Average precision (95% CI) | ROC-AUC (95% CI) |
|---|---|---:|---:|---:|---:|
| DeepMetab | CYP2D6 | 1,497 | 1.0000 | 0.2140 (0.1879–0.2457) | 0.4917 (0.4562–0.5260) |
| DeepMetab | CYP3A4 | 3,584 | 1.0000 | 0.1974 (0.1794–0.2178) | 0.4635 (0.4391–0.4877) |
| CypReact | CYP2D6 | 1,497 | 1.0000 | 0.2191 (0.1977–0.2415) | 0.5077 (0.4971–0.5177) |
| CypReact | CYP3A4 | 3,584 | 0.9997 | 0.2139 (0.1994–0.2300) | 0.5021 (0.5004–0.5037) |
| CYPMol substrate | CYP2D6 | 1,497 | 1.0000 | 0.2297 (0.2005–0.2668) | 0.5104 (0.4740–0.5455) |
| CYPMol inhibitor | CYP2D6 | 1,497 | 1.0000 | 0.1962 (0.1741–0.2251) | 0.4687 (0.4355–0.5025) |
| CYPMol substrate | CYP3A4 | 3,584 | 0.9994 | 0.1865 (0.1708–0.2035) | 0.4773 (0.4538–0.5009) |
| CYPMol inhibitor | CYP3A4 | 3,584 | 0.9994 | 0.2069 (0.1900–0.2255) | 0.5321 (0.5090–0.5541) |

The no-skill average-precision references are the label prevalences: 0.2164 for
CYP2D6 and 0.2132 for CYP3A4. CypReact called 96.3% of scored CYP2D6-labelled
molecules and 99.4% of scored CYP3A4-labelled molecules reactants. Its high
TDI-positive recall therefore comes with almost no specificity.

For tied binary CypReact scores, top-fraction enrichment uses fractional
expected inclusion at the R/N boundary rather than arbitrary molecule-ID
ordering. Expected lift was 1.013 for CYP2D6 and 1.003 for CYP3A4 at each of the
top 5%, 10%, and 20% fractions: numerically near one and not useful for
prioritization.

CYPMol CYP3A4 inhibitor shows why global ROC-AUC is insufficient: ROC-AUC is
slightly above 0.5, but enrichment is only 0.443 in the top 5%, 0.522 in the
top 10%, and 0.745 in the top 20%. Its highest scores are depleted in TDI
positives. CYPMol CYP3A4 substrate is more strongly depleted: 0.208, 0.261,
and 0.504 at the same fractions.

## Model queue

- DeepMetab substrate score: completed.
- CypReact public reactant output: completed; almost-all-positive and not useful
  for TDI prioritization.
- CYPMol substrate and inhibitor scores: completed; no useful positive TDI
  enrichment.
- DeepP450 substrate score: deferred; the approximately 58.8 GB Baidu-hosted
  bundle requires authenticated download. The inference adapter is prepared.

## Outputs

- `index.html` — compact standalone procedure and results dashboard.
- `artifacts/metrics.csv` — point estimates and confidence intervals.
- `artifacts/enrichment.csv` — top-fraction enrichment counts and intervals.
- `artifacts/bootstrap_metrics.csv` — all 2,000 scaffold-bootstrap replicates.
- `artifacts/joined_predictions.csv` — endpoint labels joined to both model scores.
- `artifacts/deepmetab_predictions.csv` — all five checkpoint scores, means, and
  variances for all nine CYP tasks.
- `artifacts/cypreact_predictions.csv` — public CYP2D6/CYP3A4 R/N calls and
  encoded binary scores.
- `artifacts/cypreact_failures.csv` — explicit one-molecule omission ledger.
- `artifacts/inference_manifest.json` — checkpoint/runtime/input provenance.
- `artifacts/cypreact_inference_manifest.json` — CypReact release, assets,
  runtime, input, output, and failure provenance.
- `artifacts/cypmol_substrate_predictions.csv` and
  `artifacts/cypmol_inhibitor_predictions.csv` — all checkpoint scores,
  ensemble means, and variances.
- `artifacts/cypmol_substrate_manifest.json` and
  `artifacts/cypmol_inhibitor_manifest.json` — revisions, checkpoint hashes,
  endpoint contexts, runtime, and output provenance.
- `artifacts/cypmol_failures.csv` — the two conformer-generation failures.
- `artifacts/benchmark_run_manifest.json` — metric and output provenance.
- `artifacts/model_registry.csv` — completed and planned comparisons.
- `artifacts/run_manifest.json` — experiment bundle hashes and source mapping.

## Reproduction

The canonical inference/evaluation implementation remains in
`tdi_public_model_benchmark/`. Regenerate this experiment bundle and verify it:

```bash
/home/dan/swr/miniconda3/envs/cheminf/bin/python \
  experiments/20260820_public-cyp-tdi-models/src/build_dashboard.py

/home/dan/swr/miniconda3/envs/cheminf/bin/python \
  experiments/20260820_public-cyp-tdi-models/src/verify_outputs.py
```

This is an experiment record, not a challenge submission.
