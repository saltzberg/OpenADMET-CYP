# Public CYP models versus OpenADMET TDI

This subproject measures whether public CYP substrate, inhibitor, and reaction scores rank OpenADMET time-dependent inhibition (TDI) positives for the matching isoform.

It does not retrain the public models and does not reinterpret their outputs as calibrated TDI probabilities. Substrate recognition, reversible inhibition, predicted metabolism, and TDI are related but distinct claims.

## Comparisons

The registry defines ten endpoint-matched comparisons:

| Public model | Score | CYP2D6 TDI | CYP3A4 TDI |
|---|---|---:|---:|
| DeepP450 | substrate | yes | yes |
| CYPMol | substrate | yes | yes |
| CYPMol | inhibitor | yes | yes |
| DeepMetab | substrate | yes | yes |
| CyProduct / CypReact | public reactant score (hard-vote fraction or R/N) | yes | yes |

`config/model_registry.csv` is the source of truth. Rows remain `planned` until an upstream model revision, score column, score direction, and prediction artifact have been verified. The evaluator runs only rows marked `ready`.

Completed comparisons:

- [`reports/deepmetab_v1.md`](reports/deepmetab_v1.md)
- [`reports/cypreact_v1.md`](reports/cypreact_v1.md)
- [`reports/cypmol_v1.md`](reports/cypmol_v1.md)

## Training labels

The source is the released OpenADMET file:

`../data/cyp-challenge-train-test/cyp-challenge-TRAIN_TDI.csv`

Prepared label coverage is:

| Endpoint | Labels | TDI positive | TDI negative | Positive prevalence |
|---|---:|---:|---:|---:|
| CYP2D6 | 1,497 | 324 | 1,173 | 21.64% |
| CYP3A4 | 3,584 | 764 | 2,820 | 21.32% |

There are 4,822 unique molecules across either endpoint. The generated label registry contains no invalid SMILES.

Prepare it with the existing `cheminf` environment:

```bash
cd /home/dan/projects/ADMET-CYP/tdi_public_model_benchmark

PYTHONPATH=src /home/dan/swr/miniconda3/envs/cheminf/bin/python \
  -m tdi_public_model_benchmark prepare-labels \
  --source ../data/cyp-challenge-train-test/cyp-challenge-TRAIN_TDI.csv \
  --output-dir data
```

This writes ignored, reproducible artifacts:

- `data/tdi_labels.csv`: one row per labelled molecule × endpoint;
- `data/model_input.csv`: one row per unique molecule for model inference;
- `data/label_manifest.json`: source/output hashes and label/scaffold counts.

## Prediction contract

Each registry row points to a CSV and names:

- `id_column`: normally `Molecule_Name`;
- `score_column`: the exact upstream score to assess;
- `score_direction`: `higher` if larger means more likely substrate/inhibitor/reaction, otherwise `lower`;
- `source_revision`: immutable commit, release, DOI/version, or model-asset hash;
- `status`: `planned` or `ready`.

Prediction files may cover only part of an endpoint. Blank scores are retained as missing and reduce coverage. Scores are never imputed. Duplicate molecule IDs, nonnumeric nonblank values, infinite values, unknown endpoints, and duplicate registry keys fail explicitly.

Upstream output should be normalized by joining back to `data/model_input.csv` through a stable row identifier. Do not rely on output row order alone.

## Metrics

### Primary: precision–recall AUC

The primary statistic is average precision, computed with scikit-learn's `average_precision_score` and reported as `pr_auc_average_precision`. Naming the implementation avoids ambiguity with trapezoidal interpolation of a precision–recall curve.

### Secondary

- `roc_auc`: scikit-learn ROC-AUC;
- top 5%, 10%, and 20% enrichment;
- prediction coverage.

For top-fraction enrichment:

```text
k = ceil(fraction × number of scored labelled compounds)
enrichment = positive rate among the top k / positive prevalence among scored compounds
```

Scores are oriented so larger always means stronger model evidence, then sorted descending. If the top-k boundary cuts through an exact score tie, tied compounds receive fractional expected inclusion; this avoids arbitrary enrichment driven by molecule-ID ordering for coarse scores such as CypReact R/N. The output reports `k`, expected positive hits, boundary tie size/fraction, top-set positive rate, scored-set prevalence, and enrichment. Coverage is finite scored labelled compounds divided by all labelled compounds for that endpoint.

Coverage is operational, not evidence of model quality. PR-AUC/ROC-AUC/enrichment are conditional on the compounds the public model successfully scores, while coverage shows how broad that conditional result is.

## Scaffold-bootstrap confidence intervals

All requested metrics receive percentile confidence intervals from a cluster bootstrap over the full labelled endpoint universe:

1. assign each valid cyclic molecule its achiral Bemis–Murcko scaffold;
2. assign acyclic molecules deterministic molecular-identity singleton groups rather than putting all acyclic chemistry into one artificial cluster;
3. sample scaffold groups with replacement;
4. retain all molecules in each sampled group, including missing model predictions;
5. recompute coverage and score metrics for each replicate.

The default is 2,000 replicates, 95% intervals, and seed `20260820`. Replicates lacking both classes among scored compounds remain undefined for PR-AUC/ROC-AUC and are excluded only from that metric's percentile calculation. `valid_bootstrap_replicates` makes this visible.

## Run the benchmark

### DeepMetab inference

DeepMetab runs in an isolated Python 3.8 / ChemProp 1.5.2 environment. The
adapter preserves molecule IDs, all five checkpoint scores, the ensemble mean
and population variance for all nine CYP tasks, plus a failure ledger and
hash-bound manifest:

```bash
PYTHONPATH=src /home/dan/swr/miniconda3/envs/deepmetab-cyp/bin/python \
  -m tdi_public_model_benchmark.deepmetab \
  --input data/model_input.csv \
  --checkpoint-dir vendor/DeepMetab/Model/Substrate/Multi \
  --output predictions/raw/deepmetab.csv \
  --manifest predictions/raw/deepmetab.manifest.json \
  --failures predictions/raw/deepmetab.failures.csv \
  --batch-size 64
```

The public DeepMetab checkpoints name their output head `readout.*`, while
compatible ChemProp 1.x calls the same layers `ffn.*`. The upstream
`load_checkpoint` helper silently skips those unmatched head weights. This
adapter maps `readout.*` to `ffn.*`, verifies every key and tensor shape, and
requires a strict full-state load before inference.

### Registry evaluation

After model outputs have been verified and corresponding registry rows changed to `ready`:

```bash
PYTHONPATH=src /home/dan/swr/miniconda3/envs/cheminf/bin/python \
  -m tdi_public_model_benchmark evaluate-registry \
  --labels data/tdi_labels.csv \
  --registry config/model_registry.csv \
  --output-dir outputs/public_models_v1 \
  --bootstrap-replicates 2000 \
  --seed 20260820
```

The output directory contains:

- `metrics.csv`: point estimates and confidence intervals in tidy form;
- `enrichment.csv`: detailed top-fraction counts and enrichment intervals;
- `bootstrap_metrics.csv`: every bootstrap replicate;
- `joined_predictions.csv`: full labels joined to raw and oriented scores, including missingness;
- `registry_status.csv`: evaluated and skipped rows;
- `run_manifest.json`: input hashes, model provenance, metric definitions, and output hashes.

Existing nonempty output directories are protected unless `--overwrite` is supplied.

### CypReact inference

CypReact v1.2 is run from the public historical Java bundle. The adapter writes
an ID-labelled SDF, invokes the unmodified JAR for CYP2D6 and CYP3A4, and maps
the released `R`/`N` calls to 1/0 scores without claiming probability
calibration:

```bash
PYTHONPATH=src /home/dan/swr/miniconda3/envs/cheminf/bin/python \
  -m tdi_public_model_benchmark.cypreact \
  --input data/model_input.csv \
  --bundle-dir vendor/CypReact/CypReactBundle \
  --output predictions/raw/cypreact.csv \
  --manifest predictions/raw/cypreact.manifest.json \
  --failures predictions/raw/cypreact.failures.csv \
  --work-dir predictions/raw/cypreact_work
```

The public CLI exports only the final reactant call. The resulting PR-AUC,
ROC-AUC and enrichment comparisons are therefore coarse binary-score tests
with many ties, not evaluations of a continuous calibrated score.

### CYPMol inference

CYPMol runs in an isolated Blackwell-compatible environment. Substrate and
inhibitor tasks are run separately so only one approximately 2.81 GB
checkpoint is resident on the GPU at a time:

```bash
PYTHONPATH=src /home/dan/swr/miniconda3/envs/cypmol/bin/python \
  -m tdi_public_model_benchmark.cypmol \
  --task substrate \
  --input data/model_input.csv \
  --source-repo vendor/CYPMol \
  --assets-dir vendor/CYPMol_assets \
  --dictionary vendor/CYPMol_assets/token_list.txt \
  --output predictions/raw/cypmol_substrate.csv \
  --manifest predictions/raw/cypmol_substrate.manifest.json \
  --failures predictions/raw/cypmol.failures.csv \
  --work-dir predictions/raw/cypmol_work \
  --batch-size 6

PYTHONPATH=src /home/dan/swr/miniconda3/envs/cypmol/bin/python \
  -m tdi_public_model_benchmark.cypmol \
  --task inhibitor \
  --input data/model_input.csv \
  --source-repo vendor/CYPMol \
  --assets-dir vendor/CYPMol_assets \
  --dictionary vendor/CYPMol_assets/token_list.txt \
  --output predictions/raw/cypmol_inhibitor.csv \
  --manifest predictions/raw/cypmol_inhibitor.manifest.json \
  --failures predictions/raw/cypmol.failures.csv \
  --work-dir predictions/raw/cypmol_work \
  --batch-size 6
```

The adapter uses strict complete state loading, preserves all ten checkpoint
scores plus ensemble means/variances, and caches deterministic conformers.
Two published-recipe conformer failures remain missing rather than receiving a
fallback representation.

## Verification

```bash
PYTHONPATH=src /home/dan/swr/miniconda3/envs/cheminf/bin/python -m pytest -q
```

The test suite covers score orientation, average precision, ROC-AUC, enrichment, missing-score coverage, duplicate-ID failure, scaffold bootstrap, label preparation, and an end-to-end registry run. Synthetic predictions exist only in temporary test directories and are never represented as public-model results.

## Interpretation boundaries

- A strong substrate score association may reflect metabolism but does not establish mechanism-based inactivation.
- An inhibitor score may be closer to the assay phenotype but can predominantly encode reversible inhibition.
- A reaction or site-of-metabolism score may capture turnover without capturing a reactive intermediate or durable enzyme inactivation.
- Different public models may have been trained on overlapping databases. Their results are comparisons, not independent experimental evidence.
- Missing scores can be chemically structured. Always read ranking metrics together with coverage and the joined missingness table.

See `docs/model_sources.md` for model-specific provenance and execution notes.
