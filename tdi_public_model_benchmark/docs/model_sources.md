# Public model source and execution audit

Reviewed 2026-08-20. Repositories were cloned read-only to a temporary directory and inspected at the revisions below. No upstream source or model artifact has been copied into this project.

## Status summary

| Model | Authoritative source | Inspected revision | Continuous score available? | Current benchmark status |
|---|---|---|---|---|
| DeepP450 | https://github.com/CjmTH/DeepP450 | `dea269a228177491aa02b02f5f3c396f415778ac` | Yes in code: class-1 softmax, intended 20-model mean | Blocked on ~58.8 GB Baidu weight retrieval, repaired inference, and license clarification |
| CYPMol | https://github.com/CjmTH/CYPMol and https://huggingface.co/CJM1111/CYPMOl | Git `0c6657dc882f23293ff9a43977c8f31818f8a6f8`; HF `3810a93e04d1dbab4cde73792a1e3830b1f17e41` | Yes: class-1 softmax for substrate and inhibitor heads | Completed for CYP2D6/CYP3A4 with both ten-checkpoint ensembles |
| DeepMetab | https://github.com/YilingZhou/DeepMetab | `78c7511327f1a4042b61a64c44d35abb4c4b6b9c` | Yes internally: five-model mean score and variance | Runnable after Git-LFS retrieval and a ChemProp 1.x environment |
| CypReact used by CyProduct | https://github.com/Le0nT1/CyProduct and https://github.com/Le0nT1/CypReact | CyProduct `3fb35435cbdb2a9e1e039b59492ebc0e53dc564a`; CypReact `eb5925f364ecb3b3a4dbdf603b4abdb11a4c8c3f` | Public output is hard R/N; CYP2D6 internally has a five-model hard-vote fraction | Runnable historical bundle, but score resolution is coarse and CyProduct licensing is unresolved |

“Public” here means the source or assets are publicly reachable. It does not by itself resolve software licensing or guarantee that the published command is reproducible on a current environment.

## DeepP450

Purpose: protein-sequence-conditioned prediction of whether a molecule is a substrate/activity-positive for a P450.

Paper: Chang, Fan, and Tian, *J. Chem. Inf. Model.* 2024, 64, 3149–3160, https://doi.org/10.1021/acs.jcim.4c00115.

Relevant behavior found in source:

- input rows contain molecule SMILES, P450 sequence, label, and dataset type;
- the model combines Uni-Mol and ESM-2 (`esm2_t33_650M_UR50D`);
- `main/model.py` applies a two-class softmax and writes class-1 probability as `label_predict`;
- `main/metric.py` expects an ensemble `mean_prob` across 20 submodels.

For this benchmark, each OpenADMET molecule should be paired separately with the fixed human CYP2D6 and CYP3A4 sequences. The score of interest is the unrounded class-1 ensemble mean.

Blockers and cautions:

- the Git repository does not contain the 20 trained `.pt` files;
- the README points to a live Baidu Pan bundle containing checkpoints `1.pt`–`20.pt` plus `esm2-finetune.pt`, approximately 58.80 GB in total;
- several paths are machine-specific and require an adapter or configuration layer;
- no license file or license declaration was found in the inspected Git revision.

Do not mark the DeepP450 registry rows `ready` until the full asset set is retrieved, hashed, a representative inference smoke succeeds, and reuse terms are clarified.

## CYPMol

Purpose: one architecture with substrate, inhibitor, and bond-of-metabolism tasks conditioned on CYP sequence/pocket features.

Paper: Chang et al., *J. Chem. Inf. Model.* 2026, 66, 8980–8997, https://doi.org/10.1021/acs.jcim.6c01507.

Authoritative assets:

- Git source: https://github.com/CjmTH/CYPMol
- Hugging Face model repository: https://huggingface.co/CJM1111/CYPMOl
- inspected Hugging Face revision: `3810a93e04d1dbab4cde73792a1e3830b1f17e41`
- the Hugging Face repository is public, ungated, and tagged `license:apache-2.0`.

The inspected HF revision contains ten substrate checkpoints and ten inhibitor checkpoints (about 28.1 GB per task), plus one BoM checkpoint. The substrate and inhibitor scripts apply two-class softmax and write `label_predict[:, 1]`; this class-1 value is the required continuous score. The example input explicitly contains CYP2D6 and CYP3A4 sequence and pocket/SRS fields.

The public web UI exposes single-SMILES substrate and inhibitor endpoints at `https://tianlab-tsinghua.cn/cypmol/api/predict_substrate/` and `/api/predict_inhibitor/`, returning streamed JSON with `protein_name`, `score`, and `active`. Both inference endpoints returned a backend connection-refused error during the 2026-08-20 audit, so they are not a dependable batch route at present.

Required adapter work:

1. create two endpoint-expanded input tables from `data/model_input.csv`, one row per molecule × CYP;
2. populate the exact published CYP sequence, pocket-site, and SRS fields from the examples/source;
3. replace hard-coded `/data/cjm/...` paths with command-line paths without changing model math;
4. preserve molecule IDs in outputs rather than joining by row order;
5. ensemble the published checkpoints exactly as specified by the authors, retaining unrounded class-1 probabilities;
6. emit four columns: CYP2D6/CYP3A4 substrate and inhibitor probabilities.

This should use a dedicated environment because the repository inherits DeepP450's older Uni-Core/Uni-Mol/ESM stack and hard-codes CUDA behavior. Do not mutate the stable `cheminf` environment.

Completed benchmark runtime: Python 3.10.20, PyTorch 2.11.0+cu128,
fair-esm 2.0.0, RDKit 2026.03.2, and source-pinned Uni-Core 0.0.1 / Uni-Mol
v0.1 on an RTX 5070 Ti. Each checkpoint contains the complete Uni-Mol, ESM-2,
fusion and classifier state, so external base/fine-tuned initialization files
were not needed. Strict key/shape parity was required before loading.

The adapter applies the authors' AddHs, random-seed-42 embedding, MMFF,
heavy-atom and coordinate-centering preparation. `OCNT-0453746` and
`OCNT-0453782` failed the published conformer recipe and remain missing. The
public web API was retested on 2026-08-21 and still returned HTTP 500 because
its backend connection was refused.

## DeepMetab

Purpose: multitask prediction of CYP substrate status, sites of metabolism, and products for nine CYP enzymes.

Paper: Zhou et al., *Chemical Science* 2025, 16, 18884–18902, https://doi.org/10.1039/d5sc04631a; open full text at https://pmc.ncbi.nlm.nih.gov/articles/PMC12439208/.

Relevant behavior found in source:

- the README specifies Python 3.8, RDKit, and ChemProp;
- `Substrates/predict.py` defines nine target columns including CYP2D6 and CYP3A4;
- five multitask checkpoints are listed in `Substrates/load_model.py`;
- `ensemble_pred()` returns the mean and variance of continuous model predictions;
- `predict_main.py::substrate_result()` immediately thresholds the mean at 0.5 and returns only binary labels.

For ranking metrics, the adapter must call or preserve `ensemble_pred()` output before `get_pred_label()`. The required score is the unthresholded five-model mean probability. The variance may be retained as auxiliary provenance but is not one of the requested comparison scores.

The repository uses Git LFS for model files. All five public checkpoint objects were confirmed downloadable and are 3,184,162 bytes each; the inspected `model_1.pt` pointer identifies SHA-256 `d81e206cacabc076c8e83b32083053b56bde3048d8d18a05e4969dc6374473b5`. Pull and hash all five real checkpoints before inference. The repository includes an MIT license. The code uses the ChemProp 1.x API; `chemprop==1.5.2` is a sensible compatibility pin but was not specified by the authors.

Checkpoint compatibility audit: the public state dictionaries call the output
layers `readout.1` and `readout.4`, whereas ChemProp 1.5.2 constructs equivalent
layers as `ffn.1` and `ffn.4`. ChemProp's standard `load_checkpoint()` only warns
and skips the unmatched output layers, leaving a random prediction head. The
benchmark adapter performs the explicit name mapping and then requires exact
key and tensor-shape parity before loading; using unmodified
`load_checkpoint()` would not reproduce the trained substrate classifier.

## CypReact from CyProduct

Purpose distinction:

- CypReact predicts whether a molecule is a reactant/substrate for a selected CYP;
- CyProduct uses CypReact as a filter, then predicts bonds/sites and generated metabolites;
- CyProduct's documented `Score = (prob-threshold)/(1-threshold)` is a site-of-metabolism score from the CypBoM stage, not a molecule-level CypReact probability.

The requested benchmark should therefore use the molecule-level CypReact output, not the maximum or sum of CyProduct site scores.

The published APIs do not expose a continuous CypReact probability:

- for most CYPs, CypReact calls Weka `classifyInstance()` and returns a boolean;
- for CYP2D6, it averages five classifier hard labels and then thresholds the vote at 0.5;
- CyProduct's API receives only that boolean filter decision.

A faithful public-output benchmark should therefore use:

- CYP3A4: the public R/N call encoded as 1/0;
- CYP2D6: the mean of the five hard classifier votes, taking values 0, 0.2, ..., 1.

These are coarse ranking scores, not calibrated probabilities. Weka `distributionForInstance(oneSample)[1]` could be evaluated later as a separate exploratory adapter, but the released code never calls it for CypReact; its class-index and calibration semantics must first be checked against the serialized model header and original hard decisions. It must not silently replace the public score. The normalized public-output columns are `CYP2D6_cypreact_public_score` and `CYP3A4_cypreact_public_score`.

The clearest runnable public artifact is the historical CypReact v1.2 bundle at https://github.com/Le0nT1/CypReact_old/releases/tag/1.2. It contains the JAR, all nine model directories, and normalization files. The current source repository does not publish a runnable release JAR.

The exact v1.2 release ZIP used by this benchmark has SHA-256
`2778b65dc786125a9f6af489d599f6de6bbe81952f3d484285d0dc73370da48f`.
The adapter preserves molecule IDs through SDF titles and uses the unmodified
public CLI output: `N=0`, `R=1` for both endpoints. Although CYP2D6 internally
ensembles five classifiers, the historical CLI exports only the thresholded
final R/N decision, not the vote fraction.

CypReact contains an MIT license. No license file was found in the inspected CyProduct repository, so use and redistribution of CyProduct-specific code/assets need separate clarification. The historical v1.2 bundle was executed unchanged under Java 8 through an ID-preserving SDF/CSV adapter. It returned both endpoint calls for 4,821 of 4,822 molecules; `OCNT-2328840` was omitted by the public tool and is retained in the failure ledger rather than imputed.

## Provenance required before a registry row becomes ready

For every comparison record:

- exact source URL and immutable Git/Hugging Face revision;
- upstream license and any asset-specific license;
- model/checkpoint SHA-256 values;
- environment lock or container digest;
- command and adapter revision;
- raw input/output hashes;
- number of input rows, successful predictions, parse failures, and model failures;
- confirmation that class-1 means substrate/inhibitor/reactant positive;
- confirmation that no OpenADMET TDI labels were used in inference or tuning.

Only after these checks should `status` change from `planned` to `ready` in `config/model_registry.csv`.
