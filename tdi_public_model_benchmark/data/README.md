# Data contract

Generated files in this directory are deliberately not committed.

`prepare-labels` reads the canonical OpenADMET `cyp-challenge-TRAIN_TDI.csv` and writes:

- `tdi_labels.csv`: one row per labelled molecule × endpoint, including `is_tdi` and `scaffold_id`;
- `model_input.csv`: one row per unique labelled molecule for public-model inference;
- `label_manifest.json`: source/output hashes, counts, class balance, scaffold counts, and parse failures.

The source assay file remains in the parent repository under `data/cyp-challenge-train-test/`. The benchmark does not copy or edit it.
