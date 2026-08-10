# Publication boundaries

This project uses GitHub and Hugging Face for different artifacts. The split is intentional: GitHub remains reviewable as source code, while Hugging Face provides the versioned scientific dataset.

## Commit to GitHub

Commit:

- reusable Python and shell source;
- tests and CI configuration;
- scientific protocols and methodological documentation;
- the Bayesian framework preregistration;
- canonical target FASTA files and small source manifests;
- environment specifications;
- licenses, citations, and links to immutable dataset versions.

Before committing a script from `runs/` or `analysis/`, move or copy it into a source-owned directory and remove machine-specific assumptions. In particular, replace hard-coded home directories, container paths, IP addresses, SSH-key paths, and output roots with command-line arguments or documented configuration.

## Publish to Hugging Face

Publish as a dataset:

- aligned experimental and predicted coordinate files;
- `data/structures.parquet` and other row-level scientific tables;
- coordinate and file checksums;
- alignment transforms, anchor definitions, coverage tables, field groups, and exclusion ledgers;
- exact release-building and verification scripts copied into the dataset snapshot;
- dataset card, citation, notices, and component licensing.

The current dataset is:

https://huggingface.co/datasets/dargason/ADMET-CYP-cofolding

The local release root is `release/huggingface_cyp_cofold_v1/`. Its `.cache/` directory is Hugging Face client state and is never part of either publication.

## Keep local or in archival object storage

Do not publish as ordinary Git history:

- raw `runs/` trees and generated engine inputs;
- model checkpoints or third-party model weights;
- full MSAs and duplicated per-run MSA copies;
- logs, PID files, heartbeats, temporary status files, and upload caches;
- intermediate coordinate trees superseded by the governed Hugging Face release;
- editor history and Python bytecode;
- credentials, tokens, SSH keys, private endpoints, or personal infrastructure details.

Raw run artifacts should remain immutable locally or in archival object storage. Their durable public representation is a compact protocol and provenance record, not a dump of the working directory.

## Release gate

Before a GitHub push:

1. Run `python scripts/check_repository_payload.py --mode staged`.
2. Inspect `git diff --cached --stat` and `git diff --cached`.
3. Confirm that no generated coordinates, MSAs, Parquet files, checkpoints, logs, or caches are staged.
4. Run a secret scanner if one is available.

Before a Hugging Face update:

1. Build into a new versioned release root.
2. Run the full release verifier without `--skip-coordinate-parse`.
3. Confirm expected versus observed rows, methods, PDB IDs, exclusions, and coordinate counts.
4. Preserve the previous dataset revision; do not overwrite scientific meaning silently.
5. Tag or document the exact dataset commit used by downstream analyses.
