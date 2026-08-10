# OpenADMET CYP

Code, protocols, and scientific notes for the OpenADMET cytochrome P450 inhibition challenge.

This repository deliberately separates source code from large scientific artifacts:

- GitHub contains code, documentation, small metadata, and canonical target sequences.
- Hugging Face contains the cofolded and experimental structures, the row-level Parquet table, QC results, checksums, and release provenance.
- Raw run directories, MSAs, caches, logs, and intermediate analyses remain outside version control.

## Project site

The project summary, methodology, and data index are published at:

https://aetherark.com/sites/admet-cyp/

The Bayesian preregistration is under **Methodology**. The cofolding release is under **Data**.

## Cofolding dataset

The published structure release is available at:

https://huggingface.co/datasets/dargason/ADMET-CYP-cofolding

Version 1 contains 16,739 aligned structures: 149 experimental PDB structures and 16,590 predictions from Boltz2, Chai-1, ESMFold2, OpenFold3, and Protenix v1. The dataset card documents coverage, alignment, ligand handling, ProLIF and PoseBusters QC, limitations, and component licensing.

To download and verify it:

```bash
conda create -n admet-cyp-release python=3.12 -y
conda activate admet-cyp-release
python -m pip install -r requirements-release.txt

hf download dargason/ADMET-CYP-cofolding \
  --repo-type dataset \
  --local-dir data/ADMET-CYP-cofolding

python scripts/verify_hf_release.py data/ADMET-CYP-cofolding
```

The verifier checks the table schema, row and method counts, coordinate paths, coordinate hashes, mmCIF parsing, and the release checksum manifest.

## Repository layout

- `intro/` — challenge notes, assay interpretation, canonical sequences, and the Bayesian claim-evidence preregistration.
- `scripts/` — scripts that build, QC, finalize, and verify the Hugging Face structure release.
- `licenses/` — license text consumed by the release builder.
- `docs/PUBLICATION_BOUNDARIES.md` — what belongs on GitHub, Hugging Face, or local storage.

The current `analysis/`, `runs/`, `inputs/`, `work/`, `release/`, and editor directories are intentionally excluded from Git. They mix large generated artifacts with machine-specific paths. Reusable campaign and analysis code should be extracted into clean source directories before it is committed.

## Scientific scope

The cofolded structures are structural hypotheses and feature sources. They are not direct measurements of CYP inhibition, universal binding affinities, or complete representations of time-dependent inhibition. Failed or ligand-omitting cofolds are not treated as apo fallbacks.

## License

Repository code and OpenADMET-generated documentation are licensed under Apache License 2.0. Experimental PDB coordinates retain their upstream CC0 status; see the Hugging Face dataset's `LICENSES.md` for component details.
