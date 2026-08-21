# OpenADMET CYP Challenge

Code, protocols, and scientific notes for the [OpenADMET cytochrome P450 inhibition challenge](https://huggingface.co/spaces/openadmet/cyp-challenge).  
My handle in the challenge is **dargason**.

This repository is under heavy development and frequent updating. (8/21/2026) 

## Project site

A prettier front-end to this project is at:  https://admet-cyp.aetherark.com/


## Potentially interesting artifacts:

I've generated models and datasets along the way and below are some that may be of particular use:

* [**Cofolding dataset**](https://huggingface.co/datasets/dargason/ADMET-CYP-cofolding): Multi-method cofolding (16,590) of CYP with PDB (149) ligands using along with PoseBusters, ProLIF and other structural features.  Cofolding methods:  Boltz2, Chai-1, ESMFold2, OpenFold3, and Protenix v1.



## Repository layout

- `intro/` - challenge notes, assay interpretation, canonical sequences, and the Bayesian claim-evidence preregistration.
- `scripts/` - scripts that build, QC, finalize, and verify the Hugging Face structure release.
- `experiments/` - logged experiments: model training, etc...
- `submissions/` - submission folders to the regression and classification tracks. Eacvh 
- `licenses/` - license text consumed by the release builder.


## License

This repository's code and OpenADMET-generated documentation are licensed under Apache License 2.0. Experimental PDB coordinates retain their upstream CC0 status; see the Hugging Face dataset's `LICENSES.md` for component details.
