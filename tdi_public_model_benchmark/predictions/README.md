# Prediction artifacts

Place immutable raw or minimally normalized public-model outputs in `predictions/raw/`. CSVs are ignored by Git; preserve their hashes and upstream provenance in the model registry or a companion manifest.

Each registry row identifies one score column for one CYP endpoint. At minimum, its CSV must contain:

- a stable molecule ID column, normally `Molecule_Name`;
- the numeric score column named by `score_column`.

Blank scores are allowed and count against coverage. Duplicate molecule IDs, infinities, nonnumeric score text, and unsupported endpoints fail evaluation. Scores are never imputed.

Do not put invented, randomly generated, or retrained surrogate scores here under a public model's name. Tests use clearly synthetic temporary predictions only.
