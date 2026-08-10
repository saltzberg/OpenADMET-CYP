# Experimental PDB structures for the CYP panel

Retrieved from UniProtKB PDB cross-references and downloaded from RCSB PDB on 2026-07-30.

## Coverage

| Target | UniProtKB | Mapped PDB entries | Legacy `.pdb` files | Authoritative `.cif` fallback files |
|---|---|---:|---:|---:|
| CYP3A4 | P08684 | 119 | 101 | 18 |
| CYP2C9 | P11712 | 15 | 15 | 0 |
| CYP2D6 | P10635 | 14 | 14 | 0 |
| CYP1A2 | P05177 | 1 | 1 | 0 |
| **Total** |  | **149** | **131** | **18** |

All UniProt-mapped experimental entries were retrieved. RCSB did not provide legacy PDB-format files for 18 recent CYP3A4 entries and returned HTTP 404 for those `.pdb` URLs. Their authoritative RCSB mmCIF coordinate files were retained instead of performing a silent format conversion.

## Layout

- `CYP3A4/`
- `CYP2C9/`
- `CYP2D6/`
- `CYP1A2/`
- `manifest.csv` — one row per mapped PDB entry, including UniProt accession, method, resolution, mapped chains/residue range, coordinate format, source URL, relative path, byte size, SHA-256 checksum, and format note.
- `failures.json` — retrieval-failure ledger; currently empty.

## Provenance and interpretation

- Mapping authority: reviewed UniProtKB records P08684, P11712, P10635, and P05177.
- Coordinate authority: RCSB PDB deposited coordinate files.
- Files are raw deposited structures, not biological-assembly expansions or prepared modeling receptors.
- No chains, waters, ligands, heme groups, alternate conformers, mutations, tags, or coordinates were removed or changed.
- Experimental constructs commonly omit or replace native N-terminal membrane-anchor residues. Use the `uniprot_chains` field in `manifest.csv` as the initial construct-coverage reference, then inspect each deposition before modeling.
- A PDB cross-reference does not by itself establish a suitable ligand-bound template, inhibitor mechanism, oxidation state, complete active site, or preferred biological assembly.
