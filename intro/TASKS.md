# Tasks

## Build an agentic crystal-structure assessment workflow

Develop an agentic workflow so the full CYP crystal-structure set can be assessed without manual inspection of every entry.

The workflow should:

- retrieve and preserve the deposited coordinates, structure factors, and electron-density maps with explicit provenance;
- validate the experimental model and identify the relevant CYP chains, heme, ligands, waters, alternate conformers, missing residues/atoms, mutations, engineered termini, and crystallographic contacts;
- assess atom-level agreement with the experimental electron density, with particular attention to bound ligands;
- test whether each ligand's identity, occupancy, conformation, stereochemistry, alternate states, and modeled atom count are supported by the density;
- prefer the simplest ligand interpretation consistent with the density and chemistry, avoiding unsupported atoms, overfitting, or unjustified alternate conformers;
- distinguish deposited-model evidence from any automated reinterpretation or rebuilt model;
- flag ambiguous, weak-density, partially occupied, covalent, heme-coordinating, multi-ligand, and symmetry-related cases for human review rather than forcing a pass/fail decision;
- produce auditable per-entry and per-ligand evidence, metrics, visualizations, confidence levels, and a prioritized review queue;
- preserve a no-winner/no-reliable-ligand option when the experimental evidence does not support a clear interpretation.

The intended result is a reproducible triage system, not an automated claim that every deposited ligand pose is correct.
