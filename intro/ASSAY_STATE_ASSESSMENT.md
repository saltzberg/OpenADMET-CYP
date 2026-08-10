# Assay-specific CYP states for structural feature generation

Assessed 2026-07-30T20:48:34Z. This separates facts stated for the challenge from provisional evidence in OpenADMET's earlier Octant CYP protocol.

## Decision summary

Use the cofolded parent–CYP complexes as one input to a feature system, not as literal models of a single native CYP state. The challenge contains two materially different biochemical states:

1. **Direct-inhibition arm:** parent compound is preincubated with recombinant CYP **without NADPH**, then the assay probe is added and residual catalytic activity is measured. The relevant structural hypotheses are parent binding to the resting enzyme and parent–probe competition or co-occupancy.
2. **TDI arm:** parent compound is preincubated with recombinant CYP **with NADPH before probe addition**. The parent can be turned over to metabolites or intermediates that bind more strongly or inactivate protein/heme. Parent binding alone is insufficient; catalytic orientation and metabolic/reactive chemistry are relevant.

The present CYP3A4 + Fe2+ heme + one parent-ligand model is useful for pocket and orientation features, but it is not an exact representation of either assay arm.

## Confirmed by current challenge sources

Sources:

- `intro/ANNOUNCEMENT.md`, derived from the OpenADMET announcement published 2026-07-29.
- Current public announcement: <https://openadmet.ghost.io/announcing-openadmets-cyp-inhibition-blind-challenge/>
- Current public Space source at commit `6af515996f59a276872a409e117b344a8903006e`, especially `config.py`.

Confirmed assay facts:

- Biochemical recombinant-enzyme assays, not cells.
- CYP3A4 uses DBOMF; CYP1A2 and CYP2C9 use EOMCC. These are masked fluorogenic probes.
- CYP2D6 uses dextromethorphan parent depletion measured by Echo-MS/ZenoTOF without chromatography.
- Each DRC has matched preincubation arms. The direct arm lacks NADPH during preincubation; the TDI arm contains NADPH. Probe is added after preincubation.
- Direct regression targets are four probe-dependent apparent pIC50 values, not binding free energies or Ki values.
- TDI labels are based on a >2-fold IC50 shift and scored only for CYP3A4 and CYP2D6.
- CYP active sites can bind compounds in multiple orientations and, especially for CYP3A4, multiple copies.
- The full assay recipe, exact constructs, concentrations, incubation times, and final data schema are explicitly deferred to launch or a subsequent assay-development post.

## Strong but provisional context from the precursor Octant protocol

OpenADMET's earlier public CYP3A4 inhibition protocol is:

<https://github.com/OpenADMET/Octant_CYP_blog_post/blob/main/protocols/cyp_inhibition_assay.md>

It used:

- DLS Gentest CYP3A4 Supersomes, catalog 456202.
- The catalog product is CYP3A4 + P450 oxidoreductase + cytochrome b5 in a microsomal membrane preparation.
- 5 nM CYP3A4, 100 µM NADP+, regeneration system, 2 µM DBOMF, 100 mM potassium phosphate at pH 8.
- 30-minute active preincubation and 30-minute probe reaction.

This strongly suggests that a membrane/POR/b5-supported recombinant system may underlie the challenge platform, but the challenge announcement does not yet confirm that the same product, concentrations, buffer, or timing were used for each isoform. Do not freeze these as challenge constants until the launch protocol confirms them.

## Structural states worth representing

### A. Direct-inhibition state

Primary representation:

- Parent compound in its plausible assay-pH protomer/tautomer states.
- Isoform-specific catalytic domain ensemble.
- Resting ferric heme as the chemically relevant preincubation reference, including explicit axial Cys thiolate and a governed distal-water option.

Probe-aware representation:

- Generate a baseline structure/pose ensemble for DBOMF–CYP3A4, EOMCC–CYP1A2, EOMCC–CYP2C9, and dextromethorphan–CYP2D6.
- Compare each test-compound ensemble to its assay probe: spatial overlap, contact displacement, channel occupancy, pocket-volume competition, and possible simultaneous accommodation.
- For CYP3A4, retain a multiple-occupancy/cooperative-state branch rather than assuming simple one-site competition.

The measured IC50 is probe- and condition-dependent. Without probe concentration and Km, it cannot be interpreted as Ki or direct affinity.

### B. TDI/metabolism-competent state

Use the parent pose as a substrate-orientation hypothesis, then calculate:

- Candidate site-of-metabolism atom distances and approach angles to Fe and to a virtual distal oxo.
- Heme-coordinating N/S/O geometry that could produce reversible or quasi-irreversible type-II binding.
- Number and persistence of catalytically accessible oxidation sites across poses.
- Reactive-metabolite and mechanism-based-inactivation alerts.
- Potential protein/heme adduct targets and metabolite re-binding features.

Do not interpret a static Fe2+ cofold as a full TDI state. NADPH drives a catalytic ensemble spanning ferric substrate-bound, reduced/oxygenated intermediates, Compound I-like chemistry, products, and possible covalent or metabolic-intermediate complexes. Geometry can inform this branch, but oxidation chemistry must come from separate descriptors or simulations.

### C. Membrane/cofactor context

If the challenge confirms Supersomes, the experimental entity includes a microsomal membrane, POR, and possibly cytochrome b5—not an isolated soluble catalytic domain. The present crystal-core constructs remain reasonable for active-site features, but they omit:

- membrane orientation and access channels;
- CYP–POR/b5 coupling and electron-transfer context;
- full N-terminal membrane anchors;
- membrane partitioning of hydrophobic compounds.

Do not add full membrane complexes to the initial feature pipeline unless validation shows value. Instead retain 2D lipophilicity/charge descriptors and add membrane-access/channel features as a separate, testable block.

## Recommended feature blocks

### Shared structural QC and uncertainty

- Ligand topology/stereochemistry pass flags.
- Per-engine and per-seed ensemble dispersion; no silent averaging of failed poses.
- Protein/heme geometry and pocket-confidence flags.
- No-winner flag when no chemically valid pose survives.

### Direct-inhibition features

- Pocket contacts and interaction fingerprints.
- Ligand burial, solvent exposure, pocket volume, shape complementarity, and channel occupancy.
- Fe distance and heme-coordination geometry.
- Assay-probe overlap/displacement and possible co-occupancy.
- Isoform-relative feature differences for the same compound.

### TDI features

- Everything above, plus site-of-metabolism/virtual-oxo geometry.
- Accessible C–H and heteroatom oxidation-site counts.
- Heme-ligation and metabolic-intermediate-complex potential.
- Reactive-metabolite/covalent-warhead alerts and likely metabolite states.
- Pose multiplicity: tight nonproductive binding versus catalytically oriented binding.

### Assay-context features

- Probe identity: DBOMF, EOMCC, or dextromethorphan.
- Readout modality: fluorescence versus parent-depletion Echo-MS.
- Preincubation arm and isoform.
- Fluorescence-interference/autofluorescence flags for CYP3A4/2C9/1A2; ionization/parent-detection flags for CYP2D6 where available.

These assay variables should be explicit model inputs rather than left confounded with isoform.

## What to request/check at data release

Do not lock the final structural protocol until these fields are known:

1. Exact recombinant enzyme products/lots and sequence constructs for each isoform.
2. POR and cytochrome b5 inclusion and ratios.
3. CYP, probe, NADPH/NADP+, regeneration-system, and compound concentrations.
4. Buffer composition, pH, temperature, DMSO, preincubation time, and reaction time.
5. Exact order of additions in both arms.
6. Probe Km or assay substrate concentration relative to Km.
7. Whether the direct and TDI arms share all conditions except NADPH preincubation.
8. Curve-censoring and credible-interval fields.
9. Raw fluorescence or Echo-MS interference/QC flags.
10. Whether CYP2D6 reports dextromethorphan disappearance only or also metabolite formation.

## Practical conclusion

The highest-value structural design is not “one native CYP state.” It is a small assay-aware ensemble:

- **resting/probe-competition branch** for direct pIC50;
- **substrate/reactive-orientation branch** for TDI;
- **isoform- and probe-specific context** for all endpoints.

Our current cofolding pipeline supplies a useful parent-bound geometry block. It should be augmented with ligand microstates, a ferric/resting heme reference, probe-aware comparisons, virtual-oxo/site-of-metabolism geometry, and explicit assay metadata. The resulting features must then be tested as an orthogonal increment over the clean 2D baseline after the activity drop.
