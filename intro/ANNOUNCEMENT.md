# OpenADMET CYP Inhibition Blind Challenge

## Source and status

- Announcement: [Announcing OpenADMET’s CYP inhibition blind challenge](https://openadmet.ghost.io/announcing-openadmets-cyp-inhibition-blind-challenge/)
- Published: 2026-07-29
- Challenge platform: [OpenADMET CYP Challenge on Hugging Face](https://huggingface.co/spaces/openadmet/cyp-challenge)
- Support: [OpenADMET Discord challenge channel](https://discord.com/channels/1412827471488745545/1480419832787763412)
- Status described in the announcement: challenge announced; launch and training-data release planned for 2026-08-17. No challenge dataset has been added to this project yet.

## Scope

The challenge targets inhibition of four human cytochrome P450 isoforms central to small-molecule drug metabolism:

1. CYP3A4
2. CYP2C9
3. CYP2D6
4. CYP1A2

OpenADMET says CYP2C19, CYP2C8, CYP2B6, and other isoforms may be addressed later; they are not part of this challenge.

## Assays

These are recombinant-enzyme biochemical assays, not cell-based assays.

- CYP3A4, CYP2C9, and CYP1A2: fluorescence-based assays adapted from ThermoFisher Vivid kits and miniaturized to 1536-well plates.
  - CYP3A4 probe: DBOMF.
  - CYP1A2 and CYP2C9 probe: EOMCC.
- CYP2D6: label-free parent-depletion assay because the tested fluorogenic probes had inadequate signal-to-noise and dynamic range.
  - Substrate: dextromethorphan.
  - Readout: acoustic ejection mass spectrometry on a SCIEX Echo-MS / ZenoTOF 7600, without chromatography.

Each dose-response experiment has two preincubation arms:

- Direct inhibition: preincubation without NADPH. This measures inhibition by the parent compound without metabolism-dependent activation.
- Time-dependent inhibition (TDI): preincubation with NADPH. This permits turnover and captures direct inhibition plus additional inhibition caused by metabolites or mechanism-based inactivation.

The TDI signal is the leftward IC50 shift in the +NADPH arm relative to the -NADPH arm.

## Dataset construction announced

1. OpenADMET screened Enamine DDS10 and FDA-approved-compound (FDAA) libraries at one concentration against each CYP.
2. The primary screen used the +NADPH/TDI condition to maximize positive detection.
3. Processing included variance stabilization, spatial-artifact correction, a per-plate linear model against negative controls, and Benjamini-Hochberg false-discovery-rate correction. Negative log-fold change denotes inhibition.
4. Selected hits received 12-point dose-response curves, with Bayesian curve fitting and both direct and TDI arms. The target was approximately 1,500 DRCs per CYP.
5. For the blinded test set, the top 25 hits from each of CYP1A2, CYP2C9, and CYP3A4 were expanded with the 10 nearest Enamine US in-stock chemisimilars by Tanimoto similarity.
6. This produced 750 test compounds. All 750 were assayed densely against all four CYPs, including CYP2D6. The announced training DRC matrix is sparse across CYPs; the test matrix is dense.

CYP2D6 did not seed hit-expansion series because its Echo-MS assay was still being developed at that point.

The announced training pack will include the full-library primary screen in the TDI condition and about 1,500 DRCs per isoform with both preincubation arms.

## Prediction tracks

### 1. Direct inhibition regression

Predict direct-inhibition pIC50 for CYP3A4, CYP2D6, CYP2C9, and CYP1A2 for every test compound: four regression targets per compound.

- Primary metric: modified macro-averaged Soft-Threshold Relative Absolute Error (MA-ST-RAE), averaged over the four endpoints.
- Low-activity compounds are downweighted.
- The metric uses fitted pIC50 credible intervals: predictions inside the interval incur zero error; otherwise error is measured to the nearest interval bound.
- Additional metrics will also be reported.

### 2. Time-dependent inhibition classification

Predict whether the IC50 shift exceeds two-fold. Predictions are requested for every test compound, but only confidently assigned labels contribute to scoring.

- Evaluated isoforms: CYP3A4 and CYP2D6 only.
- Metric: Matthews correlation coefficient (MCC).
- For direct-inhibition pIC50 > 4, the label is determined by whether the TDI shift exceeds two-fold.
- If direct-inhibition pIC50 is below 4 and TDI-arm pIC50 exceeds 4.301, the compound is an inferred positive because the shift must exceed two-fold (`log10(2) = 0.301`).
- If both arms are below pIC50 4, the compound is an assigned negative because the shift cannot be resolved reliably at that activity level.
- Positive class: measured positives plus inferred positives.
- Negative class: measured negatives plus assigned negatives.

## Leaderboard and challenge structure

- One challenge stage; unlike the PXR challenge, no half-test-set data release midway through the competition.
- Half of the test set contributes to the live leaderboard.
- The live/blinded split is by chemisimilar series so compounds descending from one parent remain together.
- At the intermediate deadline, full-test performance is revealed once.
- Further implementation details are deferred to launch day.

## Rules and awards announced

- Standard leaderboard awards plus an award for innovative machine-learning approaches. Innovation may be recognized independently of leaderboard rank.
- Only one leaderboard submission per cooperating team/lab.
- Use of proprietary data must be disclosed.
- Open/reproducible code will be tracked through a submission checkbox; a public repository with substantive code is expected for an open-code claim.
- OpenADMET plans shorter challenge-result blog posts rather than committing to a preprint for every new challenge. Earlier ExpansionRx and PXR preprint commitments remain.

## Timeline

Times are one minute before midnight UTC unless otherwise specified.

| Event | Date |
|---|---:|
| Challenge announced | 2026-07-29 |
| Challenge launch | 2026-08-17 |
| Intermediate-leaderboard submission deadline | 2026-09-24 |
| Intermediate leaderboard release | 2026-09-25 |
| Final submission deadline / challenge close | 2026-11-03 |
| Webinars, blog, and wrap-up | After 2026-11-03 |

## Notes for modeling

- Direct inhibition and TDI are distinct targets generated from matched -NADPH and +NADPH arms.
- CYP2D6 uses a different readout and probe system from the other three isoforms, so assay modality is confounded with isoform.
- Test-series construction is based on potent hits from CYP1A2, CYP2C9, and CYP3A4, not CYP2D6.
- The public primary screen is announced as +NADPH/TDI-condition data, whereas the regression target is the -NADPH direct-inhibition pIC50.
- Exact data schema, final metric implementation, submission format, and detailed assay protocol were not specified in this announcement and are expected at or after launch.
