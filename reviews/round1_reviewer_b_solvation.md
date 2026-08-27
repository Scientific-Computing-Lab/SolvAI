# Round 1 — Reviewer B: molecular simulation and solvation

## Critical issues

- PIMD8 must not appear as a retained teacher. The final selected feature block
  contains zero PIMD-trained coordinates. Figure 1 and the text now isolate
  PIMD8 as an accuracy reference; the PIMD2/PIMD8 supervision experiments are
  labelled non-retained ablations.
- Endpoint experimental supervision must be visible. Figure 1 now shows the
  experimental hydration labels entering endpoint training separately from
  response-surrogate training.

## Important issues

- Distinguish short PIMD2 response fingerprints from converged free energies.
- State that direct response integration failed because component errors of
  1.27–5.20 kcal mol−1 accumulate along the alchemical path.
- Avoid putting SMD, conformer records and molecular identities on one common
  count axis.

## Optional polish

- Use “nuclear delocalization” before the abbreviation NQE.
- Describe the published ARROW timing as total-work context because the
  hardware differs from the released inference benchmark.

## Revision made

The physical training/inference boundary, lambda-probe status and source units
were corrected in the manuscript, Figure 1, Extended Data Figures 2 and 5, and
Supplementary Methods.
