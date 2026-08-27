# Round 1 — Reviewer B: molecular simulation and solvation

## Critical issues

1. Distinguish SMD/COSMO-RS response labels, explicit-solvent alchemical labels
   and PIMD. They are not interchangeable physical calculations.
2. Do not imply that a three-point, 5-ps PIMD2 response is a converged free
   energy.
3. The ARROW/PIMD8 comparison must use the reconstructed same-set value and
   must not reuse the ARROW-NN water-dimer error as a hydration result.

## Important issues

- State that response quantities at inference are surrogate predictions.
- Describe $\mathrm{d}H/\mathrm{d}\lambda$ component errors in their exported
  energy units and explain why integration failed.
- Separate total PIMD work from ideal parallel wall time in runtime discussion.

## Optional polish

- Use “comparable to” or “PIMD8-level” rather than statistical equivalence.

## Revision made

The Methods now enumerate all 15 priors and their sources, identify PIMD8 as the
comparator rather than a selected label, define the short PIMD2 calculations as
response fingerprints, and state the published PIMD work protocol separately
from SolvAI latency.
