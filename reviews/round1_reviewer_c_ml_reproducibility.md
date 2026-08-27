# Round 1 — Reviewer C: machine learning and reproducibility

## Critical issues

1. Teacher leakage must be audited at connectivity level, not only by raw
   SMILES.
2. Feature-block choice must be nested inside outer folds.
3. The deployment refit on all 85 labels must not be confused with OOF scoring.
4. The schematic must depict a model stack, not a fictitious monolithic
   teacher–student network.

## Important issues

- Freeze every manuscript number from molecule-level outputs.
- Hash all executable artifacts and test from raw SMILES.
- Preserve seeds, split assignments, dependency versions and non-improving
  ablations.

## Optional polish

- CI should fail on metric drift, artifact drift, leakage or LaTeX errors.

## Revision made

The release now has one metric module with numerical assertions, independent
identity and runtime-schema audits, a nested-selection record, frozen repeat
predictions, artifact SHA-256 values and end-to-end smoke tests. The paper states
that deployment is a post-evaluation refit and never supplies its predictions as
OOF evidence.
