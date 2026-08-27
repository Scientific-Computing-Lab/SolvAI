# Round 1 — Reviewer C: machine learning and reproducibility

## Critical issues

- Clarify that the public deployment refit uses all 85 labels only after the
  evaluation is frozen. It cannot be the source of any reported held-out
  metric.
- Make the two training stages explicit: response-surrogate supervision and
  experimental endpoint supervision. A classic one-teacher/one-student diagram
  would misrepresent the implementation.

## Important issues

- Publish exact feature counts, endpoint hyperparameters, D-MPNN settings,
  teacher-source exclusions, repeat seeds and fold assignments.
- Do not hide that 1,280 benchmark-disjoint experimental labels train the
  endpoint.
- Preserve candidate non-improving experiments in a machine-readable ledger.

## Optional polish

- Keep “model stack” or “inference system” available as implementation terms;
  do not call SolvAI one monolithic neural network.

## Revision made

Supplementary Methods now specify every teacher and endpoint stage. Four
Supplementary Data packages contain the experiment ledger, predictions, fold
assignments and teacher manifests. The main architecture shows both kinds of
training supervision and a molecule-only deployment path.
