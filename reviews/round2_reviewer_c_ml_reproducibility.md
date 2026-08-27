# Round 2 — Reviewer C: machine learning and reproducibility

## Critical issues

None.

## Important issues

None unresolved. The evaluation artifact, deployment refit, teacher labels and
endpoint labels are now explicitly separated. Exact split assignments,
molecule predictions, experiment ledger, teacher manifests and artifact hashes
are distributed in machine-readable form.

## Optional polish

The current subprocess-based D-MPNN inference could be engineered for lower
single-query latency without changing weights, but this is not required for the
scientific release and the measured latency is reported honestly.
