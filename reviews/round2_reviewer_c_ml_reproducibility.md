# Round 2 — Reviewer C: machine learning and reproducibility

## Critical issues

None remaining. Frozen OOF predictions, nested selection, source exclusion,
artifact hashes and inference-only dependencies are separated coherently.

## Important issues

CombiSolv-QM source rows must remain outside Git because redistribution terms
are not explicit. The checkpoint, hash, DOI and opt-in download path preserve
the maximum reproducibility compatible with that constraint.

## Optional polish

Keep Git LFS enabled for the model binaries and processed CC BY tables.
