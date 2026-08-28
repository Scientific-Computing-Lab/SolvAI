# Confirmatory round 1 — machine learning and reproducibility reviewer

## Critical issues

- The old connectivity-only audit is insufficient after standardized equivalents were
  found. Make the corrected refits and split-preservation checks part of `make verify`.
- The matched baseline must differ only by the 15 prior columns.

## Important issues

- Publish shuffled-prior predictions, all repeat predictions and global-separation
  training counts.
- Make the all-data deployment refit distinct from evaluation artifacts.

## Optional polish

- Drive smoke expectations from the model manifest rather than copied constants.
