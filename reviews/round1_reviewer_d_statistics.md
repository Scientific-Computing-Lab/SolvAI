# Round 1 — Reviewer D: statistics and skeptical generalist

## Critical issues

1. A single 0.197 estimate from 85 molecules cannot establish robust sub-0.20
   population performance.
2. The PIMD8 and SolvAI point estimates should not be called equivalent without
   an equivalence design.
3. Hard chemical splits must not be hidden behind random OOF.

## Important issues

- Report all five repeats, mean and sample standard deviation.
- Use paired molecule resampling for method differences.
- Present wide bootstrap intervals as a consequence of sample size and
  heterogeneity, not as a model-selection device.

## Optional polish

- Avoid treating 0.20 as a natural discontinuity.

## Revision made

The manuscript reports 0.197 together with 0.204 ± 0.005 across repeats, uses
“reaches PIMD8-level accuracy,” gives paired bootstrap intervals, reports family
and scaffold MAEs near 0.24, and labels the 0.20 line as a visual reference.
