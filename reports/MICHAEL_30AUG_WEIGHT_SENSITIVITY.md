# ARROW endpoint-weight sensitivity

This one-factor sensitivity was frozen in
`release/MICHAEL_30AUG_SENSITIVITY_FREEZE.md` and committed as `9851ed2` before
fitting. It changes only the sample weight of available ARROW outer-training rows
from 3 to 1. The 1,280 external labels, fixed five folds, feature schema, teachers,
endpoint architecture and seeds are unchanged.

| Method | N | MAE | RMSE | Median absolute error |
|---|---:|---:|---:|---:|
| Matched structure-only | 85 | 0.30977 | 0.59585 | 0.14972 |
| Full SolvAI | 85 | 0.20642 | 0.32627 | 0.11653 |

The paired difference, defined as SolvAI minus structure-only absolute error, is
**-0.10335 kcal mol-1** (95% 100,000-resample molecule-bootstrap interval,
**-0.21883 to -0.01995**). SolvAI improves 52 of 85 molecules (61.2%). Under the
prospectively frozen rule, this is a positive result.

The response-prior advantage therefore does not depend on upweighting available
ARROW outer-training labels. The separate zero-ARROW-label analysis is stronger with
respect to reference-set adaptation (0.385 versus 0.257 kcal mol-1), although it is
not an independent external benchmark.

Machine-readable results are in `results/michael_30aug_sensitivity/`.
