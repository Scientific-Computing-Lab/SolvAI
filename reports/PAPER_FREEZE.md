# SolvAI paper freeze

This file freezes the quantitative state used by the manuscript. It is generated
only by `solv_ai.paper_metrics`; hand-edited numerical values are not authoritative.

## Scientific decision

SolvAI is a SMILES-only inference stack. The strict fixed five-fold OOF result is
**0.19705 kcal/mol**, and nested
feature-block selection is **0.19931**.
The fixed model averages **0.20375 ±
0.00493** across five independent split repeats.
Accordingly, the single-partition threshold crossing is valid but robust sub-0.20
generalization is not claimed.

## Frozen values

| Quantity | Value |
|---|---:|
| Molecules | 85 |
| Classical ARROW MAE | 0.78465 kcal/mol |
| ARROW/PIMD8 MAE | 0.20484 kcal/mol |
| Previous structure-only MAE | 0.23861 kcal/mol |
| Narrow-response MAE | 0.21362 kcal/mol |
| + MolSolv SMD(water) MAE | 0.20161 kcal/mol |
| + ConfSolv response MAE | 0.19705 kcal/mol |
| Nested-selection MAE | 0.19931 kcal/mol |
| Five-repeat fixed mean ± SD | 0.20375 ± 0.00493 kcal/mol |
| Five-repeat nested mean ± SD | 0.20860 ± 0.00275 kcal/mol |
| Family-held-out MAE | 0.23957 kcal/mol |
| Scaffold-held-out MAE | 0.24128 kcal/mol |
| MolSolv source / retained | 1,729,545 / 350,391 |
| ConfSolv source / retained | 5,392,567 / 39,878 |
| ARROW-85 connectivities also in FreeSolv | 80/85 |
| Simulation at inference | No |

## Integrity

- Every OOF method has exactly 85 unique molecule predictions.
- Every repeated method has five complete 85-molecule partitions.
- Stored absolute errors equal values recomputed from predictions and targets.
- External count and model-schema assertions pass before files are written.
- Bootstrap interval: molecule-level resampling, 100,000 replicates,
  seed 20260827.
