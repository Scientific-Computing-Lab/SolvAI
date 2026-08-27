# Structure-only lambda-response distillation

Updated: 2026-08-27 03:54 UTC

## Question

Can a structure-only student learn the short-PIMD2 alchemical response at
lambda = 0.1, 0.5, and 0.9, then use that learned response to improve hydration
free-energy prediction without any simulation at inference?

All PIMD2 observables are training-only privileged labels. In every outer fold,
the held-out molecule's response and classical/PIMD hierarchy labels are absent
from training. Every evaluated prediction is generated from molecular structure
and deterministic structure-derived features only.

## Teacher coverage

| Lambda | Attempted | Successful |
|---:|---:|---:|
| 0.1 | 76 | 72 |
| 0.5 | 76 | 74 |
| 0.9 | 76 | 73 |

There are 72 molecules with complete three-point curves. The fixed protocol used
PIMD2, one lambda state per run, and 5 ps per state. These short runs are not
treated as converged free energies.

## Locked three-seed random-OOF result

| Structure-only method | MAE (kcal/mol) |
|---|---:|
| A: matched SMD + ConfSolv response baseline | **0.19592** |
| B2: A + predicted PIMD2 lambda response | 0.20136 |
| B1: A + predicted classical to NQE to PIMD hierarchy | 0.21474 |
| B: A + both physics blocks | 0.21901 |
| C: integrate predicted dH/dlambda, then outer-train affine calibration | 1.51356 |

The response block was specified before the final teacher collection and was run
once at full coverage. It was not tuned after seeing these results.

The corresponding plot is `figures/structure_only_multilambda_ablation.png`.

## Response-head diagnostics

| Lambda | Component | N | MAE | Pearson r | MAE / target SD |
|---:|---|---:|---:|---:|---:|
| 0.1 | total dH/dlambda | 72 | 3.565 | 0.903 | 0.319 |
| 0.1 | Coulomb | 72 | 4.119 | 0.867 | 0.395 |
| 0.1 | van der Waals | 72 | 5.199 | 0.912 | 0.308 |
| 0.5 | total dH/dlambda | 74 | 2.350 | 0.542 | 0.694 |
| 0.5 | Coulomb | 74 | 4.028 | 0.774 | 0.495 |
| 0.5 | van der Waals | 74 | 3.524 | 0.903 | 0.331 |
| 0.9 | total dH/dlambda | 73 | 3.481 | 0.793 | 0.466 |
| 0.9 | Coulomb | 73 | 3.917 | 0.873 | 0.362 |
| 0.9 | van der Waals | 73 | 1.275 | 0.940 | 0.257 |

Energy-like response quantities are in ARBALEST's reported kcal/mol units. The
polarization column is identically zero in this probe export and therefore
contains no learnable signal.

## Conclusion

The student captures broad response ordering but not response values precisely
enough. Absolute component errors of 1.27-5.20 kcal/mol compound under integration
and are far larger than the final 0.20-kcal/mol target. Adding these predicted
observables therefore degrades rather than improves the matched experimental
head.

This identifies the missing information specifically: a larger, chemically
diverse, protocol-aligned set of full lambda curves with component-resolved
dH/dlambda and paired classical/PIMD8 labels. The highest-value acquisition is at
least 50 diverse benchmark-disjoint molecules in each current high-error family:
amides, aromatics, ethers, acids, and alkanes (roughly 250-500 total). More generic
structure-model tuning does not address the observed response-surrogate error.
