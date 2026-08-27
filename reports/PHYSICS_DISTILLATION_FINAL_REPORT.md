# Structure-only ARROW/PIMD physics distillation: final report

## Decision

The structure-only model crosses 0.20 on one leakage-safe five-fold partition:
the fixed SMD+ConfSolv student reaches **0.19705 kcal/mol**, and nested
feature-block selection reaches **0.19931**. The deployed artifact consumes
SMILES only and runs no simulation.

The crossing is **not robustly confirmed**. Across five independent outer-CV
partitions, the fixed student averages **0.20375 ±
0.00493**, nested selection averages **0.20860**, and
the release-time molecule bootstrap 95% interval is **[0.14810,
0.25226]**. Family and scaffold
holdouts remain **0.23957** and **0.24128**. The scientifically honest
answer is therefore: a sub-0.20 point estimate exists, but robust sub-0.20
generalization has not been demonstrated.

## Headline comparison

| Method | Simulation at inference? | Fixed random OOF | Five-repeat mean | Family | Scaffold |
|---|---:|---:|---:|---:|---:|
| ARROW/PIMD8 | Yes | 0.20484 | — | — | — |
| PIMD8 + nested residual | Yes | 0.18682 | — | 0.19644 | 0.20056 |
| Previous structure-only baseline | No | 0.23861 | — | 0.31698 | 0.32892 |
| Structure + narrow response | No | 0.21362 | 0.21545 | — | — |
| + MolSolv SMD teacher | No | 0.20161 | 0.20396 | — | — |
| + ConfSolv response (fixed) | No | **0.19705** | **0.20375** | 0.23957 | 0.24128 |
| Nested SMD/ConfSolv selection | No | **0.19931** | 0.20860 | — | — |

## What transferred

MolSolv contributes 350,391 benchmark-disjoint SMD(water) structures; adding its
water-response teacher improves the fixed model from 0.21362 to 0.20161. ConfSolv
adds 39,878 benchmark-disjoint water structures with conformer-response moments and
improves the fixed result further to 0.19705. This is the strongest direct evidence
that physics-rich supervision transfers into a structure-only student.

Broader physical blocks did not help this already strong representation. Predicted
OpenFE diagnostics, the MLFF/force-field hierarchy, and DES370K water/SAPT response
score 0.20575, 0.20109, and 0.20993 in matched one-seed screens versus 0.19191 for
their unchanged base. They were stopped rather than combined or tuned.

## Lambda-response experiment

- Multi-lambda physics distillation A: structure/response baseline: **0.19592 kcal/mol**
- Multi-lambda physics distillation B2: +distilled PIMD2 lambda response: **0.20136 kcal/mol**
- Multi-lambda physics distillation B1: +distilled classical-NQE-PIMD hierarchy: **0.21474 kcal/mol**
- Multi-lambda physics distillation B: +full distilled physics hierarchy: **0.21901 kcal/mol**
- Multi-lambda physics distillation C: integrated predicted dH/dlambda (affine calibration): **1.51356 kcal/mol**

Nonconstant response-head MAE ranges from 1.27 to 5.20 kcal/mol. All measured PIMD2 values are training-only labels. For every
outer-test molecule, both its response curve and classical/PIMD hierarchy labels
are withheld; the final prediction uses only its structure.

## Leakage and deployment

The external sources used by the final artifact have zero connectivity overlap
with all 85 benchmark molecules. The packaged feature audit finds zero forbidden
test-time fields. End-to-end inference reconstructs all 15 learned physics features
from SMILES within 1.9e-6 of their cached values. The artifact is under
`models/final_structure_only/`; the CLI is `scripts/predict_structure_only.py`.

## Remaining blocker and exact next data

The limiting step is not the final regressor. It is structure-to-response transfer
from sparse, chemically narrow high-fidelity labels. The largest errors are amides,
aromatics, ethers, acids, and alkanes, while current independent classical/PIMD or
lambda-response sets have low nearest-neighbor similarity and only tens of examples
per family. The next acquisition should be a protocol-matched, benchmark-disjoint
set containing full lambda-resolved dH/dlambda plus electrostatic, polarization,
dispersion/repulsion, and classical/PIMD8 pairs for at least 50 diverse molecules
in each of those five families (roughly 250-500 molecules total). That directly
targets the observed response-surrogate error; more generic SMILES tuning does not.

## Required conclusions

1. Best zero-simulation experimental MAE: 0.19705 fixed OOF; 0.19931 selection-adjusted nested OOF.
2. Best family-held-out MAE: 0.23957.
3. Raw ARROW/PIMD8 MAE on comparable data: 0.20484.
4. PIMD8 + ML residual MAE: 0.18682 (not eligible for deployment).
5. PIMD-distilled MAE: 0.29169 for the earlier PIMD-only student.
6. Final ensemble MAE: 0.19705 on the original split; 0.20375 five-repeat mean.
7. Number of PIMD8 labels used: 85 in candidate distillation experiments; 0 in the selected artifact because that block worsened OOF validation.
8. Number of NEW PIMD8 simulations performed: 0.
9. Estimated inference speedup vs PIMD8: not used as a claim; inference is structure-model evaluation only.
10. Does the project beat 0.20 legitimately? NO — one strict split crosses, but repeated/hard validation does not.
11. Does it beat 0.20 without simulation at inference? NO robustly; YES only for the original strict nested point estimate.
12. Most important model component: benchmark-disjoint SMD(water) plus ConfSolv conformer-response distillation.
13. Main failure family: amides.
14. Best next experiment if given one additional day: acquire protocol-matched full lambda/component and paired classical/PIMD labels for the five high-error families above.
