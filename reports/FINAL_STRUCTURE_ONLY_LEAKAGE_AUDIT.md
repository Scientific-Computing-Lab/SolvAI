# Final structure-only leakage audit

## Verdict

**PASS.** All external supervision sources used by the confirmed student have
zero connectivity overlap with the 85-molecule ARROW benchmark. The headline
prediction table has exactly one outer-fold prediction per molecule and no row
requires a simulated test-molecule property.

## Identity audit

| Source | Rows | Unique connectivities | Benchmark overlap |
|---|---:|---:|---:|
| Expanded public hydration | 5,075 | 5,075 | 0 |
| CombiSolv-QM water | 3,963 | 3,961 | 0 |
| MolSolv SMD(water) | 350,391 | 348,587 | 0 |
| ConfSolv water response | 39,878 | 39,878 | 0 |
| SoluteML Abraham | 8,098 | 8,098 | 0 |
| OpenFF explicit alchemical | 520 | 520 | 0 |
| GBn2/learned implicit | 550 | 550 | 0 |

Identity exclusion uses the first block of the standard InChIKey, deliberately
collapsing stereochemical aliases and thereby applying a stricter criterion than
exact canonical-SMILES matching.

## Inference boundary

Allowed inputs are SMILES and deterministic RDKit/Morgan features. The SMD,
CombiSolv-QM, Abraham, OpenFF, implicit-solvent, and ConfSolv quantities are
predictions from benchmark-excluded structure models. No experimental label,
classical ARROW free energy, PIMD value, trajectory, measured short-probe
observable, family label, or scaffold label is consumed for a query molecule.

The post-evaluation deployable refit is audited separately from the OOF result;
its presence never substitutes a training-set score for the reported OOF MAE.
