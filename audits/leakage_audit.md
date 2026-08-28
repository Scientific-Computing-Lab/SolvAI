# Independent SolvAI leakage and inference audit

**Status: PASS.** Exact-connectivity, fragment-parent, uncharged-parent and
canonical-tautomer equivalence have been audited. The raw source audit found
2 CombiSolv-QM, 32 MolSolv and 22 ConfSolv standardized equivalents; each was
removed before the released teachers were refitted. Retained split membership
was independently verified. The 1,280-label endpoint pool has zero standardized
benchmark equivalent.

| Source | Rows | Unique structures | SMILES overlap | Full-key overlap | Connectivity overlap | Verification |
|---|---:|---:|---:|---:|---:|---|
| Expanded public hydration | 5,075 | 5,075 | 0 | 0 | 0 | frozen clean-room audit; source not redistributed |
| Legacy public hydration | 1,147 | 1,147 | 0 | 0 | 0 | frozen clean-room audit; source not redistributed |
| CombiSolv-QM water | 3,963 | 3,963 | 0 | 0 | 0 | frozen clean-room audit; source not redistributed |
| SoluteML Abraham | 8,098 | 8,098 | 0 | 0 | 0 | independent release-time recanonicalization |
| OpenFF explicit alchemical | 520 | 520 | 0 | 0 | 0 | independent release-time recanonicalization |
| GBn2 implicit-solvent | 550 | 550 | 0 | 0 | 0 | independent release-time recanonicalization |
| MolSolv SMD(water) | 350,391 | 350,391 | 0 | 0 | 0 | independent release-time recanonicalization |
| ConfSolv water response | 39,878 | 39,878 | 0 | 0 | 0 | independent release-time recanonicalization |

Morgan-similarity-1 cases with different identities are enumerated rather than
silently removed; global similarity-exclusion controls are reported separately.

The released endpoint contains 2,265 deterministic RDKit/Morgan features and
15 structure-predicted physical-response priors. Its schema contains no
experimental, ARROW, PIMD, trajectory, probe, family, scaffold or fold field.
The runtime code does not open benchmark or prediction tables.

CombiSolv-QM and mixed public-hydration source rows are not redistributed
where publisher supplements have no standalone data licence. Their final
clean-room audit records, counts and source-table SHA-256 values are frozen
in this release; all included tables are recanonicalized when the audit is rerun.
