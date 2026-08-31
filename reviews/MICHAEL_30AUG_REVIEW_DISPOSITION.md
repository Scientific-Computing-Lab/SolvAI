# Disposition of Michael Levitt's 30 August 2026 review

**Audit date:** 31 August 2026

**Repository HEAD audited:** `1a8931bcec439e67c1f1cf3844887434827a9085`

**Review source:** `Review_Gal_draft_30Aug2026.docx` was not found in the local
workspace. This disposition therefore covers every substantive point supplied with
the execution request and compares it with the current manuscript rather than an
older draft.

The title remains **“Structure-predicted solvent responses enable simulation-free
hydration free-energy prediction”** by fixed author decision.

| Review point | Disposition | Current evidence or required action |
|---|---|---|
| Expanded identity and standardized-equivalence audit | **ALREADY RESOLVED** | Current Results, Methods, Supplementary Methods and `audits/confirmatory/` report exact, fragment-parent, uncharged-parent and canonical-tautomer checks. The affected teachers were refitted with original split membership preserved. Do not duplicate this text. |
| Zero-ARROW-label transfer | **ALREADY RESOLVED** | The frozen analysis gives 0.385 MAE for matched structure-only and 0.257 for SolvAI, paired change -0.128 kcal mol\(^{-1}\) with 95% CI [-0.246, -0.037]. It is correctly described as transport without reference-label adaptation, not an independent external benchmark. |
| Family, scaffold, cluster and nearest-neighbour separation | **ALREADY RESOLVED** | Current Results and Fig. 3 apply separation to all endpoint-supervised rows, not only ARROW rows. The relative response-prior advantage survives all predeclared regimes, while absolute family/scaffold extrapolation remains weaker. |
| ARROW outer-training weight 3 | **NEW SENSITIVITY — RESOLVED** | Weight 3 was frozen prospectively in `release/CONFIRMATORY_FREEZE.md`, not chosen post hoc. Under the committed weight-1 sensitivity, matched structure-only gives 0.30977 MAE and SolvAI 0.20642; paired change -0.10335, 95% CI [-0.21883, -0.01995]. The advantage survives. |
| PIMD8 role | **ALREADY RESOLVED** | PIMD-derived features were not retained; PIMD8 is an accuracy comparator. The Abstract still needs one sentence-level clarification that SolvAI's endpoint is experimentally supervised whereas PIMD8 is a simulation comparator. |
| Corrected primary and repeat results | **ALREADY RESOLVED** | The canonical values are 0.20223 kcal mol\(^{-1}\) on the fixed partition and 0.20737 ± 0.00444 across five complete partitions. Historical 0.197/0.204 values are not the confirmatory headline. |
| Abstract comparison asymmetry | **MANUSCRIPT FIX** | Add minimal wording that SolvAI reaches the PIMD8 accuracy scale on this chemistry despite using experimental endpoint supervision, and do not imply equivalent information or PIMD distillation. |
| Runtime in the main paper | **MANUSCRIPT FIX** | Add the machine-readable CPU runtime: 15.294 s warm single-molecule median; 15.824 s for batch 32, or 0.494 s per molecule. Any PIMD work comparison must retain the published-2080-Ti versus measured-CPU caveat. |
| SMILES preprocessing | **MANUSCRIPT FIX** | `solv_ai.features.canonicalize` parses with RDKit and converts each query to canonical isomeric SMILES before descriptors and teacher inference. It does not claim tautomer, protonation or salt invariance. Add one precise Methods sentence. |
| Per-query applicability/reliability score | **MANUSCRIPT FIX** | The release returns ensemble spread but has no calibrated applicability-domain or reliability score. State this once as a limitation; do not build a new model. |
| Tier-A 221 external cohort | **NEW EXTERNAL VALIDATION — RESOLVED** | Qualification was completed without model errors and committed before prediction. One pre-existing registry conflict was excluded; all 220 retained molecules are endpoint-disjoint and 97 are also disjoint from all teacher-source tables. MAE changes from 1.53165 to 1.15255 on N=220 and from 2.13830 to 1.53560 on strict N=97, with both paired intervals below zero. Absolute accuracy is substantially worse than on ARROW-85 and must be reported. |
| GAFF/AM1-BCC, ASP-19/GB and AMOEBA baseline ladder | **NOT ACTIONABLE** | Classical ARROW and PIMD8 are reconstructed on the same 85 molecules and remain valid context. The proposed GAFF/AM1-BCC and AMOEBA values use different subsets (reported n=80 and n=33), while ASP-19/GB is evaluated on a different cohort and target reconstruction. None can be placed on a common same-set ladder or support a “5x” statement without conflating subset and protocol effects. |
| Thiol/sulfide anecdote | **NOT ACTIONABLE** | A favourable small post-hoc family must not replace the complete family analysis. No new claim or selected example will be added. |
| Charged/ionic and peptide stress tests | **NOT ACTIONABLE** for the current validation; future work | The current scope is neutral organic hydration. Numerical outputs for ions, zwitterions or peptides are not validated applicability evidence. |
| Broad retuning or a new applicability model | **NOT ACTIONABLE** | Both would reopen model development and conflict with the frozen confirmatory design. |

## Governing interpretation

The bounded new analyses are complete. Both preserve the matched response-prior
advantage, but Tier-A also exposes a sharp scope limit: the relative gain transfers
while absolute MAE rises above 1 kcal mol\(^{-1}\). The manuscript should therefore
add the external causal contrast without extending the PIMD8-level absolute-accuracy
claim beyond ARROW-85. Teachers, feature selection, endpoint hyperparameters and the
title remain unchanged.

## Baseline-ladder audit

The proposed ladder is not suitable for a common quantitative axis. The valid
same-chemistry comparisons remain classical ARROW (85 solutes, 0.78465 kcal
mol\(^{-1}\)) and ARROW/PIMD8 (85 solutes, 0.20484 kcal mol\(^{-1}\)), both
reconstructed molecule by molecule from the ARROW source table. The proposed
GAFF/AM1-BCC value is reported for only the 80 ARROW identities that overlap the
FreeSolv calculation table, and the proposed AMOEBA value covers 33 molecules; those
subsets cannot be treated as the same 85-solute test. The ASP-19/GB archive reports
an author aggregate of 0.73257 kcal mol\(^{-1}\) on 85 names, but the independently
reconstructed disclosed protocol gives 0.67724--0.67959 and the archive does not
contain the row-level structure/prediction file needed to resolve that discrepancy
(`docs/V1.2_BENCHMARK_RECONSTRUCTION.md` in the separately audited ASP repository).
It is therefore omitted rather than
approximated. None of these values is used to replace the matched 0.30335
structure-only control, which remains the causal baseline for the response layer.
