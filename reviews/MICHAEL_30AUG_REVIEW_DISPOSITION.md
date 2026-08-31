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
| ARROW outer-training weight 3 | **NEW SENSITIVITY** | Weight 3 was frozen prospectively in `release/CONFIRMATORY_FREEZE.md`, not chosen post hoc. A single predeclared weight-1 sensitivity will change only this weight and retain all folds, rows, teachers, features, model settings and seeds. |
| PIMD8 role | **ALREADY RESOLVED** | PIMD-derived features were not retained; PIMD8 is an accuracy comparator. The Abstract still needs one sentence-level clarification that SolvAI's endpoint is experimentally supervised whereas PIMD8 is a simulation comparator. |
| Corrected primary and repeat results | **ALREADY RESOLVED** | The canonical values are 0.20223 kcal mol\(^{-1}\) on the fixed partition and 0.20737 ± 0.00444 across five complete partitions. Historical 0.197/0.204 values are not the confirmatory headline. |
| Abstract comparison asymmetry | **MANUSCRIPT FIX** | Add minimal wording that SolvAI reaches the PIMD8 accuracy scale on this chemistry despite using experimental endpoint supervision, and do not imply equivalent information or PIMD distillation. |
| Runtime in the main paper | **MANUSCRIPT FIX** | Add the machine-readable CPU runtime: 15.294 s warm single-molecule median; 15.824 s for batch 32, or 0.494 s per molecule. Any PIMD work comparison must retain the published-2080-Ti versus measured-CPU caveat. |
| SMILES preprocessing | **MANUSCRIPT FIX** | `solv_ai.features.canonicalize` parses with RDKit and converts each query to canonical isomeric SMILES before descriptors and teacher inference. It does not claim tautomer, protonation or salt invariance. Add one precise Methods sentence. |
| Per-query applicability/reliability score | **MANUSCRIPT FIX** | The release returns ensemble spread but has no calibrated applicability-domain or reliability score. State this once as a limitation; do not build a new model. |
| Tier-A 221 external cohort | **NEW EXTERNAL VALIDATION** | First audit endpoint compatibility, standardized identity, experimental-source provenance and exposure to all six teacher sources without inspecting SolvAI errors. Freeze the eligible cohort and terminology before evaluating matched A versus F. |
| GAFF/AM1-BCC, ASP-19/GB and AMOEBA baseline ladder | **NOT ACTIONABLE** unless exact target/subset comparability is established | Classical ARROW and PIMD8 are already reconstructed on the same 85 molecules. Literature values with n=80 or n=33 and ASP-19/GB on a different cohort must not be plotted as if they were same-set causal baselines. Audit provenance; omit an unfair ladder. |
| Thiol/sulfide anecdote | **NOT ACTIONABLE** | A favourable small post-hoc family must not replace the complete family analysis. No new claim or selected example will be added. |
| Charged/ionic and peptide stress tests | **NOT ACTIONABLE** for the current validation; future work | The current scope is neutral organic hydration. Numerical outputs for ions, zwitterions or peptides are not validated applicability evidence. |
| Broad retuning or a new applicability model | **NOT ACTIONABLE** | Both would reopen model development and conflict with the frozen confirmatory design. |

## Governing interpretation

The remaining work is bounded: one weight sensitivity, one independently qualified
external-cohort evaluation if Tier-A permits it, and targeted transparency edits.
Neither analysis may alter teachers, feature selection, endpoint hyperparameters or
the current title. Negative or neutral results will be reported without rescue
tuning.
