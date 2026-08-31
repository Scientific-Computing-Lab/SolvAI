# Tier-A external-cohort validation

## Result in one paragraph

The pre-existing Tier-A cohort from Sander's Henry-law compilation contained 221
neutral organic solutes and was assembled without SolvAI predictions. A previously
documented registry name/CAS/structure conflict excluded one row before SolvAI
evaluation. All remaining 220 molecules are absent, exactly and under the frozen
standardized-equivalence rules, from the 1,280 external and 85 ARROW endpoint labels.
Full SolvAI reduces MAE from **1.53165 to 1.15255 kcal mol-1** (paired change,
**-0.37909**; 95% interval, **-0.51692 to -0.24817**). Ninety-seven molecules are
also absent from all six supervised response-source tables; on this strict subset,
MAE falls from **2.13830 to 1.53560 kcal mol-1** (paired change, **-0.60271**; 95%
interval, **-0.86277 to -0.35858**). The response-prior advantage therefore transfers,
but absolute accuracy degrades substantially outside the compact ARROW reference
chemistry.

## Cohort qualification

The source is Sander's Henry-law compilation v5.0.0 (DOI
`10.5194/acp-23-10901-2023`). Tier-A uses original-measurement (`M`) records in pure
water, normalized to 298.15 K. The frozen conversion is
`Kc = Hscp R T` and `delta_G = -R T ln(Kc)`, which yields the 1 M ideal-gas to 1 M
ideal-dilute aqueous convention with negative values favouring hydration. The source
prints Henry constants to two significant figures; this limits target precision.

Eligibility was based only on chemistry, identity, target convention and provenance.
It was fixed in `release/TIER_A_EXTERNAL_VALIDATION_FREEZE.md` at commit `b01351a`
before any SolvAI or matched-baseline prediction was generated.

| Stage | N | Decision |
|---|---:|---|
| Original ASP Tier-A freeze | 221 | Pre-existing error-blind cohort |
| Compatible neutral molecular solutes | 221 | All passed target and molecular-scope checks |
| Registry-consistent, endpoint-disjoint cohort | 220 | One pre-existing name/CAS/structure conflict excluded; zero endpoint overlaps |
| Strict response-source-disjoint subset | 97 | Additionally absent from every supervised teacher-source table |

Identity checks used canonical isomeric SMILES, full InChIKey, connectivity block,
fragment parent, uncharged parent and canonical tautomer. Every row and every match is
retained in `results/tier_a_external/qualification/`. No molecule was removed based
on prediction error.

## Response-source exposure

Teacher-source membership is not endpoint-label leakage, but it can make response
prediction easier for a molecule. It is therefore reported explicitly.

| Response source | Exact overlap among N=220 | Any standardized overlap |
|---|---:|---:|
| CombiSolv-QM | 24 | 25 |
| Abraham/SoluteML | 92 | 92 |
| OpenFF | 0 | 0 |
| GBn2 | 0 | 0 |
| MolSolv SMD(water) | 53 | 54 |
| ConfSolv H2O | 14 | 14 |

In total, 123 endpoint-disjoint molecules occur in at least one teacher source and
97 occur in none. Because the source compilation is distinct but complete
measurement-level provenance cannot rule out every historical database relationship,
the appropriate term is **external molecule-disjoint cohort**, not proof of wholly
independent experimental provenance.

## Frozen matched comparison

Both endpoints were refitted on exactly the same 1,280 external experimental labels
plus all 85 ARROW labels, with weights 1 and 3 respectively. They use the same three
360-tree ExtraTrees members, seeds, imputation and structure representation. The sole
difference is the 15 frozen response coordinates. Tier-A labels never entered fitting,
feature selection, teacher selection, thresholding or model choice. The refitted full
endpoint reproduced the released artifact's Tier-A predictions to a maximum absolute
difference of `3.553e-15` kcal mol-1.

| Cohort | Method | N | MAE | RMSE | Median absolute error |
|---|---|---:|---:|---:|---:|
| Endpoint-disjoint | Matched structure-only | 220 | 1.53165 | 2.30974 | 0.96182 |
| Endpoint-disjoint | Full SolvAI | 220 | 1.15255 | 1.57723 | 0.84508 |
| Strict source-disjoint | Matched structure-only | 97 | 2.13830 | 3.08051 | 1.10237 |
| Strict source-disjoint | Full SolvAI | 97 | 1.53560 | 1.98792 | 1.24392 |

| Cohort | Paired MAE change | 95% paired bootstrap interval | Molecules improved |
|---|---:|---:|---:|
| Endpoint-disjoint | -0.37909 | [-0.51692, -0.24817] | 60.9% |
| Strict source-disjoint | -0.60271 | [-0.86277, -0.35858] | 58.8% |

The strict subset's median absolute error is slightly worse with SolvAI even though
its mean absolute error and RMSE improve substantially. The gain is therefore broad
enough to improve 58.8% of molecules but is amplified by correction of several large
structure-only errors. This nuance is retained rather than summarized as universal
per-molecule improvement.

## Chemical-distance diagnostics

The endpoint-disjoint cohort is chemically more difficult than ARROW-85. Its median
maximum Morgan similarity to the 1,365 endpoint-training molecules is 0.441. SolvAI's
largest average improvement occurs in the lowest-similarity quartile (MAE 3.116 to
2.057) and for the largest heavy-atom quartile (3.287 to 2.166). These are descriptive
strata; no threshold was selected.

Two molecules have Morgan similarity 1.0 to a different endpoint-training molecule:
ethyl octanoate and ethyl decanoate collide with a longer homolog under the fixed
radius-2 fingerprint. Their identities remain distinct under every exact and
standardized check, so they are retained and the collision is disclosed.

Large remaining SolvAI errors concentrate in chemistry that is sparse in the endpoint
training set, notably perfluorinated chains, PFBHA derivatives and sulfonic-acid or
sulfonamide structures. This is a post-hoc diagnostic, not an eligibility rule.

## Interpretation

The external result strengthens the causal claim that structure-predicted response
coordinates add information beyond the same endpoint trained on structure alone. It
does not extend the PIMD8-level absolute-accuracy claim beyond ARROW-85. On the broad
Tier-A chemistry, absolute errors of 1.15--1.54 kcal mol-1 remain too large for the
precision established on the compact reference set. The appropriate conclusion is
therefore two-part: the response advantage transfers to new molecules and survives
complete response-source disjointness, while calibrated broad-domain prediction
remains an open limitation.

Machine-readable predictions, metrics, response features, descriptive strata and
qualification records are under `results/tier_a_external/`.
