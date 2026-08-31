# Prospective freeze: Tier-A external-cohort validation

Frozen on 31 August 2026 after cohort qualification and before calculating any
SolvAI or matched-baseline prediction for Tier-A. The cohort was defined without
access to SolvAI errors. This analysis does not alter the frozen ARROW-85 model,
teachers, feature selection or confirmatory evaluation.

## Scientific question

Does the advantage of the 15 predicted solvent-response coordinates over an
otherwise identical structure-only endpoint survive on a pre-existing external
cohort whose experimental endpoint labels were not used to train either endpoint?

## Source and target

- Original cohort: 221 Tier-A molecules from the independent ASP-19/GB validation
  project, frozen before that project's predictions.
- Primary source: Sander Henry-law compilation v5.0.0, DOI
  `10.5194/acp-23-10901-2023`, original-measurement (`M`) records only.
- Solvent and condition: pure water, normalized to 298.15 K.
- Conversion: `Kc = Hscp R T` and `delta_G = -R T ln(Kc)`, giving the 1 M ideal-gas
  to 1 M ideal-dilute aqueous convention in kcal mol-1; negative values favour
  hydration.
- Original Tier-A input SHA-256:
  `7b7be03c5124559551c5380b0429d5c0156c13060af9902d26266510637fcf45`.
- ASP source commit: `8f1a0c6c5429c27e9fda96422763fbf35becc2ca`.

The source is a distinct compilation, but shared underlying measurements cannot be
excluded from provenance alone. Accordingly, the strongest allowed term is
**external molecule-disjoint cohort**, not proof of experimental-source
independence.

## Prospective eligibility and exclusions

Eligibility was fixed from chemistry, identity, target compatibility and provenance,
never from SolvAI error:

1. valid RDKit molecule; one connected neutral molecular solute;
2. H, C, N, O, F, P, S, Cl, Br or I only;
3. finite original-measurement water target at 298.15 K with the convention above;
4. one qualified target per connectivity;
5. no exact or standardized-equivalent identity in the 1,280 external endpoint
   labels or 85 ARROW endpoint labels;
6. exclusion of `MKERQGLKSFEKAE`, whose source name/CAS/structure conflict was
   documented in ASP protocol amendment 002 before this SolvAI analysis.

Identity comparison uses canonical isomeric SMILES, full InChIKey and connectivity,
then RDKit FragmentParent, Uncharger and canonical-tautomer connectivity. No error-
based or performance-based exclusion is permitted.

## Frozen cohorts

| Cohort | Definition | N | SHA-256 (CSV) |
|---|---|---:|---|
| Original | ASP Tier-A freeze | 221 | `7b7be03c5124559551c5380b0429d5c0156c13060af9902d26266510637fcf45` |
| Endpoint-disjoint | Scientifically eligible and absent, exactly or by standardized equivalence, from every endpoint label used in deployment fitting | 220 | `967a9794a5d0e3f131dc5bb6921fecd234809ca21143a5a3ce9ed7895d18b273` |
| Strict response-source-disjoint | Endpoint-disjoint and additionally absent, exactly or by standardized equivalence, from all six supervised response-source tables | 97 | `486040f54082f606ca21f461eb76075648b885c712823c72f7f69317960463bf` |

No scientifically eligible row overlapped an endpoint label. Among the 220
endpoint-disjoint rows, 123 occur in at least one response-teacher source and 97 do
not. Source-specific exposure is frozen in
`results/tier_a_external/qualification/tier_a_teacher_exposure_summary.csv`.

Frozen qualification artifacts:

- complete audit CSV:
  `a32f40de3befe18a1bacfc299a93fbc9afb5972e2828b1c1ff2bd52b94339252`;
- all identity matches:
  `f18a33eb404b2c9b61ef9d5b60a35448d242e6a7fc7fd7fd25060f9137b4feae`;
- endpoint-disjoint Parquet:
  `a364da2a9fbf4603162cb40a58e0ccf233f073e5ed3d98ca79f9fccaa391eaac`;
- strict Parquet:
  `13ee38d907d9b9afb82bb6a390aa96e1b475d9ef214ab8c96a7129fb4553d414`.

## Frozen model comparison

Both endpoints are fitted deployment-style on exactly the same 1,280 external
experimental labels plus all 85 ARROW labels. External rows have weight 1 and ARROW
rows weight 3, as frozen before the confirmatory campaign. Tier-A labels never enter
training, feature selection, thresholding, teacher selection or model choice.

- **Matched structure-only:** 2,048-bit radius-2 Morgan fingerprint plus 217 RDKit
  descriptors.
- **Full SolvAI:** the identical 2,265 structure features plus the frozen 15
  structure-predicted response coordinates.
- **Endpoint:** three identically configured ExtraTrees pipelines, 360 trees each,
  `max_features=0.7`, `min_samples_leaf=2`, squared-error criterion, no bootstrap,
  unlimited depth, minimum split 2, seeds 11, 29 and 47; mean prediction.
- **Preprocessing:** median imputation with missingness indicators, identical between
  comparators.
- **Teachers:** the six released standardized-exclusion artifacts. They remain
  frozen and are not refitted for Tier-A. The strict subset addresses direct
  molecule exposure to teacher-source data.

Released teacher SHA-256 values:

| Artifact | SHA-256 |
|---|---|
| `combisolv_qm.pt` | `2950be1232e2f868136e3388cf139e671d4196922ece0d3c7043d9dadece8b40` |
| `molsolv_smd.pt` | `35601030f030af72a9e6990bf8003df9d481858b27994b3fc70c81d605d130be` |
| `abraham_teacher.joblib` | `8c8a6365e7162dbe6ae445ce4004928ed86b2669d3e7a70d3d054f143ca97554` |
| `openff_teacher.joblib` | `55b7d31a293dd61fe57f5505e015775d7d86f53dca69980d398843706d58bb9b` |
| `implicit_teacher.joblib` | `484d5e0c1efa759c8e7acbfcbbd87e218b93ae239e7f35303fe4cd7e310910f8` |
| `confsolv_teacher.joblib` | `a71487e53d5085e897dad02e62e95562933bc47262145fbf44286aae04b70b0d` |

The evaluated code and frozen artifact state are based on SolvAI commit
`9851ed2154be5a8be1f41f3d07975e20fec9a900`, plus this prospective freeze and its
qualification artifacts. No scientific model change is permitted after the result
is observed.

## Metrics and descriptive analyses

For both the endpoint-disjoint cohort (N=220) and strict subset (N=97), report:

- MAE (primary), RMSE and median absolute error, kcal mol-1;
- fraction of molecules for which full SolvAI has smaller absolute error;
- paired mean absolute-error difference, `abs(error_full) - abs(error_structure)`;
- two-sided 95% interval from 100,000 molecule-level paired bootstrap resamples,
  seed 20260828.

Molecule-level predictions and all qualified rows will be published. Error versus
heavy-atom count and maximum Morgan similarity to the endpoint-training pool may be
summarized descriptively without selecting thresholds or excluding rows.

## Interpretation fixed in advance

- **Strong transfer evidence:** the paired interval is below zero in both the N=220
  endpoint-disjoint cohort and N=97 strict response-source-disjoint subset.
- **Endpoint-transfer only:** N=220 is positive but the strict interval includes
  zero. Claims must be limited to new endpoint labels, with teacher exposure noted.
- **Neutral:** the relevant interval includes zero.
- **Negative:** the relevant interval lies above zero; structure-only is better.

Absolute MAE may be worse than on ARROW-85 without invalidating the paired question.
Every outcome will be reported. This external analysis cannot establish universal
chemical generalization, and it does not alter the frozen ARROW-85 estimate.

The frozen evaluation command will be:

```bash
uv run python scripts/run_tier_a_external_validation.py
```
