# SolvAI Phase 1 confirmatory-analysis freeze

Frozen on 2026-08-28 before executing any analysis described below.

This protocol is prospective with respect to the Phase 1 confirmatory results. It
does not change the previously frozen SolvAI teachers, endpoint architecture,
benchmark targets, or exploratory results. Analysis definitions in this file will
not be changed in response to the observed results. Corrections required by an
implementation error will be documented as amendments, with the original text
retained in Git history, before the corrected analysis is run.

## Amendment 1: required execution environment

Added on 2026-08-28 before any valid confirmatory result was produced. The first
diagnostic invocation used the workstation base interpreter, which contains
scikit-learn 1.6.1 rather than the frozen SolvAI environment. Its outputs were
quarantined and are not confirmatory results. All analyses governed by this protocol
must be invoked with `uv run` from the parent `Freecurve_AI_Solvation` environment,
which provides Python 3.11.15, RDKit 2026.03.5 and scikit-learn 1.7.2. This amendment
does not alter a dataset, feature block, split, model hyperparameter, metric or
interpretation rule.

## Scientific questions

Phase 1 asks five primary questions.

1. Does the final response-prior stack improve over a genuinely matched molecular-
   structure-only endpoint?
2. Which predeclared response-source blocks contribute useful endpoint signal?
3. Does molecule--response alignment matter, or do arbitrary additional columns
   produce a similar result?
4. Does any response-prior advantage survive chemical separation applied to every
   endpoint-supervised molecule, including the public experimental-label pool?
5. Does the response representation retain value when no ARROW experimental label
   is used for endpoint fitting?

The analysis will not search architectures, tune hyperparameters, alter teachers,
launch simulations, or select a preferred split or feature combination after the
results are known.

## Immutable inputs

The canonical input files and SHA-256 hashes are:

| Input | Rows used or available | SHA-256 |
| --- | ---: | --- |
| `data/processed/arrow_solvation_master.parquet` | 85 unique water-solvation molecules after the frozen water/identity selection | `2b7928f162d094e7ee10d197e66636ba4ae09b0f76d626136b79c0975d3b0310` |
| `data/processed/expanded_public_hydration_nonbenchmark.parquet` | 5,075 candidate external records; the frozen source-quality mask retains exactly 1,280 | `603ed02b6be25d9a3057e321f2c6ea135b012666cfdb8a1b160e37f347951ec4` |
| `data/processed/expanded_public_hydration_rdkit_morgan_features.parquet` | 5,075 | `f6d9cd37a90bfc0718261f7251c70100be15947f73d7db12c85659ffc05b28e9` |
| `data/processed/rdkit_morgan_features.parquet` | 85 | `39877c3938616445d9996e093f8f37744a9c13d56202db9e294466396df298b4` |
| `data/processed/combisolv_qm_teacher_predictions.parquet` | 5,160 benchmark/public predictions | `b96472d73c889c53b190c5836beb2613745021f530026324dc4b90c7f576c3c0` |
| `data/processed/soluteml_abraham_teacher_predictions.parquet` | 5,160 | `d3bbcc28c4893e886564a8dda46d651ec83f2e202cd541e5038231ef91c34033` |
| `data/processed/openff_alchemical_teacher_predictions.parquet` | 5,160 | `b4a2f7ba810321d5daa112e019d52b6deae9923261c867fc3fa9a09f0549e841` |
| `data/processed/implicit_solvent_teacher_predictions.parquet` | 5,160 | `278a3071f95acddace4513d0794b7eed5a66b181c88adb523dd794a74a86af95` |
| `data/processed/molsolv_smd_teacher_predictions.parquet` | 5,160 | `a0412a23ab92b1036d9a7574d3c566d6a5065d868703505ed90fb81cad5fdd83` |
| `data/processed/confsolv_water_teacher_predictions.parquet` | 5,160 | `a30dc4a57be79720c58738c191f9764676754ffc9b8263df5b79ecd62803be5b` |

The 1,280-row endpoint mask is the frozen rule already used by SolvAI: retain every
connectivity in `public_hydration_nonbenchmark.parquet`, plus SoluteML-only rows in
the expanded table with `source_measurement_count >= 2`. No row will be selected or
removed using its prediction error.

The physics teachers and their predictions are frozen. All teacher source sets
already exclude every exact ARROW connectivity. Phase 1 does not refit or select a
teacher.

## Structure representation and endpoint model

Every matched endpoint uses the same deterministic structure representation:

- 2,048-bit radius-2 Morgan fingerprint;
- 217 RDKit two-dimensional descriptors;
- no static ARROW parameters;
- no conformer, trajectory, PIMD, MD, SMD, or probe calculation at inference.

Every endpoint member is `sklearn.ensemble.ExtraTreesRegressor` with:

- `n_estimators=360`;
- `max_features=0.7`;
- `min_samples_leaf=2`;
- `criterion="squared_error"`;
- `bootstrap=False`;
- `max_depth=None`;
- `min_samples_split=2`;
- `n_jobs` may vary without changing predictions;
- model seeds `11`, `29`, and `47`.

The three predictions are averaged. External experimental rows have sample weight
1.0. Available ARROW outer-training rows have sample weight 3.0. The test molecule's
experimental label is always absent. No block-specific tuning or nested selection is
performed in the primary confirmatory comparison.

## Predeclared feature blocks

The following six feature sets will be evaluated. Each is structure plus the listed
block; these are not selected after evaluation.

| ID | Feature set | Added predicted response coordinates |
| --- | --- | --- |
| A | Structure only | none |
| B | Empirical/residual-corrected | five Abraham `E/S/A/B/L` axes; corrected OpenFF response defined as predicted OpenFF hydration free energy plus predicted experimental residual; corrected GBn2 response defined analogously (7 coordinates) |
| C | Computation-derived core | CombiSolv-QM/COSMOtherm water response; raw predicted OpenFF hydration free energy; raw predicted GBn2 alchemical hydration free energy (3 coordinates) |
| D | SMD(water) | predicted MolSolv SMD(water) response only (1 coordinate) |
| E | ConfSolv | the six retained gas/solution/hydration conformer-correction and water-response distribution summaries (6 coordinates) |
| F | Full retained SolvAI | CombiSolv-QM (1), Abraham axes (5), corrected OpenFF (1), corrected GBn2 (1), SMD(water) (1), and ConfSolv summaries (6): 15 coordinates total |

For continuity with the frozen campaign, two cumulative stacks will also be
recomputed without selection: `narrow` (CombiSolv-QM, Abraham, corrected OpenFF,
corrected GBn2; 8 coordinates) and `narrow+SMD` (9 coordinates). They are secondary
descriptive results and do not replace comparisons A--F.

## Random-fold and repeated evaluation

The primary partition is the existing `fold_random` column in the 85-molecule
benchmark. Five complete repeat partitions use shuffled five-fold `KFold` with
seeds:

`314159, 271828, 161803, 141421, 173205`.

The benchmark row order in the immutable parquet file is retained. All six feature
sets use identical folds and training rows. The public 1,280-row pool is present in
every random-fold fit. These repeated partitions measure partition sensitivity; they
are not an independent external test and will not be described as one.

## Shuffled-prior negative control

The complete 15-column response block is shuffled as a row-aligned block, preserving
the joint distribution and inter-prior correlations while destroying alignment to
the molecule.

For every outer fold:

- the 1,280 public response rows are permuted among public training rows;
- response rows for ARROW outer-training molecules are permuted only among those
  outer-training molecules;
- response rows for outer-test molecules are permuted only among outer-test
  molecules;
- no experimental target participates in a permutation.

Five independent permutation replicates use base seeds
`88001, 88002, 88003, 88004, 88005`; the deterministic per-fold seed is
`base_seed + 1000 * repeat_index + outer_fold`, with `repeat_index=0` for the frozen
primary partition and 1--5 for the five repeat partitions. The shuffled control uses
the same endpoint model seeds and hyperparameters as A--F.

## Global chemical-separation evaluation

Chemical separation applies to all endpoint-supervised molecules. The frozen
physics teachers are not refitted: their scientific role is external response
pretraining, and exact ARROW connectivity was already excluded globally. For each
outer test partition, public experimental rows and ARROW training rows matching the
held-out chemistry are excluded before endpoint fitting.

### Functional-family separation

A deterministic primary family is assigned from RDKit structure to both public and
ARROW molecules using the following precedence:

1. carboxylic acid;
2. amide;
3. ester;
4. aldehyde;
5. ketone;
6. alcohol/phenol;
7. ether;
8. amine;
9. thiol;
10. sulfide/disulfide;
11. aromatic heterocycle;
12. carbocyclic aromatic;
13. alkene/diene;
14. saturated hydrocarbon;
15. other.

The implementation will store the exact SMARTS and assigned family for every row.
Five-fold `GroupKFold` is constructed on the 85 ARROW primary-family labels. For an
outer fold, every public label and every nominal ARROW training row assigned to a
family present in the outer test fold is removed.

### Bemis--Murcko scaffold separation

RDKit Bemis--Murcko scaffolds are canonicalized to isomeric SMILES. Molecules with an
empty Murcko scaffold receive `ACYCLIC::<primary-family>`. Five-fold `GroupKFold` is
constructed from the 85 ARROW scaffold keys. For each fold, every public label and
nominal ARROW training row with a test scaffold key is removed.

### Molecular-cluster separation

The combined 1,280-public-plus-85-ARROW endpoint pool is sorted by connectivity key
and clustered once using Taylor--Butina clustering of radius-2, 2,048-bit Morgan
fingerprints at Tanimoto similarity 0.70 (distance cutoff 0.30). Five-fold
`GroupKFold` is constructed from the ARROW cluster identifiers. Every endpoint row
sharing a held-out cluster is excluded from training.

### Strict nearest-neighbour exclusion

Using the frozen five random ARROW folds, every public or nominal ARROW training row
with Tanimoto similarity greater than or equal to 0.70 to any outer-test molecule is
removed. Prespecified sensitivity analyses use thresholds 0.50, 0.60 and 0.80. The
0.70 result is primary; no threshold is selected by endpoint performance.

All chemical-separation analyses compare A (structure only) with F (full SolvAI) on
identical retained training rows. Training counts, removed counts, and the maximum
remaining train--test similarity are reported per fold.

## Chemical-distance and leakage audit

The audit covers all endpoint labels and every final teacher source. It will report:

- canonical isomeric SMILES, full InChIKey, and connectivity-block matches;
- RDKit MolStandardize cleanup, largest-fragment selection, uncharging, canonical
  tautomer generation, and resulting standardized-identity matches;
- formal-charge and protonation relationships where a common standardized parent is
  detected;
- Bemis--Murcko scaffold matches;
- maximum radius-2 Morgan Tanimoto similarity for each ARROW molecule against each
  endpoint and teacher source;
- every similarity-1.0 non-identical pair, including structures and identifiers.

The audit is descriptive. No molecule is removed after inspecting its error. Any
newly discovered exact, tautomeric, or protonation-equivalent supervised overlap is
reported first; affected primary analyses will then be repeated using the
predeclared standardized-identity exclusion without changing any other definition.

## Zero-ARROW-label transfer

Models A and F are fitted once using only the same 1,280 external experimental
labels, weight 1.0, and endpoint seeds 11, 29 and 47. All 85 ARROW molecules are
predicted without any ARROW experimental label in training. This is a frozen-method
transport analysis, not an independent model-development test, because the method
was historically developed with reference-set feedback.

## Optional analyses, executed only after the primary package

### Experimental-label learning curve

If runtime remains modest, public endpoint-label sizes are
`80, 160, 320, 640, 1280`. For each size, five samples without replacement use seeds
`314159, 271828, 161803, 141421, 173205`. The primary frozen five-fold ARROW
partition and all available ARROW outer-training labels remain fixed. A and F are
compared using identical sampled rows. No curve point is used for model selection.

### Teacher fidelity versus downstream benefit

Existing frozen teacher test predictions only will be used; no teacher is refitted.
For each available scalar target, report MAE, RMSE, Pearson correlation, Spearman
correlation, and MAE divided by the target 5th--95th percentile range. Downstream
effect is the predeclared A--F endpoint change. Because sources, scales, and model
classes differ, this is descriptive and will not be presented as a causal regression
or universal law.

## Metrics and uncertainty

Primary metric: MAE in kcal mol\(^{-1}\).

Secondary metrics: RMSE, median absolute error, coefficient of determination, and
the fraction of molecules with lower absolute error than A.

Paired uncertainty uses 100,000 molecule-level bootstrap resamples with seed
`20260828`. For repeated partitions, each molecule's absolute error is first averaged
across the five repeats, then the 85 molecule-level paired differences are
bootstrapped. Repeat MAEs are also reported individually and as mean plus sample
standard deviation. Chemical-separation comparisons use paired bootstrap resampling
of the 85 held-out predictions.

No correction for multiple testing will be used to manufacture binary discoveries.
Every predeclared effect and interval will be reported.

## Interpretation rules

For F versus A and for each source block versus A:

- **positive:** the paired 95% bootstrap interval for the MAE difference is below
  zero;
- **negative:** the paired 95% bootstrap interval is above zero;
- **neutral:** the interval includes zero.

An effect is additionally called **material** only if the absolute MAE difference is
at least 0.010 kcal mol\(^{-1}\). Repeat consistency is reported independently as the
number of the five partitions in which the response model has lower MAE.

For the shuffled control, support for molecule-aligned response information requires
F to outperform the mean shuffled-prior control with a paired 95% interval below
zero. A shuffled result statistically indistinguishable from F is evidence against
the response-alignment claim.

For globally separated evaluations, the thesis is strengthened if F retains a
negative paired MAE difference from A; it is weakened if the interval includes zero;
and the transfer claim is falsified for that regime if the interval is above zero.
The absolute MAE is reported but is not required to remain below 0.20.

For zero-ARROW transfer, the analysis supports transport of the learned response
representation only if F improves on A. Because the method was developed using the
85-molecule reference set, this result will not be called an independent external
validation.

## Conclusions these analyses cannot support

Regardless of outcome, Phase 1 will not establish:

- universal hydration accuracy;
- state-of-the-art performance on FreeSolv or another external benchmark;
- superiority over PIMD8;
- that PIMD8 was distilled into the final model;
- that training was simulation-free;
- that the 15 priors contain information not mathematically determined by molecular
  structure;
- that five split repeats are five independent experimental datasets;
- transfer to ions, salts, proteins, arbitrary solvents, or substantially larger
  molecules.

## Outputs

All new artifacts are written under clearly separated paths:

- `results/confirmatory/` for predictions, metrics, configurations and bootstrap
  outputs;
- `audits/confirmatory/` for identity and chemical-distance audits;
- `reports/CONFIRMATORY_ANALYSIS.md` for the final Phase 1 interpretation;
- `figures/confirmatory/` for diagnostic figures only.

The existing frozen campaign artifacts are read-only inputs and will not be
overwritten.
