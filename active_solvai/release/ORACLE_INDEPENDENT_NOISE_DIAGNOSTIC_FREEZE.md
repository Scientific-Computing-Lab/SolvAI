# Post-hoc independent-noise oracle diagnostic freeze

**Frozen before any cross-block performance calculation:** 2026-09-03 UTC  
**Immutable no-go source:** `8fb984c2eb26d016c6b81cf488f88dc667ca9cd3`  
**Hostile-audit commit:** `19f1a584b51be8a63a90d37e8ca29d710201abf6`  
**Diagnostic branch:** `active-solvai-v2-diagnostics`

## Status and question

This is a pre-specified **post-hoc diagnostic**, not new prospective validation. The eight molecules, all raw trajectories, the aggregate replay results and the same-curve oracle results are development-exposed. No new simulation, model fitting, endpoint scoring or Tier-B access is permitted.

The sole question is whether the apparent non-deployable-oracle headroom survives when lambda locations are selected with one non-overlapping part of each trajectory and evaluated against response values and a dense integral built only from held-out frames.

## Frozen inputs

Primary molecules are the eight originally prospective dense sentinels:

`Cyclohexane, Octane, Octanol, AceticAcid, Acetamide, DiMethylEther, EthylAcetate, Anthracene`.

Each molecule has the registered 15 lambda windows and 51 stored production frames per window. The four calibration molecules are not included in the primary performance estimate because they selected the frozen Gaussian hyperparameters.

| Input | SHA-256 |
|---|---|
| `active_solvai/results/phase2/dense_responses_prospective.parquet` | `4ce191a2ebfabb10dfcf5e11f98fef91b0bd5b2993f066ab47cfa9cc261c1097` |
| `active_solvai/release/DENSE_SENTINEL_CALIBRATION_LOCK.json` | `3d1f411c5320bf17b0bf83d7f6280645848136c88fed7d7be1a54e9e0fd4975a` |
| `active_solvai/configs/dense_sentinel_v1.json` | `fc6afd2a4f8a255a641228df624769e6ed437dc4cc69dc244480a33c2936594d` |
| `active_solvai/results/phase1/phase1_response_predictions.parquet` | `9d9a514e2cabec3ea6f43600d9987d886ca11f2c96e2cc4f5f0b18f79397c91b` |
| `active_solvai/simulations/dense_pimd2/manifest.csv` | `9c46ad5b78af027f3c82d84ff300c1ea839fbd13650eba704bbfef95f5871358` |
| Sorted prospective `(molecule_id, lambda, raw-energy SHA-256)` records (120 rows) | `e0585cc8e1f0d19754b43eec6eaa54c760fcbf106fe76315e46fbf3b7374451e` |

The implementation must assert these hashes and every raw-energy hash before analysis.

## Frame blocking and all assignments

For each 51-frame window, frames are divided in temporal order into four contiguous blocks with `numpy.array_split`, producing fixed block sizes 13, 13, 13 and 12. No equilibration frame, lambda, molecule or block may be dropped.

All six ordered choices of two training blocks from four are evaluated; the held-out blocks are the exact complement:

| Split | Selection blocks | Evaluation blocks |
|---:|---|---|
| 0 | 0,1 | 2,3 |
| 1 | 0,2 | 1,3 |
| 2 | 0,3 | 1,2 |
| 3 | 1,2 | 0,3 |
| 4 | 1,3 | 0,2 |
| 5 | 2,3 | 0,1 |

Thus every split uses non-overlapping frames, every block appears equally often on each side, and splits 0/5, 1/4 and 2/3 are exact reversals. Some block means can remain temporally correlated because these are divisions of one short trajectory, not independent replicas; the report must state this limitation.

For each window and side, the response is the mean of all frames in its two assigned blocks. Its standard error is the sample standard deviation of the two constituent block means divided by square root two. No result from one side is used to estimate the other side's response or dense integral.

## Frozen schedules and estimators

The lambda grid, initial indices `(2, 6, 12)` corresponding to 0.1, 0.5 and 0.9, Gaussian kernels, hyperparameters, structure priors, generic prior, expected calibration SEMs, sign convention and trapezoidal weights remain exactly those in the immutable campaign.

The following are evaluated at total budgets 5 and 7:

1. **Cross-fitted oracle BQ:** construct its greedy oracle order using only the selection-side response curve, selection-side response errors and selection-side dense integral. Apply the selected lambda locations to evaluation-side response values; condition the frozen SolvAI-prior Gaussian with evaluation-side errors; score only against the evaluation-side dense integral.
2. **Active SolvAI BQ:** use the original frozen molecule-conditioned variance-reduction order. It never uses either selection- or evaluation-side unobserved responses.
3. **Uniform SolvAI BQ:** use the frozen maximin order and SolvAI-prior posterior.
4. **Fixed SolvAI BQ:** use the frozen population order and SolvAI-prior posterior.
5. **Generic BQ:** use the frozen generic calibration prior, covariance and acquisition order.
6. **Uniform direct** and **fixed direct** are retained secondary conventional controls because uniform direct was the strongest original comparator.

At a given split and budget, every method receives evaluation-side values from exactly the same number of selected windows. The cross-fitted oracle alone may use selection-side values to choose locations; it never uses selection-side values in its final integral estimate or target.

No experimental hydration label is read anywhere in this diagnostic.

## Frozen outcomes and aggregation

For every molecule, split, method and budget, save selected lambdas, evaluation-side dense target, estimate, signed error, absolute error, posterior standard deviation and interval inclusion where defined.

Primary performance aggregation is:

1. average absolute error over all six splits within each molecule;
2. average the eight resulting molecule means.

This prevents 48 split-molecule rows from being treated as independent molecules. Paired differences use the same molecule-level aggregation. Confidence intervals use 100,000 molecule-clustered bootstrap resamples with seed `20260903`; report 90% intervals for continuity with the frozen replay and 95% intervals secondarily.

Report:

- MAE at five and seven windows for every method;
- molecule-level mean errors and paired effects;
- evaluation-side dense-integral reproducibility across the three unique complementary partitions, including MAE, RMSE, median absolute difference and 90%/95% molecule-clustered intervals;
- cross-fitted oracle lambda frequencies;
- Jaccard similarity of **added** oracle lambdas (the common inherited three are excluded) for the three exact reversal pairs, plus all-pairs Jaccard as a secondary statistic;
- fraction of molecules improved by the cross-fitted oracle over uniform direct;
- original same-curve headroom `uniform_direct MAE - oracle MAE` and cross-fitted headroom under the same definition;
- headroom survival ratio `cross-fitted headroom / original same-curve headroom`. A negative ratio means the cross-fitted oracle is worse than uniform direct.

No family, molecule, block split, lambda subset or integration convention may be selected or omitted after results are visible.

## Frozen conclusion rule

Let `delta_b = cross-fitted oracle MAE - uniform-direct MAE` at budget `b`; negative values favor the oracle. Let `S_b` be the headroom survival ratio. Let `J_7` be mean reversal-pair Jaccard similarity of added oracle lambdas at seven windows.

Return exactly one conclusion, applied hierarchically:

### B. STABLE BUT UNLEARNED HEADROOM

Return B only if, at both five and seven windows:

- `delta_b <= -0.10 kcal mol^-1`;
- the molecule-clustered 90% upper confidence bound is below zero;
- at least six of eight molecules have lower split-mean oracle error than uniform direct; and
- `S_b >= 0.50`;

and additionally `J_7 >= 0.40`.

### A. NO-GO CONFIRMED

If B fails, return A only if, at both budgets:

- `delta_b >= -0.10 kcal mol^-1`;
- the molecule-clustered 90% lower confidence bound is above `-0.10 kcal mol^-1`; and
- `S_b <= 0.25`.

These conditions rule out a material, stable oracle advantage at the same 0.10-kcal mol^-1 scale used by the original protocol.

### C. INCONCLUSIVE DUE TO RESPONSE NOISE

Return C if neither A nor B is established. Also return C, unless A or B meets all of its strong effect/interval conditions, if the complementary-half dense-integral reproducibility MAE exceeds 0.50 kcal mol^-1 or its median absolute difference exceeds 0.30 kcal mol^-1. These thresholds are 2.5 times and 1.5 times the original 0.20-kcal mol^-1 reconstruction target and are fixed before calculation.

If C is returned, the report must specify the smallest replicated or longer-trajectory experiment capable of separating placement instability from finite-trajectory noise, but no such simulation may be launched.

## What this diagnostic can and cannot support

It can determine whether the striking same-curve oracle values survive a non-overlapping held-out-frame test and whether the original headroom statement should be retained, narrowed or withdrawn.

It cannot become prospective molecular validation, reopen the failed experimental-endpoint direction, validate PIMD8 accuracy, justify multi-fidelity escalation, authorize Tier-B, or overturn the immutable no-go for the deployed policy.
