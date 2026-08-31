# Prospective freeze: ARROW endpoint-weight sensitivity

Frozen on 31 August 2026 before fitting or evaluating the sensitivity described
below. This is a one-factor robustness analysis requested after the primary
confirmatory campaign. It does not replace or modify the preregistered primary
analysis in `release/CONFIRMATORY_FREEZE.md`.

## Question

Does the held-out advantage of the full 15-response SolvAI endpoint over its matched
structure-only counterpart persist when each available ARROW outer-training label
has sample weight 1 rather than the prospectively frozen primary weight 3?

## Immutable analysis

- Evaluation set: the same 85 neutral water-solvation molecules.
- Partition: the existing fixed five-fold `fold_random` assignment.
- External endpoint pool: the same frozen 1,280 benchmark-disjoint experimental
  hydration labels, each with weight 1.
- ARROW endpoint rows: only outer-training rows, each changed from weight 3 to
  weight 1. The outer-test experimental label remains absent.
- Teachers: the same standardized-equivalence-exclusion teacher predictions used by
  the current confirmatory endpoint. No teacher is refitted.
- Models compared:
  - **A**, 2,048-bit radius-2 Morgan fingerprint plus 217 RDKit descriptors;
  - **F**, the identical structure representation plus the same 15 response priors.
- Endpoint: the same three `ExtraTreesRegressor` pipelines, each with 360 trees,
  `max_features=0.7`, `min_samples_leaf=2`, `criterion="squared_error"`,
  `bootstrap=False`, no maximum depth, minimum split size 2, and seeds 11, 29 and
  47. Predictions are averaged.
- No architecture, feature, fold, teacher, preprocessing or hyperparameter change is
  permitted after results are observed.

## Frozen inputs

| Input | SHA-256 |
|---|---|
| ARROW benchmark | `2b7928f162d094e7ee10d197e66636ba4ae09b0f76d626136b79c0975d3b0310` |
| Expanded public hydration table | `603ed02b6be25d9a3057e321f2c6ea135b012666cfdb8a1b160e37f347951ec4` |
| Public structure features | `f6d9cd37a90bfc0718261f7251c70100be15947f73d7db12c85659ffc05b28e9` |
| ARROW structure features | `39877c3938616445d9996e093f8f37744a9c13d56202db9e294466396df298b4` |
| Abraham predictions | `d3bbcc28c4893e886564a8dda46d651ec83f2e202cd541e5038231ef91c34033` |
| OpenFF predictions | `b4a2f7ba810321d5daa112e019d52b6deae9923261c867fc3fa9a09f0549e841` |
| GBn2 predictions | `278a3071f95acddace4513d0794b7eed5a66b181c88adb523dd794a74a86af95` |
| Standardized-exclusion CombiSolv-QM predictions | `c63e37268fbd6d557cb66abf6d5dbcbec5deeb31a72403f970ad330560592d52` |
| Standardized-exclusion MolSolv predictions | `e3219006795cbb930db14421365b5d34a947029c986a08e8dec53d8058869533` |
| Standardized-exclusion ConfSolv predictions | `a755f999d40bb8ab5db012fd2a598aa19a8162537e395189d5281beff5455825` |
| Frozen training configuration | `349af4d244c5f17c6728ee098d11e732dd033de224723498983f0cb1a1063a37` |

## Outputs and metrics

The implementation command will be:

```bash
uv run python scripts/run_michael_weight1_sensitivity.py
```

It will write molecule-level predictions, metrics, paired comparison and metadata
under `results/michael_30aug_sensitivity/`. The primary metric is MAE in kcal
mol\(^{-1}\). Secondary metrics are RMSE, median absolute error, \(R^2\) and fraction
of molecules improved. The paired estimand is

\[
|e_i^{F}|-|e_i^{A}|,
\]

so negative values favour the response priors. Its 95% interval uses 100,000
molecule-level bootstrap resamples with seed 20260828, matching the confirmatory
convention.

## Interpretation fixed in advance

- **Positive:** the paired 95% interval lies entirely below zero.
- **Neutral:** the interval includes zero.
- **Negative:** the interval lies entirely above zero.

All outcomes will be retained. A positive result supports robustness to equal
endpoint weighting; it does not establish external validation. A neutral or negative
result weakens the weight-robustness claim but does not erase the separately frozen
zero-ARROW-label transport result. No result licenses tuning of the weight.
