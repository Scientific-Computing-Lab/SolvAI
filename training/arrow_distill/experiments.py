"""Leakage-safe cross-validation utilities for the ARROW distillation sprint."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

REGIMES = {
    "random_oof": "fold_random",
    "family_holdout": "fold_family",
    "scaffold_holdout": "fold_scaffold",
}


def mae(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(mean_absolute_error(y_true, y_pred))


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "mae": mae(y_true, y_pred),
        "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
        "r2": float(r2_score(y_true, y_pred)),
    }


def ridge(alpha: float = 100.0):
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), Ridge(alpha=alpha))


def extra_trees(
    seed: int = 11,
    n_estimators: int = 300,
    max_features: float = 0.7,
    min_samples_leaf: int = 2,
):
    return make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        ExtraTreesRegressor(
            n_estimators=n_estimators,
            max_features=max_features,
            min_samples_leaf=min_samples_leaf,
            random_state=seed,
            n_jobs=-1,
        ),
    )


def random_forest(seed: int = 11):
    return make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        RandomForestRegressor(
            n_estimators=300,
            max_features=0.7,
            min_samples_leaf=2,
            random_state=seed,
            n_jobs=-1,
        ),
    )


def xgboost(seed: int = 11):
    return make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        XGBRegressor(
            n_estimators=300,
            max_depth=3,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_lambda=5.0,
            objective="reg:absoluteerror",
            random_state=seed,
            n_jobs=-1,
        ),
    )


def fit_with_optional_weights(model, x, y, sample_weight=None):
    if sample_weight is None:
        return model.fit(x, y)
    final_name = next(reversed(model.named_steps))
    return model.fit(x, y, **{f"{final_name}__sample_weight": sample_weight})


def inner_splits(
    indices: np.ndarray,
    benchmark: pd.DataFrame,
    regime: str,
    n_splits: int = 4,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return absolute-index inner splits that preserve the outer regime."""
    if regime == "random_oof":
        splitter = KFold(n_splits=n_splits, shuffle=True, random_state=1969)
        pairs = splitter.split(indices)
    else:
        group_column = "functional_group_family" if regime == "family_holdout" else "scaffold"
        groups = benchmark.iloc[indices][group_column].astype(str).to_numpy()
        count = min(n_splits, len(np.unique(groups)))
        splitter = GroupKFold(n_splits=count)
        pairs = splitter.split(indices, groups=groups)
    return [(indices[train], indices[valid]) for train, valid in pairs]


def make_prediction_rows(
    benchmark: pd.DataFrame,
    prediction: np.ndarray,
    method: str,
    regime: str,
    requires_pimd: bool,
    uncertainty: np.ndarray | None = None,
    details: list[str] | None = None,
) -> pd.DataFrame:
    result = benchmark[
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "functional_group_family",
            "scaffold",
            REGIMES[regime],
            "delta_g_exp",
            "delta_g_pimd8",
            "delta_g_classical_arrow",
        ]
    ].copy()
    result = result.rename(columns={REGIMES[regime]: "fold", "delta_g_exp": "y_true"})
    result["method"] = method
    result["regime"] = regime
    result["requires_pimd_at_inference"] = requires_pimd
    result["y_pred"] = prediction
    result["residual"] = result.y_true - result.y_pred
    result["absolute_error"] = result.residual.abs()
    result["uncertainty"] = np.nan if uncertainty is None else uncertainty
    result["fit_details"] = "" if details is None else details
    return result


def oof_single_source(
    benchmark: pd.DataFrame,
    x: np.ndarray,
    y: np.ndarray,
    regime: str,
    factory: Callable[[], object],
) -> np.ndarray:
    prediction = np.full(len(benchmark), np.nan)
    folds = benchmark[REGIMES[regime]].to_numpy()
    for fold in sorted(np.unique(folds)):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        model = factory()
        model.fit(x[train], y[train])
        prediction[test] = model.predict(x[test]).reshape(-1)
    return prediction


@dataclass(frozen=True)
class PublicExtraConfig:
    benchmark_weight: float
    max_features: float
    min_samples_leaf: int

    @property
    def label(self) -> str:
        return f"w={self.benchmark_weight:g};mf={self.max_features:g};leaf={self.min_samples_leaf}"


PUBLIC_EXTRA_CONFIGS = (
    PublicExtraConfig(1.0, 0.7, 2),
    PublicExtraConfig(3.0, 0.7, 2),
    PublicExtraConfig(10.0, 0.7, 2),
    PublicExtraConfig(3.0, 1.0, 2),
)


def _fit_public_extra(
    x_public: np.ndarray,
    y_public: np.ndarray,
    x_benchmark: np.ndarray,
    y_benchmark: np.ndarray,
    config: PublicExtraConfig,
    seed: int,
    n_estimators: int,
):
    model = extra_trees(
        seed=seed,
        n_estimators=n_estimators,
        max_features=config.max_features,
        min_samples_leaf=config.min_samples_leaf,
    )
    x_fit = np.vstack([x_public, x_benchmark])
    y_fit = np.concatenate([y_public, y_benchmark])
    weights = np.concatenate(
        [np.ones(len(y_public)), np.full(len(y_benchmark), config.benchmark_weight)]
    )
    return fit_with_optional_weights(model, x_fit, y_fit, weights)


def nested_public_extra(
    benchmark: pd.DataFrame,
    x: np.ndarray,
    y: np.ndarray,
    x_public: np.ndarray,
    y_public: np.ndarray,
    regime: str,
    seeds: tuple[int, ...] = (11, 29, 47),
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Nested selection of public/benchmark weighting, then seed ensemble."""
    predictions = np.full((len(seeds), len(benchmark)), np.nan)
    details = [""] * len(benchmark)
    folds = benchmark[REGIMES[regime]].to_numpy()
    for outer_fold in sorted(np.unique(folds)):
        test = np.flatnonzero(folds == outer_fold)
        train = np.flatnonzero(folds != outer_fold)
        candidates: list[tuple[float, PublicExtraConfig]] = []
        splits = inner_splits(train, benchmark, regime)
        for config in PUBLIC_EXTRA_CONFIGS:
            inner_prediction = np.full(len(benchmark), np.nan)
            for inner_train, inner_valid in splits:
                model = _fit_public_extra(
                    x_public,
                    y_public,
                    x[inner_train],
                    y[inner_train],
                    config,
                    seed=11,
                    n_estimators=120,
                )
                inner_prediction[inner_valid] = model.predict(x[inner_valid])
            valid = np.concatenate([item[1] for item in splits])
            candidates.append((mae(y[valid], inner_prediction[valid]), config))
        _, selected = min(candidates, key=lambda item: item[0])
        for seed_index, seed in enumerate(seeds):
            model = _fit_public_extra(
                x_public,
                y_public,
                x[train],
                y[train],
                selected,
                seed=seed,
                n_estimators=300,
            )
            predictions[seed_index, test] = model.predict(x[test])
        for index in test:
            details[index] = selected.label
    return predictions.mean(axis=0), predictions.std(axis=0), details


def nested_ceiling_residual(
    benchmark: pd.DataFrame,
    x: np.ndarray,
    regime: str,
) -> tuple[np.ndarray, list[str]]:
    """Nested selection of residual estimator and shrinkage."""
    residual = benchmark.experimental_residual.to_numpy()
    pimd = benchmark.delta_g_pimd8.to_numpy()
    truth = benchmark.delta_g_exp.to_numpy()
    folds = benchmark[REGIMES[regime]].to_numpy()
    prediction = np.full(len(benchmark), np.nan)
    details = [""] * len(benchmark)
    factories: dict[str, Callable[[], object]] = {
        "ridge_10": lambda: ridge(10.0),
        "ridge_100": lambda: ridge(100.0),
        "extra": lambda: extra_trees(seed=11, n_estimators=180),
    }
    alphas = (0.25, 0.5, 0.75, 1.0)
    for outer_fold in sorted(np.unique(folds)):
        test = np.flatnonzero(folds == outer_fold)
        train = np.flatnonzero(folds != outer_fold)
        splits = inner_splits(train, benchmark, regime)
        choices: list[tuple[float, str, float]] = []
        for name, factory in factories.items():
            residual_prediction = np.full(len(benchmark), np.nan)
            for inner_train, inner_valid in splits:
                model = factory()
                model.fit(x[inner_train], residual[inner_train])
                residual_prediction[inner_valid] = model.predict(x[inner_valid]).reshape(-1)
            valid = np.concatenate([item[1] for item in splits])
            for alpha in alphas:
                corrected = pimd[valid] + alpha * residual_prediction[valid]
                choices.append((mae(truth[valid], corrected), name, alpha))
        _, selected_name, selected_alpha = min(choices, key=lambda item: item[0])
        model = factories[selected_name]()
        model.fit(x[train], residual[train])
        prediction[test] = pimd[test] + selected_alpha * model.predict(x[test]).reshape(-1)
        for index in test:
            details[index] = f"{selected_name};alpha={selected_alpha:g}"
    return prediction, details


def multitask_pls_oof(
    benchmark: pd.DataFrame,
    x: np.ndarray,
    target_columns: list[str],
    regime: str,
    n_components: int = 8,
) -> np.ndarray:
    """Small-data shared-latent encoder with one linear head per task."""
    y = benchmark[target_columns].to_numpy(dtype=float)
    folds = benchmark[REGIMES[regime]].to_numpy()
    prediction = np.full(len(benchmark), np.nan)
    for fold in sorted(np.unique(folds)):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        imputer = SimpleImputer(strategy="median")
        x_train = imputer.fit_transform(x[train])
        x_test = imputer.transform(x[test])
        model = PLSRegression(
            n_components=min(n_components, len(train) - 1, x_train.shape[1]),
            scale=True,
            max_iter=1000,
        )
        model.fit(x_train, y[train])
        predicted = np.asarray(model.predict(x_test))
        if predicted.ndim == 1:
            prediction[test] = predicted
        else:
            prediction[test] = predicted[:, 0]
    return prediction


def constrained_weights(base: np.ndarray, truth: np.ndarray) -> np.ndarray:
    weights, _ = nnls(base, truth)
    if weights.sum() <= 1e-12:
        return np.full(base.shape[1], 1.0 / base.shape[1])
    return weights / weights.sum()
