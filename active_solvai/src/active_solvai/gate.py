"""Fold-safe helpers for the Active SolvAI actual-observation gate."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from sklearn.covariance import LedoitWolf
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

MODEL_SEEDS = (11, 29, 47)
RIDGE_ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
LAMBDA_VALUES = np.array([0.1, 0.5, 0.9], dtype=float)


def response_model(seed: int) -> object:
    """Return the preregistered multi-output response estimator."""
    return make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        ExtraTreesRegressor(
            n_estimators=360,
            max_features=0.7,
            min_samples_leaf=3,
            min_samples_split=2,
            max_depth=None,
            bootstrap=False,
            random_state=seed,
            n_jobs=-1,
        ),
    )


def ensemble_fit_predict(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
) -> np.ndarray:
    """Fit the frozen response ensemble and average member predictions."""
    members = []
    for seed in MODEL_SEEDS:
        model = response_model(seed)
        model.fit(x_train, y_train)
        members.append(np.asarray(model.predict(x_test), dtype=float))
    return np.mean(members, axis=0)


def integration_proxy(response: np.ndarray) -> np.ndarray:
    """Integrate three annihilation derivatives using the frozen endpoint rule."""
    values = np.asarray(response, dtype=float)
    if values.shape[-1] != 3:
        raise ValueError("The three-point integration proxy requires three lambda values")
    integral = (
        0.1 * values[..., 0]
        + (0.8 / 6.0) * (values[..., 0] + 4.0 * values[..., 1] + values[..., 2])
        + 0.1 * values[..., 2]
    )
    return -integral


@dataclass(frozen=True)
class GaussianCurve:
    mean: np.ndarray
    covariance: np.ndarray


def fit_gaussian_curve(samples: np.ndarray) -> GaussianCurve:
    """Fit the fixed shrinkage Gaussian model used by the sparse posterior."""
    values = np.asarray(samples, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("Curve samples must have shape (n, 3)")
    fit = LedoitWolf().fit(values)
    return GaussianCurve(mean=values.mean(axis=0), covariance=fit.covariance_)


def condition_gaussian(
    prior_mean: np.ndarray,
    covariance: np.ndarray,
    observed_indices: Iterable[int],
    observed_values: np.ndarray,
    observed_noise_variance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Condition one three-point Gaussian curve on noisy observations."""
    mean = np.asarray(prior_mean, dtype=float).reshape(3)
    cov = np.asarray(covariance, dtype=float).reshape(3, 3)
    indices = np.asarray(tuple(observed_indices), dtype=int)
    values = np.asarray(observed_values, dtype=float).reshape(len(indices))
    noise = np.asarray(observed_noise_variance, dtype=float).reshape(len(indices))
    if not len(indices):
        return mean.copy(), cov.copy()
    scale = max(float(np.trace(cov) / 3.0), 1.0)
    cov = cov + np.eye(3) * (1e-8 * scale)
    k_oo = cov[np.ix_(indices, indices)] + np.diag(np.maximum(noise, 0.0))
    k_all_o = cov[:, indices]
    solve = np.linalg.solve(k_oo, values - mean[indices])
    posterior_mean = mean + k_all_o @ solve
    posterior_cov = cov - k_all_o @ np.linalg.solve(k_oo, k_all_o.T)
    posterior_cov = 0.5 * (posterior_cov + posterior_cov.T)
    return posterior_mean, posterior_cov


def fit_affine(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Fit y = intercept + slope*x by least squares."""
    design = np.column_stack([np.ones(len(x)), np.asarray(x, dtype=float)])
    coefficients, *_ = np.linalg.lstsq(design, np.asarray(y, dtype=float), rcond=None)
    return float(coefficients[0]), float(coefficients[1])


def apply_affine(x: np.ndarray, coefficients: tuple[float, float]) -> np.ndarray:
    return coefficients[0] + coefficients[1] * np.asarray(x, dtype=float)


def ridge_correct(
    train_features: np.ndarray,
    train_residual: np.ndarray,
    test_features: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """Fit and apply the frozen standardized ridge correction."""
    model = make_pipeline(StandardScaler(), Ridge(alpha=alpha, fit_intercept=True))
    model.fit(np.asarray(train_features, dtype=float), np.asarray(train_residual, dtype=float))
    return np.asarray(model.predict(np.asarray(test_features, dtype=float)), dtype=float)


def choose_ridge_alpha(
    features: np.ndarray,
    residual: np.ndarray,
    baseline: np.ndarray,
    truth: np.ndarray,
    folds: np.ndarray,
) -> float:
    """Select alpha using only the supplied outer-training cross-fit rows."""
    scores: list[tuple[float, float]] = []
    for alpha in RIDGE_ALPHAS:
        predictions = np.full(len(truth), np.nan, dtype=float)
        for fold in sorted(np.unique(folds)):
            train = folds != fold
            test = folds == fold
            predictions[test] = baseline[test] + ridge_correct(
                features[train], residual[train], features[test], alpha
            )
        mae = float(np.mean(np.abs(np.asarray(truth) - predictions)))
        scores.append((mae, -float(alpha)))
    # The negative alpha key makes the larger alpha win an exact MAE tie.
    return -min(scores)[1]
