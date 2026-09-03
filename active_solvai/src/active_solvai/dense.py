"""Transparent dense-curve posterior and replay utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator


def trapezoid_weights(grid: np.ndarray) -> np.ndarray:
    """Return linear weights whose dot product is the trapezoidal integral."""
    x = np.asarray(grid, dtype=float)
    if x.ndim != 1 or len(x) < 2 or np.any(np.diff(x) <= 0):
        raise ValueError("grid must be one-dimensional and strictly increasing")
    weights = np.zeros_like(x)
    delta = np.diff(x)
    weights[:-1] += delta / 2.0
    weights[1:] += delta / 2.0
    return weights


def interpolate_three_point_prior(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Interpolate the frozen 0.1/0.5/0.9 response prior onto the dense grid.

    PCHIP is used within the observed range. Endpoint extrapolation is linear
    using the adjacent PCHIP derivative, which is fixed before dense data are
    inspected.
    """
    anchors = np.array([0.1, 0.5, 0.9], dtype=float)
    y = np.asarray(values, dtype=float).reshape(3)
    x = np.asarray(grid, dtype=float)
    interpolator = PchipInterpolator(anchors, y, extrapolate=True)
    return np.asarray(interpolator(x), dtype=float)


def rbf_covariance(
    grid: np.ndarray,
    amplitude: float,
    lengthscale: float,
    local_scale: np.ndarray | None = None,
) -> np.ndarray:
    """Squared-exponential covariance, optionally modulated over lambda."""
    x = np.asarray(grid, dtype=float)
    distance = x[:, None] - x[None, :]
    covariance = amplitude**2 * np.exp(-0.5 * (distance / lengthscale) ** 2)
    if local_scale is not None:
        scale = np.asarray(local_scale, dtype=float).reshape(len(x))
        covariance = scale[:, None] * covariance * scale[None, :]
    covariance += np.eye(len(x)) * max(amplitude**2, 1.0) * 1e-10
    return covariance


def curvature_scale(prior_mean: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Return a bounded, molecule-specific variance scale from prior curvature."""
    mean = np.asarray(prior_mean, dtype=float)
    x = np.asarray(grid, dtype=float)
    first = np.gradient(mean, x)
    second = np.abs(np.gradient(first, x))
    robust = float(np.median(second[second > 0])) if np.any(second > 0) else 1.0
    normalized = second / max(robust, 1e-8)
    return np.clip(0.75 + 0.25 * normalized, 0.75, 2.5)


@dataclass(frozen=True)
class Posterior:
    mean: np.ndarray
    covariance: np.ndarray

    def integral(self, weights: np.ndarray) -> tuple[float, float]:
        w = np.asarray(weights, dtype=float)
        variance = max(float(w @ self.covariance @ w), 0.0)
        return float(w @ self.mean), float(np.sqrt(variance))


def condition_curve(
    prior_mean: np.ndarray,
    prior_covariance: np.ndarray,
    observed_indices: np.ndarray,
    observed_values: np.ndarray,
    observed_variance: np.ndarray,
) -> Posterior:
    """Condition a finite-grid Gaussian response curve on noisy observations."""
    mean = np.asarray(prior_mean, dtype=float)
    covariance = np.asarray(prior_covariance, dtype=float)
    indices = np.asarray(observed_indices, dtype=int)
    values = np.asarray(observed_values, dtype=float)
    noise = np.maximum(np.asarray(observed_variance, dtype=float), 0.0)
    if not len(indices):
        return Posterior(mean.copy(), covariance.copy())
    k_oo = covariance[np.ix_(indices, indices)] + np.diag(noise)
    k_all_o = covariance[:, indices]
    solve = np.linalg.solve(k_oo, values - mean[indices])
    posterior_mean = mean + k_all_o @ solve
    posterior_covariance = covariance - k_all_o @ np.linalg.solve(k_oo, k_all_o.T)
    posterior_covariance = 0.5 * (posterior_covariance + posterior_covariance.T)
    return Posterior(posterior_mean, posterior_covariance)


def integral_variance_reduction(
    covariance: np.ndarray,
    observed_indices: np.ndarray,
    candidate_index: int,
    noise_variance: np.ndarray,
    weights: np.ndarray,
) -> float:
    """Expected reduction in integral variance from one additional point."""
    dummy = np.zeros(covariance.shape[0], dtype=float)
    current = condition_curve(
        dummy,
        covariance,
        observed_indices,
        np.zeros(len(observed_indices)),
        noise_variance[observed_indices],
    )
    current_variance = float(weights @ current.covariance @ weights)
    updated_indices = np.append(observed_indices, candidate_index)
    updated_noise = np.append(
        noise_variance[observed_indices], max(noise_variance[candidate_index], 0.0)
    )
    updated = condition_curve(
        dummy,
        covariance,
        updated_indices,
        np.zeros(len(updated_indices)),
        updated_noise,
    )
    updated_variance = float(weights @ updated.covariance @ weights)
    return max(current_variance - updated_variance, 0.0)


def maximin_order(grid: np.ndarray, initial: tuple[int, ...]) -> list[int]:
    """Deterministic uniform refinement order with lower-lambda tie breaking."""
    chosen = list(initial)
    remaining = set(range(len(grid))) - set(chosen)
    order: list[int] = []
    while remaining:
        candidate = max(
            remaining,
            key=lambda index: (min(abs(grid[index] - grid[seen]) for seen in chosen), -index),
        )
        order.append(candidate)
        chosen.append(candidate)
        remaining.remove(candidate)
    return order


def variance_reduction_order(
    covariance: np.ndarray,
    weights: np.ndarray,
    initial: tuple[int, ...],
    noise_variance: np.ndarray,
) -> list[int]:
    """Sequential Bayesian-quadrature acquisition order."""
    chosen = list(initial)
    remaining = set(range(len(weights))) - set(chosen)
    order: list[int] = []
    while remaining:
        candidate = max(
            remaining,
            key=lambda index: (
                integral_variance_reduction(
                    covariance,
                    np.asarray(chosen, dtype=int),
                    index,
                    noise_variance,
                    weights,
                ),
                -index,
            ),
        )
        order.append(candidate)
        chosen.append(candidate)
        remaining.remove(candidate)
    return order


def fixed_order(grid: np.ndarray) -> list[int]:
    """Frozen population schedule after the inherited 0.1/0.5/0.9 points."""
    requested = (0.0, 1.0, 0.3, 0.7, 0.05, 0.95, 0.2, 0.4, 0.6, 0.8, 0.75, 0.85)
    return [int(np.flatnonzero(np.isclose(grid, value))[0]) for value in requested]


def curvature_order(
    prior_mean: np.ndarray, grid: np.ndarray, initial: tuple[int, ...]
) -> list[int]:
    """Greedy non-probabilistic schedule from current interpolant curvature."""
    chosen = list(initial)
    remaining = set(range(len(grid))) - set(chosen)
    order: list[int] = []
    while remaining:
        interpolation = PchipInterpolator(grid[chosen], prior_mean[chosen], extrapolate=True)
        dense = interpolation(grid)
        curvature = np.abs(np.gradient(np.gradient(dense, grid), grid))
        candidate = max(
            remaining,
            key=lambda index: (
                curvature[index] * min(abs(grid[index] - grid[seen]) for seen in chosen),
                -index,
            ),
        )
        order.append(candidate)
        chosen.append(candidate)
        chosen.sort()
        remaining.remove(candidate)
    return order


def observed_pchip(
    grid: np.ndarray, observed_indices: np.ndarray, observed_values: np.ndarray
) -> np.ndarray:
    """Conventional PCHIP reconstruction from only the revealed observations."""
    indices = np.asarray(observed_indices, dtype=int)
    order = np.argsort(indices)
    interpolator = PchipInterpolator(
        np.asarray(grid, dtype=float)[indices[order]],
        np.asarray(observed_values, dtype=float)[order],
        extrapolate=True,
    )
    return np.asarray(interpolator(np.asarray(grid, dtype=float)), dtype=float)
