from __future__ import annotations

import numpy as np

from active_solvai.dense import (
    condition_curve,
    curvature_order,
    fixed_order,
    interpolate_three_point_prior,
    maximin_order,
    rbf_covariance,
    trapezoid_weights,
    variance_reduction_order,
)

GRID = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
INITIAL = (2, 6, 12)


def test_trapezoid_weights_integrate_linear_curve() -> None:
    weights = trapezoid_weights(GRID)
    assert np.isclose(weights.sum(), 1.0)
    assert np.isclose(weights @ (2.0 + 3.0 * GRID), 3.5)


def test_conditioning_interpolates_nearly_noise_free_observations() -> None:
    mean = interpolate_three_point_prior(np.array([1.0, 0.0, -1.0]), GRID)
    covariance = rbf_covariance(GRID, amplitude=4.0, lengthscale=0.2)
    observed = np.array(INITIAL)
    values = np.array([2.0, -1.0, 0.5])
    posterior = condition_curve(mean, covariance, observed, values, np.full(3, 1e-12))
    assert np.allclose(posterior.mean[observed], values, atol=1e-8)
    assert np.all(np.diag(posterior.covariance)[observed] < 1e-7)


def test_every_schedule_visits_each_missing_window_once() -> None:
    prior = interpolate_three_point_prior(np.array([2.0, 0.0, -2.0]), GRID)
    covariance = rbf_covariance(GRID, amplitude=3.0, lengthscale=0.2)
    noise = np.full(len(GRID), 0.1)
    schedules = [
        fixed_order(GRID),
        maximin_order(GRID, INITIAL),
        curvature_order(prior, GRID, INITIAL),
        variance_reduction_order(covariance, trapezoid_weights(GRID), INITIAL, noise),
    ]
    expected = set(range(len(GRID))) - set(INITIAL)
    for schedule in schedules:
        assert len(schedule) == len(expected)
        assert set(schedule) == expected
