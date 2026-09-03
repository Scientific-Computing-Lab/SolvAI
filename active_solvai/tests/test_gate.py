import numpy as np

from active_solvai.gate import (
    apply_affine,
    condition_gaussian,
    fit_affine,
    integration_proxy,
)


def test_integration_proxy_constant_curve():
    response = np.array([[2.0, 2.0, 2.0]])
    assert np.allclose(integration_proxy(response), [-2.0])


def test_gaussian_conditioning_moves_observed_coordinate():
    mean = np.zeros(3)
    covariance = np.eye(3)
    posterior_mean, posterior_covariance = condition_gaussian(
        mean, covariance, [1], np.array([2.0]), np.array([0.0])
    )
    assert np.isclose(posterior_mean[1], 2.0, atol=1e-6)
    assert posterior_covariance[1, 1] < 1e-6


def test_affine_roundtrip():
    x = np.arange(5, dtype=float)
    y = 3.0 + 2.0 * x
    assert np.allclose(apply_affine(x, fit_affine(x, y)), y)
