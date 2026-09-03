"""Time-series diagnostics for fixed-grid trajectory-effort allocation."""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import ks_2samp, rankdata


def _finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def batch_variance_rate(values: np.ndarray, dt_ps: float, batches: int) -> float:
    """Estimate the asymptotic variance rate from non-overlapping batch means."""
    x = _finite(values)
    if batches < 2 or len(x) < 2 * batches:
        return math.nan
    chunks = [chunk for chunk in np.array_split(x, batches) if len(chunk)]
    means = np.asarray([chunk.mean() for chunk in chunks], dtype=float)
    mean_batch_length = float(np.mean([len(chunk) for chunk in chunks]))
    return float(means.var(ddof=1) * mean_batch_length * dt_ps)


def overlapping_batch_variance_rate(values: np.ndarray, dt_ps: float) -> float:
    """Estimate variance rate with deterministic sqrt(n)-length overlapping batches."""
    x = _finite(values)
    n = len(x)
    batch = max(2, int(np.floor(np.sqrt(n))))
    if n < 2 * batch:
        return math.nan
    means = np.convolve(x, np.ones(batch) / batch, mode="valid")
    return float(batch * dt_ps * np.sum((means - x.mean()) ** 2) / len(means))


def newey_west_variance_rate(values: np.ndarray, dt_ps: float) -> float:
    """Bartlett/Newey-West long-run variance rate with a fixed length bandwidth."""
    x = _finite(values)
    n = len(x)
    if n < 4:
        return math.nan
    centered = x - x.mean()
    bandwidth = max(1, int(np.floor(4 * (n / 100) ** (2 / 9))))
    gamma0 = float(centered @ centered / n)
    long_run = gamma0
    for lag in range(1, min(bandwidth, n - 1) + 1):
        gamma = float(centered[lag:] @ centered[:-lag] / n)
        long_run += 2 * (1 - lag / (bandwidth + 1)) * gamma
    return float(max(long_run, 0.0) * dt_ps)


def initial_positive_iat(values: np.ndarray) -> float:
    """Initial-positive-sequence estimate of integrated autocorrelation time."""
    x = _finite(values)
    n = len(x)
    if n < 4 or np.isclose(x.var(), 0.0):
        return math.nan
    centered = x - x.mean()
    gamma0 = float(centered @ centered / n)
    correlations: list[float] = []
    for lag in range(1, n):
        gamma = float(centered[lag:] @ centered[:-lag] / n)
        correlations.append(gamma / gamma0)
    total = 1.0
    for index in range(0, len(correlations) - 1, 2):
        pair = correlations[index] + correlations[index + 1]
        if pair <= 0:
            break
        total += 2 * pair
    return float(max(total, 1.0))


def prefix_diagnostics(values: np.ndarray, times_ps: np.ndarray) -> dict[str, float | bool]:
    """Compute only diagnostics available in a revealed trajectory prefix."""
    x = _finite(values)
    t = np.asarray(times_ps, dtype=float)[: len(x)]
    n = len(x)
    if n < 2:
        raise ValueError("a prefix requires at least two finite observations")
    dt = float(np.median(np.diff(t)))
    split = n // 2
    first, second = x[:split], x[split:]
    slope = float(np.polyfit(t, x, 1)[0]) if n >= 3 else math.nan
    rank_correlation = float(np.corrcoef(rankdata(t), rankdata(x))[0, 1]) if n >= 3 else math.nan
    variance = float(x.var(ddof=1))
    naive_sem = float(np.sqrt(variance / n))
    iat = initial_positive_iat(x)
    ess = float(n / iat) if np.isfinite(iat) and iat > 0 else math.nan
    first_sem = float(first.std(ddof=1) / np.sqrt(len(first))) if len(first) > 1 else math.nan
    second_sem = float(second.std(ddof=1) / np.sqrt(len(second))) if len(second) > 1 else math.nan
    pooled = np.sqrt(first_sem**2 + second_sem**2)
    half_difference = float(second.mean() - first.mean())
    lag1 = float(np.corrcoef(x[:-1], x[1:])[0, 1]) if n >= 3 and variance > 0 else math.nan
    return {
        "prefix_frames": n,
        "prefix_mean": float(x.mean()),
        "prefix_variance": variance,
        "naive_sem": naive_sem,
        "batch_variance_rate_2": batch_variance_rate(x, dt, 2),
        "batch_variance_rate_5": batch_variance_rate(x, dt, 5),
        "overlap_batch_variance_rate": overlapping_batch_variance_rate(x, dt),
        "newey_west_variance_rate": newey_west_variance_rate(x, dt),
        "lag1_autocorrelation": lag1,
        "iat_initial_positive": iat,
        "effective_sample_size": ess,
        "half_mean_difference": half_difference,
        "half_variance_ratio_log": float(
            np.log((second.var(ddof=1) + 1e-8) / (first.var(ddof=1) + 1e-8))
        )
        if len(first) > 1 and len(second) > 1
        else math.nan,
        "half_ks_distance": float(ks_2samp(first, second).statistic),
        "linear_drift_per_ps": slope,
        "rank_drift": rank_correlation,
        "unresolved_equilibration": bool(np.isfinite(pooled) and abs(half_difference) > 2 * pooled),
    }


def complementary_log_difficulty(
    block_a: np.ndarray, block_b: np.ndarray, stabilizer: float = 0.25
) -> float:
    """Return the frozen realized complementary-block log-difficulty proxy."""
    difference = float(np.mean(_finite(block_a)) - np.mean(_finite(block_b)))
    return float(np.log(0.5 * difference**2 + stabilizer**2))
