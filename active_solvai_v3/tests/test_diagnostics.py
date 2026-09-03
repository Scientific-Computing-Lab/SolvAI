import numpy as np

from active_solvai_v3.diagnostics import complementary_log_difficulty, prefix_diagnostics


def test_prefix_diagnostics_are_prefix_local() -> None:
    times = np.arange(0.1, 2.1, 0.1)
    prefix = np.sin(times) + 0.1 * times
    first = prefix_diagnostics(prefix[:10], times[:10])
    changed_future = prefix.copy()
    changed_future[10:] = 1e9
    second = prefix_diagnostics(changed_future[:10], times[:10])
    assert first == second


def test_complementary_difficulty_is_symmetric() -> None:
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([3.0, 4.0, 5.0])
    assert complementary_log_difficulty(a, b) == complementary_log_difficulty(b, a)
