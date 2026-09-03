import numpy as np
import pandas as pd

from active_solvai_v3.replay import (
    ChunkLibrary,
    EqualTimePolicy,
    MaximumSemPolicy,
    RandomChunkPolicy,
    run_replay,
)


def synthetic_library() -> ChunkLibrary:
    rng = np.random.default_rng(42)
    chunks = rng.normal(size=(3, 4, 20))
    chunks[2] *= 8.0
    return ChunkLibrary(
        lambdas=np.array([0.0, 0.5, 1.0]),
        chunks=chunks,
        chunk_costs=np.ones((3, 4)),
    )


def test_equal_time_is_balanced_and_deterministic() -> None:
    library = synthetic_library()
    left = run_replay(library, EqualTimePolicy(), initial_chunks=1, total_additional_chunks=6)
    right = run_replay(library, EqualTimePolicy(), initial_chunks=1, total_additional_chunks=6)
    pd.testing.assert_frame_equal(left, right)
    assert left.chosen_window.tolist() == [0, 1, 2, 0, 1, 2]


def test_random_policy_replays_from_seed() -> None:
    library = synthetic_library()
    left = run_replay(library, RandomChunkPolicy(9), initial_chunks=1, total_additional_chunks=5)
    right = run_replay(library, RandomChunkPolicy(9), initial_chunks=1, total_additional_chunks=5)
    pd.testing.assert_frame_equal(left, right)


def test_maximum_sem_uses_only_visible_high_variance_window() -> None:
    replay = run_replay(
        synthetic_library(), MaximumSemPolicy(), initial_chunks=1, total_additional_chunks=1
    )
    assert replay.iloc[0].chosen_window == 2


def test_invalid_policy_cannot_read_or_choose_future() -> None:
    class BadPolicy:
        name = "bad"

        def choose(self, state, eligible):
            assert not hasattr(state, "chunks")
            return int(max(eligible) + 1)

    try:
        run_replay(synthetic_library(), BadPolicy(), initial_chunks=1, total_additional_chunks=1)
    except ValueError as error:
        assert "ineligible" in str(error)
    else:  # pragma: no cover
        raise AssertionError("future-data-invalid policy was accepted")
