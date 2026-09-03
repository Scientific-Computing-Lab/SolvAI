"""Deterministic, future-data-safe sequential replay over immutable chunks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class RevealedState:
    """The complete and only policy-visible state at one replay decision."""

    lambdas: np.ndarray
    counts: np.ndarray
    means: np.ndarray
    sems: np.ndarray
    cumulative_cost: float
    step: int


class AllocationPolicy(Protocol):
    name: str

    def choose(self, state: RevealedState, eligible: np.ndarray) -> int: ...


@dataclass(frozen=True)
class ChunkLibrary:
    """An immutable rectangular molecule stream: window x chunk x observation."""

    lambdas: np.ndarray
    chunks: np.ndarray
    chunk_costs: np.ndarray

    def __post_init__(self) -> None:
        if self.chunks.ndim != 3:
            raise ValueError("chunks must have shape (windows, chunks, observations)")
        if self.chunks.shape[0] != len(self.lambdas):
            raise ValueError("lambda and window dimensions differ")
        if self.chunk_costs.shape != self.chunks.shape[:2]:
            raise ValueError("chunk_costs must have shape (windows, chunks)")
        if np.any(~np.isfinite(self.chunks)) or np.any(self.chunk_costs < 0):
            raise ValueError("chunks and costs must be finite; costs must be nonnegative")


class EqualTimePolicy:
    name = "equal_time"

    def choose(self, state: RevealedState, eligible: np.ndarray) -> int:
        return int(min(eligible, key=lambda index: (state.counts[index], index)))


class RandomChunkPolicy:
    name = "random_chunk"

    def __init__(self, seed: int) -> None:
        self._rng = np.random.default_rng(seed)

    def choose(self, state: RevealedState, eligible: np.ndarray) -> int:
        del state
        return int(self._rng.choice(np.asarray(eligible, dtype=int)))


class MaximumSemPolicy:
    name = "maximum_sem"

    def choose(self, state: RevealedState, eligible: np.ndarray) -> int:
        return int(max(eligible, key=lambda index: (state.sems[index], -index)))


def _state(library: ChunkLibrary, counts: np.ndarray, step: int) -> RevealedState:
    means = np.empty(len(library.lambdas), dtype=float)
    sems = np.empty(len(library.lambdas), dtype=float)
    for index, count in enumerate(counts):
        values = library.chunks[index, :count].reshape(-1)
        means[index] = values.mean()
        sems[index] = values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 else np.inf
    cost = sum(
        float(library.chunk_costs[index, :count].sum()) for index, count in enumerate(counts)
    )
    return RevealedState(
        lambdas=library.lambdas.copy(),
        counts=counts.copy(),
        means=means,
        sems=sems,
        cumulative_cost=cost,
        step=step,
    )


def run_replay(
    library: ChunkLibrary,
    policy: AllocationPolicy,
    *,
    initial_chunks: int,
    total_additional_chunks: int,
) -> pd.DataFrame:
    """Reveal a common floor, then sequentially expose only chosen next chunks."""
    if initial_chunks < 1 or initial_chunks >= library.chunks.shape[1]:
        raise ValueError("initial_chunks must leave at least one additional chunk")
    counts = np.full(len(library.lambdas), initial_chunks, dtype=int)
    records: list[dict[str, object]] = []
    for step in range(total_additional_chunks):
        state = _state(library, counts, step)
        eligible = np.flatnonzero(counts < library.chunks.shape[1])
        if not len(eligible):
            break
        chosen = int(policy.choose(state, eligible.copy()))
        if chosen not in eligible:
            raise ValueError(f"policy selected ineligible window {chosen}")
        records.append(
            {
                "policy": policy.name,
                "step": step,
                "chosen_window": chosen,
                "chosen_lambda": float(library.lambdas[chosen]),
                "count_before": int(counts[chosen]),
                "cumulative_cost_before": state.cumulative_cost,
                "visible_means": state.means.tolist(),
                "visible_sems": state.sems.tolist(),
                "visible_counts": state.counts.tolist(),
            }
        )
        counts[chosen] += 1
    return pd.DataFrame(records)
