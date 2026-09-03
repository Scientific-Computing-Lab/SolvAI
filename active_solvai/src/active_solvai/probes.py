"""Readers and summaries for native Arbalest response files."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


def energy_group(path: Path) -> str:
    match = re.search(r"__([^.]*)\.ene$", path.name)
    raw = match.group(1) if match else path.stem
    return re.sub(r"[^A-Za-z0-9]+", "_", raw).strip("_").lower()


def read_energy(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [str(column).strip() for column in frame.columns]
    frame = frame.loc[:, [bool(column) for column in frame.columns]]
    return frame.apply(pd.to_numeric, errors="coerce")


def contiguous_block_sem(values: np.ndarray, blocks: int = 5) -> tuple[float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    chunks = [chunk for chunk in np.array_split(values, blocks) if len(chunk)]
    means = np.array([chunk.mean() for chunk in chunks], dtype=float)
    if len(means) < 2:
        return math.nan, len(means)
    return float(means.std(ddof=1) / np.sqrt(len(means))), len(means)


def prefix_summary(frame: pd.DataFrame, fractions: tuple[float, ...]) -> list[dict[str, float]]:
    """Summarize sequential trajectory prefixes without future-frame selection."""
    rows: list[dict[str, float]] = []
    response_columns = [column for column in frame if column.lower().startswith("dhdl")]
    for fraction in fractions:
        count = max(1, math.ceil(len(frame) * fraction))
        prefix = frame.iloc[:count]
        for column in response_columns:
            values = prefix[column].dropna().to_numpy(float)
            if not len(values):
                continue
            block_sem, block_count = contiguous_block_sem(values)
            rows.append(
                {
                    "trajectory_fraction": float(fraction),
                    "frames": len(values),
                    "component": column,
                    "mean_kcal_mol": float(values.mean()),
                    "sd_kcal_mol": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
                    "naive_sem_kcal_mol": (
                        float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0
                    ),
                    "five_block_sem_kcal_mol": block_sem,
                    "block_count": int(block_count),
                }
            )
    return rows
