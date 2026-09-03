from __future__ import annotations

import numpy as np
import pandas as pd

from active_solvai.probes import contiguous_block_sem, prefix_summary


def test_prefix_summary_is_sequential():
    frame = pd.DataFrame({"Time": np.arange(10), "dHdL": np.arange(10, dtype=float)})
    rows = prefix_summary(frame, (0.2, 1.0))
    assert rows[0]["frames"] == 2
    assert rows[0]["mean_kcal_mol"] == 0.5
    assert rows[1]["mean_kcal_mol"] == 4.5


def test_contiguous_block_sem_constant_is_zero():
    sem, blocks = contiguous_block_sem(np.ones(20), blocks=5)
    assert sem == 0.0
    assert blocks == 5
