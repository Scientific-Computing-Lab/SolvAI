#!/usr/bin/env python3
"""Benchmark the released SolvAI artifact without fitting or model selection."""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import chemprop
import joblib
import lightgbm
import numpy
import pandas
import pandas as pd
import rdkit
import sklearn
import torch

from solv_ai import predict_smiles

ROOT = Path(__file__).resolve().parents[1]


def percentile(values: list[float], percentile_value: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentile_value / 100.0
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def timed_prediction(smiles: list[str]) -> float:
    start = time.perf_counter()
    predict_smiles(smiles)
    return time.perf_counter() - start


def cpu_model() -> str:
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
    return platform.processor() or "unknown"


def gpu_model() -> str:
    try:
        return subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], text=True
        ).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "not used / unavailable"


def measured_peak_rss_kib(smiles: str) -> int | None:
    executable = Path("/usr/bin/time")
    if not executable.is_file():
        return None
    with tempfile.TemporaryDirectory(prefix="solvai_benchmark_") as directory:
        output = Path(directory) / "time.txt"
        subprocess.run(
            [
                str(executable),
                "-v",
                "-o",
                str(output),
                sys.executable,
                "-m",
                "solv_ai.cli",
                smiles,
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        for line in output.read_text().splitlines():
            if "Maximum resident set size" in line:
                return int(line.rsplit(":", 1)[1].strip())
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--warm-runs", type=int, default=5)
    parser.add_argument("--batch-runs", type=int, default=3)
    parser.add_argument("--batch-size", type=int, default=32)
    args = parser.parse_args()
    if args.warm_runs < 3 or args.batch_runs < 3:
        raise ValueError("Use at least three repeated measurements")

    benchmark = pd.read_parquet(ROOT / "data/benchmark/arrow_solvation_master.parquet")
    smiles = benchmark.loc[benchmark.solvent.eq("water"), "canonical_smiles"].tolist()
    batch = smiles[: args.batch_size]

    cold = timed_prediction(["CCO"])
    warm = [timed_prediction(["CCO"]) for _ in range(args.warm_runs)]
    batch_times = [timed_prediction(batch) for _ in range(args.batch_runs)]
    batch_per_molecule = [value / len(batch) for value in batch_times]
    peak_rss = measured_peak_rss_kib("CCO")

    # ARROW Supplementary Table 2: 16–18 h per 1 ns PIMD8 window on one
    # RTX 2080 Ti + two CPU cores; the published hydration protocol uses 15 windows.
    published_total_hours = [15 * 16, 15 * 18]
    median_single_seconds = statistics.median(warm)
    median_batch_per_molecule_seconds = statistics.median(batch_per_molecule)
    result = {
        "schema_version": 1,
        "artifact": "models/final",
        "input": "SMILES",
        "hardware": {
            "cpu": cpu_model(),
            "logical_cpus": os.cpu_count(),
            "gpu_present_but_not_used": gpu_model(),
            "platform": platform.platform(),
        },
        "software": {
            "python": platform.python_version(),
            "chemprop": chemprop.__version__,
            "torch": torch.__version__,
            "rdkit": rdkit.__version__,
            "numpy": numpy.__version__,
            "pandas": pandas.__version__,
            "scikit_learn": sklearn.__version__,
            "lightgbm": lightgbm.__version__,
            "joblib": joblib.__version__,
        },
        "single_molecule": {
            "cold_seconds": cold,
            "warm_runs_seconds": warm,
            "warm_median_seconds": median_single_seconds,
            "warm_p95_seconds": percentile(warm, 95),
            "peak_rss_kib": peak_rss,
        },
        "batch": {
            "size": len(batch),
            "runs_seconds": batch_times,
            "median_seconds": statistics.median(batch_times),
            "p95_seconds": percentile(batch_times, 95),
            "median_seconds_per_molecule": median_batch_per_molecule_seconds,
            "molecules_per_second": 1.0 / median_batch_per_molecule_seconds,
        },
        "published_arrow_pimd8_reference": {
            "hardware": "NVIDIA GeForce RTX 2080 Ti + 2 CPU cores",
            "hours_per_1ns_window_range": [16, 18],
            "windows": 15,
            "total_serial_gpu_hours_range": published_total_hours,
            "ideal_parallel_wall_hours_range_on_15_gpus": [16, 18],
            "source": "Pereyaslavets et al., Nature Communications 13, 414 (2022), Supplementary Table 2",
            "doi": "10.1038/s41467-022-28041-0",
        },
        "derived_comparison": {
            "single_warm_total_work_speedup_range": [
                published_total_hours[0] * 3600 / median_single_seconds,
                published_total_hours[1] * 3600 / median_single_seconds,
            ],
            "batch_amortized_total_work_speedup_range": [
                published_total_hours[0] * 3600 / median_batch_per_molecule_seconds,
                published_total_hours[1] * 3600 / median_batch_per_molecule_seconds,
            ],
            "caveat": (
                "The numerator is a published RTX 2080 Ti simulation estimate and the "
                "denominator is measured CPU inference on the release host; this compares "
                "total per-molecule work, not same-hardware kernel throughput."
            ),
        },
    }
    output = ROOT / "results/runtime/runtime_benchmark.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
