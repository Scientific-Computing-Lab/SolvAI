#!/usr/bin/env python3
"""Run the frozen v3 simulation-based power and reference-reliability gate."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai_v3/results/power_reliability"
FIGURES = ROOT / "active_solvai_v3/figures"
TRIALS = 50_000
SEED = 20260903
THROUGHPUT_PS_GPU_H = 369.902
MOLECULE_COUNTS = (12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 256)
DURATIONS_PS = (20, 50, 100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def effect_patterns(path: Path) -> tuple[dict[int, np.ndarray], dict[int, float], pd.DataFrame]:
    frame = pd.read_parquet(path)
    keep = frame[
        frame.method.isin(["active_solvai_bq", "generic_bq"])
        & frame.total_windows.isin([5, 7])
    ]
    pivot = keep.pivot_table(
        index=["molecule_id", "molecule_name", "total_windows"],
        columns="method",
        values="absolute_integral_error_kcal_mol",
        aggfunc="mean",
    ).reset_index()
    pivot["observed_effect"] = pivot.active_solvai_bq - pivot.generic_bq
    patterns: dict[int, np.ndarray] = {}
    target_effects: dict[int, float] = {}
    for budget in (5, 7):
        part = pivot[pivot.total_windows == budget]
        observed = part.observed_effect.to_numpy(float)
        patterns[budget] = observed - observed.mean()
        generic_mae = float(part.generic_bq.mean())
        target_effects[budget] = -0.20 * generic_mae
        pivot.loc[pivot.total_windows == budget, "generic_reference_mae"] = generic_mae
        pivot.loc[pivot.total_windows == budget, "minimum_useful_effect"] = target_effects[
            budget
        ]
    return patterns, target_effects, pivot


def power_table(
    patterns: dict[int, np.ndarray], target_effects: dict[int, float]
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for n_molecules in MOLECULE_COUNTS:
        for budget in (5, 7):
            rng = np.random.default_rng(SEED + n_molecules * 100 + budget)
            indices = rng.integers(0, len(patterns[budget]), size=(TRIALS, n_molecules))
            centered = patterns[budget][indices]
            critical = stats.t.ppf(0.975, n_molecules - 1)
            for alternative, effect in (
                ("null", 0.0),
                ("minimum_20pct", target_effects[budget]),
                ("larger_30pct", 1.5 * target_effects[budget]),
            ):
                observations = centered + effect
                means = observations.mean(axis=1)
                standard_errors = observations.std(axis=1, ddof=1) / math.sqrt(n_molecules)
                upper = means + critical * standard_errors
                rows.append(
                    {
                        "n_molecules": n_molecules,
                        "budget_windows": budget,
                        "alternative": alternative,
                        "true_mean_effect_kcal_mol": effect,
                        "detection_probability": float(np.mean(upper < 0)),
                        "mean_ci_width_kcal_mol": float(
                            np.mean(2 * critical * standard_errors)
                        ),
                    }
                )
    return pd.DataFrame(rows)


def reliability_table(path: Path, threshold: float) -> pd.DataFrame:
    frame = pd.read_csv(path)
    unique = frame[frame.split.isin([0, 1, 2])]
    base = unique.absolute_dense_difference_kcal_mol.to_numpy(float)
    rows: list[dict[str, float | int]] = []
    for n_molecules in MOLECULE_COUNTS:
        rng = np.random.default_rng(SEED + 1_000_000 + n_molecules)
        samples = base[rng.integers(0, len(base), size=(TRIALS, n_molecules))]
        q90_at_2p5 = np.quantile(samples, 0.90, axis=1)
        for duration in DURATIONS_PS:
            scaled = q90_at_2p5 * math.sqrt(2.5 / duration)
            rows.append(
                {
                    "n_molecules": n_molecules,
                    "production_ps_per_stream_window": duration,
                    "reference_threshold_kcal_mol": threshold,
                    "reliability_probability": float(np.mean(scaled <= threshold)),
                    "median_cohort_q90_kcal_mol": float(np.median(scaled)),
                    "q90_of_cohort_q90_kcal_mol": float(np.quantile(scaled, 0.90)),
                }
            )
    return pd.DataFrame(rows)


def build_design_grid(
    power: pd.DataFrame, reliability: pd.DataFrame
) -> pd.DataFrame:
    minimum = power[power.alternative == "minimum_20pct"].pivot(
        index="n_molecules", columns="budget_windows", values="detection_probability"
    )
    null = power[power.alternative == "null"].pivot(
        index="n_molecules", columns="budget_windows", values="detection_probability"
    )
    rows: list[dict[str, float | int | bool]] = []
    for item in reliability.itertuples(index=False):
        n_molecules = int(item.n_molecules)
        duration = int(item.production_ps_per_stream_window)
        production_ps = n_molecules * 15 * 2 * duration
        gpu_hours = production_ps / THROUGHPUT_PS_GPU_H
        rows.append(
            {
                "n_molecules": n_molecules,
                "production_ps_per_stream_window": duration,
                "streams": 2,
                "lambda_windows": 15,
                "power_budget5": float(minimum.loc[n_molecules, 5]),
                "power_budget7": float(minimum.loc[n_molecules, 7]),
                "type1_budget5": float(null.loc[n_molecules, 5]),
                "type1_budget7": float(null.loc[n_molecules, 7]),
                "reference_reliability_probability": float(item.reliability_probability),
                "production_total_ps": production_ps,
                "projected_production_gpu_hours": gpu_hours,
                "operational_reservation_gpu_hours": gpu_hours * 1.11,
                "numerically_adequate": bool(
                    minimum.loc[n_molecules, 5] >= 0.80
                    and minimum.loc[n_molecules, 7] >= 0.80
                    and null.loc[n_molecules, 5] <= 0.075
                    and null.loc[n_molecules, 7] <= 0.075
                    and item.reliability_probability >= 0.80
                ),
            }
        )
    return pd.DataFrame(rows)


def make_figure(power: pd.DataFrame, reliability: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.hashsalt": "active-solvai-v3-power",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    minimum = power[power.alternative == "minimum_20pct"]
    for budget, color in ((5, "#0072B2"), (7, "#009E73")):
        part = minimum[minimum.budget_windows == budget]
        axes[0].plot(
            part.n_molecules,
            part.detection_probability,
            marker="o",
            color=color,
            label=f"{budget}-window effect",
        )
    axes[0].axhline(0.8, color="#555555", linestyle="--", linewidth=0.8)
    axes[0].set(xlabel="Molecules", ylabel="Detection probability", ylim=(0, 1.02))
    axes[0].legend(frameon=False)
    for n_molecules, color in ((20, "#E69F00"), (64, "#0072B2"), (128, "#009E73")):
        part = reliability[reliability.n_molecules == n_molecules]
        axes[1].plot(
            part.production_ps_per_stream_window,
            part.reliability_probability,
            marker="o",
            color=color,
            label=f"n={n_molecules}",
        )
    axes[1].axhline(0.8, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].set_xscale("log")
    axes[1].set(
        xlabel="Production per stream-window (ps)",
        ylabel="Reference-gate probability",
        ylim=(0, 1.02),
    )
    axes[1].legend(frameon=False)
    for extension, metadata in (
        ("pdf", {"CreationDate": None, "ModDate": None}),
        ("svg", {"Date": None}),
        ("png", {"Software": "Active SolvAI v3 deterministic figure build"}),
    ):
        fig.savefig(
            FIGURES / f"v3_power_reliability.{extension}",
            dpi=300,
            metadata=metadata,
        )
    svg_path = FIGURES / "v3_power_reliability.svg"
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n"
    )
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prediction_path = ROOT / "active_solvai/results/phase2/dense_replay_predictions.parquet"
    reproducibility_path = (
        ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/"
        "dense_integral_reproducibility_by_split.csv"
    )
    patterns, target_effects, effect_inputs = effect_patterns(prediction_path)
    power = power_table(patterns, target_effects)
    reference_threshold = 0.5 * min(abs(value) for value in target_effects.values())
    reliability = reliability_table(reproducibility_path, reference_threshold)
    design = build_design_grid(power, reliability)
    adequate = design[design.numerically_adequate].sort_values(
        ["operational_reservation_gpu_hours", "n_molecules"]
    )
    smallest = None if adequate.empty else adequate.iloc[0].to_dict()

    effect_inputs.to_csv(OUT / "inherited_effect_patterns.csv", index=False)
    power.to_csv(OUT / "power_by_molecule_count.csv", index=False)
    reliability.to_csv(OUT / "reference_reliability_grid.csv", index=False)
    design.to_csv(OUT / "prospective_design_grid.csv", index=False)
    make_figure(power, reliability)

    canonical = {
        "schema_version": 1,
        "protocol_commit": "9610813",
        "monte_carlo_trials": TRIALS,
        "seed": SEED,
        "minimum_effects_kcal_mol": {str(key): value for key, value in target_effects.items()},
        "reference_threshold_kcal_mol": reference_threshold,
        "smallest_numerically_adequate_design": smallest,
        "launch_authorized": False,
        "launch_blockers": [
            "Gate-1 SolvAI identifiability criterion failed.",
            "The planning calculation transfers heterogeneity from the failed v1 placement policy, not an observed v3 allocation policy.",
            "Inverse-square-root reference scaling is optimistic and unvalidated at the projected duration.",
        ],
    }
    canonical_path = OUT / "power_reliability_canonical.json"
    canonical_path.write_text(json.dumps(canonical, indent=2) + "\n")
    outputs = [
        OUT / "inherited_effect_patterns.csv",
        OUT / "power_by_molecule_count.csv",
        OUT / "reference_reliability_grid.csv",
        OUT / "prospective_design_grid.csv",
        canonical_path,
        FIGURES / "v3_power_reliability.pdf",
        FIGURES / "v3_power_reliability.svg",
        FIGURES / "v3_power_reliability.png",
    ]
    manifest = {
        "schema_version": 1,
        "inputs": {
            str(prediction_path): sha256(prediction_path),
            str(reproducibility_path): sha256(reproducibility_path),
        },
        "outputs": {str(path): sha256(path) for path in outputs},
    }
    (OUT / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(canonical, indent=2))


if __name__ == "__main__":
    main()
