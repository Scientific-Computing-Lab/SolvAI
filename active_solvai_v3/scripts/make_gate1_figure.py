#!/usr/bin/env python3
"""Create the deterministic Gate-1 diagnostic figure from frozen outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "active_solvai_v3/results/gate1"
FIGURES = ROOT / "active_solvai_v3/figures"


def main() -> None:
    metrics = pd.read_csv(RESULTS / "gate1_model_metrics.csv")
    predictions = pd.read_parquet(RESULTS / "gate1_oof_predictions.parquet")
    metrics = metrics[np.isclose(metrics.stabilizer, 0.25)]
    primary = predictions[
        np.isclose(predictions.prefix_ps, 1.0) & np.isclose(predictions.stabilizer, 0.25)
    ]

    palette = {
        "lambda_protocol": "#666666",
        "generic_observed": "#E69F00",
        "solvai_conditioned": "#0072B2",
    }
    labels = {
        "lambda_protocol": "lambda/protocol",
        "generic_observed": "generic diagnostics",
        "solvai_conditioned": "SolvAI-conditioned",
    }
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.hashsalt": "active-solvai-v3-gate1",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), constrained_layout=True)

    for model, color in palette.items():
        part = metrics[metrics.model == model].sort_values("prefix_ps")
        axes[0].plot(
            part.prefix_ps,
            part.mae,
            marker="o",
            color=color,
            label=labels[model],
        )
    axes[0].set(
        xlabel="Revealed prefix (ps)",
        ylabel="Log-difficulty MAE",
        xticks=[0.5, 1, 2, 3],
    )
    axes[0].legend(frameon=False, fontsize=7)

    errors = {}
    for model in ("generic_observed", "solvai_conditioned"):
        part = primary[primary.model == model].copy()
        part["absolute_error"] = abs(part.prediction - part.target)
        errors[model] = part.groupby(["molecule_id", "molecule_name"]).absolute_error.mean()
    differences = (errors["solvai_conditioned"] - errors["generic_observed"]).sort_values()
    colors = np.where(differences < 0, "#009E73", "#D55E00")
    axes[1].barh(range(len(differences)), differences, color=colors, height=0.72)
    axes[1].axvline(0, color="#333333", linewidth=0.8)
    axes[1].set_yticks(range(len(differences)))
    axes[1].set_yticklabels([name for _, name in differences.index], fontsize=6)
    axes[1].set_xlabel("SolvAI - generic MAE")

    aligned = float(metrics.loc[metrics.model == "solvai_conditioned", "mae"].iloc[0])
    shuffled = metrics[metrics.model.str.startswith("solvai_conditioned_shuffled")].mae
    axes[2].scatter(np.zeros(len(shuffled)), shuffled, color="#999999", s=25, label="shuffles")
    axes[2].scatter([1], [aligned], color="#0072B2", s=35, zorder=3, label="aligned")
    axes[2].set_xticks([0, 1], ["5 shuffled\ncontrols", "aligned"])
    axes[2].set_ylabel("Log-difficulty MAE")
    axes[2].set_xlim(-0.5, 1.5)
    axes[2].legend(frameon=False, fontsize=7)

    for label, axis in zip("abc", axes, strict=True):
        axis.text(-0.17, 1.05, label, transform=axis.transAxes, fontweight="bold", fontsize=10)
    FIGURES.mkdir(parents=True, exist_ok=True)
    for extension, metadata in (
        ("pdf", {"CreationDate": None, "ModDate": None}),
        ("svg", {"Date": None}),
        ("png", {"Software": "Active SolvAI v3 deterministic figure build"}),
    ):
        fig.savefig(
            FIGURES / f"gate1_identifiability.{extension}",
            dpi=300,
            metadata=metadata,
        )
    svg_path = FIGURES / "gate1_identifiability.svg"
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n"
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
