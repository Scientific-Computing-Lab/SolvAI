#!/usr/bin/env python3
"""Regenerate all main and Extended Data figures from frozen release artifacts."""

from __future__ import annotations

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "figures/main"
EXTENDED = ROOT / "figures/extended_data"

BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILION = "#D55E00"
PURPLE = "#CC79A7"
GREY = "#6B7280"
LIGHT = "#F4F6F8"
DARK = "#17212B"


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "solvai-1.0.0",
        }
    )


def save(figure: plt.Figure, directory: Path, stem: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        if suffix == "pdf":
            fixed_date = datetime(2026, 8, 27, tzinfo=UTC)
            metadata = {
                "Creator": "SolvAI 1.0.0",
                "CreationDate": fixed_date,
                "ModDate": fixed_date,
            }
        elif suffix == "svg":
            metadata = {"Creator": "SolvAI 1.0.0", "Date": "2026-08-27"}
        else:
            metadata = {"Software": "SolvAI 1.0.0"}
        figure.savefig(
            directory / f"{stem}.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
            metadata=metadata,
        )
    plt.close(figure)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.08, 1.04, label, transform=axis.transAxes, weight="bold", fontsize=11)


def rounded_box(
    axis: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    *,
    face: str,
    edge: str,
    size: float = 9,
    weight: str = "normal",
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.2,
        facecolor=face,
        edgecolor=edge,
    )
    axis.add_patch(patch)
    axis.text(
        x + width / 2, y + height / 2, text, ha="center", va="center", size=size, weight=weight
    )


def arrow(
    axis: plt.Axes, start: tuple[float, float], end: tuple[float, float], *, dashed: bool = False
) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.4,
            linestyle="--" if dashed else "-",
            color=DARK,
        )
    )


def figure_1(metrics: dict) -> None:
    data = metrics["data_counts"]
    fig, axis = plt.subplots(figsize=(10.8, 4.4))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.add_patch(
        FancyBboxPatch(
            (0.01, 0.08),
            0.31,
            0.84,
            boxstyle="round,pad=0.02",
            facecolor="#FFF7E6",
            edgecolor=ORANGE,
            linewidth=1.5,
        )
    )
    axis.add_patch(
        FancyBboxPatch(
            (0.36, 0.08),
            0.34,
            0.84,
            boxstyle="round,pad=0.02",
            facecolor="#EEF7FB",
            edgecolor=BLUE,
            linewidth=1.5,
        )
    )
    axis.add_patch(
        FancyBboxPatch(
            (0.74, 0.08),
            0.25,
            0.84,
            boxstyle="round,pad=0.02",
            facecolor="#EDF8F4",
            edgecolor=GREEN,
            linewidth=1.5,
        )
    )
    axis.text(0.165, 0.86, "TRAINING-TIME PHYSICS", ha="center", weight="bold", color="#7C4A03")
    rounded_box(
        axis,
        (0.035, 0.61),
        0.125,
        0.16,
        f"SMD(water)\n{data['molsolv_training_structures']:,}",
        face="white",
        edge=ORANGE,
        size=8,
    )
    rounded_box(
        axis,
        (0.185, 0.61),
        0.125,
        0.16,
        f"Conformer\nresponse\n{data['confsolv_training_connectivities']:,}",
        face="white",
        edge=ORANGE,
        size=7.5,
    )
    rounded_box(
        axis,
        (0.035, 0.38),
        0.125,
        0.16,
        "QM / implicit\nresponse",
        face="white",
        edge=GREY,
    )
    rounded_box(
        axis,
        (0.185, 0.38),
        0.125,
        0.16,
        "Alchemical /\ninteraction\nresponse",
        face="white",
        edge=GREY,
        size=7.5,
    )
    axis.text(
        0.165,
        0.20,
        "Physics calculations supply labels once",
        ha="center",
        color=GREY,
        style="italic",
    )

    axis.text(0.53, 0.86, "PHYSICS DISTILLATION", ha="center", weight="bold", color=BLUE)
    rounded_box(
        axis,
        (0.40, 0.61),
        0.26,
        0.15,
        "Structure → response surrogates\n15 predicted physical priors",
        face="white",
        edge=BLUE,
        weight="bold",
    )
    rounded_box(
        axis,
        (0.40, 0.33),
        0.26,
        0.15,
        "Molecular descriptors\n+ response vector",
        face="white",
        edge=SKY,
    )
    rounded_box(
        axis,
        (0.43, 0.14),
        0.20,
        0.10,
        "Hydration endpoint model",
        face="white",
        edge=BLUE,
        weight="bold",
    )
    arrow(axis, (0.53, 0.60), (0.53, 0.49))
    arrow(axis, (0.53, 0.32), (0.53, 0.25))
    arrow(axis, (0.32, 0.58), (0.395, 0.67))

    axis.text(0.865, 0.86, "DEPLOYMENT", ha="center", weight="bold", color=GREEN)
    rounded_box(axis, (0.785, 0.67), 0.16, 0.10, "SMILES", face="white", edge=GREEN, weight="bold")
    rounded_box(
        axis, (0.785, 0.46), 0.16, 0.10, "SolvAI", face="white", edge=GREEN, weight="bold", size=11
    )
    rounded_box(
        axis,
        (0.785, 0.25),
        0.16,
        0.10,
        "Hydration free energy",
        face="white",
        edge=GREEN,
        weight="bold",
    )
    arrow(axis, (0.865, 0.66), (0.865, 0.57))
    arrow(axis, (0.865, 0.45), (0.865, 0.36))
    axis.text(
        0.865,
        0.15,
        "NO MD  •  NO PIMD  •  NO PROBE",
        ha="center",
        color=GREEN,
        weight="bold",
        size=8,
    )
    arrow(axis, (0.70, 0.50), (0.78, 0.51))
    save(fig, MAIN, "figure_1_concept")


def figure_2(metrics: dict) -> None:
    data = metrics["data_counts"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})
    names = ["MolSolv SMD", "ConfSolv H$_2$O"]
    raw = [data["molsolv_source_conformers"], data["confsolv_source_water_conformers"]]
    retained = [data["molsolv_training_structures"], data["confsolv_training_connectivities"]]
    y = np.arange(2)
    axes[0].barh(y + 0.16, raw, height=0.28, color="#C7CDD4", label="Source records")
    axes[0].barh(y - 0.16, retained, height=0.28, color=BLUE, label="Strict training structures")
    axes[0].set_xscale("log")
    axes[0].set_yticks(y, names)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Records or structures (log scale)")
    axes[0].legend(frameon=False, loc="lower right")
    for index, (source, keep) in enumerate(zip(raw, retained, strict=True)):
        axes[0].text(source * 1.05, index + 0.16, f"{source / 1e6:.2f}M", va="center", size=8)
        axes[0].text(keep * 1.05, index - 0.16, f"{keep:,}", va="center", size=8, color=BLUE)
    panel_label(axes[0], "a")

    axes[1].axis("off")
    panel_label(axes[1], "b")
    axes[1].text(0.5, 0.95, "Identity firewall", ha="center", weight="bold", size=10)
    rounded_box(
        axes[1], (0.12, 0.72), 0.76, 0.13, "External molecular structures", face=LIGHT, edge=GREY
    )
    rounded_box(
        axes[1],
        (0.12, 0.48),
        0.76,
        0.13,
        "Canonical SMILES + full InChIKey\n+ connectivity block",
        face="#EEF7FB",
        edge=BLUE,
    )
    rounded_box(
        axes[1],
        (0.12, 0.24),
        0.76,
        0.13,
        "Remove every ARROW-85 connectivity",
        face="#FFF2EE",
        edge=VERMILION,
        weight="bold",
    )
    axes[1].text(
        0.5,
        0.09,
        "0 benchmark overlaps in every retained source",
        ha="center",
        color=GREEN,
        weight="bold",
    )
    arrow(axes[1], (0.5, 0.71), (0.5, 0.62))
    arrow(axes[1], (0.5, 0.47), (0.5, 0.38))
    save(fig, MAIN, "figure_2_data_landscape")


def figure_3(metrics: dict) -> None:
    method = metrics["methods"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.5), gridspec_kw={"width_ratios": [1.7, 1]})
    labels = [
        "Previous\nstructure-only",
        "Narrow\nresponse",
        "+ SMD\nwater",
        "+ ConfSolv\nresponse",
    ]
    values = [
        method[key]["mae_kcal_mol"]
        for key in ("previous_structure_only", "narrow_response", "smd_water", "smd_confsolv_fixed")
    ]
    x = np.arange(len(values))
    axes[0].plot(x, values, color=BLUE, marker="o", linewidth=2.4, markersize=8)
    axes[0].fill_between(x, values, 0.25, color=SKY, alpha=0.12)
    axes[0].axhline(
        method["arrow_pimd8"]["mae_kcal_mol"],
        color=PURPLE,
        linestyle="--",
        linewidth=1.6,
        label="ARROW/PIMD8",
    )
    axes[0].axhline(0.20, color=DARK, linestyle=":", linewidth=1.2, label="0.20")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylim(0.185, 0.248)
    axes[0].set_ylabel("Strict five-fold OOF MAE (kcal/mol)")
    axes[0].legend(frameon=False, loc="upper right")
    for xx, value in zip(x, values, strict=True):
        axes[0].text(xx, value + 0.0025, f"{value:.3f}", ha="center", weight="bold", color=BLUE)
    panel_label(axes[0], "a")

    overview_labels = ["Classical\nARROW", "ARROW/\nPIMD8", "SolvAI"]
    overview_values = [
        method["classical_arrow"]["mae_kcal_mol"],
        method["arrow_pimd8"]["mae_kcal_mol"],
        method["smd_confsolv_fixed"]["mae_kcal_mol"],
    ]
    bars = axes[1].bar(overview_labels, overview_values, color=[GREY, PURPLE, GREEN], width=0.62)
    axes[1].set_ylim(0, 0.86)
    axes[1].set_ylabel("MAE (kcal/mol)")
    for bar, value in zip(bars, overview_values, strict=True):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.025,
            f"{value:.3f}",
            ha="center",
            weight="bold",
        )
    axes[1].text(2, 0.09, "SMILES only", ha="center", color=GREEN, weight="bold")
    panel_label(axes[1], "b")
    fig.tight_layout(w_pad=3)
    save(fig, MAIN, "figure_3_headline_result")


def figure_4(metrics: dict) -> None:
    alternatives = metrics["alternative_supervision"]
    base = alternatives["matched_one_seed_base"]
    labels = [
        "MLFF hierarchy",
        "OpenFE diagnostics",
        "GNNIS force",
        "Lambda-aware implicit",
        "DES370K/SAPT",
        "ConfSolv graph latent",
        "ConfSolv FFN latent",
    ]
    keys = [
        "mlff_hierarchy",
        "openfe_diagnostics",
        "gnnis_force_response",
        "lambda_aware_implicit",
        "des370k_water_response",
        "confsolv_graph_embedding",
        "confsolv_ffn_embedding",
    ]
    changes = np.asarray([alternatives[key] - base for key in keys])
    response = pd.DataFrame(metrics["multilambda"]["response_head_metrics"])
    response = response.loc[response.target_std.astype(float).gt(0)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), gridspec_kw={"width_ratios": [1.15, 1]})
    y = np.arange(len(labels))
    axes[0].barh(y, changes, color=[ORANGE if value > 0 else GREEN for value in changes])
    axes[0].axvline(0, color=DARK, linewidth=1)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Change from matched base MAE (kcal/mol)")
    axes[0].set_title("Additional response representations")
    panel_label(axes[0], "a")

    component_label = {
        "lig_slv__dhdl_mean": "total",
        "lig_slv__dhdl_coul_mean": "electrostatic",
        "lig_slv__dhdl_vdw_mean": "vdW",
    }
    for component, group in response.groupby("component", sort=False):
        axes[1].plot(
            group["lambda"], group.mae, marker="o", linewidth=1.8, label=component_label[component]
        )
    axes[1].axhline(0.20, color=DARK, linestyle="--", linewidth=1, label="endpoint target")
    axes[1].set_xticks([0.1, 0.5, 0.9])
    axes[1].set_xlabel("Alchemical coupling λ")
    axes[1].set_ylabel("Response-head MAE (kcal/mol)")
    axes[1].set_title("Structure → PIMD2 response remains limiting")
    axes[1].legend(frameon=False)
    panel_label(axes[1], "b")
    fig.tight_layout(w_pad=3)
    save(fig, MAIN, "figure_4_what_transfers")


def figure_5(metrics: dict) -> None:
    repeated = metrics["repeated_splits"]
    methods = metrics["methods"]
    family = pd.DataFrame(metrics["chemistry_family"]).sort_values("mae_kcal_mol")
    fig, axes = plt.subplots(
        1, 3, figsize=(12.2, 4.7), gridspec_kw={"width_ratios": [1, 0.8, 1.25]}
    )
    repeats = np.arange(1, 6)
    axes[0].plot(
        repeats, repeated["fixed"]["values_kcal_mol"], marker="o", color=BLUE, label="Fixed"
    )
    axes[0].plot(
        repeats, repeated["nested"]["values_kcal_mol"], marker="s", color=PURPLE, label="Nested"
    )
    axes[0].axhline(0.20, color=DARK, linestyle="--", linewidth=1)
    axes[0].set_xticks(repeats)
    axes[0].set_xlabel("Independent split")
    axes[0].set_ylabel("OOF MAE (kcal/mol)")
    axes[0].set_ylim(0.193, 0.214)
    axes[0].legend(frameon=False)
    panel_label(axes[0], "a")

    holdout_labels = ["Random\nOOF", "Family\nholdout", "Scaffold\nholdout"]
    holdout_values = [
        methods["smd_confsolv_fixed"]["mae_kcal_mol"],
        methods["family_holdout"]["mae_kcal_mol"],
        methods["scaffold_holdout"]["mae_kcal_mol"],
    ]
    bars = axes[1].bar(holdout_labels, holdout_values, color=[GREEN, ORANGE, ORANGE])
    axes[1].axhline(0.20, color=DARK, linestyle="--", linewidth=1)
    axes[1].set_ylim(0.18, 0.255)
    axes[1].set_ylabel("MAE (kcal/mol)")
    for bar, value in zip(bars, holdout_values, strict=True):
        axes[1].text(
            bar.get_x() + bar.get_width() / 2, value + 0.003, f"{value:.3f}", ha="center", size=8
        )
    panel_label(axes[1], "b")

    colors = [VERMILION if value > 0.20 else GREEN for value in family.mae_kcal_mol]
    axes[2].barh(np.arange(len(family)), family.mae_kcal_mol, color=colors)
    axes[2].axvline(0.20, color=DARK, linestyle="--", linewidth=1)
    axes[2].set_yticks(
        np.arange(len(family)),
        [f"{name} (n={n})" for name, n in zip(family.family, family.n, strict=True)],
    )
    axes[2].set_xlabel("Random-OOF MAE (kcal/mol)")
    axes[2].set_title("Chemical families")
    panel_label(axes[2], "c")
    fig.tight_layout(w_pad=2.4)
    save(fig, MAIN, "figure_5_robustness")


def extended_figures(metrics: dict) -> None:
    methods = metrics["methods"]
    benchmark = pd.read_parquet(ROOT / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark.solvent.eq("water")]
    headline = pd.read_parquet(ROOT / "results/predictions/headline_oof.parquet")
    final = headline.loc[headline.method.eq("Fixed narrow response + SMD + ConfSolv response")]

    fig, axis = plt.subplots(figsize=(6.2, 5.6))
    axis.scatter(
        final.y_true, final.y_pred, color=BLUE, alpha=0.78, s=28, edgecolor="white", linewidth=0.3
    )
    low = min(final.y_true.min(), final.y_pred.min()) - 0.4
    high = max(final.y_true.max(), final.y_pred.max()) + 0.4
    axis.plot([low, high], [low, high], color=DARK, linestyle="--", linewidth=1)
    axis.set(
        xlabel="Experimental hydration free energy (kcal/mol)",
        ylabel="SolvAI OOF prediction (kcal/mol)",
        xlim=(low, high),
        ylim=(low, high),
    )
    axis.text(
        0.04,
        0.96,
        f"MAE = {methods['smd_confsolv_fixed']['mae_kcal_mol']:.3f}\nn = 85",
        transform=axis.transAxes,
        va="top",
        weight="bold",
    )
    save(fig, EXTENDED, "extended_data_1_parity")

    data = metrics["data_counts"]
    labels = ["FreeSolv identity overlap", "MolSolv excluded", "ConfSolv excluded"]
    values = [
        data["arrow85_connectivities_also_in_freesolv"],
        data["molsolv_benchmark_structure_matches_removed"],
        data["confsolv_benchmark_connectivities_removed"],
    ]
    fig, axis = plt.subplots(figsize=(6.5, 4.3))
    bars = axis.bar(labels, values, color=[GREY, BLUE, SKY])
    axis.set_ylabel("ARROW-85 identities detected")
    axis.set_title("Identity matching precedes supervised training")
    axis.tick_params(axis="x", rotation=15)
    for bar, value in zip(bars, values, strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2, value + 2, str(value), ha="center", weight="bold"
        )
    axis.text(
        0.98,
        0.95,
        "Retained overlap = 0 for every source",
        transform=axis.transAxes,
        ha="right",
        va="top",
        color=GREEN,
        weight="bold",
    )
    fig.tight_layout()
    save(fig, EXTENDED, "extended_data_2_identity_audit")

    alternatives = metrics["alternative_supervision"]
    selected = [
        (key, value) for key, value in alternatives.items() if key != "matched_one_seed_base"
    ]
    selected.sort(key=lambda item: item[1])
    fig, axis = plt.subplots(figsize=(7.4, 5.2))
    labels = [key.replace("_", " ") for key, _ in selected]
    values = [value for _, value in selected]
    axis.barh(np.arange(len(values)), values, color=ORANGE)
    axis.axvline(
        alternatives["matched_one_seed_base"],
        color=BLUE,
        linewidth=1.5,
        label="Matched one-seed base",
    )
    axis.axvline(0.20, color=DARK, linestyle="--", linewidth=1, label="0.20")
    axis.set_yticks(np.arange(len(labels)), labels)
    axis.invert_yaxis()
    axis.set_xlabel("Random-OOF MAE (kcal/mol)")
    axis.legend(frameon=False)
    fig.tight_layout()
    save(fig, EXTENDED, "extended_data_3_alternative_supervision")

    selective = pd.read_csv(ROOT / "results/ablations/selective_compute_pareto.csv")
    policies = {
        "Fixed risk-rank selective PIMD": "fixed routing",
        "Nested learned selective PIMD": "nested learned routing",
        "Oracle selective PIMD": "oracle lower bound",
    }
    fig, axis = plt.subplots(figsize=(6.5, 4.5))
    selected = selective.loc[selective.regime.eq("random_oof") & selective.policy.isin(policies)]
    for policy, group in selected.groupby("policy", sort=False):
        group = group.sort_values("full_pimd_fraction")
        axis.plot(
            group.full_pimd_fraction * 100,
            group.mae,
            marker="o",
            label=policies[policy],
        )
    zero_sim = selective.loc[
        selective.regime.eq("random_oof") & selective.policy.eq("Zero-simulation fast model")
    ].iloc[0]
    axis.scatter([0], [zero_sim.mae], color=GREY, marker="s", zorder=4, label="fast baseline")
    axis.axhline(0.20, color=DARK, linestyle="--", linewidth=1)
    axis.set_xlabel("Full PIMD8 fallback (%)")
    axis.set_ylabel("MAE (kcal/mol)")
    axis.set_title("Non-deployable simulation-assisted comparison")
    axis.legend(frameon=False)
    fig.tight_layout()
    save(fig, EXTENDED, "extended_data_4_simulation_assisted")

    multilambda = metrics["multilambda"]["method_mae_kcal_mol"]
    labels = [
        "Matched base",
        "+ predicted response",
        "+ fidelity hierarchy",
        "+ both",
        "Integrated curve",
    ]
    keys = list(multilambda)
    values = [multilambda[key] for key in keys]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    axes[0].bar(labels[:4], values[:4], color=[GREEN, ORANGE, ORANGE, VERMILION])
    axes[0].axhline(0.20, color=DARK, linestyle="--", linewidth=1)
    axes[0].set_ylim(0.19, 0.225)
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].set_ylabel("OOF MAE (kcal/mol)")
    axes[1].bar(["Integrated\npredicted curve"], [values[4]], color=VERMILION, width=0.5)
    axes[1].set_ylim(0, 1.65)
    axes[1].set_ylabel("OOF MAE (kcal/mol)")
    fig.tight_layout(w_pad=3)
    save(fig, EXTENDED, "extended_data_5_multilambda")

    repeats = metrics["repeated_splits"]
    fig, axis = plt.subplots(figsize=(6.5, 4.3))
    positions = [0, 1]
    means = [repeats["fixed"]["mean_kcal_mol"], repeats["nested"]["mean_kcal_mol"]]
    errors = [repeats["fixed"]["sd_kcal_mol"], repeats["nested"]["sd_kcal_mol"]]
    axis.errorbar(positions, means, yerr=errors, fmt="o", color=BLUE, capsize=5, markersize=8)
    for index, values_for_method in enumerate(
        (repeats["fixed"]["values_kcal_mol"], repeats["nested"]["values_kcal_mol"])
    ):
        jitter = np.linspace(-0.07, 0.07, 5)
        axis.scatter(index + jitter, values_for_method, color=GREY, s=22, zorder=3)
    axis.axhline(0.20, color=DARK, linestyle="--", linewidth=1)
    axis.set_xticks(positions, ["Fixed SolvAI", "Nested selection"])
    axis.set_ylabel("MAE (kcal/mol)")
    axis.set_title("Five independent split repeats; mean ± SD")
    fig.tight_layout()
    save(fig, EXTENDED, "extended_data_6_repeat_statistics")

    bootstrap = metrics["bootstrap"]
    labels = ["Fixed SolvAI", "Nested selection"]
    estimates = [
        methods["smd_confsolv_fixed"]["mae_kcal_mol"],
        methods["nested_selection"]["mae_kcal_mol"],
    ]
    intervals = [bootstrap["fixed"]["ci95_kcal_mol"], bootstrap["nested"]["ci95_kcal_mol"]]
    lower = [
        estimate - interval[0] for estimate, interval in zip(estimates, intervals, strict=True)
    ]
    upper = [
        interval[1] - estimate for estimate, interval in zip(estimates, intervals, strict=True)
    ]
    fig, axis = plt.subplots(figsize=(6.2, 3.8))
    axis.errorbar(estimates, [0, 1], xerr=[lower, upper], fmt="o", color=BLUE, capsize=5)
    axis.axvline(0.20, color=DARK, linestyle="--", linewidth=1)
    axis.set_yticks([0, 1], labels)
    axis.set_xlabel("MAE with molecule-bootstrap 95% interval (kcal/mol)")
    axis.set_ylim(-0.6, 1.6)
    fig.tight_layout()
    save(fig, EXTENDED, "extended_data_7_bootstrap")

    fig, axis = plt.subplots(figsize=(8.4, 3.8))
    axis.axis("off")
    rounded_box(axis, (0.04, 0.39), 0.18, 0.20, "SMILES", face=LIGHT, edge=GREY, weight="bold")
    rounded_box(
        axis,
        (0.30, 0.57),
        0.25,
        0.18,
        "2,265 deterministic\nstructure descriptors",
        face="#EEF7FB",
        edge=BLUE,
    )
    rounded_box(
        axis, (0.30, 0.24), 0.25, 0.18, "15 predicted\nresponse priors", face="#FFF7E6", edge=ORANGE
    )
    rounded_box(axis, (0.64, 0.39), 0.15, 0.20, "3-model\nensemble", face="#EDF8F4", edge=GREEN)
    rounded_box(
        axis, (0.85, 0.39), 0.12, 0.20, "ΔG$_{hyd}$", face="white", edge=GREEN, weight="bold"
    )
    arrow(axis, (0.22, 0.49), (0.29, 0.66))
    arrow(axis, (0.22, 0.49), (0.29, 0.33))
    arrow(axis, (0.55, 0.66), (0.64, 0.53))
    arrow(axis, (0.55, 0.33), (0.64, 0.45))
    arrow(axis, (0.79, 0.49), (0.84, 0.49))
    axis.text(0.5, 0.90, "Released inference graph", ha="center", weight="bold", size=11)
    axis.text(
        0.5,
        0.08,
        "No benchmark table, trajectory, experimental value or simulation executable is read",
        ha="center",
        color=GREEN,
        weight="bold",
    )
    save(fig, EXTENDED, "extended_data_8_artifact_audit")


def main() -> None:
    configure()
    metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    figure_1(metrics)
    figure_2(metrics)
    figure_3(metrics)
    figure_4(metrics)
    figure_5(metrics)
    extended_figures(metrics)
    submission = ROOT / "paper/figures"
    submission.mkdir(parents=True, exist_ok=True)
    for path in [*MAIN.glob("*.pdf"), *EXTENDED.glob("*.pdf")]:
        shutil.copy2(path, submission / path.name)
    print("Generated 5 main and 8 Extended Data figures in PDF, SVG and PNG.")


if __name__ == "__main__":
    main()
