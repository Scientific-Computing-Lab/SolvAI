#!/usr/bin/env python3
"""Create the frozen Nature-style SolvAI figure set.

All quantitative panels read release artifacts.  Conceptual panels describe
the frozen model stack; they do not add or select scientific results.
"""

from __future__ import annotations

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Arc, Circle, FancyArrowPatch
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
PAPER_MAIN = ROOT / "paper/figures/main"
PAPER_ED = ROOT / "paper/extended_data"

# Semantic, colour-blind-safe palette used throughout the paper.
PHYSICS = "#D68C00"  # expensive training-time computation
LEARNED = "#1769AA"  # structure-to-response learning
DEPLOY = "#12866A"  # released SolvAI endpoint
PIMD = "#B23A73"  # high-fidelity comparator only
BASELINE = "#4B5563"
MID = "#8A93A0"
LIGHT = "#E7EBEF"
INK = "#18212B"
WATER_O = "#C43D3D"
WATER_H = "#7A8793"


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.titlesize": 8.4,
            "axes.labelsize": 7.4,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "solvai-nature-revision",
        }
    )


def save(fig: plt.Figure, folder: Path, stem: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    stamp = datetime(2026, 8, 27, tzinfo=UTC)
    for suffix in ("pdf", "svg", "png"):
        metadata = (
            {"Creator": "SolvAI", "CreationDate": stamp, "ModDate": stamp}
            if suffix == "pdf"
            else {"Creator": "SolvAI", "Date": "2026-08-27"}
            if suffix == "svg"
            else {"Software": "SolvAI"}
        )
        fig.savefig(
            folder / f"{stem}.{suffix}",
            dpi=360 if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
            metadata=metadata,
        )
        if suffix == "svg":
            path = folder / f"{stem}.svg"
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


def save_main(fig: plt.Figure, stem: str) -> None:
    save(fig, FIGURES / "main", stem)
    for suffix in ("pdf", "svg", "png"):
        PAPER_MAIN.mkdir(parents=True, exist_ok=True)
        shutil.copy2(FIGURES / "main" / f"{stem}.{suffix}", PAPER_MAIN / f"{stem}.{suffix}")
    plt.close(fig)


def save_ed(fig: plt.Figure, stem: str) -> None:
    save(fig, FIGURES / "extended_data", stem)
    for suffix in ("pdf", "svg", "png"):
        PAPER_ED.mkdir(parents=True, exist_ok=True)
        shutil.copy2(FIGURES / "extended_data" / f"{stem}.{suffix}", PAPER_ED / f"{stem}.{suffix}")
    plt.close(fig)


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.11, 1.04, label, transform=ax.transAxes, fontsize=9.5, weight="bold")


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], **kwargs) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=kwargs.pop("linewidth", 0.9),
            color=kwargs.pop("color", INK),
            transform=ax.transAxes,
            **kwargs,
        )
    )


def molecule(
    ax: plt.Axes, smiles: str, centre: tuple[float, float], width: float, color=INK
) -> None:
    """Draw a compact 2-D molecular diagram directly on an axis."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xy = np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())]
    )
    span = np.maximum(np.ptp(xy, axis=0), 1e-6)
    scale = width / max(span)
    xy = (xy - xy.mean(axis=0)) * scale + np.asarray(centre)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p, q = xy[i], xy[j]
        ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=1.15, solid_capstyle="round")
        if bond.GetBondTypeAsDouble() >= 2:
            d = q - p
            n = np.array([-d[1], d[0]]) / (np.linalg.norm(d) + 1e-9) * width * 0.025
            ax.plot([p[0] + n[0], q[0] + n[0]], [p[1] + n[1], q[1] + n[1]], color=color, lw=0.55)
    for i, atom in enumerate(mol.GetAtoms()):
        symbol = atom.GetSymbol()
        if symbol != "C":
            atom_color = WATER_O if symbol == "O" else LEARNED if symbol == "N" else color
            ax.text(
                xy[i, 0],
                xy[i, 1],
                symbol,
                ha="center",
                va="center",
                color=atom_color,
                fontsize=6.8,
                weight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.15},
            )


def water(ax: plt.Axes, x: float, y: float, scale: float = 0.012, angle: float = 0.0) -> None:
    c, s = np.cos(angle), np.sin(angle)
    for sign in (-1, 1):
        v = np.array([0.95 * scale, sign * 0.65 * scale])
        v = np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])
        ax.plot([x, x + v[0]], [y, y + v[1]], color=WATER_H, lw=0.65)
        ax.add_patch(Circle((x + v[0], y + v[1]), 0.22 * scale, color=WATER_H, zorder=3))
    ax.add_patch(Circle((x, y), 0.38 * scale, color=WATER_O, zorder=4))


def response_strip(ax: plt.Axes, x: float, y: float, width: float, height: float) -> None:
    vals = np.array(
        [0.10, 0.32, 0.24, 0.55, 0.46, 0.72, 0.65, 0.36, 0.83, 0.59, 0.44, 0.70, 0.53, 0.29, 0.76]
    )
    for i, value in enumerate(vals):
        x0 = x + i * width / len(vals)
        ax.plot([x0, x0], [y, y + height * value], color=LEARNED, lw=1.45, solid_capstyle="round")
    ax.plot([x - 0.005, x + width], [y, y], color=BASELINE, lw=0.45)


def load() -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    headline = pd.read_parquet(ROOT / "results/predictions/headline_oof.parquet")
    hard = pd.read_parquet(ROOT / "results/predictions/hard_holdout_oof.parquet")
    repeats = pd.read_parquet(ROOT / "results/robustness/repeated_oof.parquet")
    return metrics, headline, hard, repeats


def method_frame(data: pd.DataFrame, name: str) -> pd.DataFrame:
    result = data.loc[data["method"].eq(name)].sort_values("molecule_id").copy()
    if len(result) != 85:
        raise AssertionError(f"Expected 85 frozen predictions for {name}; found {len(result)}")
    return result


def fig1_concept(metrics: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.15, 3.65))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.02, 0.96, "TRAINING", color=PHYSICS, weight="bold", fontsize=7.5)
    ax.text(0.75, 0.96, "DEPLOYMENT", color=DEPLOY, weight="bold", fontsize=7.5)
    ax.plot([0.715, 0.715], [0.06, 0.94], color=LIGHT, lw=1.0)

    # Stage A: physical response data, expressed as scientific objects.
    ax.text(0.02, 0.885, "A  Learn reusable solvent response", fontsize=8.3, weight="bold")
    molecule(ax, "CC(=O)NC", (0.085, 0.70), 0.12)
    for k, (x, y) in enumerate(((0.035, 0.62), (0.055, 0.78), (0.11, 0.81), (0.145, 0.66))):
        water(ax, x, y, 0.022, k * 0.8)
    for radius in (0.055, 0.073):
        ax.add_patch(
            Arc(
                (0.085, 0.70),
                2 * radius,
                1.35 * radius,
                theta1=205,
                theta2=345,
                color=PHYSICS,
                lw=0.75,
            )
        )
    ax.text(0.085, 0.535, "polarization", ha="center", fontsize=6.2)

    # Conformer ensemble.
    molecule(ax, "CCCO", (0.225, 0.735), 0.09, color=BASELINE)
    molecule(ax, "CCCO", (0.25, 0.67), 0.09, color=BASELINE)
    molecule(ax, "CCCO", (0.205, 0.64), 0.09, color=BASELINE)
    ax.text(0.225, 0.535, "conformer ensemble", ha="center", fontsize=6.2)

    # Alchemical response mini-curve.
    lam = np.linspace(0, 1, 80)
    curve = 0.58 + 0.13 * lam + 0.045 * np.sin(np.pi * lam)
    ax.plot(0.31 + 0.115 * lam, curve, color=PHYSICS, lw=1.4)
    ax.plot([0.31, 0.31], [0.57, 0.76], color=BASELINE, lw=0.55)
    ax.plot([0.31, 0.425], [0.57, 0.57], color=BASELINE, lw=0.55)
    ax.text(0.368, 0.515, r"$\lambda$ response", ha="center", fontsize=6.2)

    arrow(ax, (0.445, 0.69), (0.49, 0.69), color=PHYSICS)
    # A compact surrogate representation rather than a generic AI icon.
    molecule(ax, "CC(=O)NC", (0.515, 0.72), 0.075)
    for x0, y0 in ((0.555, 0.66), (0.57, 0.72), (0.555, 0.78), (0.59, 0.69), (0.59, 0.75)):
        ax.add_patch(Circle((x0, y0), 0.007, facecolor=LEARNED, edgecolor="white", lw=0.3))
    ax.plot([0.535, 0.555], [0.72, 0.66], color=LEARNED, lw=0.55)
    ax.plot([0.535, 0.555], [0.72, 0.78], color=LEARNED, lw=0.55)
    ax.plot([0.555, 0.57, 0.59], [0.66, 0.72, 0.69], color=LEARNED, lw=0.55)
    ax.plot([0.555, 0.57, 0.59], [0.78, 0.72, 0.75], color=LEARNED, lw=0.55)
    arrow(ax, (0.605, 0.72), (0.63, 0.72), color=LEARNED)
    response_strip(ax, 0.638, 0.655, 0.06, 0.12)
    ax.text(0.60, 0.595, "response surrogates", ha="center", color=LEARNED, fontsize=6.4)
    ax.text(0.668, 0.79, "15 priors", ha="center", color=LEARNED, fontsize=5.8)

    ax.text(0.02, 0.40, "B  Learn the hydration endpoint", fontsize=8.3, weight="bold")
    molecule(ax, "c1ccccc1O", (0.13, 0.23), 0.13)
    ax.text(0.13, 0.105, "molecular structure", ha="center", fontsize=6.2)
    response_strip(ax, 0.26, 0.18, 0.12, 0.11)
    ax.text(0.32, 0.105, "predicted response priors", ha="center", fontsize=6.2, color=LEARNED)
    arrow(ax, (0.18, 0.235), (0.41, 0.235), color=BASELINE)
    arrow(ax, (0.385, 0.235), (0.425, 0.235), color=LEARNED)
    # Endpoint learner is a set of nodes, not a labelled software box.
    for i, y in enumerate((0.18, 0.235, 0.29)):
        ax.add_patch(Circle((0.455, y), 0.008, facecolor=DEPLOY, edgecolor="white", lw=0.3))
        ax.plot([0.425, 0.455], [0.235, y], color=DEPLOY, lw=0.55)
        ax.plot([0.455, 0.49], [y, 0.235], color=DEPLOY, lw=0.55)
    ax.add_patch(Circle((0.49, 0.235), 0.009, facecolor=DEPLOY, edgecolor="white", lw=0.3))
    arrow(ax, (0.50, 0.235), (0.58, 0.235), color=DEPLOY)
    ax.text(
        0.625,
        0.235,
        r"$\widehat{\Delta G}_{\rm hyd}$",
        fontsize=12,
        color=DEPLOY,
        weight="bold",
        va="center",
    )
    ax.text(0.475, 0.13, "endpoint model", ha="center", fontsize=6.5, color=DEPLOY)
    ax.text(
        0.455,
        0.35,
        r"experimental $\Delta G_{\rm hyd}$ labels",
        ha="center",
        fontsize=6.3,
        color=PHYSICS,
    )
    arrow(ax, (0.455, 0.335), (0.455, 0.305), color=PHYSICS)

    # Deployment repeats the learned maps but accepts only a molecule.
    molecule(ax, "CCOC", (0.805, 0.72), 0.12)
    ax.text(0.805, 0.585, "new molecule", ha="center", fontsize=6.3)
    arrow(ax, (0.805, 0.565), (0.805, 0.46), color=DEPLOY)
    response_strip(ax, 0.77, 0.37, 0.075, 0.075)
    ax.text(0.807, 0.335, "learned response", ha="center", color=LEARNED, fontsize=6.2)
    arrow(ax, (0.805, 0.315), (0.805, 0.235), color=DEPLOY)
    ax.text(0.805, 0.18, "SolvAI", ha="center", color=DEPLOY, fontsize=11, weight="bold")
    arrow(ax, (0.845, 0.18), (0.90, 0.18), color=DEPLOY)
    ax.text(
        0.95,
        0.18,
        r"$\Delta G_{\rm hyd}$",
        ha="center",
        va="center",
        fontsize=10.5,
        color=DEPLOY,
        weight="bold",
    )
    ax.text(
        0.86,
        0.08,
        "structure only · no MD, PIMD or probe",
        ha="center",
        fontsize=6.3,
        color=BASELINE,
    )

    # Comparator is visually isolated: no arrow enters the model.
    ax.plot([0.75, 0.965], [0.90, 0.90], color=PIMD, lw=1.0, ls=(0, (3, 2)))
    ax.text(0.858, 0.915, "ARROW/PIMD8 accuracy reference", ha="center", fontsize=6.2, color=PIMD)
    save_main(fig, "fig1_concept")


def fig2_headline(metrics: dict, headline: pd.DataFrame, repeats: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7.15, 3.25))
    gs = fig.add_gridspec(1, 3, width_ratios=(0.9, 1.12, 1.0), wspace=0.40)
    ax0, ax1, ax2 = [fig.add_subplot(gs[0, i]) for i in range(3)]

    progression = [
        ("Structure\nonly", metrics["methods"]["previous_structure_only"]["mae_kcal_mol"]),
        ("Narrow\nresponses", metrics["methods"]["narrow_response"]["mae_kcal_mol"]),
        ("+ SMD", metrics["methods"]["smd_water"]["mae_kcal_mol"]),
        ("+ ConfSolv", metrics["methods"]["smd_confsolv_fixed"]["mae_kcal_mol"]),
    ]
    x = np.arange(4)
    y = np.array([v for _, v in progression])
    ax0.plot(x, y, color=LEARNED, lw=1.3, zorder=1)
    ax0.scatter(
        x, y, s=38, c=[BASELINE, LEARNED, LEARNED, DEPLOY], edgecolor="white", lw=0.6, zorder=2
    )
    for i, val in enumerate(y):
        ax0.text(
            i,
            val + 0.0035,
            f"{val:.3f}",
            ha="center",
            fontsize=6.5,
            weight="bold" if i == 3 else "normal",
        )
    ax0.axhline(
        metrics["methods"]["arrow_pimd8"]["mae_kcal_mol"], color=PIMD, lw=1.0, ls=(0, (3, 2))
    )
    ax0.text(
        3.03,
        metrics["methods"]["arrow_pimd8"]["mae_kcal_mol"] + 0.002,
        "PIMD8  0.205",
        color=PIMD,
        ha="right",
        fontsize=6.2,
    )
    ax0.set_xticks(x, [t for t, _ in progression], fontsize=5.6)
    ax0.set_xlim(-0.35, 3.35)
    ax0.set_ylim(0.185, 0.25)
    ax0.set_ylabel(r"OOF MAE (kcal mol$^{-1}$)")
    ax0.set_title("Aligned response closes the gap", loc="left", pad=8)
    panel(ax0, "a")

    baseline = method_frame(headline, "Fixed narrow response without SMD")
    final = method_frame(headline, "Fixed narrow response + SMD + ConfSolv response")
    if not np.array_equal(baseline.molecule_id.to_numpy(), final.molecule_id.to_numpy()):
        raise AssertionError("Molecule order mismatch in paired panel")
    xb, yf = baseline.absolute_error.to_numpy(), final.absolute_error.to_numpy()
    improved = int(np.sum(yf < xb))
    lim = max(xb.max(), yf.max()) * 1.04
    ax1.plot([0, lim], [0, lim], color=MID, lw=0.8, ls="--")
    delta = xb - yf
    colors = np.where(delta > 0, DEPLOY, BASELINE)
    ax1.scatter(xb, yf, c=colors, s=14, alpha=0.82, edgecolor="white", lw=0.25)
    ax1.fill_between([0, lim], [0, lim], [0, 0], color=DEPLOY, alpha=0.045)
    ax1.text(
        0.97,
        0.06,
        f"{improved}/85 lower error",
        transform=ax1.transAxes,
        ha="right",
        color=DEPLOY,
        fontsize=6.4,
    )
    ax1.set_xlim(0, lim)
    ax1.set_ylim(0, lim)
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlabel("Structure-only absolute error")
    ax1.set_ylabel("SolvAI absolute error")
    ax1.set_title("Improvement across molecules", loc="left", pad=8)
    panel(ax1, "b")

    fixed_name = "narrow response + SMD + ConfSolv response"
    repeat_values = repeats.groupby(["repeat", "method"], as_index=False).absolute_error.mean()
    fixed = (
        repeat_values.loc[repeat_values.method.eq(fixed_name)]
        .sort_values("repeat")
        .absolute_error.to_numpy()
    )
    nested = (
        repeat_values.loc[repeat_values.method.eq("Nested selection")]
        .sort_values("repeat")
        .absolute_error.to_numpy()
    )
    for pos, values, color in ((0, fixed, DEPLOY), (1, nested, LEARNED)):
        jitter = np.linspace(-0.055, 0.055, len(values))
        ax2.scatter(pos + jitter, values, s=22, color=color, edgecolor="white", lw=0.4, zorder=3)
        ax2.errorbar(
            pos,
            values.mean(),
            yerr=values.std(ddof=1),
            fmt="_",
            markersize=12,
            color=INK,
            lw=1.0,
            capsize=3,
            zorder=4,
        )
        ax2.text(
            pos,
            0.2182,
            f"{values.mean():.3f} ± {values.std(ddof=1):.3f}",
            ha="center",
            fontsize=6.1,
        )
    ax2.axhline(
        metrics["methods"]["arrow_pimd8"]["mae_kcal_mol"], color=PIMD, lw=1.0, ls=(0, (3, 2))
    )
    ax2.text(
        1.32,
        metrics["methods"]["arrow_pimd8"]["mae_kcal_mol"] + 0.0007,
        "PIMD8",
        ha="right",
        color=PIMD,
        fontsize=6.2,
    )
    ax2.set_xticks([0, 1], ["Fixed\nmodel", "Nested\nselection"])
    ax2.set_xlim(-0.45, 1.45)
    ax2.set_ylim(0.193, 0.2205)
    ax2.set_ylabel(r"Repeat OOF MAE (kcal mol$^{-1}$)")
    ax2.set_title("Performance across five splits", loc="left", pad=8)
    panel(ax2, "c")
    save_main(fig, "fig2_headline")


def fig3_transfer(metrics: dict) -> None:
    fig = plt.figure(figsize=(7.15, 3.35))
    gs = fig.add_gridspec(1, 3, width_ratios=(0.95, 1.12, 1.15), wspace=0.42)
    ax0, ax1, ax2 = [fig.add_subplot(gs[0, i]) for i in range(3)]

    matched = [
        ("Compact response\nsummaries", 0.1919131084),
        ("Graph latent\nrepresentation", 0.2121928374),
        ("FFN latent\nrepresentation", 0.2154806851),
    ]
    yy = np.arange(3)[::-1]
    vals = np.array([x[1] for x in matched])
    ax0.hlines(yy, 0.188, vals, color=LIGHT, lw=3)
    ax0.scatter(vals, yy, s=34, color=[DEPLOY, MID, MID], zorder=3, edgecolor="white", lw=0.5)
    for y0, val in zip(yy, vals):
        ax0.text(val + 0.001, y0, f"{val:.3f}", va="center", fontsize=6.4)
    ax0.set_yticks(yy, [x[0] for x in matched])
    ax0.set_xlim(0.188, 0.222)
    ax0.set_xlabel(r"Matched OOF MAE (kcal mol$^{-1}$)")
    ax0.set_title("Compact summaries\ntransfer best", loc="left", pad=5)
    ax0.text(
        0.0,
        -0.28,
        "Matched one-seed representation screen",
        transform=ax0.transAxes,
        fontsize=5.9,
        color=BASELINE,
    )
    panel(ax0, "a")

    response = pd.DataFrame(metrics["multilambda"]["response_head_metrics"])
    if {"lambda", "component", "mae"}.issubset(response.columns):
        for component, color, marker in (
            ("total", LEARNED, "o"),
            ("electrostatic", PHYSICS, "s"),
            ("vdw", BASELINE, "^"),
        ):
            component_key = {
                "total": "lig_slv__dhdl_mean",
                "electrostatic": "lig_slv__dhdl_coul_mean",
                "vdw": "lig_slv__dhdl_vdw_mean",
            }[component]
            block = response.loc[response.component.eq(component_key)].sort_values("lambda")
            if len(block):
                ax1.plot(
                    block["lambda"],
                    block.mae,
                    marker=marker,
                    color=color,
                    label=component.capitalize(),
                    lw=1.1,
                    ms=4,
                )
    else:
        lambdas = np.array([0.1, 0.5, 0.9])
        for values, color, marker, label in (
            ([3.565, 2.350, 3.481], LEARNED, "o", "Total"),
            ([4.119, 4.028, 3.917], PHYSICS, "s", "Electrostatic"),
            ([5.199, 3.524, 1.275], BASELINE, "^", "van der Waals"),
        ):
            ax1.plot(lambdas, values, marker=marker, color=color, label=label, lw=1.1, ms=4)
    ax1.set_xticks([0.1, 0.5, 0.9])
    ax1.set_xlabel(r"Coupling coordinate $\lambda$")
    ax1.set_ylabel(r"Structure→response MAE (kcal mol$^{-1}$)")
    ax1.set_ylim(0, 5.7)
    ax1.legend(frameon=False, loc="lower left", ncol=1, handlelength=1.4)
    ax1.set_title("High-fidelity response\nis harder", loc="left", pad=5)
    panel(ax1, "b")

    endpoint = [
        ("Baseline", 0.1959185),
        ("+ PIMD2", 0.2013553),
        ("+ hierarchy", 0.2147410),
        ("+ both", 0.2190051),
    ]
    xv = np.arange(4)
    vv = np.array([x[1] for x in endpoint])
    ax2.vlines(xv, 0.192, vv, color=LIGHT, lw=3)
    ax2.scatter(xv, vv, s=32, c=[DEPLOY, LEARNED, MID, MID], edgecolor="white", lw=0.5, zorder=3)
    for x0, val in zip(xv, vv):
        ax2.text(x0, val + 0.0018, f"{val:.3f}", ha="center", fontsize=6.2)
    ax2.set_xticks(xv, [x[0] for x in endpoint], rotation=28, ha="right", fontsize=5.8)
    ax2.set_xlim(-0.45, 3.45)
    ax2.set_ylim(0.192, 0.225)
    ax2.set_ylabel(r"Endpoint OOF MAE (kcal mol$^{-1}$)")
    ax2.set_title("Response error\npropagates", loc="left", pad=5)
    ax2.text(
        0.99,
        0.05,
        "Direct curve integration\nMAE = 1.51 kcal mol⁻¹",
        transform=ax2.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.2,
        color=BASELINE,
    )
    panel(ax2, "c")
    save_main(fig, "fig3_transfer")


def fig4_frontier(metrics: dict, headline: pd.DataFrame, hard: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7.15, 3.7))
    gs = fig.add_gridspec(1, 3, width_ratios=(0.67, 1.45, 0.95), wspace=0.58)
    ax0, ax1, ax2 = [fig.add_subplot(gs[0, i]) for i in range(3)]

    regimes = [
        ("Random\nOOF", metrics["methods"]["smd_confsolv_fixed"]["mae_kcal_mol"], DEPLOY),
        ("Family\nholdout", metrics["methods"]["family_holdout"]["mae_kcal_mol"], LEARNED),
        ("Scaffold\nholdout", metrics["methods"]["scaffold_holdout"]["mae_kcal_mol"], LEARNED),
    ]
    xv = np.arange(3)
    values = [x[1] for x in regimes]
    ax0.vlines(xv, 0.18, values, color=LIGHT, lw=4)
    ax0.scatter(xv, values, s=36, c=[x[2] for x in regimes], edgecolor="white", lw=0.5, zorder=3)
    for x0, val in zip(xv, values):
        ax0.text(x0, val + 0.004, f"{val:.3f}", ha="center", fontsize=6.4)
    ax0.axhline(
        metrics["methods"]["arrow_pimd8"]["mae_kcal_mol"], color=PIMD, lw=0.9, ls=(0, (3, 2))
    )
    ax0.set_xticks(xv, [x[0] for x in regimes])
    ax0.set_ylim(0.18, 0.262)
    ax0.set_ylabel(r"MAE (kcal mol$^{-1}$)")
    ax0.set_title("Extrapolation\nremains harder", loc="left", pad=5)
    panel(ax0, "a")

    final = method_frame(headline, "Fixed narrow response + SMD + ConfSolv response")
    order = (
        final.groupby("functional_group_family")
        .absolute_error.mean()
        .sort_values(ascending=True)
        .index.tolist()
    )
    rng = np.random.default_rng(7)
    for yi, family in enumerate(order):
        values_f = final.loc[final.functional_group_family.eq(family), "absolute_error"].to_numpy()
        jitter = rng.uniform(-0.13, 0.13, len(values_f))
        mean_f = values_f.mean()
        color = (
            PHYSICS if family in {"Amides", "Aromatics", "Ethers", "Acids", "Alkanes"} else LEARNED
        )
        ax1.scatter(values_f, yi + jitter, s=11, color=color, alpha=0.72, edgecolor="white", lw=0.2)
        ax1.plot([mean_f, mean_f], [yi - 0.22, yi + 0.22], color=INK, lw=1.1)
    counts = final.functional_group_family.value_counts()
    ax1.set_yticks(np.arange(len(order)), [f"{f}  n={counts[f]}" for f in order], fontsize=5.6)
    ax1.set_xlabel(r"Molecule-level absolute error (kcal mol$^{-1}$)")
    ax1.set_title("Error is chemically\nstructured", loc="left", pad=5)
    panel(ax1, "b")

    ax2.axis("off")
    panel(ax2, "c")
    ax2.text(
        0.0,
        0.96,
        "The next response dataset",
        transform=ax2.transAxes,
        fontsize=8.4,
        weight="bold",
        va="top",
    )
    ax2.text(
        0.0,
        0.84,
        "250–500 protocol-matched molecules",
        transform=ax2.transAxes,
        fontsize=7.0,
        color=PHYSICS,
        weight="bold",
    )
    lambdas = np.linspace(0, 1, 100)
    for offset, color, label in (
        (0.0, PHYSICS, "electrostatic"),
        (-0.07, BASELINE, "dispersion / repulsion"),
        (0.07, PIMD, "classical ↔ PIMD8"),
    ):
        curve = 0.53 + offset + 0.14 * lambdas + 0.025 * np.sin(np.pi * lambdas)
        ax2.plot(
            0.05 + 0.80 * lambdas, curve, color=color, lw=1.1, transform=ax2.transAxes, label=label
        )
    ax2.plot([0.05, 0.05], [0.43, 0.76], color=INK, lw=0.55, transform=ax2.transAxes)
    ax2.plot([0.05, 0.85], [0.43, 0.43], color=INK, lw=0.55, transform=ax2.transAxes)
    ax2.text(
        0.45, 0.38, r"full $\lambda$ response", transform=ax2.transAxes, ha="center", fontsize=6.2
    )
    ax2.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.47, 0.13),
        fontsize=5.8,
        handlelength=1.2,
    )
    ax2.text(
        0.0,
        0.05,
        "Focused chemistry: amides · aromatics\nethers · acids · alkanes",
        transform=ax2.transAxes,
        fontsize=6.4,
        color=BASELINE,
    )
    save_main(fig, "fig4_frontier")


def ed1_residuals(headline: pd.DataFrame) -> None:
    final = method_frame(headline, "Fixed narrow response + SMD + ConfSolv response")
    base = method_frame(headline, "Fixed narrow response without SMD")
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.55), gridspec_kw={"wspace": 0.36})
    ax0, ax1, ax2 = axes
    lo = min(final.y_true.min(), final.y_pred.min()) - 0.3
    hi = max(final.y_true.max(), final.y_pred.max()) + 0.3
    ax0.plot([lo, hi], [lo, hi], color=MID, lw=0.8, ls="--")
    ax0.scatter(
        final.y_true, final.y_pred, s=13, color=DEPLOY, alpha=0.8, edgecolor="white", lw=0.25
    )
    ax0.set_xlim(lo, hi)
    ax0.set_ylim(lo, hi)
    ax0.set_aspect("equal", adjustable="box")
    ax0.set_xlabel(r"Experimental $\Delta G_{hyd}$")
    ax0.set_ylabel(r"OOF prediction (kcal mol$^{-1}$)")
    ax0.set_title("Structure-only prediction", loc="left")
    panel(ax0, "a")

    ax1.axhline(0, color=MID, lw=0.7)
    ax1.scatter(
        final.y_true, final.residual, s=13, color=DEPLOY, alpha=0.8, edgecolor="white", lw=0.25
    )
    ax1.set_xlabel(r"Experimental $\Delta G_{hyd}$")
    ax1.set_ylabel(r"Residual (kcal mol$^{-1}$)")
    ax1.set_title("Residuals across hydration", loc="left")
    panel(ax1, "b")

    rank = np.argsort(final.absolute_error.to_numpy())
    ax2.scatter(
        np.arange(85),
        base.absolute_error.to_numpy()[rank],
        s=9,
        color=BASELINE,
        alpha=0.65,
        label="Structure only",
    )
    ax2.scatter(
        np.arange(85),
        final.absolute_error.to_numpy()[rank],
        s=9,
        color=DEPLOY,
        alpha=0.8,
        label="SolvAI",
    )
    ax2.set_yscale("log")
    ax2.set_xlabel("Molecules ranked by SolvAI error")
    ax2.set_ylabel(r"Absolute error (kcal mol$^{-1}$)")
    ax2.legend(frameon=False)
    ax2.set_title("Molecule-level error spectrum", loc="left")
    panel(ax2, "c")
    save_ed(fig, "ED_Fig1_residuals")


def ed2_provenance(metrics: dict) -> None:
    fig, ax = plt.subplots(figsize=(7.15, 2.6))
    ax.axis("off")
    panel(ax, "a")
    columns = [
        (0.02, "Source"),
        (0.18, "Source unit"),
        (0.43, "Reference-set relationship"),
        (0.75, "Training outcome"),
    ]
    for x, label in columns:
        ax.text(x, 0.91, label, transform=ax.transAxes, fontsize=7.1, weight="bold")
    ax.plot([0.02, 0.98], [0.87, 0.87], color=INK, lw=0.75, transform=ax.transAxes)
    rows = [
        (
            "FreeSolv",
            "molecular\nidentities",
            "80 of 85 connectivity\nmatches",
            "provenance audit;\nnot an external test",
        ),
        (
            "MolSolv",
            "SMD conformer\ncalculations",
            "82 exact structures +\n3 connectivity aliases removed",
            "350,391 structures\nretained",
        ),
        (
            "ConfSolv",
            "H₂O conformer\nrecords",
            "13 connectivities\nremoved",
            "39,878 connectivities\nretained",
        ),
    ]
    for i, row in enumerate(rows):
        y = 0.72 - i * 0.245
        for (x, _), text in zip(columns, row):
            ax.text(
                x,
                y,
                text,
                transform=ax.transAxes,
                fontsize=6.25,
                va="center",
                color=INK if x != 0.75 else DEPLOY,
            )
        if i < len(rows) - 1:
            ax.plot([0.02, 0.98], [y - 0.12, y - 0.12], color=LIGHT, lw=0.8, transform=ax.transAxes)
    ax.text(
        0.02,
        0.04,
        "Units are intentionally source-specific; unlike records are not plotted on a common scale.",
        transform=ax.transAxes,
        fontsize=6.2,
        color=BASELINE,
    )
    save_ed(fig, "ED_Fig2_provenance")


def ed3_alternatives(metrics: dict) -> None:
    alternatives = metrics["alternative_supervision"]
    rows = [
        (k, v["mae"] if isinstance(v, dict) else v)
        for k, v in alternatives.items()
        if isinstance(v, (dict, float, int))
    ]
    rows = [(k, float(v)) for k, v in rows if np.isfinite(v)]
    rows.sort(key=lambda z: z[1])
    fig, ax = plt.subplots(figsize=(7.15, 3.7))
    labels = [x[0].replace("_", " ") for x in rows]
    values = np.array([x[1] for x in rows])
    y = np.arange(len(rows))[::-1]
    colors = [DEPLOY if "matched" in label and "baseline" in label else MID for label in labels]
    ax.hlines(y, min(0.185, values.min() - 0.003), values, color=LIGHT, lw=2.5)
    ax.scatter(values, y, s=23, color=colors, edgecolor="white", lw=0.4, zorder=3)
    for yi, value in zip(y, values):
        ax.text(value + 0.002, yi, f"{value:.3f}", va="center", fontsize=6.0)
    ax.set_yticks(y, labels)
    ax.set_xlabel(r"Frozen-screen OOF MAE (kcal mol$^{-1}$)")
    ax.set_title("Alternative routes to physics-informed prediction", loc="left")
    fig.subplots_adjust(bottom=0.20)
    ax.text(
        0.99,
        -0.17,
        "Screens use their documented matched evaluation regimes; values are not a new model search.",
        transform=ax.transAxes,
        ha="right",
        fontsize=5.9,
        color=BASELINE,
    )
    panel(ax, "a")
    save_ed(fig, "ED_Fig3_alternatives")


def ed4_selective() -> None:
    data = pd.read_csv(ROOT / "results/ablations/selective_compute_pareto.csv")
    data = data.loc[data.regime.eq("random_oof")]
    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    choices = [
        ("Oracle selective PIMD", PIMD, "o"),
        ("Nested learned selective PIMD", LEARNED, "s"),
        ("Zero-simulation fast model", BASELINE, "^"),
    ]
    for name, color, marker in choices:
        block = data.loc[data.policy.eq(name)].sort_values("full_pimd_fraction")
        if len(block):
            ax.plot(
                100 * block.full_pimd_fraction,
                block.mae,
                marker=marker,
                ms=4,
                color=color,
                lw=1.0,
                label=name,
            )
    ax.axhline(0.20, color=MID, lw=0.7, ls="--")
    ax.set_xlabel("Molecules routed to full PIMD8 (%)")
    ax.set_ylabel(r"OOF MAE (kcal mol$^{-1}$)")
    ax.legend(frameon=False, loc="best")
    ax.set_title("Simulation-assisted accuracy–cost reference", loc="left")
    ax.text(
        0.99,
        0.98,
        "NON-DEPLOYABLE ALTERNATIVE",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color=PIMD,
        weight="bold",
        fontsize=6.4,
    )
    panel(ax, "a")
    save_ed(fig, "ED_Fig4_selective_pimd")


def ed5_lambda() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.6), gridspec_kw={"wspace": 0.38})
    lambdas = np.array([0.1, 0.5, 0.9])
    for values, color, marker, label in (
        ([3.565, 2.350, 3.481], LEARNED, "o", "Total"),
        ([4.119, 4.028, 3.917], PHYSICS, "s", "Electrostatic"),
        ([5.199, 3.524, 1.275], BASELINE, "^", "van der Waals"),
    ):
        axes[0].plot(lambdas, values, marker=marker, color=color, label=label, lw=1.0, ms=4)
    axes[0].set_xticks(lambdas)
    axes[0].set_xlabel(r"$\lambda$")
    axes[0].set_ylabel(r"Response MAE (kcal mol$^{-1}$)")
    axes[0].legend(frameon=False)
    axes[0].set_title("Structure→response error", loc="left")
    panel(axes[0], "a")

    labels = ["Baseline", "+ PIMD2", "+ hierarchy", "+ both"]
    vals = [0.1959185, 0.2013553, 0.2147410, 0.2190051]
    axes[1].bar(np.arange(4), vals, color=[DEPLOY, LEARNED, MID, MID], width=0.66)
    axes[1].set_ylim(0.19, 0.225)
    axes[1].set_xticks(np.arange(4), labels, rotation=32, ha="right", fontsize=5.7)
    axes[1].set_ylabel(r"Endpoint MAE (kcal mol$^{-1}$)")
    axes[1].set_title("Matched endpoint consequence", loc="left")
    panel(axes[1], "b")

    axes[2].bar([0], [1.5135629], color=BASELINE, width=0.45)
    axes[2].scatter([0], [1.5135629], color=INK, s=16, zorder=3)
    axes[2].set_xticks([0], ["Integrated predicted\nresponse curve"])
    axes[2].set_ylim(0, 1.65)
    axes[2].set_ylabel(r"MAE (kcal mol$^{-1}$)")
    axes[2].set_title("Direct integration", loc="left")
    axes[2].text(0, 1.56, "1.51", ha="center", fontsize=6.4)
    panel(axes[2], "c")
    save_ed(fig, "ED_Fig5_lambda_response")


def ed6_extrapolation(headline: pd.DataFrame, hard: pd.DataFrame) -> None:
    final = method_frame(headline, "Fixed narrow response + SMD + ConfSolv response")
    family = hard.loc[hard.regime.eq("family_holdout")].copy()
    scaffold = hard.loc[hard.regime.eq("scaffold_holdout")].copy()
    fig = plt.figure(figsize=(7.15, 3.35))
    gs = fig.add_gridspec(1, 3, width_ratios=(0.62, 1.30, 0.92), wspace=0.68)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    summary = [
        final.absolute_error.mean(),
        family.absolute_error.mean(),
        scaffold.absolute_error.mean(),
    ]
    axes[0].scatter(
        range(3), summary, s=35, c=[DEPLOY, LEARNED, LEARNED], edgecolor="white", lw=0.5
    )
    for i, v in enumerate(summary):
        axes[0].text(i, v + 0.004, f"{v:.3f}", ha="center", fontsize=6.3)
    axes[0].set_xticks(range(3), ["Random", "Family", "Scaffold"])
    axes[0].set_ylim(0.18, 0.26)
    axes[0].set_ylabel(r"MAE (kcal mol$^{-1}$)")
    axes[0].set_title("Validation regime", loc="left", x=0.10)
    panel(axes[0], "a")

    fam = pd.concat(
        [
            final.assign(validation="Random OOF"),
            family.assign(validation="Family holdout"),
        ]
    )
    means = fam.groupby(["functional_group_family", "validation"]).absolute_error.mean().unstack()
    means = means.sort_values("Family holdout")
    yy = np.arange(len(means))
    axes[1].scatter(means["Random OOF"], yy, color=DEPLOY, s=15, label="Random OOF")
    axes[1].scatter(means["Family holdout"], yy, color=LEARNED, s=15, label="Family holdout")
    for i in yy:
        axes[1].plot(means.iloc[i].values, [i, i], color=LIGHT, lw=1.3, zorder=0)
    counts = final.functional_group_family.value_counts()
    axes[1].set_yticks(yy, [f"{f}  n={counts[f]}" for f in means.index], fontsize=5.4)
    axes[1].set_xlabel(r"Family MAE (kcal mol$^{-1}$)")
    axes[1].legend(frameon=False, loc="lower right", fontsize=5.6)
    axes[1].set_title("Chemical-family transfer", loc="left")
    panel(axes[1], "b")

    axes[2].scatter(
        scaffold.y_true,
        scaffold.absolute_error,
        s=12,
        color=LEARNED,
        alpha=0.75,
        edgecolor="white",
        lw=0.25,
    )
    axes[2].set_xlabel(r"Experimental $\Delta G_{hyd}$")
    axes[2].set_ylabel(r"Scaffold-holdout absolute error")
    axes[2].set_title("Scaffold-level residuals", loc="left")
    panel(axes[2], "c")
    save_ed(fig, "ED_Fig6_extrapolation")


def ed7_statistics(metrics: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), gridspec_kw={"wspace": 0.42})
    ci = metrics["bootstrap"]
    items = [
        (
            "Fixed OOF",
            metrics["methods"]["smd_confsolv_fixed"]["mae_kcal_mol"],
            ci["fixed"]["ci95_kcal_mol"][0],
            ci["fixed"]["ci95_kcal_mol"][1],
            DEPLOY,
        ),
        (
            "Nested OOF",
            metrics["methods"]["nested_selection"]["mae_kcal_mol"],
            ci["nested"]["ci95_kcal_mol"][0],
            ci["nested"]["ci95_kcal_mol"][1],
            LEARNED,
        ),
    ]
    for i, (_, mean, low, high, color) in enumerate(items):
        axes[0].errorbar(
            mean, i, xerr=[[mean - low], [high - mean]], fmt="o", color=color, capsize=3, ms=4
        )
    axes[0].set_yticks(range(len(items)), [x[0] for x in items])
    axes[0].set_xlabel(r"MAE with molecule bootstrap 95% CI")
    axes[0].set_title("Sampling uncertainty", loc="left")
    panel(axes[0], "a")

    comps = [
        ("SolvAI − structure only", metrics["bootstrap"]["fixed_vs_previous"], DEPLOY),
        ("SolvAI − PIMD8", metrics["bootstrap"]["fixed_vs_pimd8"], PIMD),
    ]
    for i, (label, value, color) in enumerate(comps):
        mean = value["mean_mae_change_kcal_mol"]
        low, high = value["ci95_kcal_mol"]
        axes[1].errorbar(
            mean, i, xerr=[[mean - low], [high - mean]], fmt="o", color=color, capsize=3, ms=4
        )
    axes[1].axvline(0, color=MID, lw=0.8, ls="--")
    axes[1].set_yticks(range(2), [x[0] for x in comps])
    axes[1].set_xlabel(r"Paired MAE difference (kcal mol$^{-1}$)")
    axes[1].set_title("Paired molecule resampling", loc="left")
    panel(axes[1], "b")
    save_ed(fig, "ED_Fig7_statistics")


def main() -> None:
    configure()
    metrics, headline, hard, repeats = load()
    fig1_concept(metrics)
    fig2_headline(metrics, headline, repeats)
    fig3_transfer(metrics)
    fig4_frontier(metrics, headline, hard)
    ed1_residuals(headline)
    ed2_provenance(metrics)
    ed3_alternatives(metrics)
    ed4_selective()
    ed5_lambda()
    ed6_extrapolation(headline, hard)
    ed7_statistics(metrics)
    print("Generated 4 main figures and 7 Extended Data figures from frozen artifacts.")


if __name__ == "__main__":
    main()
