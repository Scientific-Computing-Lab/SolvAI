#!/usr/bin/env python3
"""Create the Nature Communications figure set from frozen result tables."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "figures/main"
ED = ROOT / "figures/extended_data"
PAPER_MAIN = ROOT / "paper/figures/main"
PAPER_ED = ROOT / "paper/extended_data"

INK = "#17212B"
MID = "#68737D"
GRID = "#D9DEE2"
PHYSICS = "#D68C2E"
PHYSICS_LIGHT = "#F4D9B4"
LEARNED = "#2878B5"
LEARNED_LIGHT = "#C8E0F1"
DEPLOY = "#14927D"
DEPLOY_LIGHT = "#BFE4DC"
PIMD = "#B43B75"
NEGATIVE = "#B95C50"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "axes.titlesize": 8.0,
        "axes.titleweight": "bold",
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "solvai-publication-release",
        "savefig.dpi": 450,
    }
)


def clean(ax: plt.Axes, *, grid: str | None = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.55, zorder=0)
    ax.tick_params(color=MID, labelcolor=INK)


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.07,
        label,
        transform=ax.transAxes,
        fontsize=9.5,
        fontweight="bold",
        va="top",
    )


def save(fig: plt.Figure, name: str, *, extended: bool = False) -> None:
    targets = (ED, PAPER_ED) if extended else (MAIN, PAPER_MAIN)
    for directory in targets:
        directory.mkdir(parents=True, exist_ok=True)
        for suffix in ("pdf", "svg", "png"):
            output = directory / f"{name}.{suffix}"
            metadata = {
                "pdf": {"CreationDate": None, "ModDate": None, "Creator": "SolvAI"},
                "svg": {"Date": None, "Creator": "SolvAI"},
                "png": {"Software": "SolvAI"},
            }[suffix]
            fig.savefig(
                output,
                bbox_inches="tight",
                pad_inches=0.035,
                dpi=450,
                transparent=False,
                metadata=metadata,
            )
            # Matplotlib writes path data with line-ending spaces. Normalizing
            # generated SVG text keeps repository whitespace checks useful.
            if suffix == "svg":
                lines = output.read_text().splitlines()
                output.write_text("\n".join(line.rstrip() for line in lines) + "\n")
    plt.close(fig)


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = INK,
    width: float = 1.0,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            lw=width,
            color=color,
            shrinkA=0,
            shrinkB=0,
        )
    )


def draw_molecule(
    ax: plt.Axes,
    smiles: str,
    centre: tuple[float, float],
    scale: float,
    *,
    alpha: float = 1.0,
) -> None:
    molecule = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(molecule)
    coordinates = molecule.GetConformer().GetPositions()[:, :2]
    span = np.ptp(coordinates, axis=0)
    normalization = max(float(span.max()), 1.0)
    coordinates = (coordinates - coordinates.mean(axis=0)) / normalization * scale
    coordinates[:, 0] += centre[0]
    coordinates[:, 1] += centre[1]
    for bond in molecule.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ax.plot(
            coordinates[[i, j], 0],
            coordinates[[i, j], 1],
            color=INK,
            lw=1.25,
            alpha=alpha,
            solid_capstyle="round",
            zorder=4,
        )
    atom_colors = {"O": NEGATIVE, "N": LEARNED, "S": PHYSICS, "Cl": DEPLOY}
    for atom, (x, y) in zip(molecule.GetAtoms(), coordinates, strict=True):
        symbol = atom.GetSymbol()
        if symbol != "C":
            ax.text(
                x,
                y,
                symbol,
                ha="center",
                va="center",
                fontsize=6.6,
                color=atom_colors.get(symbol, INK),
                weight="bold",
                alpha=alpha,
                zorder=5,
            )


def draw_water(ax: plt.Axes, x: float, y: float, angle: float = 0.0) -> None:
    angle = np.deg2rad(angle)
    oxygen = np.array([x, y])
    offsets = np.array([[-0.016, 0.019], [0.016, 0.019]])
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    hydrogens = oxygen + offsets @ rotation.T
    for hydrogen in hydrogens:
        ax.plot(
            [oxygen[0], hydrogen[0]],
            [oxygen[1], hydrogen[1]],
            color=MID,
            lw=0.65,
            zorder=1,
        )
        ax.add_patch(Circle(hydrogen, 0.006, color="#E9EEF1", ec=MID, lw=0.25))
    ax.add_patch(Circle(oxygen, 0.010, color=NEGATIVE, ec="white", lw=0.3, zorder=2))


def response_strip(ax: plt.Axes, x: float, y: float, width: float, height: float) -> None:
    values = np.array([0.35, 0.72, 0.48, 0.86, 0.58, 0.28, 0.67, 0.42])
    for index, value in enumerate(values):
        bar_x = x + width * index / len(values)
        ax.plot(
            [bar_x, bar_x],
            [y, y + height * value],
            color=LEARNED,
            lw=2.0,
            solid_capstyle="round",
        )


def fig1_concept() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.65))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.02,
        0.96,
        "TRAINING: physical response is generated once",
        color=PHYSICS,
        weight="bold",
        fontsize=8,
    )
    ax.text(0.69, 0.96, "DEPLOYMENT", color=DEPLOY, weight="bold", fontsize=8)
    ax.plot([0.655, 0.655], [0.06, 0.96], color=GRID, lw=0.8)

    # Three physical response vignettes, without workflow boxes.
    draw_molecule(ax, "CC(=O)NC", (0.105, 0.73), 0.10)
    for x, y, angle in ((0.03, 0.80, 25), (0.17, 0.80, -25), (0.03, 0.65, 150), (0.18, 0.65, 210)):
        draw_water(ax, x, y, angle)
    ax.add_patch(Circle((0.105, 0.73), 0.115, fill=False, ec=PHYSICS, lw=0.8, ls=(0, (2, 2))))
    ax.text(0.105, 0.56, "water response", ha="center", color=PHYSICS, fontsize=6.8)

    for offset, alpha in ((-0.028, 0.35), (0.0, 0.65), (0.028, 1.0)):
        draw_molecule(ax, "CCCO", (0.32 + offset, 0.73), 0.085, alpha=alpha)
    ax.annotate(
        "",
        xy=(0.37, 0.79),
        xytext=(0.27, 0.66),
        arrowprops={"arrowstyle": "<->", "color": PHYSICS, "lw": 0.8},
    )
    ax.text(0.32, 0.56, "conformer response", ha="center", color=PHYSICS, fontsize=6.8)

    lam = np.linspace(0.0, 1.0, 100)
    response = 0.64 + 0.085 * np.sin(np.pi * lam) - 0.13 * lam
    ax.plot(0.45 + 0.14 * lam, response, color=PHYSICS, lw=1.3)
    ax.fill_between(0.45 + 0.14 * lam, 0.60, response, color=PHYSICS_LIGHT, alpha=0.75)
    ax.text(0.45, 0.59, r"$\lambda=0$", ha="center", fontsize=5.8, color=MID)
    ax.text(0.59, 0.59, r"$\lambda=1$", ha="center", fontsize=5.8, color=MID)
    ax.text(0.52, 0.56, "alchemical response", ha="center", color=PHYSICS, fontsize=6.8)

    ax.text(
        0.31, 0.88, "external calculations and measurements", ha="center", color=MID, fontsize=6.3
    )
    arrow(ax, (0.31, 0.52), (0.31, 0.42), PHYSICS)
    ax.text(
        0.31,
        0.445,
        "structure → response surrogates",
        ha="center",
        color=LEARNED,
        weight="bold",
        fontsize=7.2,
    )
    response_strip(ax, 0.255, 0.31, 0.13, 0.08)
    ax.text(0.32, 0.285, "15 predicted response priors", ha="center", color=LEARNED, fontsize=6.6)

    draw_molecule(ax, "CCOC", (0.105, 0.25), 0.085)
    ax.text(0.105, 0.14, "molecular structure", ha="center", fontsize=6.5)
    arrow(ax, (0.17, 0.25), (0.235, 0.25), LEARNED)
    arrow(ax, (0.39, 0.35), (0.44, 0.29), LEARNED)
    ax.text(0.43, 0.36, "+", fontsize=11, weight="bold", color=MID)
    ax.text(
        0.53, 0.405, r"experimental $\Delta G_{\rm hyd}$", ha="center", color=PHYSICS, fontsize=6.5
    )
    arrow(ax, (0.53, 0.38), (0.53, 0.30), PHYSICS)
    for y in (0.20, 0.25, 0.30):
        ax.add_patch(Circle((0.48, y), 0.007, color=DEPLOY))
        ax.plot([0.44, 0.48], [0.25, y], color=DEPLOY, lw=0.6)
        ax.plot([0.48, 0.55], [y, 0.25], color=DEPLOY, lw=0.6)
    ax.text(
        0.50,
        0.14,
        "hydration endpoint model",
        ha="center",
        color=DEPLOY,
        weight="bold",
        fontsize=6.8,
    )
    ax.text(
        0.31,
        0.06,
        "Stage 1: learn response coordinates        Stage 2: learn the endpoint",
        ha="center",
        color=INK,
        fontsize=6.3,
    )

    # Deployment side.
    ax.text(0.74, 0.80, "CC(=O)NC", ha="center", family="monospace", fontsize=7.0, color=INK)
    draw_molecule(ax, "CC(=O)NC", (0.82, 0.77), 0.095)
    arrow(ax, (0.79, 0.69), (0.79, 0.57), DEPLOY)
    response_strip(ax, 0.735, 0.48, 0.125, 0.075)
    ax.text(0.80, 0.445, "structure-predicted response", ha="center", fontsize=6.3, color=LEARNED)
    arrow(ax, (0.80, 0.425), (0.80, 0.32), DEPLOY)
    ax.text(0.80, 0.27, "SolvAI", ha="center", color=DEPLOY, weight="bold", fontsize=12)
    arrow(ax, (0.85, 0.27), (0.92, 0.27), DEPLOY)
    ax.text(
        0.96,
        0.27,
        r"$\Delta G_{\rm hyd}$",
        ha="center",
        va="center",
        color=DEPLOY,
        weight="bold",
        fontsize=9,
    )
    ax.text(0.82, 0.12, "SMILES only", ha="center", color=INK, fontsize=7, weight="bold")
    ax.text(0.82, 0.07, "no MD · no PIMD · no probe", ha="center", color=MID, fontsize=6.3)
    ax.plot([0.69, 0.96], [0.88, 0.88], color=PIMD, lw=0.8, ls=(0, (3, 2)))
    ax.text(
        0.825,
        0.895,
        "PIMD8: high-fidelity accuracy reference—not a retained teacher",
        ha="center",
        color=PIMD,
        fontsize=5.8,
    )
    save(fig, "fig1_concept")


def fig2_headline(
    metrics: dict, primary: pd.DataFrame, repeats: pd.DataFrame, paired: pd.DataFrame
) -> None:
    fig = plt.figure(figsize=(7.2, 5.35))
    grid = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.38)
    ax0, ax1, ax2, ax3 = [fig.add_subplot(grid[i, j]) for i in range(2) for j in range(2)]

    progression = [
        ("Structure\nonly", metrics["methods"]["matched_structure_only"]["mae_kcal_mol"]),
        ("+ compact\nresponse", metrics["methods"]["narrow_response"]["mae_kcal_mol"]),
        ("+ SMD\nwater", metrics["methods"]["narrow_plus_smd"]["mae_kcal_mol"]),
        ("+ conformer\nresponse", metrics["methods"]["full_solvai"]["mae_kcal_mol"]),
    ]
    x = np.arange(len(progression))
    y = np.array([value for _, value in progression])
    ax0.plot(x, y, color=LEARNED, lw=1.4, zorder=1)
    ax0.scatter(
        x, y, s=48, color=[MID, LEARNED, LEARNED, DEPLOY], edgecolor="white", lw=0.6, zorder=2
    )
    for position, value in zip(x, y, strict=True):
        ax0.text(
            position,
            value + 0.006,
            f"{value:.3f}",
            ha="center",
            fontsize=6.7,
            weight="bold" if position == 3 else "normal",
        )
    pimd = metrics["methods"]["arrow_pimd8"]["mae_kcal_mol"]
    ax0.axhline(pimd, color=PIMD, ls=(0, (3, 2)), lw=1)
    ax0.text(3.25, pimd + 0.002, "PIMD8  0.205", ha="right", color=PIMD, fontsize=6.3)
    ax0.set_xticks(x, [label for label, _ in progression])
    ax0.set_ylabel(r"OOF MAE (kcal mol$^{-1}$)")
    ax0.set_ylim(0.18, 0.325)
    ax0.set_title("Response priors close the accuracy gap", loc="left")
    clean(ax0)
    panel(ax0, "a")

    baseline = primary.loc[primary.method.eq("A_structure_only")].sort_values("molecule_id")
    full = primary.loc[primary.method.eq("F_full_solvai")].sort_values("molecule_id")
    lim = max(baseline.absolute_error.max(), full.absolute_error.max()) * 1.05
    improved = full.absolute_error.to_numpy() < baseline.absolute_error.to_numpy()
    ax1.plot([0, lim], [0, lim], color=MID, lw=0.8, ls="--")
    ax1.scatter(
        baseline.absolute_error,
        full.absolute_error,
        c=np.where(improved, DEPLOY, MID),
        s=17,
        edgecolor="white",
        lw=0.3,
        alpha=0.9,
    )
    ax1.fill_between([0, lim], [0, lim], [0, 0], color=DEPLOY_LIGHT, alpha=0.35)
    ax1.set(
        xlim=(0, lim),
        ylim=(0, lim),
        xlabel="Structure-only absolute error",
        ylabel="SolvAI absolute error",
    )
    ax1.set_aspect("equal", adjustable="box")
    ax1.text(
        0.97,
        0.06,
        f"{int(improved.sum())}/85 lower error",
        transform=ax1.transAxes,
        ha="right",
        color=DEPLOY,
        fontsize=6.5,
    )
    ax1.set_title("The gain is resolved molecule by molecule", loc="left")
    clean(ax1, grid=None)
    panel(ax1, "b")

    block_order = [
        ("Empirical +\ncorrected", "primary_B_empirical_residual"),
        ("Computation\ncore", "primary_C_computation_core"),
        ("SMD water", "primary_D_smd_water"),
        ("ConfSolv", "primary_E_confsolv"),
        ("Full\nSolvAI", "primary_F_full_solvai"),
    ]
    block = paired.set_index("analysis")
    for position, (label, key) in enumerate(block_order):
        row = block.loc[key]
        color = DEPLOY if key.endswith("F_full_solvai") else LEARNED
        ax2.errorbar(
            row.difference,
            position,
            xerr=[[row.difference - row.ci_low_95], [row.ci_high_95 - row.difference]],
            fmt="o",
            color=color,
            ecolor=color,
            capsize=2.5,
            ms=4.5,
            lw=1.0,
        )
    ax2.axvline(0, color=MID, lw=0.8)
    ax2.set_yticks(range(len(block_order)), [label for label, _ in block_order])
    ax2.invert_yaxis()
    ax2.set_xlabel(r"ΔMAE vs structure-only (kcal mol$^{-1}$)")
    ax2.set_title("Predeclared source blocks", loc="left")
    clean(ax2, grid="x")
    panel(ax2, "c")

    repeat_metrics = (
        repeats.groupby(["repeat", "method"], as_index=False)
        .absolute_error.mean()
        .rename(columns={"absolute_error": "mae"})
    )
    a = repeat_metrics.loc[repeat_metrics.method.eq("A_structure_only")].sort_values("repeat")
    f = repeat_metrics.loc[repeat_metrics.method.eq("F_full_solvai")].sort_values("repeat")
    for index in range(5):
        ax3.plot([0, 1], [a.mae.iloc[index], f.mae.iloc[index]], color=GRID, lw=0.8, zorder=1)
    ax3.scatter(np.zeros(5), a.mae, color=MID, s=27, edgecolor="white", lw=0.4, zorder=2)
    ax3.scatter(np.ones(5), f.mae, color=DEPLOY, s=27, edgecolor="white", lw=0.4, zorder=2)
    ax3.errorbar(
        [0, 1],
        [a.mae.mean(), f.mae.mean()],
        yerr=[a.mae.std(ddof=1), f.mae.std(ddof=1)],
        fmt="_",
        color=INK,
        ms=14,
        capsize=3,
        lw=1.0,
        zorder=3,
    )
    ax3.axhline(pimd, color=PIMD, ls=(0, (3, 2)), lw=1)
    ax3.set_xticks([0, 1], ["Structure\nonly", "SolvAI"])
    ax3.set_xlim(-0.35, 1.35)
    ax3.set_ylim(0.18, 0.33)
    ax3.set_ylabel(r"OOF MAE (kcal mol$^{-1}$)")
    ax3.text(1.32, pimd + 0.002, "PIMD8", ha="right", color=PIMD, fontsize=6.2)
    ax3.set_title("Five complete split repeats", loc="left")
    clean(ax3)
    panel(ax3, "d")
    save(fig, "fig2_headline")


def fig3_transfer(metrics: dict, separation: pd.DataFrame, zero: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7.2, 3.25))
    grid = fig.add_gridspec(1, 3, width_ratios=(1.05, 1.0, 1.0), wspace=0.45)
    ax0, ax1, ax2 = [fig.add_subplot(grid[0, index]) for index in range(3)]

    regimes = ["global_butina_0_70", "global_scaffold", "global_family"]
    labels = ["Molecular\nclusters", "Scaffolds", "Functional\nfamilies"]
    for y, regime in enumerate(regimes):
        values = separation.loc[separation.regime.eq(regime)].set_index("method").mae
        ax0.plot([values["F_full_solvai"], values["A_structure_only"]], [y, y], color=GRID, lw=2)
        ax0.scatter(
            values["A_structure_only"],
            y,
            color=MID,
            s=29,
            zorder=2,
            label="Structure only" if y == 0 else None,
        )
        ax0.scatter(
            values["F_full_solvai"],
            y,
            color=DEPLOY,
            s=29,
            zorder=2,
            label="SolvAI" if y == 0 else None,
        )
    ax0.set_yticks(range(3), labels)
    ax0.invert_yaxis()
    ax0.set_xlabel(r"MAE (kcal mol$^{-1}$)")
    ax0.set_title("Globally separated chemistry", loc="left")
    ax0.legend(frameon=False, fontsize=6.1, loc="lower right")
    clean(ax0, grid="x")
    panel(ax0, "a")

    thresholds = [0.5, 0.6, 0.7, 0.8]
    for method, label, color in (
        ("A_structure_only", "Structure only", MID),
        ("F_full_solvai", "SolvAI", DEPLOY),
    ):
        values = [
            separation.loc[
                separation.regime.eq(f"global_nn_{threshold:.2f}") & separation.method.eq(method),
                "mae",
            ].iloc[0]
            for threshold in thresholds
        ]
        ax1.plot(thresholds, values, marker="o", ms=4, lw=1.2, color=color, label=label)
    ax1.set_xlabel("Maximum allowed train–test similarity")
    ax1.set_ylabel(r"MAE (kcal mol$^{-1}$)")
    ax1.set_xticks(thresholds)
    ax1.set_ylim(0.18, 0.39)
    ax1.set_title("Nearest-neighbour exclusion", loc="left")
    ax1.legend(frameon=False, fontsize=6.1)
    clean(ax1)
    panel(ax1, "b")

    baseline = zero.loc[zero.method.eq("A_structure_only")].sort_values("molecule_id")
    full = zero.loc[zero.method.eq("F_full_solvai")].sort_values("molecule_id")
    differences = full.absolute_error.to_numpy() - baseline.absolute_error.to_numpy()
    ordered = np.sort(differences)
    ax2.axhline(0, color=MID, lw=0.8)
    ax2.scatter(
        np.arange(1, 86), ordered, s=12, c=np.where(ordered < 0, DEPLOY, MID), edgecolor="none"
    )
    ax2.axhline(differences.mean(), color=DEPLOY, lw=1.2, ls=(0, (3, 2)))
    ax2.text(
        84,
        differences.mean() - 0.025,
        f"mean {differences.mean():.3f}",
        ha="right",
        color=DEPLOY,
        fontsize=6.2,
    )
    ax2.set_xlabel("Molecules ranked by paired change")
    ax2.set_ylabel(r"SolvAI − structure-only error (kcal mol$^{-1}$)")
    ax2.set_title("No ARROW labels in training", loc="left")
    clean(ax2)
    panel(ax2, "c")
    save(fig, "fig3_transfer")


def fig4_frontier(metrics: dict) -> None:
    fig = plt.figure(figsize=(7.2, 3.35))
    grid = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.05, 1.0), wspace=0.62)
    ax0, ax1, ax2 = [fig.add_subplot(grid[0, index]) for index in range(3)]

    representation = pd.DataFrame(
        [
            ("Compact response\nsummaries", 0.1919131084),
            ("+ graph latent", 0.2121932664),
            ("+ feed-forward\nlatent", 0.2154807450),
        ],
        columns=["representation", "mae"],
    )
    ax0.barh(
        np.arange(3), representation.mae, color=[DEPLOY, LEARNED_LIGHT, LEARNED_LIGHT], height=0.56
    )
    ax0.set_yticks(np.arange(3), representation.representation)
    ax0.invert_yaxis()
    ax0.set_xlim(0.18, 0.225)
    ax0.set_xlabel(r"Exploratory OOF MAE (kcal mol$^{-1}$)")
    ax0.set_title("Compact coordinates transfer", loc="left", fontsize=7.4)
    for y, value in enumerate(representation.mae):
        ax0.text(value + 0.001, y, f"{value:.3f}", va="center", fontsize=6.2)
    clean(ax0, grid="x")
    panel(ax0, "a")

    response = pd.DataFrame(metrics["multilambda"]["response_head_metrics"])
    response = response.loc[~response.component.eq("lig__dhdl_pol_mean")].copy()
    response["component"] = response.component.map(
        {
            "lig_slv__dhdl_mean": "total",
            "lig_slv__dhdl_coul_mean": "electrostatic",
            "lig_slv__dhdl_vdw_mean": "van der Waals",
        }
    )
    pivot = response.pivot(index="component", columns="lambda", values="mae").loc[
        ["total", "electrostatic", "van der Waals"]
    ]
    image = ax1.imshow(pivot.to_numpy(), cmap="YlOrBr", vmin=0, vmax=5.5, aspect="auto")
    for row in range(pivot.shape[0]):
        for column in range(pivot.shape[1]):
            ax1.text(
                column,
                row,
                f"{pivot.iloc[row, column]:.2f}",
                ha="center",
                va="center",
                fontsize=6.2,
                color=INK,
            )
    ax1.set_xticks(range(3), [r"$\lambda=0.1$", r"$0.5$", r"$0.9$"])
    ax1.set_yticks(range(3), pivot.index)
    ax1.set_title("Response error", loc="left", fontsize=7.4)
    colorbar = fig.colorbar(image, ax=ax1, fraction=0.046, pad=0.03)
    colorbar.ax.set_title("MAE", fontsize=5.8, pad=2)
    colorbar.ax.tick_params(labelsize=5.8)
    panel(ax1, "b")

    lambda_metrics = metrics["multilambda"]["method_mae_kcal_mol"]
    downstream = [
        (
            "Base",
            lambda_metrics["Multi-lambda physics distillation A: structure/response baseline"],
        ),
        (
            "+ PIMD2",
            lambda_metrics[
                "Multi-lambda physics distillation B2: +distilled PIMD2 lambda response"
            ],
        ),
        (
            "+ hierarchy",
            lambda_metrics[
                "Multi-lambda physics distillation B1: +distilled classical-NQE-PIMD hierarchy"
            ],
        ),
        (
            "+ both",
            lambda_metrics[
                "Multi-lambda physics distillation B: +full distilled physics hierarchy"
            ],
        ),
    ]
    x = np.arange(len(downstream))
    values = [value for _, value in downstream]
    ax2.bar(x, values, color=[DEPLOY, PHYSICS_LIGHT, PHYSICS_LIGHT, NEGATIVE], width=0.62)
    ax2.set_xticks(x, [label for label, _ in downstream], rotation=24, ha="right")
    ax2.set_ylabel(r"OOF MAE (kcal mol$^{-1}$)")
    ax2.set_ylim(0.18, 0.23)
    for position, value in zip(x, values, strict=True):
        ax2.text(position, value + 0.0013, f"{value:.3f}", ha="center", fontsize=6.1)
    ax2.set_title("Endpoint effect", loc="left", fontsize=7.4)
    clean(ax2)
    panel(ax2, "c")
    save(fig, "fig4_frontier")


def ed_fig1_residuals(primary: pd.DataFrame) -> None:
    full = primary.loc[primary.method.eq("F_full_solvai")].sort_values("molecule_id")
    baseline = primary.loc[primary.method.eq("A_structure_only")].sort_values("molecule_id")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax0, ax1, ax2 = axes
    limits = [
        min(full.y_true.min(), full.y_pred.min()) - 0.4,
        max(full.y_true.max(), full.y_pred.max()) + 0.4,
    ]
    ax0.plot(limits, limits, color=MID, ls="--", lw=0.8)
    ax0.scatter(full.y_true, full.y_pred, s=15, color=DEPLOY, edgecolor="white", lw=0.3)
    ax0.set(
        xlabel=r"Experimental $\Delta G_{\rm hyd}$",
        ylabel=r"OOF prediction",
        xlim=limits,
        ylim=limits,
    )
    ax0.set_aspect("equal", adjustable="box")
    clean(ax0, grid=None)
    panel(ax0, "a")
    ax1.axhline(0, color=MID, lw=0.8)
    ax1.scatter(
        full.y_true, full.y_pred - full.y_true, s=15, color=DEPLOY, edgecolor="white", lw=0.3
    )
    ax1.set(xlabel=r"Experimental $\Delta G_{\rm hyd}$", ylabel="Prediction residual")
    clean(ax1)
    panel(ax1, "b")
    delta = full.absolute_error.to_numpy() - baseline.absolute_error.to_numpy()
    ax2.hist(delta, bins=np.linspace(delta.min(), delta.max(), 17), color=DEPLOY, alpha=0.85)
    ax2.axvline(0, color=MID, lw=0.8)
    ax2.axvline(delta.mean(), color=INK, lw=1.1, ls=(0, (3, 2)))
    ax2.set(xlabel="Paired absolute-error change", ylabel="Molecules")
    clean(ax2)
    panel(ax2, "c")
    save(fig, "ED_Fig1_residuals", extended=True)


def ed_fig2_provenance() -> None:
    source = pd.DataFrame(
        [
            ("CombiSolv-QM", "structures", 3961, 2, 3959),
            ("MolSolv", "SMD calculations", 350391, 32, 350359),
            ("ConfSolv", "usable connectivities", 17851, 22, 17829),
            ("Endpoint labels", "connectivities", 1280, 0, 1280),
        ],
        columns=["source", "unit", "before", "removed", "retained"],
    )
    fig, ax = plt.subplots(figsize=(7.2, 2.65))
    ax.axis("off")
    y_positions = np.linspace(0.78, 0.18, len(source))
    ax.text(0.03, 0.94, "Source", weight="bold")
    ax.text(0.31, 0.94, "Source-specific unit", weight="bold")
    ax.text(0.60, 0.94, "Standardized-equivalent exclusion", weight="bold")
    ax.text(0.90, 0.94, "Retained", weight="bold", ha="right")
    for row, y in zip(source.itertuples(), y_positions, strict=True):
        ax.text(0.03, y, row.source, weight="bold", color=INK)
        ax.text(0.31, y, row.unit, color=MID)
        ax.plot([0.57, 0.82], [y, y], color=GRID, lw=5, solid_capstyle="round")
        width = max(0.008, 0.25 * row.removed / max(row.before, 1))
        ax.plot([0.57, 0.57 + width], [y, y], color=NEGATIVE, lw=5, solid_capstyle="round")
        ax.text(
            0.695,
            y + 0.045,
            f"{row.removed:,} removed from {row.before:,}",
            ha="center",
            fontsize=6.2,
            color=NEGATIVE if row.removed else MID,
        )
        ax.text(0.90, y, f"{row.retained:,}", ha="right", color=DEPLOY, weight="bold")
    ax.text(
        0.03,
        0.05,
        "Source units are reported separately. Bar lengths are normalized within source,\nnot compared across calculations, structures and connectivities.",
        color=MID,
        fontsize=6.2,
    )
    save(fig, "ED_Fig2_provenance", extended=True)


def ed_fig3_alternatives(metrics: dict) -> None:
    alternatives = {
        "OpenFE diagnostics": metrics["alternative_supervision"]["openfe_diagnostics"],
        "MLFF hierarchy": metrics["alternative_supervision"]["mlff_hierarchy"],
        "DES370K response": metrics["alternative_supervision"]["des370k_water_response"],
        "ConfSolv graph latent": 0.2121932664,
        "ConfSolv FFN latent": 0.2154807450,
        "PIMD2 lambda response": metrics["multilambda"]["method_mae_kcal_mol"][
            "Multi-lambda physics distillation B2: +distilled PIMD2 lambda response"
        ],
        "Classical/NQE/PIMD": metrics["multilambda"]["method_mae_kcal_mol"][
            "Multi-lambda physics distillation B1: +distilled classical-NQE-PIMD hierarchy"
        ],
    }
    frame = pd.Series(alternatives).sort_values()
    fig, ax = plt.subplots(figsize=(6.2, 3.15))
    ax.scatter(frame.values, np.arange(len(frame)), color=LEARNED, s=28)
    ax.axvline(
        metrics["multilambda"]["method_mae_kcal_mol"][
            "Multi-lambda physics distillation A: structure/response baseline"
        ],
        color=DEPLOY,
        lw=1.0,
        ls=(0, (3, 2)),
        label="matched campaign base",
    )
    ax.set_yticks(np.arange(len(frame)), frame.index)
    ax.invert_yaxis()
    ax.set_xlabel(r"Exploratory OOF MAE (kcal mol$^{-1}$)")
    ax.legend(frameon=False, fontsize=6.2)
    clean(ax, grid="x")
    save(fig, "ED_Fig3_alternatives", extended=True)


def ed_fig4_selective() -> None:
    frontier = pd.read_csv(ROOT / "results/ablations/selective_compute_pareto.csv")
    frontier = frontier.loc[
        frontier.regime.eq("random_oof")
        & frontier.policy.isin(["Oracle selective PIMD", "Nested learned selective PIMD"])
    ]
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    for label, group in frontier.groupby("policy"):
        if "oracle" in label.lower():
            color, style = MID, "--"
        else:
            color, style = LEARNED, "-"
        group = group.sort_values("full_pimd_fraction")
        ax.plot(
            group["full_pimd_fraction"] * 100,
            group["mae"],
            marker="o",
            ms=3.5,
            lw=1.0,
            ls=style,
            color=color,
            label=label,
        )
    ax.set(xlabel="Full PIMD8 fallback (%)", ylabel=r"MAE (kcal mol$^{-1}$)")
    ax.set_title("Simulation-assisted reference (not SolvAI)", loc="left")
    ax.legend(frameon=False, fontsize=5.7)
    clean(ax)
    save(fig, "ED_Fig4_selective_pimd", extended=True)


def ed_fig5_lambda(metrics: dict) -> None:
    response = pd.DataFrame(metrics["multilambda"]["response_head_metrics"])
    response = response.loc[~response.component.eq("lig__dhdl_pol_mean")].copy()
    components = ["lig_slv__dhdl_mean", "lig_slv__dhdl_coul_mean", "lig_slv__dhdl_vdw_mean"]
    labels = ["total", "electrostatic", "van der Waals"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.75), gridspec_kw={"width_ratios": [1.3, 1.0]})
    for component, label, color in zip(components, labels, [INK, LEARNED, PHYSICS], strict=True):
        group = response.loc[response.component.eq(component)].sort_values("lambda")
        axes[0].plot(
            group["lambda"], group.mae, marker="o", ms=4, lw=1.15, color=color, label=label
        )
    axes[0].set(xlabel=r"Coupling coordinate $\lambda$", ylabel=r"Response MAE (kcal mol$^{-1}$)")
    axes[0].legend(frameon=False, fontsize=6.2)
    clean(axes[0])
    panel(axes[0], "a")
    values = metrics["multilambda"]["method_mae_kcal_mol"]
    names = ["base", "+ PIMD2", "+ hierarchy", "+ both", "integrated curve"]
    vals = [values[key] for key in values]
    axes[1].bar(np.arange(4), vals[:4], color=[DEPLOY, PHYSICS_LIGHT, PHYSICS_LIGHT, NEGATIVE])
    axes[1].set_xticks(np.arange(4), names[:4], rotation=25, ha="right")
    axes[1].set_ylabel(r"Endpoint OOF MAE (kcal mol$^{-1}$)")
    twin = axes[1].twinx()
    twin.scatter([3.8], [vals[4]], marker="D", color=NEGATIVE, s=28)
    twin.set_ylim(0, 1.7)
    twin.set_ylabel("Integrated-curve MAE", color=NEGATIVE)
    axes[1].text(
        3.8,
        0.07,
        f"{vals[4]:.2f}",
        transform=twin.get_xaxis_transform(),
        ha="center",
        color=NEGATIVE,
        fontsize=6.2,
    )
    clean(axes[1])
    panel(axes[1], "b")
    save(fig, "ED_Fig5_lambda_response", extended=True)


def ed_fig6_extrapolation(primary: pd.DataFrame, separation: pd.DataFrame) -> None:
    full = primary.loc[primary.method.eq("F_full_solvai")].copy()
    family_counts = (
        full.groupby("functional_group_family", sort=True)
        .size()
        .rename("count")
        .reset_index()
        .sort_values(
            ["count", "functional_group_family"],
            ascending=[False, True],
            kind="mergesort",
        )
    )
    families = family_counts["functional_group_family"].to_list()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), gridspec_kw={"width_ratios": [1.25, 1.0]})
    rng = np.random.default_rng(20260828)
    for position, family in enumerate(families):
        values = (
            full.loc[full.functional_group_family.eq(family)]
            .sort_values("molecule_id", kind="mergesort")["absolute_error"]
            .to_numpy()
        )
        jitter = rng.uniform(-0.12, 0.12, len(values))
        axes[0].scatter(values, position + jitter, s=12, color=DEPLOY, alpha=0.8)
        axes[0].plot(
            [values.mean(), values.mean()], [position - 0.18, position + 0.18], color=INK, lw=1
        )
    axes[0].set_yticks(
        np.arange(len(families)),
        [
            f"{family} (n={len(full.loc[full.functional_group_family.eq(family)])})"
            for family in families
        ],
    )
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Absolute OOF error")
    clean(axes[0], grid="x")
    panel(axes[0], "a")
    regimes = ["global_nn_0.70", "global_butina_0_70", "global_scaffold", "global_family"]
    labels = ["NN ≤ 0.70", "clusters", "scaffolds", "families"]
    label_offsets = [-10, 10, 0, 0]
    for position, (regime, label) in enumerate(zip(regimes, labels, strict=True)):
        group = separation.loc[separation.regime.eq(regime)].set_index("method")
        axes[1].plot(
            [0, 1],
            [group.loc["A_structure_only", "mae"], group.loc["F_full_solvai", "mae"]],
            color=GRID,
            lw=1,
        )
        axes[1].scatter([0], [group.loc["A_structure_only", "mae"]], color=MID, s=23)
        axes[1].scatter([1], [group.loc["F_full_solvai", "mae"]], color=DEPLOY, s=23)
        axes[1].annotate(
            label,
            (1, group.loc["F_full_solvai", "mae"]),
            xytext=(8, label_offsets[position]),
            textcoords="offset points",
            va="center",
            fontsize=6.1,
        )
    axes[1].set_xticks([0, 1], ["Structure only", "SolvAI"])
    axes[1].set_xlim(-0.25, 1.55)
    axes[1].set_ylabel(r"MAE (kcal mol$^{-1}$)")
    clean(axes[1])
    panel(axes[1], "b")
    save(fig, "ED_Fig6_extrapolation", extended=True)


def main() -> None:
    metrics = json.loads((ROOT / "results/paper_metrics.json").read_text())
    endpoint = pd.read_parquet(
        ROOT / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    )
    primary = endpoint.loc[endpoint.partition.eq("standardized_exclusion_primary")]
    repeats = endpoint.loc[endpoint.partition.eq("standardized_exclusion_repeat")]
    zero = endpoint.loc[endpoint.partition.eq("standardized_exclusion_zero_arrow")]
    paired = pd.read_csv(ROOT / "results/confirmatory/confirmatory_paired_comparisons.csv")
    separation = pd.read_csv(
        ROOT / "results/confirmatory/standardized_exclusion_global_separation_metrics.csv"
    )

    fig1_concept()
    fig2_headline(metrics, primary, repeats, paired)
    fig3_transfer(metrics, separation, zero)
    fig4_frontier(metrics)
    ed_fig1_residuals(primary)
    ed_fig2_provenance()
    ed_fig3_alternatives(metrics)
    ed_fig4_selective()
    ed_fig5_lambda(metrics)
    ed_fig6_extrapolation(primary, separation)

    obsolete = [
        "ED_Fig7_statistics",
    ]
    for stem in obsolete:
        for directory in (ED, PAPER_ED):
            for suffix in ("pdf", "svg", "png"):
                path = directory / f"{stem}.{suffix}"
                if path.exists():
                    path.unlink()
    print("Rendered four main and six Extended Data figures from frozen results.")


if __name__ == "__main__":
    main()
