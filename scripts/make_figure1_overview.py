#!/usr/bin/env python3
"""Build the full-page SolvAI Figure 1 as a reproducible hybrid SVG.

The scientific composition, labels, molecular structures and quantitative marks
are deterministic.  Curated, text-free ImageGen assets add visual depth; they are
embedded inside the SVG while text, connectors and data remain editable.
"""

from __future__ import annotations

import base64
import csv
import re
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
TARGETS = (ROOT / "figures/main", ROOT / "paper/figures/main")
ASSET_DIR = ROOT / "assets/figure1/generated"

# Color semantics follow the existing SolvAI figure system.
INK = "#17212B"
MID = "#66737D"
LIGHT = "#D9DEE2"
PALE = "#F6F8F9"
PHYSICS = "#C98222"
PHYSICS_LIGHT = "#F7E8D1"
LEARNED = "#2A78AD"
LEARNED_LIGHT = "#DDECF5"
DEPLOY = "#138B78"
PIMD = "#B33B75"
NEGATIVE = "#B44E46"
WHITE = "#FFFFFF"

# Match the paper's \usepackage{lmodern} typography exactly.  Register the TeX
# OpenType files explicitly because they are not exposed through macOS Font Book.
LM_FONT_DIR = Path("/usr/local/texlive/2025/texmf-dist/fonts/opentype/public/lm")
LM_FONTS = {
    "roman_regular": LM_FONT_DIR / "lmroman10-regular.otf",
    "roman_bold": LM_FONT_DIR / "lmroman10-bold.otf",
    "mono_regular": LM_FONT_DIR / "lmmono10-regular.otf",
}
for font_path in LM_FONT_DIR.glob("lmroman*.otf"):
    font_manager.fontManager.addfont(font_path)
for font_path in LM_FONT_DIR.glob("lmmono*.otf"):
    font_manager.fontManager.addfont(font_path)

mpl.rcParams.update(
    {
        "font.family": "Latin Modern Roman",
        "font.serif": ["Latin Modern Roman"],
        "font.monospace": ["Latin Modern Mono"],
        "font.size": 6.0,
        "axes.linewidth": 0.55,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.image_inline": True,
        "svg.hashsalt": "solvai-figure1-overview-v1",
        "savefig.dpi": 450,
    }
)


def rounded_box(
    ax,
    xy,
    width,
    height,
    *,
    facecolor=WHITE,
    edgecolor=LIGHT,
    linewidth=0.7,
    radius=0.02,
    zorder=1,
):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, *, color=INK, width=0.9, mutation=7, style="-|>", zorder=4):
    patch = FancyArrowPatch(
        start,
        end,
        transform=ax.transAxes,
        arrowstyle=style,
        mutation_scale=mutation,
        linewidth=width,
        color=color,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def lock_asset(parent, x, y, *, color=LEARNED, scale=1.0):
    """Place the Azure-generated lock as a clean, recolorable alpha mask."""

    path = ASSET_DIR / "lock_icon_flat_azure_v1.png"
    if not path.is_file():
        raise FileNotFoundError(path)
    pixels = mpimg.imread(path)
    # The generated source intentionally has generous whitespace.  Crop around
    # the icon, then discard the near-white background and retain antialiasing.
    pixels = pixels[307:645, 369:655]
    rgb = pixels[..., :3]
    luminance = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
    alpha = np.clip((0.97 - luminance) / 0.22, 0.0, 1.0)
    if pixels.shape[-1] == 4:
        alpha *= pixels[..., 3]
    tinted = np.empty((*pixels.shape[:2], 4), dtype=float)
    tinted[..., :3] = mpl.colors.to_rgb(color)
    tinted[..., 3] = alpha

    width = 0.029 * scale
    height = 0.090 * scale
    inset = parent.inset_axes(
        [x - width / 2, y - height / 2, width, height],
        zorder=8,
    )
    inset.imshow(tinted, interpolation="lanczos", aspect="equal")
    inset.axis("off")
    return inset


def panel_axes(fig, rect, label, title, *, facecolor=None, title_size=7.0):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    if facecolor:
        rounded_box(ax, (0.0, 0.0), 1.0, 1.0, facecolor=facecolor, edgecolor=LIGHT, radius=0.025)
    ax.text(0.01, 0.985, label, ha="left", va="top", fontsize=10.2, weight="bold", color=INK)
    ax.text(0.075, 0.973, title, ha="left", va="top", fontsize=title_size, weight="bold", color=INK)
    return ax


def draw_molecule(ax, smiles, *, line_color=INK):
    """Draw a chemically explicit 2D molecule into its own equal-aspect axes."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    molecule = Chem.Mol(molecule)
    try:
        Chem.Kekulize(molecule, clearAromaticFlags=True)
    except Chem.KekulizeException:
        pass
    AllChem.Compute2DCoords(molecule)
    coordinates = molecule.GetConformer().GetPositions()[:, :2]
    center = coordinates.mean(axis=0)
    coordinates = coordinates - center
    span = np.ptp(coordinates, axis=0)
    scale = max(float(span.max()), 1.0)
    coordinates /= scale

    for bond in molecule.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        p0, p1 = coordinates[i], coordinates[j]
        order = int(round(bond.GetBondTypeAsDouble()))
        direction = p1 - p0
        normal = np.array([-direction[1], direction[0]])
        norm = np.linalg.norm(normal)
        normal = normal / norm if norm else normal
        offsets = {1: [0.0], 2: [-0.018, 0.018], 3: [-0.027, 0.0, 0.027]}.get(order, [0.0])
        for offset in offsets:
            shifted = normal * offset
            ax.plot(
                [p0[0] + shifted[0], p1[0] + shifted[0]],
                [p0[1] + shifted[1], p1[1] + shifted[1]],
                color=line_color,
                linewidth=1.25 if order == 1 else 0.85,
                solid_capstyle="round",
                zorder=2,
            )

    atom_colors = {"O": NEGATIVE, "N": LEARNED, "S": PHYSICS, "Cl": DEPLOY, "F": DEPLOY}
    for atom, (x, y) in zip(molecule.GetAtoms(), coordinates, strict=True):
        symbol = atom.GetSymbol()
        if symbol != "C":
            label = symbol
            hydrogens = atom.GetTotalNumHs()
            if hydrogens:
                label += "H" if hydrogens == 1 else f"H{hydrogens}"
            ax.text(
                x,
                y,
                label,
                ha="center",
                va="center",
                fontsize=6.4,
                color=atom_colors.get(symbol, INK),
                weight="bold",
                bbox={"facecolor": WHITE, "edgecolor": "none", "pad": 0.4},
                zorder=3,
            )
    ax.set_xlim(-0.62, 0.62)
    ax.set_ylim(-0.52, 0.52)
    ax.set_aspect("equal")
    ax.axis("off")


def molecule_inset(parent, bounds, smiles):
    inset = parent.inset_axes(bounds)
    draw_molecule(inset, smiles)
    return inset


def image_asset(
    parent,
    bounds,
    filename,
    *,
    alpha=1.0,
    zorder=0,
    crop=None,
    preserve_aspect=False,
):
    """Place a curated generated bitmap inside the deterministic SVG layout."""

    path = ASSET_DIR / filename
    if not path.is_file():
        raise FileNotFoundError(path)
    pixels = mpimg.imread(path)
    if crop is not None:
        left, top, right, bottom = crop
        height, width = pixels.shape[:2]
        pixels = pixels[
            int(top * height) : int(bottom * height),
            int(left * width) : int(right * width),
        ]
    inset = parent.inset_axes(bounds, zorder=zorder)
    inset.imshow(
        pixels,
        interpolation="lanczos",
        alpha=alpha,
        aspect="equal" if preserve_aspect else "auto",
    )
    if preserve_aspect:
        inset.set_anchor("C")
    else:
        inset.set_aspect("auto")
    inset.axis("off")
    return inset


def draw_water(ax, x, y, angle=0.0, size=18):
    theta = np.deg2rad(angle)
    oxygen = np.array([x, y])
    raw = np.array([[-0.022, 0.028], [0.022, 0.028]])
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    hydrogens = oxygen + raw @ rotation.T
    for hydrogen in hydrogens:
        ax.plot(
            [oxygen[0], hydrogen[0]],
            [oxygen[1], hydrogen[1]],
            transform=ax.transAxes,
            color=MID,
            linewidth=0.55,
            zorder=2,
        )
    ax.scatter([oxygen[0]], [oxygen[1]], transform=ax.transAxes, s=size, c=NEGATIVE, edgecolors=WHITE, linewidths=0.4, zorder=3)
    ax.scatter(hydrogens[:, 0], hydrogens[:, 1], transform=ax.transAxes, s=size * 0.38, c=WHITE, edgecolors=MID, linewidths=0.35, zorder=3)


def response_vector(ax, x, y, width, height, *, compact=False):
    """Draw the exact six-family, 15-coordinate concatenated response vector."""
    groups = [
        (1, "#376B91"),
        (5, "#6FA3C5"),
        (1, "#D08B2D"),
        (1, "#A5AAB0"),
        (1, "#3E89B8"),
        (6, DEPLOY),
    ]
    total = sum(count for count, _ in groups)
    gap = width * (0.012 if compact else 0.017)
    cell_width = (width - gap * (total + len(groups) - 1)) / total
    cursor = x
    centers = []
    for group_index, (count, color) in enumerate(groups):
        start = cursor
        for _ in range(count):
            ax.add_patch(
                Rectangle(
                    (cursor, y),
                    cell_width,
                    height,
                    transform=ax.transAxes,
                    facecolor=color,
                    edgecolor=WHITE,
                    linewidth=0.45,
                    zorder=3,
                )
            )
            cursor += cell_width + gap
        end = cursor - gap
        centers.append((start + end) / 2)
        if group_index < len(groups) - 1:
            cursor += gap
    return centers


def response_group(ax, x, y, width, height, count, color):
    """Draw one family-owned response group without implying a shared output."""
    gap = width * 0.07
    cell_width = (width - gap * (count - 1)) / count
    for index in range(count):
        ax.add_patch(
            Rectangle(
                (x + index * (cell_width + gap), y),
                cell_width,
                height,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor=WHITE,
                linewidth=0.38,
                zorder=4,
            )
        )


def draw_source_icon(ax, kind, x, y):
    if kind == "continuum":
        ax.add_patch(Circle((x, y), 0.036, transform=ax.transAxes, fill=False, ec=PHYSICS, lw=0.9))
        ax.add_patch(Circle((x, y), 0.021, transform=ax.transAxes, fill=False, ec=PHYSICS, lw=0.7, ls=(0, (2, 1))))
    elif kind == "axes":
        ax.plot([x - 0.035, x - 0.035, x + 0.035], [y + 0.035, y - 0.035, y - 0.035], transform=ax.transAxes, color=PHYSICS, lw=0.8)
        for dx, dy in ((-0.018, 0.006), (0.0, 0.022), (0.020, -0.010), (0.030, 0.017), (-0.006, -0.016)):
            ax.scatter([x + dx], [y + dy], transform=ax.transAxes, s=5, color=PHYSICS)
    elif kind == "water":
        draw_water(ax, x, y, angle=180, size=13)
    elif kind == "shell":
        ax.add_patch(Circle((x, y), 0.034, transform=ax.transAxes, fc=PHYSICS_LIGHT, ec=PHYSICS, lw=0.7))
        ax.text(x, y, "ε", transform=ax.transAxes, ha="center", va="center", fontsize=6.5, color=PHYSICS, weight="bold")
    elif kind == "conformers":
        for dx, alpha in ((-0.018, 0.35), (0.0, 0.65), (0.018, 1.0)):
            ax.plot([x - 0.025 + dx, x + 0.025 + dx], [y - 0.012, y + 0.014], transform=ax.transAxes, color=PHYSICS, lw=0.9, alpha=alpha)
            ax.plot([x - 0.025 + dx, x + 0.010 + dx], [y - 0.012, y - 0.028], transform=ax.transAxes, color=PHYSICS, lw=0.9, alpha=alpha)


def source_tile(ax, x, y, width, title, subtitle, count, asset, *, subtitle_size=4.15):
    rounded_box(ax, (x, y), width, 0.19, facecolor=WHITE, edgecolor="#E2C18F", linewidth=0.75, radius=0.018)
    image_asset(
        ax,
        [x + 0.010, y + 0.022, 0.095, 0.145],
        asset,
        crop=(0.10, 0.10, 0.90, 0.90),
        zorder=2,
    )
    ax.text(
        x + 0.105,
        y + 0.137,
        title,
        transform=ax.transAxes,
        fontsize=4.85,
        weight="bold",
        color=INK,
        va="center",
        linespacing=0.94,
    )
    ax.text(
        x + 0.105,
        y + 0.087,
        subtitle,
        transform=ax.transAxes,
        fontsize=subtitle_size,
        color=MID,
        va="center",
        linespacing=0.94,
    )
    if count:
        ax.text(
            x + width - 0.014,
            y + 0.028,
            count,
            transform=ax.transAxes,
            fontsize=3.45,
            weight="bold",
            color=PHYSICS,
            ha="right",
            va="center",
            zorder=5,
        )


def fingerprint_bins(smiles, bins=22):
    molecule = Chem.MolFromSmiles(smiles)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    bits = np.fromiter((int(char) for char in generator.GetFingerprint(molecule).ToBitString()), dtype=float)
    return np.array([chunk.sum() for chunk in np.array_split(bits, bins)])


def draw_fingerprint(ax, x, y, width, height, smiles):
    values = fingerprint_bins(smiles)
    maximum = max(values.max(), 1.0)
    gap = width * 0.012
    cell_width = (width - gap * (len(values) - 1)) / len(values)
    for index, value in enumerate(values):
        alpha = 0.10 + 0.80 * value / maximum
        ax.add_patch(
            Rectangle(
                (x + index * (cell_width + gap), y),
                cell_width,
                height,
                transform=ax.transAxes,
                facecolor=LEARNED,
                alpha=alpha,
                edgecolor=LEARNED,
                linewidth=0.25,
            )
        )


def load_metrics():
    rows = {}
    with (ROOT / "paper/tables/main_comparison.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows[row["Method"]] = float(row["Fixed OOF MAE"])
    macro_text = (ROOT / "paper/tables/metrics_macros.tex").read_text(encoding="utf-8")
    shuffled_match = re.search(r"\\newcommand\{\\ShuffledMAE\}\{([0-9.]+)\}", macro_text)
    if not shuffled_match:
        raise RuntimeError("ShuffledMAE was not found in the frozen paper macros")
    return {
        "structure": rows["Matched structure-only"],
        "solvai": rows["SolvAI"],
        "pimd8": rows["ARROW/PIMD8"],
        "shuffled": float(shuffled_match.group(1)),
    }


def build_panel_a(fig, ax):
    image_asset(
        ax,
        [0.025, 0.34, 0.23, 0.43],
        "acetamide_3d_v2.png",
        crop=(0.10, 0.00, 0.90, 1.00),
        preserve_aspect=True,
    )
    ax.text(0.145, 0.79, "new solute", transform=ax.transAxes, ha="center", fontsize=5.0, weight="bold", color=INK)
    ax.text(0.145, 0.23, "gas", transform=ax.transAxes, ha="center", fontsize=5.2, color=MID)
    # Physical-computation arrows are orange and terminate in whitespace before assets.
    arrow(ax, (0.255, 0.54), (0.305, 0.54), color=PHYSICS, width=0.9)
    image_asset(ax, [0.31, 0.30, 0.29, 0.50], "acetamide_hydrated_v1.png")
    ax.text(0.46, 0.23, "water", transform=ax.transAxes, ha="center", fontsize=5.2, color=MID)

    rounded_box(ax, (0.61, 0.30), 0.36, 0.48, facecolor=WHITE, edgecolor=LIGHT, radius=0.018)
    image_asset(
        ax,
        [0.615, 0.305, 0.35, 0.43],
        "sampling_frames_v1.png",
        crop=(0.02, 0.14, 0.98, 0.90),
        zorder=2,
    )
    ax.text(0.79, 0.705, "repeat per new solute", transform=ax.transAxes, ha="center", fontsize=4.65, weight="bold", color=INK)
    ax.text(0.79, 0.25, "sampling + free-energy calculation", transform=ax.transAxes, ha="center", fontsize=4.15, color=MID)
    ax.text(0.50, 0.10, "accurate physical response, but costly per molecule", transform=ax.transAxes, ha="center", fontsize=5.65, weight="bold", color=PHYSICS)


def build_panel_b(ax):
    image_asset(ax, [0.13, 0.00, 0.74, 0.53], "response_manifold_v1.png", alpha=0.20, zorder=0)
    source_tile(ax, 0.03, 0.58, 0.28, "CombiSolv-QM", "COSMOtherm water", "1 coordinate", "cosmo_continuum_v1.png")
    source_tile(ax, 0.36, 0.58, 0.28, "Abraham", "empirical axes", "5 coordinates", "abraham_axes_v1.png")
    source_tile(ax, 0.69, 0.58, 0.28, "MolSolv", "SMD(water) continuum", "1 coordinate", "smd_continuum_v1.png", subtitle_size=3.75)
    source_tile(ax, 0.03, 0.28, 0.28, "OpenFF", "explicit-water correction", "1 corrected", "openff_explicit_v1.png", subtitle_size=3.70)
    source_tile(ax, 0.36, 0.28, 0.28, "GBn2", "implicit-water correction", "1 corrected", "gbn2_implicit_v1.png", subtitle_size=3.70)
    source_tile(
        ax,
        0.69,
        0.28,
        0.28,
        "ConfSolv",
        "conformer/water summaries",
        "6 coordinates",
        "confsolv_conformers_v1.png",
        subtitle_size=3.45,
    )
    ax.text(
        0.50,
        0.180,
        "six source families · 15 response coordinates",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=4.15,
        color=PHYSICS,
        weight="bold",
    )
    arrow(ax, (0.50, 0.154), (0.50, 0.095), color=PHYSICS, width=0.85, mutation=5)
    ax.text(0.50, 0.055, "complementary solvent-response information", transform=ax.transAxes, ha="center", fontsize=5.35, color=DEPLOY, weight="bold")
    ax.text(0.50, 0.020, "heterogeneous scales/approximations; no shared fidelity axis", transform=ax.transAxes, ha="center", fontsize=4.05, color=MID)


def build_panel_c(ax):
    ax.text(
        0.50,
        0.865,
        "STAGE 1 · ONE QUERY MOLECULE → SIX FROZEN SURROGATE FAMILIES",
        transform=ax.transAxes,
        ha="center",
        fontsize=4.45,
        weight="bold",
        color=LEARNED,
    )

    # The query is deliberately shown once.  Each branch therefore reads as a
    # prediction for the same neutral organic molecule, never as a different
    # molecule or a repeated source calculation.
    image_asset(
        ax,
        [0.012, 0.385, 0.150, 0.235],
        "panelc_query_acetophenone_v1.png",
        crop=(0.07, 0.07, 0.93, 0.93),
        zorder=2,
        preserve_aspect=True,
    )
    ax.text(0.088, 0.355, "one neutral query", transform=ax.transAxes, ha="center", fontsize=3.85, weight="bold", color=INK)
    ax.text(0.088, 0.315, "same molecule\nin every branch", transform=ax.transAxes, ha="center", va="top", fontsize=3.15, linespacing=0.92, color=MID)
    # Join the query to the collector with a plain line: the six branch arrows
    # carry the direction.  The collector stops exactly at the outer branches,
    # so no stray vertical stubs extend past the top or bottom arrow.
    arrow(ax, (0.158, 0.505), (0.190, 0.505), color=INK, width=0.72, style="-")
    ax.plot([0.190, 0.190], [0.275, 0.725], transform=ax.transAxes, color=INK, lw=0.7, zorder=3)

    rows = (
        (0.725, "model_icon_dmpnn_azure_v1.png", (0.245, 0.215, 0.830, 0.725), (0.225, 0.685, 0.072, 0.080), "CombiSolv-QM", "CHEMELEON · D-MPNN", 1, "#376B91"),
        (0.635, "model_icon_extratrees_azure_v1.png", (0.125, 0.335, 0.880, 0.600), (0.222, 0.602, 0.130, 0.066), "Abraham", "ExtraTrees", 5, "#6FA3C5"),
        (0.545, "model_icon_extratrees_azure_v1.png", (0.125, 0.335, 0.880, 0.600), (0.222, 0.512, 0.130, 0.066), "OpenFF corrected", "ExtraTrees", 1, "#D08B2D"),
        (0.455, "model_icon_extratrees_azure_v1.png", (0.125, 0.335, 0.880, 0.600), (0.222, 0.422, 0.130, 0.066), "GBn2 corrected", "ExtraTrees", 1, "#A5AAB0"),
        (0.365, "model_icon_dmpnn_azure_v1.png", (0.245, 0.215, 0.830, 0.725), (0.225, 0.325, 0.072, 0.080), "SMD(water)", "CHEMELEON · D-MPNN", 1, "#3E89B8"),
        (0.275, "model_icon_lightgbm_azure_v1.png", (0.085, 0.385, 0.925, 0.575), (0.220, 0.248, 0.145, 0.054), "ConfSolv", "LightGBM", 6, DEPLOY),
    )
    for center_y, asset, crop, bounds, family, model, count, color in rows:
        arrow(ax, (0.190, center_y), (0.216, center_y), color=INK, width=0.58, mutation=4.5)
        image_asset(ax, bounds, asset, crop=crop, zorder=2, preserve_aspect=True)
        ax.text(0.375, center_y + 0.014, family, transform=ax.transAxes, ha="left", va="center", fontsize=3.85, weight="bold", color=INK)
        ax.text(0.375, center_y - 0.020, model, transform=ax.transAxes, ha="left", va="center", fontsize=3.25, color=LEARNED)
        arrow(ax, (0.620, center_y), (0.658, center_y), color=color, width=0.72, mutation=4.5)
        group_width = 0.030 if count == 1 else 0.128
        response_group(ax, 0.665, center_y - 0.023, group_width, 0.046, count, color)
        ax.text(0.806, center_y, f"{count}", transform=ax.transAxes, ha="center", va="center", fontsize=3.7, weight="bold", color=color)
        ax.plot([0.820, 0.850], [center_y, center_y], transform=ax.transAxes, color=color, lw=0.62, zorder=3)

    # A single collector makes concatenation explicit; no branch appears to
    # generate the full vector by itself.
    ax.plot([0.850, 0.850], [0.122, 0.725], transform=ax.transAxes, color=DEPLOY, lw=0.82, zorder=3)
    ax.text(0.875, 0.500, "concatenate", transform=ax.transAxes, ha="center", va="center", rotation=90, fontsize=3.55, weight="bold", color=DEPLOY)

    ax.text(0.525, 0.195, "15 molecule-aligned response coordinates", transform=ax.transAxes, ha="center", fontsize=4.35, weight="bold", color=DEPLOY)
    ax.text(0.525, 0.157, "1 | 5 | 1 | 1 | 1 | 6 = 15", transform=ax.transAxes, ha="center", fontsize=4.15, weight="bold", color=LEARNED)
    response_vector(ax, 0.280, 0.105, 0.490, 0.035, compact=True)
    arrow(ax, (0.850, 0.122), (0.778, 0.122), color=DEPLOY, width=0.85, mutation=5)
    ax.text(0.50, 0.066, "predicted per molecule by frozen structure→response surrogates", transform=ax.transAxes, ha="center", fontsize=3.65, color=MID)
    ax.text(0.50, 0.035, "source calculations are not run at inference", transform=ax.transAxes, ha="center", fontsize=3.65, weight="bold", color=MID)
    ax.text(0.50, 0.009, "benchmark-equivalent molecules excluded from supervised response teachers", transform=ax.transAxes, ha="center", fontsize=3.25, color=MID)


def build_panel_d(ax):
    ax.text(0.50, 0.855, "STAGE 2 · RESPONSES + STRUCTURE → ENDPOINT", transform=ax.transAxes, ha="center", fontsize=4.75, weight="bold", color=DEPLOY)
    inputs = [
        (0.61, "response_channels_15_v1.png", "15 predicted response coordinates", LEARNED, (0.02, 0.20, 0.98, 0.82)),
        (0.38, "fingerprint_v1.png", "2,048-bit Morgan fingerprint", INK, (0.02, 0.28, 0.98, 0.74)),
        (0.15, "descriptors_v1.png", "217 RDKit 2D descriptors", INK, (0.03, 0.08, 0.97, 0.90)),
    ]
    for y, filename, label, color, crop in inputs:
        ax.text(0.16, y + 0.145, label, transform=ax.transAxes, ha="center", fontsize=4.15, color=color)
        image_asset(ax, [0.015, y, 0.30, 0.13], filename, crop=crop, zorder=2)
        # The input connectors meet at one explicit hub without arrowheads
        # colliding with each other or implying one forest per input block.
        arrow(ax, (0.315, y + 0.065), (0.365, 0.500), color=color, width=0.70, style="-")

    ax.scatter([0.365], [0.500], transform=ax.transAxes, s=17, color=DEPLOY, edgecolor=WHITE, linewidth=0.45, zorder=6)
    arrow(ax, (0.374, 0.500), (0.397, 0.500), color=DEPLOY, width=0.90, mutation=4.5)

    # A compact segmented vector replaces the former text-heavy box.  Segment
    # widths are qualitative so the small 15-coordinate block remains visible.
    rounded_box(ax, (0.397, 0.450), 0.148, 0.100, facecolor=WHITE, edgecolor=LIGHT, linewidth=0.65, radius=0.014, zorder=2)
    segments = (
        (0.404, 0.024, LEARNED, "15"),
        (0.430, 0.072, "#405E74", "2,048"),
        (0.504, 0.034, DEPLOY, "217"),
    )
    for x, width, color, label in segments:
        ax.add_patch(Rectangle((x, 0.463), width, 0.074, transform=ax.transAxes, facecolor=color, edgecolor=WHITE, linewidth=0.45, zorder=4))
        ax.text(x + width / 2, 0.500, label, transform=ax.transAxes, ha="center", va="center", fontsize=2.55, weight="bold", color=WHITE, zorder=5)
    ax.text(0.471, 0.590, "CONCATENATE", transform=ax.transAxes, ha="center", va="center", fontsize=3.35, weight="bold", color=DEPLOY)
    ax.text(0.471, 0.410, "2,280-D", transform=ax.transAxes, ha="center", va="center", fontsize=3.25, weight="bold", color=INK)
    ax.text(0.471, 0.375, "molecular representation", transform=ax.transAxes, ha="center", va="center", fontsize=2.85, color=INK)
    arrow(ax, (0.545, 0.500), (0.565, 0.500), color=DEPLOY, width=0.95, mutation=4.5)

    image_asset(
        ax,
        [0.565, 0.315, 0.255, 0.375],
        "endpoint_extratrees_forest_flat_v3.png",
        crop=(0.03, 0.04, 0.97, 0.96),
        zorder=2,
        preserve_aspect=True,
    )
    ax.text(0.692, 0.745, "3 × ExtraTrees ensembles", transform=ax.transAxes, ha="center", fontsize=4.65, weight="bold", color=DEPLOY)
    ax.text(0.692, 0.705, "360 trees each · seeds 11 / 29 / 47", transform=ax.transAxes, ha="center", fontsize=3.25, color=MID)
    ax.text(0.692, 0.055, "experimental hydration labels", transform=ax.transAxes, ha="center", fontsize=3.55, color=PHYSICS)
    image_asset(ax, [0.595, 0.090, 0.195, 0.065], "experimental_dg_strip_v1.png", crop=(0.03, 0.28, 0.97, 0.72), zorder=2)
    # Experimental labels supervise all three endpoint regressors.  A centered,
    # vertical arrow terminates at the ensemble's lower edge instead of pointing
    # ambiguously toward one internal tree glyph.
    arrow(ax, (0.692, 0.158), (0.692, 0.315), color=PHYSICS, width=0.75, mutation=5)
    arrow(ax, (0.820, 0.500), (0.865, 0.500), color=DEPLOY, width=1.0)
    ax.text(0.842, 0.435, "mean prediction", transform=ax.transAxes, fontsize=3.15, color=DEPLOY, ha="center")
    ax.text(0.910, 0.54, "ΔG", transform=ax.transAxes, fontsize=10.5, weight="bold", color=DEPLOY, ha="center")
    ax.text(0.962, 0.50, "hyd", transform=ax.transAxes, fontsize=4.35, color=DEPLOY, ha="left")
    ax.text(0.185, 0.055, "No PIMD-derived feature in the final stack", transform=ax.transAxes, ha="center", fontsize=3.55, color=MID)


def build_panel_e(ax):
    ax.text(0.09, 0.49, "SMILES input", transform=ax.transAxes, ha="center", fontsize=5.1, weight="bold", color=INK)
    ax.text(0.09, 0.36, "CC(=O)N", transform=ax.transAxes, ha="center", fontsize=5.9, family="monospace", color=INK)
    ax.text(0.09, 0.285, "acetamide", transform=ax.transAxes, ha="center", fontsize=3.7, color=MID)
    pipeline_y = 0.44
    arrow(ax, (0.145, pipeline_y), (0.187, pipeline_y), color=INK, width=0.8, mutation=5)
    image_asset(
        ax,
        [0.190, 0.270, 0.115, 0.34],
        "acetamide_3d_v2.png",
        crop=(0.10, 0.00, 0.90, 1.00),
        zorder=2,
        preserve_aspect=True,
    )
    arrow(ax, (0.305, pipeline_y), (0.332, pipeline_y), color=LEARNED, width=0.9, mutation=5)
    image_asset(
        ax,
        [0.335, 0.315, 0.15, 0.25],
        "frozen_surrogate_wedge_flat_v2.png",
        crop=(0.08, 0.10, 0.92, 0.90),
        zorder=2,
        preserve_aspect=True,
    )
    ax.text(0.410, 0.625, "six frozen surrogates", transform=ax.transAxes, ha="center", fontsize=4.0, weight="bold", color=LEARNED)
    ax.text(0.410, 0.590, "2 D-MPNN · 3 ExtraTrees · 1 LightGBM", transform=ax.transAxes, ha="center", fontsize=2.85, color=MID)
    lock_asset(ax, 0.470, 0.485, color=LEARNED, scale=0.90)

    arrow(ax, (0.485, pipeline_y), (0.512, pipeline_y), color=LEARNED, width=0.9, mutation=5)
    image_asset(ax, [0.515, 0.35, 0.14, 0.18], "response_channels_15_v1.png", crop=(0.02, 0.20, 0.98, 0.82), zorder=2)
    ax.text(0.585, 0.625, "15 response priors", transform=ax.transAxes, ha="center", fontsize=4.0, color=LEARNED, weight="bold")

    ax.text(0.585, 0.325, "fingerprint + descriptors", transform=ax.transAxes, ha="center", fontsize=3.35, color=INK, weight="bold")
    image_asset(ax, [0.515, 0.225, 0.14, 0.075], "fingerprint_v1.png", crop=(0.02, 0.28, 0.98, 0.74), zorder=2)
    image_asset(ax, [0.515, 0.125, 0.14, 0.085], "descriptors_v1.png", crop=(0.03, 0.10, 0.97, 0.90), zorder=2)
    arrow(ax, (0.305, 0.34), (0.502, 0.25), color=INK, width=0.7, mutation=5)

    # Both branches terminate at a visible '+' junction.  Their arrowheads point
    # to the node itself; a separate, longer green arrow then enters the endpoint
    # model, so the black connector cannot be mistaken for pointing into another
    # arrow or into empty space.
    merge_x = 0.690
    arrow(ax, (0.655, pipeline_y), (merge_x - 0.016, pipeline_y), color=DEPLOY, width=0.85, mutation=5)
    # Route the structure branch orthogonally into the bottom of the merge node.
    # Keeping the arrowhead vertical removes the former premature rightward jog.
    ax.plot([0.665, merge_x], [0.225, 0.225], transform=ax.transAxes, color=INK, lw=0.72, zorder=4)
    arrow(ax, (merge_x, 0.225), (merge_x, pipeline_y - 0.016), color=INK, width=0.72, mutation=4.5)
    # A scatter marker stays perfectly circular under the panel's non-square
    # transform, unlike an axes-coordinate Circle patch.
    ax.scatter([merge_x], [pipeline_y], transform=ax.transAxes, s=31, facecolor=WHITE, edgecolor=DEPLOY, linewidth=0.85, zorder=6)
    ax.text(merge_x, pipeline_y - 0.001, "+", transform=ax.transAxes, ha="center", va="center", fontsize=4.1, weight="bold", color=DEPLOY, zorder=7)
    ax.text(0.714, 0.390, "2,280-D", transform=ax.transAxes, ha="center", va="center", fontsize=2.85, weight="bold", color=DEPLOY)
    arrow(ax, (merge_x + 0.016, pipeline_y), (0.728, pipeline_y), color=DEPLOY, width=0.9, mutation=5)
    image_asset(
        ax,
        [0.728, 0.330, 0.140, 0.225],
        "endpoint_extratrees_forest_flat_v3.png",
        crop=(0.03, 0.04, 0.97, 0.96),
        zorder=2,
        preserve_aspect=True,
    )
    ax.text(0.825, 0.250, "3 × frozen ExtraTrees ensembles", transform=ax.transAxes, ha="center", fontsize=3.65, weight="bold", color=DEPLOY)
    ax.text(0.825, 0.220, "360 trees each · predictions averaged", transform=ax.transAxes, ha="center", fontsize=2.65, color=MID)
    lock_asset(ax, 0.850, 0.502, color=DEPLOY, scale=0.90)
    arrow(ax, (0.867, pipeline_y), (0.88, pipeline_y), color=DEPLOY, width=1.05, mutation=4.5)
    ax.text(0.910, pipeline_y, "ΔG", transform=ax.transAxes, fontsize=11.0, weight="bold", color=DEPLOY, ha="center", va="center")
    ax.text(0.953, 0.415, "hyd", transform=ax.transAxes, fontsize=5.0, color=DEPLOY, ha="left", va="center")
    ax.text(0.915, 0.31, "predicted hydration\nfree energy", transform=ax.transAxes, fontsize=3.45, color=MID, ha="center", linespacing=0.9)
    ax.text(0.50, 0.72, "learned response information is reused", transform=ax.transAxes, ha="center", fontsize=5.55, weight="bold", color=DEPLOY)
    ax.text(0.50, 0.055, "NO MD   ·   NO PIMD   ·   NO PROBE CALCULATION", transform=ax.transAxes, ha="center", fontsize=6.0, weight="bold", color=INK)


def build_panel_f(fig, ax, metrics):
    chart = ax.inset_axes([0.34, 0.13, 0.62, 0.51])
    labels = ["Matched structure-only", "Shuffled response priors", "SolvAI", "ARROW/PIMD8"]
    values = [metrics["structure"], metrics["shuffled"], metrics["solvai"], metrics["pimd8"]]
    colors = [MID, "#9AA3AA", DEPLOY, PIMD]
    y = np.arange(4)[::-1]
    chart.set_xlim(0, 0.35)
    chart.set_ylim(-0.55, 3.55)
    chart.set_yticks(y, labels)
    chart.set_xticks([0.0, 0.1, 0.2, 0.3])
    chart.tick_params(axis="both", labelsize=4.25, length=2, width=0.5, colors=INK)
    chart.spines[["top", "right", "left"]].set_visible(False)
    chart.spines["bottom"].set_color(MID)
    chart.spines["bottom"].set_linewidth(0.55)
    chart.grid(axis="x", color=LIGHT, linewidth=0.45, zorder=0)
    for yi, value, color in zip(y, values, colors, strict=True):
        chart.plot([0, value], [yi, yi], color=color, lw=1.4, solid_capstyle="round", zorder=2)
        chart.scatter([value], [yi], s=24, color=color, edgecolor=WHITE, linewidth=0.5, zorder=3)
        chart.text(value + 0.010, yi, f"{value:.3f}", va="center", fontsize=5.0, color=INK)
    chart.set_xlabel("MAE (kcal/mol) - lower is better", fontsize=5.0, color=INK, labelpad=3)
    chart.axhline(0.5, color=LIGHT, lw=0.65)
    ax.text(0.97, 0.875, "matched 85-solute ARROW chemistry", transform=ax.transAxes, ha="right", fontsize=4.7, color=MID)
    ax.text(0.97, 0.825, "response priors reduce matched MAE · shuffling abolishes the gain", transform=ax.transAxes, ha="right", fontsize=4.05, color=MID)
    ax.text(0.97, 0.775, "PIMD8-level accuracy on this reference chemistry", transform=ax.transAxes, ha="right", fontsize=4.35, color=PIMD, weight="bold")
    ax.text(0.97, 0.725, "accuracy reference; not a teacher · never used as a SolvAI input", transform=ax.transAxes, ha="right", fontsize=3.85, color=PIMD)


def embed_svg_fonts(path):
    """Keep SVG text editable while making Latin Modern rendering portable."""

    faces = []
    specifications = (
        ("Latin Modern Roman", 400, LM_FONTS["roman_regular"]),
        ("Latin Modern Roman", 700, LM_FONTS["roman_bold"]),
        ("Latin Modern Mono", 400, LM_FONTS["mono_regular"]),
    )
    for family, weight, font_path in specifications:
        if not font_path.is_file():
            raise FileNotFoundError(font_path)
        encoded = base64.b64encode(font_path.read_bytes()).decode("ascii")
        faces.append(
            "@font-face {"
            f"font-family:'{family}';font-style:normal;font-weight:{weight};"
            f"src:url(data:font/otf;base64,{encoded}) format('opentype');"
            "}"
        )
    svg = path.read_text(encoding="utf-8")
    font_style = '<style type="text/css"><![CDATA[' + "".join(faces) + "]]></style>"
    if "<defs>" not in svg:
        raise RuntimeError(f"SVG definitions block not found in {path}")
    svg = svg.replace("<defs>", "<defs>\n  " + font_style, 1)
    path.write_text(svg, encoding="utf-8")


def validate_panel_text_geometry(fig, panels):
    """Reject overlapping or clipped panel text before any artifact is saved."""

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    problems = []
    for panel in panels:
        panel_bounds = panel.get_window_extent(renderer)
        texts = [text for text in panel.texts if text.get_visible() and text.get_text().strip()]
        extents = [(text, text.get_window_extent(renderer)) for text in texts]
        for text, extent in extents:
            if (
                extent.x0 < panel_bounds.x0 - 0.5
                or extent.x1 > panel_bounds.x1 + 0.5
                or extent.y0 < panel_bounds.y0 - 0.5
                or extent.y1 > panel_bounds.y1 + 0.5
            ):
                overflow = (
                    max(panel_bounds.x0 - extent.x0, 0.0),
                    max(extent.x1 - panel_bounds.x1, 0.0),
                    max(panel_bounds.y0 - extent.y0, 0.0),
                    max(extent.y1 - panel_bounds.y1, 0.0),
                )
                problems.append(
                    f"out of bounds: {text.get_text()!r}; "
                    f"overflow L/R/B/T={tuple(round(value, 2) for value in overflow)} px"
                )
        for index, (left_text, left_extent) in enumerate(extents):
            for right_text, right_extent in extents[index + 1 :]:
                x_overlap = min(left_extent.x1, right_extent.x1) - max(left_extent.x0, right_extent.x0)
                y_overlap = min(left_extent.y1, right_extent.y1) - max(left_extent.y0, right_extent.y0)
                if x_overlap > 0.5 and y_overlap > 0.5:
                    problems.append(
                        f"overlap: {left_text.get_text()!r} with {right_text.get_text()!r}"
                    )
    if problems:
        raise RuntimeError("Panel text geometry validation failed:\n- " + "\n- ".join(problems))


def save_figure(fig, stem):
    primary = ROOT / "paper/figures/main"
    primary.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "pdf", "png"):
        path = primary / f"{stem}.{suffix}"
        metadata = {
            "svg": {"Date": None, "Creator": "SolvAI deterministic figure builder"},
            "pdf": {"CreationDate": None, "ModDate": None, "Creator": "SolvAI deterministic figure builder"},
            "png": {"Software": "SolvAI deterministic figure builder"},
        }[suffix]
        fig.savefig(path, dpi=450, facecolor=WHITE, metadata=metadata)
        if suffix == "svg":
            lines = path.read_text(encoding="utf-8").splitlines()
            path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")
            embed_svg_fonts(path)

    mirror = ROOT / "figures/main"
    mirror.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "pdf", "png"):
        shutil.copy2(primary / f"{stem}.{suffix}", mirror / f"{stem}.{suffix}")


def main():
    metrics = load_metrics()
    # 7.2 inches = 183 mm Nature double-column width; 5.9 inches = 150 mm.
    fig = plt.figure(figsize=(7.2, 5.9), facecolor=WHITE)

    panel_a = panel_axes(fig, [0.018, 0.700, 0.390, 0.282], "A", "Why physical response is expensive")
    panel_b = panel_axes(fig, [0.424, 0.700, 0.558, 0.282], "B", "Compute a complementary response vocabulary once")
    panel_c = panel_axes(fig, [0.018, 0.360, 0.462, 0.305], "C", "Distil response into structure-predictable coordinates", title_size=6.4)
    panel_d = panel_axes(fig, [0.502, 0.360, 0.480, 0.305], "D", "Learn the experimental endpoint separately")
    panel_e = panel_axes(fig, [0.018, 0.035, 0.622, 0.290], "E", "Simulation-free deployment")
    panel_f = panel_axes(fig, [0.660, 0.035, 0.322, 0.290], "F", "Matched evidence")

    build_panel_a(fig, panel_a)
    build_panel_b(panel_b)
    build_panel_c(panel_c)
    build_panel_d(panel_d)
    build_panel_e(panel_e)
    build_panel_f(fig, panel_f, metrics)

    validate_panel_text_geometry(
        fig,
        (panel_a, panel_b, panel_c, panel_d, panel_e, panel_f),
    )

    # Subtle separators establish reading order without enclosing each panel.
    overlay = fig.add_axes([0, 0, 1, 1], zorder=-1)
    overlay.axis("off")
    overlay.plot([0.018, 0.982], [0.682, 0.682], transform=overlay.transAxes, color=LIGHT, lw=0.65)
    overlay.plot([0.018, 0.982], [0.342, 0.342], transform=overlay.transAxes, color=LIGHT, lw=0.65)

    save_figure(fig, "fig1_concept")
    plt.close(fig)
    print("Rendered deterministic Figure 1 SVG/PDF/PNG at 183 × 150 mm.")


if __name__ == "__main__":
    main()
