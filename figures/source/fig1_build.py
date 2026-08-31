#!/usr/bin/env python3
"""Build the vector-first SolvAI concept figure from frozen scientific sources.

The script deliberately separates scientific source assets from final composition:
RDKit generates molecular depictions and conformers; Matplotlib renders the cavity
and the repository's measured lambda response; SVGWrite assembles the architecture;
CairoSVG exports publication PDF and the inspection PNG.
"""

from __future__ import annotations

import base64
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import svgwrite
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdFingerprintGenerator, rdMolAlign
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "figures/source/fig1_assets"
MAIN = ROOT / "figures/main"
PAPER_MAIN = ROOT / "paper/figures/main"

CANVAS_W = 1800
CANVAS_H = 1050
PNG_W = 3600

INK = "#19242E"
MID = "#66727D"
LIGHT = "#D9E0E5"
PALE = "#F5F7F8"
ORANGE = "#D78525"
ORANGE_LIGHT = "#F6E7D3"
BLUE = "#2678B8"
BLUE_DARK = "#195B8E"
BLUE_LIGHT = "#E4F0F8"
TEAL = "#118F7A"
TEAL_LIGHT = "#DDF2ED"
MAGENTA = "#B23A74"
RED = "#C65D4C"

FONT = "Arial, Liberation Sans, DejaVu Sans, sans-serif"
MONO = "DejaVu Sans Mono, Liberation Mono, monospace"

REFERENCE_SMILES = "CNC(C)=O"
REFERENCE_NAME = "N-methylacetamide"
CONFORMER_SMILES = "COCCOC"
CONFORMER_NAME = "1,2-dimethoxyethane"

PRIORS = [
    (1, "Calculated solvation", "COSMOtherm water", "combisolv_qm"),
    (9, "Calculated solvation", "SMD(water)", "smd_water"),
    (2, "Polarity / H-bonding", "Abraham E · excess molar refraction", "abraham_e"),
    (3, "Polarity / H-bonding", "Abraham S · dipolarity / polarizability", "abraham_s"),
    (4, "Polarity / H-bonding", "Abraham A · H-bond acidity", "abraham_a"),
    (5, "Polarity / H-bonding", "Abraham B · H-bond basicity", "abraham_b"),
    (6, "Polarity / H-bonding", "Abraham L · hexadecane–air partition", "abraham_l"),
    (7, "Explicit / implicit water", "corrected OpenFF explicit-water ΔG", "openff_corrected"),
    (8, "Explicit / implicit water", "corrected GBn2 implicit-solvent ΔG", "gbn2_corrected"),
    (10, "Conformational response", "gas conformer correction", "conf_gas_corr"),
    (11, "Conformational response", "solution conformer correction", "conf_solution_corr"),
    (12, "Conformational response", "hydration conformer correction", "conf_hydration_corr"),
    (13, "Conformational response", "conformer solvation-energy spread", "conf_gsolv_sd"),
    (14, "Conformational response", "mean conformer solvent response", "conf_response_mean"),
    (15, "Conformational response", "conformer response spread", "conf_response_sd"),
]

CONFIG_RESPONSE_FEATURES = [
    "combisolv_qm_teacher",
    "abraham_e_teacher",
    "abraham_s_teacher",
    "abraham_a_teacher",
    "abraham_b_teacher",
    "abraham_l_teacher",
    "openff_corrected_teacher",
    "gbn2_corrected_teacher",
    "molsolv_smd_teacher",
    "confsolv_gas_conformer_correction_teacher",
    "confsolv_solution_conformer_correction_teacher",
    "confsolv_hydration_conformer_correction_teacher",
    "confsolv_water_gsolv_std_teacher",
    "confsolv_water_response_mean_teacher",
    "confsolv_water_response_std_teacher",
]


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Liberation Sans",
            "font.size": 9.0,
            "axes.labelsize": 9.0,
            "axes.titlesize": 9.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "solvai-fig1-scientific-assets",
            "savefig.transparent": True,
        }
    )


def normalize_svg(path: Path) -> None:
    lines = path.read_text().splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n")


def normalize_pdf_metadata(path: Path) -> None:
    """Fix Cairo's timestamp without changing byte offsets in the PDF."""
    payload = path.read_bytes()
    normalized, count = re.subn(
        rb"/CreationDate \(D:\d{14}Z\)",
        b"/CreationDate (D:20000101000000Z)",
        payload,
    )
    if count != 1:
        raise AssertionError(f"Expected one Cairo CreationDate, found {count}")
    path.write_bytes(normalized)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def svg_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def molecule_svg(smiles: str, output: Path, width: int, height: int) -> tuple[Chem.Mol, dict[int, tuple[float, float]]]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Could not parse {smiles}")
    rdDepictor.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.08
    options.bondLineWidth = 2.1
    options.fixedBondLength = 32
    options.setAtomPalette(
        {
            6: (0.10, 0.14, 0.18),
            7: (0.12, 0.39, 0.66),
            8: (0.76, 0.25, 0.20),
            9: (0.10, 0.55, 0.42),
            16: (0.84, 0.54, 0.15),
            17: (0.10, 0.55, 0.42),
        }
    )
    drawer.DrawMolecule(molecule)
    coordinates = {
        index: (float(drawer.GetDrawCoords(index).x), float(drawer.GetDrawCoords(index).y))
        for index in range(molecule.GetNumAtoms())
    }
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    output.write_text(svg)
    normalize_svg(output)
    return molecule, coordinates


def indigo_molecule_svg(smiles: str, output: Path) -> None:
    """Create the renderer-comparison asset requested for art-direction review."""
    from indigo import Indigo
    from indigo.renderer import IndigoRenderer

    indigo = Indigo()
    renderer = IndigoRenderer(indigo)
    molecule = indigo.loadMolecule(smiles)
    molecule.layout()
    indigo.setOption("render-output-format", "svg")
    indigo.setOption("render-background-color", 1.0, 1.0, 1.0)
    indigo.setOption("render-coloring", True)
    indigo.setOption("render-bond-length", 38)
    indigo.setOption("render-margins", 12, 12)
    renderer.renderToFile(molecule, str(output))
    normalize_svg(output)


def water_shell_vector_svg(pdb_path: Path, output: Path) -> None:
    """Project a Packmol solute--water configuration into an editable SVG."""
    atoms: list[tuple[str, str, np.ndarray]] = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        residue = line[17:20].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        atoms.append((element, residue, xyz))
    solute = atoms[:12]
    waters = [atoms[index : index + 3] for index in range(12, len(atoms), 3)]
    solute_heavy = np.array([xyz for element, _, xyz in solute if element != "H"])
    ranked_waters = sorted(
        waters,
        key=lambda water: float(
            np.min(np.linalg.norm(solute_heavy - water[0][2][None, :], axis=1))
        ),
    )[:18]

    centre = solute_heavy.mean(axis=0)
    _, _, vh = np.linalg.svd(solute_heavy - centre, full_matrices=False)
    rotation = np.vstack([vh[0], vh[1], vh[2]])

    def project(xyz: np.ndarray) -> np.ndarray:
        return (xyz - centre) @ rotation.T

    projected_solute = [(element, project(xyz)) for element, _, xyz in solute]
    projected_waters = [[(element, project(xyz)) for element, _, xyz in water] for water in ranked_waters]
    all_xy = np.array([item[1][:2] for item in projected_solute] + [item[1][:2] for water in projected_waters for item in water])
    span = np.ptp(all_xy, axis=0)
    scale = min(380 / max(span[0], 1), 275 / max(span[1], 1))
    offset = np.array([220.0, 155.0]) - (all_xy.min(axis=0) + span / 2) * scale

    drawing = svgwrite.Drawing(output, size=(440, 310), viewBox="0 0 440 310")

    def xy(point: np.ndarray) -> tuple[float, float]:
        value = point[:2] * scale + offset
        return float(value[0]), float(310 - value[1])

    for water in sorted(projected_waters, key=lambda value: float(value[0][1][2])):
        oxygen = water[0][1]
        depth = float(np.clip((oxygen[2] + 6) / 12, 0, 1))
        opacity = 0.28 + 0.38 * depth
        o_xy = xy(oxygen)
        for _, hydrogen in water[1:]:
            drawing.add(
                drawing.line(
                    start=o_xy,
                    end=xy(hydrogen),
                    stroke="#8A969F",
                    stroke_width=1.1,
                    opacity=opacity,
                )
            )
        drawing.add(
            drawing.circle(
                center=o_xy,
                r=4.0 + 1.5 * depth,
                fill=RED,
                opacity=opacity,
            )
        )

    topology = Chem.AddHs(Chem.MolFromSmiles(REFERENCE_SMILES))
    for bond in topology.GetBonds():
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        start = np.array(xy(projected_solute[first][1]))
        end = np.array(xy(projected_solute[second][1]))
        vector = end - start
        normal = np.array([-vector[1], vector[0]])
        normal /= max(float(np.linalg.norm(normal)), 1.0)
        offsets = [0.0] if bond.GetBondTypeAsDouble() < 1.5 else [-2.2, 2.2]
        for displacement in offsets:
            delta = normal * displacement
            drawing.add(
                drawing.line(
                    start=tuple(start + delta),
                    end=tuple(end + delta),
                    stroke=INK,
                    stroke_width=3.1 if len(offsets) == 2 else 3.8,
                )
            )
    atom_colors = {"C": INK, "N": BLUE, "O": RED, "H": "#AEB8BF"}
    for element, point in projected_solute:
        drawing.add(
            drawing.circle(
                center=xy(point),
                r={"H": 3.2, "C": 7.0, "N": 7.6, "O": 7.8}.get(element, 6.5),
                fill=atom_colors.get(element, MID),
                stroke="white",
                stroke_width=1.0,
            )
        )
    drawing.save(pretty=True)
    normalize_svg(output)


def fingerprint_svg(smiles: str, output: Path) -> None:
    """Render the actual 2,048-bit Morgan fingerprint as a sparse vector matrix."""
    molecule = Chem.MolFromSmiles(smiles)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    vector = np.asarray(generator.GetFingerprintAsNumPy(molecule), dtype=np.uint8)
    matrix = vector.reshape(32, 64)
    drawing = svgwrite.Drawing(output, size=(256, 128), viewBox="0 0 256 128")
    for row in range(32):
        for column in range(64):
            if matrix[row, column]:
                drawing.add(
                    drawing.rect(
                        insert=(column * 4 + 0.6, row * 4 + 0.6),
                        size=(2.8, 2.8),
                        fill=INK,
                    )
                )
    drawing.save(pretty=True)
    normalize_svg(output)


def crop_rdkit_svg(source: Path, output: Path, coordinates: dict[int, tuple[float, float]]) -> None:
    """Crop an RDKit-authored SVG while leaving every molecular path untouched."""
    pad = 30.0
    xs = [point[0] for point in coordinates.values()]
    ys = [point[1] for point in coordinates.values()]
    x0, y0 = min(xs) - pad, min(ys) - pad
    width, height = max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad
    svg = source.read_text()
    svg = re.sub(
        r"width='[^']+' height='[^']+' viewBox='[^']+'",
        f"width='{width:.1f}px' height='{height:.1f}px' viewBox='0 0 {width:.1f} {height:.1f}'",
        svg,
        count=1,
    )
    svg = svg.replace(
        "<!-- END OF HEADER -->",
        f"<!-- END OF HEADER -->\n<g transform='translate({-x0:.1f},{-y0:.1f})'>",
        1,
    )
    svg = svg.rsplit("</svg>", 1)[0] + "</g>\n</svg>\n"
    output.write_text(svg)
    normalize_svg(output)


def cavity_svg(
    molecule: Chem.Mol,
    coordinates: dict[int, tuple[float, float]],
    molecule_asset: Path,
    output: Path,
    width: int = 440,
    height: int = 250,
) -> None:
    """Draw a deterministic 2D van der Waals envelope around the RDKit depiction."""
    table = Chem.GetPeriodicTable()
    bond_lengths = []
    for bond in molecule.GetBonds():
        x1, y1 = coordinates[bond.GetBeginAtomIdx()]
        x2, y2 = coordinates[bond.GetEndAtomIdx()]
        bond_lengths.append(math.hypot(x2 - x1, y2 - y1))
    pixels_per_angstrom = float(np.median(bond_lengths)) / 1.5
    xs = np.linspace(0, width, 360)
    ys = np.linspace(0, height, 220)
    xx, yy = np.meshgrid(xs, ys)
    distance_ratio = np.full_like(xx, np.inf, dtype=float)
    for atom in molecule.GetAtoms():
        x, y = coordinates[atom.GetIdx()]
        radius = table.GetRvdw(atom.GetAtomicNum()) * pixels_per_angstrom
        ratio = np.sqrt((xx - x) ** 2 + (yy - y) ** 2) / radius
        distance_ratio = np.minimum(distance_ratio, ratio)

    fig, ax = plt.subplots(figsize=(4.4, 2.5), dpi=100)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.contourf(
        xx,
        yy,
        distance_ratio,
        levels=[0.0, 1.0, 1.30, 1.62],
        colors=["#F9EDD9", "#F8E4C7", "#FBEFDF"],
        alpha=0.95,
    )
    ax.contour(
        xx,
        yy,
        distance_ratio,
        levels=[1.0, 1.30, 1.62],
        colors=[ORANGE, "#E4A75E", "#EDC994"],
        linewidths=[1.4, 1.0, 0.8],
    )
    ax.set(xlim=(0, width), ylim=(height, 0))
    ax.set_aspect("equal")
    ax.axis("off")
    base = ASSETS / "continuum_cavity_field.svg"
    fig.savefig(
        base,
        format="svg",
        bbox_inches=None,
        pad_inches=0,
        metadata={"Date": None, "Creator": "SolvAI Figure 1 build"},
    )
    plt.close(fig)
    normalize_svg(base)

    drawing = svgwrite.Drawing(output, size=(width, height), viewBox=f"0 0 {width} {height}")
    drawing.add(drawing.image(href=svg_uri(base), insert=(0, 0), size=(width, height)))
    drawing.add(
        drawing.image(href=svg_uri(molecule_asset), insert=(0, 0), size=(width, height))
    )
    drawing.save(pretty=True)
    normalize_svg(output)


def abraham_axes_svg(output: Path, width: int = 420, height: int = 250) -> None:
    drawing = svgwrite.Drawing(output, size=(width, height), viewBox=f"0 0 {width} {height}")
    centre = np.array([width / 2, height / 2 + 5])
    radius = 78.0
    labels = ["E", "S", "A", "B", "L"]
    angles = np.deg2rad(np.array([-90, -18, 54, 126, 198]))
    points = [centre + radius * np.array([math.cos(a), math.sin(a)]) for a in angles]
    drawing.add(
        drawing.polygon(
            points=[tuple(point) for point in points],
            fill="none",
            stroke="#BCD7EA",
            stroke_width=2,
        )
    )
    drawing.add(drawing.circle(center=tuple(centre), r=23, fill=BLUE_DARK))
    drawing.add(
        drawing.text(
            "solute",
            insert=(centre[0], centre[1] + 5),
            text_anchor="middle",
            font_family=FONT,
            font_size=17,
            fill="white",
            font_weight="bold",
        )
    )
    for point, symbol in zip(points, labels, strict=True):
        drawing.add(
            drawing.line(
                start=tuple(centre), end=tuple(point), stroke="#80B1D3", stroke_width=1.3
            )
        )
        drawing.add(drawing.circle(center=tuple(point), r=24, fill=BLUE_LIGHT, stroke=BLUE, stroke_width=1.5))
        drawing.add(
            drawing.text(
                symbol,
                insert=(point[0], point[1] + 8),
                text_anchor="middle",
                font_family=FONT,
                font_size=24,
                fill=BLUE_DARK,
                font_weight="bold",
            )
        )
    drawing.save(pretty=True)
    normalize_svg(output)


def select_conformers(molecule: Chem.Mol, count: int = 3) -> tuple[list[int], dict[int, float]]:
    params = AllChem.ETKDGv3()
    params.randomSeed = 0x5A17
    params.numThreads = 1
    params.pruneRmsThresh = 0.15
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=30, params=params))
    if not conformer_ids:
        raise RuntimeError("ETKDG did not generate conformers")
    optimisation = AllChem.MMFFOptimizeMoleculeConfs(molecule, numThreads=1, maxIters=1000)
    energies = {cid: float(optimisation[index][1]) for index, cid in enumerate(conformer_ids)}
    ranked = sorted(conformer_ids, key=energies.get)
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    selected: list[int] = []
    for cid in ranked:
        if not selected:
            selected.append(cid)
        else:
            distinct = all(
                rdMolAlign.GetBestRMS(
                    molecule,
                    molecule,
                    prbId=cid,
                    refId=other,
                    map=[[(index, index) for index in heavy]],
                )
                >= 0.42
                for other in selected
            )
            if distinct:
                selected.append(cid)
        if len(selected) == count:
            break
    if len(selected) < count:
        selected = ranked[:count]
    return selected, energies


def conformer_svg(output: Path, sdf_output: Path, metadata_output: Path) -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles(CONFORMER_SMILES))
    selected, energies = select_conformers(molecule)
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    reference = selected[0]
    atom_map = [(index, index) for index in heavy]
    for cid in selected[1:]:
        rdMolAlign.AlignMol(molecule, molecule, prbCid=cid, refCid=reference, atomMap=atom_map)

    writer = Chem.SDWriter(str(sdf_output))
    minimum = min(energies[cid] for cid in selected)
    for order, cid in enumerate(selected, start=1):
        molecule.SetProp("_Name", f"{CONFORMER_NAME} conformer {order}")
        molecule.SetProp("conformer_id", str(cid))
        molecule.SetProp("MMFF_relative_energy_kcal_mol", f"{energies[cid] - minimum:.6f}")
        writer.write(molecule, confId=cid)
    writer.close()

    fig, axes = plt.subplots(1, 3, figsize=(4.5, 2.2), constrained_layout=True)
    atom_colors = {6: INK, 8: RED}
    records = []
    for axis, cid in zip(axes, selected, strict=True):
        conformer = molecule.GetConformer(cid)
        coords = np.array([conformer.GetAtomPosition(index) for index in heavy], dtype=float)
        coords -= coords.mean(axis=0)
        _, _, vh = np.linalg.svd(coords, full_matrices=False)
        xy = coords @ vh[:2].T
        mapping = {atom_index: xy[position] for position, atom_index in enumerate(heavy)}
        for bond in molecule.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i in mapping and j in mapping:
                axis.plot(
                    [mapping[i][0], mapping[j][0]],
                    [mapping[i][1], mapping[j][1]],
                    color=INK,
                    lw=1.4,
                    solid_capstyle="round",
                    zorder=1,
                )
        for atom_index in heavy:
            atom = molecule.GetAtomWithIdx(atom_index)
            x, y = mapping[atom_index]
            axis.scatter(
                [x],
                [y],
                s=42 if atom.GetAtomicNum() == 8 else 25,
                color=atom_colors.get(atom.GetAtomicNum(), MID),
                edgecolor="white",
                linewidth=0.5,
                zorder=2,
            )
            if atom.GetAtomicNum() == 8:
                axis.text(x, y, "O", ha="center", va="center", color="white", fontsize=6.5, weight="bold")
        relative = energies[cid] - minimum
        axis.text(
            0.08,
            0.90,
            str(len(records) + 1),
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=9.0,
            weight="bold",
            color=BLUE_DARK,
        )
        axis.set_aspect("equal")
        axis.axis("off")
        margin = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1]))) * 0.65
        axis.set(xlim=(-margin * 0.88, margin * 0.88), ylim=(-margin * 0.88, margin * 0.88))
        records.append({"conformer_id": cid, "relative_mmff_energy_kcal_mol": relative})
    fig.savefig(
        output,
        format="svg",
        bbox_inches="tight",
        pad_inches=0.02,
        metadata={"Date": None, "Creator": "SolvAI Figure 1 build"},
    )
    plt.close(fig)
    normalize_svg(output)
    metadata_output.write_text(
        json.dumps(
            {
                "molecule": CONFORMER_NAME,
                "smiles": CONFORMER_SMILES,
                "method": "RDKit ETKDGv3; MMFF optimisation; fixed seed 0x5A17",
                "selected": records,
            },
            indent=2,
        )
        + "\n"
    )


def lambda_response_svg(output: Path, data_output: Path) -> None:
    source = ROOT / "results/ablations/pimd2_multilambda_teacher.parquet"
    table = pd.read_parquet(source)
    row = table.loc[table.molecule_name.str.casefold().eq(REFERENCE_NAME.casefold())]
    if len(row) != 1:
        raise AssertionError(f"Expected one {REFERENCE_NAME} lambda record, found {len(row)}")
    row = row.iloc[0]
    lambdas = np.array([0.1, 0.5, 0.9])
    component_columns = {
        "total": "lig_slv__dhdl_mean",
        "Coulomb": "lig_slv__dhdl_coul_mean",
        "van der Waals": "lig_slv__dhdl_vdw_mean",
    }
    colors = {"total": ORANGE, "Coulomb": BLUE, "van der Waals": MID}
    records = []
    fig, ax = plt.subplots(figsize=(3.6, 2.15))
    for label, prefix in component_columns.items():
        values = np.array([float(row[f"{prefix}__lambda{suffix}"]) for suffix in ("01", "05", "09")])
        ax.plot(
            lambdas,
            values,
            marker="o",
            ms=3.8,
            lw=1.25,
            color=colors[label],
            label=label,
        )
        records.extend(
            {"molecule": REFERENCE_NAME, "lambda": float(lam), "component": label, "dhdl": float(value)}
            for lam, value in zip(lambdas, values, strict=True)
        )
    ax.axhline(0, color=LIGHT, lw=0.7, zorder=0)
    ax.set(
        xlim=(0.05, 0.95),
        xticks=lambdas,
        xlabel=r"coupling coordinate $\lambda$",
        ylabel=r"$\langle\mathrm{d}H/\mathrm{d}\lambda\rangle$ (kcal mol$^{-1}$)",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, color=MID)
    ax.legend(frameon=False, fontsize=7.8, ncol=3, loc="upper right", handlelength=1.4, columnspacing=0.8)
    fig.savefig(
        output,
        format="svg",
        bbox_inches="tight",
        pad_inches=0.03,
        metadata={"Date": None, "Creator": "SolvAI Figure 1 build"},
    )
    plt.close(fig)
    normalize_svg(output)
    pd.DataFrame(records).to_csv(data_output, index=False)


def math_label_svg(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(2.5, 0.9))
    ax.text(
        0.5,
        0.5,
        r"$\Delta G_{\mathrm{hyd}}$",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=34,
        color=TEAL,
        weight="bold",
    )
    ax.axis("off")
    fig.savefig(
        output,
        format="svg",
        bbox_inches="tight",
        pad_inches=0,
        metadata={"Date": None, "Creator": "SolvAI Figure 1 build"},
    )
    plt.close(fig)
    normalize_svg(output)


def response_priors_svg(output: Path) -> None:
    drawing = svgwrite.Drawing(output, size=(800, 560), viewBox="0 0 800 560")
    placements = [
        (
            18,
            20,
            365,
            "Calculated solvation",
            [item for item in PRIORS if item[1] == "Calculated solvation"],
            "#DCECF7",
        ),
        (
            417,
            20,
            365,
            "Polarity / H-bonding",
            [item for item in PRIORS if item[1] == "Polarity / H-bonding"],
            "#E8F2F8",
        ),
        (
            18,
            284,
            365,
            "Explicit / implicit water",
            [item for item in PRIORS if item[1] == "Explicit / implicit water"],
            "#D8EAF6",
        ),
        (
            417,
            284,
            365,
            "Conformational response",
            [item for item in PRIORS if item[1] == "Conformational response"],
            "#E3F0F6",
        ),
    ]
    for x, y, width, title, items, fill in placements:
        header_height = 36
        drawing.add(drawing.rect(insert=(x, y), size=(width, header_height), fill=fill))
        drawing.add(
            drawing.text(
                f"{title}  ({len(items)})",
                insert=(x + 12, y + 24),
                font_family=FONT,
                font_size=19,
                font_weight="bold",
                fill=BLUE_DARK,
            )
        )
        row_height = 36 if len(items) >= 5 else 47
        for index, (number, _, label, _) in enumerate(items):
            ry = y + header_height + 12 + index * row_height
            drawing.add(drawing.circle(center=(x + 19, ry + 12), r=10, fill=BLUE))
            drawing.add(
                drawing.text(
                    str(number),
                    insert=(x + 19, ry + 17),
                    text_anchor="middle",
                    font_family=FONT,
                    font_size=12,
                    font_weight="bold",
                    fill="white",
                )
            )
            drawing.add(
                drawing.text(
                    label,
                    insert=(x + 38, ry + 17),
                    font_family=FONT,
                    font_size=16 if len(label) > 34 else 18,
                    fill=INK,
                )
            )
    drawing.add(
        drawing.text(
            "15 named scalar targets — not latent embeddings",
            insert=(400, 550),
            text_anchor="middle",
            font_family=FONT,
            font_size=18,
            font_weight="bold",
            fill=BLUE,
        )
    )
    drawing.save(pretty=True)
    normalize_svg(output)


def add_text(
    drawing: svgwrite.Drawing,
    text: str,
    x: float,
    y: float,
    *,
    size: float = 24,
    color: str = INK,
    weight: str = "normal",
    anchor: str = "start",
    family: str = FONT,
    italic: bool = False,
) -> None:
    drawing.add(
        drawing.text(
            text,
            insert=(x, y),
            text_anchor=anchor,
            font_family=family,
            font_size=size,
            font_weight=weight,
            font_style="italic" if italic else "normal",
            fill=color,
        )
    )


def add_multiline(
    drawing: svgwrite.Drawing,
    lines: list[str],
    x: float,
    y: float,
    *,
    size: float = 22,
    leading: float = 1.25,
    color: str = INK,
    weight: str = "normal",
    anchor: str = "start",
    family: str = FONT,
) -> None:
    text = drawing.text(
        "",
        insert=(x, y),
        text_anchor=anchor,
        font_family=family,
        font_size=size,
        font_weight=weight,
        fill=color,
    )
    for index, line in enumerate(lines):
        text.add(drawing.tspan(line, x=[x], dy=[0 if index == 0 else size * leading]))
    drawing.add(text)


def markers(drawing: svgwrite.Drawing) -> dict[str, svgwrite.container.Marker]:
    result = {}
    for name, color in (("orange", ORANGE), ("blue", BLUE), ("teal", TEAL), ("ink", INK)):
        marker = drawing.marker(
            insert=(12, 6),
            size=(12, 12),
            orient="auto",
            markerUnits="userSpaceOnUse",
            id=f"arrow-{name}",
        )
        marker.add(drawing.path(d="M 0 0 L 12 6 L 0 12 z", fill=color))
        drawing.defs.add(marker)
        result[name] = marker
    return result


def arrow(
    drawing: svgwrite.Drawing,
    marker: svgwrite.container.Marker,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    width: float = 4,
    dash: str | None = None,
) -> None:
    line = drawing.line(start=start, end=end, stroke=color, stroke_width=width)
    line["marker-end"] = marker.get_funciri()
    if dash:
        line["stroke-dasharray"] = dash
    drawing.add(line)


def embed(
    drawing: svgwrite.Drawing,
    path: Path,
    x: float,
    y: float,
    width: float,
    height: float,
) -> None:
    drawing.add(
        drawing.image(
            href=svg_uri(path),
            insert=(x, y),
            size=(width, height),
            preserveAspectRatio="xMidYMid meet",
        )
    )


def source_card(
    drawing: svgwrite.Drawing,
    x: float,
    y: float,
    width: float,
    title: str,
    subtitle: str,
) -> None:
    add_text(drawing, title, x, y, size=22, color=INK, weight="bold")
    add_text(drawing, subtitle, x, y + 25, size=16, color=MID)
    drawing.add(
        drawing.line(
            start=(x, y + 36),
            end=(x + width, y + 36),
            stroke=ORANGE,
            stroke_width=2,
        )
    )


def final_figure(assets: dict[str, Path], output: Path) -> None:
    drawing = svgwrite.Drawing(
        output,
        size=("180mm", "105mm"),
        viewBox=f"0 0 {CANVAS_W} {CANVAS_H}",
        profile="full",
    )
    drawing.add(drawing.rect(insert=(0, 0), size=(CANVAS_W, CANVAS_H), fill="white"))
    arrowheads = markers(drawing)

    # Panel boundaries and top-level hierarchy.
    drawing.add(drawing.line(start=(560, 58), end=(560, 1015), stroke=LIGHT, stroke_width=2))
    drawing.add(drawing.line(start=(1250, 58), end=(1250, 1015), stroke=LIGHT, stroke_width=2))
    add_text(drawing, "a", 28, 38, size=30, weight="bold")
    add_text(drawing, "SOLVATION-RESPONSE SOURCES", 70, 38, size=27, color=ORANGE, weight="bold")
    add_text(drawing, "TRAINING · generated once", 70, 65, size=16, color=ORANGE, weight="bold")
    add_text(drawing, "b", 582, 38, size=30, weight="bold")
    add_text(drawing, "REUSABLE RESPONSE LAYER", 624, 38, size=27, color=BLUE, weight="bold")
    add_text(drawing, "STAGE 1 → STAGE 2", 624, 65, size=16, color=BLUE, weight="bold")
    add_text(drawing, "c", 1273, 38, size=30, weight="bold")
    add_text(drawing, "STRUCTURE-ONLY DEPLOYMENT", 1315, 38, size=27, color=TEAL, weight="bold")
    add_text(drawing, "SMILES input", 1315, 65, size=16, color=TEAL, weight="bold")

    # Panel a: four scientifically grounded response-source vignettes.
    source_card(drawing, 35, 92, 235, "Calculated solvation", "COSMOtherm water · SMD(water)")
    embed(drawing, assets["cavity"], 43, 136, 220, 155)
    add_text(drawing, "RDKit van der Waals cavity", 152, 311, size=15, color=MID, anchor="middle")

    source_card(drawing, 300, 92, 225, "Polarity / H-bonding", "Abraham E · S · A · B · L")
    embed(drawing, assets["abraham"], 304, 133, 216, 168)
    add_text(drawing, "empirical solute coordinates", 412, 311, size=15, color=MID, anchor="middle")

    source_card(drawing, 35, 365, 235, "Water-model response", "OpenFF explicit · GBn2 implicit")
    embed(drawing, assets["lambda"], 43, 410, 220, 147)
    drawing.add(
        drawing.rect(
            insert=(43, 568),
            size=(220, 38),
            fill="white",
            stroke=MID,
            stroke_width=1.2,
            stroke_dasharray="5,4",
        )
    )
    add_text(drawing, "PIMD2 λ trace: explored, not retained", 153, 592, size=14, color=MID, anchor="middle")

    source_card(drawing, 300, 365, 225, "Conformational response", "ConfSolv H2O · six summaries")
    embed(drawing, assets["conformers"], 300, 420, 225, 150)
    add_text(drawing, "actual ETKDG/MMFF conformers", 412, 592, size=14, color=MID, anchor="middle")

    add_text(drawing, "STAGE 1", 35, 681, size=18, color=BLUE, weight="bold")
    add_text(drawing, "learn structure → response", 130, 681, size=18, color=BLUE, weight="bold")
    embed(drawing, assets["molecule"], 42, 712, 115, 105)
    add_text(drawing, REFERENCE_SMILES, 99, 839, size=17, color=MID, anchor="middle", family=MONO)
    arrow(drawing, arrowheads["blue"], (170, 765), (220, 765), color=BLUE, width=4)
    drawing.add(drawing.rect(insert=(232, 716), size=(285, 102), fill=BLUE_LIGHT, stroke=BLUE, stroke_width=2))
    add_text(drawing, "6 frozen response surrogates", 374, 753, size=19, color=BLUE_DARK, weight="bold", anchor="middle")
    add_text(drawing, "2 D-MPNN · 3 ExtraTrees · 1 LightGBM", 374, 785, size=15, color=MID, anchor="middle")
    add_text(drawing, "trained on benchmark-disjoint source data", 374, 807, size=14, color=MID, anchor="middle")
    stage_one_path = drawing.path(
        d="M 518 765 L 542 765 L 542 535 L 585 535",
        fill="none",
        stroke=BLUE,
        stroke_width=4,
    )
    stage_one_path["marker-end"] = arrowheads["blue"].get_funciri()
    drawing.add(stage_one_path)
    add_multiline(
        drawing,
        ["Source calculations", "and measurements", "do not run again"],
        280,
        904,
        size=20,
        leading=1.18,
        color=ORANGE,
        weight="bold",
        anchor="middle",
    )
    drawing.add(drawing.line(start=(84, 858), end=(476, 858), stroke=ORANGE, stroke_width=2))

    # Panel b: exact 15-prior vocabulary and matched endpoint architecture.
    add_text(drawing, "15 interpretable predicted response coordinates", 900, 93, size=23, color=BLUE_DARK, weight="bold", anchor="middle")
    embed(drawing, assets["priors"], 585, 108, 635, 445)
    arrow(drawing, arrowheads["blue"], (900, 556), (900, 592), color=BLUE, width=4)

    add_text(drawing, "STAGE 2", 590, 613, size=18, color=TEAL, weight="bold")
    add_text(drawing, "fit the hydration endpoint", 685, 613, size=18, color=TEAL, weight="bold")

    drawing.add(drawing.rect(insert=(590, 638), size=(270, 130), fill=PALE, stroke=INK, stroke_width=1.6))
    add_text(drawing, "MOLECULAR STRUCTURE", 725, 667, size=15, color=MID, weight="bold", anchor="middle")
    add_text(drawing, "2,048-bit Morgan", 725, 704, size=20, color=INK, weight="bold", anchor="middle")
    add_text(drawing, "+ 217 RDKit descriptors", 725, 734, size=18, color=INK, anchor="middle")
    add_text(drawing, "= 2,265 deterministic features", 725, 759, size=15, color=MID, anchor="middle")

    add_text(drawing, "+", 884, 716, size=45, color=MID, weight="bold", anchor="middle")
    drawing.add(drawing.rect(insert=(915, 638), size=(290, 130), fill=BLUE_LIGHT, stroke=BLUE, stroke_width=1.8))
    add_text(drawing, "PREDICTED RESPONSE", 1060, 667, size=15, color=BLUE_DARK, weight="bold", anchor="middle")
    add_text(drawing, "15 named scalar priors", 1060, 710, size=22, color=BLUE_DARK, weight="bold", anchor="middle")
    add_text(drawing, "from the frozen surrogates", 1060, 741, size=17, color=MID, anchor="middle")

    arrow(drawing, arrowheads["teal"], (725, 788), (820, 846), color=TEAL, width=3.5)
    arrow(drawing, arrowheads["teal"], (1060, 788), (970, 846), color=TEAL, width=3.5)
    drawing.add(drawing.rect(insert=(760, 852), size=(275, 94), fill=TEAL_LIGHT, stroke=TEAL, stroke_width=2))
    add_text(drawing, "ExtraTrees endpoint ensemble", 898, 887, size=18, color=TEAL, weight="bold", anchor="middle")
    add_text(drawing, "3 seeds × 360 trees", 898, 918, size=16, color=MID, anchor="middle")

    # Experimental labels supervise Stage 2 but are not deployment inputs.
    drawing.add(drawing.line(start=(1090, 839), end=(1172, 839), stroke=TEAL, stroke_width=2))
    drawing.add(drawing.line(start=(1131, 822), end=(1131, 856), stroke=TEAL, stroke_width=2))
    drawing.add(drawing.circle(center=(1131, 839), r=7, fill=TEAL))
    add_multiline(
        drawing,
        ["experimental hydration", "free-energy labels"],
        1131,
        784,
        size=15,
        leading=1.15,
        color=TEAL,
        weight="bold",
        anchor="middle",
    )
    arrow(drawing, arrowheads["teal"], (1131, 860), (1019, 884), color=TEAL, width=3)
    add_text(drawing, "frozen after evaluation", 898, 975, size=16, color=MID, anchor="middle")

    # Panel c: deployment repeats learned mappings only.
    add_text(drawing, "NEW MOLECULE", 1510, 101, size=17, color=MID, weight="bold", anchor="middle")
    add_text(drawing, REFERENCE_SMILES, 1510, 137, size=23, color=INK, anchor="middle", family=MONO)
    embed(drawing, assets["molecule"], 1420, 158, 180, 150)
    arrow(drawing, arrowheads["teal"], (1510, 318), (1510, 365), color=TEAL, width=4)

    drawing.add(drawing.line(start=(1325, 401), end=(1695, 401), stroke=LIGHT, stroke_width=1.5))
    add_text(drawing, "deterministic structure", 1390, 429, size=17, color=INK, weight="bold", anchor="middle")
    add_text(drawing, "2,265 features", 1390, 456, size=17, color=MID, anchor="middle")
    add_text(drawing, "+", 1510, 446, size=30, color=MID, weight="bold", anchor="middle")
    add_text(drawing, "frozen response surrogates", 1630, 429, size=17, color=BLUE, weight="bold", anchor="middle")
    add_text(drawing, "→ 15 priors", 1630, 456, size=17, color=BLUE, anchor="middle")
    arrow(drawing, arrowheads["teal"], (1510, 478), (1510, 535), color=TEAL, width=4)

    drawing.add(drawing.rect(insert=(1378, 545), size=(264, 104), fill=TEAL_LIGHT, stroke=TEAL, stroke_width=2.5))
    add_text(drawing, "SolvAI", 1510, 590, size=38, color=TEAL, weight="bold", anchor="middle")
    add_text(drawing, "frozen endpoint ensemble", 1510, 624, size=16, color=MID, anchor="middle")
    arrow(drawing, arrowheads["teal"], (1510, 660), (1510, 713), color=TEAL, width=4.5)
    embed(drawing, assets["delta_g"], 1415, 712, 190, 74)
    add_text(drawing, "kcal mol−1", 1510, 808, size=18, color=MID, anchor="middle")

    drawing.add(drawing.rect(insert=(1323, 842), size=(374, 54), fill=TEAL_LIGHT))
    add_text(drawing, "NO MD  ·  NO PIMD  ·  NO PROBE", 1510, 877, size=20, color=TEAL, weight="bold", anchor="middle")

    drawing.add(drawing.line(start=(1300, 939), end=(1720, 939), stroke=MAGENTA, stroke_width=2, stroke_dasharray="7,5"))
    add_text(drawing, "ARROW/PIMD8  ≈ 0.205 kcal mol−1", 1510, 970, size=18, color=MAGENTA, weight="bold", anchor="middle")
    add_text(drawing, "accuracy reference only · not a retained teacher", 1510, 996, size=15, color=MID, anchor="middle")

    drawing.save(pretty=True)
    normalize_svg(output)


def validate_science() -> None:
    config = json.loads((ROOT / "data/manifests/final_training_config.json").read_text())
    model_card = json.loads((ROOT / "models/final/model_card.json").read_text())
    if len(config["response_features"]) != 15 or len(PRIORS) != 15:
        raise AssertionError("Figure must contain exactly 15 retained priors")
    if config["response_features"] != CONFIG_RESPONSE_FEATURES:
        raise AssertionError("Figure response-prior order disagrees with the frozen model")
    if config["structure_features"] != {
        "rdkit_2d_descriptors": 217,
        "morgan_radius": 2,
        "morgan_bits": 2048,
    }:
        raise AssertionError("Unexpected deterministic structure-feature schema")
    prohibited = ("pimd", "lambda", "dhdl", "nqe")
    if any(any(term in name.casefold() for term in prohibited) for name in config["response_features"]):
        raise AssertionError("A prohibited PIMD/lambda/NQE feature entered the retained vector")
    if config["selected_pimd8_training_labels"] != 0 or config["simulation_at_inference"]:
        raise AssertionError("Final model boundary disagrees with figure")
    if config["endpoint"] != {
        "algorithm": "ExtraTreesRegressor",
        "trees": 360,
        "max_features": 0.7,
        "min_samples_leaf": 2,
        "seeds": [11, 29, 47],
        "external_endpoint_weight": 1,
        "arrow_outer_training_weight": 3,
    }:
        raise AssertionError("Endpoint architecture disagrees with the frozen configuration")
    if model_card["pimd8_labels_in_selected_artifact"] != 0:
        raise AssertionError("PIMD8 must remain a comparator")
    table = (ROOT / "paper/supplementary/tables/response_priors.tex").read_text()
    expected_names = [item[3] for item in PRIORS]
    if not all(name.replace("_", r"\_") in table for name in expected_names):
        missing = [name for name in expected_names if name.replace("_", r"\_") not in table]
        raise AssertionError(f"Figure prior name not present in Supplementary Table 1: {missing}")


def write_manifest(paths: dict[str, Path]) -> None:
    manifest = {
        "schema_version": 1,
        "reference_molecule": {"name": REFERENCE_NAME, "smiles": REFERENCE_SMILES},
        "conformer_molecule": {"name": CONFORMER_NAME, "smiles": CONFORMER_SMILES},
        "lambda_source": "results/ablations/pimd2_multilambda_teacher.parquet",
        "lambda_status": "explored training-time supervision; not retained in SolvAI",
        "retained_prior_count": len(PRIORS),
        "assets": {
            name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for name, path in sorted(paths.items())
        },
    }
    (ASSETS / "asset_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    configure_matplotlib()
    validate_science()
    ASSETS.mkdir(parents=True, exist_ok=True)
    MAIN.mkdir(parents=True, exist_ok=True)
    PAPER_MAIN.mkdir(parents=True, exist_ok=True)

    molecule_canvas_path = ASSETS / "n_methylacetamide_rdkit_canvas.svg"
    molecule_path = ASSETS / "n_methylacetamide_rdkit.svg"
    molecule, coordinates = molecule_svg(REFERENCE_SMILES, molecule_canvas_path, 440, 250)
    crop_rdkit_svg(molecule_canvas_path, molecule_path, coordinates)
    indigo_path = ASSETS / "n_methylacetamide_indigo.svg"
    indigo_molecule_svg(REFERENCE_SMILES, indigo_path)
    cavity_path = ASSETS / "continuum_cavity.svg"
    cavity_svg(molecule, coordinates, molecule_canvas_path, cavity_path)
    water_shell_path = ASSETS / "nma_water_shell_vector.svg"
    water_shell_vector_svg(ASSETS / "nma_water_shell.pdb", water_shell_path)

    abraham_path = ASSETS / "abraham_axes.svg"
    abraham_axes_svg(abraham_path)

    conformer_path = ASSETS / "dimethoxyethane_conformers.svg"
    conformer_svg(
        conformer_path,
        ASSETS / "dimethoxyethane_selected_conformers.sdf",
        ASSETS / "dimethoxyethane_conformer_metadata.json",
    )

    lambda_path = ASSETS / "nmethylacetamide_lambda_response.svg"
    lambda_response_svg(lambda_path, ASSETS / "nmethylacetamide_lambda_response.csv")

    priors_path = ASSETS / "response_priors_15.svg"
    response_priors_svg(priors_path)

    delta_g_path = ASSETS / "delta_g_hyd.svg"
    math_label_svg(delta_g_path)

    fingerprint_path = ASSETS / "n_methylacetamide_morgan2048.svg"
    fingerprint_svg(REFERENCE_SMILES, fingerprint_path)

    assets = {
        "molecule": molecule_path,
        "molecule_indigo": indigo_path,
        "cavity": cavity_path,
        "water_shell": water_shell_path,
        "abraham": abraham_path,
        "conformers": conformer_path,
        "lambda": lambda_path,
        "priors": priors_path,
        "delta_g": delta_g_path,
        "fingerprint": fingerprint_path,
    }

    from fig1_compositions import build_all_variants

    build_all_variants(assets, PRIORS, MAIN, PAPER_MAIN)
    svg_output = MAIN / "fig1_concept.svg"

    paths = {
        **assets,
        "molecule_canvas": molecule_canvas_path,
        "conformer_sdf": ASSETS / "dimethoxyethane_selected_conformers.sdf",
        "conformer_metadata": ASSETS / "dimethoxyethane_conformer_metadata.json",
        "lambda_values": ASSETS / "nmethylacetamide_lambda_response.csv",
        "water_shell_packmol_input": ASSETS / "water_shell.packmol.inp",
        "water_shell_pdb": ASSETS / "nma_water_shell.pdb",
        "water_shell_pymol_script": ASSETS / "water_shell.pml",
        "water_shell_pymol_preview": ASSETS / "nma_water_shell_pymol.png",
        "figure_svg": svg_output,
        "figure_pdf": MAIN / "fig1_concept.pdf",
        "figure_png": MAIN / "fig1_concept.png",
    }
    for stem, directory in (
        ("fig1_variant_A_minimal", MAIN.parent / "alternatives"),
        ("fig1_variant_B_molecular", MAIN.parent / "alternatives"),
        ("fig1_variant_C_balanced", MAIN),
    ):
        for suffix in ("svg", "pdf", "png"):
            paths[f"{stem}_{suffix}"] = directory / f"{stem}.{suffix}"
        paths[f"{stem}_print_180mm"] = directory / f"{stem}_print_180mm.png"
    write_manifest(paths)
    print(f"Built {svg_output.relative_to(ROOT)} from {len(PRIORS)} verified response priors")


if __name__ == "__main__":
    main()
