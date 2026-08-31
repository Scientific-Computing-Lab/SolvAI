"""Art-directed SVG compositions for SolvAI Figure 1.

All chemistry and quantitative panels are embedded from deterministic assets made
by ``fig1_build.py``. This module contains only typography, layout and abstract
architecture marks; it does not manufacture molecular or numerical content.
"""

from __future__ import annotations

import base64
import re
import shutil
from pathlib import Path

import cairosvg
import svgwrite

WIDTH = 1800
HEIGHT = 1000
PNG_WIDTH = 3600

INK = "#202A33"
MID = "#67737D"
LIGHT = "#D9E0E5"
PALE = "#F6F8F9"
ORANGE = "#D78324"
ORANGE_LIGHT = "#F7E8D5"
BLUE = "#2278B5"
BLUE_DARK = "#185986"
BLUE_LIGHT = "#E4F0F7"
TEAL = "#078D78"
TEAL_LIGHT = "#DDF1EC"
MAGENTA = "#B33A74"
FONT = "Arial, Liberation Sans, DejaVu Sans, sans-serif"
MONO = "DejaVu Sans Mono, Liberation Mono, monospace"


def uri(path: Path) -> str:
    return "data:image/svg+xml;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def text(
    drawing: svgwrite.Drawing,
    value: str,
    x: float,
    y: float,
    *,
    size: float = 17,
    color: str = INK,
    bold: bool = False,
    anchor: str = "start",
    family: str = FONT,
) -> None:
    drawing.add(
        drawing.text(
            value,
            insert=(x, y),
            text_anchor=anchor,
            font_family=family,
            font_size=size,
            font_weight="bold" if bold else "normal",
            fill=color,
        )
    )


def multiline(
    drawing: svgwrite.Drawing,
    lines: list[str],
    x: float,
    y: float,
    *,
    size: float = 16,
    color: str = INK,
    bold: bool = False,
    anchor: str = "start",
    leading: float = 1.2,
) -> None:
    node = drawing.text(
        "",
        insert=(x, y),
        text_anchor=anchor,
        font_family=FONT,
        font_size=size,
        font_weight="bold" if bold else "normal",
        fill=color,
    )
    for index, line in enumerate(lines):
        node.add(drawing.tspan(line, x=[x], dy=[0 if index == 0 else size * leading]))
    drawing.add(node)


def embed(
    drawing: svgwrite.Drawing,
    path: Path,
    x: float,
    y: float,
    width: float,
    height: float,
) -> None:
    drawing.add(drawing.image(href=uri(path), insert=(x, y), size=(width, height)))


def marker(drawing: svgwrite.Drawing, color: str) -> svgwrite.container.Marker:
    arrowhead = drawing.marker(insert=(7, 4), size=(8, 8), orient="auto", markerUnits="strokeWidth")
    arrowhead.add(drawing.path(d="M 0 0 L 8 4 L 0 8 z", fill=color))
    drawing.defs.add(arrowhead)
    return arrowhead


def arrow(
    drawing: svgwrite.Drawing,
    arrowhead: svgwrite.container.Marker,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    width: float = 3.2,
) -> None:
    line = drawing.line(start=start, end=end, stroke=color, stroke_width=width)
    line["marker-end"] = arrowhead.get_funciri()
    drawing.add(line)


def base_canvas(dividers: tuple[int, int]) -> tuple[svgwrite.Drawing, dict[str, svgwrite.container.Marker]]:
    drawing = svgwrite.Drawing(size=("180mm", "100mm"), viewBox=f"0 0 {WIDTH} {HEIGHT}", profile="full")
    drawing.add(drawing.rect(insert=(0, 0), size=(WIDTH, HEIGHT), fill="white"))
    for x in dividers:
        drawing.add(drawing.line(start=(x, 54), end=(x, 958), stroke=LIGHT, stroke_width=1.5))
    return drawing, {
        "orange": marker(drawing, ORANGE),
        "blue": marker(drawing, BLUE),
        "teal": marker(drawing, TEAL),
    }


def panel_heading(
    drawing: svgwrite.Drawing,
    label: str,
    title: str,
    subtitle: str,
    x: float,
    color: str,
) -> None:
    text(drawing, label, x, 38, size=30, bold=True)
    text(drawing, title, x + 42, 38, size=27, color=color, bold=True)
    text(drawing, subtitle, x + 42, 66, size=16, color=color, bold=True)


def source_title(drawing: svgwrite.Drawing, x: float, y: float, title: str, subtitle: str, width: float) -> None:
    text(drawing, title, x, y, size=20, bold=True)
    text(drawing, subtitle, x, y + 25, size=15, color=MID)
    drawing.add(drawing.line(start=(x, y + 36), end=(x + width, y + 36), stroke=ORANGE, stroke_width=2))


def surrogate_bank(drawing: svgwrite.Drawing, x: float, y: float, *, compact: bool = False) -> None:
    rows = [(0, 0), (42, 0), (84, 0), (0, 42), (42, 42), (84, 42)]
    for index, (dx, dy) in enumerate(rows):
        drawing.add(
            drawing.circle(
                center=(x + dx, y + dy),
                r=10 if compact else 12,
                fill=BLUE if index < 2 else BLUE_LIGHT,
                stroke=BLUE,
                stroke_width=1.5,
            )
        )


def response_dots(drawing: svgwrite.Drawing, x: float, y: float, *, scale: float = 1.0) -> None:
    counts = (2, 5, 2, 6)
    offsets = (0, 52, 151, 203)
    for group, (count, offset) in enumerate(zip(counts, offsets, strict=True)):
        color = ["#4F9BC8", "#76A9C8", "#2478B5", "#2D678F"][group]
        for index in range(count):
            drawing.add(
                drawing.circle(
                    center=(x + (offset + index * 15) * scale, y),
                    r=4.8 * scale,
                    fill=color,
                )
            )


def forest(drawing: svgwrite.Drawing, x: float, y: float, *, scale: float = 1.0) -> None:
    for offset in (-34, 0, 34):
        x0 = x + offset * scale
        drawing.add(drawing.line(start=(x0, y + 32 * scale), end=(x0, y), stroke=TEAL, stroke_width=2.4 * scale))
        drawing.add(drawing.line(start=(x0, y + 9 * scale), end=(x0 - 11 * scale, y - 6 * scale), stroke=TEAL, stroke_width=2.0 * scale))
        drawing.add(drawing.line(start=(x0, y + 9 * scale), end=(x0 + 11 * scale, y - 6 * scale), stroke=TEAL, stroke_width=2.0 * scale))
        for cx, cy in ((x0, y - 2 * scale), (x0 - 11 * scale, y - 8 * scale), (x0 + 11 * scale, y - 8 * scale)):
            drawing.add(drawing.circle(center=(cx, cy), r=3.4 * scale, fill=TEAL))


def experimental_symbol(drawing: svgwrite.Drawing, x: float, y: float) -> None:
    drawing.add(drawing.line(start=(x - 34, y), end=(x + 34, y), stroke=TEAL, stroke_width=2))
    drawing.add(drawing.line(start=(x, y - 20), end=(x, y + 20), stroke=TEAL, stroke_width=2))
    drawing.add(drawing.circle(center=(x, y), r=7, fill=TEAL))


def prior_atlas(
    drawing: svgwrite.Drawing,
    priors: list[tuple[int, str, str, str]],
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    compact: bool = False,
) -> None:
    wrapped = {
        "Abraham E · excess molar refraction": ["Abraham E · excess molar", "refraction"],
        "Abraham S · dipolarity / polarizability": ["Abraham S · dipolarity /", "polarizability"],
        "Abraham L · hexadecane–air partition": ["Abraham L · hexadecane–air", "partition"],
    }
    groups = [
        ("CALCULATED SOLVATION", "Calculated solvation", (x, y), width * 0.43),
        ("POLARITY / H-BONDING", "Polarity / H-bonding", (x + width * 0.52, y), width * 0.48),
        ("EXPLICIT / IMPLICIT WATER", "Explicit / implicit water", (x, y + height * 0.51), width * 0.43),
        ("CONFORMATIONAL RESPONSE", "Conformational response", (x + width * 0.52, y + height * 0.51), width * 0.48),
    ]
    centre = (x + width * 0.49, y + height * 0.48)
    drawing.add(drawing.circle(center=centre, r=39 if not compact else 31, fill="white", stroke=BLUE, stroke_width=2.2))
    text(drawing, "15", centre[0], centre[1] - 2, size=25 if not compact else 21, color=BLUE_DARK, bold=True, anchor="middle")
    text(drawing, "scalars", centre[0], centre[1] + 18, size=13, color=MID, anchor="middle")
    for title, key, (gx, gy), group_width in groups:
        items = [item for item in priors if item[1] == key]
        anchor_x = gx + (group_width if gx < centre[0] else 0)
        anchor_y = gy + 35
        drawing.add(drawing.line(start=centre, end=(anchor_x, anchor_y), stroke="#B9D4E6", stroke_width=1.4))
        drawing.add(drawing.line(start=(gx, gy + 23), end=(gx + group_width, gy + 23), stroke=BLUE, stroke_width=1.8))
        text(drawing, title, gx, gy + 15, size=14 if compact else 15, color=BLUE_DARK, bold=True)
        row_step = 36 if len(items) <= 2 else (34 if len(items) == 5 else 29)
        for index, (number, _, label, _) in enumerate(items):
            row_y = gy + 51 + index * row_step
            drawing.add(drawing.circle(center=(gx + 10, row_y - 5), r=8, fill=BLUE))
            text(drawing, str(number), gx + 10, row_y - 1, size=10.5, color="white", bold=True, anchor="middle")
            lines = wrapped.get(label, [label])
            multiline(
                drawing,
                lines,
                gx + 25,
                row_y,
                size=14.5 if compact else 15.5,
                leading=1.0,
                color=INK,
            )


def stage_two(
    drawing: svgwrite.Drawing,
    arrows: dict[str, svgwrite.container.Marker],
    assets: dict[str, Path],
    *,
    x: float,
    y: float,
    width: float,
    compact: bool = False,
) -> None:
    text(drawing, "STAGE 2", x, y, size=17, color=TEAL, bold=True)
    text(drawing, "experimentally supervised endpoint", x + 90, y, size=17, color=TEAL, bold=True)
    centres = (x + 105, x + width * 0.52, x + width - 78)
    embed(drawing, assets["fingerprint"], centres[0] - 64, y + 30, 128, 64)
    multiline(drawing, ["2,265 structure features", "2,048 Morgan + 217 RDKit"], centres[0], y + 115, size=14, color=INK, bold=True, anchor="middle")
    response_dots(drawing, centres[1] - 82, y + 64, scale=0.70)
    multiline(drawing, ["15 predicted responses", "named scalar coordinates"], centres[1], y + 115, size=14, color=BLUE_DARK, bold=True, anchor="middle")
    experimental_symbol(drawing, centres[2], y + 64)
    multiline(drawing, ["experimental ΔGhyd", "training labels"], centres[2], y + 115, size=14, color=TEAL, bold=True, anchor="middle")
    join_y = y + 205
    endpoint_x = x + width / 2
    for centre, target_x in zip(
        centres,
        (endpoint_x - 38, endpoint_x, endpoint_x + 38),
        strict=True,
    ):
        arrow(drawing, arrows["teal"], (centre, y + 150), (target_x, join_y), color=TEAL, width=2.5)
    forest(drawing, endpoint_x, y + 230, scale=0.9 if compact else 1.0)
    text(drawing, "ExtraTrees endpoint ensemble", endpoint_x, y + 285, size=17, color=TEAL, bold=True, anchor="middle")
    text(drawing, "3 × 360 trees · frozen after evaluation", endpoint_x, y + 311, size=14, color=MID, anchor="middle")


def deployment(
    drawing: svgwrite.Drawing,
    arrows: dict[str, svgwrite.container.Marker],
    assets: dict[str, Path],
    *,
    x: float,
    width: float,
    compact: bool = False,
) -> None:
    centre = x + width / 2
    text(drawing, "NEW MOLECULE", centre, 103, size=16, color=MID, bold=True, anchor="middle")
    text(drawing, "CNC(C)=O", centre, 137, size=22, color=INK, anchor="middle", family=MONO)
    embed(drawing, assets["molecule"], centre - 85, 157, 170, 140)
    arrow(drawing, arrows["teal"], (centre, 306), (centre, 350), color=TEAL)

    branch_y = 405
    drawing.add(drawing.line(start=(centre - 160, branch_y), end=(centre + 160, branch_y), stroke=LIGHT, stroke_width=1.5))
    embed(drawing, assets["fingerprint"], centre - 170, branch_y + 32, 110, 55)
    text(drawing, "structure", centre - 115, branch_y + 112, size=15, color=INK, bold=True, anchor="middle")
    text(drawing, "2,265 features", centre - 115, branch_y + 135, size=14, color=MID, anchor="middle")
    text(drawing, "+", centre, branch_y + 78, size=31, color=MID, bold=True, anchor="middle")
    response_dots(drawing, centre + 44, branch_y + 61, scale=0.66)
    text(drawing, "frozen surrogates", centre + 120, branch_y + 112, size=15, color=BLUE, bold=True, anchor="middle")
    text(drawing, "15 responses", centre + 120, branch_y + 135, size=14, color=BLUE, anchor="middle")
    arrow(drawing, arrows["teal"], (centre, branch_y + 155), (centre, 610), color=TEAL)

    forest(drawing, centre, 649, scale=1.05)
    text(drawing, "SolvAI", centre, 720, size=35, color=TEAL, bold=True, anchor="middle")
    text(drawing, "frozen endpoint ensemble", centre, 746, size=14, color=MID, anchor="middle")
    arrow(drawing, arrows["teal"], (centre, 765), (centre, 806), color=TEAL)
    embed(drawing, assets["delta_g"], centre - 88, 806, 176, 66)
    text(drawing, "kcal mol−1", centre, 890, size=16, color=MID, anchor="middle")
    text(drawing, "STRUCTURE ONLY", centre, 907, size=17, color=TEAL, bold=True, anchor="middle")
    text(drawing, "no MD · no PIMD · no probe", centre, 931, size=14, color=MID, anchor="middle")


def comparator_footer(drawing: svgwrite.Drawing, x: float, width: float) -> None:
    drawing.add(drawing.line(start=(x + 30, 943), end=(x + width - 30, 943), stroke=MAGENTA, stroke_width=1.8, stroke_dasharray="7,5"))
    text(drawing, "ARROW/PIMD8 ≈ 0.205 kcal mol−1", x + width / 2, 966, size=14, color=MAGENTA, bold=True, anchor="middle")
    text(drawing, "accuracy reference · not a retained teacher", x + width / 2, 987, size=12.5, color=MID, anchor="middle")


def minimal_variant(assets: dict[str, Path], priors: list[tuple[int, str, str, str]], output: Path) -> None:
    drawing, arrows = base_canvas((500, 1280))
    panel_heading(drawing, "a", "RESPONSE SUPERVISION", "paid once", 24, ORANGE)
    embed(drawing, assets["cavity"], 120, 105, 255, 165)
    source_title(drawing, 52, 300, "Four complementary families", "calculated · empirical · water-model · conformational", 390)
    embed(drawing, assets["abraham"], 55, 350, 130, 90)
    embed(drawing, assets["lambda"], 195, 345, 135, 90)
    embed(drawing, assets["conformers"], 340, 350, 125, 86)
    drawing.add(drawing.line(start=(55, 477), end=(445, 477), stroke=ORANGE, stroke_width=2))
    text(drawing, "external response datasets", 250, 510, size=17, color=ORANGE, bold=True, anchor="middle")
    arrow(drawing, arrows["blue"], (250, 535), (250, 612), color=BLUE)
    embed(drawing, assets["molecule"], 72, 620, 112, 95)
    surrogate_bank(drawing, 270, 646)
    text(drawing, "structure → response", 312, 730, size=18, color=BLUE_DARK, bold=True, anchor="middle")
    text(drawing, "six frozen surrogate models", 312, 756, size=15, color=MID, anchor="middle")
    arrow(drawing, arrows["blue"], (385, 680), (497, 680), color=BLUE)
    multiline(drawing, ["The costly information is generated once", "and reused through learned molecular responses."], 250, 850, size=18, color=ORANGE, bold=True, anchor="middle")

    panel_heading(drawing, "b", "INTERPRETABLE RESPONSE LAYER", "15 named scalar coordinates", 522, BLUE)
    prior_atlas(drawing, priors, x=535, y=100, width=720, height=455, compact=True)
    drawing.add(drawing.line(start=(540, 575), end=(1245, 575), stroke=LIGHT, stroke_width=1.5))
    stage_two(drawing, arrows, assets, x=550, y=610, width=675, compact=True)

    panel_heading(drawing, "c", "STRUCTURE-ONLY DEPLOYMENT", "frozen mappings", 1300, TEAL)
    deployment(drawing, arrows, assets, x=1285, width=500, compact=True)
    comparator_footer(drawing, 1285, 500)
    drawing.saveas(output, pretty=True)


def molecular_variant(assets: dict[str, Path], priors: list[tuple[int, str, str, str]], output: Path) -> None:
    drawing, arrows = base_canvas((620, 1260))
    panel_heading(drawing, "a", "MOLECULE IN ITS RESPONSE SPACE", "real coordinates and frozen response data", 24, ORANGE)
    embed(drawing, assets["water_shell"], 55, 100, 510, 330)
    text(drawing, "Packmol solute–water configuration · vector projection", 310, 446, size=15, color=MID, anchor="middle")
    source_title(drawing, 45, 495, "Conformation", "ETKDGv3 + MMFF", 245)
    embed(drawing, assets["conformers"], 50, 545, 240, 145)
    source_title(drawing, 330, 495, "Alchemical response", "measured PIMD2 trace", 245)
    embed(drawing, assets["lambda"], 338, 538, 228, 152)
    text(drawing, "explored · not retained", 452, 714, size=14, color=MID, anchor="middle")
    drawing.add(drawing.line(start=(55, 758), end=(565, 758), stroke=ORANGE, stroke_width=2))
    embed(drawing, assets["molecule"], 72, 780, 110, 90)
    arrow(drawing, arrows["blue"], (190, 825), (255, 825), color=BLUE)
    surrogate_bank(drawing, 315, 805)
    multiline(drawing, ["structure→response", "surrogates"], 410, 817, size=18, color=BLUE_DARK, bold=True)
    arrow(drawing, arrows["blue"], (510, 825), (617, 825), color=BLUE)
    text(drawing, "calculation is paid once", 310, 924, size=18, color=ORANGE, bold=True, anchor="middle")

    panel_heading(drawing, "b", "RESPONSE COORDINATE ATLAS", "interpretable, structure-predicted scalars", 642, BLUE)
    prior_atlas(drawing, priors, x=655, y=105, width=575, height=500, compact=True)
    stage_two(drawing, arrows, assets, x=650, y=635, width=570, compact=True)

    panel_heading(drawing, "c", "DEPLOYMENT", "SMILES only", 1282, TEAL)
    deployment(drawing, arrows, assets, x=1270, width=515)
    comparator_footer(drawing, 1270, 515)
    drawing.saveas(output, pretty=True)


def balanced_variant(assets: dict[str, Path], priors: list[tuple[int, str, str, str]], output: Path) -> None:
    drawing, arrows = base_canvas((540, 1260))
    panel_heading(drawing, "a", "SOLVATION RESPONSE", "training information · generated once", 24, ORANGE)
    source_title(drawing, 38, 105, "Calculated solvation", "COSMOtherm water · SMD(water)", 220)
    embed(drawing, assets["cavity"], 45, 150, 205, 145)
    source_title(drawing, 285, 105, "Polarity / H-bonding", "Abraham E · S · A · B · L", 220)
    embed(drawing, assets["abraham"], 295, 146, 200, 150)
    source_title(drawing, 38, 345, "Water-model response", "OpenFF explicit · GBn2 implicit", 220)
    embed(drawing, assets["lambda"], 45, 390, 205, 137)
    text(drawing, "PIMD2 λ diagnostic · not retained", 147, 550, size=14, color=MID, anchor="middle")
    source_title(drawing, 285, 345, "Conformational response", "ConfSolv H2O · six summaries", 220)
    embed(drawing, assets["conformers"], 288, 399, 215, 135)
    text(drawing, "ETKDG/MMFF conformers", 395, 550, size=14, color=MID, anchor="middle")
    drawing.add(drawing.line(start=(42, 606), end=(500, 606), stroke=ORANGE, stroke_width=2))
    text(drawing, "STAGE 1", 42, 644, size=17, color=BLUE, bold=True)
    text(drawing, "learn structure → response", 137, 644, size=17, color=BLUE, bold=True)
    embed(drawing, assets["molecule"], 47, 685, 108, 92)
    arrow(drawing, arrows["blue"], (165, 731), (222, 731), color=BLUE)
    surrogate_bank(drawing, 274, 710)
    multiline(drawing, ["six frozen", "response surrogates"], 414, 721, size=16, color=BLUE_DARK, bold=True, anchor="middle")
    arrow(drawing, arrows["blue"], (485, 731), (537, 731), color=BLUE)
    multiline(drawing, ["CALCULATIONS + MEASUREMENTS", "do not run again"], 270, 876, size=17, color=ORANGE, bold=True, anchor="middle")

    panel_heading(drawing, "b", "REUSABLE RESPONSE LAYER", "15 interpretable scalar coordinates", 562, BLUE)
    prior_atlas(drawing, priors, x=575, y=105, width=655, height=455, compact=False)
    drawing.add(drawing.line(start=(570, 579), end=(1230, 579), stroke=LIGHT, stroke_width=1.5))
    stage_two(drawing, arrows, assets, x=580, y=615, width=635, compact=False)

    panel_heading(drawing, "c", "STRUCTURE-ONLY DEPLOYMENT", "SMILES input", 1282, TEAL)
    deployment(drawing, arrows, assets, x=1270, width=515)
    comparator_footer(drawing, 1270, 515)
    drawing.saveas(output, pretty=True)


def normalize_svg(path: Path) -> None:
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


def normalize_pdf(path: Path) -> None:
    payload = path.read_bytes()
    payload, count = re.subn(
        rb"/CreationDate \(D:\d{14}Z\)",
        b"/CreationDate (D:20000101000000Z)",
        payload,
    )
    if count != 1:
        raise AssertionError(f"Expected one PDF CreationDate, found {count}")
    path.write_bytes(payload)


def export(svg_path: Path) -> None:
    normalize_svg(svg_path)
    pdf_path = svg_path.with_suffix(".pdf")
    png_path = svg_path.with_suffix(".png")
    cairosvg.svg2pdf(bytestring=svg_path.read_bytes(), write_to=str(pdf_path))
    normalize_pdf(pdf_path)
    cairosvg.svg2png(
        bytestring=svg_path.read_bytes(),
        write_to=str(png_path),
        output_width=PNG_WIDTH,
        output_height=round(PNG_WIDTH * HEIGHT / WIDTH),
    )
    cairosvg.svg2png(
        bytestring=svg_path.read_bytes(),
        write_to=str(svg_path.with_name(f"{svg_path.stem}_print_180mm.png")),
        output_width=round(180 / 25.4 * 300),
        output_height=round(100 / 25.4 * 300),
    )


def build_all_variants(
    assets: dict[str, Path],
    priors: list[tuple[int, str, str, str]],
    main_dir: Path,
    paper_main_dir: Path,
) -> None:
    alternatives = main_dir.parent / "alternatives"
    alternatives.mkdir(parents=True, exist_ok=True)
    builders = [
        (alternatives, "fig1_variant_A_minimal", minimal_variant),
        (alternatives, "fig1_variant_B_molecular", molecular_variant),
        (main_dir, "fig1_variant_C_balanced", balanced_variant),
    ]
    for directory, stem, builder in builders:
        svg_path = directory / f"{stem}.svg"
        builder(assets, priors, svg_path)
        export(svg_path)

    for suffix in ("svg", "pdf", "png"):
        winner = main_dir / f"fig1_variant_C_balanced.{suffix}"
        shutil.copy2(winner, main_dir / f"fig1_concept.{suffix}")
        shutil.copy2(winner, paper_main_dir / f"fig1_concept.{suffix}")
