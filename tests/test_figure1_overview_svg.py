from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
SVG = ROOT / "paper/figures/main/fig1_concept.svg"
GENERATOR = ROOT / "scripts/make_figure1_overview.py"
FIGURE_DRIVER = ROOT / "scripts/make_figures.py"
PENPOT_INSTALLER = ROOT / "figures/penpot/install_exports.py"
NS = "{http://www.w3.org/2000/svg}"


def parsed_svg():
    root = ET.parse(SVG).getroot()
    text = " ".join("".join(root.itertext()).split())
    return root, text


def test_svg_is_editable_hybrid_with_only_curated_embedded_assets():
    root, _ = parsed_svg()
    raw = SVG.read_text(encoding="utf-8")
    assert root.attrib["width"] == "518.4pt"
    assert root.attrib["height"] == "424.8pt"
    images = root.findall(f".//{NS}image")
    assert len(images) == 30
    assert all(
        image.attrib["{http://www.w3.org/1999/xlink}href"].startswith("data:image/png;base64,")
        for image in images
    )
    assert len(root.findall(f".//{NS}text")) >= 80
    assert raw.count("data:font/otf;base64,") == 3
    assert "font-family: 'Latin Modern Roman'" in raw
    assert "font-family: 'Latin Modern Mono'" in raw
    assert "DejaVu" not in raw


def test_svg_preserves_core_scientific_constraints():
    _, text = parsed_svg()
    required = (
        "ONE QUERY MOLECULE → SIX FROZEN SURROGATE FAMILIES",
        "same molecule in every branch",
        "CombiSolv-QM",
        "Abraham",
        "OpenFF corrected",
        "GBn2 corrected",
        "SMD(water)",
        "ConfSolv",
        "CHEMELEON · D-MPNN",
        "ExtraTrees",
        "LightGBM",
        "15 molecule-aligned response coordinates",
        "1 | 5 | 1 | 1 | 1 | 6 = 15",
        "predicted per molecule by frozen structure→response surrogates",
        "source calculations are not run at inference",
        "CONCATENATE",
        "2,280-D molecular representation",
        "3 × ExtraTrees ensembles",
        "360 trees each · seeds 11 / 29 / 47",
        "mean prediction",
        "six frozen surrogates",
        "2 D-MPNN · 3 ExtraTrees · 1 LightGBM",
        "3 × frozen ExtraTrees ensembles",
        "360 trees each · predictions averaged",
        "No PIMD-derived feature in the final stack",
        "NO MD · NO PIMD · NO PROBE CALCULATION",
        "SMILES input",
        "CC(=O)N",
        "acetamide",
        "shuffling abolishes the gain",
        "PIMD8-level accuracy on this reference chemistry",
        "accuracy reference; not a teacher · never used as a SolvAI input",
        "0.303",
        "0.306",
        "0.202",
        "0.205",
    )
    for phrase in required:
        assert phrase in text


def test_panels_d_and_e_use_the_preferred_flat_forest_asset():
    source = GENERATOR.read_text(encoding="utf-8")
    assert '"endpoint_ensemble_v1.png"' not in source
    assert '"frozen_model_v1.png"' not in source
    assert '"endpoint_extratrees_triplet_azure_v1.png"' not in source
    assert source.count('"endpoint_extratrees_forest_flat_v3.png"') == 2
    assert source.count('"frozen_surrogate_wedge_flat_v2.png"') == 1
    assert '"endpoint_extratrees_wedge_flat_v2.png"' not in source
    assert '"endpoint_extratrees_wedge_v1.png"' not in source
    assert '"frozen_surrogate_wedge_v1.png"' not in source
    assert source.count('"lock_icon_flat_azure_v1.png"') == 1
    assert source.count('"model_icon_dmpnn_azure_v1.png"') == 2
    assert source.count('"model_icon_extratrees_azure_v1.png"') == 3
    assert source.count('"model_icon_lightgbm_azure_v1.png"') == 1
    assert "The input connectors meet at one explicit hub" in source
    assert "A compact segmented vector replaces the former text-heavy box" in source
    assert "ax.scatter([merge_x], [pipeline_y]" in source


def test_penpot_variants_cannot_overwrite_the_reviewed_hybrid_figure():
    installer = PENPOT_INSTALLER.read_text(encoding="utf-8")
    driver = FIGURE_DRIVER.read_text(encoding="utf-8")
    assert "fig1_concept" not in installer
    assert "make_figure1_overview.py" in driver
