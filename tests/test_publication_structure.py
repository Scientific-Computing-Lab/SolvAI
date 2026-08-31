from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "paper/main.tex"
SUPPLEMENTARY = ROOT / "paper/supplementary/supplementary.tex"
MAKEFILE = ROOT / "Makefile"

SUPPLEMENTARY_FIGURES = (
    "Supp_Fig1_residuals",
    "Supp_Fig2_provenance",
    "Supp_Fig3_alternatives",
    "Supp_Fig4_lambda_response",
    "Supp_Fig5_extrapolation",
)


def test_supporting_figures_are_embedded_and_cited() -> None:
    main = MAIN.read_text(encoding="utf-8")
    supplementary = SUPPLEMENTARY.read_text(encoding="utf-8")

    assert "Extended Data" not in main
    for number, stem in enumerate(SUPPLEMENTARY_FIGURES, start=1):
        assert f"Supplementary Fig.~{number}" in main
        assert f"Supplementary Figure {number} |" in supplementary
        assert f"{{{stem}.pdf}}" in supplementary
        for root in (ROOT / "figures/supplementary", ROOT / "paper/supplementary/figures"):
            for suffix in ("pdf", "svg", "png"):
                assert (root / f"{stem}.{suffix}").is_file()


def test_review_bundle_is_main_followed_by_complete_supplementary_information() -> None:
    makefile = MAKEFILE.read_text(encoding="utf-8")

    assert "paper/extended_data" not in makefile
    assert (
        "--pages paper/main.pdf paper/supplementary/supplementary.pdf "
        "-- paper/review_combined.pdf"
    ) in makefile


def test_discontinued_selective_pimd_is_not_in_submission_sources() -> None:
    publication_text = MAIN.read_text(encoding="utf-8") + SUPPLEMENTARY.read_text(
        encoding="utf-8"
    )

    assert "selective PIMD" not in publication_text
    assert not (ROOT / "paper/extended_data").exists()
    assert (ROOT / "figures/diagnostics/selective_pimd_reference.pdf").is_file()
