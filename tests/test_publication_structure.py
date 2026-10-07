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
            for suffix in ("pdf", "svg"):
                assert (root / f"{stem}.{suffix}").is_file()


def test_discontinued_selective_pimd_is_not_in_submission_sources() -> None:
    publication_text = MAIN.read_text(encoding="utf-8") + SUPPLEMENTARY.read_text(
        encoding="utf-8"
    )

    assert "selective PIMD" not in publication_text
    assert not (ROOT / "paper/extended_data").exists()
    assert (ROOT / "figures/diagnostics/selective_pimd_reference.pdf").is_file()


def test_journal_has_one_compilation_root() -> None:
    main = MAIN.read_text(encoding="utf-8")
    supplementary = SUPPLEMENTARY.read_text(encoding="utf-8")
    makefile = MAKEFILE.read_text(encoding="utf-8")
    assert r"\input{supplementary/supplementary.tex}" in main
    assert main.index(r"\bibliography{references}") < main.index(
        r"\input{supplementary/supplementary.tex}"
    ) < main.index(r"\end{document}")
    assert r"\documentclass" not in supplementary
    assert r"\begin{document}" not in supplementary
    assert r"\end{document}" not in supplementary
    assert r"\setcounter{page}" not in main + supplementary
    assert r"\pagenumbering" not in supplementary
    assert "-cd paper/supplementary/supplementary.tex" not in makefile
    assert supplementary.count(r"\input{supplementary/tables/") == 16
    assert supplementary.count(r"\includegraphics[") == 12


def test_response_framework_precedes_individual_endpoint_results() -> None:
    main = MAIN.read_text(encoding="utf-8")
    endpoint = (ROOT / "paper/endpoint_results.tex").read_text(encoding="utf-8")
    size = (ROOT / "paper/endpoint_size_results.tex").read_text(encoding="utf-8")
    assert "SolvAI is a solvent-response representation coupled to a learned hydration endpoint" in main
    assert "We evaluate the same response representation with tree ensembles" in main
    assert main.index(r"\input{endpoint_results.tex}") < main.index(r"\subsection*{Molecule alignment and matched controls}")
    assert main.index(r"\input{endpoint_size_results.tex}") > main.index(r"\subsection*{Source composition changes transfer behavior}")
    assert r"{F3_endpoint_comparison.pdf}" in endpoint
    assert r"{F8_size_and_stress.pdf}" in size
    assert "No endpoint has the lowest error in every evaluation" in main
