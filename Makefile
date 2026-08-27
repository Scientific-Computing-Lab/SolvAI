.PHONY: setup test verify metrics tables figures paper security claims clean

setup:
	uv sync --extra dev

test:
	uv run pytest -q

metrics:
	uv run python scripts/reproduce_metrics.py

verify: metrics
	uv run python scripts/audit_leakage.py
	uv run python scripts/verify_artifact.py
	uv run python scripts/check_claims.py
	uv run python scripts/security_scan.py
	uv run python scripts/verify_release_manifest.py

tables: metrics
	uv run python scripts/make_tables.py

figures: tables
	uv run python scripts/make_figures.py

paper: figures
	SOURCE_DATE_EPOCH=1787788800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/main.tex
	SOURCE_DATE_EPOCH=1787788800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/supplementary/supplementary.tex
	SOURCE_DATE_EPOCH=1787788800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=ED_Table1 -cd paper/extended_data/ED_Table1_standalone.tex
	qpdf --empty --static-id --pages paper/main.pdf paper/extended_data/ED_Fig1_residuals.pdf paper/extended_data/ED_Fig2_provenance.pdf paper/extended_data/ED_Fig3_alternatives.pdf paper/extended_data/ED_Fig4_selective_pimd.pdf paper/extended_data/ED_Fig5_lambda_response.pdf paper/extended_data/ED_Fig6_extrapolation.pdf paper/extended_data/ED_Fig7_statistics.pdf paper/extended_data/ED_Table1.pdf -- paper/review_combined.pdf

security:
	uv run python scripts/security_scan.py

claims:
	uv run python scripts/check_claims.py

clean:
	latexmk -C -cd paper/main.tex
	latexmk -C -cd paper/supplementary/supplementary.tex
	latexmk -C -jobname=ED_Table1 -cd paper/extended_data/ED_Table1_standalone.tex
