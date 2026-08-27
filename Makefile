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
	latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/main.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/supplementary.tex

security:
	uv run python scripts/security_scan.py

claims:
	uv run python scripts/check_claims.py

clean:
	latexmk -C -cd paper/main.tex
	latexmk -C -cd paper/supplementary.tex
