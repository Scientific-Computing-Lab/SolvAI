.PHONY: setup test verify metrics tables figures paper journal-paper clean

setup:
	uv sync --extra dev

test:
	uv run pytest -q

metrics:
	uv run python scripts/reproduce_metrics.py

verify: metrics
	uv run python scripts/audit_leakage.py
	uv run python scripts/verify_artifact.py
	uv run python scripts/build_release_manifest.py
	uv run python scripts/verify_release_manifest.py

tables: metrics
	uv run python scripts/make_tables.py

figures: tables
	uv run python scripts/make_figures.py

paper: figures
	SOURCE_DATE_EPOCH=1787788800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/main.tex

clean:
	latexmk -C -cd paper/main.tex

# One canonical PDF: main text, references and Supplementary Information.
# Compile from verified frozen figures and tables, without retraining.
journal-paper:
	latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/main.tex
