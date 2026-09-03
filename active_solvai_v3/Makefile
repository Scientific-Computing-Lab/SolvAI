.PHONY: setup test reproduce figures metrics paper verify clean

setup:
	uv sync --project . --extra dev

test:
	PYTHONPATH=src uv run --project . pytest -q tests
	PYTHONPATH=src uv run --project . ruff check src scripts tests

reproduce:
	uv run --project ../active_solvai python scripts/reproduce_inherited.py
	uv run --project ../active_solvai python scripts/inventory_trajectories.py
	PYTHONPATH=src uv run --project .. python scripts/build_gate1_dataset.py
	PYTHONPATH=src uv run --project . python scripts/run_gate1_identifiability.py
	PYTHONPATH=src uv run --project . python scripts/run_power_reliability_gate.py

figures:
	PYTHONPATH=src uv run --project . python scripts/make_gate1_figure.py

metrics:
	PYTHONPATH=src uv run --project . python scripts/finalize_v3.py >/dev/null

paper: metrics figures
	cd paper && SOURCE_DATE_EPOCH=1788393600 FORCE_SOURCE_DATE=1 latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex
	cd paper/supplementary && SOURCE_DATE_EPOCH=1788393600 FORCE_SOURCE_DATE=1 latexmk -g -pdf -interaction=nonstopmode -halt-on-error supplementary.tex

verify: test metrics figures paper
	PYTHONPATH=src uv run --project . python scripts/verify_v3_release.py
	PYTHONPATH=src uv run --project . python scripts/build_release_manifest.py

clean:
	cd paper && latexmk -C main.tex
	cd paper/supplementary && latexmk -C supplementary.tex
