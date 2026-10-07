"""Keep best-result emphasis aligned with the displayed, comparable errors."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "paper"


def _rows(path: str, data_row) -> list[list[str]]:
    rows = []
    for line in (TABLES / path).read_text().splitlines():
        if " & " not in line or not line.endswith(r"\\"):
            continue
        cells = [part.strip() for part in line[:-2].split(" & ")]
        if data_row(cells):
            rows.append(cells)
    assert rows, path
    return rows


def _score(cell: str) -> float | None:
    plain = re.sub(r"\\(?:textbf|mathbf|ensuremath)", "", cell)
    plain = plain.replace("{", "").replace("}", "").replace("$", "")
    return float(plain) if re.fullmatch(r"-?\d+\.\d+", plain) else None


def _check(path: str, data_row, columns: tuple[int, ...], group: int | None = None) -> None:
    rows = _rows(path, data_row)
    by_group = defaultdict(list)
    for row in rows:
        by_group[row[group] if group is not None else "all"].append(row)
    for cohort in by_group.values():
        for column in columns:
            values = [(row[column], _score(row[column])) for row in cohort]
            values = [(cell, score) for cell, score in values if score is not None]
            assert values, (path, column)
            minimum = min(score for _, score in values)
            for cell, score in values:
                assert (r"\textbf{" in cell or r"\mathbf{" in cell) == (score == minimum), (
                    path, column, cell, minimum
                )


def test_main_endpoint_table_bolds_lowest_full_response_mae_by_cohort():
    endpoints = {
        "ExtraTrees", "Residual MLP", "Dual branch", "Ridge + neural residual",
        "Message passing", "MoLFormer + neural head", "TabM", "TabPFN-3.5",
        "MoLFormer + ExtraTrees",
    }
    _check("tables/endpoint_summary.tex", lambda cells: len(cells) == 6 and cells[0] in endpoints,
           (2, 3, 4, 5))


def test_supplementary_matched_comparisons_bold_only_within_matching_set():
    _check("supplementary/tables/repeat_values.tex", lambda cells: len(cells) == 4 and cells[0].isdigit(), (3,), 0)
    _check("supplementary/tables/global_separation.tex", lambda cells: len(cells) == 4 and cells[2].isdigit(), (3,), 0)
    _check("supplementary/tables/weight_one_sensitivity.tex", lambda cells: len(cells) == 5 and cells[1].isdigit(), (2, 3, 4))
    _check("supplementary/tables/tier_a_external.tex", lambda cells: len(cells) == 6 and cells[2].isdigit(), (3, 4, 5), 0)


def test_supplementary_endpoint_tables_bold_descriptive_minima():
    cohort = lambda cells: cells[0] in {"ARROW-85", "External-220", "Strict-97"}
    _check("supplementary/tables/journal_baselines.tex", cohort, (3, 4, 5), 0)
    _check("supplementary/tables/endpoint_metrics.tex", cohort, (3, 4, 5, 6), 0)
    _check("supplementary/tables/endpoint_comparisons.tex", cohort, (2,), 0)
    _check("supplementary/tables/endpoint_size.tex", lambda cells: len(cells) == 4 and _score(cells[1]) is not None, (1, 2, 3))
