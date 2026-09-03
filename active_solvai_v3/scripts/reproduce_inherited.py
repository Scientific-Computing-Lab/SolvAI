#!/usr/bin/env python3
"""Recompute inherited SolvAI and Active SolvAI headline metrics from row data."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai_v3/results/inherited"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mae(frame: pd.DataFrame, method: str) -> float:
    rows = frame[frame["method"] == method]
    if rows.empty:
        raise AssertionError(f"missing method {method}")
    return float(rows["absolute_error"].mean())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sources = {
        "parent_primary": ROOT
        / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet",
        "phase1": ROOT / "active_solvai/results/phase1/phase1_endpoint_predictions.parquet",
        "phase2": ROOT / "active_solvai/results/phase2/dense_replay_predictions.parquet",
        "v2_crossfit": ROOT
        / "active_solvai/results/v2_diagnostics/oracle_independent_noise/cross_block_predictions.parquet",
        "v2_power": ROOT
        / "active_solvai/results/v2_independent_replicas/power/power_analysis.json",
        "tier_a": ROOT
        / "results/tier_a_external/evaluation/tier_a_external_predictions.parquet",
    }
    for path in sources.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    rows: list[dict[str, object]] = []

    parent = pd.read_parquet(sources["parent_primary"])
    primary = parent[parent["partition"] == "standardized_exclusion_primary"]
    rows.extend(
        [
            {"program": "SolvAI", "result": "matched_structure_mae", "value": mae(primary, "A_structure_only")},
            {"program": "SolvAI", "result": "full_solvai_mae", "value": mae(primary, "F_full_solvai")},
        ]
    )
    repeats = parent[parent["partition"] == "standardized_exclusion_repeat"]
    repeat_metrics = (
        repeats.groupby(["repeat", "method"], as_index=False)["absolute_error"].mean()
    )
    for method, label in (("A_structure_only", "matched_repeat"), ("F_full_solvai", "solvai_repeat")):
        values = repeat_metrics[repeat_metrics.method == method]["absolute_error"].to_numpy(float)
        rows.append({"program": "SolvAI", "result": f"{label}_mean", "value": float(values.mean())})
        rows.append({"program": "SolvAI", "result": f"{label}_sd", "value": float(values.std(ddof=1))})

    zero = parent[parent["partition"] == "standardized_exclusion_zero_arrow"]
    rows.extend(
        [
            {"program": "SolvAI", "result": "zero_arrow_structure_mae", "value": mae(zero, "A_structure_only")},
            {"program": "SolvAI", "result": "zero_arrow_solvai_mae", "value": mae(zero, "F_full_solvai")},
        ]
    )

    tier_a = pd.read_parquet(sources["tier_a"])
    for cohort_flag, label in (
        ("endpoint_disjoint_eligible", "tier_a_endpoint_disjoint"),
        ("strict_response_source_disjoint", "tier_a_strict_response_source_disjoint"),
    ):
        part = tier_a[tier_a[cohort_flag].astype(bool)]
        for error_column, suffix in (
            ("structure_only_absolute_error", "structure_mae"),
            ("solvai_absolute_error", "solvai_mae"),
        ):
            rows.append(
                {
                    "program": "SolvAI",
                    "result": f"{label}_{suffix}",
                    "value": float(part[error_column].mean()),
                }
            )
        rows.append(
            {
                "program": "SolvAI",
                "result": f"{label}_n",
                "value": int(part["candidate_id"].nunique()),
            }
        )

    phase1 = pd.read_parquet(sources["phase1"])
    phase1 = phase1[phase1["partition"] == "standardized_exclusion_repeat"]
    baseline = phase1[phase1.method == "P1-A_frozen_solvai"]
    candidate = phase1[
        (phase1.method == "P1-D_actual_minus_predicted")
        & (phase1.lambda_subset == "0p1_0p5_0p9")
        & (phase1.response_scope == "total")
        & np.isclose(phase1.trajectory_fraction, 1.0)
    ]
    rows.extend(
        [
            {"program": "Active_v1", "result": "frozen_solvai_repeat_mae", "value": float(baseline.absolute_error.mean())},
            {"program": "Active_v1", "result": "actual_residual_repeat_mae", "value": float(candidate.absolute_error.mean())},
        ]
    )

    phase2 = pd.read_parquet(sources["phase2"])
    for budget in (5, 7):
        part = phase2[phase2.total_windows == budget]
        for method in ("active_solvai_bq", "uniform_direct", "oracle_non_deployable"):
            values = part[part.method == method]
            if method.startswith("random"):
                values = values.groupby("molecule_id", as_index=False)["absolute_integral_error_kcal_mol"].mean()
            rows.append(
                {
                    "program": "Active_v1",
                    "result": f"phase2_{method}_mae_budget{budget}",
                    "value": float(values.absolute_integral_error_kcal_mol.mean()),
                }
            )

    crossfit = pd.read_parquet(sources["v2_crossfit"])
    molecule = (
        crossfit.groupby(["molecule_id", "method", "total_windows"], as_index=False)["absolute_error_kcal_mol"]
        .mean()
    )
    for budget in (5, 7):
        part = molecule[molecule.total_windows == budget]
        for method in ("crossfit_oracle_bq", "uniform_direct"):
            rows.append(
                {
                    "program": "Active_v2",
                    "result": f"{method}_mae_budget{budget}",
                    "value": float(part[part.method == method].absolute_error_kcal_mol.mean()),
                }
            )

    power = json.loads(sources["v2_power"].read_text())
    proposed = power["proposed_design"]
    for key in ("power_budget_5", "power_budget_7", "dense_reliability_probability"):
        rows.append({"program": "Active_v2", "result": key, "value": float(proposed[key])})

    result = pd.DataFrame(rows)
    result.to_csv(OUT / "canonical_reproduction.csv", index=False)
    payload = {
        "schema_version": 1,
        "branch_parent": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_hashes": {name: sha256(path) for name, path in sources.items()},
        "metrics": {row["result"]: row["value"] for row in rows},
    }
    (OUT / "canonical_reproduction.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
