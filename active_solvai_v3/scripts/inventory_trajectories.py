#!/usr/bin/env python3
"""Inventory inherited trajectories without interpreting them as independent replicas."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
WORKSPACE = ROOT.parents[1]
OUT = ROOT / "active_solvai_v3/data/manifests"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_energy_metadata(path: Path) -> dict[str, object]:
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [str(column).strip() for column in frame.columns]
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    times = numeric["Time"].dropna().to_numpy(float) if "Time" in numeric else np.array([])
    return {
        "frames": len(frame),
        "observed_time_ps": float(times.max() - times.min()) if len(times) else np.nan,
        "has_dhdl": bool("dHdL" in numeric and numeric["dHdL"].notna().any()),
        "finite_dhdl": bool("dHdL" in numeric and np.isfinite(numeric["dHdL"]).all()),
    }


def parse_param(text: str, title: str) -> str | None:
    pattern = rf'<Param\s+Title="{re.escape(title)}"[^>]*>(.*?)</Param>'
    match = re.search(pattern, text, flags=re.DOTALL)
    return match.group(1).strip() if match else None


def infer_config(energy: Path) -> Path | None:
    case = energy.parent.parent if energy.parent.name.lower().startswith("output") else energy.parent
    direct = case / "conf.xml"
    if direct.is_file():
        return direct
    candidates = sorted(case.glob("conf*.xml"))
    if not candidates:
        candidates = sorted(energy.parent.glob("conf*.xml"))
    match = re.search(r"Annihilation_(\d+)", energy.name)
    if match:
        index = int(match.group(1))
        exact = [path for path in candidates if path.stem in {f"conf{index}", f"conf{index:02d}"}]
        if exact:
            return exact[0]
    return candidates[0] if len(candidates) == 1 else None


def config_metadata(path: Path | None) -> dict[str, object]:
    if path is None:
        return {
            "config": None,
            "config_sha256": None,
            "beads": np.nan,
            "ti_point": np.nan,
            "lambda": np.nan,
            "time_step_ps": np.nan,
            "random_seed_fields": None,
        }
    text = path.read_text(errors="replace")
    bead_raw = parse_param(text, "NumPIMDReplicas")
    point_raw = parse_param(text, "TIPoint") or parse_param(text, "LambdaPointNumber")
    lambda_raw = parse_param(text, "LambdaValues")
    step_raw = parse_param(text, "TimeStep")
    seeds = re.findall(r'<Param\s+Title="[^"]*Seed[^"]*"[^>]*>(.*?)</Param>', text, re.DOTALL)
    ti_point = int(point_raw) if point_raw and point_raw.strip().isdigit() else np.nan
    lambdas = [float(value) for value in lambda_raw.split()] if lambda_raw else []
    lambda_value = lambdas[int(ti_point)] if lambdas and np.isfinite(ti_point) else np.nan
    return {
        "config": str(path),
        "config_sha256": sha256(path),
        "beads": int(bead_raw) if bead_raw and bead_raw.strip().isdigit() else 1,
        "ti_point": ti_point,
        "lambda": lambda_value,
        "time_step_ps": float(step_raw) if step_raw else np.nan,
        "random_seed_fields": ";".join(seed.strip() for seed in seeds) or None,
    }


def source_class(path: Path) -> str:
    value = str(path)
    if "/results/physics_probes/pimd2_lambda" in value:
        return "historical_three_point_pimd2"
    if "/active_solvai/simulations/dense_pimd2/" in value:
        return "active_v1_dense_pimd2"
    if "/results/rtx3090_benchmark/" in value:
        return "rtx3090_protocol_benchmark"
    if "/repositories/Simulations/" in value or "/repositories_recovered/Simulations/" in value:
        return "repository_test_fixture"
    if "/repositories_recovered/scripts/arbalest/ene-temp/" in value:
        return "legacy_unqualified_output"
    return "ancillary_unqualified_output"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    probe_cases = pd.read_csv(ROOT / "active_solvai/data/manifests/response_case_inventory.csv")
    probe_files = pd.read_csv(ROOT / "active_solvai/data/manifests/response_file_manifest.csv")
    probe_success = probe_cases[["campaign", "lambda", "molecule_id", "success"]]
    probe_system = probe_files[
        (probe_files["role"] == "energy") & (probe_files["energy_group"] == "system")
    ].merge(probe_success, on=["campaign", "lambda", "molecule_id"], how="left")
    successful_probe_paths = set(probe_system.loc[probe_system.success.fillna(False), "path"])
    system_files = sorted(
        path
        for path in WORKSPACE.rglob("*.ene")
        if path.name.endswith("SYSTEM.ene") and ".git" not in path.parts and ".venv" not in path.parts
    )
    rows: list[dict[str, object]] = []
    for path in system_files:
        source = source_class(path)
        config = infer_config(path)
        metadata = config_metadata(config)
        energy = read_energy_metadata(path)
        value = str(path)
        eligible = source == "active_v1_dense_pimd2" or (
            source == "historical_three_point_pimd2" and value in successful_probe_paths
        )
        # Active-v1 dense paths cover only the twelve newly generated lambda windows;
        # the inherited 0.1/0.5/0.9 files remain in the historical source tree.
        rows.append(
            {
                "source_class": source,
                "path": value,
                "energy_sha256": sha256(path),
                **metadata,
                **energy,
                "stream_id": "single_inherited_stream",
                "independence_status": "one trajectory; contiguous blocks are not replicas",
                "development_exposed": True,
                "eligible_v3_development": eligible,
                "eligible_v3_prospective_sentinel": False,
            }
        )

    frame = pd.DataFrame(rows).sort_values(["source_class", "path"]).reset_index(drop=True)
    frame.to_csv(OUT / "trajectory_inventory.csv", index=False)
    frame.to_parquet(OUT / "trajectory_inventory.parquet", index=False)

    dense_manifest = pd.read_csv(ROOT / "active_solvai/simulations/dense_pimd2/manifest.csv")
    parent = pd.read_parquet(ROOT / "data/benchmark/arrow_solvation_master.parquet")
    summary = {
        "schema_version": 1,
        "all_system_energy_files": len(frame),
        "source_counts": frame.groupby("source_class").size().astype(int).to_dict(),
        "qualified_pimd2_system_files": int(frame.eligible_v3_development.sum()),
        "independent_replica_sets": 0,
        "historical_probe_attempts": len(probe_cases),
        "historical_probe_successes": int(probe_cases.success.sum()),
        "historical_probe_unique_molecules": int(probe_cases.molecule_id.nunique()),
        "historical_probe_complete_three_point_molecules": int(
            probe_cases[probe_cases.success]
            .groupby("molecule_id")["lambda"]
            .nunique()
            .eq(3)
            .sum()
        ),
        "dense_manifest_rows": len(dense_manifest),
        "dense_molecules": int(dense_manifest.molecule_id.nunique()),
        "dense_lambda_values": sorted(dense_manifest["lambda"].unique().tolist()),
        "dense_existing_three_point_rows": int(dense_manifest.existing_observation.sum()),
        "dense_new_window_rows": int((~dense_manifest.existing_observation).sum()),
        "dense_streams_per_molecule_window": 1,
        "dense_production_ps_per_window": sorted(dense_manifest.production_ps.unique().tolist()),
        "parent_scalar_rows": len(parent),
        "parent_scalar_fields": [
            column for column in parent.columns if "classical" in column.lower() or "pimd" in column.lower()
        ],
        "scientific_limit": (
            "The inherited full-grid library contains one 5-ps trajectory per molecule-window. "
            "It permits complementary-block development diagnostics but no independent-replica "
            "policy validation or prospective sentinel claim."
        ),
    }
    (OUT / "trajectory_inventory_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
