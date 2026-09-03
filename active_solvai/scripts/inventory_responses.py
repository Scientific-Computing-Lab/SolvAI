#!/usr/bin/env python3
"""Inventory existing response observations without running simulation."""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from active_solvai.probes import energy_group, prefix_summary, read_energy

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = REPO.parent.parent
ACTIVE = REPO / "active_solvai"
PROBE_ROOT = WORKSPACE / "results/physics_probes"
CAMPAIGNS = {
    "lambda01": (PROBE_ROOT / "pimd2_lambda01_5ps", 0.1),
    "lambda05": (PROBE_ROOT / "pimd2_lambda05_5ps", 0.5),
    "lambda09": (PROBE_ROOT / "pimd2_lambda09_5ps", 0.9),
}
FRACTIONS = (0.10, 0.20, 0.40, 0.70, 1.00)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def succeeded(case: Path) -> bool:
    console = case / "console.txt"
    return console.is_file() and "Simulation has successfully finished!" in console.read_text(
        errors="ignore"
    )


def main() -> None:
    benchmark = pd.read_parquet(REPO / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark["solvent"].eq("water")].copy()
    identity = benchmark[
        [
            "molecule_id",
            "molecule_name",
            "canonical_smiles",
            "inchi_key",
            "inchi_connectivity_key",
            "functional_group_family",
            "scaffold",
        ]
    ].drop_duplicates("molecule_id")
    name_to_identity = identity.set_index("molecule_name").to_dict(orient="index")

    cases: list[dict[str, object]] = []
    files: list[dict[str, object]] = []
    blocks: list[dict[str, object]] = []
    for campaign, (root, lambda_value) in CAMPAIGNS.items():
        manifest = pd.read_csv(root / "manifest.csv")
        for manifest_row in manifest.to_dict(orient="records"):
            name = str(manifest_row["molecule_name"])
            case = Path(str(manifest_row["config"])).parent
            row: dict[str, object] = {
                "campaign": campaign,
                "lambda": lambda_value,
                "molecule_name": name,
                "case_directory": str(case),
                "success": succeeded(case),
                "beads": int(manifest_row["beads"]),
                "ti_point": int(manifest_row["ti_point"]),
                "production_steps": int(manifest_row["production_steps"]),
                "production_ps": int(manifest_row["production_steps"]) * 0.002,
            }
            row.update(name_to_identity.get(name, {}))
            timing = case / "timing.json"
            if timing.exists():
                timing_data = json.loads(timing.read_text())
                row["elapsed_seconds"] = timing_data.get("elapsed_seconds")
                row["exit_status"] = timing_data.get("exit_status")
                row["timed_out"] = timing_data.get("timed_out")
            energy_paths = sorted((case / "output").glob("*.ene"))
            row["energy_file_count"] = len(energy_paths)
            cases.append(row)

            candidate_files = [
                case / "conf.xml",
                case / "console.txt",
                case / "timing.json",
                *energy_paths,
            ]
            for path in candidate_files:
                if not path.is_file():
                    continue
                files.append(
                    {
                        "campaign": campaign,
                        "lambda": lambda_value,
                        "molecule_name": name,
                        "molecule_id": row.get("molecule_id"),
                        "role": "energy" if path.suffix == ".ene" else path.name,
                        "energy_group": energy_group(path) if path.suffix == ".ene" else None,
                        "path": str(path),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
                if path.suffix != ".ene" or not row["success"]:
                    continue
                frame = read_energy(path)
                for summary in prefix_summary(frame, FRACTIONS):
                    blocks.append(
                        {
                            "campaign": campaign,
                            "lambda": lambda_value,
                            "molecule_name": name,
                            "molecule_id": row.get("molecule_id"),
                            "functional_group_family": row.get("functional_group_family"),
                            "energy_group": energy_group(path),
                            **summary,
                        }
                    )

    case_frame = pd.DataFrame(cases)
    file_frame = pd.DataFrame(files)
    block_frame = pd.DataFrame(blocks)
    manifests = ACTIVE / "data/manifests"
    results = ACTIVE / "results/phase0"
    manifests.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    case_frame.to_csv(manifests / "response_case_inventory.csv", index=False)
    file_frame.to_csv(manifests / "response_file_manifest.csv", index=False)
    block_frame.to_parquet(results / "response_prefix_blocks.parquet", index=False)

    historical_wide = WORKSPACE / "results/pimd2_multilambda_teacher.parquet"
    historical_long = WORKSPACE / "results/pimd2_multilambda_teacher_long.parquet"
    wide = pd.read_parquet(historical_wide)
    long = pd.read_parquet(historical_long)
    complete_ids = set(wide.loc[wide["complete_three_point_curve"], "molecule_id"].astype(str))
    if len(complete_ids) != 72:
        raise AssertionError(f"Expected 72 complete three-point curves, found {len(complete_ids)}")
    if case_frame["molecule_id"].isna().any():
        missing = case_frame.loc[case_frame["molecule_id"].isna(), "molecule_name"].unique()
        raise AssertionError(f"Probe cases lack benchmark identities: {missing.tolist()}")

    availability = (
        case_frame.pivot_table(
            index="molecule_id",
            columns="campaign",
            values="success",
            aggfunc="max",
            fill_value=False,
        )
        .rename(columns=lambda value: f"has_probe_{value}")
        .reset_index()
    )
    identity_manifest = identity.merge(availability, on="molecule_id", how="left")
    for campaign in CAMPAIGNS:
        column = f"has_probe_{campaign}"
        identity_manifest[column] = identity_manifest[column].fillna(False).astype(bool)
    identity_manifest["complete_three_point_curve"] = (
        identity_manifest["molecule_id"].astype(str).isin(complete_ids)
    )
    fold_columns = [
        "molecule_id",
        "fold_random",
        "fold_family",
        "fold_scaffold",
        "delta_g_exp",
        "delta_g_classical_arrow",
        "delta_g_pimd4",
        "delta_g_pimd8",
    ]
    identity_manifest = identity_manifest.merge(
        benchmark[fold_columns].drop_duplicates("molecule_id"), on="molecule_id", how="left"
    )
    identity_manifest["active_exposure_tier"] = "ARROW-development"
    identity_dir = ACTIVE / "data/identity"
    identity_dir.mkdir(parents=True, exist_ok=True)
    identity_manifest.to_parquet(identity_dir / "probe_identity_manifest.parquet", index=False)
    identity_manifest.to_csv(identity_dir / "probe_identity_manifest.csv", index=False)

    repeat_path = REPO / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet"
    repeated = pd.read_parquet(repeat_path)
    split_assignments = (
        repeated.loc[
            repeated["method"].eq("F_full_solvai"),
            ["molecule_id", "partition", "repeat", "split_seed", "fold"],
        ]
        .drop_duplicates()
        .merge(identity_manifest[["molecule_id", "complete_three_point_curve"]], on="molecule_id")
    )
    split_assignments.to_csv(manifests / "probe_split_assignments.csv", index=False)

    # Dense local data are suitable for software tests or single-molecule diagnostics,
    # not population replay. Duplicated recovered-repository copies are counted once.
    dense_candidates = [
        {
            "dataset": "toluene_11_window_test",
            "path": str(
                WORKSPACE
                / "repositories/Simulations/MD_TOOLS_PY/test/test_runscripts/test_analyze/data/286227.InWater_C07H08_Toluene/Output"
            ),
            "molecules": 1,
            "windows": 11,
            "population_replay_eligible": False,
            "reason": "single software-test molecule; protocol compatibility requires audit",
        },
        {
            "dataset": "legacy_21_window_ene_temp",
            "path": str(WORKSPACE / "repositories_recovered/scripts/arbalest/ene-temp"),
            "molecules": None,
            "windows": 21,
            "population_replay_eligible": False,
            "reason": "molecular identity and protocol provenance unresolved",
        },
    ]
    (manifests / "dense_response_candidates.json").write_text(
        json.dumps(dense_candidates, indent=2) + "\n"
    )
    summary = {
        "generated_utc": datetime.now(UTC).isoformat(),
        "campaigns": {
            campaign: {
                "lambda": lambda_value,
                "attempted": int((case_frame["campaign"] == campaign).sum()),
                "successful": int(
                    case_frame.loc[case_frame["campaign"].eq(campaign), "success"].sum()
                ),
            }
            for campaign, (_, lambda_value) in CAMPAIGNS.items()
        },
        "unique_probe_molecules": int(case_frame["molecule_id"].nunique()),
        "complete_three_point_molecules": len(complete_ids),
        "raw_files_hashed": len(file_frame),
        "raw_bytes_hashed": int(file_frame["bytes"].sum()),
        "prefix_summary_rows": len(block_frame),
        "identity_manifest_rows": len(identity_manifest),
        "trajectory_fractions": list(FRACTIONS),
        "historical_derived": {
            "wide_path": str(historical_wide),
            "wide_sha256": sha256(historical_wide),
            "wide_shape": list(wide.shape),
            "long_path": str(historical_long),
            "long_sha256": sha256(historical_long),
            "long_shape": list(long.shape),
        },
        "retrospective_replay": {
            "three_point_endpoint_gate": True,
            "trajectory_prefix_gate": True,
            "population_dense_hidden_window_replay": False,
            "reason": "No compatible dense multi-window population was located in authorized local storage.",
        },
    }
    (results / "response_inventory_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
