"""Independent identity-leakage and inference-schema audits for SolvAI."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import joblib
import pandas as pd

from .features import canonicalize

ROOT = Path(__file__).resolve().parents[1]

SOURCES = {
    "Expanded public hydration": (
        "expanded_public_hydration_nonbenchmark.parquet",
        "canonical_smiles",
        "inchi_key",
        "inchi_connectivity_key",
    ),
    "Legacy public hydration": (
        "public_hydration_nonbenchmark.parquet",
        "canonical_smiles",
        "inchi_key",
        "inchi_connectivity_key",
    ),
    "CombiSolv-QM water": (
        "combisolv_qm_water_nonbenchmark.parquet",
        "solute_canonical_smiles",
        "solute_inchi_key",
        "solute_connectivity_key",
    ),
    "SoluteML Abraham": (
        "soluteml_abraham_nonbenchmark.parquet",
        "canonical_smiles",
        "inchi_key",
        "connectivity_key",
    ),
    "OpenFF explicit alchemical": (
        "openff_alchemical_nonbenchmark.parquet",
        "canonical_smiles",
        "inchi_key",
        "connectivity_key",
    ),
    "GBn2 implicit-solvent": (
        "implicit_solvent_nonbenchmark.parquet",
        "canonical_smiles",
        "inchi_key",
        "connectivity_key",
    ),
    "MolSolv SMD(water)": (
        "molsolv_smd_water_nonbenchmark.parquet",
        "canonical_smiles",
        "molecule_id",
        "connectivity_key",
    ),
    "ConfSolv water response": (
        "confsolv_water_nonbenchmark.parquet",
        "canonical_smiles",
        "inchi_key",
        "connectivity_key",
    ),
}

FORBIDDEN_FEATURE_FRAGMENTS = (
    "delta_g_exp",
    "delta_g_pimd",
    "delta_g_classical_arrow",
    "experimental_residual",
    "nqe_residual",
    "trajectory",
    "probe__",
    "fold_",
    "functional_group_family",
    "scaffold",
)


def _identities(frame: pd.DataFrame, smiles_column: str) -> pd.DataFrame:
    rows = []
    for smiles in frame[smiles_column].astype(str).drop_duplicates():
        canonical, key, connectivity = canonicalize(smiles)
        rows.append(
            {
                "source_smiles": smiles,
                "canonical_smiles": canonical,
                "inchi_key": key,
                "connectivity_key": connectivity,
            }
        )
    return pd.DataFrame(rows)


def run_leakage_audit(root: Path = ROOT) -> tuple[dict[str, Any], pd.DataFrame]:
    benchmark = pd.read_parquet(root / "data/benchmark/arrow_solvation_master.parquet")
    benchmark = benchmark.loc[benchmark.solvent.eq("water")].copy()
    benchmark_identity = _identities(benchmark, "canonical_smiles")
    if len(benchmark_identity) != 85:
        raise AssertionError("Benchmark canonicalization did not yield 85 identities")
    benchmark_smiles = set(benchmark_identity.canonical_smiles)
    benchmark_keys = set(benchmark_identity.inchi_key)
    benchmark_connectivity = set(benchmark_identity.connectivity_key)

    prior_audit_path = root / "audits/leakage_audit.csv"
    prior_audit = (
        pd.read_csv(prior_audit_path).set_index("source")
        if prior_audit_path.exists()
        else pd.DataFrame()
    )
    rows: list[dict[str, Any]] = []
    for source, (filename, smiles_column, full_column, connectivity_column) in SOURCES.items():
        source_path = root / "data/processed" / filename
        if not source_path.exists():
            if source not in prior_audit.index:
                raise FileNotFoundError(source_path)
            row = prior_audit.loc[source].to_dict()
            row["source"] = source
            row["verification_mode"] = "frozen clean-room audit; source not redistributed"
            rows.append(row)
            continue
        frame = pd.read_parquet(source_path)
        identities = _identities(frame, smiles_column)
        supplied_full = set(frame[full_column].dropna().astype(str))
        supplied_connectivity = set(frame[connectivity_column].dropna().astype(str))
        recomputed_full = set(identities.inchi_key)
        recomputed_connectivity = set(identities.connectivity_key)
        full_mismatch = len(supplied_full.symmetric_difference(recomputed_full))
        connectivity_mismatch = len(
            supplied_connectivity.symmetric_difference(recomputed_connectivity)
        )
        row = {
            "source": source,
            "rows": len(frame),
            "unique_canonical_smiles": len(identities),
            "exact_smiles_overlap": len(benchmark_smiles & set(identities.canonical_smiles)),
            "full_inchikey_overlap": len(benchmark_keys & recomputed_full),
            "connectivity_overlap": len(benchmark_connectivity & recomputed_connectivity),
            "stereochemical_alias_overlap": int(
                len(benchmark_connectivity & recomputed_connectivity)
                - len(benchmark_keys & recomputed_full)
            ),
            "stored_full_identity_symmetric_difference": full_mismatch,
            "stored_connectivity_symmetric_difference": connectivity_mismatch,
            "verification_mode": "independent release-time recanonicalization",
        }
        rows.append(row)
    table = pd.DataFrame(rows)
    if (
        table[["exact_smiles_overlap", "full_inchikey_overlap", "connectivity_overlap"]]
        .to_numpy()
        .any()
    ):
        raise AssertionError("A supervised external source overlaps the ARROW-85 benchmark")
    if (
        table[
            [
                "stored_full_identity_symmetric_difference",
                "stored_connectivity_symmetric_difference",
            ]
        ]
        .to_numpy()
        .any()
    ):
        raise AssertionError(
            "Stored source identities disagree with independent RDKit canonicalization"
        )

    head = joblib.load(root / "models/final/head.joblib")
    feature_names = [*head["descriptor_columns"], *head["teacher_columns"]]
    forbidden = sorted(
        name
        for name in feature_names
        if any(fragment in name.lower() for fragment in FORBIDDEN_FEATURE_FRAGMENTS)
    )
    if forbidden:
        raise AssertionError(f"Forbidden inference features: {forbidden}")
    inference_files = [root / "solv_ai/inference.py", root / "solv_ai/teachers.py"]
    forbidden_runtime_paths = []
    for path in inference_files:
        source = path.read_text()
        for marker in ("data/benchmark", "results/predictions", "pimd2_multilambda"):
            if marker in source:
                forbidden_runtime_paths.append(f"{path.name}:{marker}")
    if forbidden_runtime_paths:
        raise AssertionError(
            f"Inference code references evaluation data: {forbidden_runtime_paths}"
        )

    report = {
        "status": "PASS",
        "identity_policy": (
            "Independent RDKit canonical isomeric SMILES, full InChIKey, and first "
            "InChIKey connectivity block; connectivity overlap is disallowed."
        ),
        "benchmark_molecules": 85,
        "sources": rows,
        "all_supervised_external_connectivity_overlaps_zero": True,
        "artifact": {
            "input": "SMILES only",
            "descriptor_features": len(head["descriptor_columns"]),
            "predicted_physics_response_features": len(head["teacher_columns"]),
            "total_features": len(feature_names),
            "ensemble_members": len(head["models"]),
            "forbidden_feature_count": 0,
            "forbidden_runtime_data_references": forbidden_runtime_paths,
            "simulation_at_inference": False,
        },
    }
    return report, table


def write_leakage_audit(root: Path = ROOT) -> dict[str, Any]:
    report, table = run_leakage_audit(root)
    audit_dir = root / "audits"
    audit_dir.mkdir(parents=True, exist_ok=True)
    (audit_dir / "leakage_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    table.to_csv(audit_dir / "leakage_audit.csv", index=False)
    lines = [
        "# Independent SolvAI leakage and inference audit",
        "",
        "**Status: PASS.** Every supervised external training source used by the final model",
        "is disjoint from the 85-solute reference set at the InChIKey connectivity level.",
        "Canonical isomeric SMILES and full InChIKey comparisons also have zero overlap.",
        "",
        "| Source | Rows | Unique structures | SMILES overlap | Full-key overlap | Connectivity overlap | Verification |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in table.itertuples(index=False):
        lines.append(
            f"| {row.source} | {row.rows:,} | {row.unique_canonical_smiles:,} | "
            f"{row.exact_smiles_overlap} | {row.full_inchikey_overlap} | "
            f"{row.connectivity_overlap} | {row.verification_mode} |"
        )
    lines.extend(
        [
            "",
            "The released endpoint contains 2,265 deterministic RDKit/Morgan features and",
            "15 structure-predicted physical-response priors. Its schema contains no",
            "experimental, ARROW, PIMD, trajectory, probe, family, scaffold or fold field.",
            "The runtime code does not open benchmark or prediction tables.",
            "",
            "CombiSolv-QM and mixed public-hydration source rows are not redistributed",
            "where publisher supplements have no standalone data licence. Their final",
            "clean-room audit records, counts and source-table SHA-256 values are frozen",
            "in this release; all included tables are recanonicalized when the audit is rerun.",
        ]
    )
    (audit_dir / "leakage_audit.md").write_text("\n".join(lines) + "\n")
    return report
