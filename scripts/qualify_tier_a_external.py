#!/usr/bin/env python3
"""Qualify the pre-existing Tier-A cohort without evaluating SolvAI predictions.

This script is deliberately limited to target compatibility, molecular scope,
identity equivalence and source exposure. It does not load endpoint models,
teacher predictions or Tier-A prediction errors.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

RELEASE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = RELEASE_ROOT.parents[1]
ASP_ROOT = WORKSPACE_ROOT.parent / "ASP"
SOURCE_COHORT = ASP_ROOT / "data/processed/tier_a_frozen.csv"
SOURCE_IDENTITIES = (
    WORKSPACE_ROOT
    / "release/large_confirmatory_intermediates/supervision_standardized_identities.parquet"
)
BENCHMARK_IDENTITIES = RELEASE_ROOT / "audits/confirmatory/benchmark_standardized_identities.parquet"
OUT = RELEASE_ROOT / "results/tier_a_external/qualification"

EXPECTED_HASHES = {
    SOURCE_COHORT: "7b7be03c5124559551c5380b0429d5c0156c13060af9902d26266510637fcf45",
    SOURCE_IDENTITIES: "6d35c5ea58192d8fe3d0b5c022effcfa93fc8c76cdbd19ed424d0b7b50d2b85f",
}
KNOWN_REGISTRY_CONFLICT = "MKERQGLKSFEKAE"
ALLOWED_ATOMIC_NUMBERS = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
IDENTITY_COLUMNS = (
    "full_inchi_key",
    "connectivity_key",
    "fragment_parent_key",
    "uncharged_parent_key",
    "canonical_tautomer_key",
)
TEACHER_SOURCES = (
    "combisolv_qm",
    "abraham",
    "openff",
    "gbn2",
    "molsolv_smd",
    "confsolv",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable_path(path: Path) -> str:
    """Describe an input without embedding a release-host home directory."""
    for prefix, root in (
        ("repository", RELEASE_ROOT),
        ("workspace", WORKSPACE_ROOT),
        ("external/ASP", ASP_ROOT),
    ):
        try:
            return f"{prefix}/{path.relative_to(root)}"
        except ValueError:
            continue
    return path.name


def connectivity(key: str) -> str:
    return str(key).split("-")[0] if key else ""


def identity(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return {"parse_status": "parse_error"}
    fragment = rdMolStandardize.FragmentParent(mol)
    uncharged = rdMolStandardize.Uncharger().uncharge(fragment)
    tautomer = rdMolStandardize.TautomerEnumerator().Canonicalize(uncharged)
    full_key = Chem.MolToInchiKey(mol)
    return {
        "parse_status": "ok",
        "canonical_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "full_inchi_key": full_key,
        "connectivity_key": connectivity(full_key),
        "fragment_parent_key": connectivity(Chem.MolToInchiKey(fragment)),
        "uncharged_parent_key": connectivity(Chem.MolToInchiKey(uncharged)),
        "canonical_tautomer_key": connectivity(Chem.MolToInchiKey(tautomer)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "fragment_count": len(Chem.GetMolFrags(mol)),
        "supported_elements": all(
            atom.GetAtomicNum() in ALLOWED_ATOMIC_NUMBERS for atom in mol.GetAtoms()
        ),
    }


def match_records(
    candidates: pd.DataFrame, source: pd.DataFrame, source_name: str
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    lookups: dict[str, dict[str, list[int]]] = {}
    for column in IDENTITY_COLUMNS:
        mapping: dict[str, list[int]] = {}
        for index, value in source[column].fillna("").astype(str).items():
            if value:
                mapping.setdefault(value, []).append(index)
        lookups[column] = mapping

    exact_flags: list[bool] = []
    standardized_flags: list[bool] = []
    details: list[dict[str, object]] = []
    for row in candidates.itertuples():
        matched_exact = False
        matched_standardized = False
        seen: set[tuple[str, int]] = set()
        for column in IDENTITY_COLUMNS:
            value = str(getattr(row, column, "") or "")
            if not value:
                continue
            for source_index in lookups[column].get(value, []):
                token = (column, int(source_index))
                if token in seen:
                    continue
                seen.add(token)
                exact = column in {"full_inchi_key", "connectivity_key"}
                matched_exact |= exact
                matched_standardized |= not exact
                source_row = source.loc[source_index]
                details.append(
                    {
                        "candidate_id": row.candidate_id,
                        "candidate_name": row.name,
                        "candidate_smiles": row.smiles,
                        "source": source_name,
                        "match_type": column,
                        "match_value": value,
                        "source_row": int(source_row.get("source_row", source_index)),
                        "source_smiles": source_row.get("source_smiles", ""),
                        "source_connectivity_key": source_row.get("connectivity_key", ""),
                    }
                )
        exact_flags.append(matched_exact)
        standardized_flags.append(matched_standardized)
    result = pd.DataFrame(
        {
            f"overlap_{source_name}_exact": exact_flags,
            f"overlap_{source_name}_standardized": standardized_flags,
        }
    )
    return result, details


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path, expected in EXPECTED_HASHES.items():
        observed = sha256(path)
        if observed != expected:
            raise AssertionError(f"Input hash changed for {path}: {observed}")

    cohort = pd.read_csv(SOURCE_COHORT)
    if len(cohort) != 221:
        raise AssertionError(f"Expected original Tier-A N=221, found {len(cohort)}")
    identities = pd.DataFrame([identity(value) for value in cohort.smiles.astype(str)])
    candidates = pd.concat([cohort.reset_index(drop=True), identities], axis=1)

    compatible = (
        candidates.parse_status.eq("ok")
        & candidates.source.eq("Sander Henry database v5.0.0")
        & candidates.source_type.eq("M")
        & np.isclose(candidates.temperature_k.astype(float), 298.15)
        & candidates.hscp_units.eq("mol m^-3 Pa^-1")
        & candidates.standard_state.eq("1 M ideal gas to 1 M ideal-dilute aqueous")
        & candidates.sign_convention.eq("negative is favorable hydration")
        & np.isfinite(candidates.experimental_dg_kcal_mol.astype(float))
        & candidates.formal_charge.eq(0)
        & candidates.fragment_count.eq(1)
        & candidates.supported_elements
    )
    candidates["target_scope_compatible"] = compatible
    candidates["known_registry_conflict"] = candidates.candidate_id.eq(KNOWN_REGISTRY_CONFLICT)
    candidates["scientifically_eligible"] = compatible & ~candidates.known_registry_conflict

    source_identities = pd.read_parquet(SOURCE_IDENTITIES).reset_index(drop=True)
    benchmark_identities = pd.read_parquet(BENCHMARK_IDENTITIES).reset_index(drop=True)
    benchmark_source = benchmark_identities.copy()
    benchmark_source["source_row"] = np.arange(len(benchmark_source))
    benchmark_source["source_smiles"] = benchmark_source.source_smiles.astype(str)

    all_match_details: list[dict[str, object]] = []
    endpoint_public = source_identities.loc[
        source_identities.source.eq("endpoint_experimental")
    ].reset_index(drop=True)
    flags, details = match_records(candidates, endpoint_public, "endpoint_public_1280")
    candidates = pd.concat([candidates, flags], axis=1)
    all_match_details.extend(details)
    flags, details = match_records(candidates, benchmark_source, "endpoint_arrow_85")
    candidates = pd.concat([candidates, flags], axis=1)
    all_match_details.extend(details)

    for source_name in TEACHER_SOURCES:
        source = source_identities.loc[source_identities.source.eq(source_name)].reset_index(
            drop=True
        )
        flags, details = match_records(candidates, source, f"teacher_{source_name}")
        candidates = pd.concat([candidates, flags], axis=1)
        all_match_details.extend(details)

    endpoint_columns = [
        "overlap_endpoint_public_1280_exact",
        "overlap_endpoint_public_1280_standardized",
        "overlap_endpoint_arrow_85_exact",
        "overlap_endpoint_arrow_85_standardized",
    ]
    teacher_columns = [
        column
        for column in candidates
        if column.startswith("overlap_teacher_")
    ]
    candidates["any_endpoint_overlap"] = candidates[endpoint_columns].any(axis=1)
    candidates["any_teacher_source_overlap"] = candidates[teacher_columns].any(axis=1)
    candidates["endpoint_disjoint_eligible"] = (
        candidates.scientifically_eligible & ~candidates.any_endpoint_overlap
    )
    candidates["strict_response_source_disjoint"] = (
        candidates.endpoint_disjoint_eligible & ~candidates.any_teacher_source_overlap
    )

    reason_columns = []
    for row in candidates.itertuples():
        reasons: list[str] = []
        if not row.target_scope_compatible:
            reasons.append("target_or_neutral_scope_incompatible")
        if row.known_registry_conflict:
            reasons.append("preexisting_registry_name_structure_conflict")
        if row.any_endpoint_overlap:
            reasons.append("experimental_endpoint_identity_overlap")
        reason_columns.append(";".join(reasons))
    candidates["endpoint_cohort_exclusion_reason"] = reason_columns

    output_columns = [
        "candidate_id",
        "name",
        "smiles",
        "canonical_isomeric_smiles",
        "full_inchi_key",
        "connectivity_key",
        "fragment_parent_key",
        "uncharged_parent_key",
        "canonical_tautomer_key",
        "experimental_dg_kcal_mol",
        "temperature_k",
        "hscp_units",
        "standard_state",
        "sign_convention",
        "heavy_atom_count",
        "molecular_weight",
        "n_eligible_measurement_rows",
        "n_independent_source_refs",
        "source_record_ids",
        "source_refs",
        "source_dois",
        "target_scope_compatible",
        "known_registry_conflict",
        "scientifically_eligible",
        *endpoint_columns,
        *teacher_columns,
        "any_endpoint_overlap",
        "any_teacher_source_overlap",
        "endpoint_disjoint_eligible",
        "strict_response_source_disjoint",
        "endpoint_cohort_exclusion_reason",
    ]
    audit = candidates[output_columns].copy()
    audit.to_csv(OUT / "tier_a_qualification_audit.csv", index=False)
    audit.to_parquet(OUT / "tier_a_qualification_audit.parquet", index=False)
    pd.DataFrame(all_match_details).to_csv(OUT / "tier_a_identity_matches.csv", index=False)
    audit.loc[~audit.endpoint_disjoint_eligible].to_csv(
        OUT / "tier_a_exclusions.csv", index=False
    )
    endpoint = audit.loc[audit.endpoint_disjoint_eligible].reset_index(drop=True)
    strict = audit.loc[audit.strict_response_source_disjoint].reset_index(drop=True)
    endpoint.to_csv(OUT / "tier_a_endpoint_disjoint.csv", index=False)
    endpoint.to_parquet(OUT / "tier_a_endpoint_disjoint.parquet", index=False)
    strict.to_csv(OUT / "tier_a_strict_response_source_disjoint.csv", index=False)
    strict.to_parquet(OUT / "tier_a_strict_response_source_disjoint.parquet", index=False)

    source_exposure = []
    for source_name in TEACHER_SOURCES:
        exact = f"overlap_teacher_{source_name}_exact"
        standardized = f"overlap_teacher_{source_name}_standardized"
        source_exposure.append(
            {
                "teacher_source": source_name,
                "endpoint_disjoint_rows": len(endpoint),
                "exact_overlap_rows": int(endpoint[exact].sum()),
                "standardized_overlap_rows": int(endpoint[standardized].sum()),
                "any_overlap_rows": int(endpoint[[exact, standardized]].any(axis=1).sum()),
            }
        )
    exposure_frame = pd.DataFrame(source_exposure)
    exposure_frame.to_csv(OUT / "tier_a_teacher_exposure_summary.csv", index=False)

    files = [
        "tier_a_qualification_audit.csv",
        "tier_a_qualification_audit.parquet",
        "tier_a_identity_matches.csv",
        "tier_a_exclusions.csv",
        "tier_a_endpoint_disjoint.csv",
        "tier_a_endpoint_disjoint.parquet",
        "tier_a_strict_response_source_disjoint.csv",
        "tier_a_strict_response_source_disjoint.parquet",
        "tier_a_teacher_exposure_summary.csv",
    ]
    summary = {
        "qualification_only": True,
        "prediction_errors_loaded": False,
        "original_cohort_rows": len(cohort),
        "target_scope_compatible_rows": int(candidates.target_scope_compatible.sum()),
        "scientifically_eligible_rows": int(candidates.scientifically_eligible.sum()),
        "known_registry_conflicts": int(candidates.known_registry_conflict.sum()),
        "endpoint_overlap_rows": int(
            (candidates.scientifically_eligible & candidates.any_endpoint_overlap).sum()
        ),
        "endpoint_disjoint_eligible_rows": len(endpoint),
        "teacher_exposed_endpoint_disjoint_rows": int(endpoint.any_teacher_source_overlap.sum()),
        "strict_response_source_disjoint_rows": len(strict),
        "source_exposure": source_exposure,
        "rdkit_version": rdBase.rdkitVersion,
        "inputs": {
            "tier_a_source": portable_path(SOURCE_COHORT),
            "tier_a_source_sha256": sha256(SOURCE_COHORT),
            "source_identity_index": portable_path(SOURCE_IDENTITIES),
            "source_identity_index_sha256": sha256(SOURCE_IDENTITIES),
            "benchmark_identity_index": portable_path(BENCHMARK_IDENTITIES),
            "benchmark_identity_index_sha256": sha256(BENCHMARK_IDENTITIES),
        },
        "outputs_sha256": {name: sha256(OUT / name) for name in files},
    }
    (OUT / "tier_a_qualification_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
