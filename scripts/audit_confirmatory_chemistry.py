#!/usr/bin/env python3
"""Audit exact, standardized and chemical-distance overlap for final supervision."""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
TAUTOMER_ENUMERATOR = None
UNCHARGER = None


def initialize_worker() -> None:
    global TAUTOMER_ENUMERATOR, UNCHARGER
    TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()
    UNCHARGER = rdMolStandardize.Uncharger()


def connectivity(full_key: str) -> str:
    return full_key.split("-")[0] if full_key else ""


def identity_record(smiles: str):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return (str(smiles), "", "", "", "", "", "", "", b"", "parse_error")
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        full_key = Chem.MolToInchiKey(mol)
        fragment = rdMolStandardize.FragmentParent(mol)
        fragment_key = connectivity(Chem.MolToInchiKey(fragment))
        uncharged = UNCHARGER.uncharge(fragment)
        uncharged_key = connectivity(Chem.MolToInchiKey(uncharged))
        tautomer = TAUTOMER_ENUMERATOR.Canonicalize(uncharged)
        tautomer_key = connectivity(Chem.MolToInchiKey(tautomer))
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
        fingerprint = DataStructs.BitVectToBinaryText(FP_GENERATOR.GetFingerprint(mol))
        return (
            str(smiles),
            canonical,
            full_key,
            connectivity(full_key),
            fragment_key,
            uncharged_key,
            tautomer_key,
            scaffold,
            fingerprint,
            "ok",
        )
    except Exception as exc:  # noqa: BLE001  # pragma: no cover - audit records it
        return (str(smiles), "", "", "", "", "", "", "", b"", type(exc).__name__)


def process_smiles(smiles: list[str], processes: int) -> pd.DataFrame:
    if processes == 1:
        initialize_worker()
        records = [identity_record(value) for value in smiles]
    else:
        with mp.get_context("spawn").Pool(
            processes=processes, initializer=initialize_worker
        ) as pool:
            records = list(pool.imap(identity_record, smiles, chunksize=256))
    return pd.DataFrame(
        records,
        columns=[
            "source_smiles",
            "canonical_isomeric_smiles",
            "full_inchi_key",
            "connectivity_key",
            "fragment_parent_key",
            "uncharged_parent_key",
            "canonical_tautomer_key",
            "murcko_scaffold",
            "fingerprint_binary",
            "parse_status",
        ],
    )


def source_specs(processed: Path, endpoint: pd.DataFrame):
    return (
        ("endpoint_experimental", endpoint, "canonical_smiles", "inchi_connectivity_key"),
        (
            "combisolv_qm",
            pd.read_parquet(processed / "combisolv_qm_water_nonbenchmark.parquet"),
            "solute_canonical_smiles",
            "solute_connectivity_key",
        ),
        (
            "abraham",
            pd.read_parquet(processed / "soluteml_abraham_nonbenchmark.parquet"),
            "canonical_smiles",
            "connectivity_key",
        ),
        (
            "openff",
            pd.read_parquet(processed / "openff_alchemical_nonbenchmark.parquet"),
            "canonical_smiles",
            "connectivity_key",
        ),
        (
            "gbn2",
            pd.read_parquet(processed / "implicit_solvent_nonbenchmark.parquet"),
            "canonical_smiles",
            "connectivity_key",
        ),
        (
            "molsolv_smd",
            pd.read_parquet(processed / "molsolv_smd_water_nonbenchmark.parquet"),
            "canonical_smiles",
            "connectivity_key",
        ),
        (
            "confsolv",
            pd.read_parquet(processed / "confsolv_water_nonbenchmark.parquet"),
            "canonical_smiles",
            "connectivity_key",
        ),
    )


def set_lookup(frame: pd.DataFrame, column: str) -> dict[str, list[int]]:
    lookup: dict[str, list[int]] = {}
    for index, value in enumerate(frame[column].astype(str)):
        if value:
            lookup.setdefault(value, []).append(index)
    return lookup


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--release-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--processes", type=int, default=min(20, mp.cpu_count()))
    args = parser.parse_args()
    processed = args.workspace_root / "data" / "processed"
    out = args.release_root / "audits" / "confirmatory"
    out.mkdir(parents=True, exist_ok=True)

    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = (
        benchmark.loc[benchmark.solvent.eq("water")]
        .drop_duplicates("molecule_id")
        .reset_index(drop=True)
    )
    public_all = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    public_old = pd.read_parquet(processed / "public_hydration_nonbenchmark.parquet")
    old_keys = set(public_old.inchi_connectivity_key.astype(str))
    endpoint = public_all.loc[
        public_all.inchi_connectivity_key.astype(str).isin(old_keys)
        | public_all.source_measurement_count.fillna(0).ge(2)
    ].reset_index(drop=True)
    if len(benchmark) != 85 or len(endpoint) != 1280:
        raise AssertionError("Frozen benchmark or endpoint pool changed")

    benchmark_identity = process_smiles(
        benchmark.canonical_smiles.astype(str).tolist(), min(args.processes, 8)
    )
    benchmark_identity.insert(0, "molecule_id", benchmark.molecule_id.astype(str))
    benchmark_identity.insert(1, "molecule_name", benchmark.molecule_name.astype(str))
    benchmark_identity.to_parquet(out / "benchmark_standardized_identities.parquet", index=False)
    benchmark_fps = [
        DataStructs.CreateFromBinaryText(value) for value in benchmark_identity.fingerprint_binary
    ]
    benchmark_lookups = {
        column: set_lookup(benchmark_identity, column)
        for column in (
            "full_inchi_key",
            "connectivity_key",
            "fragment_parent_key",
            "uncharged_parent_key",
            "canonical_tautomer_key",
        )
    }
    benchmark_scaffolds = set(benchmark_identity.murcko_scaffold) - {""}

    summary_rows: list[dict] = []
    distance_rows: list[dict] = []
    match_rows: list[dict] = []
    for source, table, smiles_column, provided_key_column in source_specs(processed, endpoint):
        print(f"auditing {source}: {len(table):,} rows", flush=True)
        unique = (
            table[[smiles_column, provided_key_column]]
            .dropna(subset=[smiles_column])
            .drop_duplicates(smiles_column)
            .reset_index(drop=True)
        )
        identities = process_smiles(unique[smiles_column].astype(str).tolist(), args.processes)
        identities.insert(0, "source", source)
        identities.insert(1, "source_row", np.arange(len(identities)))
        identities.insert(
            2,
            "provided_connectivity_key",
            unique[provided_key_column].astype(str).to_numpy(),
        )
        parse_failures = int((~identities.parse_status.eq("ok")).sum())
        source_fps = [
            DataStructs.CreateFromBinaryText(value)
            for value in identities.loc[identities.parse_status.eq("ok"), "fingerprint_binary"]
        ]
        valid_indices = identities.index[identities.parse_status.eq("ok")].to_numpy()
        maximum_similarities = np.zeros(len(benchmark), dtype=float)
        nearest_indices = np.full(len(benchmark), -1, dtype=int)
        for benchmark_index, benchmark_fp in enumerate(benchmark_fps):
            similarities = np.asarray(
                DataStructs.BulkTanimotoSimilarity(benchmark_fp, source_fps),
                dtype=float,
            )
            nearest_local = int(np.argmax(similarities))
            nearest_index = int(valid_indices[nearest_local])
            maximum_similarities[benchmark_index] = similarities[nearest_local]
            nearest_indices[benchmark_index] = nearest_index
            nearest = identities.loc[nearest_index]
            distance_rows.append(
                {
                    "source": source,
                    "benchmark_molecule_id": benchmark.loc[benchmark_index, "molecule_id"],
                    "benchmark_molecule_name": benchmark.loc[benchmark_index, "molecule_name"],
                    "benchmark_smiles": benchmark.loc[benchmark_index, "canonical_smiles"],
                    "maximum_morgan_tanimoto": similarities[nearest_local],
                    "nearest_source_smiles": nearest.source_smiles,
                    "nearest_source_connectivity_key": nearest.connectivity_key,
                    "same_murcko_scaffold": bool(
                        nearest.murcko_scaffold
                        and nearest.murcko_scaffold
                        == benchmark_identity.loc[benchmark_index, "murcko_scaffold"]
                    ),
                }
            )
            if similarities[nearest_local] == 1.0 and (
                nearest.connectivity_key
                != benchmark_identity.loc[benchmark_index, "connectivity_key"]
            ):
                match_rows.append(
                    {
                        "source": source,
                        "match_type": "morgan_similarity_1_nonidentical",
                        "benchmark_molecule_id": benchmark.loc[benchmark_index, "molecule_id"],
                        "benchmark_smiles": benchmark.loc[benchmark_index, "canonical_smiles"],
                        "source_smiles": nearest.source_smiles,
                        "source_connectivity_key": nearest.connectivity_key,
                        "value": 1.0,
                    }
                )

        overlap_counts = {}
        for column, lookup in benchmark_lookups.items():
            count = 0
            for source_index, key in enumerate(identities[column].astype(str)):
                if key and key in lookup:
                    count += 1
                    for benchmark_index in lookup[key]:
                        match_rows.append(
                            {
                                "source": source,
                                "match_type": column,
                                "benchmark_molecule_id": benchmark.loc[
                                    benchmark_index, "molecule_id"
                                ],
                                "benchmark_smiles": benchmark.loc[
                                    benchmark_index, "canonical_smiles"
                                ],
                                "source_smiles": identities.loc[source_index, "source_smiles"],
                                "source_connectivity_key": identities.loc[
                                    source_index, "connectivity_key"
                                ],
                                "value": key,
                            }
                        )
            overlap_counts[column] = count
        scaffold_match_rows = int(identities.murcko_scaffold.isin(benchmark_scaffolds).sum())
        summary_rows.append(
            {
                "source": source,
                "source_rows": len(table),
                "unique_smiles_audited": len(unique),
                "parse_failures": parse_failures,
                "full_inchi_key_matches": overlap_counts["full_inchi_key"],
                "connectivity_matches": overlap_counts["connectivity_key"],
                "fragment_parent_matches": overlap_counts["fragment_parent_key"],
                "uncharged_parent_matches": overlap_counts["uncharged_parent_key"],
                "canonical_tautomer_matches": overlap_counts["canonical_tautomer_key"],
                "rows_sharing_benchmark_scaffold": scaffold_match_rows,
                "maximum_similarity": float(maximum_similarities.max()),
                "median_of_benchmark_max_similarity": float(np.median(maximum_similarities)),
            }
        )
        print(
            source,
            summary_rows[-1],
            flush=True,
        )

    summary = pd.DataFrame(summary_rows)
    distances = pd.DataFrame(distance_rows)
    matches = pd.DataFrame(match_rows)
    summary.to_csv(out / "chemical_distance_summary.csv", index=False)
    distances.to_csv(out / "chemical_distance_by_benchmark.csv", index=False)
    matches.to_csv(out / "chemical_identity_matches.csv", index=False)
    excluded = matches.loc[
        matches.match_type.isin(
            [
                "fragment_parent_key",
                "uncharged_parent_key",
                "canonical_tautomer_key",
            ]
        ),
        ["source", "source_smiles", "source_connectivity_key"],
    ].drop_duplicates()
    excluded["source_smiles_sha256"] = excluded.source_smiles.map(
        lambda value: hashlib.sha256(value.encode()).hexdigest()
    )
    excluded.to_csv(out / "standardized_exclusion_records.csv", index=False)

    exact_clean = bool(
        summary.full_inchi_key_matches.eq(0).all() and summary.connectivity_matches.eq(0).all()
    )
    standardized_clean = bool(
        summary.fragment_parent_matches.eq(0).all()
        and summary.uncharged_parent_matches.eq(0).all()
        and summary.canonical_tautomer_matches.eq(0).all()
    )
    audit = {
        "benchmark_rows": len(benchmark),
        "endpoint_rows": len(endpoint),
        "sources": summary.to_dict(orient="records"),
        "exact_identity_clean": exact_clean,
        "standardized_identity_clean": standardized_clean,
        "requires_standardized_exclusion_rerun": not standardized_clean,
        "fingerprint": "Morgan radius 2, 2,048 bits",
        "standardization": [
            "RDKit parse and canonical isomeric SMILES",
            "MolStandardize FragmentParent",
            "MolStandardize Uncharger",
            "MolStandardize canonical tautomer",
        ],
    }
    (out / "chemical_distance_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    lines = [
        "# Confirmatory chemical-distance and leakage audit",
        "",
        f"- Benchmark molecules: **{len(benchmark)}**",
        f"- Endpoint experimental training rows: **{len(endpoint):,}**",
        f"- Exact identity/connectivity firewall clean: **{'YES' if exact_clean else 'NO'}**",
        f"- Standardized parent/tautomer firewall clean: **{'YES' if standardized_clean else 'NO'}**",
        "",
        "| Source | Rows | Exact connectivity | Uncharged parent | Canonical tautomer | Shared scaffold rows | Maximum similarity |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary.itertuples():
        lines.append(
            f"| {row.source} | {row.source_rows:,} | {row.connectivity_matches} | "
            f"{row.uncharged_parent_matches} | {row.canonical_tautomer_matches} | "
            f"{row.rows_sharing_benchmark_scaffold:,} | {row.maximum_similarity:.3f} |"
        )
    lines.extend(
        [
            "",
            "The scaffold and similarity results are proximity diagnostics, not identity leakage. Every identity-equivalent match is listed in `chemical_identity_matches.csv`; nearest neighbours for every benchmark/source pair are listed in `chemical_distance_by_benchmark.csv`.",
            "",
        ]
    )
    (out / "chemical_distance_audit.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
