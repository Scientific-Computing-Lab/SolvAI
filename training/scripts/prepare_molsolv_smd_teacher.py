"""Extract a leakage-safe, water-specific SMD teacher set from MolSolv.

The Zenodo archive contains 1.73 million M06-2X/6-31G* SMD(water) records.
To keep one-night training tractable, this script retains every neutral supported
molecule with at most ten heavy atoms and a deterministic one-in-eight sample of
larger (11--18 heavy atom) molecules.  Selection is independent of the ARROW
targets and structures.  All benchmark connectivities are then removed globally.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from rdkit import Chem, RDLogger

SOURCE = ROOT / "data/external/zenodo/7262826/smd_solvation_energy_dataset.sdf"
PROCESSED = ROOT / "data/processed"
ALLOWED_ATOMIC_NUMBERS = {1, 6, 7, 8, 9, 16, 17}


def identity_smiles(smiles_or_mol: str | Chem.Mol) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
        if mol is None:
            return ""
        mol = Chem.RemoveHs(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except (RuntimeError, ValueError):
        return ""


def inchi_key(smiles: str) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        value = Chem.MolToInchiKey(mol)
        return value or None
    except (RuntimeError, ValueError):
        return None


def selected_by_structure(smiles: str, heavy_atoms: int) -> bool:
    if heavy_atoms <= 10:
        return True
    digest = int.from_bytes(hashlib.sha256(smiles.encode()).digest()[:8], "little")
    return digest % 8 == 0


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    benchmark = pd.read_parquet(PROCESSED / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    public = pd.read_parquet(PROCESSED / "expanded_public_hydration_nonbenchmark.parquet")
    benchmark_by_smiles: dict[str, list[str]] = defaultdict(list)
    public_by_smiles: dict[str, list[str]] = defaultdict(list)
    for row in benchmark.itertuples(index=False):
        benchmark_by_smiles[identity_smiles(row.canonical_smiles)].append(str(row.molecule_id))
    for row in public.itertuples(index=False):
        public_by_smiles[identity_smiles(row.canonical_smiles)].append(str(row.molecule_id))

    rows: list[dict[str, object]] = []
    benchmark_matches: dict[str, list[float]] = defaultdict(list)
    public_matches: dict[str, list[float]] = defaultdict(list)
    total = valid = eligible = 0
    with SOURCE.open("rb") as stream:
        supplier = Chem.ForwardSDMolSupplier(
            stream,
            sanitize=True,
            removeHs=True,
            strictParsing=False,
        )
        for source_index, mol in enumerate(supplier):
            total += 1
            if mol is None or not mol.HasProp("solvdG"):
                continue
            try:
                target = float(mol.GetProp("solvdG"))
            except ValueError:
                continue
            if not np.isfinite(target):
                continue
            valid += 1
            if Chem.GetFormalCharge(mol) != 0 or len(Chem.GetMolFrags(mol)) != 1:
                continue
            if any(atom.GetAtomicNum() not in ALLOWED_ATOMIC_NUMBERS for atom in mol.GetAtoms()):
                continue
            heavy_atoms = mol.GetNumHeavyAtoms()
            if heavy_atoms > 18:
                continue
            smiles = identity_smiles(mol)
            if not smiles:
                continue
            eligible += 1
            for molecule_id in benchmark_by_smiles.get(smiles, []):
                benchmark_matches[molecule_id].append(target)
            for molecule_id in public_by_smiles.get(smiles, []):
                public_matches[molecule_id].append(target)
            if smiles in benchmark_by_smiles or not selected_by_structure(smiles, heavy_atoms):
                continue
            rows.append(
                {
                    "canonical_smiles": smiles,
                    "smd_water_dg": target,
                    "heavy_atom_count": heavy_atoms,
                    "source_index": source_index,
                }
            )
            if total % 100_000 == 0:
                print(
                    f"records={total} valid={valid} eligible={eligible} retained={len(rows)}",
                    flush=True,
                )

    source = pd.DataFrame(rows)
    source = source.groupby("canonical_smiles", as_index=False).agg(
        smd_water_dg=("smd_water_dg", "median"),
        smd_water_dg_std=("smd_water_dg", "std"),
        source_conformer_count=("smd_water_dg", "size"),
        heavy_atom_count=("heavy_atom_count", "first"),
        source_index=("source_index", "first"),
    )
    source["molecule_id"] = source.canonical_smiles.map(inchi_key)
    source = source.dropna(subset=["molecule_id"]).reset_index(drop=True)
    source["connectivity_key"] = source.molecule_id.str.split("-").str[0]
    benchmark_connectivity = set(benchmark.inchi_connectivity_key.astype(str))
    # A different tautomeric/protonation spelling can have a different canonical
    # SMILES while still resolving to a benchmark connectivity.  Exclude at the
    # final InChI connectivity level as a second, stricter leakage barrier.
    connectivity_overlap = source.connectivity_key.isin(benchmark_connectivity)
    connectivity_matches = source.loc[
        connectivity_overlap, ["molecule_id", "canonical_smiles", "connectivity_key"]
    ].to_dict("records")
    source = source.loc[~connectivity_overlap].reset_index(drop=True)
    if set(source.connectivity_key) & benchmark_connectivity:
        raise AssertionError("Benchmark connectivity survived MolSolv exclusion")
    source.to_parquet(PROCESSED / "molsolv_smd_water_nonbenchmark.parquet", index=False)

    exact_rows = []
    for scope, mapping in (("benchmark", benchmark_matches), ("public", public_matches)):
        for molecule_id, values in mapping.items():
            exact_rows.append(
                {
                    "prediction_scope": scope,
                    "molecule_id": molecule_id,
                    "smd_water_dg": float(np.median(values)),
                    "smd_water_dg_std": float(np.std(values)),
                    "source_conformer_count": len(values),
                }
            )
    exact = pd.DataFrame(exact_rows)
    exact.to_parquet(PROCESSED / "molsolv_smd_exact_overlaps.parquet", index=False)
    metadata = {
        "source": "MolSolv SMD water dataset",
        "doi": "10.5281/zenodo.7262826",
        "license": "CC-BY-4.0 data; GPL-2.0 code/checkpoint",
        "method": "M06-2X/6-31G* SMD(water)",
        "records_total": total,
        "records_valid": valid,
        "records_eligible_neutral_supported_heavy_le_18": eligible,
        "training_structures_after_sampling_and_deduplication": len(source),
        "benchmark_exact_structure_matches_removed": len(benchmark_matches),
        "benchmark_exact_match_ids": sorted(benchmark_matches),
        "benchmark_connectivity_alias_rows_removed": len(connectivity_matches),
        "benchmark_connectivity_aliases_removed": connectivity_matches,
        "public_exact_structure_matches": len(public_matches),
        "selection": "all heavy<=10; SHA256 one-in-eight for 11<=heavy<=18",
        "selection_uses_benchmark_structure_or_target": False,
        "benchmark_connectivity_overlap_after_filter": 0,
        "inference_simulation": False,
    }
    (PROCESSED / "molsolv_smd_teacher_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
