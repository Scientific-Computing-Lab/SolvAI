"""Dataset construction and static molecular/ARROW feature extraction."""

from __future__ import annotations

import json
import math
import re
import shlex
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold, KFold

from .structures import BENCHMARK_SMILES, HIN_BY_NAME

ROOT = Path(__file__).resolve().parents[2]
PAPER_XLSX = (
    ROOT / "data/raw/freecurve_paper/PMC8776904_supplementary/41467_2022_28041_MOESM3_ESM.xlsx"
)
FREESOLV_JSON = ROOT / "data/raw/public/freesolv_database.json"
COMBISOLV_EXP = ROOT / "data/raw/public/CombiSolv-Exp-8780.csv"
QMPFF_ATOM = (
    ROOT / "repositories/FF-NN-2/src/nn_tools/resources/Input/QMPFF/933/output/QMPFFATOM.PAR"
)
ARROW_NN_DIR = ROOT / "repositories/ARROW_NN_PARAM"
HYDRATION_PIMD4 = {"Water": -6.20, "Ethanol": -5.18}


def canonicalize(smiles: str) -> tuple[str, str, str, bool]:
    """Return canonical SMILES, InChIKey, connectivity key, and stereo flag."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    key = Chem.MolToInchiKey(mol)
    stereo = any(a.HasProp("_CIPCode") for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )
    return canonical, key, key.split("-")[0], stereo


def _scaffold(mol: Chem.Mol, family: str) -> str:
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    return scaffold or f"ACYCLIC::{family.lower()}"


def basic_descriptors(mol: Chem.Mol) -> dict[str, float | int]:
    return {
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "molecular_weight": float(Descriptors.MolWt(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
    }


def rdkit_descriptor_frame(records: pd.DataFrame) -> pd.DataFrame:
    """Compute all stable scalar RDKit descriptors plus Morgan fingerprints."""
    descriptor_list = Descriptors._descList
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    rows: list[dict[str, Any]] = []
    for row in records.itertuples(index=False):
        mol = Chem.MolFromSmiles(row.canonical_smiles)
        assert mol is not None
        out: dict[str, Any] = {
            "molecule_id": row.molecule_id,
            "canonical_smiles": row.canonical_smiles,
        }
        for name, function in descriptor_list:
            try:
                value = float(function(mol))
                out[f"rdkit__{name}"] = value if math.isfinite(value) else np.nan
            except (ValueError, OverflowError, RuntimeError, ZeroDivisionError):
                out[f"rdkit__{name}"] = np.nan
        fp = fpgen.GetFingerprintAsNumPy(mol)
        out.update({f"morgan2__{i:04d}": int(v) for i, v in enumerate(fp)})
        rows.append(out)
    return pd.DataFrame(rows)


def load_benchmark() -> pd.DataFrame:
    raw = pd.read_excel(PAPER_XLSX, sheet_name="Solvation_Data", header=2)
    family: str | None = None
    records: list[dict[str, Any]] = []
    for _, row in raw.iterrows():
        name = str(row.iloc[0]).strip()
        exp = pd.to_numeric(row.iloc[1], errors="coerce")
        if pd.isna(exp):
            family = name
            continue
        if name not in BENCHMARK_SMILES:
            raise KeyError(f"No curated structure for benchmark molecule {name!r}")
        assert family is not None
        smiles, key, connectivity, stereo = canonicalize(BENCHMARK_SMILES[name])
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        common = {
            "molecule_id": key,
            "canonical_smiles": smiles,
            "inchi_key": key,
            "inchi_connectivity_key": connectivity,
            "molecule_name": name,
            "temperature": 298.0,
            "source_repo": "freecurve/interx_solvation_suite (paper-linked)",
            "source_file": (
                "data/raw/freecurve_paper/PMC8776904_supplementary/"
                "41467_2022_28041_MOESM3_ESM.xlsx#Solvation_Data"
            ),
            "source_commit": None,
            "source_publication": "doi:10.1038/s41467-022-28041-0",
            "functional_group_family": family,
            "scaffold": _scaffold(mol, family),
            "arrow_parameterization_supported": True,
            "stereochemistry_specified": stereo,
            **basic_descriptors(mol),
        }
        hydration = {
            **common,
            "condition_id": f"{key}::water::298K",
            "solvent": "water",
            "delta_g_exp": float(exp),
            "delta_g_exp_uncertainty": np.nan,
            "delta_g_classical_arrow": float(row.iloc[2]),
            "delta_g_classical_uncertainty": np.nan,
            "delta_g_pimd4": HYDRATION_PIMD4.get(name, np.nan),
            "delta_g_pimd8": float(row.iloc[3]),
            "delta_g_pimd8_uncertainty": np.nan,
            "nqe_residual": float(row.iloc[3] - row.iloc[2]),
            "experimental_residual": float(row.iloc[1] - row.iloc[3]),
        }
        if name in HYDRATION_PIMD4:
            hydration["source_file"] += ";41467_2022_28041_MOESM1_ESM.pdf#Supplementary_Table_4"
        records.append(hydration)
        cyclohexane_exp = pd.to_numeric(row.iloc[6], errors="coerce")
        if pd.notna(cyclohexane_exp):
            records.append(
                {
                    **common,
                    "condition_id": f"{key}::cyclohexane::298K",
                    "solvent": "cyclohexane",
                    "delta_g_exp": float(cyclohexane_exp),
                    "delta_g_exp_uncertainty": np.nan,
                    "delta_g_classical_arrow": float(row.iloc[7]),
                    "delta_g_classical_uncertainty": np.nan,
                    "delta_g_pimd4": float(row.iloc[8]),
                    "delta_g_pimd8": np.nan,
                    "delta_g_pimd8_uncertainty": np.nan,
                    "nqe_residual": np.nan,
                    "experimental_residual": np.nan,
                }
            )
    result = pd.DataFrame(records)
    hydration = result[result.solvent == "water"]
    if len(hydration) != 85 or not hydration.molecule_id.is_unique:
        raise AssertionError(f"Expected 85 unique hydration molecules, got {len(hydration)}")
    if result.duplicated(["molecule_id", "solvent", "temperature"]).any():
        raise AssertionError("Molecule-condition rows are not unique")
    return result


def load_freesolv() -> pd.DataFrame:
    payload = json.loads(FREESOLV_JSON.read_text())
    records: list[dict[str, Any]] = []
    for identifier, item in payload.items():
        try:
            smiles, key, connectivity, stereo = canonicalize(item["smiles"])
        except (ValueError, RuntimeError):
            continue
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        records.append(
            {
                "public_id": identifier,
                "molecule_id": key,
                "canonical_smiles": smiles,
                "inchi_key": key,
                "inchi_connectivity_key": connectivity,
                "molecule_name": str(
                    item.get("iupac") or item.get("nickname") or identifier
                ).strip(),
                "delta_g_exp": float(item["expt"]),
                "delta_g_exp_uncertainty": float(item["d_expt"]),
                "formal_charge": int(Chem.GetFormalCharge(mol)),
                "fragment_count": len(Chem.GetMolFrags(mol)),
                "stereochemistry_specified": stereo,
                "source_publication": "FreeSolv v0.52; doi:10.1021/acs.jced.7b00104",
                "source_file": "data/raw/public/freesolv_database.json",
                **basic_descriptors(mol),
            }
        )
    frame = pd.DataFrame(records)
    # The source has unique IDs but can contain stereoisomers or aliases sharing
    # connectivity. Preserve all here; fold-specific code uses connectivity keys.
    return frame.sort_values("public_id").reset_index(drop=True)


def load_combisolv_water() -> pd.DataFrame:
    """Load the water-solvent subset of the public CombiSolv-Exp-8780 copy."""
    raw = pd.read_csv(COMBISOLV_EXP)
    pairs = raw.ssid.str.split(".", n=1, expand=True, regex=False)
    raw = raw.loc[pairs[0].eq("O")].copy()
    raw["solute_smiles"] = pairs.loc[pairs[0].eq("O"), 1]
    records: list[dict[str, Any]] = []
    for index, row in raw.iterrows():
        try:
            smiles, key, connectivity, stereo = canonicalize(str(row.solute_smiles))
        except (ValueError, RuntimeError):
            continue
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        records.append(
            {
                "public_id": f"combisolv_exp_{index}",
                "molecule_id": key,
                "canonical_smiles": smiles,
                "inchi_key": key,
                "inchi_connectivity_key": connectivity,
                "molecule_name": None,
                "delta_g_exp": float(row.dgsolv),
                "delta_g_exp_uncertainty": np.nan,
                "formal_charge": int(Chem.GetFormalCharge(mol)),
                "fragment_count": len(Chem.GetMolFrags(mol)),
                "stereochemistry_specified": stereo,
                "source_publication": "CombiSolv-Exp-8780 via su-group/SolvBERT",
                "source_file": "data/raw/public/CombiSolv-Exp-8780.csv",
                **basic_descriptors(mol),
            }
        )
    return pd.DataFrame(records)


def combine_public_hydration(
    freesolv: pd.DataFrame, combisolv: pd.DataFrame, benchmark: pd.DataFrame
) -> pd.DataFrame:
    """Make a deduplicated, benchmark-free experimental hydration training set."""
    source = pd.concat([freesolv, combisolv], ignore_index=True, sort=False)
    benchmark_connectivity = set(benchmark.inchi_connectivity_key)
    source = source[
        (source.formal_charge == 0)
        & (source.fragment_count == 1)
        & (source.heavy_atom_count > 0)
        & (~source.inchi_connectivity_key.isin(benchmark_connectivity))
    ].copy()
    rows: list[dict[str, Any]] = []
    # Collapse stereochemical aliases by connectivity because the benchmark is
    # mostly stereo-unspecified and hydration labels are not consistently
    # stereospecific across the contributing historical databases.
    for _, group in source.groupby("inchi_connectivity_key", sort=True):
        representative = group.iloc[0].to_dict()
        values = group.delta_g_exp.astype(float).to_numpy()
        representative["delta_g_exp"] = float(np.median(values))
        representative["public_label_count"] = len(values)
        representative["public_label_range"] = float(values.max() - values.min())
        uncertainties = pd.to_numeric(group.delta_g_exp_uncertainty, errors="coerce")
        representative["delta_g_exp_uncertainty"] = (
            float(np.nanmedian(uncertainties)) if uncertainties.notna().any() else np.nan
        )
        representative["public_sources"] = ";".join(sorted(set(group.source_file.astype(str))))
        representative["public_ids"] = ";".join(sorted(group.public_id.astype(str)))
        rows.append(representative)
    return pd.DataFrame(rows).reset_index(drop=True)


def attach_freesolv_metadata(benchmark: pd.DataFrame, freesolv: pd.DataFrame) -> pd.DataFrame:
    """Attach uncertainty/reference metadata without replacing paper targets."""
    result = benchmark.copy()
    by_connectivity: dict[str, pd.DataFrame] = {
        key: group for key, group in freesolv.groupby("inchi_connectivity_key", sort=False)
    }
    matches: list[str | None] = []
    uncertainties: list[float] = []
    match_kinds: list[str] = []
    external_values: list[float] = []
    for row in result.itertuples():
        if row.solvent != "water":
            matches.append(None)
            uncertainties.append(np.nan)
            match_kinds.append("not_hydration")
            external_values.append(np.nan)
            continue
        candidates = by_connectivity.get(row.inchi_connectivity_key)
        if candidates is None:
            matches.append(None)
            uncertainties.append(np.nan)
            match_kinds.append("none")
            external_values.append(np.nan)
            continue
        exact = candidates[candidates.inchi_key == row.inchi_key]
        pool = exact if len(exact) else candidates
        chosen = pool.iloc[(pool.delta_g_exp - row.delta_g_exp).abs().argmin()]
        matches.append(str(chosen.public_id))
        uncertainties.append(float(chosen.delta_g_exp_uncertainty))
        match_kinds.append("exact_inchikey" if len(exact) else "connectivity_only")
        external_values.append(float(chosen.delta_g_exp))
    result["freesolv_match_id"] = matches
    result["freesolv_match_method"] = match_kinds
    result["freesolv_delta_g_exp"] = external_values
    result["delta_g_exp_uncertainty"] = uncertainties
    return result


def _load_qmpff_atom_parameters(path: Path = QMPFF_ATOM) -> dict[str, dict[str, float]]:
    header: list[str] | None = None
    result: dict[str, dict[str, float]] = {}
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith(";#!"):
            header = stripped[3:].split()
            continue
        if not stripped or stripped.startswith(";") or header is None:
            continue
        parts = stripped.split()
        if len(parts) < len(header):
            continue
        values = dict(zip(header, parts, strict=False))
        atom_type = values["ATOM"]
        result[atom_type] = {
            key: float(value) for key, value in values.items() if key not in {"ID", "ATOM"}
        }
    if not result:
        raise ValueError(f"No atom parameters parsed from {path}")
    return result


def _parse_hin(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    in_first_molecule = False
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("mol "):
            if in_first_molecule:
                break
            in_first_molecule = True
        elif stripped.startswith("endmol") and in_first_molecule:
            break
        elif stripped.startswith("atom ") and in_first_molecule:
            parts = shlex.split(stripped)
            if len(parts) < 11:
                continue
            try:
                atoms.append(
                    {
                        "element": parts[3],
                        "atom_type": parts[4],
                        "charge": float(parts[6]),
                        "x": float(parts[7]),
                        "y": float(parts[8]),
                        "z": float(parts[9]),
                    }
                )
            except ValueError:
                continue
    return atoms


def _nn_supported_atom_types() -> set[str]:
    types: set[str] = set()
    for path in ARROW_NN_DIR.glob("*.pb"):
        left, _, right = path.stem.partition("_")
        if right in {"OW", "HW"}:
            types.add(left)
    return types


def arrow_static_features(benchmark: pd.DataFrame) -> pd.DataFrame:
    params = _load_qmpff_atom_parameters()
    nn_types = _nn_supported_atom_types()
    preliminary: list[dict[str, Any]] = []
    all_types: set[str] = set()
    for row in benchmark.itertuples(index=False):
        relpath = HIN_BY_NAME.get(row.molecule_name)
        out: dict[str, Any] = {
            "molecule_id": row.molecule_id,
            "canonical_smiles": row.canonical_smiles,
            "arrow_hin_path": relpath,
            "arrow_static_available": False,
        }
        if relpath is None or not (ROOT / relpath).is_file():
            preliminary.append(out)
            continue
        atoms = _parse_hin(ROOT / relpath)
        if not atoms:
            preliminary.append(out)
            continue
        out["arrow_static_available"] = True
        counts = Counter(atom["atom_type"] for atom in atoms)
        all_types.update(counts)
        out["_type_counts"] = counts
        out["arrow_atom_count"] = len(atoms)
        out["arrow_unique_atom_types"] = len(counts)
        charges = np.asarray([atom["charge"] for atom in atoms], dtype=float)
        out["arrow_charge_sum"] = float(charges.sum())
        out["arrow_abs_charge_sum"] = float(np.abs(charges).sum())
        out["arrow_charge_std"] = float(charges.std())
        coords = np.asarray([[atom["x"], atom["y"], atom["z"]] for atom in atoms], dtype=float)
        centered = coords - coords.mean(axis=0)
        out["arrow_radius_gyration"] = float(np.sqrt(np.mean(np.sum(centered**2, axis=1))))
        if len(coords) > 1:
            distances = np.sqrt(np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=2))
            out["arrow_max_interatomic_distance"] = float(distances.max())
        else:
            out["arrow_max_interatomic_distance"] = 0.0
        weights = np.asarray(list(counts.values()), dtype=float)
        probabilities = weights / weights.sum()
        out["arrow_atom_type_entropy"] = float(-np.sum(probabilities * np.log(probabilities)))
        out["arrow_nn_supported_atom_fraction"] = float(
            np.mean([atom["atom_type"] in nn_types for atom in atoms])
        )
        parameter_rows = [
            params[atom["atom_type"]] for atom in atoms if atom["atom_type"] in params
        ]
        out["arrow_parameter_match_fraction"] = len(parameter_rows) / len(atoms)
        if parameter_rows:
            for key in next(iter(parameter_rows)):
                values = np.asarray([item[key] for item in parameter_rows], dtype=float)
                for statistic, value in (
                    ("sum", values.sum()),
                    ("mean", values.mean()),
                    ("std", values.std()),
                    ("min", values.min()),
                    ("max", values.max()),
                ):
                    out[f"arrow_{key.lower()}_{statistic}"] = float(value)
        preliminary.append(out)
    rows: list[dict[str, Any]] = []
    for out in preliminary:
        counts = out.pop("_type_counts", None)
        for atom_type in sorted(all_types):
            clean = re.sub(r"[^A-Za-z0-9]+", "_", atom_type).strip("_").lower()
            out[f"arrow_type_count__{clean}"] = (
                float(counts.get(atom_type, 0)) if counts else np.nan
            )
        rows.append(out)
    return pd.DataFrame(rows)


def assign_folds(benchmark: pd.DataFrame, n_splits: int = 5) -> pd.DataFrame:
    result = benchmark.copy()
    molecules = result.drop_duplicates("molecule_id").reset_index(drop=True)
    random_folds = np.full(len(molecules), -1, dtype=int)
    for fold, (_, test) in enumerate(
        KFold(n_splits=n_splits, shuffle=True, random_state=20260826).split(molecules)
    ):
        random_folds[test] = fold
    random_map = dict(zip(molecules.molecule_id, random_folds))
    result["fold_random"] = result.molecule_id.map(random_map).astype(int)
    for column, output in (
        ("functional_group_family", "fold_family"),
        ("scaffold", "fold_scaffold"),
    ):
        assignments = np.full(len(molecules), -1, dtype=int)
        groups = molecules[column].astype(str).to_numpy()
        for fold, (_, test) in enumerate(
            GroupKFold(n_splits=n_splits).split(molecules, groups=groups)
        ):
            assignments[test] = fold
        fold_map = dict(zip(molecules.molecule_id, assignments))
        result[output] = result.molecule_id.map(fold_map).astype(int)
    if (result[["fold_random", "fold_family", "fold_scaffold"]] < 0).any().any():
        raise AssertionError("Not all rows received fold assignments")
    return result


def build_all() -> dict[str, pd.DataFrame]:
    benchmark = load_benchmark()
    freesolv = load_freesolv()
    combisolv = load_combisolv_water()
    benchmark = attach_freesolv_metadata(benchmark, freesolv)
    benchmark = assign_folds(benchmark)
    hydration = benchmark[benchmark.solvent == "water"].reset_index(drop=True)
    static = arrow_static_features(hydration)
    features = rdkit_descriptor_frame(hydration)
    public_hydration = combine_public_hydration(freesolv, combisolv, benchmark)
    public_features = rdkit_descriptor_frame(public_hydration)
    benchmark_connectivity = set(benchmark.inchi_connectivity_key)
    freesolv["benchmark_connectivity_match"] = freesolv.inchi_connectivity_key.isin(
        benchmark_connectivity
    )
    public_nonbenchmark = freesolv[
        (freesolv.formal_charge == 0)
        & (freesolv.fragment_count == 1)
        & (~freesolv.benchmark_connectivity_match)
    ].copy()
    return {
        "benchmark": benchmark,
        "freesolv": freesolv,
        "freesolv_nonbenchmark": public_nonbenchmark,
        "combisolv_water": combisolv,
        "public_hydration": public_hydration,
        "static": static,
        "features": features,
        "public_features": public_features,
    }
