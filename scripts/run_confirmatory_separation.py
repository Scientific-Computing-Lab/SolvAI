#!/usr/bin/env python3
"""Run the frozen all-label chemical-separation evaluations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from confirmatory_common import (
    bootstrap_difference,
    fit_predict,
    load_confirmatory_data,
    metric_record,
)
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.model_selection import GroupKFold

FAMILY_RULES = (
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("amide", "[NX3][CX3](=O)"),
    ("ester", "[CX3](=O)[OX2][#6]"),
    ("aldehyde", "[CX3H1](=O)"),
    ("ketone", "[#6][CX3](=O)[#6]"),
    ("alcohol_or_phenol", "[OX2H1][#6]"),
    ("ether", "[#6][OD2][#6]"),
    ("amine", "[NX3;!$(N-C=O);!$(N-S=O)]"),
    ("thiol", "[SX2H1]"),
    ("sulfide_or_disulfide", "[#6][SX2,SX1][#6,S]"),
)
COMPILED_FAMILY_RULES = tuple((name, Chem.MolFromSmarts(smarts)) for name, smarts in FAMILY_RULES)
FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def molecule(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"RDKit could not parse {smiles!r}")
    return mol


def primary_family(mol: Chem.Mol) -> str:
    for name, query in COMPILED_FAMILY_RULES:
        if mol.HasSubstructMatch(query):
            return name
    if any(atom.GetIsAromatic() and atom.GetAtomicNum() != 6 for atom in mol.GetAtoms()):
        return "aromatic_heterocycle"
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return "carbocyclic_aromatic"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[CX3]")):
        return "alkene_or_diene"
    if all(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return "saturated_hydrocarbon"
    return "other"


def scaffold_key(mol: Chem.Mol, family: str) -> str:
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
    return scaffold or f"ACYCLIC::{family}"


def butina_assignments(fingerprints: list, cutoff: float = 0.30) -> np.ndarray:
    distances: list[float] = []
    for index in range(1, len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[index], fingerprints[:index])
        distances.extend(1.0 - value for value in similarities)
    clusters = Butina.ClusterData(
        distances,
        len(fingerprints),
        cutoff,
        isDistData=True,
        reordering=True,
    )
    assignments = np.full(len(fingerprints), -1, dtype=int)
    for cluster_id, members in enumerate(clusters):
        assignments[list(members)] = cluster_id
    if np.any(assignments < 0):
        raise AssertionError("Incomplete Butina assignments")
    return assignments


def max_similarity_to_test(candidate_fps: list, test_fps: list) -> np.ndarray:
    maxima = np.zeros(len(candidate_fps), dtype=float)
    for test_fp in test_fps:
        maxima = np.maximum(
            maxima,
            np.asarray(DataStructs.BulkTanimotoSimilarity(test_fp, candidate_fps)),
        )
    return maxima


def build_assignments(data: object) -> tuple[pd.DataFrame, list, list]:
    benchmark_mols = [molecule(value) for value in data.benchmark.canonical_smiles]
    public_mols = [molecule(value) for value in data.public.canonical_smiles]
    benchmark_families = [primary_family(value) for value in benchmark_mols]
    public_families = [primary_family(value) for value in public_mols]
    benchmark_scaffolds = [
        scaffold_key(mol, family)
        for mol, family in zip(benchmark_mols, benchmark_families, strict=True)
    ]
    public_scaffolds = [
        scaffold_key(mol, family) for mol, family in zip(public_mols, public_families, strict=True)
    ]
    benchmark_fps = [FP_GENERATOR.GetFingerprint(value) for value in benchmark_mols]
    public_fps = [FP_GENERATOR.GetFingerprint(value) for value in public_mols]

    combined = pd.concat(
        [
            pd.DataFrame(
                {
                    "pool": "benchmark",
                    "pool_index": np.arange(len(data.benchmark)),
                    "molecule_id": data.benchmark.molecule_id.astype(str),
                    "connectivity_key": data.benchmark.inchi_connectivity_key.astype(str),
                    "canonical_smiles": data.benchmark.canonical_smiles.astype(str),
                    "primary_family": benchmark_families,
                    "scaffold_key": benchmark_scaffolds,
                }
            ),
            pd.DataFrame(
                {
                    "pool": "public",
                    "pool_index": np.arange(len(data.public)),
                    "molecule_id": data.public.molecule_id.astype(str),
                    "connectivity_key": data.public.inchi_connectivity_key.astype(str),
                    "canonical_smiles": data.public.canonical_smiles.astype(str),
                    "primary_family": public_families,
                    "scaffold_key": public_scaffolds,
                }
            ),
        ],
        ignore_index=True,
    )
    combined["_original_index"] = np.arange(len(combined))
    order = np.argsort(combined.connectivity_key.to_numpy(), kind="stable")
    sorted_fps = benchmark_fps + public_fps
    sorted_fps = [sorted_fps[index] for index in order]
    sorted_clusters = butina_assignments(sorted_fps)
    clusters = np.empty(len(combined), dtype=int)
    clusters[order] = sorted_clusters
    combined["butina_cluster_0_70"] = clusters
    combined = combined.drop(columns="_original_index")
    return combined, benchmark_fps, public_fps


def group_folds(groups: np.ndarray):
    splitter = GroupKFold(n_splits=5)
    dummy = np.zeros((len(groups), 1))
    for fold, (train, test) in enumerate(splitter.split(dummy, groups=groups)):
        yield fold, train, test


def random_folds(data: object):
    fold_ids = data.benchmark.fold_random.to_numpy(dtype=int)
    for fold in sorted(np.unique(fold_ids)):
        yield int(fold), np.flatnonzero(fold_ids != fold), np.flatnonzero(fold_ids == fold)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--release-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--standardized-exclusion", action="store_true")
    args = parser.parse_args()
    out = args.release_root / "results" / "confirmatory"
    out.mkdir(parents=True, exist_ok=True)
    suffix = "standardized_exclusion_" if args.standardized_exclusion else ""
    teacher_overrides = None
    if args.standardized_exclusion:
        root = out / "teacher_refits"
        teacher_overrides = {
            "combisolv_qm": root / "combisolv_qm/teacher_predictions.parquet",
            "molsolv_smd": root / "molsolv_smd/teacher_predictions.parquet",
            "confsolv": root / "confsolv/teacher_predictions.parquet",
        }
    data = load_confirmatory_data(args.workspace_root, teacher_overrides)
    assignments, benchmark_fps, public_fps = build_assignments(data)
    assignments.to_parquet(out / f"{suffix}chemical_separation_assignments.parquet", index=False)
    benchmark_assignment = assignments.loc[assignments.pool.eq("benchmark")].sort_values(
        "pool_index"
    )
    public_assignment = assignments.loc[assignments.pool.eq("public")].sort_values("pool_index")

    truth = data.benchmark.delta_g_exp.to_numpy(dtype=float)
    public_truth = data.public.delta_g_exp.to_numpy(dtype=float)
    prediction_rows: list[pd.DataFrame] = []
    fold_rows: list[dict] = []
    regimes: list[tuple[str, object]] = [
        ("global_family", benchmark_assignment.primary_family.to_numpy()),
        ("global_scaffold", benchmark_assignment.scaffold_key.to_numpy()),
        ("global_butina_0_70", benchmark_assignment.butina_cluster_0_70.to_numpy()),
    ]

    for regime, groups in regimes:
        for fold, nominal_train, test in group_folds(groups):
            test_groups = set(groups[test])
            benchmark_keep = np.asarray(
                [index for index in nominal_train if groups[index] not in test_groups],
                dtype=int,
            )
            public_group_column = {
                "global_family": "primary_family",
                "global_scaffold": "scaffold_key",
                "global_butina_0_70": "butina_cluster_0_70",
            }[regime]
            public_groups = public_assignment[public_group_column].to_numpy()
            public_keep = np.flatnonzero(~np.isin(public_groups, list(test_groups)))
            for method in ("A_structure_only", "F_full_solvai"):
                benchmark_x, public_x = data.feature_sets[method]
                pred = fit_predict(
                    public_x[public_keep],
                    public_truth[public_keep],
                    benchmark_x,
                    truth,
                    benchmark_keep,
                    test,
                )
                for row_index, value in zip(test, pred, strict=True):
                    prediction_rows.append(
                        pd.DataFrame(
                            {
                                "regime": [regime],
                                "threshold": [np.nan],
                                "fold": [fold],
                                "method": [method],
                                "molecule_id": [data.benchmark.loc[row_index, "molecule_id"]],
                                "molecule_name": [data.benchmark.loc[row_index, "molecule_name"]],
                                "canonical_smiles": [
                                    data.benchmark.loc[row_index, "canonical_smiles"]
                                ],
                                "y_true": [truth[row_index]],
                                "y_pred": [value],
                            }
                        )
                    )
            fold_rows.append(
                {
                    "regime": regime,
                    "threshold": np.nan,
                    "fold": fold,
                    "n_test": len(test),
                    "n_arrow_train": len(benchmark_keep),
                    "n_public_train": len(public_keep),
                    "n_arrow_removed": len(nominal_train) - len(benchmark_keep),
                    "n_public_removed": len(data.public) - len(public_keep),
                    "test_groups": "|".join(map(str, sorted(test_groups))),
                    "maximum_remaining_similarity": np.nan,
                }
            )

    for threshold in (0.50, 0.60, 0.70, 0.80):
        regime = f"global_nn_{threshold:.2f}"
        for fold, nominal_train, test in random_folds(data):
            test_fps = [benchmark_fps[index] for index in test]
            public_similarity = max_similarity_to_test(public_fps, test_fps)
            arrow_train_fps = [benchmark_fps[index] for index in nominal_train]
            arrow_similarity = max_similarity_to_test(arrow_train_fps, test_fps)
            public_keep = np.flatnonzero(public_similarity < threshold)
            benchmark_keep = nominal_train[arrow_similarity < threshold]
            remaining = np.concatenate(
                [public_similarity[public_keep], arrow_similarity[arrow_similarity < threshold]]
            )
            for method in ("A_structure_only", "F_full_solvai"):
                benchmark_x, public_x = data.feature_sets[method]
                pred = fit_predict(
                    public_x[public_keep],
                    public_truth[public_keep],
                    benchmark_x,
                    truth,
                    benchmark_keep,
                    test,
                )
                for row_index, value in zip(test, pred, strict=True):
                    prediction_rows.append(
                        pd.DataFrame(
                            {
                                "regime": [regime],
                                "threshold": [threshold],
                                "fold": [fold],
                                "method": [method],
                                "molecule_id": [data.benchmark.loc[row_index, "molecule_id"]],
                                "molecule_name": [data.benchmark.loc[row_index, "molecule_name"]],
                                "canonical_smiles": [
                                    data.benchmark.loc[row_index, "canonical_smiles"]
                                ],
                                "y_true": [truth[row_index]],
                                "y_pred": [value],
                            }
                        )
                    )
            fold_rows.append(
                {
                    "regime": regime,
                    "threshold": threshold,
                    "fold": fold,
                    "n_test": len(test),
                    "n_arrow_train": len(benchmark_keep),
                    "n_public_train": len(public_keep),
                    "n_arrow_removed": len(nominal_train) - len(benchmark_keep),
                    "n_public_removed": len(data.public) - len(public_keep),
                    "test_groups": "",
                    "maximum_remaining_similarity": float(remaining.max()),
                }
            )

    predictions = pd.concat(prediction_rows, ignore_index=True)
    predictions["residual"] = predictions.y_true - predictions.y_pred
    predictions["absolute_error"] = predictions.residual.abs()
    predictions.to_parquet(out / f"{suffix}global_separation_predictions.parquet", index=False)
    fold_table = pd.DataFrame(fold_rows)
    fold_table.to_csv(out / f"{suffix}global_separation_training_counts.csv", index=False)

    metric_rows: list[dict] = []
    comparison_rows: list[dict] = []
    for regime, group in predictions.groupby("regime", sort=False):
        methods = {}
        for method, method_group in group.groupby("method"):
            if len(method_group) != 85 or method_group.molecule_id.nunique() != 85:
                raise AssertionError(f"Incomplete predictions for {regime}/{method}")
            method_group = method_group.set_index("molecule_id").loc[data.benchmark.molecule_id]
            methods[method] = method_group
            metric_rows.append(
                {
                    "regime": regime,
                    "threshold": method_group.threshold.iloc[0],
                    "method": method,
                    "n": len(method_group),
                    **metric_record(method_group.y_true, method_group.y_pred),
                }
            )
        baseline = methods["A_structure_only"]
        candidate = methods["F_full_solvai"]
        comparison_rows.append(
            {
                "regime": regime,
                "threshold": candidate.threshold.iloc[0],
                **bootstrap_difference(
                    candidate.y_true.to_numpy(),
                    candidate.y_pred.to_numpy(),
                    baseline.y_pred.to_numpy(),
                ),
            }
        )
        print(
            regime,
            "A",
            metric_record(baseline.y_true, baseline.y_pred),
            "F",
            metric_record(candidate.y_true, candidate.y_pred),
            flush=True,
        )
    pd.DataFrame(metric_rows).to_csv(out / f"{suffix}global_separation_metrics.csv", index=False)
    pd.DataFrame(comparison_rows).to_csv(
        out / f"{suffix}global_separation_paired_comparisons.csv", index=False
    )
    metadata = {
        "family_rules": [{"name": name, "smarts": smarts} for name, smarts in FAMILY_RULES],
        "family_precedence": [name for name, _ in FAMILY_RULES]
        + [
            "aromatic_heterocycle",
            "carbocyclic_aromatic",
            "alkene_or_diene",
            "saturated_hydrocarbon",
            "other",
        ],
        "butina_similarity": 0.70,
        "nearest_neighbor_thresholds": [0.50, 0.60, 0.70, 0.80],
        "benchmark_rows": len(data.benchmark),
        "public_rows": len(data.public),
        "standardized_exclusion_teachers": args.standardized_exclusion,
    }
    (out / f"{suffix}global_separation_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
