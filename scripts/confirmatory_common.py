#!/usr/bin/env python3
"""Shared, frozen utilities for the SolvAI Phase 1 confirmation package."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline

MODEL_SEEDS = (11, 29, 47)
REPEAT_SEEDS = (314159, 271828, 161803, 141421, 173205)
SHUFFLE_SEEDS = (88001, 88002, 88003, 88004, 88005)
BOOTSTRAP_SEED = 20260828


@dataclass
class ConfirmatoryData:
    benchmark: pd.DataFrame
    public: pd.DataFrame
    benchmark_structure: np.ndarray
    public_structure: np.ndarray
    benchmark_responses: dict[str, np.ndarray]
    public_responses: dict[str, np.ndarray]
    feature_sets: dict[str, tuple[np.ndarray, np.ndarray]]
    response_names: list[str]


def _align(table: pd.DataFrame, ids: pd.Series, columns: list[str]) -> np.ndarray:
    indexed = table.drop_duplicates("molecule_id").set_index("molecule_id")
    values = indexed.reindex(ids.astype(str))[columns].to_numpy(dtype=np.float32)
    if values.shape != (len(ids), len(columns)) or not np.isfinite(values).all():
        raise AssertionError(f"Invalid alignment for {columns}")
    return values


def load_confirmatory_data(
    workspace_root: Path, teacher_overrides: dict[str, Path] | None = None
) -> ConfirmatoryData:
    processed = workspace_root / "data" / "processed"
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = (
        benchmark.loc[benchmark.solvent.eq("water")]
        .drop_duplicates("molecule_id")
        .reset_index(drop=True)
    )
    if len(benchmark) != 85:
        raise AssertionError(f"Expected 85 benchmark rows, found {len(benchmark)}")

    public_all = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    old_public = pd.read_parquet(processed / "public_hydration_nonbenchmark.parquet")
    old_keys = set(old_public.inchi_connectivity_key.astype(str))
    source_mask = (
        public_all.inchi_connectivity_key.astype(str).isin(old_keys)
        | public_all.source_measurement_count.fillna(0).ge(2)
    ).to_numpy()
    public = public_all.loc[source_mask].reset_index(drop=True)
    if len(public) != 1280:
        raise AssertionError(f"Expected 1,280 public endpoint rows, found {len(public)}")
    if set(benchmark.inchi_connectivity_key.astype(str)) & set(
        public.inchi_connectivity_key.astype(str)
    ):
        raise AssertionError("Exact benchmark connectivity in endpoint training pool")

    benchmark_features = pd.read_parquet(processed / "rdkit_morgan_features.parquet")
    public_features_all = pd.read_parquet(
        processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    )
    if not benchmark.molecule_id.equals(benchmark_features.molecule_id):
        raise AssertionError("Benchmark structure-feature order changed")
    if not public_all.molecule_id.equals(public_features_all.molecule_id):
        raise AssertionError("Public structure-feature order changed")
    public_features = public_features_all.loc[source_mask].reset_index(drop=True)
    descriptor_columns = [
        column for column in benchmark_features if column.startswith(("rdkit__", "morgan2__"))
    ]
    if len(descriptor_columns) != 2265:
        raise AssertionError(f"Expected 2,265 structure columns, found {len(descriptor_columns)}")
    benchmark_structure = benchmark_features[descriptor_columns].to_numpy(dtype=np.float32)
    public_structure = public_features[descriptor_columns].to_numpy(dtype=np.float32)

    teacher_overrides = teacher_overrides or {}

    def teacher_table(name: str, default_name: str) -> pd.DataFrame:
        return pd.read_parquet(teacher_overrides.get(name, processed / default_name))

    qm = teacher_table("combisolv_qm", "combisolv_qm_teacher_predictions.parquet")
    abraham = pd.read_parquet(processed / "soluteml_abraham_teacher_predictions.parquet")
    openff = pd.read_parquet(processed / "openff_alchemical_teacher_predictions.parquet")
    implicit = pd.read_parquet(processed / "implicit_solvent_teacher_predictions.parquet")
    smd = teacher_table("molsolv_smd", "molsolv_smd_teacher_predictions.parquet")
    confsolv = teacher_table("confsolv", "confsolv_water_teacher_predictions.parquet")

    qcols = ["combisolv_qm_teacher"]
    acols = [f"abraham_{name}_teacher" for name in "esabl"]
    ocols = ["openff23_dg_teacher", "openff23_exp_residual_teacher"]
    icols = ["gbn2_alchemical_dg_teacher", "gbn2_exp_residual_teacher"]
    scols = ["molsolv_smd_teacher"]
    ccols = [
        "confsolv_gas_conformer_correction_teacher",
        "confsolv_solution_conformer_correction_teacher",
        "confsolv_hydration_conformer_correction_teacher",
        "confsolv_water_gsolv_std_teacher",
        "confsolv_water_response_mean_teacher",
        "confsolv_water_response_std_teacher",
    ]

    def response_parts(ids: pd.Series) -> dict[str, np.ndarray]:
        qm_values = _align(qm, ids, qcols)
        abraham_values = _align(abraham, ids, acols)
        openff_values = _align(openff, ids, ocols)
        implicit_values = _align(implicit, ids, icols)
        smd_values = _align(smd, ids, scols)
        confsolv_values = _align(confsolv, ids, ccols)
        openff_corrected = openff_values.sum(axis=1, keepdims=True)
        gbn2_corrected = implicit_values.sum(axis=1, keepdims=True)
        empirical = np.column_stack([abraham_values, openff_corrected, gbn2_corrected]).astype(
            np.float32
        )
        computation_core = np.column_stack(
            [qm_values, openff_values[:, [0]], implicit_values[:, [0]]]
        ).astype(np.float32)
        narrow = np.column_stack(
            [qm_values, abraham_values, openff_corrected, gbn2_corrected]
        ).astype(np.float32)
        full = np.column_stack([narrow, smd_values, confsolv_values]).astype(np.float32)
        return {
            "empirical_residual": empirical,
            "computation_core": computation_core,
            "smd": smd_values,
            "confsolv": confsolv_values,
            "narrow": narrow,
            "narrow_smd": np.column_stack([narrow, smd_values]).astype(np.float32),
            "full": full,
        }

    benchmark_responses = response_parts(benchmark.molecule_id)
    public_responses_all = response_parts(public_all.molecule_id)
    public_responses = {name: values[source_mask] for name, values in public_responses_all.items()}

    def add(base: np.ndarray, response: np.ndarray | None) -> np.ndarray:
        return base if response is None else np.column_stack([base, response])

    feature_spec = {
        "A_structure_only": None,
        "B_empirical_residual": "empirical_residual",
        "C_computation_core": "computation_core",
        "D_smd_water": "smd",
        "E_confsolv": "confsolv",
        "F_full_solvai": "full",
        "G_narrow_reference": "narrow",
        "H_narrow_smd_reference": "narrow_smd",
    }
    feature_sets = {
        name: (
            add(benchmark_structure, benchmark_responses[key] if key else None),
            add(public_structure, public_responses[key] if key else None),
        )
        for name, key in feature_spec.items()
    }
    response_names = [
        "combisolv_qm",
        "abraham_e",
        "abraham_s",
        "abraham_a",
        "abraham_b",
        "abraham_l",
        "openff_corrected",
        "gbn2_corrected",
        "smd_water",
        "confsolv_gas_conformer_correction",
        "confsolv_solution_conformer_correction",
        "confsolv_hydration_conformer_correction",
        "confsolv_water_gsolv_std",
        "confsolv_water_response_mean",
        "confsolv_water_response_std",
    ]
    if benchmark_responses["full"].shape[1] != 15:
        raise AssertionError("Full response block does not contain 15 coordinates")
    return ConfirmatoryData(
        benchmark=benchmark,
        public=public,
        benchmark_structure=benchmark_structure,
        public_structure=public_structure,
        benchmark_responses=benchmark_responses,
        public_responses=public_responses,
        feature_sets=feature_sets,
        response_names=response_names,
    )


def endpoint_model(seed: int):
    return make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        ExtraTreesRegressor(
            n_estimators=360,
            max_features=0.7,
            min_samples_leaf=2,
            criterion="squared_error",
            bootstrap=False,
            max_depth=None,
            min_samples_split=2,
            random_state=seed,
            n_jobs=-1,
        ),
    )


def fit_predict(
    public_x: np.ndarray,
    public_y: np.ndarray,
    benchmark_x: np.ndarray,
    benchmark_y: np.ndarray,
    train_indices: np.ndarray,
    test_indices: np.ndarray,
) -> np.ndarray:
    x_fit = np.vstack([public_x, benchmark_x[train_indices]])
    y_fit = np.concatenate([public_y, benchmark_y[train_indices]])
    weights = np.concatenate([np.ones(len(public_y)), np.full(len(train_indices), 3.0)])
    predictions = []
    for seed in MODEL_SEEDS:
        model = endpoint_model(seed)
        model.fit(x_fit, y_fit, extratreesregressor__sample_weight=weights)
        predictions.append(model.predict(benchmark_x[test_indices]))
    return np.mean(predictions, axis=0)


def metric_record(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    residual = y_true - y_pred
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
        "median_absolute_error": float(np.median(np.abs(residual))),
        "r2": float(r2_score(y_true, y_pred)),
    }


def bootstrap_difference(
    truth: np.ndarray,
    candidate: np.ndarray,
    baseline: np.ndarray,
    *,
    seed: int = BOOTSTRAP_SEED,
    draws: int = 100_000,
) -> dict[str, float | str]:
    differences = np.abs(truth - candidate) - np.abs(truth - baseline)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(truth), size=(draws, len(truth)))
    sampled = differences[indices].mean(axis=1)
    low, high = np.quantile(sampled, [0.025, 0.975])
    mean = float(differences.mean())
    if high < 0:
        outcome = "positive"
    elif low > 0:
        outcome = "negative"
    else:
        outcome = "neutral"
    return {
        "mae_difference_vs_structure": mean,
        "ci_low_95": float(low),
        "ci_high_95": float(high),
        "probability_lower_error": float(np.mean(sampled < 0)),
        "fraction_molecules_improved": float(np.mean(differences < 0)),
        "material": bool(abs(mean) >= 0.010),
        "outcome": outcome,
    }
