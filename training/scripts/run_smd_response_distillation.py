"""Distill a benchmark-excluded SMD(water) teacher into hydration prediction.

The SMD values used here are predictions from a graph encoder pretrained after
removing every ARROW benchmark connectivity.  Consequently all evaluated inputs
remain deterministic functions of molecular structure.  The script compares a
direct response model with a physically anchored residual form and selects the
residual shrinkage only on inner outer-training predictions.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from arrow_distill.data import ROOT
from arrow_distill.experiments import (
    REGIMES,
    extra_trees,
    fit_with_optional_weights,
    make_prediction_rows,
)


def align(table: pd.DataFrame, ids: pd.Series, columns: list[str]) -> np.ndarray:
    indexed = table.drop_duplicates("molecule_id").set_index("molecule_id")
    values = indexed.reindex(ids.astype(str))[columns].to_numpy(dtype=np.float32)
    if not np.isfinite(values).all():
        raise AssertionError(f"Non-finite aligned values in {columns}")
    return values


def align_imputed(table: pd.DataFrame, ids: pd.Series, columns: list[str]) -> np.ndarray:
    indexed = table.drop_duplicates("molecule_id").set_index("molecule_id")
    numeric = indexed[columns].apply(pd.to_numeric, errors="coerce")
    medians = numeric.median()
    numeric = numeric.reindex(ids.astype(str)).fillna(medians).fillna(0.0)
    values = numeric.to_numpy(dtype=np.float32)
    if not np.isfinite(values).all():
        raise AssertionError(f"Non-finite aligned values in {columns}")
    return values


def train_predict(
    public_x: np.ndarray,
    public_y: np.ndarray,
    benchmark_x: np.ndarray,
    benchmark_y: np.ndarray,
    train: np.ndarray,
    test: np.ndarray,
    seed: int,
    n_estimators: int,
) -> np.ndarray:
    model = extra_trees(
        seed=seed,
        n_estimators=n_estimators,
        max_features=0.7,
        min_samples_leaf=2,
    )
    weights = np.concatenate([np.ones(len(public_y)), np.full(len(train), 3.0)])
    fit_with_optional_weights(
        model,
        np.vstack([public_x, benchmark_x[train]]),
        np.concatenate([public_y, benchmark_y[train]]),
        weights,
    )
    return model.predict(benchmark_x[test])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regimes", nargs="+", choices=tuple(REGIMES), default=["random_oof"])
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--seeds", type=int, nargs="+", default=[11, 29, 47])
    parser.add_argument(
        "--no-static-arrow",
        action="store_true",
        help="Use only SMILES-derived RDKit/Morgan features and distilled teachers.",
    )
    parser.add_argument("--output-stem", default="smd_response_distillation_oof")
    parser.add_argument(
        "--direct-only",
        action="store_true",
        help="Skip residual/shrinkage variants for fast fixed-model confirmation.",
    )
    parser.add_argument(
        "--confsolv",
        choices=("none", "scalar", "response", "hierarchy"),
        default="none",
        help="Optional structure-predicted ConfSolv water-response block.",
    )
    parser.add_argument(
        "--smd-feature",
        choices=("raw", "calibrated", "both"),
        default="raw",
        help="MolSolv structure-only teacher representation.",
    )
    parser.add_argument(
        "--gnnis",
        choices=("none", "energy", "force", "summary"),
        default="none",
        help="Optional static explicit-water neural-force teacher block.",
    )
    parser.add_argument(
        "--lambda-potential",
        action="store_true",
        help="Add the static lambda-response neural-potential teacher block.",
    )
    parser.add_argument(
        "--smd-embedding",
        choices=("none", "graph", "ffn", "both"),
        default="none",
        help="Add frozen latent coordinates from the benchmark-excluded SMD encoder.",
    )
    parser.add_argument(
        "--confsolv-embedding",
        choices=("none", "graph", "ffn", "both"),
        default="none",
        help="Add a frozen latent trained on the ConfSolv water-response hierarchy.",
    )
    parser.add_argument(
        "--openfe-diagnostics",
        action="store_true",
        help="Add structure-distilled OpenFE MBAR/convergence diagnostics.",
    )
    parser.add_argument(
        "--mlff-hfe",
        action="store_true",
        help="Add the narrow structure-distilled MLFF/force-field hydration hierarchy.",
    )
    parser.add_argument(
        "--des-water",
        action="store_true",
        help="Add the narrow DES370K water-dimer interaction/SAPT response block.",
    )
    args = parser.parse_args()

    processed = ROOT / "data/processed"
    benchmark = pd.read_parquet(processed / "arrow_solvation_master.parquet")
    benchmark = benchmark[benchmark.solvent.eq("water")].reset_index(drop=True)
    public = pd.read_parquet(processed / "expanded_public_hydration_nonbenchmark.parquet")
    base = pd.read_parquet(processed / "rdkit_morgan_features.parquet")
    public_base = pd.read_parquet(
        processed / "expanded_public_hydration_rdkit_morgan_features.parquet"
    )
    static = pd.read_parquet(processed / "arrow_static_features.parquet")
    qm = pd.read_parquet(processed / "combisolv_qm_teacher_predictions.parquet")
    abraham = pd.read_parquet(processed / "soluteml_abraham_teacher_predictions.parquet")
    openff = pd.read_parquet(processed / "openff_alchemical_teacher_predictions.parquet")
    implicit = pd.read_parquet(processed / "implicit_solvent_teacher_predictions.parquet")
    smd = pd.read_parquet(processed / "molsolv_smd_teacher_predictions.parquet")
    calibrated_smd_path = processed / "molsolv_smd_calibrated_teacher_predictions.parquet"
    calibrated_smd = pd.read_parquet(calibrated_smd_path) if calibrated_smd_path.is_file() else None
    confsolv_path = processed / "confsolv_water_teacher_predictions.parquet"
    confsolv = pd.read_parquet(confsolv_path) if confsolv_path.is_file() else None
    gnnis = pd.read_parquet(processed / "gnnis_static_response_features.parquet")
    lambda_potential = pd.read_parquet(processed / "lambda_potential_static_features.parquet")
    smd_embedding_path = processed / "molsolv_smd_embeddings.parquet"
    smd_embeddings = (
        pd.read_parquet(smd_embedding_path)
        if args.smd_embedding != "none" and smd_embedding_path.is_file()
        else None
    )
    confsolv_embedding_path = processed / "confsolv_physics_embeddings.parquet"
    confsolv_embeddings = (
        pd.read_parquet(confsolv_embedding_path)
        if args.confsolv_embedding != "none" and confsolv_embedding_path.is_file()
        else None
    )
    openfe_diagnostics_path = processed / "openfe_diagnostics_teacher_predictions.parquet"
    openfe_diagnostics = (
        pd.read_parquet(openfe_diagnostics_path)
        if args.openfe_diagnostics and openfe_diagnostics_path.is_file()
        else None
    )
    mlff_hfe_path = processed / "mlff_hfe_teacher_predictions.parquet"
    mlff_hfe = pd.read_parquet(mlff_hfe_path) if args.mlff_hfe and mlff_hfe_path.is_file() else None
    des_water_path = processed / "des370k_water_teacher_predictions.parquet"
    des_water = (
        pd.read_parquet(des_water_path) if args.des_water and des_water_path.is_file() else None
    )

    if set(benchmark.inchi_connectivity_key.astype(str)) & set(
        public.inchi_connectivity_key.astype(str)
    ):
        raise AssertionError("Benchmark identity leaked into public hydration data")
    descriptors = [column for column in base if column.startswith(("rdkit__", "morgan2__"))]
    arrow_columns = [
        column
        for column in static
        if column.startswith("arrow_")
        and column not in {"arrow_hin_path", "arrow_static_available"}
        and pd.api.types.is_numeric_dtype(static[column])
    ]
    if args.no_static_arrow:
        x_structure = base[descriptors].to_numpy(dtype=np.float32)
        p_structure = public_base[descriptors].to_numpy(dtype=np.float32)
    else:
        x_structure = np.column_stack(
            [
                base[descriptors].to_numpy(dtype=np.float32),
                static.arrow_static_available.astype(float).to_numpy(),
                static[arrow_columns].to_numpy(dtype=np.float32),
            ]
        )
        p_structure = np.column_stack(
            [
                public_base[descriptors].to_numpy(dtype=np.float32),
                np.zeros((len(public), 1 + len(arrow_columns)), dtype=np.float32),
            ]
        )
    qcols = ["combisolv_qm_teacher"]
    acols = [f"abraham_{name}_teacher" for name in "esabl"]
    ocols = ["openff23_dg_teacher", "openff23_exp_residual_teacher"]
    icols = ["gbn2_alchemical_dg_teacher", "gbn2_exp_residual_teacher"]
    x_openff = align(openff, benchmark.molecule_id, ocols)
    p_openff = align(openff, public.molecule_id, ocols)
    x_implicit = align(implicit, benchmark.molecule_id, icols)
    p_implicit = align(implicit, public.molecule_id, icols)
    raw_x_smd = align(smd, benchmark.molecule_id, ["molsolv_smd_teacher"])
    raw_p_smd = align(smd, public.molecule_id, ["molsolv_smd_teacher"])
    if args.smd_feature == "raw":
        x_smd_block = raw_x_smd
        p_smd_block = raw_p_smd
    else:
        if calibrated_smd is None:
            raise FileNotFoundError(calibrated_smd_path)
        calibrated_columns = ["molsolv_smd_calibrated_teacher"]
        calibrated_x = align(calibrated_smd, benchmark.molecule_id, calibrated_columns)
        calibrated_p = align(calibrated_smd, public.molecule_id, calibrated_columns)
        if args.smd_feature == "calibrated":
            x_smd_block = calibrated_x
            p_smd_block = calibrated_p
        else:
            x_smd_block = np.column_stack([raw_x_smd, calibrated_x])
            p_smd_block = np.column_stack([raw_p_smd, calibrated_p])
    # The residual formulation needs a single physical anchor.  Use the
    # calibrated SMD value when requested and the original graph prediction otherwise.
    x_smd = x_smd_block[:, -1]
    p_smd = p_smd_block[:, -1]
    confsolv_columns = [
        "confsolv_water_dg_ensemble_teacher",
        "confsolv_water_gsolv_gas_lec_teacher",
        "confsolv_water_gsolv_solution_lec_teacher",
        "confsolv_water_gsolv_gas_weighted_teacher",
        "confsolv_water_gsolv_solution_weighted_teacher",
        "confsolv_gas_conformer_correction_teacher",
        "confsolv_solution_conformer_correction_teacher",
        "confsolv_hydration_conformer_correction_teacher",
        "confsolv_water_gsolv_std_teacher",
        "confsolv_water_response_mean_teacher",
        "confsolv_water_response_std_teacher",
    ]
    if args.confsolv != "none":
        if confsolv is None:
            raise FileNotFoundError(confsolv_path)
        if args.confsolv == "scalar":
            selected_confsolv = confsolv_columns[:1]
        elif args.confsolv == "response":
            selected_confsolv = confsolv_columns[5:]
        else:
            selected_confsolv = confsolv_columns
        x_confsolv = align(confsolv, benchmark.molecule_id, selected_confsolv)
        p_confsolv = align(confsolv, public.molecule_id, selected_confsolv)
    else:
        x_confsolv = np.empty((len(benchmark), 0), dtype=np.float32)
        p_confsolv = np.empty((len(public), 0), dtype=np.float32)
    if args.gnnis == "none":
        x_gnnis = np.empty((len(benchmark), 0), dtype=np.float32)
        p_gnnis = np.empty((len(public), 0), dtype=np.float32)
    else:
        prefixes = {
            "energy": ("gnnis_energy_",),
            "force": ("gnnis_force_", "gnnis_response_"),
            "summary": (
                "gnnis_conformer_",
                "gnnis_energy_",
                "gnnis_force_",
                "gnnis_response_",
            ),
        }[args.gnnis]
        gnnis_columns = [
            column
            for column in gnnis
            if column.startswith(prefixes) and pd.api.types.is_numeric_dtype(gnnis[column])
        ]
        x_gnnis = align_imputed(gnnis, benchmark.molecule_id, gnnis_columns)
        p_gnnis = align_imputed(gnnis, public.molecule_id, gnnis_columns)
    if args.lambda_potential:
        lambda_columns = [
            column
            for column in lambda_potential
            if column.startswith("lambda_potential_")
            and column not in {"lambda_potential_available"}
            and pd.api.types.is_numeric_dtype(lambda_potential[column])
        ]
        x_lambda = align_imputed(lambda_potential, benchmark.molecule_id, lambda_columns)
        p_lambda = align_imputed(lambda_potential, public.molecule_id, lambda_columns)
    else:
        x_lambda = np.empty((len(benchmark), 0), dtype=np.float32)
        p_lambda = np.empty((len(public), 0), dtype=np.float32)
    if args.smd_embedding == "none":
        x_embedding = np.empty((len(benchmark), 0), dtype=np.float32)
        p_embedding = np.empty((len(public), 0), dtype=np.float32)
    else:
        if smd_embeddings is None:
            raise FileNotFoundError(smd_embedding_path)
        embedding_prefixes = {
            "graph": ("molsolv_smd_graph_pca_",),
            "ffn": ("molsolv_smd_ffn_",),
            "both": ("molsolv_smd_graph_pca_", "molsolv_smd_ffn_"),
        }[args.smd_embedding]
        embedding_columns = [
            column for column in smd_embeddings if column.startswith(embedding_prefixes)
        ]
        x_embedding = align(smd_embeddings, benchmark.molecule_id, embedding_columns)
        p_embedding = align(smd_embeddings, public.molecule_id, embedding_columns)
    if args.confsolv_embedding == "none":
        x_confsolv_embedding = np.empty((len(benchmark), 0), dtype=np.float32)
        p_confsolv_embedding = np.empty((len(public), 0), dtype=np.float32)
    else:
        if confsolv_embeddings is None:
            raise FileNotFoundError(confsolv_embedding_path)
        confsolv_prefixes = {
            "graph": ("confsolv_graph_pca_",),
            "ffn": ("confsolv_ffn_",),
            "both": ("confsolv_graph_pca_", "confsolv_ffn_"),
        }[args.confsolv_embedding]
        confsolv_embedding_columns = [
            column for column in confsolv_embeddings if column.startswith(confsolv_prefixes)
        ]
        x_confsolv_embedding = align(
            confsolv_embeddings, benchmark.molecule_id, confsolv_embedding_columns
        )
        p_confsolv_embedding = align(
            confsolv_embeddings, public.molecule_id, confsolv_embedding_columns
        )
    if args.openfe_diagnostics:
        if openfe_diagnostics is None:
            raise FileNotFoundError(openfe_diagnostics_path)
        # Deliberately omit final dG and experimental-residual targets here: the
        # existing OpenFF block already supplies that information.  This block
        # asks whether structure-distilled sampling/response diagnostics add
        # genuinely new information to the strong SMD+ConfSolv student.
        diagnostic_targets = [
            "openfe_archive_repeat_sd",
            "openfe_solvent_estimate_error_mean",
            "openfe_solvent_overlap_scalar_mean",
            "openfe_solvent_neighbor_overlap_mean_mean",
            "openfe_solvent_exchange_spectral_gap_mean",
            "openfe_solvent_forward_reverse_disagreement_mean",
            "openfe_solvent_early_to_final_drift_mean",
            "openfe_solvent_convergence_range_mean",
            "openfe_vacuum_estimate_error_mean",
            "openfe_vacuum_overlap_scalar_mean",
            "openfe_vacuum_exchange_spectral_gap_mean",
            "openfe_vacuum_forward_reverse_disagreement_mean",
            "nagl_charge_abs_sum",
            "nagl_charge_sq_sum",
            "nagl_charge_std",
            "nagl_heavy_charge_abs_max",
        ]
        diagnostic_columns = [f"{target}_teacher" for target in diagnostic_targets]
        x_openfe_diagnostics = align(openfe_diagnostics, benchmark.molecule_id, diagnostic_columns)
        p_openfe_diagnostics = align(openfe_diagnostics, public.molecule_id, diagnostic_columns)
    else:
        x_openfe_diagnostics = np.empty((len(benchmark), 0), dtype=np.float32)
        p_openfe_diagnostics = np.empty((len(public), 0), dtype=np.float32)
    if args.mlff_hfe:
        if mlff_hfe is None:
            raise FileNotFoundError(mlff_hfe_path)
        mlff_columns = [
            "mlff_hfe_organic_mpnice_teacher",
            "mlff_hfe_gaff_teacher",
            "mlff_hfe_opls4_teacher",
            "mlff_hfe_esol_teacher",
            "mlff_hfe_dft_pbf_teacher",
            "mlff_hfe_mpnice_minus_gaff_teacher",
            "mlff_hfe_mpnice_minus_dft_pbf_teacher",
        ]
        x_mlff_hfe = align(mlff_hfe, benchmark.molecule_id, mlff_columns)
        p_mlff_hfe = align(mlff_hfe, public.molecule_id, mlff_columns)
    else:
        x_mlff_hfe = np.empty((len(benchmark), 0), dtype=np.float32)
        p_mlff_hfe = np.empty((len(public), 0), dtype=np.float32)
    if args.des_water:
        if des_water is None:
            raise FileNotFoundError(des_water_path)
        des_columns = [
            "des_md_solvation_cbs_mean_teacher",
            "des_md_solvation_cbs_std_teacher",
            "des_md_solvation_cbs_min_teacher",
            "des_md_solvation_cbs_q10_teacher",
            "des_md_solvation_cbs_q50_teacher",
            "des_md_solvation_cbs_q90_teacher",
            "des_md_solvation_sapt_es_mean_teacher",
            "des_md_solvation_sapt_ex_mean_teacher",
            "des_md_solvation_sapt_ind_mean_teacher",
            "des_md_solvation_sapt_disp_mean_teacher",
            "des_md_solvation_sapt_all_mean_teacher",
        ]
        x_des_water = align(des_water, benchmark.molecule_id, des_columns)
        p_des_water = align(des_water, public.molecule_id, des_columns)
    else:
        x_des_water = np.empty((len(benchmark), 0), dtype=np.float32)
        p_des_water = np.empty((len(public), 0), dtype=np.float32)
    x = np.column_stack(
        [
            x_structure,
            align(qm, benchmark.molecule_id, qcols),
            align(abraham, benchmark.molecule_id, acols),
            x_openff[:, 0] + x_openff[:, 1],
            x_implicit[:, 0] + x_implicit[:, 1],
            x_smd_block,
            x_confsolv,
            x_gnnis,
            x_lambda,
            x_embedding,
            x_confsolv_embedding,
            x_openfe_diagnostics,
            x_mlff_hfe,
            x_des_water,
        ]
    )
    p = np.column_stack(
        [
            p_structure,
            align(qm, public.molecule_id, qcols),
            align(abraham, public.molecule_id, acols),
            p_openff[:, 0] + p_openff[:, 1],
            p_implicit[:, 0] + p_implicit[:, 1],
            p_smd_block,
            p_confsolv,
            p_gnnis,
            p_lambda,
            p_embedding,
            p_confsolv_embedding,
            p_openfe_diagnostics,
            p_mlff_hfe,
            p_des_water,
        ]
    )
    old_public = pd.read_parquet(processed / "public_hydration_nonbenchmark.parquet")
    old_keys = set(old_public.inchi_connectivity_key.astype(str))
    source_mask = (
        public.inchi_connectivity_key.astype(str).isin(old_keys)
        | public.source_measurement_count.fillna(0).ge(2)
    ).to_numpy()
    p = p[source_mask]
    public_truth = public.delta_g_exp.to_numpy(dtype=float)[source_mask]
    p_smd = p_smd[source_mask]
    truth = benchmark.delta_g_exp.to_numpy(dtype=float)
    alphas = np.asarray([0.25, 0.5, 0.75, 1.0])
    rows: list[pd.DataFrame] = []
    selections: list[dict[str, object]] = []

    for regime in args.regimes:
        folds = benchmark[REGIMES[regime]].to_numpy()
        direct = np.full(len(benchmark), np.nan)
        residual_by_alpha = {alpha: np.full(len(benchmark), np.nan) for alpha in alphas}
        nested = np.full(len(benchmark), np.nan)
        for fold in sorted(np.unique(folds)):
            outer_train = np.flatnonzero(folds != fold)
            outer_test = np.flatnonzero(folds == fold)
            direct_seed_predictions = []
            residual_seed_predictions = []
            for seed in args.seeds:
                direct_seed_predictions.append(
                    train_predict(
                        p,
                        public_truth,
                        x,
                        truth,
                        outer_train,
                        outer_test,
                        seed,
                        args.n_estimators,
                    )
                )
                if not args.direct_only:
                    residual_seed_predictions.append(
                        train_predict(
                            p,
                            public_truth - p_smd,
                            x,
                            truth - x_smd,
                            outer_train,
                            outer_test,
                            seed + 1000,
                            args.n_estimators,
                        )
                    )
            direct[outer_test] = np.mean(direct_seed_predictions, axis=0)
            if not args.direct_only:
                residual = np.mean(residual_seed_predictions, axis=0)
                for alpha in alphas:
                    residual_by_alpha[alpha][outer_test] = x_smd[outer_test] + alpha * residual

                inner_prediction = {alpha: np.full(len(outer_train), np.nan) for alpha in alphas}
                for inner_fold in sorted(np.unique(folds[outer_train])):
                    inner_test = outer_train[folds[outer_train] == inner_fold]
                    inner_train = outer_train[folds[outer_train] != inner_fold]
                    values = []
                    for seed in args.seeds:
                        values.append(
                            train_predict(
                                p,
                                public_truth - p_smd,
                                x,
                                truth - x_smd,
                                inner_train,
                                inner_test,
                                seed + 2000,
                                args.n_estimators,
                            )
                        )
                    inner_residual = np.mean(values, axis=0)
                    locations = np.searchsorted(outer_train, inner_test)
                    for alpha in alphas:
                        inner_prediction[alpha][locations] = (
                            x_smd[inner_test] + alpha * inner_residual
                        )
                inner_mae = {
                    alpha: float(np.mean(np.abs(truth[outer_train] - values)))
                    for alpha, values in inner_prediction.items()
                }
                selected_alpha = min(inner_mae, key=inner_mae.get)
                nested[outer_test] = residual_by_alpha[selected_alpha][outer_test]
                selections.append(
                    {
                        "regime": regime,
                        "outer_fold": int(fold),
                        "selected_alpha": float(selected_alpha),
                        **{f"inner_mae_alpha_{alpha:g}": mae for alpha, mae in inner_mae.items()},
                    }
                )

        parts = []
        if args.smd_feature != "raw":
            parts.append(f"SMD {args.smd_feature}")
        if args.confsolv != "none":
            parts.append(f"ConfSolv {args.confsolv}")
        if args.gnnis != "none":
            parts.append(f"GNNIS {args.gnnis}")
        if args.lambda_potential:
            parts.append("lambda potential")
        if args.smd_embedding != "none":
            parts.append(f"SMD {args.smd_embedding} embedding")
        if args.confsolv_embedding != "none":
            parts.append(f"ConfSolv {args.confsolv_embedding} embedding")
        if args.openfe_diagnostics:
            parts.append("OpenFE diagnostics")
        if args.mlff_hfe:
            parts.append("MLFF HFE hierarchy")
        if args.des_water:
            parts.append("DES370K water response")
        suffix = "" if not parts else " + " + " + ".join(parts)
        candidates = {
            f"SMD response direct{suffix}": direct,
        }
        if not args.direct_only:
            candidates[f"SMD residual nested shrinkage{suffix}"] = nested
            candidates.update(
                {
                    f"SMD residual alpha={alpha:g}{suffix}": prediction
                    for alpha, prediction in residual_by_alpha.items()
                }
            )
        for method, prediction in candidates.items():
            rows.append(
                make_prediction_rows(
                    benchmark,
                    prediction,
                    method,
                    regime,
                    False,
                    details=["SMD teacher was pretrained after global benchmark-identity exclusion"]
                    * len(benchmark),
                )
            )
            print(regime, method, np.mean(np.abs(truth - prediction)), flush=True)

    output = pd.concat(rows, ignore_index=True)
    output_path = ROOT / f"results/{args.output_stem}.parquet"
    output.to_parquet(output_path, index=False)
    metrics = (
        output.groupby(["method", "regime"], as_index=False)
        .agg(n=("absolute_error", "size"), mae=("absolute_error", "mean"))
        .sort_values(["regime", "mae"])
    )
    metrics.to_csv(output_path.with_suffix(".csv"), index=False)
    pd.DataFrame(selections).to_csv(
        ROOT / f"results/{args.output_stem}_selections.csv", index=False
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
