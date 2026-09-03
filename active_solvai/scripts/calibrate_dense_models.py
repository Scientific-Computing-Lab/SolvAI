"""Calibrate frozen dense-curve GP hyperparameters on four designated molecules."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from active_solvai.dense import (
    condition_curve,
    curvature_scale,
    interpolate_three_point_prior,
    rbf_covariance,
)
from active_solvai.ledger import sha256

RELEASE_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"
CONFIG = ACTIVE_ROOT / "configs/dense_sentinel_v1.json"
RESPONSES = ACTIVE_ROOT / "results/phase2/dense_responses_calibration.parquet"
PHASE1_PRIORS = ACTIVE_ROOT / "results/phase1/phase1_response_predictions.parquet"
OUT = ACTIVE_ROOT / "release/DENSE_SENTINEL_CALIBRATION_LOCK.json"
OUT_MD = ACTIVE_ROOT / "release/DENSE_SENTINEL_CALIBRATION_LOCK.md"


def load_priors(names: list[str], grid: np.ndarray) -> dict[str, np.ndarray]:
    frame = pd.read_parquet(PHASE1_PRIORS)
    frame = frame.loc[
        frame.partition.eq("standardized_exclusion_primary")
        & frame.repeat.eq(-1)
        & np.isclose(frame.trajectory_fraction, 1.0)
        & frame.component.eq("total")
        & frame.molecule_name.isin(names)
    ]
    priors: dict[str, np.ndarray] = {}
    for name, group in frame.groupby("molecule_name"):
        group = group.sort_values("lambda")
        if len(group) != 3:
            raise AssertionError(f"Expected three prior coordinates for {name}")
        priors[name] = interpolate_three_point_prior(
            group.predicted_structure_only.to_numpy(float), grid
        )
    if set(priors) != set(names):
        raise AssertionError("Missing frozen structure priors")
    return priors


def nlpd(y: np.ndarray, mean: np.ndarray, variance: np.ndarray) -> float:
    variance = np.maximum(np.asarray(variance, dtype=float), 1e-8)
    residual = np.asarray(y, dtype=float) - np.asarray(mean, dtype=float)
    return float(np.mean(0.5 * (np.log(2.0 * np.pi * variance) + residual**2 / variance)))


def main() -> None:
    prospective_energy = list(
        (ACTIVE_ROOT / "simulations/dense_pimd2/prospective").glob("**/output/*__SYSTEM.ene")
    )
    if prospective_energy:
        raise AssertionError(
            "Prospective sentinel responses already exist; calibration lock must precede them"
        )
    config = json.loads(CONFIG.read_text())
    grid = np.asarray(config["lambda_grid"], dtype=float)
    initial = np.array([2, 6, 12], dtype=int)
    hidden = np.array([index for index in range(len(grid)) if index not in initial])
    frame = pd.read_parquet(RESPONSES)
    names = config["calibration_molecules"]
    counts = frame.groupby("molecule_name").size().reindex(names)
    if not (counts == len(grid)).all():
        raise AssertionError(f"Calibration curves incomplete: {counts.to_dict()}")
    curves = {
        name: group.sort_values("lambda").mean_dhdl_kcal_mol.to_numpy(float)
        for name, group in frame.groupby("molecule_name")
    }
    sems = {
        name: group.sort_values("lambda").five_block_sem_kcal_mol.to_numpy(float)
        for name, group in frame.groupby("molecule_name")
    }
    structure_priors = load_priors(names, grid)
    rows: list[dict[str, float | str]] = []
    for prior_kind in ("generic", "solvai"):
        for lengthscale in config["calibration"]["lengthscale_candidates"]:
            for noise_multiplier in config["calibration"]["noise_inflation_candidates"]:
                scores = []
                for target in names:
                    training = [name for name in names if name != target]
                    if prior_kind == "generic":
                        prior = np.mean([curves[name] for name in training], axis=0)
                        train_residuals = np.concatenate(
                            [curves[name] - prior for name in training]
                        )
                        local = None
                    else:
                        prior = structure_priors[target]
                        train_residuals = np.concatenate(
                            [curves[name] - structure_priors[name] for name in training]
                        )
                        local = curvature_scale(prior, grid)
                    amplitude = float(np.clip(np.sqrt(np.mean(train_residuals**2)), 1.0, 20.0))
                    covariance = rbf_covariance(grid, amplitude, lengthscale, local)
                    posterior = condition_curve(
                        prior,
                        covariance,
                        initial,
                        curves[target][initial],
                        noise_multiplier * sems[target][initial] ** 2,
                    )
                    scores.append(
                        nlpd(
                            curves[target][hidden],
                            posterior.mean[hidden],
                            np.diag(posterior.covariance)[hidden],
                        )
                    )
                rows.append(
                    {
                        "prior": prior_kind,
                        "lengthscale": float(lengthscale),
                        "noise_inflation": float(noise_multiplier),
                        "mean_nlpd": float(np.mean(scores)),
                    }
                )
    search = pd.DataFrame(rows).sort_values(
        ["prior", "mean_nlpd", "lengthscale", "noise_inflation"]
    )
    selected: dict[str, dict[str, float]] = {}
    all_curves = np.stack([curves[name] for name in names])
    all_sems = np.stack([sems[name] for name in names])
    generic_mean = all_curves.mean(axis=0)
    for prior_kind in ("generic", "solvai"):
        best = search.loc[search.prior.eq(prior_kind)].iloc[0]
        if prior_kind == "generic":
            residuals = all_curves - generic_mean
        else:
            residuals = np.stack([curves[name] - structure_priors[name] for name in names])
        selected[prior_kind] = {
            "lengthscale": float(best.lengthscale),
            "noise_inflation": float(best.noise_inflation),
            "amplitude": float(np.clip(np.sqrt(np.mean(residuals**2)), 1.0, 20.0)),
            "calibration_mean_nlpd": float(best.mean_nlpd),
        }
    search_path = ACTIVE_ROOT / "results/phase2/dense_calibration_search.csv"
    search.to_csv(search_path, index=False)
    payload = {
        "schema_version": 1,
        "status": "calibration_locked_before_prospective_generation",
        "calibration_molecules": names,
        "lambda_grid": grid.tolist(),
        "initial_indices": initial.tolist(),
        "initial_lambda_values": grid[initial].tolist(),
        "selection": selected,
        "generic_prior_mean": generic_mean.tolist(),
        "expected_sem_by_lambda": np.median(all_sems, axis=0).tolist(),
        "inputs": {
            "configuration": {"path": str(CONFIG), "sha256": sha256(CONFIG)},
            "calibration_responses": {"path": str(RESPONSES), "sha256": sha256(RESPONSES)},
            "phase1_priors": {"path": str(PHASE1_PRIORS), "sha256": sha256(PHASE1_PRIORS)},
            "search_table": {"path": str(search_path), "sha256": sha256(search_path)},
            "calibration_script": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256(Path(__file__).resolve()),
            },
        },
        "note": "No prospective sentinel dense response was generated or read before this lock.",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    OUT_MD.write_text(
        "# Dense sentinel calibration lock\n\n"
        "This file records the automatically selected Gaussian-process settings before any "
        "prospective-sentinel dense response was generated. Selection followed the committed "
        "protocol in `DENSE_SENTINEL_FREEZE.md`; no setting was selected using endpoint error.\n\n"
        f"- Calibration response SHA-256: `{payload['inputs']['calibration_responses']['sha256']}`\n"
        f"- Phase 1 prior SHA-256: `{payload['inputs']['phase1_priors']['sha256']}`\n"
        f"- Generic prior: length scale `{selected['generic']['lengthscale']}`, noise inflation "
        f"`{selected['generic']['noise_inflation']}`, amplitude `{selected['generic']['amplitude']:.8f}` kcal mol^-1.\n"
        f"- SolvAI-conditioned prior: length scale `{selected['solvai']['lengthscale']}`, noise inflation "
        f"`{selected['solvai']['noise_inflation']}`, amplitude `{selected['solvai']['amplitude']:.8f}` kcal mol^-1.\n"
        "- Prospective dense responses available at lock time: `0`.\n"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
