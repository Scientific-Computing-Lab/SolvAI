#!/usr/bin/env python3
"""Post-review comparisons; test labels never select fits or hyperparameters.

Run with the original data workspace and the release root supplied explicitly.
Stages are restartable. Protocol and input hashes precede model fitting.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

SEEDS = (11, 29, 47)
PROTOCOL = {
    "version": "2026-10-02-postreview-v1",
    "status": "retrospective extension; not a new preregistration or untouched test",
    "test_sets": ["original ARROW fixed five folds", "original external 220; strict subset 97"],
    "endpoint_labels": "same 1280 public labels; eligible ARROW labels weighted 3, public 1",
    "selection": "no test-dependent tuning or result-dependent removal of models",
    "seeds": SEEDS,
    "tree": {"n_estimators": 360, "min_samples_leaf": 2, "max_features": 0.7},
    "representations": ["structure", "SolvAI", "MoLFormer + structure", "MoLFormer + SolvAI", "Uni-Mol + structure", "Uni-Mol + SolvAI"],
    "graph": {"models": ["D-MPNN", "CheMeleon"], "epoch_cap": 80,
              "inner_validation_fraction_per_source": 0.15, "selection_metric": "unweighted validation MAE",
              "early_stopping_patience": 12, "batch_size": 64,
              "optimizer": "Adam constant learning rate", "scratch_lr": 0.0003,
              "pretrained_lr": 0.0001, "scratch_message_dim": 300,
              "scratch_depth": 3, "ffn_hidden": 300, "ffn_layers": 2, "dropout": 0.1,
              "refit": "fresh initialization on ALL outer-training labels for selected epoch count"},
    "limitations": ["foundation-model pretraining membership is not an audited disjoint source",
                    "external cohort labels have been used in earlier reporting, not these fits",
                    "no universal state-of-the-art claim; no unfiltered FreeSolv comparison"],
    "strata": "external nearest-training Morgan Tanimoto: [0,.3), [.3,.5), [.5,.7), [.7,1]",
    "statistics": "paired molecule bootstrap, 100000 resamples, seed 20260828; descriptive, no multiplicity adjustment",
}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + "\n")


def prepare(args):
    from confirmatory_common import load_confirmatory_data
    sys.path.insert(0, str(args.release))
    from solv_ai.features import descriptor_frame
    overrides = {k: args.release / "results/confirmatory/teacher_refits" / k / "teacher_predictions.parquet"
                 for k in ("combisolv_qm", "molsolv_smd", "confsolv")}
    d = load_confirmatory_data(args.workspace, overrides)
    extpath = args.release / "results/tier_a_external/evaluation/tier_a_external_predictions.parquet"
    e = pd.read_parquet(extpath)
    rows = []
    for frame, part in ((d.public, "public"), (d.benchmark, "arrow")):
        rows.append(pd.DataFrame({"id": frame.molecule_id, "smiles": frame.canonical_smiles,
                                 "y": frame.delta_g_exp, "part": part,
                                 "fold": frame.fold_random if part == "arrow" else -1,
                                 "weight": 3.0 if part == "arrow" else 1.0,
                                 "strict": False, "similarity": np.nan}))
    rows.append(pd.DataFrame({"id": e.candidate_id, "smiles": e.canonical_smiles_evaluated,
                             "y": e.y_true, "part": "external", "fold": -1, "weight": 0.,
                             "strict": e.strict_response_source_disjoint,
                             "similarity": e.maximum_endpoint_training_morgan_tanimoto}))
    frame = pd.concat(rows, ignore_index=True)
    assert frame.id.is_unique and len(frame) == 1585
    descriptor_columns = [c for c in pd.read_parquet(args.workspace / "data/processed/rdkit_morgan_features.parquet")
                          if c.startswith(("rdkit__", "morgan2__"))]
    ex = descriptor_frame(e.canonical_smiles_evaluated.tolist())[descriptor_columns].to_numpy(np.float32)
    er = pd.read_parquet(args.release / "results/tier_a_external/evaluation/tier_a_external_response_features.parquet")
    er = er.set_index("candidate_id").loc[e.candidate_id, d.response_names].to_numpy(np.float32)
    x = np.vstack([d.public_structure, d.benchmark_structure, ex])
    r = np.vstack([d.public_responses["full"], d.benchmark_responses["full"], er])
    assert x.shape == (1585, 2265) and r.shape == (1585, 15)
    frame.to_csv(args.output / "molecules.csv", index=False)
    np.savez_compressed(args.output / "features.npz", structure=x, response=r)
    dump(args.output / "protocol.json", PROTOCOL)
    paths = [args.output / "molecules.csv", args.output / "features.npz", Path(__file__), extpath, *overrides.values()]
    dump(args.output / "input_hashes.json", {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def embeddings(args):
    frame = pd.read_csv(args.output / "molecules.csv", float_precision="round_trip")
    name = args.encoder
    output = args.output / f"{name}.npy"
    if output.exists():
        return
    old = pd.read_parquet(args.workspace / "data/processed" / f"{name}_embeddings.parquet")
    cols = [c for c in old if c.startswith(name + "__")]
    old = old.drop_duplicates("canonical_smiles").set_index("canonical_smiles")
    x = old.reindex(frame.smiles)[cols].to_numpy(np.float32, copy=True)
    missing = np.flatnonzero(~np.isfinite(x).all(axis=1))
    smiles = frame.smiles.iloc[missing].tolist()
    print(name, "cached", len(frame) - len(missing), "missing", len(missing), flush=True)
    if name == "molformer":
        import torch
        from transformers import AutoModel, AutoTokenizer
        model_id = "ibm-research/MoLFormer-XL-both-10pct"
        revision = "7b12d946c181a37f6012b9dc3b002275de070314"
        tok = AutoTokenizer.from_pretrained(model_id, revision=revision, trust_remote_code=True, local_files_only=True)
        model = AutoModel.from_pretrained(model_id, revision=revision, trust_remote_code=True, local_files_only=True).eval().cuda()
        with torch.inference_mode():
            for start in range(0, len(smiles), 32):
                batch = tok(smiles[start:start+32], padding=True, truncation=True, max_length=256, return_tensors="pt")
                batch = {k: v.cuda() for k, v in batch.items()}
                h = model(**batch).last_hidden_state
                mask = batch["attention_mask"].unsqueeze(-1).to(h.dtype)
                x[missing[start:start+32]] = ((h*mask).sum(1)/mask.sum(1).clamp_min(1)).cpu().numpy()
    else:
        from unimol_tools import UniMolRepr
        encoder = UniMolRepr(data_type="molecule", batch_size=32, remove_hs=False,
                             model_name="unimolv1", model_size="84m", use_cuda=True)
        x[missing] = np.asarray(encoder.get_repr(smiles))
    assert np.isfinite(x).all()
    np.save(output, x)
    dump(args.output / f"{name}_metadata.json", {"shape": x.shape, "cached_rows": len(frame)-len(missing),
        "new_rows": len(missing), "labels_used": False, "sha256": hashlib.sha256(output.read_bytes()).hexdigest()})


def partitions(frame):
    for fold in sorted(frame.loc[frame.part.eq("arrow"), "fold"].unique()):
        test = np.flatnonzero(frame.part.eq("arrow") & frame.fold.eq(fold))
        train = np.flatnonzero(frame.part.eq("public") | (frame.part.eq("arrow") & frame.fold.ne(fold)))
        yield f"arrow_{int(fold)}", train, test
    yield "external", np.flatnonzero(frame.part.ne("external")), np.flatnonzero(frame.part.eq("external"))


def computed(args):
    """Separately declared source-attribution control with verified component replay."""
    from confirmatory_common import load_confirmatory_data
    sys.path.insert(0, str(args.release))
    from solv_ai.features import descriptor_frame
    from solv_ai.teachers import _tree_predictions
    overrides = {k: args.release / "results/confirmatory/teacher_refits" / k / "teacher_predictions.parquet"
                 for k in ("combisolv_qm", "molsolv_smd", "confsolv")}
    d = load_confirmatory_data(args.workspace, overrides)
    parts = []
    for r in [d.public_responses, d.benchmark_responses]:
        parts.append(np.column_stack([r["computation_core"], r["smd"], r["confsolv"]]))
    x = np.vstack(parts)
    assert x.shape == (1365, 10) and np.isfinite(x).all()
    if not (args.output / "computed_response_protocol.json").exists():
        dump(args.output / "computed_response_protocol.json", {
        "declared": "2026-10-02 before fitting this contrast; post-review exploratory source attribution",
        "question": "Do solely computed source targets improve the same structural learner?",
        "features": "COSMOtherm; uncorrected OpenFF and GBn2; SMD; six ConfSolv summaries",
        "excluded": "all five empirical Abraham descriptors and both experimental residual corrections",
        "evaluation": "same fixed ARROW five folds, same labels, weights and seeds; no external extrapolation claimed",
        "estimator": PROTOCOL["tree"], "seeds": SEEDS,
        "feature_sha256": hashlib.sha256(x.tobytes()).hexdigest()})
    frame = pd.read_csv(args.output / "molecules.csv", float_precision="round_trip")
    ext = frame.loc[frame.part.eq("external")]
    desc = descriptor_frame(ext.smiles.tolist())
    md = args.release / "models/final"
    op = _tree_predictions(desc, md / "openff_teacher.joblib", ["openff23_dg", "openff23_exp_residual"])
    gb = _tree_predictions(desc, md / "implicit_teacher.joblib", ["gbn2_alchemical_dg", "gbn2_exp_residual"])
    cached = pd.read_parquet(args.release / "results/tier_a_external/evaluation/tier_a_external_response_features.parquet")
    cached = cached.set_index("candidate_id").loc[ext.id]
    # The frozen response cache was explicitly stored in float32; replay that conversion.
    replay = {"openff": float(np.max(np.abs((op["openff23_dg"]+op["openff23_exp_residual"]).astype(np.float32)-cached.openff_corrected))),
              "gbn2": float(np.max(np.abs((gb["gbn2_alchemical_dg"]+gb["gbn2_exp_residual"]).astype(np.float32)-cached.gbn2_corrected)))}
    assert max(replay.values()) < 1e-10, replay
    ex = np.column_stack([cached.combisolv_qm, op["openff23_dg"], gb["gbn2_alchemical_dg"], cached.smd_water,
                          cached[d.response_names[-6:]]]).astype(np.float32)
    np.save(args.output / "computed_response.npy", np.vstack([x, ex]))
    dump(args.output / "computed_response_external_protocol.json", {
        "declared": "2026-10-02 before fitting the external control, after ARROW source-control evaluation",
        "question": "Does the same computed-only descriptor set transfer to the frozen external cohort?",
        "changes_to_model_or_training": "none; same ten descriptors, all 1365 training labels and weights",
        "components": "recovered from original packaged surrogates; corrected sums replay before fitting",
        "component_sum_max_replay_difference": replay,
        "feature_sha256": hashlib.sha256((args.output / "computed_response.npy").read_bytes()).hexdigest()})


def trees(args):
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import make_pipeline
    frame = pd.read_csv(args.output / "molecules.csv", float_precision="round_trip")
    z = np.load(args.output / "features.npz")
    xs = {"structure": z["structure"], "solvai": np.column_stack([z["structure"], z["response"]])}
    if (args.output / "computed_response.npy").exists():
        xs["computed_only"] = np.column_stack([z["structure"], np.load(args.output / "computed_response.npy")])
    for name in ("molformer", "unimol"):
        enc = np.load(args.output / f"{name}.npy")
        xs[name + "_structure"] = np.column_stack([z["structure"], enc])
        xs[name + "_solvai"] = np.column_stack([z["structure"], enc, z["response"]])
    for name, x in xs.items():
        if args.tree_variants and name not in args.tree_variants:
            continue
        for part, train, test in partitions(frame):
            path = args.output / "predictions" / f"{name}_{part}.csv"
            if path.exists():
                continue
            pred = []
            for seed in SEEDS:
                model = make_pipeline(SimpleImputer(strategy="median", add_indicator=True),
                    ExtraTreesRegressor(**PROTOCOL["tree"], random_state=seed, n_jobs=6))
                model.fit(x[train], frame.y.iloc[train], extratreesregressor__sample_weight=frame.weight.iloc[train])
                pred.append(model.predict(x[test]))
            out = frame.iloc[test].copy()
            out["model"] = name
            out["prediction"] = np.mean(pred, axis=0)
            for seed, values in zip(SEEDS, pred):
                out[f"seed_{seed}"] = values
            out.to_csv(path, index=False)
            print("done", name, part, flush=True)


def graphs(args):
    import torch
    from lightning import pytorch as pl
    from chemprop import data, models, nn
    torch.set_num_threads(6)
    torch.set_float32_matmul_precision("high")
    frame = pd.read_csv(args.output / "molecules.csv", float_precision="round_trip")
    cfg = PROTOCOL["graph"]

    class ConstantAdamMPNN(models.MPNN):
        def configure_optimizers(self):
            return torch.optim.Adam(self.parameters(), lr=self.max_lr)

    class SelectEpoch(pl.Callback):
        def __init__(self, target_scale):
            self.best = float("inf")
            self.epoch = 1
            self.trace = []
            self.target_scale = target_scale
        def on_validation_epoch_end(self, trainer, module):
            if trainer.sanity_checking:
                return
            value = float(trainer.callback_metrics["val/mae"]) * self.target_scale
            epoch = trainer.current_epoch + 1
            self.trace.append([epoch, value])
            if value < self.best:
                self.best, self.epoch = value, epoch
            if epoch - self.epoch >= cfg["early_stopping_patience"]:
                trainer.should_stop = True

    def fit(train, valid, seed, variant, epochs, select):
        pl.seed_everything(seed, workers=True)
        def dataset(idx):
            result = data.MoleculeDataset([data.MoleculeDatapoint.from_smi(
                frame.smiles.iloc[i], y=np.array([frame.y.iloc[i]], np.float32),
                weight=float(frame.weight.iloc[i])) for i in idx])
            result.cache = True  # identical deterministic graphs, built once instead of each epoch
            return result
        trainset = dataset(train)
        scaler = trainset.normalize_targets()
        validset = dataset(valid) if len(valid) else None
        # Chemprop disables output unscaling during validation; use the TRAIN scaler.
        if validset is not None:
            validset.normalize_targets(scaler)
        if variant == "chemeleon":
            cp = torch.load(Path.home() / ".chemprop/chemeleon_mp.pt", weights_only=True)
            mp = nn.BondMessagePassing(**cp["hyper_parameters"])
            mp.load_state_dict(cp["state_dict"])
            lr = cfg["pretrained_lr"]
        else:
            mp = nn.BondMessagePassing(d_h=cfg["scratch_message_dim"], depth=cfg["scratch_depth"], dropout=cfg["dropout"])
            lr = cfg["scratch_lr"]
        ffn = nn.RegressionFFN(input_dim=mp.output_dim, hidden_dim=cfg["ffn_hidden"], n_layers=cfg["ffn_layers"],
                               dropout=cfg["dropout"], output_transform=nn.UnscaleTransform.from_standard_scaler(scaler))
        model = ConstantAdamMPNN(mp, nn.MeanAggregation(), ffn, metrics=[nn.MAE()], max_lr=lr)
        callback = SelectEpoch(float(scaler.scale_[0]))
        trainer = pl.Trainer(accelerator="gpu", devices=1, max_epochs=epochs, logger=False,
                             enable_checkpointing=False, enable_progress_bar=False, enable_model_summary=False,
                             num_sanity_val_steps=0, callbacks=[callback] if select else [],
                             limit_val_batches=1.0 if select else 0)
        loader = data.build_dataloader(trainset, batch_size=cfg["batch_size"], num_workers=0, shuffle=True, seed=seed)
        val_loader = data.build_dataloader(validset, batch_size=64, shuffle=False) if select else None
        trainer.fit(model, loader, val_loader)
        return model, trainer, callback

    for variant in args.graph_variants:
        for part, train, test in partitions(frame):
            path = args.output / "predictions" / f"{variant}_{part}.csv"
            if path.exists():
                continue
            preds = []
            for seed in SEEDS:
                run_path = args.output / "graph_runs" / f"{variant}_{part}_{seed}.json"
                if run_path.exists():
                    result = json.loads(run_path.read_text())
                    preds.append(np.array(result["predictions"]))
                    continue
                rng = np.random.default_rng(20261002 + seed)
                valid = np.concatenate([rng.permutation(train[frame.part.iloc[train].to_numpy() == s])[:max(1,round((frame.part.iloc[train] == s).sum()*.15))]
                                         for s in ("public", "arrow")])
                inner = np.setdiff1d(train, valid)
                assert not set(test) & (set(inner) | set(valid))
                model, trainer, cb = fit(inner, valid, seed, variant, cfg["epoch_cap"], True)
                epochs, trace = cb.epoch, cb.trace
                del model, trainer
                torch.cuda.empty_cache()
                model, trainer, _ = fit(train, [], seed, variant, epochs, False)
                testset = data.MoleculeDataset([data.MoleculeDatapoint.from_smi(frame.smiles.iloc[i]) for i in test])
                testset.cache = True
                pred = torch.cat(trainer.predict(model, data.build_dataloader(testset, batch_size=64, shuffle=False))).cpu().numpy().ravel()
                assert pred.shape == (len(test),) and np.isfinite(pred).all()
                dump(run_path, {"variant": variant, "partition": part, "seed": seed, "selected_epochs": epochs,
                    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    "inner_train_ids": frame.id.iloc[inner].tolist(), "inner_validation_ids": frame.id.iloc[valid].tolist(),
                    "outer_train_ids": frame.id.iloc[train].tolist(), "test_ids": frame.id.iloc[test].tolist(),
                    "validation_trace": trace, "predictions": pred.tolist()})
                preds.append(pred)
                print("done", variant, part, seed, "epochs", epochs, flush=True)
                del model, trainer
                torch.cuda.empty_cache()
            out = frame.iloc[test].copy()
            out["model"] = variant
            out["prediction"] = np.mean(preds, axis=0, dtype=np.float64)
            for seed, values in zip(SEEDS, preds):
                out[f"seed_{seed}"] = values.astype(np.float64)
            out.to_csv(path, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prepare", "embeddings", "computed", "trees", "graphs"])
    ap.add_argument("--workspace", type=Path, required=True)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--encoder", choices=["molformer", "unimol"])
    ap.add_argument("--graph-variants", nargs="+", choices=["dmpnn", "chemeleon"], default=["dmpnn", "chemeleon"])
    ap.add_argument("--tree-variants", nargs="+")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "predictions").mkdir(exist_ok=True)
    (args.output / "graph_runs").mkdir(exist_ok=True)
    globals()[args.stage](args)


if __name__ == "__main__":
    main()
