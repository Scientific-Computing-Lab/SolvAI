#!/usr/bin/env python3
"""Verify complete post-review comparisons and build vector figures and tables."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

LABELS = {
    "structure": "Structural descriptors",
    "solvai": "SolvAI",
    "computed_only": "Computation-only + structure",
    "dmpnn": "D-MPNN (from scratch)",
    "chemeleon": "CheMeleon (fine-tuned)",
    "molformer_structure": "MoLFormer + structure",
    "molformer_solvai": "MoLFormer + responses",
    "unimol_structure": "Uni-Mol + structure",
    "unimol_solvai": "Uni-Mol + responses",
}


def latex_number(value):
    """Retain the reported precision while typesetting signs mathematically."""
    return rf"\ensuremath{{{value:.4f}}}"


def paired(a, b):
    d = np.asarray(a) - np.asarray(b)
    rng = np.random.default_rng(20260828)
    boot = np.concatenate([d[rng.integers(0, len(d), size=(1000, len(d)))].mean(1) for _ in range(100)])
    low, high = np.quantile(boot, [.025, .975])
    return dict(delta=float(d.mean()), low=float(low), high=float(high), n=len(d))


def canonical_graph_exports(output):
    """Aggregate authoritative seed arrays in float64, without refitting models.

    Neural predictions are float32 tensors; explicit float64 exports avoid losing
    significant digits when pandas serializes float32 columns with short decimals.
    """
    source = output / "molecules.csv"
    if source.exists():
        frame = pd.read_csv(source, float_precision="round_trip").set_index("id")
    else:
        # Frozen-artifact reproduction needs only released TEST metadata, not
        # redistributing third-party training-label tables.
        files = [output / "predictions" / f"structure_arrow_{i}.csv" for i in range(5)]
        files.append(output / "predictions/structure_external.csv")
        frame = pd.concat([pd.read_csv(p, float_precision="round_trip") for p in files], ignore_index=True)
        frame = frame.drop(columns=[c for c in frame if c == "model" or c == "prediction" or c.startswith("seed_")]).set_index("id")
        assert frame.index.is_unique and len(frame) == 305
    changes = {}
    for model in ["dmpnn", "chemeleon"]:
        for part in [f"arrow_{i}" for i in range(5)] + ["external"]:
            runs = [json.loads((output / "graph_runs" / f"{model}_{part}_{seed}.json").read_text()) for seed in (11,29,47)]
            ids = runs[0]["test_ids"]
            assert all(r["test_ids"] == ids for r in runs)
            pred = np.array([r["predictions"] for r in runs], dtype=np.float64)
            path = output / "predictions" / f"{model}_{part}.csv"
            old = pd.read_csv(path, float_precision="round_trip").set_index("id")
            difference = float(np.max(np.abs(old.loc[ids, "prediction"] - pred.mean(axis=0))))
            assert difference < 1e-5, (model, part, difference)
            out = frame.loc[ids].reset_index().copy()
            out["model"] = model
            out["prediction"] = pred.mean(axis=0)
            for seed, values in zip((11,29,47), pred):
                out[f"seed_{seed}"] = values
            out.to_csv(path, index=False)
            changes[path.name] = difference
    report = output / "neural_export_precision.json"
    if not report.exists():
        report.write_text(json.dumps({"source": "saved full-precision per-seed JSON arrays", "models_refitted": False,
                                      "maximum_mean_export_changes": changes}, indent=2)+"\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    canonical_graph_exports(args.output)
    frames = []
    for model in LABELS:
        for part in [f"arrow_{i}" for i in range(5)] + ["external"]:
            p = args.output / "predictions" / f"{model}_{part}.csv"
            if not p.is_file():
                raise FileNotFoundError(f"Incomplete campaign: {p}")
            d = pd.read_csv(p, float_precision="round_trip")
            assert (d.model == model).all() and np.isfinite(d.prediction).all()
            assert np.allclose(d.prediction, d[["seed_11", "seed_29", "seed_47"]].mean(axis=1), atol=1e-12, rtol=0)
            frames.append(d)
    allpred = pd.concat(frames, ignore_index=True)
    assert not allpred.duplicated(["id", "model"]).any()
    for model, d in allpred.groupby("model"):
        assert len(d) == 305 and (d.part == "arrow").sum() == 85
        assert (d.part.eq("external") & d.strict).sum() == 97
    original_ext = pd.read_parquet(args.release / "results/tier_a_external/evaluation/tier_a_external_predictions.parquet").set_index("candidate_id")
    original_arrow = pd.read_parquet(args.release / "results/confirmatory/standardized_exclusion_endpoint_predictions.parquet")
    original_arrow = original_arrow.loc[original_arrow.partition.eq("standardized_exclusion_primary")]
    replay = {}
    for model, col in [("structure", "structure_only_prediction"), ("solvai", "solvai_prediction")]:
        d = allpred.loc[allpred.model.eq(model) & allpred.part.eq("external")].set_index("id")
        delta = float(np.max(np.abs(d.prediction - original_ext.loc[d.index, col])))
        assert delta < 1e-10, (model, delta)
        replay[model] = delta
        code = "A_structure_only" if model == "structure" else "F_full_solvai"
        old = original_arrow.loc[original_arrow.method.eq(code)].set_index("molecule_id")
        new = allpred.loc[allpred.model.eq(model) & allpred.part.eq("arrow")].set_index("id")
        delta = float(np.max(np.abs(new.prediction - old.loc[new.index, "y_pred"])))
        assert delta < 1e-10, (model, "ARROW replay", delta)
        replay[model + "_arrow"] = delta
    allpred.to_csv(args.output / "all_predictions.csv", index=False)
    metrics, comparisons = [], []
    data = {}
    for cohort, mask in [("ARROW-85", allpred.part.eq("arrow")),
                         ("External-220", allpred.part.eq("external")),
                         ("Strict-97", allpred.part.eq("external") & allpred.strict)]:
        selected = allpred.loc[mask]
        data[cohort] = {}
        for model, d in selected.groupby("model"):
            d = d.sort_values("id").reset_index(drop=True)
            data[cohort][model] = d
            metrics.append({"cohort": cohort, "model": model, "n": len(d),
                "mae": mean_absolute_error(d.y, d.prediction),
                "rmse": mean_squared_error(d.y, d.prediction)**.5,
                "median_ae": np.median(np.abs(d.y-d.prediction)), "r2": r2_score(d.y, d.prediction)})
        pairs = [("solvai", m) for m in LABELS if m != "solvai"]
        pairs += [("molformer_solvai", "molformer_structure"), ("unimol_solvai", "unimol_structure"), ("computed_only", "structure")]
        for a, b in pairs:
            da, db = data[cohort][a], data[cohort][b]
            assert da.id.equals(db.id) and np.array_equal(da.y, db.y)
            comparisons.append({"cohort": cohort, "candidate": a, "reference": b,
                                **paired(np.abs(da.y-da.prediction), np.abs(db.y-db.prediction))})
    met, comp = pd.DataFrame(metrics), pd.DataFrame(comparisons)
    met.to_csv(args.output / "metrics.csv", index=False)
    comp.to_csv(args.output / "paired_comparisons.csv", index=False)

    computed = pd.concat([pd.read_csv(args.output / "predictions" / f"computed_only_arrow_{i}.csv",
                                      float_precision="round_trip") for i in range(5)], ignore_index=True)
    computed = computed.sort_values("id").reset_index(drop=True)
    structural = data["ARROW-85"]["structure"]
    assert computed.id.equals(structural.id) and np.array_equal(computed.y, structural.y)
    computation_control = {"n": len(computed), "mae": float(np.abs(computed.y-computed.prediction).mean()),
        **paired(np.abs(computed.y-computed.prediction), np.abs(structural.y-structural.prediction))}
    (args.output / "computed_response_comparison.json").write_text(json.dumps(computation_control, indent=2)+"\n")
    macro_names = {"mae":"ComputedOnlyMAE", "delta":"ComputedOnlyDelta", "low":"ComputedOnlyLow", "high":"ComputedOnlyHigh"}
    macros = [rf"\newcommand{{\{name}}}{{\ensuremath{{{computation_control[key]:.3f}}}}}" for key,name in macro_names.items()]
    (args.release / "paper/tables/journal_macros.tex").write_text("\n".join(macros)+"\n")

    strata = []
    bins = [0., .3, .5, .7, 1.00001]
    for left, right in zip(bins[:-1], bins[1:]):
        a, b = data["External-220"]["solvai"], data["External-220"]["structure"]
        mask = a.similarity.ge(left) & a.similarity.lt(right)
        if not mask.any():
            continue
        aa, bb = a.loc[mask], b.loc[mask]
        strata.append({"lower": left, "upper": min(right, 1.), "n": len(aa),
                       "structure_mae": np.abs(bb.y-bb.prediction).mean(),
                       "solvai_mae": np.abs(aa.y-aa.prediction).mean(),
                       **paired(np.abs(aa.y-aa.prediction), np.abs(bb.y-bb.prediction))})
    strata = pd.DataFrame(strata)
    strata.to_csv(args.output / "similarity_strata.csv", index=False)

    table = args.release / "paper/tables/journal_baselines.tex"
    lines = [r"\begin{tabular}{@{}lrrr@{}}", r"\toprule",
             r"Predictor & ARROW-85 & External-220 & Strict-97 \\", r"\midrule"]
    for model, label in LABELS.items():
        vals = [met.loc[met.model.eq(model)&met.cohort.eq(c), "mae"].item() for c in data]
        lines.append(label + " & " + " & ".join(f"{v:.3f}" for v in vals) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    table.write_text("\n".join(lines)+"\n")
    tabdir = args.release / "paper/supplementary/tables"
    lines = [r"\begin{longtable}{@{}llrrrr@{}}", r"\toprule",
             r"Set & Model & $N$ & MAE & RMSE & Median AE \\", r"\midrule\endhead"]
    for _, row in met.iterrows():
        lines.append(f"{row.cohort} & {LABELS[row.model]} & {row.n} & {row.mae:.4f} & {row.rmse:.4f} & {row.median_ae:.4f}" + r" \\")
    lines += [r"\bottomrule\end{longtable}"]
    (tabdir / "journal_baselines.tex").write_text("\n".join(lines)+"\n")
    lines = [r"\begin{longtable}{@{}llrrr@{}}", r"\toprule",
             r"Set & Contrast (candidate $-$ reference) & $\Delta$MAE & 95\% lower & 95\% upper \\", r"\midrule\endhead"]
    for _, row in comp.iterrows():
        lines.append(f"{row.cohort} & {LABELS[row.candidate]} $-$ {LABELS[row.reference]} & {latex_number(row.delta)} & {latex_number(row.low)} & {latex_number(row.high)}" + r" \\")
    lines += [r"\bottomrule\end{longtable}"]
    (tabdir / "journal_comparisons.tex").write_text("\n".join(lines)+"\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.9), layout="constrained")
    colors = {m: ("#008f80" if "solvai" in m else "#667783") for m in LABELS}
    ax = axes[0]
    records = []
    for cohort in data:
        for a,b,label in [("solvai","structure","Structural"),("molformer_solvai","molformer_structure","MoLFormer"),("unimol_solvai","unimol_structure","Uni-Mol")]:
            row = comp.loc[comp.cohort.eq(cohort)&comp.candidate.eq(a)&comp.reference.eq(b)].iloc[0]
            records.append((f"{cohort}: {label}",row))
    for i,(label,row) in enumerate(records):
        ax.errorbar(row.delta,i,xerr=[[row.delta-row.low],[row.high-row.delta]],fmt="o",color="#008f80",capsize=2)
    ax.axvline(0,color="#777777",lw=.8)
    ax.set_yticks(range(len(records)),[r[0] for r in records],fontsize=8)
    ax.invert_yaxis();ax.set_xlabel("MAE change after adding responses\n(kcal/mol)")
    ax.set_title("a  Benefit of response descriptors",loc="left",fontweight="bold",fontsize=9)
    ax=axes[1]
    x=np.arange(len(strata)); w=.34
    ax.bar(x-w/2,strata.structure_mae,w,label="Structural descriptors",color="#667783")
    ax.bar(x+w/2,strata.solvai_mae,w,label="SolvAI",color="#008f80")
    labs=[f"{r.lower:.1f}–{r.upper:.1f}\n(n={int(r.n)})" for _,r in strata.iterrows()]
    ax.set_xticks(x,labs);ax.set_xlabel("Nearest training-molecule similarity\n(Morgan Tanimoto)")
    ax.set_ylabel("External-cohort MAE (kcal/mol)")
    ax.set_title("b  Error and chemical distance",loc="left",fontweight="bold",fontsize=9)
    ax.legend(fontsize=8,frameon=False)
    fig.savefig(args.release / "paper/figures/main/journal_baselines.pdf")
    fig.savefig(args.release / "paper/figures/main/journal_baselines.svg")
    plt.close(fig)
    meta={"control_replay_max_absolute_difference":replay,"models":list(LABELS),"molecule_predictions":len(allpred),
          "protocol_sha256":hashlib.sha256((args.output/"protocol.json").read_bytes()).hexdigest(),
          "analysis_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/"validation.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(met.pivot(index="model",columns="cohort",values="mae").to_string())
    print(comp.to_string(index=False))
    print(strata.to_string(index=False))
    print("Computed-only control:", json.dumps(computation_control))


if __name__ == "__main__":
    main()
