#!/usr/bin/env python3
"""Check split integrity, epoch selection and label-free inference inputs."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frame = pd.read_csv(args.output / "molecules.csv", float_precision="round_trip")
    assert frame.id.is_unique and len(frame) == 1585
    public = set(frame.loc[frame.part.eq("public"), "id"])
    arrow = set(frame.loc[frame.part.eq("arrow"), "id"])
    external = set(frame.loc[frame.part.eq("external"), "id"])
    assert len(public) == 1280 and len(arrow) == 85 and len(external) == 220
    assert not public & arrow and not (public | arrow) & external
    assert (frame.loc[frame.part.eq("public"), "weight"] == 1).all()
    assert (frame.loc[frame.part.eq("arrow"), "weight"] == 3).all()
    trace_count = 0
    for model in ["dmpnn", "chemeleon"]:
        for part in [f"arrow_{i}" for i in range(5)] + ["external"]:
            if part == "external":
                expected_test, expected_train = external, public | arrow
            else:
                fold = int(part[-1])
                expected_test = set(frame.loc[frame.part.eq("arrow") & frame.fold.eq(fold), "id"])
                expected_train = (public | arrow) - expected_test
            predictions = pd.read_csv(args.output / "predictions" / f"{model}_{part}.csv", float_precision="round_trip").set_index("id")
            for seed in (11, 29, 47):
                path = args.output / "graph_runs" / f"{model}_{part}_{seed}.json"
                run = json.loads(path.read_text())
                train, test = set(run["outer_train_ids"]), set(run["test_ids"])
                inner, valid = set(run["inner_train_ids"]), set(run["inner_validation_ids"])
                assert train == expected_train and test == expected_test
                assert not train & test and not inner & valid and inner | valid == train
                trace = np.array(run["validation_trace"])
                assert np.isfinite(trace).all() and len(trace) <= 80
                best_epoch = int(trace[np.argmin(trace[:, 1]), 0])
                assert run["selected_epochs"] == best_epoch
                actual = predictions.loc[run["test_ids"], f"seed_{seed}"].to_numpy()
                assert np.allclose(actual, run["predictions"], rtol=0, atol=1e-12)
                trace_count += 1
    validation = json.loads((args.output / "validation.json").read_text())
    assert all(v < 1e-10 for v in validation["control_replay_max_absolute_difference"].values())
    ext = json.loads((args.output / "computed_response_external_protocol.json").read_text())
    assert all(v == 0 for v in ext["component_sum_max_replay_difference"].values())
    result = {"status": "pass", "graph_seed_runs": trace_count,
              "outer_test_labels_entered_fitting_or_epoch_selection": False,
              "complete_inner_partitions_checked": True,
              "epoch_selection_traces_checked": True,
              "control_prediction_replay_checked": True,
              "computed_external_component_replay_checked": True,
              "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output / "audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
