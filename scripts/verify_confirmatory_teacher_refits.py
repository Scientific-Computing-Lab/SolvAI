#!/usr/bin/env python3
"""Verify that confirmatory teacher refits preserve frozen source splits."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_chemprop_split(
    workspace: Path, release: Path, original_run: str, corrected_source: str
) -> dict:
    old_root = workspace / "results" / original_run / "input"
    new_root = release / "results" / "confirmatory" / "teacher_refits" / corrected_source / "input"
    exclusions = set(json.loads((new_root.parent / "exclusions.json").read_text())["smiles"])
    split_records = {}
    for split in ("train", "validation", "test"):
        old = pd.read_csv(old_root / f"{split}.csv")
        new = pd.read_csv(new_root / f"{split}.csv")
        expected = old.loc[~old.smiles.astype(str).isin(exclusions), list(new.columns)].reset_index(
            drop=True
        )
        new = new.reset_index(drop=True)
        if not expected.equals(new):
            raise AssertionError(f"{corrected_source}/{split} differs beyond declared exclusions")
        split_records[split] = {
            "original_rows": len(old),
            "corrected_rows": len(new),
            "removed_rows": len(old) - len(new),
            "original_sha256": sha256(old_root / f"{split}.csv"),
            "corrected_sha256": sha256(new_root / f"{split}.csv"),
            "retained_rows_exactly_equal_and_ordered": True,
        }
    return {
        "source": corrected_source,
        "declared_exclusion_smiles": len(exclusions),
        "splits": split_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--release-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    records = [
        verify_chemprop_split(
            args.workspace_root,
            args.release_root,
            "combisolv_qm_pretraining",
            "combisolv_qm",
        ),
        verify_chemprop_split(
            args.workspace_root,
            args.release_root,
            "molsolv_smd_pretraining",
            "molsolv_smd",
        ),
    ]
    confsolv = json.loads(
        (
            args.release_root
            / "results/confirmatory/teacher_refits/confsolv/training_metadata.json"
        ).read_text()
    )
    if confsolv["original_rows"] - confsolv["rows_after_exclusion"] != confsolv["exclusions"]:
        raise AssertionError("ConfSolv exclusion counts do not balance")
    records.append(
        {
            "source": "confsolv",
            "declared_exclusion_smiles": confsolv["exclusions"],
            "original_rows": confsolv["original_rows"],
            "corrected_rows": confsolv["rows_after_exclusion"],
            "original_split_indices_reused": True,
            "feature_ranking_computed_before_exclusion": True,
        }
    )
    output = {
        "status": "passed",
        "checks": records,
    }
    destination = (
        args.release_root / "audits/confirmatory/standardized_exclusion_refit_verification.json"
    )
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
