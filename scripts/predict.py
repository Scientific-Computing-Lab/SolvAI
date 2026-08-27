#!/usr/bin/env python3
"""Predict hydration free energies from one or more SMILES strings."""

from __future__ import annotations

import argparse

from solv_ai import predict_smiles


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="+")
    args = parser.parse_args()
    predictions, spread = predict_smiles(args.smiles)
    for smiles, value, uncertainty in zip(args.smiles, predictions, spread, strict=True):
        print(f"{smiles}\t{value:.6f}\t{uncertainty:.6f}")


if __name__ == "__main__":
    main()
