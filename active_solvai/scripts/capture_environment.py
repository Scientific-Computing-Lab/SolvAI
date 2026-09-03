#!/usr/bin/env python3
"""Capture the immutable Phase-0 host and software inventory."""

from __future__ import annotations

import importlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "active_solvai/data/manifests/environment.json"


def command(args: list[str]) -> str | None:
    try:
        return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def package_version(name: str) -> str | None:
    try:
        module = importlib.import_module(name)
    except ImportError:
        return None
    return str(getattr(module, "__version__", "installed-version-not-exposed"))


def main() -> None:
    tools = [
        "uv",
        "python3",
        "git",
        "gh",
        "cmake",
        "ninja",
        "make",
        "clang",
        "gcc",
        "nvcc",
        "docker",
        "latexmk",
        "pandoc",
        "sqlite3",
        "duckdb",
        "Arbalest",
        "gmx",
        "lmp",
        "namd3",
    ]
    packages = [
        "numpy",
        "pandas",
        "pyarrow",
        "scipy",
        "sklearn",
        "rdkit",
        "torch",
        "gpytorch",
        "openmm",
        "chemprop",
        "lightgbm",
    ]
    payload = {
        "captured_utc": datetime.now(UTC).isoformat(),
        "parent_commit": command(["git", "rev-parse", "HEAD"]),
        "branch": command(["git", "branch", "--show-current"]),
        "working_tree_porcelain": command(["git", "status", "--porcelain=v1"]),
        "host": {
            "platform": platform.platform(),
            "python": sys.version,
            "cpu_count_logical": os.cpu_count(),
            "lscpu": command(["lscpu"]),
            "memory": command(["free", "-h"]),
            "disk": command(["df", "-h", str(ROOT)]),
            "gpu": command(
                [
                    "nvidia-smi",
                    "--query-gpu=name,memory.total,memory.free,driver_version,temperature.gpu,power.limit",
                    "--format=csv,noheader",
                ]
            ),
        },
        "tools": {name: shutil.which(name) for name in tools},
        "packages": {name: package_version(name) for name in packages},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(OUT)


if __name__ == "__main__":
    main()
