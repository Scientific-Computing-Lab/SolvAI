"""Run prepared dense PIMD2 windows sequentially and account for every attempt."""

from __future__ import annotations

import argparse
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from active_solvai.ledger import append_record, sha256
from active_solvai.probes import read_energy

RELEASE_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_ROOT = RELEASE_ROOT / "active_solvai"
DEFAULT_ROOT = ACTIVE_ROOT / "simulations/dense_pimd2"
DEFAULT_BINARY = (
    RELEASE_ROOT.parents[1] / "repositories/arbalest/build_solvation_gcc10/ARBALEST/ARBALEST"
)


def slug(value: str) -> str:
    return "".join(character.lower() if character.isalnum() else "_" for character in value)


def system_energy(case: Path) -> Path | None:
    candidates = sorted((case / "output").glob("*__SYSTEM.ene"))
    return candidates[0] if len(candidates) == 1 else None


def quality_control(case: Path, return_code: int, console: str) -> dict[str, object]:
    energy = system_energy(case)
    result: dict[str, object] = {
        "success_marker": "Simulation has successfully finished!" in console,
        "system_energy_file": str(energy) if energy else None,
        "frames": 0,
        "finite_response": False,
        "mean_temperature_kelvin": None,
        "mean_density_g_cm3": None,
    }
    if energy is not None:
        frame = read_energy(energy)
        result["frames"] = len(frame)
        response = frame.get("dHdL", pd.Series(dtype=float)).to_numpy(dtype=float)
        result["finite_response"] = bool(len(response) and np.isfinite(response).all())
        if "Temp" in frame:
            result["mean_temperature_kelvin"] = float(frame.Temp.mean())
        if "Density" in frame:
            result["mean_density_g_cm3"] = float(frame.Density.mean())
        result["energy_sha256"] = sha256(energy)
    temperature = result["mean_temperature_kelvin"]
    density = result["mean_density_g_cm3"]
    result["passed"] = bool(
        return_code == 0
        and result["success_marker"]
        and result["frames"] >= 45
        and result["finite_response"]
        and temperature is not None
        and 250.0 <= temperature <= 350.0
        and density is not None
        and 0.70 <= density <= 1.30
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--role", choices=("calibration", "prospective"), required=True)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--omp", type=int, default=8)
    parser.add_argument("--timeout-minutes", type=float, default=30.0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    manifest = pd.read_csv(args.root / "manifest.csv")
    manifest = manifest.loc[(manifest.role == args.role) & ~manifest.existing_observation].copy()
    if args.limit is not None:
        manifest = manifest.head(args.limit)
    status_path = args.root / f"run_status_{args.role}.csv"
    status = pd.read_csv(status_path).to_dict("records") if status_path.exists() else []
    succeeded = {record["case_directory"] for record in status if bool(record.get("passed"))}
    ledger_path = ACTIVE_ROOT / "runs/ledger.jsonl"

    for row in manifest.to_dict("records"):
        case = Path(row["case_directory"])
        if str(case) in succeeded:
            print(
                f"SKIP {row['molecule_name']} lambda={row['lambda']}: complete",
                flush=True,
            )
            continue
        if sha256(Path(row["config"])) != row["config_sha256"]:
            raise AssertionError(f"Configuration hash changed: {row['config']}")
        prior_attempts = sum(record["case_directory"] == str(case) for record in status)
        if prior_attempts >= 2:
            print(
                f"SKIP {row['molecule_name']} lambda={row['lambda']}: two attempts used",
                flush=True,
            )
            continue
        attempt = prior_attempts + 1
        run_id = (
            f"AS-P2-{args.role.upper()}-{slug(row['molecule_name'])}-"
            f"L{int(row['ti_point']):02d}-A{attempt:02d}"
        )
        command = [
            str(args.binary.resolve()),
            "--gpu",
            "true",
            "--gpudeviceid",
            "0",
            "--omp",
            str(args.omp),
            "--outTS",
            "0",
            "--log",
            f"Arbalest_attempt_{attempt:02d}.log",
            "-C",
            "conf.xml",
        ]
        started = time.perf_counter()
        timeout = False
        console_path = case / f"console_attempt_{attempt:02d}.txt"
        with console_path.open("w") as console_handle:
            try:
                completed = subprocess.run(
                    command,
                    cwd=case,
                    stdout=console_handle,
                    stderr=subprocess.STDOUT,
                    timeout=args.timeout_minutes * 60.0,
                    check=False,
                )
                return_code = completed.returncode
            except subprocess.TimeoutExpired:
                timeout = True
                return_code = 124
        elapsed = time.perf_counter() - started
        console = console_path.read_text(errors="ignore")
        qc = quality_control(case, return_code, console)
        record = {
            "run_id": run_id,
            "timestamp_utc": datetime.now(UTC).isoformat(),
            "campaign": row["campaign"],
            "stage": "phase2_dense_calibration"
            if args.role == "calibration"
            else "phase2_dense_prospective",
            "role": args.role,
            "molecule_name": row["molecule_name"],
            "molecule_id": row["molecule_id"],
            "lambda": float(row["lambda"]),
            "ti_point": int(row["ti_point"]),
            "case_directory": str(case),
            "config_sha256": row["config_sha256"],
            "command": command,
            "attempt": attempt,
            "status": "completed" if qc["passed"] else "failed",
            "passed": bool(qc["passed"]),
            "exit_status": return_code,
            "timed_out": timeout,
            "quality_control": qc,
            "wall_seconds": elapsed,
            "gpu_hours": elapsed / 3600.0,
            "cpu_hours": elapsed * args.omp / 3600.0,
            "production_ps": float(row["production_ps"]),
            "simulated_time_ps": float(row["production_ps"]),
            "beads": int(row["beads"]),
            "bead_windows": int(row["beads"]),
            "nominal_bead_steps": 6250,
            "force_evaluations": None,
            "force_evaluation_note": "Arbalest does not expose exact fast/slow kernel call counts; nominal bead-steps are reported instead.",
            "freeze_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=RELEASE_ROOT, text=True
            ).strip(),
        }
        status.append(record)
        pd.DataFrame(status).to_csv(status_path, index=False)
        append_record(ledger_path, record)
        print(
            f"DONE {row['molecule_name']} lambda={row['lambda']}: "
            f"passed={qc['passed']} elapsed={elapsed:.1f}s",
            flush=True,
        )


if __name__ == "__main__":
    main()
