#!/usr/bin/env python3
"""Red-team manuscript claims and record every sensitive phrase in context."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = [ROOT / "paper/main.tex", ROOT / "README.md", ROOT / "submission/editorial_summary.md"]
PATTERNS = {
    "state_of_the_art": re.compile(r"state[- ]of[- ]the[- ]art", re.IGNORECASE),
    "robust_sub_threshold": re.compile(
        r"robust(?:ly)?[^.]{0,35}(?:sub-|below )0\.20", re.IGNORECASE
    ),
    "standard_arrow_benchmark": re.compile(r"standard ARROW benchmark", re.IGNORECASE),
    "universal": re.compile(r"\buniversal(?:ly)?\b", re.IGNORECASE),
    "first_ever": re.compile(r"\bfirst ever\b", re.IGNORECASE),
    "speedup_superlative": re.compile(r"(?:orders of magnitude|millions-fold)", re.IGNORECASE),
    "simulation_free_training": re.compile(r"simulation-free training", re.IGNORECASE),
    "single_neural_network": re.compile(r"single (?:monolithic )?neural network", re.IGNORECASE),
    "pimd_unnecessary": re.compile(r"PIMD is unnecessary", re.IGNORECASE),
    "exact_equivalence": re.compile(r"exactly equivalent", re.IGNORECASE),
}
SAFE_CONTEXT = (
    "not claimed",
    "not as",
    "no external",
    "not be read",
    "not a claim",
    "should be described",
    "has not been established",
    "rather than",
)


def main() -> None:
    findings = []
    unreviewed = []
    for path in FILES:
        text = path.read_text()
        for line_number, line in enumerate(text.splitlines(), 1):
            for claim, pattern in PATTERNS.items():
                if not pattern.search(line):
                    continue
                safe = any(marker in line.lower() for marker in SAFE_CONTEXT)
                record = {
                    "claim": claim,
                    "file": str(path.relative_to(ROOT)),
                    "line": line_number,
                    "text": line.strip(),
                    "disposition": "explicit limitation/negation" if safe else "unreviewed",
                }
                findings.append(record)
                if not safe:
                    unreviewed.append(record)
    report = {
        "status": "PASS" if not unreviewed else "FAIL",
        "findings": findings,
        "unreviewed": unreviewed,
    }
    (ROOT / "audits/claim_red_team.json").write_text(json.dumps(report, indent=2) + "\n")
    if unreviewed:
        raise AssertionError(f"Unreviewed high-risk claims: {len(unreviewed)}")
    print(f"Claim red-team PASS ({len(findings)} explicit limitations reviewed)")


if __name__ == "__main__":
    main()
