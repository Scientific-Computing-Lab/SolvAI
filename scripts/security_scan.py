#!/usr/bin/env python3
"""Fail when high-confidence credential material occurs in the release tree."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_PARTS = {".git", ".venv", ".pytest_cache", ".ruff_cache", "__pycache__"}
TEXT_SUFFIXES = {
    ".bib",
    ".cff",
    ".csv",
    ".json",
    ".md",
    ".py",
    ".tex",
    ".toml",
    ".txt",
    ".yml",
    ".yaml",
}
PATTERNS = {
    "github_classic_pat": re.compile("gh" + "p_[A-Za-z0-9]{20,}"),
    "github_fine_grained_pat": re.compile("github" + "_pat_[A-Za-z0-9_]{20,}"),
    "openai_secret": re.compile("sk" + "-[A-Za-z0-9_-]{20,}"),
    "authorization_header": re.compile(
        "Author" + "ization:\\s*(?:Bearer|token)\\s+\\S+", re.IGNORECASE
    ),
    "assigned_secret": re.compile(
        r"(?:API_KEY|TOKEN|PASSWORD)\s*=\s*[^\s$<{][^\s]{7,}", re.IGNORECASE
    ),
    "private_ssh_key": re.compile("BEGIN OPENSSH " + "PRIVATE KEY"),
}


def main() -> None:
    findings = []
    scanned = 0
    for path in ROOT.rglob("*"):
        if not path.is_file() or any(part in EXCLUDED_PARTS for part in path.parts):
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name not in {"Makefile", "LICENSE"}:
            continue
        scanned += 1
        text = path.read_text(errors="ignore")
        for name, pattern in PATTERNS.items():
            if pattern.search(text):
                findings.append({"type": name, "file": str(path.relative_to(ROOT))})
    report = {
        "status": "PASS" if not findings else "FAIL",
        "files_scanned": scanned,
        "findings": findings,
    }
    (ROOT / "audits/security_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    if findings:
        raise AssertionError(f"Potential credentials in {len(findings)} release files")
    print(f"Security scan PASS ({scanned} text files)")


if __name__ == "__main__":
    main()
