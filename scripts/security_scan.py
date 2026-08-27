#!/usr/bin/env python3
"""Fail when high-confidence credential material occurs in the release tree."""

from __future__ import annotations

import json
import re
import subprocess
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


def is_scannable(path: Path) -> bool:
    return path.suffix.lower() in TEXT_SUFFIXES or path.name in {"Makefile", "LICENSE"}


def matches(text: str, location: str) -> list[dict[str, str]]:
    return [
        {"type": name, "file": location}
        for name, pattern in PATTERNS.items()
        if pattern.search(text)
    ]


def scan_history() -> tuple[int, list[dict[str, str]]]:
    """Inspect text blobs in all local revisions without echoing matched content."""

    listing = subprocess.run(
        ["git", "rev-list", "--objects", "--all"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    findings: list[dict[str, str]] = []
    seen: set[str] = set()
    scanned = 0
    for line in listing:
        object_id, separator, name = line.partition(" ")
        if not separator or object_id in seen or not is_scannable(Path(name)):
            continue
        seen.add(object_id)
        blob = subprocess.run(
            ["git", "cat-file", "-p", object_id],
            cwd=ROOT,
            check=True,
            capture_output=True,
        ).stdout
        scanned += 1
        findings.extend(matches(blob.decode(errors="ignore"), f"git:{object_id[:12]}:{name}"))
    return scanned, findings


def main() -> None:
    findings = []
    scanned = 0
    for path in ROOT.rglob("*"):
        if not path.is_file() or any(part in EXCLUDED_PARTS for part in path.parts):
            continue
        if not is_scannable(path):
            continue
        scanned += 1
        text = path.read_text(errors="ignore")
        findings.extend(matches(text, str(path.relative_to(ROOT))))
    _, history_findings = scan_history()
    findings.extend(history_findings)
    report = {
        "status": "PASS" if not findings else "FAIL",
        "files_scanned": scanned,
        "git_history_scanned": True,
        "findings": findings,
    }
    (ROOT / "audits/security_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    if findings:
        raise AssertionError(f"Potential credentials in {len(findings)} release files")
    print(f"Security scan PASS ({scanned} text files)")


if __name__ == "__main__":
    main()
