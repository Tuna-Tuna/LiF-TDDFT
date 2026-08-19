#!/usr/bin/env python3
"""Build a checksum manifest and enforce release naming boundaries."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.provenance import git_commit, sha256_file

FORBIDDEN = {
    "pauli-free": re.compile(r"pauli[-_ ]free", re.I),
    "pauli-off": re.compile(r"pauli[-_ ]off", re.I),
    "force_pauli": re.compile(r"force_pauli", re.I),
    "first-principles final yield": re.compile(r"first[- ]principles final yield", re.I),
}


def audit_text(path: Path) -> list[str]:
    if path.name in {"build_release_manifest.py", "README.md", "README_revision.md"}:
        return []
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return []
    return [name for name, pattern in FORBIDDEN.items() if pattern.search(text)]


def build(root: Path) -> dict:
    files = []
    violations = []
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    relative_paths = sorted(
        item for item in result.stdout.decode("utf-8").split("\0") if item
    )
    for relative in relative_paths:
        if relative == "results/manifests/release_manifest.json":
            continue
        path = root / relative
        if not path.is_file():
            continue
        files.append({"path": relative, "sha256": sha256_file(path), "bytes": path.stat().st_size})
        hits = audit_text(path)
        if hits:
            violations.append({"path": relative, "terms": hits})
    return {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "base_git_commit": git_commit(root),
        "content_identity": "per-file SHA-256; the base commit may precede this manifest update",
        "files": files,
        "naming_violations": violations,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "manifests" / "release_manifest.json")
    args = parser.parse_args(argv)
    manifest = build(ROOT)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if manifest["naming_violations"]:
        print(json.dumps(manifest["naming_violations"], indent=2), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
