#!/usr/bin/env python3
"""Fail if Git tracks numerical datasets or generated research outputs."""

from __future__ import annotations

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_PREFIXES = (
    "data/example/",
    "data/experiment/",
    "data/processed/",
    "results/",
)


def tracked_paths() -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return [
        item.decode("utf-8", errors="surrogateescape")
        for item in completed.stdout.split(b"\0")
        if item
    ]


def main() -> int:
    forbidden = [
        path
        for path in tracked_paths()
        if path.replace("\\", "/").startswith(FORBIDDEN_PREFIXES)
    ]
    if forbidden:
        formatted = "\n".join(f"  - {path}" for path in forbidden)
        raise SystemExit(
            "Git tracks numerical data or generated outputs:\n"
            f"{formatted}\n"
            "Keep these files local and ignored."
        )
    print("repository payload contains no tracked numerical datasets or generated outputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
