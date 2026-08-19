"""Provenance helpers for auditable raw and derived results."""

from __future__ import annotations

import hashlib
import json
import platform
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(repo: str | Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def software_version(executable: str, *args: str) -> str | None:
    try:
        result = subprocess.run(
            [executable, *args], capture_output=True, text=True, check=True, timeout=20
        )
        return (result.stdout or result.stderr).strip().splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return None


def build_provenance(
    repo: str | Path,
    inputs: Iterable[str | Path],
    pseudopotentials: Iterable[str | Path],
    *,
    octopus_version: str | None,
    mpi_command: list[str],
) -> dict:
    def records(paths: Iterable[str | Path]) -> list[dict[str, str]]:
        result = []
        for item in paths:
            path = Path(item)
            if not path.is_file():
                raise FileNotFoundError(path)
            result.append({"path": str(path.resolve()), "sha256": sha256_file(path)})
        return result

    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(repo),
        "octopus_version": octopus_version,
        "mpi_command": mpi_command,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "inputs": records(inputs),
        "pseudopotentials": records(pseudopotentials),
    }


def write_json(path: str | Path, data: dict) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
