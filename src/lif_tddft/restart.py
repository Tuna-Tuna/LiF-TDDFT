"""Conservative restart validation."""

from pathlib import Path


def validate_restart(stage_dir: str | Path, *, require_coordinates: bool = True) -> dict[str, object]:
    root = Path(stage_dir)
    restart = root / "restart"
    coordinates = root / "td.general" / "coordinates"
    issues: list[str] = []
    if not restart.is_dir() or not any(restart.iterdir()):
        issues.append("restart directory is missing or empty")
    if require_coordinates and (not coordinates.is_file() or coordinates.stat().st_size == 0):
        issues.append("td.general/coordinates is missing or empty")
    return {"valid": not issues, "issues": issues}
