"""Campaign expansion with one independent history per target velocity."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable


@dataclass(frozen=True)
class RunSpec:
    campaign_id: str
    run_id: str
    height_bohr: float
    relative_velocity_au: float
    cluster_model: str
    variant: str
    projectile_velocity_au: float
    cluster_velocity_au: float
    config_path: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path)
    text = config_path.read_text(encoding="utf-8")
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml
        except ImportError as exc:
            raise RuntimeError(
                "non-JSON YAML configuration requires PyYAML; the bundled "
                "default campaign works with the Python standard library"
            ) from exc
        data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError("campaign configuration must be a mapping")
    data["_config_path"] = str(config_path.resolve())
    return data


def _number_token(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def expand_campaign(config: dict[str, Any]) -> list[RunSpec]:
    required = (
        "campaign_id",
        "heights_bohr",
        "relative_velocities_au",
        "cluster_models",
        "variants",
    )
    missing = [key for key in required if key not in config]
    if missing:
        raise ValueError(f"missing campaign keys: {', '.join(missing)}")
    heights = config['heights_bohr']
    if not isinstance(heights, list) or not heights:
        raise ValueError('heights_bohr is required: supply actual calculated nodes within 1.2-10 bohr')
    if any(not math.isfinite(float(h)) or not 1.2 <= float(h) <= 10.0 for h in heights):
        raise ValueError('PRA collision heights must lie in [1.2,10] bohr')
    if any(float(b) <= float(a) for a,b in zip(heights, heights[1:])):
        raise ValueError('height nodes must be strictly increasing')
    pf = float(config.get("projectile_fraction", -1.0))
    cf = float(config.get("cluster_fraction", 0.0))
    if abs(abs(pf - cf) - 1.0) > 1.0e-12:
        raise ValueError("projectile_fraction and cluster_fraction must give unit relative speed")
    if pf != -1.0 or cf != 0.0:
        raise ValueError("PRA flyby requires a fixed surface and projectile moving toward decreasing x")
    specs: list[RunSpec] = []
    for height in map(float, config["heights_bohr"]):
        if height <= 0:
            raise ValueError("height must be positive")
        for velocity in map(float, config["relative_velocities_au"]):
            if not math.isfinite(velocity) or not 0.1 <= velocity <= 0.5:
                raise ValueError("PRA target velocity must lie in [0.1,0.5] a.u.")
            for model in config["cluster_models"]:
                for variant in config["variants"]:
                    run_id = (
                        f"{config['campaign_id']}__{model}__{variant}"
                        f"__h{_number_token(height)}__v{_number_token(velocity)}"
                    )
                    specs.append(
                        RunSpec(
                            campaign_id=str(config["campaign_id"]),
                            run_id=run_id,
                            height_bohr=height,
                            relative_velocity_au=velocity,
                            cluster_model=str(model),
                            variant=str(variant),
                            projectile_velocity_au=pf * velocity,
                            cluster_velocity_au=cf * velocity,
                            config_path=str(config.get("_config_path", "")),
                        )
                    )
    if len({item.run_id for item in specs}) != len(specs):
        raise ValueError("campaign expansion produced duplicate run_id values")
    return specs


def ramp_velocities(target: float, initial: float, increment: float) -> list[float]:
    if not (0 < initial <= target) or increment <= 0:
        raise ValueError("ramp requires 0 < initial <= target and increment > 0")
    values: list[float] = []
    current = initial
    while current < target - 1.0e-12:
        values.append(round(current, 12))
        current += increment
    values.append(float(target))
    return values


def independent_histories(specs: Iterable[RunSpec]) -> bool:
    """Return True when each run id identifies exactly one target velocity."""
    histories: dict[str, float] = {}
    for spec in specs:
        previous = histories.setdefault(spec.run_id, spec.relative_velocity_au)
        if previous != spec.relative_velocity_au:
            return False
    return True
