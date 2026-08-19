"""Source-resolved numerical convergence summaries."""

from __future__ import annotations

import numpy as np


def relative_change(reference: float, candidate: float) -> float:
    scale = max(abs(float(reference)), 1.0e-15)
    return abs(float(candidate) - float(reference)) / scale


def summarize(reference: dict[str, float], variants: dict[str, dict[str, float]]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for source, values in variants.items():
        missing = set(reference) - set(values)
        if missing:
            raise ValueError(f"{source} lacks metrics: {sorted(missing)}")
        for metric, ref in reference.items():
            rows.append(
                {
                    "source": source,
                    "metric": metric,
                    "reference": float(ref),
                    "candidate": float(values[metric]),
                    "absolute_change": abs(float(values[metric]) - float(ref)),
                    "relative_change": relative_change(ref, values[metric]),
                }
            )
    return rows


def source_resolved_quadrature(components: dict[str, float]) -> float:
    return float(np.sqrt(sum(float(value) ** 2 for value in components.values())))
