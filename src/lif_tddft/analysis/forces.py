"""Force naming and alignment; deliberately contains no Pauli-force alias."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ForceSeries:
    coordinate: np.ndarray
    force_full_electronic: np.ndarray
    force_point_charge_reference: np.ndarray

    def __post_init__(self) -> None:
        coordinate = np.asarray(self.coordinate, dtype=float)
        full = np.asarray(self.force_full_electronic, dtype=float)
        pc = np.asarray(self.force_point_charge_reference, dtype=float)
        if coordinate.ndim != 1 or full.shape != pc.shape or full.shape[0] != coordinate.size:
            raise ValueError("coordinate and force arrays are inconsistent")

    @property
    def force_full_minus_pc(self) -> np.ndarray:
        """Non-unique electronic contribution; not a Pauli-force decomposition."""
        return np.asarray(self.force_full_electronic) - np.asarray(self.force_point_charge_reference)


def impulse(times: np.ndarray, forces: np.ndarray) -> np.ndarray:
    return np.trapz(np.asarray(forces, dtype=float), np.asarray(times, dtype=float), axis=0)
