"""Nonresonant Demkov capture of Supporting Information Eq. (S14)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DemkovCapture:
    """Evaluate the SI Demkov probability without empirical envelopes."""

    gamma_bohr_inverse: float

    def __post_init__(self) -> None:
        if self.gamma_bohr_inverse <= 0:
            raise ValueError("Demkov gamma must be positive and provenance-backed")

    @staticmethod
    def _sech_squared(value: float) -> float:
        """Evaluate sech squared with an algebraically exact stable identity."""

        magnitude = abs(float(value))
        exponential = np.exp(-2.0 * magnitude)
        return float(4.0 * exponential / (1.0 + exponential) ** 2)

    def probability(self, event: dict) -> float:
        """Return ``1/2 sech^2[pi(DeltaE+v^2/2)/(2 gamma v)]``.

        ``energy_defect_au`` and ``v_parallel`` are mandatory event fields.
        Surface height enters only through the independently evaluated energy
        defect (SI Eqs. S15--S17); no additional height envelope is permitted.
        """

        energy_defect = float(event["energy_defect_au"])
        speed = float(event["v_parallel"])
        if not np.isfinite(energy_defect):
            raise ValueError("energy_defect_au must be finite")
        if not np.isfinite(speed) or speed <= 0.0:
            raise ValueError("v_parallel must be finite and positive")
        argument = np.pi * (energy_defect + 0.5 * speed**2) / (
            2.0 * self.gamma_bohr_inverse * speed
        )
        return 0.5 * self._sech_squared(argument)

    def sensitivity(self, event: dict, relative_span: float = 0.2) -> dict[str, tuple[float, float]]:
        """Re-evaluate S14 at alternate physical inputs for diagnostics only.

        This method does not multiply, normalize, or replace a computed
        production probability and is not called by figure or yield exports.
        """

        if not 0.0 < relative_span < 1.0:
            raise ValueError("relative_span must lie between zero and one")
        gamma_values = [
            DemkovCapture(self.gamma_bohr_inverse * factor).probability(event)
            for factor in (1.0 - relative_span, 1.0 + relative_span)
        ]
        energy = float(event["energy_defect_au"])
        energy_values = []
        for factor in (1.0 - relative_span, 1.0 + relative_span):
            varied = dict(event, energy_defect_au=energy * factor)
            energy_values.append(self.probability(varied))
        return {
            "gamma_bohr_inverse": (min(gamma_values), max(gamma_values)),
            "energy_defect_au": (min(energy_values), max(energy_values)),
        }
