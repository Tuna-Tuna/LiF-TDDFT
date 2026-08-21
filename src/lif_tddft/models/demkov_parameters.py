"""Load provenance-backed Demkov inputs and derive gamma at runtime."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from .capture_demkov import DemkovCapture, gamma_from_binding_energies


@dataclass(frozen=True)
class DemkovParameters:
    """Physical inputs used to construct the Demkov capture parameter."""

    fluorine_electron_affinity_ev: float
    madelung_energy_magnitude_ev: float
    ev_per_hartree: float

    def __post_init__(self) -> None:
        for label, value in (
            ("fluorine electron affinity", self.fluorine_electron_affinity_ev),
            ("Madelung energy magnitude", self.madelung_energy_magnitude_ev),
            ("eV per Hartree", self.ev_per_hartree),
        ):
            if value <= 0.0:
                raise ValueError(f"{label} must be positive")

    @property
    def projectile_binding_energy_hartree(self) -> float:
        """Return ``Ep=d_EF`` after conversion from eV to Hartree."""

        return self.fluorine_electron_affinity_ev / self.ev_per_hartree

    @property
    def target_binding_energy_hartree(self) -> float:
        """Return ``Et=d_EF+V_Mad`` after conversion to Hartree."""

        return (
            self.fluorine_electron_affinity_ev + self.madelung_energy_magnitude_ev
        ) / self.ev_per_hartree

    @property
    def gamma_bohr_inverse(self) -> float:
        """Derive gamma; no saved MATLAB value is read or accepted."""

        return gamma_from_binding_energies(
            self.target_binding_energy_hartree,
            self.projectile_binding_energy_hartree,
        )

    def capture_model(self) -> DemkovCapture:
        return DemkovCapture.from_binding_energies(
            self.target_binding_energy_hartree,
            self.projectile_binding_energy_hartree,
        )


def load_demkov_parameters(path: Path) -> DemkovParameters:
    """Load primary energy inputs from YAML and derive all other quantities."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    parameters = payload.get("parameters", {})
    forbidden = {"gamma_bohr_inverse", "gamma", "gama"} & set(parameters)
    if forbidden:
        raise ValueError(
            "Demkov gamma must be derived from Et and Ep, not configured directly: "
            f"{sorted(forbidden)}"
        )
    required = {
        "fluorine_electron_affinity_ev",
        "madelung_energy_magnitude_ev",
        "ev_per_hartree",
    }
    missing = required - set(parameters)
    if missing:
        raise ValueError(f"Demkov configuration lacks primary inputs: {sorted(missing)}")
    return DemkovParameters(
        fluorine_electron_affinity_ev=float(
            parameters["fluorine_electron_affinity_ev"]
        ),
        madelung_energy_magnitude_ev=float(
            parameters["madelung_energy_magnitude_ev"]
        ),
        ev_per_hartree=float(parameters["ev_per_hartree"]),
    )
