"""Explicit unit conversions used at input/output boundaries."""

BOHR_PER_ANGSTROM = 1.8897261246257702
HARTREE_PER_EV = 1.0 / 27.211386245988
ATOMIC_TIME_PER_FS = 41.34137333518211


def angstrom_to_bohr(value: float) -> float:
    return float(value) * BOHR_PER_ANGSTROM


def bohr_to_angstrom(value: float) -> float:
    return float(value) / BOHR_PER_ANGSTROM


def ev_to_hartree(value: float) -> float:
    return float(value) * HARTREE_PER_EV


def fs_to_atomic_time(value: float) -> float:
    return float(value) * ATOMIC_TIME_PER_FS
