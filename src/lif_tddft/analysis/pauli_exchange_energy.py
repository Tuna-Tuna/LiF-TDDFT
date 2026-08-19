"""Naming gate for an optional frozen Pauli/exchange decomposition.

This module does not compute the decomposition.  It prevents provisional
orthogonalization proxies from being mislabeled as Pauli energy or force.
"""

REQUIRED_GATES = (
    "reference_states_defined",
    "all_energy_terms_closed",
    "far_distance_zero",
    "orthogonalization_schemes_agree",
    "conditioning_safe",
    "continuous_curve_differentiable",
    "short_range_force_repulsive",
    "aligned_with_td_force_and_loss",
)


def allowed_label(gates: dict[str, bool]) -> str:
    missing = [name for name in REQUIRED_GATES if not gates.get(name, False)]
    if missing:
        return "frozen-fragment orthogonalization descriptor"
    return "frozen Pauli/exchange"


def assert_pauli_exchange_allowed(gates: dict[str, bool]) -> None:
    missing = [name for name in REQUIRED_GATES if not gates.get(name, False)]
    if missing:
        raise ValueError("Pauli/exchange label is blocked; failed gates: " + ", ".join(missing))
