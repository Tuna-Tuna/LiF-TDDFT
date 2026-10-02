"""Ordered propagation implementing the ordered charge-state model."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


STATES = ("F_minus", "F_neutral", "F_plus")


def transition_matrix(
    p_det: float,
    p_cap: float,
    *,
    p_det_two: float = 0.0,
    f_plus_to_neutral: float = 1.0,
) -> np.ndarray:
    """Return the column-stochastic F-/F0/F+ encounter matrix.

    ``p_det`` is ``1-P0`` and ``p_det_two`` is the two-electron sector P2.
    The charge-state closure imposes unity F+ -> F0 at the next encounter.
    """

    p_det = float(p_det)
    p_cap = float(p_cap)
    p_det_two = float(p_det_two)
    p_plus_neutral = float(f_plus_to_neutral)
    values = np.asarray([p_det, p_cap, p_det_two, p_plus_neutral], dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError("transition probabilities must be finite")
    if not 0.0 <= p_det <= 1.0 or not 0.0 <= p_cap <= 1.0:
        raise ValueError("Pdet and Pcap must lie in [0,1]; clipping is forbidden")
    if not 0.0 <= p_det_two <= p_det:
        raise ValueError("P2 must lie in [0,Pdet]; clipping is forbidden")
    if not 0.0 <= p_plus_neutral <= 1.0:
        raise ValueError("P(F+->F0) must lie in [0,1]; clipping is forbidden")
    p_det_one = p_det - p_det_two
    return np.asarray([
        [1.0 - p_det, p_cap, 0.0],
        [p_det_one, 1.0 - p_cap, p_plus_neutral],
        [p_det_two, 0.0, 1.0 - p_plus_neutral],
    ])


def propagate_ordered_events(
    events: list[dict],
    detachment_probability: Callable[[dict], float] | None,
    capture_probability: Callable[[dict], float],
    initial_state: str = "F_neutral",
    *,
    detachment_sector_probabilities: Callable[[dict], tuple[float, float, float]] | None = None,
) -> dict[str, object]:
    if detachment_probability is None and detachment_sector_probabilities is None:
        raise ValueError("a detachment probability or sector callback is required")
    state = np.zeros(3, dtype=float)
    state[STATES.index(initial_state)] = 1.0
    history = []
    transient_f_plus_created = 0.0
    for event in sorted(events, key=lambda item: (item["event_time"], item["event_index"])):
        if detachment_sector_probabilities is not None:
            p0, p1, p2 = map(float, detachment_sector_probabilities(event))
            if min(p0, p1, p2) < -1.0e-12 or not np.isclose(
                p0 + p1 + p2, 1.0, atol=1.0e-12
            ):
                raise ValueError("detachment sectors must be nonnegative and sum to one")
            p_det = 1.0 - p0
            mean_loss = p1 + 2.0 * p2
        else:
            assert detachment_probability is not None
            p_det = float(detachment_probability(event))
            p0, p1, p2 = 1.0 - p_det, p_det, 0.0
            mean_loss = p_det
        p_cap = float(capture_probability(event))
        incoming = state.copy()
        transient_f_plus_created += incoming[0] * p2
        state = transition_matrix(p_det, p_cap, p_det_two=p2) @ state
        if not np.isclose(state.sum(), 1.0, atol=1.0e-12):
            raise RuntimeError("charge-state probability is not conserved")
        history.append({
            "event_index": event["event_index"],
            "mean_loss": mean_loss,
            "p_det_zero": p0,
            "p_det_one": p1,
            "p_det_two": p2,
            "p_det": p_det,
            "p_cap": p_cap,
            "state": state.tolist(),
        })
    return {
        "states": dict(zip(STATES, state.tolist())),
        "yield_final_hybrid": float(state[0]),
        "transient_f_plus_probability_created": float(transient_f_plus_created),
        "history": history,
    }
