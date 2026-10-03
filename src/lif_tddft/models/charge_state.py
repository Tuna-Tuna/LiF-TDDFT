"""PRA Eq. (10): independent encounters between neutral and negative F."""
from collections.abc import Callable
import numpy as np

STATES = ('F_minus', 'F_neutral')


def transition_matrix(p_det: float, p_cap: float) -> np.ndarray:
    values = np.asarray([p_det, p_cap], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0) or np.any(values > 1):
        raise ValueError('Pdet and Pcap must be finite and lie in [0,1]; clipping is forbidden')
    return np.asarray([[1.0-p_det, p_cap], [p_det, 1.0-p_cap]])


def propagate_ordered_events(
    events: list[dict],
    detachment_probability: Callable[[dict], float],
    capture_probability: Callable[[dict], float],
    initial_state: str = 'F_neutral',
) -> dict[str, object]:
    if initial_state not in STATES:
        raise ValueError('initial_state must be F_minus or F_neutral')
    state = np.zeros(2)
    state[STATES.index(initial_state)] = 1.0
    history = []
    for event in sorted(events, key=lambda item: (item['event_time'], item['event_index'])):
        p_det = float(detachment_probability(event))
        p_cap = float(capture_probability(event))
        state = transition_matrix(p_det, p_cap) @ state
        if not np.isclose(state.sum(), 1.0, rtol=0, atol=1e-12):
            raise RuntimeError('charge-state probability is not conserved')
        history.append(dict(event_index=event['event_index'], mean_loss=p_det,
                            p_det=p_det, p_cap=p_cap, state=state.tolist()))
    return dict(states=dict(zip(STATES, state.tolist())),
                yield_final_hybrid=float(state[0]), history=history)
