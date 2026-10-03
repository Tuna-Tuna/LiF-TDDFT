"""Load the explicitly supplied capture parameter of PRA Eq. (9)."""
from dataclasses import dataclass
from pathlib import Path
import math
import yaml
from .capture_demkov import DemkovCapture


@dataclass(frozen=True)
class DemkovParameters:
    gamma_bohr_inverse: float
    source: str

    def __post_init__(self):
        if not math.isfinite(self.gamma_bohr_inverse) or self.gamma_bohr_inverse <= 0:
            raise ValueError('gamma_c must be finite and positive')
        if not self.source.strip():
            raise ValueError('the source of gamma_c must be specified')

    def capture_model(self):
        return DemkovCapture(self.gamma_bohr_inverse)


def load_demkov_parameters(path: Path) -> DemkovParameters:
    payload = yaml.safe_load(path.read_text(encoding='utf-8'))
    parameters = payload.get('parameters', {})
    value = parameters.get('gamma_bohr_inverse')
    source = payload.get('gamma_source')
    if value is None or not isinstance(source, str) or not source.strip():
        raise ValueError('gamma_c and gamma_source are required; the PRA manuscript gives no numerical value')
    return DemkovParameters(float(value), source)
