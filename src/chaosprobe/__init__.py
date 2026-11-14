"""chaosprobe: a dephasing probe qubit as a diagnostic of many-body quantum chaos."""

from .models import ChainParams, build_hamiltonian, draw_fields, site_operator
from .otoc import otoc
from .probe import (
    backflow,
    backflow_thresholded,
    loschmidt_coherence,
    plateau_statistic,
    simulate_ramsey,
)
from .spectral import mean_gap_ratio
from .stats import auc_score, haar_state

__version__ = "1.0.0"
__all__ = [
    "ChainParams",
    "build_hamiltonian",
    "site_operator",
    "draw_fields",
    "mean_gap_ratio",
    "loschmidt_coherence",
    "backflow",
    "simulate_ramsey",
    "backflow_thresholded",
    "plateau_statistic",
    "otoc",
    "auc_score",
    "haar_state",
]
