"""Level-spacing statistics: the classical exact-diagonalisation reference."""
from __future__ import annotations

import numpy as np

R_POISSON = 0.3863  # 2 ln 2 - 1, integrable (Poisson) value
R_GOE = 0.5307  # Gaussian orthogonal ensemble value, chaotic


def mean_gap_ratio(energies: np.ndarray, keep: float = 0.5) -> float:
    """Mean adjacent-gap ratio <r> over the central ``keep`` fraction of the spectrum."""
    e = np.sort(np.asarray(energies, dtype=float))
    n = len(e)
    lo = int(n * (1 - keep) / 2)
    hi = n - lo
    gaps = np.diff(e[lo:hi])
    a, b = gaps[:-1], gaps[1:]
    big = np.maximum(a, b)
    ok = big > 1e-13
    return float(np.mean(np.minimum(a, b)[ok] / big[ok]))
