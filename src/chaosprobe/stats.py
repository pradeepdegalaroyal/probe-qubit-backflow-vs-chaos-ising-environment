"""Small statistics helpers."""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def auc_score(positives: np.ndarray, negatives: np.ndarray) -> float:
    """Area under the receiver operating characteristic curve (Mann-Whitney form).

    Probability that a random positive scores above a random negative (ties 1/2).
    """
    pos = np.asarray(positives, dtype=float)
    neg = np.asarray(negatives, dtype=float)
    ranks = rankdata(np.concatenate([pos, neg]))
    u = ranks[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2
    return float(u / (len(pos) * len(neg)))


def haar_state(dim: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-random pure state of dimension ``dim``."""
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return v / np.linalg.norm(v)
