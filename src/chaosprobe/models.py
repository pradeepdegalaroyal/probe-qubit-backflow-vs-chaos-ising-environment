"""Mixed-field Ising chain (the environment) and single-site operators.

H_E = J sum_i Z_i Z_{i+1} + sum_i hx_i X_i + hz sum_i Z_i   (open boundaries)

All matrices are real symmetric and dense; bit i of a basis index is spin i
(bit = 1 means Z = -1).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ChainParams:
    """Parameters of the environment chain."""

    n_sites: int
    j_zz: float = 1.0
    hx_mean: float = 0.9045
    hx_disorder: float = 0.1
    hz: float = 0.0


def draw_fields(params: ChainParams, rng: np.random.Generator) -> np.ndarray:
    """Transverse fields hx_i = hx_mean + uniform(-hx_disorder, hx_disorder).

    The weak disorder removes the reflection symmetry and accidental
    degeneracies but keeps the chain free-fermionic at hz = 0.
    """
    return params.hx_mean + rng.uniform(-params.hx_disorder, params.hx_disorder, params.n_sites)


def _z_signs(n_sites: int) -> np.ndarray:
    idx = np.arange(2**n_sites)
    bits = (idx[:, None] >> np.arange(n_sites)[None, :]) & 1
    return 1.0 - 2.0 * bits  # shape (D, n_sites)


def site_operator(n_sites: int, site: int, kind: str) -> np.ndarray:
    """Dense Pauli X or Z acting on ``site`` of an ``n_sites`` chain."""
    dim = 2**n_sites
    if kind == "Z":
        return np.diag(_z_signs(n_sites)[:, site])
    if kind == "X":
        op = np.zeros((dim, dim))
        idx = np.arange(dim)
        op[idx ^ (1 << site), idx] = 1.0
        return op
    raise ValueError(f"kind must be 'X' or 'Z', got {kind!r}")


def build_hamiltonian(params: ChainParams, hx: np.ndarray) -> np.ndarray:
    """Dense Hamiltonian of the mixed-field Ising chain with site fields ``hx``."""
    n = params.n_sites
    if len(hx) != n:
        raise ValueError("hx must have one entry per site")
    dim = 2**n
    z = _z_signs(n)
    diag = params.j_zz * np.sum(z[:, :-1] * z[:, 1:], axis=1) + params.hz * z.sum(axis=1)
    ham = np.diag(diag)
    idx = np.arange(dim)
    for site in range(n):
        ham[idx ^ (1 << site), idx] += hx[site]
    return ham
