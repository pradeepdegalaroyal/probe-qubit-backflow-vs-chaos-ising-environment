"""Out-of-time-order correlator: the standard quantum diagnostic of scrambling."""
from __future__ import annotations

import numpy as np
import torch

from .linalg import device, eigh, to_numpy


def otoc(
    h_env: np.ndarray, w_op: np.ndarray, v_op: np.ndarray, psi: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """F(t) = <psi| W(t) V W(t) V |psi> with W(t) = exp(iHt) W exp(-iHt)."""
    energies, vecs = eigh(h_env)
    dev = device()
    c128 = torch.complex128
    vt = vecs.T
    w_e = (vt @ torch.as_tensor(w_op, dtype=torch.float64, device=dev) @ vecs).to(c128)
    v_e = (vt @ torch.as_tensor(v_op, dtype=torch.float64, device=dev) @ vecs).to(c128)
    psi_e = vt.to(c128) @ torch.as_tensor(psi, dtype=c128, device=dev)
    t = torch.as_tensor(times, dtype=torch.float64, device=dev)
    out = torch.empty(len(times), dtype=c128, device=dev)
    chunk = 64
    for s in range(0, len(times), chunk):
        ph = torch.exp(-1j * t[s : s + chunk, None] * energies[None, :])  # (T, D)
        x = (v_e @ psi_e)[None, :].expand(ph.shape[0], -1)
        for _ in range(2):  # apply V-then-W(t) twice: W(t) V W(t) V |psi>
            x = torch.conj(ph) * ((ph * x) @ w_e.T)
            if _ == 0:
                x = x @ v_e.T
        out[s : s + chunk] = x @ torch.conj(psi_e)
    return to_numpy(out)
