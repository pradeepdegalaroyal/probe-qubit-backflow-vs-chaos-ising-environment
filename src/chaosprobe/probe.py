"""Pure-dephasing probe qubit: Loschmidt coherence, information backflow, shots.

The probe couples as  H = H_E + lam * sigma_z (probe) x B.  A probe prepared in
|+> keeps coherence  c(t) = <psi| U_-^dag(t) U_+(t) |psi>  with
U_pm = exp(-i (H_E +- lam B) t).  For pure dephasing the trace distance of the
optimal state pair equals |c(t)|, so the Breuer-Laine-Piilo backflow is the
sum of all increases of |c(t)|.
"""
from __future__ import annotations

import numpy as np
import torch

from .linalg import device, eigh, to_numpy


def loschmidt_coherence(
    h_env: np.ndarray, coupling: np.ndarray, lam: float, psi: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """Complex probe coherence c(t) for environment state ``psi`` at ``times``."""
    e_p, v_p = eigh(h_env + lam * coupling)
    e_m, v_m = eigh(h_env - lam * coupling)
    dev = device()
    psi_t = torch.as_tensor(psi, dtype=torch.complex128, device=dev)
    t = torch.as_tensor(times, dtype=torch.float64, device=dev)
    a = v_p.T.to(torch.complex128) @ psi_t
    b = v_m.T.to(torch.complex128) @ psi_t
    overlap = (v_m.T @ v_p).to(torch.complex128)  # M[m, n] = <m-|n+>
    out = torch.empty(len(times), dtype=torch.complex128, device=dev)
    chunk = 64
    for s in range(0, len(times), chunk):
        tt = t[s : s + chunk, None]
        a_t = torch.exp(-1j * tt * e_p[None, :]) * a[None, :]
        b_t = torch.exp(-1j * tt * e_m[None, :]) * b[None, :]
        out[s : s + chunk] = torch.sum(torch.conj(b_t) * (a_t @ overlap.T), dim=1)
    return to_numpy(out)


def backflow(magnitude: np.ndarray) -> float:
    """Information backflow: sum of the positive increments of |c(t)|."""
    d = np.diff(np.asarray(magnitude, dtype=float))
    return float(d[d > 0].sum())


def simulate_ramsey(
    coherence: np.ndarray, shots: int, gamma: float, times: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Finite-shot estimate of |c(t)| from sigma_x and sigma_y Ramsey readout.

    ``gamma`` is an extra Markovian probe dephasing rate (c -> c exp(-gamma t)).
    Each quadrature uses ``shots`` single-shot +-1 outcomes per time point.
    """
    c = np.asarray(coherence) * np.exp(-gamma * np.asarray(times))
    px = np.clip((1 + c.real) / 2, 0, 1)
    py = np.clip((1 + c.imag) / 2, 0, 1)
    x_hat = 2 * rng.binomial(shots, px) / shots - 1
    y_hat = 2 * rng.binomial(shots, py) / shots - 1
    return np.sqrt(x_hat**2 + y_hat**2)


def backflow_thresholded(magnitude: np.ndarray, shots: int, z: float = 2.0) -> float:
    """Backflow counting only the part of each rise above ``z`` shot-noise sigmas.

    A difference of two estimated magnitudes has standard deviation at most
    sqrt(2 / shots), which sets the threshold.
    """
    d = np.diff(np.asarray(magnitude, dtype=float)) - z * np.sqrt(2.0 / shots)
    return float(d[d > 0].sum())


def plateau_statistic(magnitude: np.ndarray, shots: int | None = None) -> float:
    """Late-time mean of |c|^2 over the second half of the record.

    With ``shots`` given, the shot-noise bias 2/shots of |c_hat|^2 is removed.
    """
    m = np.asarray(magnitude, dtype=float)
    late = m[len(m) // 2 :] ** 2
    bias = 0.0 if shots is None else 2.0 / shots
    return float(late.mean() - bias)
