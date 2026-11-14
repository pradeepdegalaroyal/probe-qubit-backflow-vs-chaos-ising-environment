"""Thin dense linear-algebra backend (PyTorch; GPU if present, else CPU)."""
from __future__ import annotations

import numpy as np
import torch

torch.set_num_threads(6)


def device() -> torch.device:
    """Return the compute device (CUDA when available)."""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def eigh(mat: np.ndarray) -> tuple[torch.Tensor, torch.Tensor]:
    """Eigen-decomposition of a real symmetric matrix; tensors stay on the device."""
    t = torch.as_tensor(mat, dtype=torch.float64, device=device())
    return torch.linalg.eigh(t)


def to_numpy(t: torch.Tensor) -> np.ndarray:
    """Move a tensor back to a NumPy array."""
    return t.detach().cpu().numpy()
