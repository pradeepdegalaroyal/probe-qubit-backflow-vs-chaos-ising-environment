import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chaosprobe import *  # noqa

n = 10
times = np.arange(0, 40.01, 0.1)
ypl = np.array([1, 1j]) / np.sqrt(2)
psi = np.array([1.0 + 0j])
for _ in range(n):
    psi = np.kron(psi, ypl)
for kind in "XZ":
    b = site_operator(n, n // 2, kind)
    for hz in (0.0, 0.1, 0.3, 0.809):
        bf, pl, early = [], [], []
        for seed in range(5):
            rng = np.random.default_rng(1000 * n + seed)
            p = ChainParams(n_sites=n, hz=hz)
            h = build_hamiltonian(p, draw_fields(p, rng))
            c = np.abs(loschmidt_coherence(h, b, 0.4, psi, times))
            bf.append(backflow(c)); pl.append(plateau_statistic(c)); early.append(backflow(c[:101]))
        print(kind, hz, "bf %.3f+-%.3f" % (np.mean(bf), np.std(bf)), "plateau %.4f+-%.4f" % (np.mean(pl), np.std(pl)),
              "bf(t<10) %.3f" % np.mean(early), flush=True)
