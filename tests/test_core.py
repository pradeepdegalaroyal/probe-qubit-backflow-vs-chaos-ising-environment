"""Fast correctness tests (small chains, run on CPU or GPU)."""
import numpy as np
import scipy.linalg as sl

from chaosprobe import (
    ChainParams,
    auc_score,
    backflow,
    backflow_thresholded,
    build_hamiltonian,
    draw_fields,
    haar_state,
    loschmidt_coherence,
    mean_gap_ratio,
    otoc,
    simulate_ramsey,
    site_operator,
)


def _setup(n=5, hz=0.4, seed=0):
    rng = np.random.default_rng(seed)
    p = ChainParams(n_sites=n, hz=hz)
    h = build_hamiltonian(p, draw_fields(p, rng))
    return h, haar_state(2**n, rng)


def test_hamiltonian_matches_kron():
    p = ChainParams(n_sites=3, hz=0.3, hx_disorder=0.0)
    h = build_hamiltonian(p, np.full(3, 0.9045))
    x = np.array([[0, 1], [1, 0.0]])
    z = np.diag([1.0, -1.0])
    i = np.eye(2)

    def k(a, b, c):  # site 0 is the least significant bit
        return np.kron(c, np.kron(b, a))

    ref = (
        k(z, z, i) + k(i, z, z)
        + 0.9045 * (k(x, i, i) + k(i, x, i) + k(i, i, x))
        + 0.3 * (k(z, i, i) + k(i, z, i) + k(i, i, z))
    )
    assert np.allclose(h, ref)
    assert np.allclose(site_operator(3, 1, "X"), k(i, x, i))
    assert np.allclose(site_operator(3, 2, "Z"), k(i, i, z))


def test_coherence_matches_expm_and_is_bounded():
    h, psi = _setup()
    b = site_operator(5, 2, "X")
    times = np.linspace(0, 3, 7)
    c = loschmidt_coherence(h, b, 0.4, psi, times)
    for t, ct in zip(times, c):
        up = sl.expm(-1j * (h + 0.4 * b) * t) @ psi
        um = sl.expm(-1j * (h - 0.4 * b) * t) @ psi
        assert abs(np.vdot(um, up) - ct) < 1e-9
    assert abs(c[0] - 1) < 1e-12
    assert np.all(np.abs(c) <= 1 + 1e-9)


def test_zero_coupling_gives_no_decay():
    h, psi = _setup()
    c = loschmidt_coherence(h, site_operator(5, 2, "X"), 0.0, psi, np.linspace(0, 5, 6))
    assert np.allclose(c, 1.0)


def test_otoc_starts_at_one_and_matches_expm():
    h, psi = _setup()
    w, v = site_operator(5, 0, "X"), site_operator(5, 4, "X")
    f = otoc(h, w, v, psi, np.array([0.0, 1.7]))
    assert abs(f[0] - 1) < 1e-10
    u = sl.expm(-1j * h * 1.7)
    wt = u.conj().T @ w @ u
    assert abs(np.vdot(psi, wt @ v @ wt @ v @ psi) - f[1]) < 1e-9


def test_backflow_and_threshold():
    assert abs(backflow(np.array([1.0, 0.5, 0.7, 0.6, 0.65])) - 0.25) < 1e-12
    assert backflow(np.exp(-np.linspace(0, 3, 20))) == 0.0
    assert backflow_thresholded(np.array([0.5, 0.7]), shots=200, z=2.0) == 0.0
    assert backflow_thresholded(np.array([0.5, 0.9]), shots=200, z=2.0) > 0.0


def test_ramsey_converges_with_shots():
    rng = np.random.default_rng(1)
    c = 0.6 * np.exp(1j * np.linspace(0, 2, 5))
    est = simulate_ramsey(c, 10**6, 0.0, np.zeros(5), rng)
    assert np.allclose(est, 0.6, atol=5e-3)


def test_gap_ratio_poisson_and_goe():
    rng = np.random.default_rng(2)
    assert abs(mean_gap_ratio(np.cumsum(rng.exponential(size=20000)), keep=1.0) - 0.386) < 0.02
    a = rng.normal(size=(400, 400))
    assert abs(mean_gap_ratio(np.linalg.eigvalsh(a + a.T)) - 0.531) < 0.04


def test_auc():
    assert auc_score(np.array([2.0, 3.0]), np.array([0.0, 1.0])) == 1.0
    assert auc_score(np.array([1.0]), np.array([1.0])) == 0.5
