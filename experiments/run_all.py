"""Reproduce every result and figure of the study with one command.

    python experiments/run_all.py           # full study (GPU: about 5 minutes)
    python experiments/run_all.py --quick   # smoke run (< 1 minute, CPU is fine)

Full results go to results/ and figures/; the smoke run writes to
results/quick/ and figures/quick/ so it never overwrites the published numbers.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chaosprobe import (  # noqa: E402
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
    plateau_statistic,
    simulate_ramsey,
    site_operator,
)
from chaosprobe import plotting  # noqa: E402
from chaosprobe.linalg import device, eigh, to_numpy  # noqa: E402

HZ_CHAOTIC = 0.809
LAM = 0.4
STATES = ("haar", "product")
ESTIMATORS = ("naive", "thresholded", "plateau")
DIAGNOSTICS = ("backflow_x_haar", "plateau_x_haar", "backflow_x_product", "plateau_x_product",
               "backflow_z_product", "plateau_z_product", "otoc_late")


def config(quick: bool) -> dict:
    """Experiment sizes for the full study and for the smoke run."""
    if quick:
        return dict(
            sweep_sizes=[6, 8], main_size=8, hz_grid=[0.0, 0.2, HZ_CHAOTIC], seeds=2,
            scaling_sizes=[6, 7, 8], scaling_seeds={}, lams=[0.2, 0.4],
            t_max=40.0, dt=0.2, shot_stride=5, shots=[100, 1000], gammas=[0.0, 0.05],
            reps=10, ecdf_shots=1000,
        )
    return dict(
        sweep_sizes=[8, 10], main_size=10,
        hz_grid=[0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, HZ_CHAOTIC], seeds=5,
        scaling_sizes=[6, 7, 8, 9, 10, 11, 12], scaling_seeds={12: 3}, lams=[0.2, 0.4, 0.8],
        t_max=40.0, dt=0.1, shot_stride=10, shots=[30, 100, 300, 1000, 3000, 10000],
        gammas=[0.0, 0.01, 0.02, 0.05, 0.1, 0.2], reps=40, ecdf_shots=1000,
    )


def instance(n: int, hz: float, seed: int) -> tuple[np.ndarray, np.random.Generator]:
    """Environment Hamiltonian of one disorder realisation and its generator."""
    rng = np.random.default_rng(1000 * n + seed)  # same fields for every hz
    params = ChainParams(n_sites=n, hz=hz)
    return build_hamiltonian(params, draw_fields(params, rng)), rng


def product_state(n: int) -> np.ndarray:
    """All spins along +y: a product state with <H_E> = Tr(H_E)/D for every hz.

    It sits at the infinite-temperature energy density, like a Haar state, but
    has no entanglement and is easy to prepare.
    """
    one = np.array([1.0, 1.0j]) / np.sqrt(2.0)
    psi = np.array([1.0 + 0.0j])
    for _ in range(n):
        psi = np.kron(psi, one)
    return psi


def env_state(name: str, n: int, rng: np.random.Generator) -> np.ndarray:
    """Initial environment state: 'haar' or 'product'."""
    return haar_state(2**n, rng) if name == "haar" else product_state(n)


def run_sweep(cfg: dict, times: np.ndarray) -> tuple[pd.DataFrame, dict]:
    """E1: integrable-to-chaotic sweep: probe statistics, OTOC and gap ratio."""
    rows, traces = [], {}
    for n in cfg["sweep_sizes"]:
        mid = n // 2
        bx, bz = site_operator(n, mid, "X"), site_operator(n, mid, "Z")
        w, v = site_operator(n, 0, "X"), site_operator(n, n - 1, "X")
        prod = product_state(n)
        for hz in cfg["hz_grid"]:
            for seed in range(cfg["seeds"]):
                h, rng = instance(n, hz, seed)
                psi = haar_state(2**n, rng)
                energies = to_numpy(eigh(h)[0])
                cxh = loschmidt_coherence(h, bx, LAM, psi, times)
                cxp = loschmidt_coherence(h, bx, LAM, prod, times)
                czp = np.abs(loschmidt_coherence(h, bz, LAM, prod, times))
                f = otoc(h, w, v, psi, times).real
                traces[("haar", n, hz, seed)] = cxh
                traces[("product", n, hz, seed)] = cxp
                rows.append(dict(
                    n_sites=n, hz=hz, seed=seed, gap_ratio=mean_gap_ratio(energies),
                    energy_width=float(energies.std()),
                    backflow_x_haar=backflow(np.abs(cxh)),
                    plateau_x_haar=plateau_statistic(np.abs(cxh)),
                    backflow_x_product=backflow(np.abs(cxp)),
                    plateau_x_product=plateau_statistic(np.abs(cxp)),
                    backflow_z_product=backflow(czp), plateau_z_product=plateau_statistic(czp),
                    otoc_late=float(f[len(f) // 2:].mean()),
                ))
    return pd.DataFrame(rows), traces


def run_scaling(cfg: dict, times: np.ndarray) -> pd.DataFrame:
    """E2: system-size scaling at the integrable and the chaotic end."""
    rows = []
    for n in cfg["scaling_sizes"]:
        bx = site_operator(n, n // 2, "X")
        for hz in (0.0, HZ_CHAOTIC):
            for seed in range(cfg["scaling_seeds"].get(n, cfg["seeds"])):
                h, rng = instance(n, hz, seed)
                for state in STATES:
                    c = np.abs(loschmidt_coherence(h, bx, LAM, env_state(state, n, rng), times))
                    rows.append(dict(n_sites=n, hz=hz, seed=seed, state=state,
                                     backflow=backflow(c), plateau=plateau_statistic(c)))
    return pd.DataFrame(rows)


def run_ablation(cfg: dict, times: np.ndarray) -> pd.DataFrame:
    """E3: coupling operator, coupling strength and initial environment state."""
    n = cfg["main_size"]
    rows = []
    for kind in ("X", "Z"):
        b = site_operator(n, n // 2, kind)
        for lam in cfg["lams"]:
            for state in STATES:
                for hz in (0.0, HZ_CHAOTIC):
                    for seed in range(cfg["seeds"]):
                        h, rng = instance(n, hz, seed)
                        c = np.abs(loschmidt_coherence(h, b, lam, env_state(state, n, rng), times))
                        rows.append(dict(coupling=kind, lam=lam, state=state, hz=hz, seed=seed,
                                         backflow=backflow(c), plateau=plateau_statistic(c)))
    return pd.DataFrame(rows)


def run_shots(cfg: dict, times: np.ndarray, traces: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """E4: finite-shot Ramsey readout with extra Markovian probe dephasing."""
    n, stride = cfg["main_size"], cfg["shot_stride"]
    t = times[::stride]
    rng = np.random.default_rng(2026)
    rows, samples = [], []
    for state in STATES:
        for gamma in cfg["gammas"]:
            for shots in cfg["shots"]:
                scores: dict[str, dict[float, list[float]]] = {
                    k: {0.0: [], HZ_CHAOTIC: []} for k in ESTIMATORS}
                for hz in (0.0, HZ_CHAOTIC):
                    for seed in range(cfg["seeds"]):
                        c = traces[(state, n, hz, seed)][::stride]
                        for _ in range(cfg["reps"]):
                            m = simulate_ramsey(c, shots, gamma, t, rng)
                            scores["naive"][hz].append(backflow(m))
                            scores["thresholded"][hz].append(backflow_thresholded(m, shots))
                            scores["plateau"][hz].append(plateau_statistic(m, shots))
                row = dict(state=state, gamma=gamma, shots=shots)
                for k, d in scores.items():
                    row[f"auc_{k}"] = auc_score(np.array(d[0.0]), np.array(d[HZ_CHAOTIC]))
                    if shots == cfg["ecdf_shots"] and gamma == 0.0 and state == "product":
                        for hz, vals in d.items():
                            samples += [dict(estimator=k, hz=hz, value=x) for x in vals]
                rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(samples)


def fit_decay(df: pd.DataFrame, col: str) -> dict:
    """Slope and standard error of ln(mean statistic) versus chain length."""
    g = df.groupby("n_sites")[col].mean()
    x, y = g.index.to_numpy(float), np.log(g.to_numpy())
    if len(x) < 4:
        return {"rate": float(-np.polyfit(x, y, 1)[0]), "stderr": None}
    coef, cov = np.polyfit(x, y, 1, cov=True)
    return {"rate": float(-coef[0]), "stderr": float(np.sqrt(cov[0, 0]))}


def _ms(x: pd.Series) -> dict:
    return {"mean": float(x.mean()), "std": float(x.std(ddof=1))}


def summarise(cfg: dict, sweep: pd.DataFrame, scaling: pd.DataFrame, ablation: pd.DataFrame,
              shots: pd.DataFrame) -> dict:
    """Headline numbers used in the report (all derived from the CSV files)."""
    n = cfg["main_size"]
    main = sweep[sweep.n_sites == n]
    integ, chaos = main[main.hz == 0.0], main[main.hz == HZ_CHAOTIC]
    out: dict = {"main_size": n, "lambda": LAM, "device": str(device()),
                 "gap_ratio": {"integrable": _ms(integ.gap_ratio), "chaotic": _ms(chaos.gap_ratio)}}
    out["diagnostics"] = {}
    for col in DIAGNOSTICS:
        rho, p = spearmanr(main[col], main["gap_ratio"])
        out["diagnostics"][col] = dict(
            spearman_vs_gap_ratio=float(rho), spearman_p=float(p),
            integrable=_ms(integ[col]), chaotic=_ms(chaos[col]),
            contrast=float(integ[col].mean() / chaos[col].mean()),
            auc_exact=auc_score(integ[col].to_numpy(), chaos[col].to_numpy()),
        )
    rho, p = spearmanr(main["backflow_x_haar"], main["energy_width"])
    out["spearman_backflow_x_haar_vs_energy_width"] = {"rho": float(rho), "p": float(p)}
    out["scaling_decay_per_site"] = {}
    for state in STATES:
        for name, hz in (("integrable", 0.0), ("chaotic", HZ_CHAOTIC)):
            d = scaling[(scaling.state == state) & (scaling.hz == hz)]
            for col in ("backflow", "plateau"):
                out["scaling_decay_per_site"][f"{col}_{state}_{name}"] = fit_decay(d, col)
    out["haar_floor_rates"] = {"amplitude_ln2_over_2": float(np.log(2) / 2),
                               "squared_ln2": float(np.log(2))}
    big = scaling[scaling.n_sites == scaling.n_sites.max()]
    out["largest_size"] = {"n_sites": int(big.n_sites.max())}
    for state in STATES:
        for name, hz in (("integrable", 0.0), ("chaotic", HZ_CHAOTIC)):
            d = big[(big.state == state) & (big.hz == hz)]
            out["largest_size"][f"{state}_{name}"] = {"plateau": _ms(d.plateau),
                                                      "backflow": _ms(d.backflow)}
        gi = big[(big.state == state) & (big.hz == 0.0)]
        gc = big[(big.state == state) & (big.hz == HZ_CHAOTIC)]
        out["largest_size"][f"{state}_plateau_contrast"] = float(gi.plateau.mean()
                                                                 / gc.plateau.mean())
    out["ablation"] = []
    for (kind, lam, state), sub in ablation.groupby(["coupling", "lam", "state"]):
        gi, gc = sub[sub.hz == 0.0], sub[sub.hz == HZ_CHAOTIC]
        out["ablation"].append(dict(
            coupling=kind, lam=float(lam), state=state,
            backflow_integrable=_ms(gi.backflow), backflow_chaotic=_ms(gc.backflow),
            plateau_integrable=_ms(gi.plateau), plateau_chaotic=_ms(gc.plateau),
            backflow_contrast=float(gi.backflow.mean() / gc.backflow.mean()),
            plateau_contrast=float(gi.plateau.mean() / gc.plateau.mean()),
            backflow_auc=auc_score(gi.backflow.to_numpy(), gc.backflow.to_numpy()),
            plateau_auc=auc_score(gi.plateau.to_numpy(), gc.plateau.to_numpy()),
        ))
    out["shot_budget"] = {}
    for state in STATES:
        d = shots[shots.state == state]
        entry = {}
        for est in ESTIMATORS:
            ok = d[(d.gamma == 0.0) & (d[f"auc_{est}"] >= 0.95)]
            best = d[d.shots == d.shots.max()]
            tol = best[best[f"auc_{est}"] >= 0.95]
            entry[est] = dict(
                min_shots_auc95_gamma0=int(ok.shots.min()) if len(ok) else None,
                max_gamma_auc95_at_max_shots=float(tol.gamma.max()) if len(tol) else None,
                auc_gamma0_max_shots=float(best[best.gamma == 0.0][f"auc_{est}"].iloc[0]),
            )
        out["shot_budget"][state] = entry
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="smoke run (< 1 minute)")
    args = parser.parse_args()
    cfg = config(args.quick)
    res = ROOT / "results" / ("quick" if args.quick else "")
    fig = ROOT / "figures" / ("quick" if args.quick else "")
    res.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    times = np.arange(0.0, cfg["t_max"] + 1e-9, cfg["dt"])
    t0 = time.time()

    sweep, traces = run_sweep(cfg, times)
    print(f"[sweep]    {time.time() - t0:6.1f} s", flush=True)
    scaling = run_scaling(cfg, times)
    print(f"[scaling]  {time.time() - t0:6.1f} s", flush=True)
    ablation = run_ablation(cfg, times)
    print(f"[ablation] {time.time() - t0:6.1f} s", flush=True)
    shots, samples = run_shots(cfg, times, traces)
    print(f"[shots]    {time.time() - t0:6.1f} s", flush=True)

    sweep.to_csv(res / "sweep.csv", index=False)
    scaling.to_csv(res / "scaling.csv", index=False)
    ablation.to_csv(res / "ablation.csv", index=False)
    shots.to_csv(res / "shots.csv", index=False)
    samples.to_csv(res / "shot_samples.csv", index=False)
    n = cfg["main_size"]
    tr = pd.DataFrame({"time": times})
    for state in STATES:
        for name, hz in (("integrable", 0.0), ("chaotic", HZ_CHAOTIC)):
            tr[f"abs_c_{state}_{name}"] = np.abs(traces[(state, n, hz, 0)])
    tr.to_csv(res / "example_traces.csv", index=False)

    summary = summarise(cfg, sweep, scaling, ablation, shots)
    summary["runtime_seconds"] = round(time.time() - t0, 1)
    summary["config"] = {k: v for k, v in cfg.items() if k != "scaling_seeds"}
    (res / "summary.json").write_text(json.dumps(summary, indent=2))

    plotting.fig_sweep_bars(sweep, fig / "fig1_sweep_bars.png")
    plotting.fig_ecdf(samples, cfg["ecdf_shots"], fig / "fig2_shot_ecdf.png")
    plotting.fig_heatmaps(shots, fig / "fig3_auc_heatmap.png")
    plotting.fig_scaling(scaling, tr, fig / "fig4_scaling.png")
    for col, d in summary["diagnostics"].items():
        print(f"{col:20s} spearman {d['spearman_vs_gap_ratio']:+.2f}  "
              f"contrast {d['contrast']:7.2f}  auc {d['auc_exact']:.2f}")
    print(f"done in {time.time() - t0:.1f} s on {device()}")


if __name__ == "__main__":
    main()
