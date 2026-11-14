"""Figures of the study (Okabe-Ito palette, 200 dpi)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OI = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000"}
DPI = 200


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def fig_sweep_bars(sweep: pd.DataFrame, path: Path) -> None:
    """Grouped, hatched bars: chaos diagnostics along the integrability sweep.

    Each diagnostic is divided by its own mean at hz = 0 so that the bars show
    how much contrast it has between the integrable and the chaotic regime.
    """
    n = int(sweep.n_sites.max())
    d = sweep[sweep.n_sites == n]
    g = d.groupby("hz")
    hz = np.array(sorted(d.hz.unique()))
    x = np.arange(len(hz))
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 6.6), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 2.4]})
    r = g.gap_ratio
    ax0.errorbar(x, r.mean(), yerr=r.std(ddof=1), color=OI["black"], marker="o", capsize=3)
    ax0.axhline(0.3863, color=OI["blue"], ls="--", lw=1)
    ax0.axhline(0.5307, color=OI["vermillion"], ls="--", lw=1)
    ax0.text(len(hz) - 0.55, 0.395, "Poisson 0.386", color=OI["blue"], ha="right", fontsize=8)
    ax0.text(-0.4, 0.542, "random matrix 0.531", color=OI["vermillion"], fontsize=8)
    ax0.set_ylim(0.34, 0.58)
    ax0.set_ylabel("mean gap ratio")
    ax0.set_title(f"Integrable-to-chaotic sweep, {n}-spin environment "
                  f"(mean ± std, {d.seed.nunique()} seeds)")
    specs = [("backflow_x_haar", "backflow, random state", OI["sky"], "//"),
             ("backflow_x_product", "backflow, product state", OI["blue"], "\\\\"),
             ("plateau_x_product", "late-time plateau, product state", OI["orange"], ".."),
             ("otoc_late", "late-time correlator (baseline)", OI["green"], "xx")]
    width = 0.2
    for k, (col, label, color, hatch) in enumerate(specs):
        ref = g[col].mean().iloc[0]
        mean, std = g[col].mean() / ref, g[col].std(ddof=1) / ref
        lower = np.minimum(std, mean * 0.999)
        ax1.bar(x + (k - 1.5) * width, mean, width, yerr=[lower, std], color=color, hatch=hatch,
                edgecolor="black", linewidth=0.6, capsize=1.5, label=label,
                error_kw={"elinewidth": 0.8})
    ax1.axhline(1.0, color="black", lw=0.6, ls=":")
    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{h:g}" for h in hz])
    ax1.set_xlabel("longitudinal field hz (units of J)")
    ax1.set_ylabel("diagnostic / its mean at hz = 0")
    ax1.legend(frameon=False, fontsize=8, ncol=4, loc="upper center",
               bbox_to_anchor=(0.5, -0.16))
    _save(fig, path)


def fig_ecdf(samples: pd.DataFrame, shots: int, path: Path) -> None:
    """Empirical CDFs of the finite-shot statistics, integrable vs chaotic."""
    names = [("naive", "naive backflow"), ("thresholded", "thresholded backflow"),
             ("plateau", "late-time plateau")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.7))
    for ax, (est, title) in zip(axes, names):
        d = samples[samples.estimator == est]
        for flag, color, ls, label in ((True, OI["blue"], "-", "integrable (hz = 0)"),
                                       (False, OI["vermillion"], "--", "chaotic (hz = 0.809)")):
            v = np.sort(d[(d.hz == 0.0) == flag].value.to_numpy())
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=color, ls=ls,
                    lw=1.8, label=label)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("measured statistic (dimensionless)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("empirical cumulative distribution")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle(f"Product-state probe at {shots} shots per quadrature and time point, "
                 "no extra dephasing", fontsize=11)
    _save(fig, path)


def fig_heatmaps(shots: pd.DataFrame, path: Path) -> None:
    """Annotated heatmaps of classification AUC over shots and probe dephasing."""
    panels = [("product", "auc_plateau", "product state, late-time plateau"),
              ("product", "auc_thresholded", "product state, thresholded backflow"),
              ("haar", "auc_plateau", "random state, late-time plateau")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    for ax, (state, col, title) in zip(axes, panels):
        piv = shots[shots.state == state].pivot(index="gamma", columns="shots",
                                                values=col).sort_index()
        im = ax.imshow(piv.to_numpy(), origin="lower", cmap="cividis", vmin=0.0, vmax=1.0,
                       aspect="auto")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(piv.shape[0]))
        ax.set_yticklabels([f"{g:g}" for g in piv.index])
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                val = piv.iloc[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8,
                        color="black" if val > 0.6 else "white")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("shots per quadrature and time point")
    axes[0].set_ylabel("extra probe dephasing rate (units of J)")
    fig.colorbar(im, ax=axes, label="AUC, integrable scored above chaotic", shrink=0.9)
    fig.suptitle("Separating integrable from chaotic environments with a noisy probe "
                 "(0.5 = chance)")
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_scaling(scaling: pd.DataFrame, traces: pd.DataFrame, path: Path) -> None:
    """Example coherence traces and system-size scaling of the plateau."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.1))
    lines = [("abs_c_product_integrable", OI["blue"], "-", "product state, integrable"),
             ("abs_c_product_chaotic", OI["vermillion"], "-", "product state, chaotic"),
             ("abs_c_haar_integrable", OI["sky"], ":", "random state, integrable"),
             ("abs_c_haar_chaotic", OI["orange"], ":", "random state, chaotic")]
    for col, color, ls, label in lines:
        ax0.plot(traces.time, traces[col], color=color, ls=ls, lw=1.2, label=label)
    ax0.set_yscale("log")
    ax0.set_ylim(1e-3, 1.5)
    ax0.set_xlabel("time (units of 1/J)")
    ax0.set_ylabel("probe coherence magnitude")
    ax0.set_title("Probe coherence, one disorder realisation")
    ax0.legend(frameon=False, fontsize=8, ncol=2, loc="lower center")
    style = [("product", True, OI["blue"], "o", "-", "product state, integrable"),
             ("product", False, OI["vermillion"], "s", "-", "product state, chaotic"),
             ("haar", True, OI["sky"], "o", ":", "random state, integrable"),
             ("haar", False, OI["orange"], "s", ":", "random state, chaotic")]
    for state, integ, color, marker, ls, label in style:
        d = scaling[(scaling.state == state) & ((scaling.hz == 0.0) == integ)]
        g = d.groupby("n_sites").plateau
        mean, std = g.mean(), g.std(ddof=1)
        ax1.errorbar(mean.index, mean, yerr=[np.minimum(std, 0.999 * mean), std], color=color,
                     marker=marker, ls=ls, capsize=3, label=label)
    sizes = np.array(sorted(scaling.n_sites.unique()), dtype=float)
    ax1.plot(sizes, 2.0**-sizes, color="black", lw=0.8, ls="--", label="1 / Hilbert dimension")
    ax1.set_yscale("log")
    ax1.set_xlabel("environment size L (spins)")
    ax1.set_ylabel("late-time plateau of squared coherence")
    ax1.set_title("System-size scaling (mean ± std over seeds)")
    ax1.legend(frameon=False, fontsize=8)
    _save(fig, path)
