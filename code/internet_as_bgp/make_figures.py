#!/usr/bin/env python3
"""Figures for the AS-network analysis from `monthly_metrics.tsv` + selected edge files.

Writes to figures/ :
  fig1_growth.png     N(t), E(t) | mean degree | LCC share
  fig2_degree.png     degree CCDFs (snapshots) | power-law alpha(t)
  fig3_turnover.png   node/edge birth-death rates | Jaccard persistence | hub overlap
  fig4_structure.png  k_nn(k) snapshots | mu(t) | Leiden Q(t) | #communities
  fig5_snapshots.png  network drawings: 1999 giant component | 2026 inner core

Run after `analyze_networks.py`:  python make_figures.py
"""
import glob
import gzip
import os
import re

import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "internet_as_bgp"))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# categorical palette (Okabe-Ito subset, CVD-validated); viridis for year ramps
BLUE, VERM, GREEN, PINK = "#0072B2", "#D55E00", "#009E73", "#CC79A7"

plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5,
    "legend.frameon": False, "figure.dpi": 200, "savefig.bbox": "tight",
})


def month_to_x(s):
    y, m = int(s[:4]), int(s[5:7])
    return y + (m - 0.5) / 12


def load_metrics():
    df = pd.read_csv(os.path.join(DATA, "monthly_metrics.tsv"), sep="\t")
    df["x"] = df["month"].map(month_to_x)
    # break lines across archive gaps: reindex on a full monthly grid
    full = pd.period_range(df["month"].iloc[0], df["month"].iloc[-1], freq="M")
    df = df.set_index(pd.PeriodIndex(df["month"], freq="M")).reindex(full)
    df["x"] = [p.year + (p.month - 0.5) / 12 for p in df.index]
    # months whose predecessor is not the previous calendar month -> NaN turnover
    prev_ok = df["prev_month"].notna() & (
        df["prev_month"] == df.index.shift(-1).astype(str))
    for c in ["new_nodes", "lost_nodes", "node_jaccard", "new_edges",
              "lost_edges", "edge_jaccard", "hub20_overlap"]:
        df.loc[~prev_ok, c] = np.nan
    return df


def degrees_of(month):
    path = os.path.join(DATA, f"as_edges_{month[:4]}_{month[5:7]}_weighted.tsv.gz")
    edges = []
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            a, b, _ = line.split("\t")
            edges.append((int(a), int(b)))
    return ig.Graph.TupleList(edges, directed=False)


def pick_snapshots(wanted):
    have = sorted(re.search(r"as_edges_(\d{4})_(\d{2})", f).group(1, 2)
                  for f in glob.glob(os.path.join(DATA, "as_edges_*.tsv.gz")))
    have = [f"{y}-{m}" for y, m in have]
    return sorted({min(have, key=lambda h: abs(month_to_x(h) - month_to_x(w)))
                   for w in wanted})


def year_axis(ax):
    ax.set_xlabel("year")
    ax.xaxis.set_major_locator(plt.MultipleLocator(5))


def fig1(df):
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.0))
    a, b, c = axes
    a.plot(df["x"], df["nodes"], color=BLUE, lw=1.6, label="nodes $N$")
    a.plot(df["x"], df["edges"], color=VERM, lw=1.6, label="edges $E$")
    a.set_yscale("log")
    a.set_ylabel("count")
    a.legend()
    b.plot(df["x"], df["mean_k"], color=GREEN, lw=1.6)
    b.set_ylabel(r"mean degree $\langle k\rangle$")
    c.plot(df["x"], 100 * df["lcc_frac"], color=BLUE, lw=1.6)
    c.set_ylabel("largest component (% of nodes)")
    c.set_ylim(99, 100.02)
    for ax in axes:
        year_axis(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_growth.png"))


def fig2(df):
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    a, b = axes
    snaps = pick_snapshots(["2000-01", "2005-01", "2010-01", "2015-01",
                            "2020-01", "2026-04"])
    cmap = plt.get_cmap("viridis")
    for i, month in enumerate(snaps):
        deg = np.sort(np.asarray(degrees_of(month).degree()))
        ccdf = 1.0 - np.arange(len(deg)) / len(deg)
        a.loglog(deg, ccdf, color=cmap(0.85 * i / max(len(snaps) - 1, 1)),
                 lw=1.4, label=month)
    a.set_xlabel("degree $k$")
    a.set_ylabel(r"CCDF  $P(K\geq k)$")
    a.legend(fontsize=7)
    b.plot(df["x"], df["alpha"], color=BLUE, lw=1.4)
    b.fill_between(df["x"], df["alpha"] - df["sigma"], df["alpha"] + df["sigma"],
                   color=BLUE, alpha=0.25, lw=0)
    b.set_ylabel(r"power-law exponent $\hat\alpha$ ($\pm\hat\sigma$)")
    year_axis(b)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_degree.png"))


def roll(s, w=12):
    return s.rolling(w, min_periods=w // 2, center=True).median()


def fig3(df):
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.0))
    a, b, c = axes
    birth = df["new_nodes"] / df["nodes"]
    death = df["lost_nodes"] / df["nodes"]
    a.plot(df["x"], birth, color=GREEN, lw=0.5, alpha=0.35)
    a.plot(df["x"], death, color=PINK, lw=0.5, alpha=0.35)
    a.plot(df["x"], roll(birth), color=GREEN, lw=1.8, label="AS births")
    a.plot(df["x"], roll(death), color=PINK, lw=1.8, label="AS deaths")
    a.set_ylabel("monthly node turnover (fraction of $N$)")
    a.legend()
    b.plot(df["x"], df["node_jaccard"], color=BLUE, lw=0.5, alpha=0.35)
    b.plot(df["x"], df["edge_jaccard"], color=VERM, lw=0.5, alpha=0.35)
    b.plot(df["x"], roll(df["node_jaccard"]), color=BLUE, lw=1.8, label="nodes")
    b.plot(df["x"], roll(df["edge_jaccard"]), color=VERM, lw=1.8, label="edges")
    b.set_ylabel("Jaccard overlap with previous month")
    b.legend()
    c.plot(df["x"], df["hub20_overlap"], color=BLUE, lw=0.5, alpha=0.35)
    c.plot(df["x"], roll(df["hub20_overlap"]), color=BLUE, lw=1.8)
    c.set_ylabel("top-20 hub overlap with previous month")
    c.set_ylim(0.5, 1.02)
    for ax in axes:
        year_axis(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_turnover.png"))


def fig4(df):
    fig, axes = plt.subplots(2, 2, figsize=(8, 6.0))
    (a, b), (c, d) = axes
    snaps = pick_snapshots(["2001-01", "2013-01", "2026-04"])
    cmap = plt.get_cmap("viridis")
    for i, month in enumerate(snaps):
        g = degrees_of(month)
        knnk = np.asarray(g.knn()[1], dtype=float)
        k = np.arange(1, len(knnk) + 1, dtype=float)
        ok = np.isfinite(knnk) & (knnk > 0)
        col = cmap(0.85 * i / max(len(snaps) - 1, 1))
        a.loglog(k[ok], knnk[ok], ".", ms=2.5, color=col, label=month)
        fit = ok & (k > 1)
        mu, c0 = np.polyfit(np.log(k[fit]), np.log(knnk[fit]), 1)
        kk = np.array([2, k[fit].max()])
        a.loglog(kk, np.exp(c0) * kk ** mu, "-", lw=1.2, color=col)
    a.set_xlabel("degree $k$")
    a.set_ylabel(r"$\langle k_{nn}\rangle(k)$")
    a.legend(fontsize=7)
    b.plot(df["x"], df["mu_knn"], color=BLUE, lw=1.4)
    b.set_ylabel(r"$k_{nn}$ exponent $\mu$")
    c.plot(df["x"], df["modularity"], color=GREEN, lw=1.4)
    c.set_ylabel("Leiden modularity $Q$")
    d.plot(df["x"], df["n_comm"], color=VERM, lw=1.4)
    d.set_ylabel("number of communities")
    for ax in (b, c, d):
        year_axis(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_structure.png"))


def fig5():
    """Draw the earliest snapshot and the dense core of the latest one."""
    import random
    from matplotlib.collections import LineCollection
    ig.set_random_number_generator(random.Random(42))
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.4))
    cmap = plt.get_cmap("viridis")
    for ax, month, mode in [(axes[0], "1999-11", "full"), (axes[1], "2026-04", "core")]:
        g = degrees_of(month)
        if mode == "full":
            g = g.connected_components().giant()
            title = f"{month}  (giant component:\n$N$={g.vcount():,}, $E$={g.ecount():,})"
        else:
            core = np.asarray(g.coreness())
            k = 15
            while (core >= k).sum() > 4000:      # keep the drawing readable
                k += 5
            g = g.induced_subgraph(np.flatnonzero(core >= k).tolist())
            title = (f"{month}  ($k$-core, $k\\geq${k}:\n"
                     f"$N$={g.vcount():,}, $E$={g.ecount():,})")
        xy = np.asarray(g.layout_fruchterman_reingold(niter=500).coords)
        deg = np.asarray(g.degree(), dtype=float)
        rel = np.log1p(deg) / np.log1p(deg.max())
        segs = xy[np.asarray(g.get_edgelist())]
        ax.add_collection(LineCollection(segs, colors="#888888", linewidths=0.12,
                                         alpha=0.18, zorder=1))
        order = np.argsort(deg)                  # draw hubs on top
        ax.scatter(xy[order, 0], xy[order, 1], s=1 + 28 * rel[order] ** 2,
                   c=cmap(0.9 * rel[order]), linewidths=0, zorder=2)
        ax.set_title(title, fontsize=8)
        ax.set_aspect("equal")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_snapshots.png"))


def main():
    df = load_metrics()
    fig1(df)
    fig2(df)
    fig3(df)
    fig4(df)
    fig5()
    print(f"figures -> {FIG}")


if __name__ == "__main__":
    main()
