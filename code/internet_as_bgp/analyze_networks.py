#!/usr/bin/env python3
"""Compute monthly network metrics for the AS-level Internet time series.

PoCN data project. Reads every `as_edges_YYYY_MM_weighted.tsv.gz` produced by
`build_edges.py` (in chronological order) and writes

  data/internet_as_bgp/monthly_metrics.tsv   one row per month:
      structure   : nodes, edges, mean degree, largest-connected-component share
      heavy tail  : discrete power-law MLE (alpha, xmin, sigma) on the degree
                    sequence + log-likelihood-ratio test vs lognormal (R, p)
      correlations: exponent mu of k_nn(k) ~ k^mu (log-log OLS on the
                    degree-averaged neighbour degree)
      communities : Leiden modularity partition (Q, number of communities)
      turnover    : node/edge births, deaths and Jaccard overlap w.r.t. the
                    previous available snapshot, top-20 hub set overlap
  data/internet_as_bgp/hubs_top20.tsv        per month, the 20 highest-degree ASNs

All metrics are computed on the *unweighted* simple graph (the weighted file
is reduced by ignoring the weight column); weights only enter the edge lists.

Usage:  python analyze_networks.py            # all months found
        python analyze_networks.py --limit 24 # first 24 months (smoke test)
"""
import argparse
import glob
import gzip
import os
import re
import warnings

import igraph as ig
import leidenalg as la
import numpy as np

DATA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "data", "internet_as_bgp"))

METRIC_COLS = ("month nodes edges mean_k lcc_frac alpha xmin sigma R_ln p_ln "
               "mu_knn modularity n_comm prev_month new_nodes lost_nodes "
               "node_jaccard new_edges lost_edges edge_jaccard hub20_overlap").split()


def read_graph(path):
    """Edge file -> (igraph.Graph, node set, edge set of (a, b) tuples)."""
    edges = []
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            a, b, _ = line.split("\t")
            edges.append((int(a), int(b)))
    g = ig.Graph.TupleList(edges, directed=False)
    return g, set(g.vs["name"]), set(edges)


def powerlaw_fit(degrees):
    """Discrete power-law MLE + comparison with lognormal. Returns 5 floats."""
    import powerlaw
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = powerlaw.Fit(degrees, discrete=True, verbose=False)
        r, p = fit.distribution_compare("power_law", "lognormal",
                                        normalized_ratio=True)
    return fit.power_law.alpha, fit.power_law.xmin, fit.power_law.sigma, r, p


def knn_exponent(g):
    """OLS slope of log <k_nn>(k) vs log k (igraph degree-averaged knn)."""
    knnk = np.asarray(g.knn()[1], dtype=float)          # index i -> degree i+1
    k = np.arange(1, len(knnk) + 1, dtype=float)
    ok = np.isfinite(knnk) & (knnk > 0) & (k > 1)       # k=1 dominated by leaves
    if ok.sum() < 5:
        return float("nan")
    return float(np.polyfit(np.log(k[ok]), np.log(knnk[ok]), 1)[0])


def leiden(g):
    part = la.find_partition(g, la.ModularityVertexPartition, seed=42,
                             n_iterations=2)
    return part.modularity, len(part)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--data-dir", default=DATA)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.data_dir, "as_edges_*_weighted.tsv.gz")))
    if args.limit:
        files = files[:args.limit]
    if not files:
        raise SystemExit(f"no edge files in {args.data_dir}")

    metrics_path = os.path.join(args.data_dir, "monthly_metrics.tsv")
    hubs_path = os.path.join(args.data_dir, "hubs_top20.tsv")
    fm = open(metrics_path, "w")
    fm.write("\t".join(METRIC_COLS) + "\n")
    fh = open(hubs_path, "w")
    fh.write("month\ttop20_asns_by_degree\n")

    prev = None                                          # (month, nodes, edges, hubs)
    for path in files:
        month = "-".join(re.search(r"as_edges_(\d{4})_(\d{2})", path).groups())
        g, nodes, edges = read_graph(path)
        n, m = g.vcount(), g.ecount()
        deg = np.asarray(g.degree())
        lcc_frac = g.connected_components().giant().vcount() / n
        alpha, xmin, sigma, r_ln, p_ln = powerlaw_fit(deg)
        mu = knn_exponent(g)
        q, n_comm = leiden(g)

        order = np.argsort(deg)[::-1][:20]
        hubs = [g.vs[int(i)]["name"] for i in order]
        fh.write(f"{month}\t{','.join(str(a) for a in hubs)}\n")

        if prev is None:
            turn = [""] + [float("nan")] * 7
        else:
            pm, pn, pe, phubs = prev
            turn = [pm,
                    len(nodes - pn), len(pn - nodes),
                    len(nodes & pn) / len(nodes | pn),
                    len(edges - pe), len(pe - edges),
                    len(edges & pe) / len(edges | pe),
                    len(set(hubs) & set(phubs)) / 20]
        row = [month, n, m, 2 * m / n, lcc_frac, alpha, xmin, sigma, r_ln, p_ln,
               mu, q, n_comm] + turn
        fm.write("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x)
                           for x in row) + "\n")
        fm.flush()
        print(f"{month}  N={n}  E={m}  lcc={lcc_frac:.3f}  alpha={alpha:.2f}  "
              f"mu={mu:.2f}  Q={q:.3f}")
        prev = (month, nodes, edges, hubs)

    fm.close()
    fh.close()
    print(f"\n-> {metrics_path}\n-> {hubs_path}")


if __name__ == "__main__":
    main()
