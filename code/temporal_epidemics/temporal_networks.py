"""
Synthetic temporal-network generators for the epidemics-on-temporal-networks study.

A temporal network is represented as a sequence of ``T`` undirected snapshots.
Snapshot ``t`` is an integer array of shape ``(E_t, 2)`` listing the edges that are
*active during time step t only* (edges last a single step, as in the
activity-driven model of Perra et al., Sci. Rep. 2, 469 (2012)).

Two structurally different generators are provided so that the epidemic models are
studied "across a variety of different synthetic topologies" (project brief):

  * activity_driven   -- no fixed backbone; heterogeneous node activity a_i ~ x^{-gamma}.
                         The canonical synthetic temporal model.
  * edge_activated    -- a fixed static backbone (Barabasi-Albert or Erdos-Renyi)
                         whose edges switch on independently at each step.

Both return a :class:`TemporalNetwork`, which also builds the time-aggregated
(static) weighted graph and the two null models used for comparison in the report:
the concurrent static sequence and the time-shuffled surrogate.

Uses python-igraph for backbone generation and for the aggregated-graph object.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import igraph as ig
import numpy as np


# --------------------------------------------------------------------------- #
#  Container
# --------------------------------------------------------------------------- #
@dataclass
class TemporalNetwork:
    """A discrete-time temporal network as a list of per-step edge arrays."""

    N: int
    snapshots: List[np.ndarray]           # length T; each (E_t, 2) int array
    kind: str = "generic"
    activities: np.ndarray | None = None  # node activities, if activity-driven
    meta: dict = field(default_factory=dict)

    @property
    def T(self) -> int:
        return len(self.snapshots)

    @property
    def n_contacts(self) -> int:
        """Total number of (edge, time) contact events over the whole horizon."""
        return int(sum(len(s) for s in self.snapshots))

    # ---- derived representations ----------------------------------------- #
    def aggregate(self) -> ig.Graph:
        """Time-aggregated static weighted graph (weight = contact frequency)."""
        weight: dict[tuple[int, int], int] = {}
        for snap in self.snapshots:
            for a, b in snap:
                key = (a, b) if a < b else (b, a)
                weight[key] = weight.get(key, 0) + 1
        edges = list(weight.keys())
        g = ig.Graph(n=self.N, edges=edges, directed=False)
        g.es["weight"] = [weight[e] for e in edges]
        return g

    def static_sequence(self) -> "TemporalNetwork":
        """
        Time-aggregated *static* network as a temporal sequence: every edge of
        the aggregate is present at every step.  This is the standard "ignore
        timing, treat all contacts as simultaneous" approximation whose epidemic
        threshold is the static mean-field <k>/<k^2> (see theory.py).  Because a
        node's whole contact neighbourhood is permanently available, spreading is
        far easier than on the temporal network -- the contrast the report makes.
        """
        base = np.array(self.aggregate().get_edgelist(), dtype=np.int64)
        snaps = [base] * self.T          # same array reused; simulate only reads it
        return TemporalNetwork(self.N, snaps, kind=self.kind + "+aggregated",
                               meta={"null": "aggregated-static"})

    def time_shuffled(self, seed: int = 0) -> "TemporalNetwork":
        """
        Null model isolating *ordering*: the multiset of contact events is kept
        but each event is reassigned to a uniformly random time step.  Destroys
        temporal correlations / causal paths while preserving the aggregate graph
        and the total contact volume exactly.
        """
        rng = np.random.default_rng(seed)
        all_edges = np.concatenate(self.snapshots, axis=0) if self.n_contacts else \
            np.empty((0, 2), dtype=np.int64)
        assign = rng.integers(0, self.T, size=len(all_edges))
        snaps = [all_edges[assign == t] for t in range(self.T)]
        return TemporalNetwork(self.N, snaps, kind=self.kind + "+shuffled",
                               meta={"null": "time-shuffled"})


# --------------------------------------------------------------------------- #
#  Generators
# --------------------------------------------------------------------------- #
def _sample_activity(N: int, gamma: float, eps: float, rng) -> np.ndarray:
    """Draw N activities from pdf f(x) ~ x^{-gamma} on [eps, 1] (inverse CDF)."""
    u = rng.random(N)
    if abs(gamma - 1.0) < 1e-9:
        return eps * (1.0 / eps) ** u
    lo, exp = eps ** (1.0 - gamma), 1.0 - gamma
    return (u * (1.0 - lo) + lo) ** (1.0 / exp)


def activity_driven(N: int = 2000, T: int = 1000, m: int = 1,
                    gamma: float = 2.1, eps: float = 1e-2, eta: float = 1.0,
                    seed: int = 0) -> TemporalNetwork:
    """
    Activity-driven temporal network (Perra et al. 2012).

    Each node i has activity a_i = eta * x_i with x_i ~ f(x) = x^{-gamma},
    x in [eps, 1].  At every step, node i is active with probability a_i and, if
    active, opens m undirected links to uniformly-chosen partners (sampled with
    replacement; accidental self-pairs are dropped).
    """
    rng = np.random.default_rng(seed)
    a = np.clip(eta * _sample_activity(N, gamma, eps, rng), 0.0, 1.0)
    snapshots: List[np.ndarray] = []
    for _ in range(T):
        active = np.nonzero(rng.random(N) < a)[0]
        if active.size == 0:
            snapshots.append(np.empty((0, 2), dtype=np.int64))
            continue
        src = np.repeat(active, m)
        dst = rng.integers(0, N, size=src.size)
        good = src != dst                       # drop accidental self-loops
        snapshots.append(np.stack([src[good], dst[good]], axis=1))
    return TemporalNetwork(N, snapshots, kind="activity-driven", activities=a,
                           meta=dict(m=m, gamma=gamma, eps=eps, eta=eta))


def edge_activated(N: int = 2000, T: int = 1000, backbone: str = "ba",
                   m_ba: int = 2, avg_k: float = 4.0, p_active: float = 0.02,
                   seed: int = 0) -> TemporalNetwork:
    """
    Edge-activated temporal network on a fixed static backbone.

    A Barabasi-Albert (``backbone="ba"``) or Erdos-Renyi (``"er"``) graph is
    generated once; at each step every backbone edge is independently active with
    probability ``p_active``.  The time-aggregated graph is (a weighted copy of)
    the backbone, giving a genuinely different topology from the activity-driven
    case, where the aggregate has no fixed skeleton.
    """
    rng = np.random.default_rng(seed)
    import random as _random
    _random.seed(seed)
    ig.set_random_number_generator(_random)
    if backbone == "ba":
        g = ig.Graph.Barabasi(n=N, m=m_ba)
    elif backbone == "er":
        g = ig.Graph.Erdos_Renyi(n=N, m=int(avg_k * N / 2))
    else:
        raise ValueError(f"unknown backbone {backbone!r}")
    base = np.array(g.get_edgelist(), dtype=np.int64)
    snapshots = []
    for _ in range(T):
        keep = rng.random(len(base)) < p_active
        snapshots.append(base[keep])
    return TemporalNetwork(N, snapshots, kind=f"edge-activated({backbone})",
                           meta=dict(backbone=backbone, p_active=p_active,
                                     backbone_edges=len(base)))


# --------------------------------------------------------------------------- #
#  Empirical temporal network (SocioPatterns)
# --------------------------------------------------------------------------- #
def load_sociopatterns(path: str, bin_seconds: int = 300) -> TemporalNetwork:
    """
    Load a SocioPatterns-style timestamped contact list into a TemporalNetwork.

    The file has one contact per line, ``source, target, time`` (comment lines
    start with ``#``); timestamps are in seconds.  Contacts are binned into
    consecutive windows of ``bin_seconds`` (one snapshot per window; duplicate
    contacts within a window are collapsed), and node ids are remapped to a
    contiguous range.  Default here: the ACM Hypertext 2009 conference network
    (113 attendees, 20 s resolution, ~59 h), retrieved from the Netzschleuder
    repository (``sp_hypertext``); see Isella et al., J. Theor. Biol. 271 (2011).
    """
    raw = np.loadtxt(path, delimiter=",", comments="#", usecols=(0, 1, 2),
                     dtype=np.int64)
    src, dst, t = raw[:, 0], raw[:, 1], raw[:, 2]
    nodes = np.unique(np.concatenate([src, dst]))
    remap = {int(n): i for i, n in enumerate(nodes)}
    src = np.array([remap[int(x)] for x in src])
    dst = np.array([remap[int(x)] for x in dst])
    step = (t - t.min()) // bin_seconds
    T = int(step.max()) + 1
    snaps = []
    for s in range(T):
        m = step == s
        if not m.any():
            snaps.append(np.empty((0, 2), dtype=np.int64))
            continue
        e = np.unique(np.stack([src[m], dst[m]], axis=1), axis=0)  # dedupe window
        snaps.append(e)
    return TemporalNetwork(len(nodes), snaps, kind="empirical(sp_hypertext)",
                           meta=dict(bin_seconds=bin_seconds, source="sp_hypertext"))


# --------------------------------------------------------------------------- #
#  I/O helpers
# --------------------------------------------------------------------------- #
def save_edge_list(g: ig.Graph, path: str) -> None:
    """Write a weighted graph as ``node_from,node_to,weight`` (weight=1 if absent)."""
    w = g.es["weight"] if "weight" in g.es.attributes() else [1] * g.ecount()
    with open(path, "w") as fh:
        fh.write("node_from,node_to,weight\n")
        for (a, b), wt in zip(g.get_edgelist(), w):
            fh.write(f"{a},{b},{int(wt)}\n")
