"""
One parametrised compartmental epidemic engine over a temporal network.

Per the project note ("if the simulation code is adequately structured, simulating
SI vs SIS vs SIR vs SEIR is a simple variation, not 4x the coding effort") a single
discrete-time stochastic core loop covers all four models; the model is just a small
transition table.

Compartments are encoded as integers::

    0 = S (susceptible)   1 = E (exposed, latent, not infectious)
    2 = I (infectious)    3 = R (recovered / removed)

Dynamics at each step, on the active edges of that snapshot:
  * infection : every S--I contact transmits with probability ``beta``.  In SEIR
                the newly infected enter E; otherwise they enter I directly.
  * progression (SEIR only): each E becomes I with probability ``sigma``.
  * recovery  : each I leaves the infectious state with probability ``mu`` -->
                back to S (SIS) or to R (SIR/SEIR).  SI has no recovery.

Transitions are evaluated on the pre-step state, and recoveries/progressions are
committed before new infections, so a node cannot be infected and recover within the
same step.  All randomness flows through one seeded Generator for reproducibility.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

S, E, I, R = 0, 1, 2, 3

# model -> (has_exposed, recovers, recovered_state)
_MODELS = {
    "SI":   (False, False, None),
    "SIS":  (False, True,  S),
    "SIR":  (False, True,  R),
    "SEIR": (True,  True,  R),
}


@dataclass
class Trajectory:
    """Compartment time series (each array length T+1, including t=0)."""
    model: str
    S: np.ndarray
    E: np.ndarray
    I: np.ndarray
    R: np.ndarray

    @property
    def T(self) -> int:
        return len(self.I) - 1

    @property
    def prevalence(self) -> np.ndarray:
        """Fraction currently infectious, I(t)/N."""
        N = self.S[0] + self.E[0] + self.I[0] + self.R[0]
        return self.I / N

    @property
    def cumulative_incidence(self) -> float:
        """Final fraction ever infected (peak of I+R+E vs the initial susceptibles)."""
        N = self.S[0] + self.E[0] + self.I[0] + self.R[0]
        return float((self.E[-1] + self.I[-1] + self.R[-1]) / N)

    @property
    def endemic_prevalence(self, frac: float = 0.2) -> float:
        """Mean prevalence over the last ``frac`` of the run (SI/SIS steady state)."""
        k = max(1, int(len(self.I) * frac))
        N = self.S[0] + self.E[0] + self.I[0] + self.R[0]
        return float(self.I[-k:].mean() / N)


def simulate(net, model: str = "SIR", beta: float = 0.5, mu: float = 0.01,
             sigma: float = 0.1, n_seeds: int = 5, seed: int = 0,
             init_infected: np.ndarray | None = None) -> Trajectory:
    """
    Run one stochastic realisation of ``model`` on temporal network ``net``.

    Parameters
    ----------
    net           : TemporalNetwork (provides ``.N`` and ``.snapshots``).
    beta          : per-contact, per-step transmission probability.
    mu            : per-step recovery probability (ignored for SI).
    sigma         : per-step E->I progression probability (SEIR only).
    n_seeds       : number of initially infectious nodes (if ``init_infected`` is None).
    """
    if model not in _MODELS:
        raise ValueError(f"unknown model {model!r}; choose from {list(_MODELS)}")
    has_E, recovers, rec_state = _MODELS[model]

    N, Tn = net.N, net.T
    rng = np.random.default_rng(seed)
    state = np.zeros(N, dtype=np.int8)
    if init_infected is None:
        init_infected = rng.choice(N, size=min(n_seeds, N), replace=False)
    state[init_infected] = I

    counts = np.zeros((Tn + 1, 4), dtype=np.int64)
    counts[0] = np.bincount(state, minlength=4)

    for t in range(Tn):
        edges = net.snapshots[t]
        # ---- 1. find susceptibles reached by an infectious contact ---------- #
        if len(edges):
            u = state[edges[:, 0]]
            v = state[edges[:, 1]]
            m_uI = (u == I) & (v == S)
            m_vI = (v == I) & (u == S)
            targets = np.concatenate([edges[m_uI, 1], edges[m_vI, 0]])
            if targets.size:
                hit = targets[rng.random(targets.size) < beta]
                newly = np.unique(hit)
            else:
                newly = np.empty(0, dtype=np.int64)
        else:
            newly = np.empty(0, dtype=np.int64)

        # ---- 2. progression + recovery on the pre-step state --------------- #
        if has_E:
            exposed = np.nonzero(state == E)[0]
            if exposed.size:
                promote = exposed[rng.random(exposed.size) < sigma]
                state[promote] = I
        if recovers:
            infectious = np.nonzero(state == I)[0]
            if infectious.size:
                rec = infectious[rng.random(infectious.size) < mu]
                state[rec] = rec_state

        # ---- 3. commit new infections -------------------------------------- #
        state[newly] = E if has_E else I

        counts[t + 1] = np.bincount(state, minlength=4)

    return Trajectory(model, counts[:, 0], counts[:, 1], counts[:, 2], counts[:, 3])


def average_runs(net, n_runs: int = 20, base_seed: int = 100, **kw) -> dict:
    """Average several stochastic realisations; returns mean trajectories + scalars."""
    prev, cum, end = [], [], []
    traj0 = None
    for r in range(n_runs):
        tr = simulate(net, seed=base_seed + r, **kw)
        if traj0 is None:
            traj0 = tr
        prev.append(tr.prevalence)
        cum.append(tr.cumulative_incidence)
        end.append(tr.endemic_prevalence)
    return dict(
        mean_prevalence=np.mean(prev, axis=0),
        std_prevalence=np.std(prev, axis=0),
        cumulative_incidence=float(np.mean(cum)),
        endemic_prevalence=float(np.mean(end)),
        model=traj0.model,
    )
