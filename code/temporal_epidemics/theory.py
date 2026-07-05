"""
Mean-field epidemic thresholds used to interpret the simulations.

We compare two predictions for the SIS/SIR critical ratio (beta/mu)_c:

1. Activity-driven temporal network (Perra et al., Sci. Rep. 2, 469 (2012); see
   also Pastor-Satorras et al., Rev. Mod. Phys. 87, 925 (2015), Sec. VIII):

       (beta/mu)_c = 1 / [ m ( <a> + sqrt(<a^2>) ) ]

   where a_i is the node activity and m the number of links fired per activation.
   The threshold is governed by the *activity* moments, not by the aggregated
   degree distribution.

2. Static (time-aggregated) network, heterogeneous mean-field / quenched
   mean-field approximation (RMP 2015):

       (beta/mu)_c^{HMF}  = <k> / <k^2>
       (beta/mu)_c^{QMF}  = 1 / lambda_max(A)

   evaluated on the *aggregated* graph, i.e. treating all of a node's contacts
   over the whole horizon as simultaneously present.

Because the aggregate accumulates many contacts, <k^2> (and lambda_max) are large,
so the static threshold is much *lower* than the temporal one: assuming
concurrency systematically overestimates spreading.  This is the quantitative
statement the report makes.

The simulator works in discrete time with per-step probabilities (beta, mu); for
small probabilities these map onto the continuous rates of the formulae above, so
the analytic (beta/mu)_c is compared against the simulated ratio (beta/mu).
"""

from __future__ import annotations

import numpy as np


def activity_moments(activities: np.ndarray) -> tuple[float, float]:
    """Return (<a>, <a^2>) of the activity distribution."""
    return float(activities.mean()), float((activities ** 2).mean())


def ad_threshold(activities: np.ndarray, m: int) -> float:
    """Activity-driven (beta/mu)_c = 1 / [ m (<a> + sqrt(<a^2>)) ]."""
    a1, a2 = activity_moments(activities)
    return 1.0 / (m * (a1 + np.sqrt(a2)))


def static_hmf_threshold(g) -> float:
    """Heterogeneous mean-field (beta/mu)_c = <k>/<k^2> on the aggregated graph."""
    k = np.asarray(g.degree(), dtype=float)
    k1, k2 = k.mean(), (k ** 2).mean()
    return float(k1 / k2) if k2 > 0 else np.inf


def static_qmf_threshold(g) -> float:
    """Quenched mean-field (beta/mu)_c = 1/lambda_max(A) on the aggregated graph."""
    # Leading adjacency eigenvalue via igraph (power iteration under the hood).
    # The aggregate can be disconnected; lambda_max is still the max over
    # components, so we silence igraph's connectivity warning.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lam = g.eigenvector_centrality(scale=False, weights=None,
                                       return_eigenvalue=True)[1]
    return float(1.0 / lam) if lam > 0 else np.inf
