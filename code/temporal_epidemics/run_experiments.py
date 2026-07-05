"""
Experiment driver: epidemic spreading on temporal networks.

Runs the full set of experiments used in the report and writes:

  figures/  fig1_models.png        SI/SIS/SIR/SEIR on the activity-driven network
            fig2_threshold.png     SIS threshold: temporal vs static vs shuffled
            fig3_temporal_vs_static.png   SIR prevalence, three protocols
            fig4_backbone.png      SIR on the edge-activated backbone topology
  data/temporal_epidemics/  activity_driven_aggregated_edges.csv   (node_from,node_to,weight)
                 backbone_aggregated_edges.csv
                 sir_trajectory_activity_driven.csv
                 thresholds.csv

All experiments are seeded, so re-running reproduces the figures.  Run from
anywhere:  ``python run_experiments.py``.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import temporal_networks as tn
import epidemic as epi
import theory as th

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "temporal_epidemics"))
os.makedirs(FIG, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

# ---- global parameters ---------------------------------------------------- #
N, T, M = 1500, 800, 1
GAMMA, EPS = 2.1, 1e-2
MU, SIGMA = 0.02, 0.10
N_RUNS = 15
SEED = 7

plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False})
COL = {"SI": "#444", "SIS": "#1b7837", "SIR": "#2166ac", "SEIR": "#b2182b",
       "temporal": "#2166ac", "concurrent": "#b2182b", "shuffled": "#f1a340"}


def banner(txt: str) -> None:
    print(f"\n{'='*66}\n{txt}\n{'='*66}")


# --------------------------------------------------------------------------- #
def main() -> None:
    banner("Building activity-driven temporal network")
    ad = tn.activity_driven(N=N, T=T, m=M, gamma=GAMMA, eps=EPS, seed=SEED)
    agg = ad.aggregate()
    a1, a2 = th.activity_moments(ad.activities)
    bc_ad = th.ad_threshold(ad.activities, M) * MU          # critical beta (temporal)
    bc_hmf = th.static_hmf_threshold(agg) * MU
    bc_qmf = th.static_qmf_threshold(agg) * MU
    print(f"N={N}  T={T}  m={M}  <a>={a1:.4f}  <a^2>={a2:.5f}")
    print(f"aggregated graph: {agg.vcount()} nodes, {agg.ecount()} edges, "
          f"<k>={np.mean(agg.degree()):.2f}, <k^2>={np.mean(np.array(agg.degree())**2):.1f}")
    print(f"contacts over horizon: {ad.n_contacts}")
    print(f"critical beta  (mu={MU}):  temporal-AD={bc_ad:.4f}  "
          f"static-HMF={bc_hmf:.5f}  static-QMF={bc_qmf:.5f}")

    tn.save_edge_list(agg, os.path.join(DATA, "activity_driven_aggregated_edges.csv"))

    # ---- Fig 1: the four models on the same network ----------------------- #
    banner("Fig 1: SI / SIS / SIR / SEIR on the activity-driven network")
    beta_hi = 4.0 * bc_ad
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for model in ("SI", "SIS", "SIR", "SEIR"):
        res = epi.average_runs(ad, n_runs=N_RUNS, model=model, beta=beta_hi,
                               mu=MU, sigma=SIGMA, n_seeds=5)
        t = np.arange(T + 1)
        ax.plot(t, res["mean_prevalence"], color=COL[model], label=model, lw=1.8)
        ax.fill_between(t, res["mean_prevalence"] - res["std_prevalence"],
                        res["mean_prevalence"] + res["std_prevalence"],
                        color=COL[model], alpha=0.15)
        print(f"  {model:4s}  peak prevalence={res['mean_prevalence'].max():.3f}  "
              f"cumulative incidence={res['cumulative_incidence']:.3f}")
    ax.set(xlabel="time step", ylabel="prevalence  I(t)/N",
           title=f"Compartmental models on activity-driven net (β={beta_hi:.3f}, μ={MU})")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_models.png")); plt.close(fig)

    # ---- Fig 2: SIS threshold on the temporal net vs analytic predictions - #
    # The temporal and time-shuffled sweeps are cheap (~<k> edges/step); the
    # fully aggregated static net is dense, so instead of sweeping it we mark its
    # analytic threshold (beta_c^HMF, far to the left) -- at any beta on this axis
    # the static network is already deeply endemic.
    banner("Fig 2: SIS endemic prevalence vs beta (threshold comparison)")
    shuffled = ad.time_shuffled(seed=SEED)
    betas = np.linspace(0.2 * bc_ad, 4.0 * bc_ad, 12)
    curves = {}
    for name, net in (("temporal", ad), ("shuffled", shuffled)):
        vals = [epi.average_runs(net, n_runs=N_RUNS, model="SIS", beta=b,
                                 mu=MU, n_seeds=10)["endemic_prevalence"] for b in betas]
        curves[name] = np.array(vals)
        print(f"  {name:11s} endemic prevalence at max beta = {vals[-1]:.3f}")

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for name in ("temporal", "shuffled"):
        ax.plot(betas, curves[name], "o-", color=COL[name], label=name, lw=1.6, ms=4)
    ax.axvline(bc_ad, ls="--", color=COL["temporal"], alpha=0.8,
               label=r"$\beta_c$ temporal (AD theory)")
    ax.axvline(bc_hmf, ls="--", color=COL["concurrent"], alpha=0.8,
               label=r"$\beta_c$ static (HMF theory)")
    ax.set(xlabel=r"transmission probability $\beta$", ylabel="endemic prevalence",
           title="SIS threshold: temporal simulation matches AD theory,\n"
                 "far above the static-network threshold")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_threshold.png")); plt.close(fig)

    # ---- Fig 3: SIR prevalence, three protocols at a single beta ---------- #
    banner("Fig 3: SIR temporal vs aggregated-static vs time-shuffled")
    aggregated = ad.static_sequence()
    beta_sir = 3.0 * bc_ad
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    traj_for_csv = None
    for name, net, nr in (("temporal", ad, N_RUNS),
                          ("concurrent", aggregated, 8),   # aggregated static (dense)
                          ("shuffled", shuffled, N_RUNS)):
        res = epi.average_runs(net, n_runs=nr, model="SIR", beta=beta_sir,
                               mu=MU, n_seeds=5)
        label = "aggregated static" if name == "concurrent" else name
        t = np.arange(T + 1)
        ax.plot(t, res["mean_prevalence"], color=COL[name], label=label, lw=1.8)
        print(f"  {label:18s} peak={res['mean_prevalence'].max():.3f}  "
              f"final size={res['cumulative_incidence']:.3f}")
        if name == "temporal":
            traj_for_csv = res["mean_prevalence"]
    ax.set(xlabel="time step", ylabel="prevalence  I(t)/N",
           title=f"SIR: aggregating the timeline inflates the outbreak "
                 f"(β={beta_sir:.3f}, μ={MU})")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_temporal_vs_static.png")); plt.close(fig)

    np.savetxt(os.path.join(DATA, "sir_trajectory_activity_driven.csv"),
               np.column_stack([np.arange(T + 1), traj_for_csv]),
               delimiter=",", header="step,prevalence", comments="")

    # ---- Fig 4: second topology (edge-activated backbone) ----------------- #
    banner("Fig 4: SIR on an edge-activated Barabasi-Albert backbone")
    ba = tn.edge_activated(N=N, T=T, backbone="ba", m_ba=2, p_active=0.03, seed=SEED)
    ba_agg = ba.aggregate()
    tn.save_edge_list(ba_agg, os.path.join(DATA, "backbone_aggregated_edges.csv"))
    ba_agg_seq = ba.static_sequence()
    print(f"backbone: {ba_agg.ecount()} edges, <k>={np.mean(ba_agg.degree()):.2f}")
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for name, net, nr in (("temporal", ba, N_RUNS), ("concurrent", ba_agg_seq, 8)):
        res = epi.average_runs(net, n_runs=nr, model="SIR", beta=0.35,
                               mu=MU, n_seeds=5)
        label = "aggregated static" if name == "concurrent" else name
        t = np.arange(T + 1)
        ax.plot(t, res["mean_prevalence"], color=COL[name], label=label, lw=1.8)
        print(f"  {label:18s} peak={res['mean_prevalence'].max():.3f}  "
              f"final size={res['cumulative_incidence']:.3f}")
    ax.set(xlabel="time step", ylabel="prevalence  I(t)/N",
           title="SIR on edge-activated BA backbone (β=0.35, μ=%.2f)" % MU)
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4_backbone.png")); plt.close(fig)

    # ---- thresholds table ------------------------------------------------- #
    with open(os.path.join(DATA, "thresholds.csv"), "w") as fh:
        fh.write("quantity,value\n")
        fh.write(f"mu,{MU}\n<a>,{a1}\n<a^2>,{a2}\n")
        fh.write(f"beta_c_temporal_AD,{bc_ad}\n")
        fh.write(f"beta_c_static_HMF,{bc_hmf}\n")
        fh.write(f"beta_c_static_QMF,{bc_qmf}\n")

    banner("Done -- figures in figures/, data in data/temporal_epidemics/")


if __name__ == "__main__":
    main()
