# Epidemic spreading on temporal networks

One parametrised SI / SIS / SIR / SEIR simulator run on synthetic **temporal**
networks, compared against the time-aggregated (static) and time-shuffled
surrogates, and against the mean-field epidemic thresholds.

## Files
| file | role |
|------|------|
| `temporal_networks.py` | generators (activity-driven; edge-activated BA/ER backbone), an empirical SocioPatterns loader, aggregation, static & time-shuffled null models. Uses **python-igraph**. |
| `epidemic.py` | single discrete-time stochastic engine; model = SI/SIS/SIR/SEIR via a small transition table. |
| `theory.py` | analytic thresholds: activity-driven `1/[m(<a>+√<a²>)]` vs static HMF `<k>/<k²>` and QMF `1/λ_max`. |
| `run_experiments.py` | driver: builds networks, runs all experiments, writes figures + data. |

## Run
```bash
pip install python-igraph numpy matplotlib     # scipy/pandas already present
python run_experiments.py                       # ~10 s, fully seeded
```
Outputs: figures in `figures/`, edge lists + trajectory + thresholds in
`../../data/temporal_epidemics/`.

## What the figures show
- **fig1** — SI saturates, SIS reaches an endemic plateau, SIR/SEIR are single
  outbreaks (SEIR delayed by the latent E stage), all on one activity-driven net.
- **fig2** — the simulated SIS threshold on the temporal net coincides with the
  activity-driven theory `β_c`, and lies **far above** the static-network HMF
  threshold: treating the timeline as static wrongly predicts an epidemic where
  the temporal dynamics have none. Time-shuffled ≈ temporal (ordering is
  irrelevant in a memoryless activity-driven model — a control).
- **fig3** — aggregating the timeline (all contacts simultaneous) turns a slow,
  bounded temporal outbreak into an explosive one (peak ≈0.98 vs ≈0.18).
- **fig4** — backbone structure at fixed ⟨k⟩=4: the scale-free BA network spreads
  more than the homogeneous ER one (peak 0.12 vs 0.06).
- **fig5** — the same analysis on the empirical SocioPatterns Hypertext-2009
  contact network; here time-shuffling *reduces* the outbreak (final 0.11 vs 0.20),
  exposing temporal correlations absent from the memoryless synthetic models.
  Data auto-loaded from `data/temporal_epidemics/empirical/` (see its `SOURCE.txt`).

## Key numbers (default seed)
`N=1500, T=800, m=1, γ=2.1, μ=0.02`: β_c temporal-AD ≈ 0.158, β_c static-HMF ≈
1.9×10⁻⁴ — a ~800× gap that is the quantitative core of the report.
