# Physics of Complex Networks — Course Project

Course project for *Physics of Complex Networks: Structure and Dynamics*
(SCQ2101383), University of Padova, a.y. 2025/26. Author: Giancarlo Venturato.

Two projects are included (numbers, names and scores as listed on the course
Moodle):

| # | Type | Title | Score |
|----|------|-------|-------|
| 22 | Theoretical | Epidemic spreading on temporal networks: SI, SIS, SIR, SEIR (mean-field) | 0.5 |
| 40 | Data | Monthly snapshots of the Internet AS-level network from BGP data | 1.2 |

## Repository layout

```
report.pdf        compiled report
latex/            LaTeX source and figures
code/temporal_epidemics/     temporal-network epidemic simulator (Python + igraph)
code/internet_as_bgp/     BGP AS-level network pipeline
data/temporal_epidemics/     output edge lists (node_from,node_to,weight) and results
data/internet_as_bgp/     output AS-level snapshot networks
```

## Epidemic spreading on temporal networks

A single parametrised SI/SIS/SIR/SEIR stochastic engine run on synthetic temporal
topologies (activity-driven, and edge-activated Barabási–Albert / Erdős–Rényi
backbones) and on an empirical SocioPatterns contact network, each compared against
the time-aggregated static and time-shuffled null models and against the mean-field
epidemic thresholds.

```bash
pip install -r requirements.txt
cd code/temporal_epidemics && python run_experiments.py     # ~15 s, seeded
```
Writes figures to `code/temporal_epidemics/figures/` and edge lists / results to
`data/temporal_epidemics/`. See `code/temporal_epidemics/README.md` for details.

## AS-level Internet from BGP

A reproducible monthly time series of the Internet AS-level topology, built from
the RIPE RIS `rrc00` routing-table archive: one `bview` dump per month from
November 1999 to April 2026 (314 usable months), parsed and cleaned into weighted
edge lists (`asn1 asn2 weight`, weight = AS-pair observation frequency), plus the
exact-file manifest and per-month cleaning statistics. Analyses cover growth,
degree distributions (power-law MLE), node/edge/hub turnover, degree correlations
and Leiden community structure over 26.5 years.

```bash
brew install bgpdump
cd code/internet_as_bgp
python build_edges.py --start 1999-11 --end 2026-04   # ~4 h, ~35 GB transferred
python analyze_networks.py                            # ~30 min
python make_figures.py
```
Monthly edge lists and analysis tables are committed under `data/internet_as_bgp/`, so
the analysis and figures can be reproduced without re-downloading the BGP dumps.
See `code/internet_as_bgp/README.md` for details.

## Building the report

```bash
cd latex && latexmk -pdf main.tex
```
