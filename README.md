# Physics of Complex Networks — Course Project

Course project for *Physics of Complex Networks: Structure and Dynamics*
(SCQ2101383), University of Padova, a.y. 2025/26. Author: Giancarlo Venturato.

Two studies are included:

| Type | Title |
|------|-------|
| Theoretical | Epidemic spreading on temporal networks: SI, SIS, SIR, SEIR (mean-field) |
| Data | Monthly snapshots of the Internet AS-level network from BGP data |

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
cd code/temporal_epidemics && python run_experiments.py     # ~10 s, seeded
```
Writes figures to `code/temporal_epidemics/figures/` and edge lists / results to
`data/temporal_epidemics/`. See `code/temporal_epidemics/README.md` for details.

## AS-level Internet from BGP

In progress.

## Building the report

```bash
cd latex && latexmk -pdf main.tex
```
