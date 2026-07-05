# Monthly AS-level Internet networks from BGP

Reconstructs a monthly time series of the Internet at the level of Autonomous
Systems (ASes) from the routing tables archived by the RIPE RIS collector
`rrc00` (<https://data.ris.ripe.net/rrc00/>), November 1999 → April 2026, and
analyzes its structural evolution.

## Scripts

| script | what it does |
|---|---|
| `fetch_bview.py` | lists the archive and picks one `bview` table dump per month (prefers day 01 @ 00:00); `--dry-run` shows the selection without downloading |
| `build_edges.py` | full pipeline: download → parse → clean AS_PATHs → aggregate → write one weighted edge list per month, then delete the raw dump |
| `analyze_networks.py` | per-month metrics (size, LCC, power-law MLE, \(k_{nn}\) exponent, Leiden modularity, node/edge/hub turnover) → `monthly_metrics.tsv` |
| `make_figures.py` | report figures from the metrics table and selected snapshots |

## Pipeline (cleaning rules)

For every RIB entry of the monthly table dump the AS_PATH is processed as
follows: consecutive duplicate ASNs are collapsed (BGP prepending); malformed
or empty paths are dropped; paths containing AS_SETs (`{...}`) are dropped
entirely (a documented simplification — they are a negligible fraction, see
`cleaning_stats.tsv`); paths with fewer than two ASNs are dropped. Every
consecutive ASN pair of a retained path becomes an undirected edge stored in
canonical `min,max` order; self-loops are discarded. The integer edge weight
counts in how many RIB entries the pair was observed that month; ignoring the
weight column gives the unweighted graph. ASNs are the integer AS numbers from
the dumps, so node labels are consistent across snapshots.

Three archive formats are handled transparently: modern MRT `TABLE_DUMP`/
`TABLE_DUMP_V2` (parsed with [`bgpdump`](https://github.com/RIPE-NCC/bgpdump)),
the pre-2001 Zebra MRT dialect (single record packing many RIB entries, 16-bit
ASNs — built-in binary parser), and the gzipped GateD *text* tables that appear
in a few months of 2000 (built-in text parser). Truncated/corrupt archive
files are detected (a real table dump always yields ≥100k usable paths) and
the pipeline falls back to the remaining files of the month, largest first.

## Outputs (in `data/internet_as_bgp/`)

- `as_edges_YYYY_MM_weighted.tsv.gz` — columns `asn1  asn2  weight`
- `bview_manifest.tsv` — the exact archive file used (or failure status) per month
- `cleaning_stats.tsv` — per month: total/malformed/AS_SET/too-short/retained
  path counts and resulting graph size
- `monthly_metrics.tsv`, `hubs_top20.tsv` — analysis results

## Requirements

`bgpdump` on PATH (`brew install bgpdump`), plus
`python-igraph`, `leidenalg`, `powerlaw`, `numpy`, `pandas`, `matplotlib`
(see repository `requirements.txt`).

## Usage

```bash
python build_edges.py --start 1999-11 --end 2026-04   # ~4 h, ~35 GB transferred
python analyze_networks.py                            # ~30 min
python make_figures.py
```

Both stages are restartable: months whose edge file already exists are skipped
(`--overwrite` forces a rebuild).
