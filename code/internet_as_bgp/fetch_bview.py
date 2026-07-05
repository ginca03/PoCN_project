#!/usr/bin/env python3
"""Download one RIPE RIS rrc00 `bview` MRT dump per month.

PoCN data project. Picks, for each month in [start, end], the first available
`bview.YYYYMMDD.HHMM.gz` (preferring day 01, time 0000), records the exact filename
used, and downloads it. Does NOT parse — that is the next stage (bgpdump / pybgpstream).

Usage:
    python fetch_bview.py --start 1999-11 --end 2026-04 --out ../../data/internet_as_bgp/raw
    python fetch_bview.py --start 2010-01 --end 2010-12 --dry-run      # list only, no download

This is a STUB to validate the download step on a small window before scaling to the
full ~318-month series. Verify behavior on a 1-year window first.
"""
import argparse, os, re, sys, urllib.request
from datetime import date

BASE = "https://data.ris.ripe.net/rrc00"
BVIEW_RE = re.compile(r'bview\.(\d{8})\.(\d{4})\.gz')


def months(start, end):
    y, m = start
    while (y, m) <= end:
        yield y, m
        m += 1
        if m > 12:
            y, m = y + 1, 1


def list_bviews(year, month, retries=3):
    url = f"{BASE}/{year:04d}.{month:02d}/"
    for attempt in range(1, retries + 1):
        try:
            html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
            return sorted(set(BVIEW_RE.findall(html)))  # list of (YYYYMMDD, HHMM)
        except Exception as e:
            print(f"  ! {year}-{month:02d}: index fetch failed ({e}), attempt {attempt}/{retries}",
                  file=sys.stderr)
            import time
            time.sleep(5 * attempt)
    return []


def pick(bviews):
    """Prefer day 01 @ 0000, else earliest (day, time)."""
    if not bviews:
        return None
    for d, t in bviews:
        if d.endswith("01") and t == "0000":
            return d, t
    return min(bviews, key=lambda dt: (dt[0], dt[1]))


def parse_ym(s):
    y, m = s.split("-")
    return int(y), int(m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="YYYY-MM")
    ap.add_argument("--end", required=True, help="YYYY-MM")
    ap.add_argument("--out", default="../../data/internet_as_bgp/raw")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--manifest", default="bview_manifest.tsv")
    args = ap.parse_args()

    start, end = parse_ym(args.start), parse_ym(args.end)
    os.makedirs(args.out, exist_ok=True)
    rows = []
    for y, m in months(start, end):
        chosen = pick(list_bviews(y, m))
        if not chosen:
            print(f"{y}-{m:02d}\tNO_BVIEW")
            rows.append((f"{y}-{m:02d}", "", "MISSING"))
            continue
        d, t = chosen
        fname = f"bview.{d}.{t}.gz"
        url = f"{BASE}/{y:04d}.{m:02d}/{fname}"
        dest = os.path.join(args.out, f"{y}_{m:02d}_{fname}")
        print(f"{y}-{m:02d}\t{fname}\t{'(dry-run)' if args.dry_run else url}")
        rows.append((f"{y}-{m:02d}", fname, "planned" if args.dry_run else "downloaded"))
        if not args.dry_run and not os.path.exists(dest):
            urllib.request.urlretrieve(url, dest)

    with open(os.path.join(args.out, args.manifest), "w") as fh:
        fh.write("month\tbview_file\tstatus\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"\nManifest → {os.path.join(args.out, args.manifest)}  ({len(rows)} months)")


if __name__ == "__main__":
    main()
