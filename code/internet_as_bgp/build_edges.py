#!/usr/bin/env python3
"""Build monthly AS-level edge lists from RIPE RIS rrc00 `bview` MRT dumps.

PoCN data project. For each month in [start, end]:
  1. pick one `bview` table dump (prefer day 01 @ 0000, else earliest that month),
  2. download it to a scratch dir,
  3. parse it with `bgpdump -m` (TABLE_DUMP v1/v2); 1999-2000 files use a legacy
     Zebra MRT dialect bgpdump cannot read and fall back to a built-in parser,
  4. clean every AS_PATH:
       - collapse consecutive duplicate ASNs (prepending),
       - drop malformed / empty paths,
       - drop paths containing AS_SETs ("{...}"; rare, documented simplification),
       - drop paths with < 2 valid ASNs,
  5. extract consecutive ASN pairs, remove self-loops, canonicalise as (min, max),
  6. aggregate pair observation counts into integer edge weights,
  7. write `as_edges_YYYY_MM_weighted.tsv.gz` (columns: asn1 asn2 weight);
     the unweighted graph is obtained by ignoring the weight column,
  8. append per-month rows to `bview_manifest.tsv` (exact file used) and
     `cleaning_stats.tsv` (processed / discarded / retained AS_PATH counts),
  9. delete the raw dump (recent dumps are ~400 MB; keeping ~318 of them is
     not feasible — the manifest makes the exact inputs reproducible).

ASN labels are the integer AS numbers from the dumps (asplain), so they are
consistent across snapshots by construction.

Usage:
    python build_edges.py --start 1999-11 --end 2026-04
    python build_edges.py --start 2011-01 --end 2011-01 --keep-raw   # single month, keep dump

Requires: bgpdump on PATH (brew install bgpdump).
"""
import argparse
import gzip
import os
import re
import struct
from collections import Counter
import subprocess
import sys
import time
import urllib.request

from fetch_bview import BASE, list_bviews, months, parse_ym, pick

MANIFEST_COLS = "month\tbview_file\tsize_bytes\tstatus\n"
STATS_COLS = ("month\tentries_total\tpath_malformed_or_empty\tpath_with_as_set\t"
              "path_too_short\tpath_retained\tnodes\tedges\n")


class Prefetcher:
    """Download next month's dump with curl while the current one is parsed."""

    def __init__(self):
        self.proc = self.url = self.dest = None

    def start(self, url, dest):
        if self.proc is not None:
            return
        self.proc = subprocess.Popen(
            ["curl", "-sS", "--fail", "--connect-timeout", "30", "-o", dest, url],
            stderr=subprocess.DEVNULL)
        self.url, self.dest = url, dest

    def claim(self, url):
        """Return the file size if `url` was prefetched successfully, else None."""
        if self.url != url:
            self.drop()
            return None
        rc, dest = self.proc.wait(), self.dest
        self.proc = self.url = self.dest = None
        if rc == 0 and os.path.exists(dest):
            return os.path.getsize(dest)
        if os.path.exists(dest):
            os.remove(dest)
        return None

    def drop(self):
        if self.proc is not None and self.proc.poll() is None:
            self.proc.kill()
            self.proc.wait()
        if self.dest and os.path.exists(self.dest):
            os.remove(self.dest)
        self.proc = self.url = self.dest = None


def head_size(url):
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            return int(r.headers.get("Content-Length", 0))
    except Exception:
        return 0


def download(url, dest, retries=3):
    # curl is an order of magnitude faster than urllib here (~60 vs ~2.5 MB/s)
    for attempt in range(1, retries + 1):
        r = subprocess.run(["curl", "-sS", "--fail", "--connect-timeout", "30",
                            "-o", dest, url], stderr=subprocess.DEVNULL)
        if r.returncode == 0:
            return os.path.getsize(dest)
        print(f"  ! download attempt {attempt}/{retries} failed (curl {r.returncode})",
              file=sys.stderr)
        if os.path.exists(dest):
            os.remove(dest)
        time.sleep(5 * attempt)
    return None


def clean_path(tokens):
    """Steps 4-7 on a tokenised AS_PATH. Returns list of ASNs or None if dropped."""
    path = []
    for tok in tokens:
        if not tok.isdigit():          # AS_SET "{...}" or garbage -> drop whole path
            return None
        asn = int(tok)
        if not path or path[-1] != asn:  # collapse consecutive duplicates
            path.append(asn)
    return path if len(path) >= 2 else []


def iter_paths_bgpdump(path):
    """Yield one raw AS_PATH string per RIB entry via bgpdump -m."""
    proc = subprocess.Popen(["bgpdump", "-m", path], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    for line in proc.stdout:
        # bgpdump -m: proto|ts|B|peer_ip|peer_asn|prefix|AS_PATH|origin|...
        fields = line.split("|")
        if len(fields) >= 7:
            yield fields[6].strip()
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"bgpdump exited with {proc.returncode} on {path}")


def _as_path_from_attrs(attrs):
    """Extract the AS_PATH (attribute type 2, 16-bit ASNs) from BGP path attributes."""
    parts, p = [], 0
    while p + 3 <= len(attrs):
        flags, atype = attrs[p], attrs[p + 1]
        if flags & 0x10:                    # extended length
            alen = struct.unpack_from(">H", attrs, p + 2)[0]
            val, p = attrs[p + 4:p + 4 + alen], p + 4 + alen
        else:
            alen = attrs[p + 2]
            val, p = attrs[p + 3:p + 3 + alen], p + 3 + alen
        if atype != 2:
            continue
        q = 0
        while q + 2 <= len(val):
            seg_type, seg_cnt = val[q], val[q + 1]
            asns = struct.unpack_from(f">{seg_cnt}H", val, q + 2)
            q += 2 + 2 * seg_cnt
            if seg_type == 2:               # AS_SEQUENCE
                parts.extend(str(a) for a in asns)
            else:                           # AS_SET & friends -> keep the marker
                parts.append("{%s}" % ",".join(str(a) for a in asns))
    return " ".join(parts)


def iter_paths_legacy(path):
    """Yield AS_PATH strings from legacy Zebra `bview` files (1999-2000).

    These predate the RFC 6396 layout and bgpdump reads only their first entry:
    each TABLE_DUMP record body is (length field + 4) bytes long and packs
    view(2) seq(2) followed by many entries
    prefix(4|16) mask(1) status(1) uptime(4) peer_ip(4|16) peer_as(2) attr_len(2) attrs.
    """
    with gzip.open(path, "rb") as fh:
        data = fh.read()
    off = 0
    while off + 12 <= len(data):
        _, typ, sub, ln = struct.unpack_from(">IHHI", data, off)
        body = data[off + 12:off + 12 + ln + 4]
        off += 12 + ln + 4
        if typ != 12 or sub not in (1, 2):  # only TABLE_DUMP, AFI IPv4/IPv6
            continue
        addr = 4 if sub == 1 else 16
        p = 4                               # skip view + sequence
        while p + 2 * addr + 10 <= len(body):
            q = p + 2 * addr + 6            # start of peer_as
            alen = struct.unpack_from(">H", body, q + 2)[0]
            yield _as_path_from_attrs(body[q + 4:q + 4 + alen])
            p = q + 4 + alen


def iter_paths_gated(path):
    """Yield AS_PATH strings from early-2000 GateD text table dumps.

    A few 2000 months are gzipped plain-text routing tables (`View #0 ip
    unicast ...` header; one route per line ending in the AS path plus an
    origin code) instead of binary MRT.
    """
    route_re = re.compile(r"^[s*>iaxdh ]+B\s+\d+\s+[\d:]+\s+\S+\s+\S+\s+\S+\s*(.*?)\s*[ie?a]$")
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            m = route_re.match(line.rstrip())
            if m:
                yield m.group(1)


def process_paths(path_iter):
    """Aggregate cleaned AS pairs from raw AS_PATH strings. Returns (stats, edges).

    The same AS_PATH string recurs across peers and prefixes, so paths are
    tallied first and each distinct one is cleaned once (~2x faster on the
    ~55M-entry recent dumps).
    """
    stats = dict(entries_total=0, malformed=0, as_set=0, too_short=0, retained=0)
    edges = {}
    for raw, n in Counter(path_iter).items():
        stats["entries_total"] += n
        if not raw:
            stats["malformed"] += n
            continue
        if "{" in raw:                      # AS_SET anywhere in the path
            stats["as_set"] += n
            continue
        path_asns = clean_path(raw.split())
        if path_asns is None:
            stats["malformed"] += n
            continue
        if not path_asns:
            stats["too_short"] += n
            continue
        stats["retained"] += n
        for a, b in zip(path_asns, path_asns[1:]):
            if a == b:                      # self-loop (cannot occur after collapsing)
                continue
            e = (a, b) if a < b else (b, a)
            edges[e] = edges.get(e, 0) + n
    return stats, edges


MIN_PATHS = 50_000

def process_dump(path):
    stats, edges = process_paths(iter_paths_bgpdump(path))
    if stats["entries_total"] < 1000:       # not modern MRT: legacy Zebra or GateD text
        with gzip.open(path, "rb") as fh:
            head = fh.read(4096)
        text = b"View #" in head or b"Destination" in head
        stats, edges = process_paths(iter_paths_gated(path) if text
                                     else iter_paths_legacy(path))
        stats["gated_text" if text else "legacy_format"] = True
    if stats["retained"] < MIN_PATHS:       # a real full table has far more entries;
        raise ValueError(                   # fewer means a truncated/partial dump
            f"only {stats['retained']} usable AS_PATHs - truncated dump?")
    return stats, edges


def write_edges(edges, out_path):
    nodes = set()
    with gzip.open(out_path, "wt") as fh:
        fh.write("asn1\tasn2\tweight\n")
        for (a, b), w in sorted(edges.items()):
            fh.write(f"{a}\t{b}\t{w}\n")
            nodes.add(a)
            nodes.add(b)
    return len(nodes)


def upsert_row(path, header, row):
    """Insert/replace the row keyed by its first column (month); keep file sorted."""
    rows = {}
    if os.path.exists(path):
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                if line.strip():
                    rows[line.split("\t", 1)[0]] = line.rstrip("\n")
    rows[str(row[0])] = "\t".join(str(x) for x in row)
    with open(path, "w") as fh:
        fh.write(header)
        for k in sorted(rows):
            fh.write(rows[k] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="YYYY-MM")
    ap.add_argument("--end", required=True, help="YYYY-MM")
    ap.add_argument("--data-dir", default="../../data/internet_as_bgp")
    ap.add_argument("--scratch", default=None, help="dir for raw dumps (default <data-dir>/raw)")
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--overwrite", action="store_true", help="rebuild months whose edge file exists")
    ap.add_argument("--min-paths", type=int, default=50_000,
                    help="reject dumps with fewer usable AS_PATHs (partial-dump guard)")
    args = ap.parse_args()
    global MIN_PATHS
    MIN_PATHS = args.min_paths

    data_dir = os.path.abspath(args.data_dir)
    scratch = os.path.abspath(args.scratch or os.path.join(data_dir, "raw"))
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(scratch, exist_ok=True)
    manifest = os.path.join(data_dir, "bview_manifest.tsv")
    stats_log = os.path.join(data_dir, "cleaning_stats.tsv")

    month_list = list(months(parse_ym(args.start), parse_ym(args.end)))
    out_for = lambda y, m: os.path.join(data_dir, f"as_edges_{y}_{m:02d}_weighted.tsv.gz")
    raw_for = lambda y, m, d, t: os.path.join(scratch, f"{y}_{m:02d}_bview.{d}.{t}.gz")
    bviews_cache = {}

    def month_bviews(y, m):
        if (y, m) not in bviews_cache:
            bviews_cache[(y, m)] = list_bviews(y, m)
        return bviews_cache[(y, m)]

    prefetch = Prefetcher()

    def prefetch_next(y, m):
        """Start downloading the next pending month's preferred dump."""
        for y2, m2 in month_list:
            if (y2, m2) <= (y, m) or (os.path.exists(out_for(y2, m2))
                                      and not args.overwrite):
                continue
            nxt = pick(month_bviews(y2, m2))
            if nxt:
                d, t = nxt
                prefetch.start(f"{BASE}/{y2:04d}.{m2:02d}/bview.{d}.{t}.gz",
                               raw_for(y2, m2, d, t))
            return

    for y, m in month_list:
        month = f"{y}-{m:02d}"
        out_path = out_for(y, m)
        if os.path.exists(out_path) and not args.overwrite:
            print(f"{month}\tSKIP (exists)")
            continue
        t0 = time.time()
        bviews = month_bviews(y, m)
        if not bviews:
            print(f"{month}\tNO_BVIEW")
            upsert_row(manifest, MANIFEST_COLS, [month, "", 0, "MISSING"])
            continue
        # preferred file first; if it is corrupt/truncated, retry with the
        # remaining files of the month ranked by size (bad dumps are small)
        url_for = lambda b: f"{BASE}/{y:04d}.{m:02d}/bview.{b[0]}.{b[1]}.gz"
        stats = edges = None
        queue, tried = [pick(bviews)], set()
        while queue and stats is None and len(tried) < 5:
            d, t = queue.pop(0)
            tried.add((d, t))
            fname = f"bview.{d}.{t}.gz"
            raw_path = raw_for(y, m, d, t)
            # a prefetched file first; then --keep-raw reuse (an untrusted
            # pre-existing file may be a partial download); else download now
            size = prefetch.claim(url_for((d, t)))
            if size is None:
                size = os.path.getsize(raw_path) if args.keep_raw and os.path.exists(raw_path) \
                    else download(url_for((d, t)), raw_path)
            if size is None:
                print(f"{month}\tDOWNLOAD_FAILED\t{fname}")
                upsert_row(manifest, MANIFEST_COLS, [month, fname, 0, "DOWNLOAD_FAILED"])
            else:
                prefetch_next(y, m)         # overlap next download with this parse
                try:
                    stats, edges = process_dump(raw_path)
                except Exception as e:
                    print(f"{month}\tPARSE_FAILED\t{fname}\t{e}")
                    upsert_row(manifest, MANIFEST_COLS, [month, fname, size, "PARSE_FAILED"])
                finally:
                    if not args.keep_raw and os.path.exists(raw_path):
                        os.remove(raw_path)
            if stats is None and not queue:
                rest = [b for b in bviews if b not in tried]
                queue = sorted(rest, key=lambda b: head_size(url_for(b)), reverse=True)
        if stats is None:
            continue
        n_nodes = write_edges(edges, out_path)
        upsert_row(manifest, MANIFEST_COLS, [month, fname, size, "OK"])
        upsert_row(stats_log, STATS_COLS,
                   [month, stats["entries_total"], stats["malformed"], stats["as_set"],
                    stats["too_short"], stats["retained"], n_nodes, len(edges)])
        print(f"{month}\t{fname}\t{size/1e6:.0f} MB\t{stats['retained']}/{stats['entries_total']} paths"
              f"\t{n_nodes} nodes\t{len(edges)} edges\t{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
