#!/usr/bin/env python3
"""Download FDA-ARGOS paired FASTQ slices for Table S5 SRR accessions.

Reads benchmark/kraken2/tables/table_s5.csv (species, sra_run). Does not invent a
run list. Headers are rewritten as ARGOS.<SRR>.<i> for scoring.

Usage:
  prepare_fda_argos_reads.py --table-s5 benchmark/kraken2/tables/table_s5.csv \\
      --out-dir /path/data/fda-argos [--pairs 10000]
"""

from __future__ import print_function

import argparse
import csv
import gzip
import io
import os
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UA = {"User-Agent": "kraken2-download_data"}
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TABLE = os.path.join(HERE, "tables", "table_s5.csv")


def ena_filereport(run):
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        "?result=read_run"
        "&accession=" + urllib.parse.quote(run)
        + "&fields=run_accession,scientific_name,fastq_ftp,read_count"
        "&format=tsv"
    )
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise RuntimeError("ENA returned no row for " + run)
    return rows[0]


def fastq_urls(row):
    ftp = (row.get("fastq_ftp") or "").strip()
    parts = [p.strip() for p in ftp.split(";") if p.strip()]
    r1 = next((p for p in parts if p.endswith("_1.fastq.gz")), None)
    r2 = next((p for p in parts if p.endswith("_2.fastq.gz")), None)
    if not r1 or not r2:
        raise RuntimeError("no paired FASTQ FTP for " + row.get("run_accession", "?"))
    if not r1.startswith("http"):
        r1 = "https://" + r1
    if not r2.startswith("http"):
        r2 = "https://" + r2
    return r1, r2


def open_gz_url(url):
    req = urllib.request.Request(url, headers=UA)
    resp = urllib.request.urlopen(req, timeout=600)
    return gzip.GzipFile(fileobj=resp), resp


def read_fastq_record(fh):
    lines = [fh.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if any(not x for x in lines):
        raise EOFError("truncated FASTQ record")
    return lines


def take_first_pairs(r1_url, r2_url, k):
    gz1 = gz2 = resp1 = resp2 = None
    try:
        gz1, resp1 = open_gz_url(r1_url)
        gz2, resp2 = open_gz_url(r2_url)
        t1 = io.TextIOWrapper(gz1, encoding="ascii", errors="replace")
        t2 = io.TextIOWrapper(gz2, encoding="ascii", errors="replace")
        kept = []
        while len(kept) < k:
            a = read_fastq_record(t1)
            b = read_fastq_record(t2)
            if a is None and b is None:
                break
            if a is None or b is None:
                raise EOFError("R1/R2 length mismatch while streaming")
            kept.append((a, b))
        return kept, len(kept)
    finally:
        for x in (gz1, gz2, resp1, resp2):
            try:
                if x is not None:
                    x.close()
            except Exception:
                pass


def take_first_pairs_retry(r1_url, r2_url, k, attempts=5):
    last = None
    for attempt in range(1, attempts + 1):
        try:
            return take_first_pairs(r1_url, r2_url, k)
        except (
            urllib.error.URLError,
            TimeoutError,
            ConnectionError,
            EOFError,
            OSError,
        ) as e:
            last = e
            wait = min(60, 2**attempt)
            print(
                "    retry {}/{} after {}: {} (sleep {}s)".format(
                    attempt, attempts, type(e).__name__, e, wait
                ),
                flush=True,
            )
            time.sleep(wait)
    raise RuntimeError("failed after {} attempts: {}".format(attempts, last))


def reservoir_pairs(r1_url, r2_url, k, seed):
    rng = random.Random(seed)
    gz1 = gz2 = resp1 = resp2 = None
    try:
        gz1, resp1 = open_gz_url(r1_url)
        gz2, resp2 = open_gz_url(r2_url)
        t1 = io.TextIOWrapper(gz1, encoding="ascii", errors="replace")
        t2 = io.TextIOWrapper(gz2, encoding="ascii", errors="replace")
        kept = []
        scanned = 0
        while True:
            a = read_fastq_record(t1)
            b = read_fastq_record(t2)
            if a is None and b is None:
                break
            if a is None or b is None:
                raise EOFError("R1/R2 length mismatch")
            scanned += 1
            if len(kept) < k:
                kept.append((a, b))
            else:
                j = rng.randrange(scanned)
                if j < k:
                    kept[j] = (a, b)
        return kept, scanned
    finally:
        for x in (gz1, gz2, resp1, resp2):
            try:
                if x is not None:
                    x.close()
            except Exception:
                pass


def rewrite_pair(a, b, run, i):
    def one(rec, mate):
        rec = list(rec)
        rest = rec[0][1:].rstrip("\n")
        rec[0] = "@ARGOS.{}.{}/{} {}\n".format(run, i, mate, rest)
        return rec

    return one(a, 1), one(b, 2)


def load_table_s5(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            run = (row.get("sra_run") or row.get("run_accession") or "").strip()
            species = (row.get("species") or "").strip()
            if not run:
                continue
            rows.append({"sra_run": run, "species": species})
    if not rows:
        sys.exit("no sra_run rows in " + path)
    return rows


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Download Table S5 FDA-ARGOS PE FASTQs from ENA."
    )
    p.add_argument("--table-s5", default=DEFAULT_TABLE)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pairs", type=int, default=10000)
    p.add_argument("--reservoir", action="store_true")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--force", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    os.makedirs(args.out_dir, exist_ok=True)
    fwd = os.path.join(args.out_dir, "R1.fastq")
    rev = os.path.join(args.out_dir, "R2.fastq")
    if (
        not args.force
        and os.path.isfile(fwd)
        and os.path.isfile(rev)
        and os.path.getsize(fwd) > 0
    ):
        print("  FASTQs already present in {} (use --force to refresh)".format(args.out_dir))
        return 0

    experiments = load_table_s5(args.table_s5)
    tmp_fwd, tmp_rev = fwd + ".partial", rev + ".partial"
    with open(tmp_fwd, "w", encoding="ascii") as o1, open(
        tmp_rev, "w", encoding="ascii"
    ) as o2:
        total = 0
        for i, exp in enumerate(experiments, 1):
            run = exp["sra_run"]
            species = exp["species"]
            print(
                "  [{}/{}] {} {}".format(i, len(experiments), run, species),
                flush=True,
            )
            row = ena_filereport(run)
            r1, r2 = fastq_urls(row)
            if args.reservoir:
                last = None
                kept = scanned = None
                for attempt in range(1, 6):
                    try:
                        kept, scanned = reservoir_pairs(r1, r2, args.pairs, args.seed)
                        break
                    except (
                        urllib.error.URLError,
                        TimeoutError,
                        ConnectionError,
                        EOFError,
                        OSError,
                    ) as e:
                        last = e
                        time.sleep(min(60, 2**attempt))
                if kept is None:
                    raise RuntimeError("reservoir failed for {}: {}".format(run, last))
            else:
                kept, scanned = take_first_pairs_retry(r1, r2, args.pairs)
            for j, (a, b) in enumerate(kept, 1):
                a2, b2 = rewrite_pair(a, b, run, j)
                o1.writelines(a2)
                o2.writelines(b2)
            total += len(kept)
            print("    kept {} pairs (scanned {})".format(len(kept), scanned), flush=True)
    os.replace(tmp_fwd, fwd)
    os.replace(tmp_rev, rev)
    print("  Wrote {} PE pairs → {} / {}".format(total, fwd, rev))
    print("  Runs from {}".format(args.table_s5))
    return 0


if __name__ == "__main__":
    sys.exit(main())
