#!/usr/bin/env python3
"""Download a mock-community SRA/ENA run and keep a deterministic read sample.

Streams both mates of one paired run from ENA and keeps every k-th pair
(k = read_count // --pairs, offset 0; no RNG), so the sample spans the whole
run (all tiles/lanes) and is identical on every machine. Nothing but the
sample is stored. Read headers are kept as ENA wrote them.

Writes in --out-dir:
  R1.fastq R2.fastq   the sampled pairs
  SOURCE.txt          run accession, ENA URLs, read count, stride

Usage:
  prepare_mock_reads.py --run SRR12324251 --out-dir data/mock [--pairs 250000]
"""

from __future__ import print_function

import argparse
import csv
import gzip
import io
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UA = {"User-Agent": "kraken2-download_data"}


def ena_filereport(run):
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        "?result=read_run"
        "&accession=" + urllib.parse.quote(run)
        + "&fields=run_accession,sample_title,instrument_model,library_strategy,"
        "fastq_ftp,read_count,base_count"
        "&format=tsv"
    )
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise SystemExit("ENA returned no row for " + run)
    return rows[0]


def mate_urls(row):
    parts = [p.strip() for p in (row.get("fastq_ftp") or "").split(";") if p.strip()]
    r1 = next((p for p in parts if p.endswith("_1.fastq.gz")), None)
    r2 = next((p for p in parts if p.endswith("_2.fastq.gz")), None)
    if not r1 or not r2:
        raise SystemExit("no paired FASTQ for " + row.get("run_accession", "?"))
    return ["https://" + u if not u.startswith("http") else u for u in (r1, r2)]


def stream_sample(url, stride, pairs, out_path, tries=4):
    """Write every stride-th record of a remote .fastq.gz; returns (n_kept, ids)."""
    last = None
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            ids = []
            with urllib.request.urlopen(req, timeout=300) as resp:
                gz = gzip.GzipFile(fileobj=resp)
                with open(out_path, "w") as out:
                    j = 0
                    while len(ids) < pairs:
                        rec = [gz.readline() for _ in range(4)]
                        if not rec[0]:
                            break
                        if j % stride == 0:
                            lines = [x.decode("ascii") for x in rec]
                            out.write("".join(lines))
                            ids.append(lines[0][1:].split()[0])
                        j += 1
            return ids
        except (OSError, EOFError, urllib.error.URLError) as e:
            last = e
            print("  retry {}/{} for {}: {}".format(attempt, tries, url, e), file=sys.stderr)
            time.sleep(5 * attempt)
    raise SystemExit("download failed: {}".format(last))


def base_id(tok):
    for suf in ("/1", "/2"):
        if tok.endswith(suf):
            return tok[:-2]
    return tok


def parse_args(argv):
    p = argparse.ArgumentParser(description="Deterministic subsample of a mock run.")
    p.add_argument("--run", required=True, help="SRR/ERR run accession")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pairs", type=int, default=250000)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    os.makedirs(args.out_dir, exist_ok=True)
    row = ena_filereport(args.run)
    n_pairs = int(row["read_count"])
    stride = max(1, n_pairs // args.pairs)
    u1, u2 = mate_urls(row)
    print("[*] {} ({}, {}): {} pairs; keeping every {}th (target {})".format(
        args.run, row.get("sample_title"), row.get("instrument_model"), n_pairs,
        stride, args.pairs))
    o1 = os.path.join(args.out_dir, "R1.fastq")
    o2 = os.path.join(args.out_dir, "R2.fastq")
    ids1 = stream_sample(u1, stride, args.pairs, o1 + ".part")
    ids2 = stream_sample(u2, stride, args.pairs, o2 + ".part")
    if len(ids1) != len(ids2) or any(base_id(a) != base_id(b) for a, b in zip(ids1, ids2)):
        raise SystemExit("mates out of sync after sampling ({} vs {})".format(len(ids1), len(ids2)))
    os.replace(o1 + ".part", o1)
    os.replace(o2 + ".part", o2)
    with open(os.path.join(args.out_dir, "SOURCE.txt"), "w") as f:
        f.write("run\t{}\n".format(args.run))
        f.write("sample_title\t{}\n".format(row.get("sample_title")))
        f.write("instrument\t{}\n".format(row.get("instrument_model")))
        f.write("url_r1\t{}\nurl_r2\t{}\n".format(u1, u2))
        f.write("run_pairs\t{}\nstride\t{}\nsampled_pairs\t{}\n".format(n_pairs, stride, len(ids1)))
    print("[*] wrote {} pairs to {}".format(len(ids1), args.out_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
