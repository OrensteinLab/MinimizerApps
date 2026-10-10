#!/usr/bin/env python3
"""Uniform random sample of N reads (minimum length L) from a FASTQ stream on
stdin, by reservoir sampling (Algorithm R) with a fixed seed -- one pass, no
need to store the input. Built for the Winnowmap paper's D3: "a small random
sample of the rel3 human CHM13 whole-genome, ultra-long ONT sequencing data"
(1,000 reads, minimum read length 1 kbp, mapped to GRCh38), whose 355 GB
FASTQ is streamed rather than downloaded:

  curl -sL https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/nanopore/rel3/rel3.fastq.gz \\
    | tee >(md5sum > rel3.md5) | pigz -dc \\
    | python3 lib/sample_reads.py --n 1000 --min-len 1000 --seed 1 --out-dir data/chm13_rel3_ont

--min-mean-q Q additionally requires a read mean quality >= Q, as basecallers
define it: -10 log10 of the mean per-base error probability (the ONT "pass"
threshold is 7; D3 uses it, see README).

Writes reads.fa (sampled reads, input order), sampled_ids.txt and
sample_stats.json (eligible/total reads seen, N50, max length).
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path


def n50(lengths: list[int]) -> int:
    total, acc = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        acc += x
        if 2 * acc >= total:
            return x
    return 0


_ERR = [10 ** (-(q - 33) / 10) for q in range(128)]


def mean_quality(qual: bytes) -> float:
    """-10 log10(mean error probability) of a Phred+33 quality string."""
    total = sum(_ERR[c] * qual.count(c) for c in set(qual))
    return -10 * math.log10(total / len(qual))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-mean-q", type=float, default=None,
                    help="keep only reads with mean quality >= this (default: no filter)")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    reservoir: list[tuple[int, str, str]] = []  # (input index, name, sequence)
    n_total = n_eligible = 0
    stdin = sys.stdin.buffer
    while True:
        header = stdin.readline()
        if not header:
            break
        seq = stdin.readline().rstrip(b"\n")
        stdin.readline()  # '+'
        qual = stdin.readline().rstrip(b"\n")
        n_total += 1
        if len(seq) < args.min_len:
            continue
        if args.min_mean_q is not None and mean_quality(qual) < args.min_mean_q:
            continue
        if n_eligible < args.n:
            reservoir.append((n_eligible, header[1:].split()[0].decode(), seq.decode()))
        else:
            j = rng.randrange(n_eligible + 1)
            if j < args.n:
                reservoir[j] = (n_eligible, header[1:].split()[0].decode(), seq.decode())
        n_eligible += 1
        if n_total % 1_000_000 == 0:
            print(f"{n_total} reads seen, {n_eligible} eligible", file=sys.stderr, flush=True)

    reservoir.sort()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "reads.fa").open("w") as f:
        for _, name, seq in reservoir:
            f.write(f">{name}\n{seq}\n")
    (args.out_dir / "sampled_ids.txt").write_text("".join(f"{name}\n" for _, name, _ in reservoir))
    lengths = [len(s) for _, _, s in reservoir]
    stats = {
        "n_reads_seen": n_total, "n_eligible": n_eligible, "min_len": args.min_len,
        "min_mean_q": args.min_mean_q, "seed": args.seed,
        "n_sampled": len(reservoir), "n50": n50(lengths), "max_len": max(lengths, default=0),
        "min_sampled_len": min(lengths, default=0),
    }
    (args.out_dir / "sample_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
