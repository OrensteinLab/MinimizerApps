#!/usr/bin/env python3
"""Count the k-mer windows minimap2 actually sketches under HPC (-H, which
`-x map-pb` turns on): the genome is homopolymer-compressed (each run of one
base collapsed to a single base) and split at every non-ACGT base, then every
k consecutive compressed bases is one candidate k-mer -- minimap2/sketch.c's
mm_sketch loop (`c < 4` branch, tiny_queue_t run skipping; `l = 0` reset on
an ambiguous base). This is the right denominator for density_factor:
minimap2 -d's own "average spacing" divides by the *raw* length instead, which
understates HPC density by the genome's compression ratio (~0.70 on human),
plus any N gaps (5.3% of GRCh38).

sketch.c also drops a compressed k-mer whose raw span is >= 256 bp; that needs
>= ~241 bp of homopolymer inside 15 runs and occurs 0 times in either
benchmark genome, so it isn't modeled here.

Cached next to the genome (genome.fa.hpc_k<k>.json, keyed on file size and
mtime) since a whole-genome scan takes minutes.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

_NON_ACGT = re.compile(rb"[^ACGT]+")
_HOMOPOLYMER = re.compile(rb"([ACGT])\1+")


def _contigs(genome: Path):
    chunks: list[bytes] = []
    with genome.open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                if chunks:
                    yield b"".join(chunks).upper()
                chunks = []
            else:
                chunks.append(line.rstrip(b"\r\n"))
    if chunks:
        yield b"".join(chunks).upper()


def count_hpc_kmers(genome: Path, k: int) -> dict:
    n_kmers = hpc_bp = 0
    for seq in _contigs(genome):
        for seg in _NON_ACGT.split(seq):
            h = len(_HOMOPOLYMER.sub(rb"\1", seg))
            hpc_bp += h
            n_kmers += max(0, h - k + 1)
    return {"k": k, "hpc_bp": hpc_bp, "n_hpc_kmers": n_kmers}


def hpc_kmers(genome: Path, k: int) -> dict:
    """count_hpc_kmers, cached next to the genome."""
    st = genome.stat()
    cache = genome.with_name(f"{genome.name}.hpc_k{k}.json")
    key = {"size": st.st_size, "mtime": int(st.st_mtime)}
    if cache.is_file():
        try:
            data = json.loads(cache.read_text())
            if data.get("genome") == key:
                return data
        except (OSError, ValueError):
            pass
    data = {"genome": key, **count_hpc_kmers(genome, k)}
    try:
        cache.write_text(json.dumps(data, indent=2) + "\n")
    except OSError:
        pass
    return data


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("genome", type=Path)
    ap.add_argument("-k", type=int, default=15)
    args = ap.parse_args(argv)
    print(json.dumps(hpc_kmers(args.genome, args.k), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
