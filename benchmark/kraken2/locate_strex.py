#!/usr/bin/env python3
"""Locate strain-exclusion FASTA, taxonomy, and FASTQs inside an extracted tree."""

from __future__ import print_function

import argparse
import os
import sys

FASTA_NAMES = (
    "strain_excluded.fna",
    "strain_exclusion.fna",
    "strex.fna",
)
NODES_NAMES = ("nodes.dmp",)
READ_PAIRS = (
    ("accuracy_1.fq", "accuracy_2.fq"),
    ("accuracy_1.fastq", "accuracy_2.fastq"),
    ("bacteria_1.fq", "bacteria_2.fq"),
)


def find_named(root, names):
    hits = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for n in filenames:
            if n in names:
                hits.append(os.path.join(dirpath, n))
    return hits


def find_pair(root, a, b):
    for dirpath, _dirnames, filenames in os.walk(root):
        if a in filenames and b in filenames:
            return os.path.join(dirpath, a), os.path.join(dirpath, b)
    return None, None


def find_per_genome(root):
    """Directories with many *_1.fq files (Selected/bacteria style)."""
    best = []
    for dirpath, _dirnames, filenames in os.walk(root):
        ones = [f for f in filenames if f.endswith("_1.fq") or f.endswith("_1.fastq")]
        if len(ones) >= 10:
            best.append((len(ones), dirpath))
    best.sort(reverse=True)
    return [p for _, p in best]


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True)
    p.add_argument("--out", required=True, help="Write KEY=value lines")
    return p.parse_args(argv)


def first(hits):
    return hits[0] if hits else ""


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    root = args.root
    fasta = first(find_named(root, set(FASTA_NAMES)))
    if not fasta:
        # any large .fna excluding selected genomes
        cands = []
        for dirpath, _, filenames in os.walk(root):
            for n in filenames:
                if n.endswith(".fna") or n.endswith(".fa"):
                    p = os.path.join(dirpath, n)
                    try:
                        sz = os.path.getsize(p)
                    except OSError:
                        continue
                    if sz > 50 * 1024 * 1024:
                        cands.append((sz, p))
        cands.sort(reverse=True)
        fasta = cands[0][1] if cands else ""
    nodes = first(find_named(root, set(NODES_NAMES)))
    names = ""
    if nodes:
        sib = os.path.join(os.path.dirname(nodes), "names.dmp")
        if os.path.isfile(sib):
            names = sib
    r1 = r2 = ""
    for a, b in READ_PAIRS:
        p1, p2 = find_pair(root, a, b)
        if p1:
            r1, r2 = p1, p2
            break
    genome_dirs = ",".join(find_per_genome(root))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        f.write("STREX_FASTA={}\n".format(fasta))
        f.write("NODES_DMP={}\n".format(nodes))
        f.write("NAMES_DMP={}\n".format(names))
        f.write("ACC_R1={}\n".format(r1))
        f.write("ACC_R2={}\n".format(r2))
        f.write("GENOME_FQ_DIRS={}\n".format(genome_dirs))
    print("  FASTA ", fasta or "(missing)")
    print("  nodes ", nodes or "(missing)")
    print("  names ", names or "(missing)")
    print("  R1    ", r1 or "(missing)")
    print("  per-genome dirs:", genome_dirs or "(none)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
