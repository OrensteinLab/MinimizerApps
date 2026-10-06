#!/usr/bin/env python3
"""Deterministically subsample a CAMI II short-read sample and build truth.

Input is one CAMI II simulated sample (interleaved anonymous_reads.fq.gz,
reads_mapping.tsv.gz, and the dataset's metadata.tsv). Every k-th read pair is
kept (k = total pairs // --pairs, offset 0; no RNG), so the subsample is the
same on every machine. Original read IDs are kept (S0R<i>/1, S0R<i>/2) so the
CAMI gold-standard mapping still applies.

Writes in --out-dir:
  R1.fastq R2.fastq   the sampled pairs
  truth.tsv           #seqid taxid domain   (one row per pair; CAMI tax_id)
  truth_extra.tsv     #seqid taxid genome_id novelty_category (CAMI metadata)

The CAMI tax_id is the NCBI taxid (2017-02 taxdump) CAMI assigned to the
source genome; for new_species / new_genus genomes it is a higher-rank taxon.
The domain column is filled from --nodes when given, else "all".

Usage:
  prepare_cami_reads.py --reads anonymous_reads.fq.gz \\
      --mapping reads_mapping.tsv.gz --metadata metadata.tsv \\
      --out-dir data/cami [--pairs 250000] [--nodes nodes.dmp]
"""

from __future__ import print_function

import argparse
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def count_pairs(mapping):
    """Pairs = mapping rows / 2 (the mapping lists both mates)."""
    n = 0
    with gzip.open(mapping, "rt") as f:
        next(f)  # header
        for _ in f:
            n += 1
    return n // 2


def read_id(header):
    tok = header[1:].split()[0]
    if tok.endswith("/1") or tok.endswith("/2"):
        tok = tok[:-2]
    return tok


def sample_pairs(reads, stride, pairs, out1, out2):
    """Write every stride-th pair; returns the ordered list of kept IDs."""
    kept = []
    with gzip.open(reads, "rt") as f, open(out1, "w") as o1, open(out2, "w") as o2:
        j = 0
        while len(kept) < pairs:
            r1 = [f.readline() for _ in range(4)]
            if not r1[0]:
                break
            r2 = [f.readline() for _ in range(4)]
            if not r2[0]:
                raise SystemExit("odd number of reads: not interleaved?")
            id1, id2 = read_id(r1[0]), read_id(r2[0])
            if id1 != id2 or not r1[0].rstrip().endswith("/1") or not r2[0].rstrip().endswith("/2"):
                raise SystemExit("reads not interleaved /1 /2 at pair {}: {} {}".format(
                    j, r1[0].strip(), r2[0].strip()))
            if j % stride == 0:
                o1.write("".join(r1))
                o2.write("".join(r2))
                kept.append(id1)
            j += 1
    return kept


def load_mapping(mapping, wanted):
    """seqid -> (genome_id, taxid) for wanted pair IDs; mates must agree."""
    out = {}
    with gzip.open(mapping, "rt") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            sid = read_id("@" + parts[0])
            if sid not in wanted:
                continue
            rec = (parts[1], int(parts[2]))
            if sid in out and out[sid] != rec:
                raise SystemExit("mates of {} disagree in the mapping".format(sid))
            out[sid] = rec
    return out


def load_novelty(metadata):
    """genome_ID -> novelty_category from CAMI metadata.tsv."""
    nov = {}
    with open(metadata) as f:
        header = f.readline().rstrip("\n").split("\t")
        gi = header.index("genome_ID") if "genome_ID" in header else 0
        ni = header.index("novelty_category")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > ni:
                nov[p[gi]] = p[ni]
    return nov


def parse_args(argv):
    p = argparse.ArgumentParser(description="Subsample a CAMI II sample + truth.")
    p.add_argument("--reads", required=True)
    p.add_argument("--mapping", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pairs", type=int, default=250000)
    p.add_argument("--nodes", default="", help="nodes.dmp to fill the domain column")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    os.makedirs(args.out_dir, exist_ok=True)
    total = count_pairs(args.mapping)
    stride = max(1, total // args.pairs)
    print("[*] {} pairs in the sample; keeping every {}th (target {})".format(
        total, stride, args.pairs))
    out1 = os.path.join(args.out_dir, "R1.fastq")
    out2 = os.path.join(args.out_dir, "R2.fastq")
    kept = sample_pairs(args.reads, stride, args.pairs, out1 + ".part", out2 + ".part")
    mapping = load_mapping(args.mapping, set(kept))
    missing = [s for s in kept if s not in mapping]
    if missing:
        raise SystemExit("{} sampled reads missing from the mapping, e.g. {}".format(
            len(missing), missing[0]))

    parent = None
    domain_of = None
    if args.nodes and os.path.isfile(args.nodes):
        from evaluate_calls import load_taxonomy
        from score_cami import domain_of
        parent, _ranks = load_taxonomy(args.nodes)
    novelty = load_novelty(args.metadata)

    with open(os.path.join(args.out_dir, "truth.tsv"), "w") as ft, \
            open(os.path.join(args.out_dir, "truth_extra.tsv"), "w") as fe:
        ft.write("#seqid\ttaxid\tdomain\n")
        fe.write("#seqid\ttaxid\tgenome_id\tnovelty_category\n")
        for sid in kept:
            genome, taxid = mapping[sid]
            dom = domain_of(parent, taxid) if parent is not None else "all"
            ft.write("{}\t{}\t{}\n".format(sid, taxid, dom))
            fe.write("{}\t{}\t{}\t{}\n".format(sid, taxid, genome, novelty.get(genome, "NA")))
    os.replace(out1 + ".part", out1)
    os.replace(out2 + ".part", out2)
    print("[*] wrote {} pairs to {}".format(len(kept), args.out_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
