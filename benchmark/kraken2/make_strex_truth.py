#!/usr/bin/env python3
"""Build truth.tsv from Mason-style or ARGOS-prefixed FASTQ headers.

Mason:  --table-s3
ARGOS:  --table-s5 --names-dmp --nodes-dmp   (seqids ARGOS.<SRR>.<i>)

Usage:
  make_strex_truth.py --fastq R1.fastq --table-s3 table_s3.csv --out truth.tsv
  make_strex_truth.py --fastq R1.fastq --table-s5 table_s5.csv \\
      --names-dmp names.dmp --nodes-dmp nodes.dmp --out truth.tsv
"""

from __future__ import print_function

import argparse
import csv
import os
import re
import sys

TAXID_RE = re.compile(r"taxid_(\d+)", re.I)
ARGOS_RE = re.compile(r"^ARGOS\.([A-Z0-9]+)\.(\d+)", re.I)
KRAKEN_TAX_RE = re.compile(r"kraken:taxid\|(\d+)")
LEADING_TAX_RE = re.compile(r"^(\d+)(?:[._]|$)")


def seq_id_from_header(line):
    tok = line[1:].split()[0] if line.startswith(">") or line.startswith("@") else line.split()[0]
    if tok.endswith("/1") or tok.endswith("/2"):
        tok = tok[:-2]
    return tok


def load_s3(path):
    by_taxid = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            by_taxid[int(row["taxid"])] = row
    return by_taxid


def write_truth(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("#seqid\ttaxid\tdomain\n")
        for seqid, taxid, domain in rows:
            f.write("{}\t{}\t{}\n".format(seqid, taxid, domain))


def parse_taxid(seqid, header):
    for rx, src in (
        (TAXID_RE, seqid),
        (TAXID_RE, header),
        (KRAKEN_TAX_RE, header),
        (LEADING_TAX_RE, seqid),
    ):
        m = rx.search(src)
        if m:
            return int(m.group(1))
    return None


def from_mason(fastq, s3):
    rows = []
    n = 0
    with open(fastq) as f:
        for i, line in enumerate(f):
            if i % 4 != 0:
                continue
            header = line.rstrip("\n")
            seqid = seq_id_from_header(header)
            taxid = parse_taxid(seqid, header)
            if taxid is None:
                sys.exit("no taxid in FASTQ header: " + header[:160])
            domain = s3[taxid]["domain"] if taxid in s3 else "unknown"
            rows.append((seqid, taxid, domain))
            n += 1
    return rows, n


def load_names_dmp(path):
    sci = {}
    other = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t|\t")
            if len(parts) < 4:
                continue
            taxid = int(parts[0])
            name = parts[1].strip()
            nclass = parts[3].rstrip("\t|").strip()
            key = name.lower()
            bucket = sci if nclass == "scientific name" else other
            bucket.setdefault(key, []).append(taxid)
    return sci, other


def load_ranks(nodes_path):
    rank = {}
    parent = {}
    with open(nodes_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t|\t")
            if len(parts) < 3:
                continue
            tid = int(parts[0])
            parent[tid] = int(parts[1])
            rank[tid] = parts[2]
    return parent, rank


def elevate(parent, rank_map, taxid, want):
    while taxid > 0:
        if rank_map.get(taxid) == want:
            return taxid
        taxid = parent.get(taxid, 0)
    return 0


def resolve_species(name, sci, other, parent, rank_map):
    key = name.lower().strip()
    candidates = sci.get(key) or other.get(key) or []
    if not candidates and " subsp." in key:
        candidates = sci.get(key.split(" subsp.")[0]) or []
    scored = []
    for tid in candidates:
        sp = elevate(parent, rank_map, tid, "species") or tid
        scored.append((0 if rank_map.get(tid) == "species" else 1, sp, tid))
    if not scored:
        return 0
    scored.sort()
    tid = scored[0][2]
    return elevate(parent, rank_map, tid, "species") or tid


def from_argos_table(fastq, table_s5, names_dmp, nodes_dmp):
    sci, other = load_names_dmp(names_dmp)
    parent, rank_map = load_ranks(nodes_dmp)
    run_to = {}
    with open(table_s5, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            run = (row.get("sra_run") or row.get("run_accession") or "").strip()
            species = (row.get("species") or "").strip()
            if not run:
                continue
            taxid = resolve_species(species, sci, other, parent, rank_map)
            run_to[run] = (taxid, "prokaryote")
    rows = []
    with open(fastq) as f:
        for i, line in enumerate(f):
            if i % 4 != 0:
                continue
            seqid = seq_id_from_header(line.rstrip("\n"))
            m = ARGOS_RE.match(seqid)
            if not m:
                sys.exit("expected ARGOS.<RUN>.<n> header, got: " + seqid)
            run = m.group(1)
            if run not in run_to:
                sys.exit("run {} not in {}".format(run, table_s5))
            taxid, domain = run_to[run]
            rows.append((seqid, taxid, domain))
    return rows, len(rows)


def from_argos(fastq, manifest):
    run_to = {}
    with open(manifest, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            run = row.get("run_accession") or row.get("sra_run")
            taxid = row.get("species_taxid") or row.get("taxid") or "0"
            domain = row.get("domain") or "prokaryote"
            run_to[run] = (int(taxid), domain)
    rows = []
    with open(fastq) as f:
        for i, line in enumerate(f):
            if i % 4 != 0:
                continue
            seqid = seq_id_from_header(line.rstrip("\n"))
            m = ARGOS_RE.match(seqid)
            if not m:
                sys.exit("expected ARGOS.<RUN>.<n> header, got: " + seqid)
            run = m.group(1)
            if run not in run_to:
                sys.exit("run {} not in manifest".format(run))
            taxid, domain = run_to[run]
            rows.append((seqid, taxid, domain))
    return rows, len(rows)


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--fastq", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--table-s3", default="")
    p.add_argument("--table-s5", default="")
    p.add_argument("--names-dmp", default="")
    p.add_argument("--nodes-dmp", default="")
    p.add_argument("--argos-manifest", default="", help="legacy TSV; prefer --table-s5")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    if args.table_s5:
        if not args.names_dmp or not args.nodes_dmp:
            sys.exit("--table-s5 requires --names-dmp and --nodes-dmp")
        rows, n = from_argos_table(
            args.fastq, args.table_s5, args.names_dmp, args.nodes_dmp
        )
    elif args.argos_manifest:
        rows, n = from_argos(args.fastq, args.argos_manifest)
    else:
        if not args.table_s3:
            sys.exit("pass --table-s3 (Mason) or --table-s5 (ARGOS)")
        rows, n = from_mason(args.fastq, load_s3(args.table_s3))
    write_truth(args.out, rows)
    print("  {} fragments → {}".format(n, args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
