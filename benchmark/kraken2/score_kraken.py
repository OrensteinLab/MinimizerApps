#!/usr/bin/env python3
"""Score Kraken 2 output at genus and species (Wood et al. 2019).

Reads standard Kraken output (C/U, seqid, taxid, ...) and a truth TSV
(seqid, taxid, [domain]). Writes a metrics TSV and optional per-rank
call files.

Usage:
  score_kraken.py --nodes nodes.dmp --truth truth.tsv --kraken out.kraken \\
      --out-tsv metrics.tsv [--expected-s1 table_s1.csv]
"""

from __future__ import print_function

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from evaluate_calls import load_taxonomy, score_calls, format_row

RANKS = ("genus", "species")


def strip_mate(seqid):
    if seqid.endswith("/1") or seqid.endswith("/2"):
        return seqid[:-2]
    return seqid


def load_truth(path):
    truth = {}
    domain = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                fields = line.split()
            seqid = strip_mate(fields[0])
            taxid = int(fields[1])
            truth[seqid] = taxid
            domain[seqid] = fields[2] if len(fields) > 2 else "all"
    return truth, domain


def load_kraken_calls(path):
    calls = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seqid = strip_mate(parts[1])
            try:
                taxid = int(parts[2])
            except ValueError:
                taxid = 0
            if parts[0] == "U":
                taxid = 0
            calls[seqid] = taxid
    return calls


def subset_map(d, keys):
    return {k: d[k] for k in keys if k in d}


def write_calls(path, calls, keys):
    with open(path, "w") as f:
        for k in keys:
            f.write("{}\t{}\n".format(k, calls.get(k, 0)))


def expected_s1_row(path, classifier="Kraken 2"):
    if not path or not os.path.isfile(path):
        return None
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("classifier") == classifier:
                return row
    return None


def parse_args(argv):
    p = argparse.ArgumentParser(description="Score Kraken output vs strain-exclusion truth.")
    p.add_argument("--nodes", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--kraken", required=True, help="Kraken 2 --output file")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--scheme", default="")
    p.add_argument("--minimum-hit-groups", type=int, default=None,
                   help="classify --minimum-hit-groups used for this output; "
                        "recorded in the minimum_hit_groups column")
    p.add_argument("--dataset", default="strain-exclusion")
    p.add_argument("--calls-dir", default="", help="Write per-domain call TSVs here")
    p.add_argument("--expected-s1", default="", help="table_s1.csv for Kraken 2 genus check")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    parent, rank_map = load_taxonomy(args.nodes)
    truth, domains = load_truth(args.truth)
    calls = load_kraken_calls(args.kraken)
    groups = {"all": list(truth)}
    for seqid, dom in domains.items():
        groups.setdefault(dom, []).append(seqid)

    rows = []
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    if args.calls_dir:
        os.makedirs(args.calls_dir, exist_ok=True)

    for domain, keys in sorted(groups.items()):
        t = subset_map(truth, keys)
        c = subset_map(calls, keys)
        if args.calls_dir:
            write_calls(
                os.path.join(args.calls_dir, "{}_{}.calls.tsv".format(args.dataset, domain)),
                c,
                keys,
            )
        for rank in RANKS:
            m = score_calls(parent, rank_map, rank, t, c)
            rec = {
                "scheme": args.scheme,
                "dataset": args.dataset,
                "domain": domain,
                "rank": rank,
                "n_truth": len(t),
                "n_calls": len(c),
            }
            rec.update(m)
            rows.append(rec)
            print(
                "{:<12} {:<12} {:<8} {}".format(domain, rank, args.scheme or "-", format_row(m))
            )

    fields = [
        "scheme",
        "minimum_hit_groups",
        "dataset",
        "domain",
        "rank",
        "n_truth",
        "n_calls",
        "tp",
        "fp",
        "fn",
        "vp",
        "no",
        "n",
        "sensitivity",
        "ppv",
        "f1",
    ]
    mhg = "" if args.minimum_hit_groups is None else args.minimum_hit_groups
    with open(args.out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for rec in rows:
            rec["minimum_hit_groups"] = mhg
            w.writerow(rec)
    print("  wrote", args.out_tsv)

    exp = expected_s1_row(args.expected_s1)
    if exp:
        print("  Table S1 Kraken 2 genus (expected):")
        print(
            "    prokaryote  sens={} ppv={} f1={}".format(
                exp.get("prok_genus_sensitivity"),
                exp.get("prok_genus_ppv"),
                exp.get("prok_genus_f1"),
            )
        )
        print(
            "    virus       sens={} ppv={} f1={}".format(
                exp.get("viral_genus_sensitivity"),
                exp.get("viral_genus_ppv"),
                exp.get("viral_genus_f1"),
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
