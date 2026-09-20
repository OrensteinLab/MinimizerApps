#!/usr/bin/env python3
"""Score FDA-ARGOS Kraken calls vs SRA-assigned taxa (Wood et al. Table S5).

Does not report sensitivity or PPV: real sequencing fragments may be
contaminated, low-quality, or distant from the reference, so the true origin
of each fragment is not known with certainty.

At genus rank, every experiment (SRA run) is scored as four percentages of
fragments whose true lineage has a genus node:
  (a) concordant    call is the SRA-assigned genus or a descendant
  (b) discordant    call is off-lineage
  (c) ancestor      call is an ancestor of the SRA-assigned genus
  (d) unclassified  no classification (taxid 0)

Paper Table S5 averages are the unweighted mean of the 25 experiments.
A pooled-fragment row is also written for scheme comparison.

Usage:
  score_argos.py --nodes nodes.dmp --truth truth.tsv --kraken out.kraken \\
      --out-tsv metrics.tsv [--expected table_s5_kraken2_expected.csv]
"""

from __future__ import print_function

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from evaluate_calls import load_taxonomy, score_calls
from score_kraken import load_kraken_calls, load_truth

ARGOS_RUN_RE = re.compile(r"^ARGOS\.([A-Z0-9]+)\.", re.I)
RANK = "genus"
FIELDS = [
    "scheme",
    "dataset",
    "rank",
    "scope",
    "sra_run",
    "species",
    "n",
    "n_scored",
    "tp",
    "fp",
    "fn",
    "vp",
    "no",
    "concordant_frac",
    "discordant_frac",
    "ancestor_frac",
    "unclassified_frac",
]
FOUR = (
    "concordant_frac",
    "discordant_frac",
    "ancestor_frac",
    "unclassified_frac",
)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Score ARGOS genus concordance (Table S5; no sens/PPV)."
    )
    p.add_argument("--nodes", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--kraken", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--scheme", default="")
    p.add_argument("--expected", default="", help="table_s5_kraken2_expected.csv")
    return p.parse_args(argv)


def load_expected(path):
    if not path or not os.path.isfile(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def run_of(seqid):
    m = ARGOS_RUN_RE.match(seqid)
    return m.group(1) if m else ""


def group_by_run(truth):
    groups = defaultdict(dict)
    for seqid, taxid in truth.items():
        run = run_of(seqid)
        groups[run][seqid] = taxid
    return groups


def expected_index(rows):
    by_run = {}
    for row in rows:
        run = (row.get("sra_run") or "").strip()
        if run:
            by_run[run] = row
    return by_run


def mean_key(rows, key):
    vals = []
    for r in rows:
        v = r.get(key)
        if v in ("", None):
            continue
        vals.append(float(v))
    return (sum(vals) / len(vals)) if vals else None


def metrics_row(scheme, m, scope, sra_run, species):
    rec = {
        "scheme": scheme,
        "dataset": "fda-argos",
        "rank": RANK,
        "scope": scope,
        "sra_run": sra_run,
        "species": species,
        "n": m["n"],
        "n_scored": m.get("n_scored", m["tp"] + m["fp"] + m["fn"] + m["vp"]),
        "tp": m["tp"],
        "fp": m["fp"],
        "fn": m["fn"],
        "vp": m["vp"],
        "no": m["no"],
    }
    for k in FOUR:
        rec[k] = m[k]
    return rec


def unweighted_mean_row(scheme, exp_rows):
    rec = {
        "scheme": scheme,
        "dataset": "fda-argos",
        "rank": RANK,
        "scope": "unweighted_mean",
        "sra_run": "",
        "species": "",
        "n": sum(int(r["n"]) for r in exp_rows),
        "n_scored": sum(int(r["n_scored"]) for r in exp_rows),
        "tp": "",
        "fp": "",
        "fn": "",
        "vp": "",
        "no": "",
    }
    for k in FOUR:
        rec[k] = mean_key(exp_rows, k) or 0.0
    return rec


def fmt_four(rec):
    return "conc={:.4f} disc={:.4f} anc={:.4f} unclass={:.4f}".format(
        float(rec["concordant_frac"]),
        float(rec["discordant_frac"]),
        float(rec["ancestor_frac"]),
        float(rec["unclassified_frac"]),
    )


def calls_for(truth_subset, calls):
    return {k: calls[k] for k in truth_subset if k in calls}


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    parent, rank_map = load_taxonomy(args.nodes)
    truth, _domains = load_truth(args.truth)
    calls = load_kraken_calls(args.kraken)
    exp_meta = expected_index(load_expected(args.expected))

    grouped = group_by_run(truth)
    if not grouped:
        print("Error: no truth records", file=sys.stderr)
        return 1
    unlabeled = grouped.pop("", None)
    if unlabeled:
        print(
            "  warning: {} truth seqids lack ARGOS.<SRR>. prefix".format(len(unlabeled)),
            file=sys.stderr,
        )

    exp_rows = []
    for run in sorted(grouped):
        t = grouped[run]
        m = score_calls(parent, rank_map, RANK, t, calls_for(t, calls))
        species = (exp_meta.get(run) or {}).get("species", "")
        rec = metrics_row(args.scheme, m, "experiment", run, species)
        exp_rows.append(rec)

    pooled = score_calls(parent, rank_map, RANK, dict(truth), dict(calls))
    rows = list(exp_rows)
    rows.append(unweighted_mean_row(args.scheme, exp_rows))
    rows.append(metrics_row(args.scheme, pooled, "pooled", "", ""))

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    with open(args.out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for rec in rows:
            w.writerow(rec)

    mean_rec = rows[-2]
    pooled_rec = rows[-1]
    print("  genus rank (no sensitivity/PPV; Wood et al. Table S5)")
    print(
        "  unweighted mean of {} experiments: {}".format(
            len(exp_rows), fmt_four(mean_rec)
        )
    )
    print("  pooled fragments:               {}".format(fmt_four(pooled_rec)))
    print("  wrote", args.out_tsv)

    if exp_meta:
        paper_mean = {
            "concordant_frac": mean_key(list(exp_meta.values()), "genus_concordant") or 0.0,
            "discordant_frac": mean_key(list(exp_meta.values()), "genus_discordant") or 0.0,
            "ancestor_frac": mean_key(list(exp_meta.values()), "genus_ancestor") or 0.0,
            "unclassified_frac": mean_key(list(exp_meta.values()), "genus_unclassified")
            or 0.0,
        }
        print("  Table S5 Kraken 2 genus averages (25 experiments, paper):")
        print("    {}".format(fmt_four(paper_mean)))
        print(
            "  Note: paper averages are unweighted per-experiment means; "
            "slices of each SRR will not match Table S5 counts."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
