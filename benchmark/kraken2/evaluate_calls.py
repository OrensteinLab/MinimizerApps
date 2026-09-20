#!/usr/bin/env python3
"""Score per-fragment calls vs truth at a taxonomic rank (Wood et al. 2019).

Matches DerrickWood/kraken2-experiment-code Utilities/evaluate_calls.cc:

  TP  call is the true taxon at RANK, or a descendant of it
  VP  call is an ancestor of the true RANK taxon (vague / above rank)
  FP  call is off-lineage
  FN  unclassified (taxid 0) or missing from calls
  NO  true lineage has no node at RANK

Sensitivity = TP / (TP + VP + FN + FP)
PPV         = TP / (TP + FP)
F1          = harmonic mean of sensitivity and PPV

Usage:
  evaluate_calls.py nodes.dmp genus truth.tsv calls.tsv
"""

from __future__ import print_function

import argparse
import sys


def load_taxonomy(nodes_path):
    parent = {}
    rank = {}
    with open(nodes_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t|\t")
            if len(parts) < 3:
                continue
            node_id = int(parts[0])
            parent_id = int(parts[1])
            rnk = parts[2].rstrip("\t|").strip()
            if node_id == 1:
                parent_id = 0
            parent[node_id] = parent_id
            rank[node_id] = rnk
    return parent, rank


def load_two_col(path):
    """Map first whitespace field → integer taxid."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.replace(",", " ").split()
            if len(fields) < 2:
                continue
            name, tax = fields[0], fields[1]
            try:
                out[name] = int(tax)
            except ValueError:
                continue
    return out


def is_ancestor(parent, a, b):
    """True if a is b or an ancestor of b."""
    while b > 0:
        if a == b:
            return True
        b = parent.get(b, 0)
    return False


def elevate_to_rank(parent, rank_map, taxid, rank):
    while taxid > 0:
        if rank_map.get(taxid) == rank:
            return taxid
        taxid = parent.get(taxid, 0)
    return 0


def concordance_fracs(tp, fp, fn, vp):
    """Wood et al. Table S5 genus quantities (fraction of fragments).

    Denominator is fragments whose true lineage has a node at the scored rank
    (excludes NO). The four values sum to 1:
      concordant   TP  labeled RANK taxon or a descendant
      discordant   FP  off-lineage
      ancestor     VP  call is an ancestor of the labeled RANK taxon
      unclassified FN  taxid 0 or missing call
    """
    n = tp + fp + fn + vp
    if not n:
        z = 0.0
        return {
            "n_scored": 0,
            "concordant_frac": z,
            "discordant_frac": z,
            "ancestor_frac": z,
            "unclassified_frac": z,
        }
    d = float(n)
    return {
        "n_scored": n,
        "concordant_frac": tp / d,
        "discordant_frac": fp / d,
        "ancestor_frac": vp / d,
        "unclassified_frac": fn / d,
    }


def score_calls(parent, rank_map, rank, truth, calls):
    truth = dict(truth)
    tp = fp = fn = vp = no = 0
    for name, taxon in calls.items():
        if name not in truth:
            continue
        correct = truth.pop(name)
        if taxon == 0:
            fn += 1
            continue
        if correct not in rank_map:
            correct = 0
        correct = elevate_to_rank(parent, rank_map, correct, rank)
        if not correct:
            no += 1
            continue
        if is_ancestor(parent, correct, taxon):
            tp += 1
        elif is_ancestor(parent, taxon, correct):
            vp += 1
        else:
            fp += 1
    fn += len(truth)
    n = tp + fp + fn + vp + no
    denom_s = tp + fn + fp + vp
    sens = (tp / float(denom_s)) if denom_s else 0.0
    prec = (tp / float(tp + fp)) if (tp + fp) else 0.0
    f1 = 0.0 if (sens * prec == 0) else (2 * sens * prec / (sens + prec))
    out = {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "vp": vp,
        "no": no,
        "n": n,
        "sensitivity": sens,
        "ppv": prec,
        "f1": f1,
    }
    out.update(concordance_fracs(tp, fp, fn, vp))
    return out


def format_row(m):
    return "{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}".format(
        m["tp"],
        m["fp"],
        m["fn"],
        m["vp"],
        m["no"],
        m["n"],
        m["sensitivity"],
        m["ppv"],
        m["f1"],
    )


def parse_args(argv):
    p = argparse.ArgumentParser(description="Score taxonomic calls at a rank.")
    p.add_argument("nodes", help="NCBI nodes.dmp")
    p.add_argument("rank", help="Rank name, e.g. genus or species")
    p.add_argument("truth", help="TSV: seqid  true_taxid")
    p.add_argument("calls", help="TSV: seqid  predicted_taxid (0 = unclassified)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    parent, rank_map = load_taxonomy(args.nodes)
    truth = load_two_col(args.truth)
    calls = load_two_col(args.calls)
    m = score_calls(parent, rank_map, args.rank, truth, calls)
    print(format_row(m))
    return 0


if __name__ == "__main__":
    sys.exit(main())
