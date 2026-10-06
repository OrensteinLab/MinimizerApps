#!/usr/bin/env python3
"""Profile-level scoring of a mock community from a Kraken 2 report.

Compares the observed taxon profile (clade read counts at genus / species rank
in the standard `kraken2 --report`) with the published composition of the mock
(tables/mock_zymo_d6300_expected.csv). No Kraken database is needed.

Per scope x rank x detection threshold it reports:
  precision / recall / F1   of the detected taxa vs the expected taxa
  l1, bray_curtis           between expected and observed relative abundances
                            (expected renormalised to 1 over the scope; observed
                            = clade reads of the taxa at the rank / sum over all
                            taxa at the rank, so spurious taxa dilute it)
  classified_frac           reads classified at any level / reads in scope
  on_target_frac            reads in expected taxa at the rank / reads in scope
  off_target_frac           reads in other taxa at the rank / reads in scope
  above_rank_frac           classified, but not at or below the rank
  unclassified_frac         reads in scope not classified
  (on + off + above + unclassified = 1)

Scopes: `all` (every expected taxon) and `bacteria` (expected restricted to
bacteria and renormalised; reads of the expected non-bacterial taxa are removed
from the denominator, which matters for libraries without fungi).

A taxon is detected when its clade reads / reads in scope >= min_frac
(min_frac 0 = any read). The expected genus of each species comes from, in
order: --nodes lineage, the report lineage (if the species is in the report),
the genus_taxid column of the composition table.

Usage:
  score_mock.py --report r.txt --composition mock_zymo_d6300_expected.csv \\
      --out-tsv metrics.tsv [--nodes nodes.dmp] [--scheme xor] \\
      [--abundance gdna_pct] [--min-fracs 0,1e-4,1e-3]
"""

from __future__ import print_function

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

RANK_CODE = {"genus": "G", "species": "S"}
SCOPES = ("all", "bacteria")
FIELDS = [
    "scheme",
    "minimum_hit_groups",
    "dataset",
    "scope",
    "rank",
    "min_frac",
    "n_reads",
    "n_expected",
    "n_detected",
    "tp_taxa",
    "fp_taxa",
    "fn_taxa",
    "precision",
    "recall",
    "f1",
    "l1",
    "bray_curtis",
    "classified_frac",
    "on_target_frac",
    "off_target_frac",
    "above_rank_frac",
    "unclassified_frac",
    "top_off_target",
]


class Node(object):
    __slots__ = ("taxid", "name", "code", "clade", "taxon", "depth", "parent")

    def __init__(self, taxid, name, code, clade, taxon, depth, parent):
        self.taxid = taxid
        self.name = name
        self.code = code
        self.clade = clade
        self.taxon = taxon
        self.depth = depth
        self.parent = parent


def parse_report(path):
    """Parse a Kraken 2 report -> (n_unclassified, nodes by taxid).

    Columns: pct, clade reads, taxon reads, rank code, taxid, name (indented
    two spaces per level). Extra trailing columns (minimizer data) are ignored.
    """
    unclassified = 0
    nodes = {}
    stack = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                clade = int(parts[1])
                taxon = int(parts[2])
                taxid = int(parts[4])
            except ValueError:
                continue
            raw = parts[5]
            depth = (len(raw) - len(raw.lstrip(" "))) // 2
            if taxid == 0 and parts[3] == "U":
                unclassified = clade
                continue
            while stack and stack[-1].depth >= depth:
                stack.pop()
            parent = stack[-1] if stack else None
            node = Node(taxid, raw.strip(), parts[3], clade, taxon, depth, parent)
            nodes[taxid] = node
            stack.append(node)
    return unclassified, nodes


def load_composition(path, abundance):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r["taxid"] = int(r["taxid"])
            r["genus_taxid"] = int(r["genus_taxid"]) if r.get("genus_taxid") else 0
            r["abundance"] = float(r[abundance])
            rows.append(r)
    if not rows:
        raise SystemExit("empty composition table " + path)
    return rows


def genus_via_report(nodes, taxid):
    node = nodes.get(taxid)
    while node is not None:
        if node.code == "G":
            return node.taxid
        node = node.parent
    return 0


def expected_taxa(comp, rank, nodes, tax):
    """taxid at rank -> abundance (summed if two species share a genus)."""
    out = {}
    for r in comp:
        if rank == "species":
            t = r["taxid"]
        else:
            t = 0
            if tax is not None:
                from evaluate_calls import elevate_to_rank
                t = elevate_to_rank(tax[0], tax[1], r["taxid"], "genus")
            if not t:
                t = genus_via_report(nodes, r["taxid"])
            if not t:
                t = r["genus_taxid"]
        out[t] = out.get(t, 0.0) + r["abundance"]
    return out


def clade_at_rank(nodes, rank):
    code = RANK_CODE[rank]
    return {t: n for t, n in nodes.items() if n.code == code}


def fmt_top(taxa, total, k=3):
    top = sorted(taxa, key=lambda x: -x[1])[:k]
    return "; ".join("{} {:.4f}".format(n, c / float(total)) for n, c in top if total)


def score_rank(scheme, dataset, scope, rank, min_fracs, comp, n_unclassified,
               nodes, tax):
    if scope == "bacteria":
        comp_s = [r for r in comp if r["domain"] == "bacteria"]
        drop = [r["taxid"] for r in comp if r["domain"] != "bacteria"]
    else:
        comp_s, drop = comp, []
    exp = expected_taxa(comp_s, rank, nodes, tax)
    exp_total = sum(exp.values())
    exp = {t: a / exp_total for t, a in exp.items()}

    total_root = nodes[1].clade if 1 in nodes else 0
    n_all = n_unclassified + total_root
    # Remove the reads of expected taxa outside the scope from the denominator.
    removed = 0
    if drop:
        drop_rows = [r for r in comp if r["taxid"] in drop]
        drop_genera = set(expected_taxa(drop_rows, "genus", nodes, tax))
        removed = sum(nodes[g].clade for g in drop_genera if g in nodes)
        removed += sum(nodes[r["taxid"]].clade for r in drop_rows
                       if r["taxid"] in nodes
                       and not (set(expected_taxa([r], "genus", nodes, tax)) & set(nodes)))
    n_reads = n_all - removed
    classified = total_root - removed
    at_rank = clade_at_rank(nodes, rank)
    dropped_at_rank = set(expected_taxa([r for r in comp if r["taxid"] in drop],
                                        rank, nodes, tax)) if drop else set()
    obs = {t: n.clade for t, n in at_rank.items() if t not in dropped_at_rank}
    obs_total = sum(obs.values())
    rows = []
    for min_frac in min_fracs:
        thr = max(1, int(min_frac * n_reads + 0.999999)) if min_frac > 0 else 1
        detected = {t for t, c in obs.items() if c >= thr}
        tp = len(detected & set(exp))
        fp = len(detected - set(exp))
        fn = len(set(exp) - detected)
        prec = tp / float(tp + fp) if (tp + fp) else 0.0
        rec = tp / float(len(exp)) if exp else 0.0
        f1 = 0.0 if prec * rec == 0 else 2 * prec * rec / (prec + rec)
        union = set(exp) | set(obs)
        if obs_total:
            l1 = sum(abs(exp.get(t, 0.0) - obs.get(t, 0) / float(obs_total)) for t in union)
        else:
            l1 = 2.0
        on = sum(c for t, c in obs.items() if t in exp)
        off_taxa = [(at_rank[t].name, c) for t, c in obs.items() if t not in exp]
        off = sum(c for _n, c in off_taxa)
        nr = float(n_reads) if n_reads else 1.0
        rows.append(
            {
                "scheme": scheme,
                "dataset": dataset,
                "scope": scope,
                "rank": rank,
                "min_frac": min_frac,
                "n_reads": n_reads,
                "n_expected": len(exp),
                "n_detected": len(detected),
                "tp_taxa": tp,
                "fp_taxa": fp,
                "fn_taxa": fn,
                "precision": prec,
                "recall": rec,
                "f1": f1,
                "l1": l1,
                "bray_curtis": l1 / 2.0,
                "classified_frac": classified / nr,
                "on_target_frac": on / nr,
                "off_target_frac": off / nr,
                "above_rank_frac": (classified - obs_total) / nr,
                "unclassified_frac": n_unclassified / nr,
                "top_off_target": fmt_top(off_taxa, nr),
            }
        )
    return rows


def score_report(report, composition, scheme="", dataset="mock", nodes_dmp="",
                 abundance="gdna_pct", min_fracs=(0.0, 1e-4, 1e-3)):
    comp = load_composition(composition, abundance)
    n_unclassified, nodes = parse_report(report)
    tax = None
    if nodes_dmp:
        from evaluate_calls import load_taxonomy
        tax = load_taxonomy(nodes_dmp)
    rows = []
    for scope in SCOPES:
        for rank in ("genus", "species"):
            rows.extend(score_rank(scheme, dataset, scope, rank, min_fracs, comp,
                                   n_unclassified, nodes, tax))
    return rows


def parse_args(argv):
    p = argparse.ArgumentParser(description="Score a mock-community Kraken 2 report.")
    p.add_argument("--report", required=True, help="kraken2 --report file")
    p.add_argument("--composition", required=True, help="expected composition CSV")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--nodes", default="", help="nodes.dmp (derive expected genera)")
    p.add_argument("--scheme", default="")
    p.add_argument("--dataset", default="mock")
    p.add_argument("--minimum-hit-groups", type=int, default=None,
                   help="classify --minimum-hit-groups used for this output; "
                        "recorded in the minimum_hit_groups column")
    p.add_argument("--abundance", default="gdna_pct",
                   help="expected column: gdna_pct (reads, default), genome_copy_pct, cell_pct")
    p.add_argument("--min-fracs", default="0,1e-4,1e-3",
                   help="comma-separated detection thresholds (fraction of reads)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    fracs = tuple(float(x) for x in args.min_fracs.split(",") if x.strip())
    rows = score_report(args.report, args.composition, args.scheme, args.dataset,
                        args.nodes, args.abundance, fracs)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r["minimum_hit_groups"] = (
                "" if args.minimum_hit_groups is None else args.minimum_hit_groups)
            w.writerow(r)
    for r in rows:
        if r["scope"] == "all" and r["min_frac"] == 1e-4:
            print("{:<8} {:<8} P={:.3f} R={:.3f} F1={:.3f} L1={:.4f} BC={:.4f} classified={:.4f} off={:.4f}".format(
                r["scheme"] or "-", r["rank"], r["precision"], r["recall"], r["f1"],
                r["l1"], r["bray_curtis"], r["classified_frac"], r["off_target_frac"]))
    print("  wrote", args.out_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
