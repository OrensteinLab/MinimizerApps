#!/usr/bin/env python3
"""Score Kraken 2 calls on CAMI-style reads, aware of what the library holds.

CAMI genomes are mostly NOT in the classification library, so a plain
Wood et al. sensitivity (which counts every unclassified read as a miss) would
mix "scheme loses reads" with "the library never had this taxon". For each
rank (genus, species) every read falls in exactly one bucket:

  no         the true taxid has no node at the rank (or is not in nodes.dmp),
             e.g. new_species reads are labelled with their genus by CAMI
  excluded   the true taxon at the rank has no reference sequence in the
             library (--library-dir); reported, never scored
  scored     true taxon at the rank is in the library: TP / VP / FP / FN as in
             evaluate_calls.py (Wood et al. 2019)

Sensitivity = TP/(TP+VP+FP+FN), PPV = TP/(TP+FP), F1 over *scored* reads only.
Excluded reads are further split by what the classifier did with them:
ex_unclassified, ex_consistent (call is on the true lineage: the true taxon, an
ancestor such as the family, or a descendant), ex_off (call is off-lineage = a spurious assignment).

Groups: all; domain (bacteria / archaea / virus / eukaryote / other, from the
superkingdom of the true taxid); novelty:<category> from truth_extra.tsv (CAMI
metadata: known_strain, new_strain, new_species, ...).

Without --kraken only the library coverage of the truth is printed
(how many reads / taxa would be excluded for that library).

Usage:
  score_cami.py --nodes nodes.dmp --truth truth.tsv --library-dir db/standard_p10 \\
      --kraken out.kraken --scheme Random --out-tsv metrics.tsv
  score_cami.py --nodes nodes.dmp --truth truth.tsv --library-dir db/strain_exclusion
"""

from __future__ import print_function

import argparse
import csv
import glob
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from evaluate_calls import elevate_to_rank, is_ancestor, load_taxonomy  # noqa: E402
from score_kraken import load_kraken_calls, load_truth  # noqa: E402

RANKS = ("genus", "species")
DOMAIN_OF_SUPERKINGDOM = {
    2: "bacteria",
    2157: "archaea",
    10239: "virus",
    2759: "eukaryote",
    12884: "virus",  # viroids
}
NOVELTY_SKIP = ("", "NA", "unknown")
FIELDS = [
    "scheme",
    "dataset",
    "group",
    "rank",
    "n_truth",
    "n_calls",
    "no",
    "excluded",
    "excluded_frac",
    "ex_unclassified",
    "ex_consistent",
    "ex_off",
    "n_scored",
    "tp",
    "fp",
    "fn",
    "vp",
    "sensitivity",
    "ppv",
    "f1",
]


def domain_of(parent, taxid):
    """Domain label from the superkingdom ancestor of taxid ('other' if none)."""
    while taxid > 0:
        if taxid in DOMAIN_OF_SUPERKINGDOM:
            return DOMAIN_OF_SUPERKINGDOM[taxid]
        taxid = parent.get(taxid, 0)
    return "other"


def library_taxids(db_dir):
    """Taxids with sequence in a Kraken library dir (prelim_map*.txt, col 3).

    Works for db/standard*/library/<kingdom>/prelim_map.txt and the strain
    exclusion library/added/prelim_map_<md5>.txt. Falls back to the DB's
    seqid2taxid.map (second column) if no prelim maps exist.
    """
    taxids = set()
    lib = os.path.join(db_dir, "library")
    maps = glob.glob(os.path.join(lib, "**", "prelim_map*.txt"), recursive=True)
    for path in maps:
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[0] == "TAXID":
                    try:
                        taxids.add(int(parts[2]))
                    except ValueError:
                        pass
    if taxids:
        return taxids
    s2t = os.path.join(db_dir, "seqid2taxid.map")
    if os.path.isfile(s2t):
        with open(s2t) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        taxids.add(int(parts[1]))
                    except ValueError:
                        pass
    return taxids


class RankIndex(object):
    """Memoised elevate-to-rank + the library's taxa at each rank."""

    def __init__(self, parent, rank_map, lib_taxids):
        self.parent = parent
        self.rank_map = rank_map
        self._memo = {}
        self.has_library = lib_taxids is not None
        self.lib = {r: set() for r in RANKS}
        for t in lib_taxids or ():
            if t not in rank_map:
                continue
            for r in RANKS:
                e = self.elevate(t, r)
                if e:
                    self.lib[r].add(e)

    def elevate(self, taxid, rank):
        key = (taxid, rank)
        if key not in self._memo:
            self._memo[key] = elevate_to_rank(self.parent, self.rank_map, taxid, rank)
        return self._memo[key]

    def in_library(self, taxon_at_rank, rank):
        return (not self.has_library) or taxon_at_rank in self.lib[rank]


def score_group(idx, rank, truth, calls):
    """Bucket and score the reads of one group at one rank."""
    c = Counter()
    for seqid, true_tax in truth.items():
        call = calls.get(seqid, 0)
        if true_tax not in idx.rank_map:
            c["no"] += 1
            continue
        at_rank = idx.elevate(true_tax, rank)
        if not at_rank:
            c["no"] += 1
            continue
        if not idx.in_library(at_rank, rank):
            c["excluded"] += 1
            if call == 0:
                c["ex_unclassified"] += 1
            elif (is_ancestor(idx.parent, call, true_tax)
                  or is_ancestor(idx.parent, true_tax, call)):
                c["ex_consistent"] += 1
            else:
                c["ex_off"] += 1
            continue
        if call == 0:
            c["fn"] += 1
        elif is_ancestor(idx.parent, at_rank, call):
            c["tp"] += 1
        elif is_ancestor(idx.parent, call, at_rank):
            c["vp"] += 1
        else:
            c["fp"] += 1
    tp, fp, fn, vp = c["tp"], c["fp"], c["fn"], c["vp"]
    n_scored = tp + fp + fn + vp
    sens = tp / float(n_scored) if n_scored else 0.0
    ppv = tp / float(tp + fp) if (tp + fp) else 0.0
    f1 = 0.0 if sens * ppv == 0 else 2 * sens * ppv / (sens + ppv)
    n = len(truth)
    return {
        "n_truth": n,
        "n_calls": sum(1 for s in truth if calls.get(s, 0) != 0),
        "no": c["no"],
        "excluded": c["excluded"],
        "excluded_frac": c["excluded"] / float(n) if n else 0.0,
        "ex_unclassified": c["ex_unclassified"],
        "ex_consistent": c["ex_consistent"],
        "ex_off": c["ex_off"],
        "n_scored": n_scored,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "vp": vp,
        "sensitivity": sens,
        "ppv": ppv,
        "f1": f1,
    }


def load_novelty(path):
    """seqid -> novelty category from truth_extra.tsv (seqid taxid genome novelty)."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                out[parts[0]] = parts[3]
    return out


def make_groups(parent, truth, novelty):
    groups = defaultdict(dict)
    for seqid, taxid in truth.items():
        groups["all"][seqid] = taxid
        groups[domain_of(parent, taxid)][seqid] = taxid
        nov = novelty.get(seqid, "")
        if nov not in NOVELTY_SKIP:
            groups["novelty:" + nov][seqid] = taxid
    return groups


def group_order(name):
    if name == "all":
        return (0, name)
    if name.startswith("novelty:"):
        return (2, name)
    return (1, name)


def coverage_report(idx, truth):
    """Print library coverage of the truth (reads and distinct taxa) per rank."""
    print("library coverage of the truth ({} reads):".format(len(truth)))
    print("{:<8} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}".format(
        "rank", "no_node", "in_lib", "excluded", "excl_%", "taxa", "taxa_excl"))
    for rank in RANKS:
        no = inl = exc = 0
        taxa = {}
        for t in truth.values():
            r = idx.elevate(t, rank) if t in idx.rank_map else 0
            if not r:
                no += 1
                continue
            taxa[r] = idx.in_library(r, rank)
            if taxa[r]:
                inl += 1
            else:
                exc += 1
        denom = float(inl + exc) or 1.0
        print("{:<8} {:>9} {:>9} {:>9} {:>8.1f}% {:>9} {:>9}".format(
            rank, no, inl, exc, 100.0 * exc / denom, len(taxa),
            sum(1 for v in taxa.values() if not v)))


def parse_args(argv):
    p = argparse.ArgumentParser(description="Score CAMI reads (library-aware).")
    p.add_argument("--nodes", required=True, help="nodes.dmp of the DB taxonomy")
    p.add_argument("--truth", required=True, help="truth.tsv (#seqid taxid domain)")
    p.add_argument("--kraken", default="", help="Kraken 2 --output file (omit: coverage only)")
    p.add_argument("--library-dir", default="",
                   help="shared DB dir (db/<id>) whose library defines 'excluded'")
    p.add_argument("--extra", default="",
                   help="truth_extra.tsv (default: next to truth.tsv if present)")
    p.add_argument("--out-tsv", default="")
    p.add_argument("--scheme", default="")
    p.add_argument("--dataset", default="cami")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    parent, rank_map = load_taxonomy(args.nodes)
    truth, _domains = load_truth(args.truth)
    lib = library_taxids(args.library_dir) if args.library_dir else None
    if args.library_dir and not lib:
        sys.exit("no library taxids found under " + args.library_dir)
    idx = RankIndex(parent, rank_map, lib)
    extra = args.extra or os.path.join(os.path.dirname(os.path.abspath(args.truth)),
                                       "truth_extra.tsv")
    groups = make_groups(parent, truth, load_novelty(extra))

    if not args.kraken:
        coverage_report(idx, truth)
        return 0
    if not args.out_tsv:
        sys.exit("--out-tsv is required with --kraken")

    calls = load_kraken_calls(args.kraken)
    rows = []
    for name in sorted(groups, key=group_order):
        for rank in RANKS:
            rec = {"scheme": args.scheme, "dataset": args.dataset, "group": name, "rank": rank}
            rec.update(score_group(idx, rank, groups[name], calls))
            rows.append(rec)
            if name == "all":
                print("{:<10} {:<8} {:<10} scored={:<8} excl={:<8} sens={:.4f} ppv={:.4f} f1={:.4f}".format(
                    name, rank, args.scheme or "-", rec["n_scored"], rec["excluded"],
                    rec["sensitivity"], rec["ppv"], rec["f1"]))
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for rec in rows:
            w.writerow(rec)
    print("  wrote", args.out_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
