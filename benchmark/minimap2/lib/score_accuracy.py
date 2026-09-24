#!/usr/bin/env python3
"""Score a PAF file against truth.tsv (see make_dataset.py): for each read,
find its primary alignment (tp:A:P) and check whether it overlaps the read's
true origin by at least `--min-overlap-frac` of the true span. Reports
overall and repeat-region-only accuracy (the paper's own headline
comparison) plus mean per-base identity (from PAF's de:f: tag, when mapped).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def load_truth(path: Path) -> dict[str, dict]:
    truth = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = {
                "start": int(row["true_start"]),
                "end": int(row["true_end"]),
                "strand": row["true_strand"],
                "in_repeat": row["in_repeat"] == "1",
            }
    return truth


def parse_paf_primary(path: Path) -> dict[str, dict]:
    """Returns {read_id: {tstart, tend, strand, de}} for each read's primary
    alignment (tp:A:P), or the first alignment seen if no tp tag is present."""
    best: dict[str, dict] = {}
    with path.open() as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qname, _qlen, _qs, _qe, strand, _tname, _tlen, tstart, tend = cols[:9]
            tags = {}
            for tag in cols[12:]:
                parts = tag.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            is_primary = tags.get("tp") == "P"
            if qname in best and not is_primary:
                continue
            if qname in best and best[qname].get("_is_primary"):
                continue
            # de:f: (gap-compressed per-base divergence) needs -c/--cs; without
            # it minimap2 still emits dv:f: (approximate divergence from
            # chaining/minimizer matches alone) -- prefer de, fall back to dv.
            div = tags.get("de", tags.get("dv"))
            best[qname] = {
                "tstart": int(tstart),
                "tend": int(tend),
                "strand": strand,
                "de": float(div) if div is not None else None,
                "_is_primary": is_primary,
            }
    return best


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paf", required=True, type=Path)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--min-overlap-frac", type=float, default=0.5)
    ap.add_argument("--json-out", type=Path)
    args = ap.parse_args(argv)

    truth = load_truth(args.truth)
    aligned = parse_paf_primary(args.paf)

    n_total = len(truth)
    n_repeat = sum(1 for t in truth.values() if t["in_repeat"])
    n_mapped = n_correct = n_correct_repeat = 0
    identities = []
    for rid, t in truth.items():
        a = aligned.get(rid)
        if a is None:
            continue
        n_mapped += 1
        if a["de"] is not None:
            identities.append(1.0 - a["de"])
        ov = overlap_len(a["tstart"], a["tend"], t["start"], t["end"])
        true_span = t["end"] - t["start"]
        correct = true_span > 0 and ov / true_span >= args.min_overlap_frac
        if correct:
            n_correct += 1
            if t["in_repeat"]:
                n_correct_repeat += 1

    result = {
        "n_reads": n_total,
        "n_reads_in_repeat": n_repeat,
        "n_mapped": n_mapped,
        "n_correct": n_correct,
        "n_correct_in_repeat": n_correct_repeat,
        "frac_mapped": round(n_mapped / n_total, 4) if n_total else None,
        "frac_correct": round(n_correct / n_total, 4) if n_total else None,
        "frac_correct_in_repeat": round(n_correct_repeat / n_repeat, 4) if n_repeat else None,
        "mean_identity": round(sum(identities) / len(identities), 4) if identities else None,
    }
    print(json.dumps(result, indent=2))
    if args.json_out:
        args.json_out.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
