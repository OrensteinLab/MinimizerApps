#!/usr/bin/env python3
"""Score a PAF file against truth.tsv (see make_dataset.py/pbsim_reads.py):
for each read, find its primary alignment (tp:A:P) and check whether it
overlaps the read's true origin, at two thresholds:

  - "correct_lit": overlap >= 10% of the UNION of the true and reported
    intervals -- the exact criterion both Winnowmap (Jain et al. 2020, "We
    followed a similar criteria previously used by Li (2018)") and minimap2
    (Li 2018, "the overlap length is >=10% of the true interval length" --
    Winnowmap's own restatement uses the union, which is what's implemented
    here) use, so this number is comparable to published minimap2/Winnowmap
    accuracy figures.
  - "correct_strict": overlap >= 50% of the true span -- this repo's
    original, stricter threshold, kept for continuity with earlier smoke-scale
    runs (see benchmark/minimap2/README.md's first sweep table).

Reports both overall and repeat-region-only (the papers' own headline
comparison -- Winnowmap's Table 3 restricts to reads sampled from the
centromeric repeat array) plus mean per-base identity (from PAF's de:f:/dv:f:
tag, when mapped).
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
                "contig": row["true_contig"],
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
            qname, _qlen, _qs, _qe, strand, tname, _tlen, tstart, tend = cols[:9]
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
                "tname": tname,
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
    ap.add_argument("--json-out", type=Path)
    args = ap.parse_args(argv)

    truth = load_truth(args.truth)
    aligned = parse_paf_primary(args.paf)

    n_total = len(truth)
    n_repeat = sum(1 for t in truth.values() if t["in_repeat"])
    n_mapped = n_correct_lit = n_correct_lit_repeat = 0
    n_correct_strict = n_correct_strict_repeat = 0
    identities = []
    for rid, t in truth.items():
        a = aligned.get(rid)
        if a is None:
            continue
        n_mapped += 1
        if a["de"] is not None:
            identities.append(1.0 - a["de"])
        # A read whose primary alignment lands on a *different* contig than
        # its true origin is wrong regardless of coordinate overlap -- e.g.
        # grch38_wg has 195 sequences, so two contigs can easily share a
        # coordinate range. Single-contig datasets (smoke, chm13_chrX) never
        # hit this (there's only one contig to land on), which is exactly
        # why it went unnoticed until a real multi-contig genome was added.
        same_contig = a["tname"] == t["contig"]
        ov = overlap_len(a["tstart"], a["tend"], t["start"], t["end"]) if same_contig else 0
        true_span = t["end"] - t["start"]
        aligned_span = a["tend"] - a["tstart"]
        union_span = true_span + aligned_span - ov

        correct_lit = same_contig and union_span > 0 and ov / union_span >= 0.10
        correct_strict = same_contig and true_span > 0 and ov / true_span >= 0.50

        if correct_lit:
            n_correct_lit += 1
            if t["in_repeat"]:
                n_correct_lit_repeat += 1
        if correct_strict:
            n_correct_strict += 1
            if t["in_repeat"]:
                n_correct_strict_repeat += 1

    def frac(n, d):
        return round(n / d, 4) if d else None

    result = {
        "n_reads": n_total,
        "n_reads_in_repeat": n_repeat,
        "n_mapped": n_mapped,
        "frac_mapped": frac(n_mapped, n_total),
        # Literature criterion (>=10% overlap of the union) -- comparable to
        # published minimap2/Winnowmap numbers.
        "n_correct": n_correct_lit,
        "n_correct_in_repeat": n_correct_lit_repeat,
        "frac_correct": frac(n_correct_lit, n_total),
        "frac_correct_in_repeat": frac(n_correct_lit_repeat, n_repeat),
        "frac_incorrectly_mapped": frac(n_mapped - n_correct_lit, n_total),
        # This repo's original stricter threshold (>=50% of the true span),
        # kept for continuity with earlier smoke-scale runs.
        "frac_correct_strict50": frac(n_correct_strict, n_total),
        "frac_correct_strict50_in_repeat": frac(n_correct_strict_repeat, n_repeat),
        "mean_identity": round(sum(identities) / len(identities), 4) if identities else None,
    }
    print(json.dumps(result, indent=2))
    if args.json_out:
        args.json_out.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
