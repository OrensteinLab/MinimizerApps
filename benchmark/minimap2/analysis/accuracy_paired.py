#!/usr/bin/env python3
"""Paired per-read accuracy of every order against native, per (dataset,
mid_occ_cap) block, from each run's map.paf and the dataset's truth.tsv.

Correctness is lib/score_accuracy.py's criterion (paftools mapeval: same
contig, same strand, overlap >= 10% of the union of true and mapped spans),
recomputed per read so
two orders can be compared on the *same* reads:

- McNemar exact test on the discordant reads (b = native right / order wrong,
  c = native wrong / order right), two-sided binomial p.
- 95% CI of the paired difference in correct rate, (c - b) / n, with the
  paired variance (b + c - (b - c)^2 / n) / n^2.
- MAPQ curve: for each threshold q, the fraction of reads with a primary at
  MAPQ >= q and the error rate among them (minimap2's own operating-point
  metric).

Writes analysis_accuracy_paired.csv and analysis_mapq_curve.csv under
--results-dir. Usage:
  python3.11 analysis/accuracy_paired.py --results-dir $MINIMAP2_RESULTS_DIR --data-dir $MINIMAP2_DATA_DIR
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from score_accuracy import load_truth, overlap_len  # noqa: E402

MAPQ_THRESHOLDS = (0, 1, 5, 10, 20, 30, 40, 50, 60)


def primary_hits(paf: Path) -> dict[str, tuple[str, int, int, int, str]]:
    """{read: (tname, tstart, tend, mapq, strand)} for the primary (tp:A:P) line."""
    out: dict[str, tuple[str, int, int, int, str]] = {}
    is_primary: dict[str, bool] = {}
    with paf.open() as f:
        for line in f:
            cols = line.split("\t", 12)
            if len(cols) < 12:
                continue
            q, prim = cols[0], "\ttp:A:P" in line
            if q in out and (is_primary[q] or not prim):
                continue
            out[q] = (cols[5], int(cols[7]), int(cols[8]), int(cols[11]), cols[4])
            is_primary[q] = prim
    return out


def correctness(truth: dict, hits: dict) -> tuple[list[str], list[int], list[int]]:
    """Per read (sorted ids): correct (1/0) and MAPQ (-1 = unmapped)."""
    ids = sorted(truth)
    ok, mq = [], []
    for rid in ids:
        t = truth[rid]
        h = hits.get(rid)
        if h is None:
            ok.append(0), mq.append(-1)
            continue
        tname, ts, te, q, strand = h
        same = tname == t["contig"] and strand == t["strand"]
        ov = overlap_len(ts, te, t["start"], t["end"]) if same else 0
        union = (t["end"] - t["start"]) + (te - ts) - ov
        ok.append(1 if same and union > 0 and ov / union >= 0.10 else 0)
        mq.append(q)
    return ids, ok, mq


def mcnemar_p(b: int, c: int) -> float:
    """Two-sided exact McNemar p (binomial, p=0.5, on b + c discordant)."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    # P(X <= k) for X ~ Bin(n, 0.5), doubled, via log-space sums
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2) for i in range(k + 1)]
    m = max(logs)
    tail = math.exp(m) * sum(math.exp(x - m) for x in logs)
    return min(1.0, 2 * tail)


def runs_by_block(results: Path) -> dict[tuple[str, str], dict[str, list[Path]]]:
    """{(dataset, cap): {order: [run dirs sorted by seed_idx]}}."""
    out: dict[tuple[str, str], dict[str, list[Path]]] = {}
    bases: dict[str, str] = {}
    for paf in sorted((results / "runs").glob("*/*/*/*/map.paf")):
        run = paf.parent
        order_dir, cond_dir, ds_dir = run.parent, run.parent.parent, run.parent.parent.parent
        if ds_dir.name == "smoke":
            continue
        cap = cond_dir.name.split("_f")[1] if "_f" in cond_dir.name else "default"
        base = cond_dir.name.split("_f")[0]  # k<k>_w<w>
        if bases.setdefault(ds_dir.name, base) != base:
            raise SystemExit(f"error: {ds_dir.name} has runs at both {bases[ds_dir.name]} and {base}; move one aside")
        out.setdefault((ds_dir.name, cap), {}).setdefault(order_dir.name, []).append(run)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--data-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    paired_rows, curve_rows = [], []
    truths: dict[str, dict] = {}
    for (ds, cap), orders in sorted(runs_by_block(args.results_dir).items()):
        if not (args.data_dir / ds / "truth.tsv").is_file():
            continue  # real reads (D3): no per-read correctness
        if ds not in truths:
            truths[ds] = load_truth(args.data_dir / ds / "truth.tsv")
        truth = truths[ds]
        rep = [truth[r]["in_repeat"] for r in sorted(truth)]
        if "native" not in orders:
            continue
        _, base_ok, _ = correctness(truth, primary_hits(orders["native"][0] / "map.paf"))
        for order, runs in sorted(orders.items()):
            for run in runs:
                _, ok, mq = correctness(truth, primary_hits(run / "map.paf"))
                for subset, mask in (("all", None), ("in_repeat", rep)):
                    idx = range(len(ok)) if mask is None else [i for i, m in enumerate(mask) if m]
                    n = len(idx)
                    b = sum(1 for i in idx if base_ok[i] and not ok[i])
                    c = sum(1 for i in idx if ok[i] and not base_ok[i])
                    diff = (c - b) / n if n else 0.0
                    var = (b + c - (b - c) ** 2 / n) / n ** 2 if n else 0.0
                    half = 1.96 * math.sqrt(max(var, 0.0))
                    paired_rows.append({
                        "dataset": ds, "mid_occ_cap": cap, "order": order, "seed_idx": run.name,
                        "subset": subset, "n_reads": n,
                        "correct": sum(ok[i] for i in idx), "native_correct": sum(base_ok[i] for i in idx),
                        "b_native_only": b, "c_order_only": c,
                        "diff_correct_rate": round(diff, 6),
                        "ci95_lo": round(diff - half, 6), "ci95_hi": round(diff + half, 6),
                        "mcnemar_p": round(mcnemar_p(b, c), 6),
                    })
                n = len(ok)
                for q in MAPQ_THRESHOLDS:
                    sel = [i for i in range(n) if mq[i] >= q]
                    wrong = sum(1 for i in sel if not ok[i])
                    curve_rows.append({
                        "dataset": ds, "mid_occ_cap": cap, "order": order, "seed_idx": run.name,
                        "mapq_ge": q, "frac_reads": round(len(sel) / n, 6),
                        "error_rate": round(wrong / len(sel), 6) if sel else None, "n_wrong": wrong,
                    })
            print(f"{ds} {cap} {order}: {len(runs)} run(s)", file=sys.stderr)

    for name, rows in (("analysis_accuracy_paired.csv", paired_rows), ("analysis_mapq_curve.csv", curve_rows)):
        path = args.results_dir / name
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {path} ({len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
