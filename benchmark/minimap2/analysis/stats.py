#!/usr/bin/env python3
"""Statistics for the scheme analysis (no scipy needed): reads the outputs of
run_seed_experiments.sh, accuracy_paired.py and lib/export.py under
--results-dir and writes analysis_*.csv tables next to them.

- mechanism: per (dataset, cap, order), from the --dbg-seed-occ pass: read
  minimizers in the index per kb, fraction dropped by mid_occ, anchors per
  read (chaining input), reference occurrences before filtering, thread-local
  memory per read.
- timing: per (dataset, cap, order), total sample time per replicate; ratio to
  native; Welch t-test on log totals (Holm over the orders of a block); 90% CI
  of the ratio and TOST equivalence at +-5% (equivalent iff the 90% CI lies in
  [0.95, 1/0.95]); paired per-read median log2 ratio with a bootstrap CI.
- per-read mechanism: within each order, Pearson r between log anchors and log
  per-read time across reads.
- accuracy: per (dataset, cap, order), paired correct-rate difference vs
  native (all seeds), Holm-adjusted McNemar p, and TOST at +-0.1 percentage
  point (90% CI of the difference inside [-0.001, 0.001]).
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
from pathlib import Path

ORDERS = ["native", "xor", "xor-pcg-hash", "greedymini", "greedymini-at1", "mod-minimizer-xor"]
TIME_MARGIN = 0.05
ACC_MARGIN = 0.001


# -- t distribution (regularized incomplete beta, Numerical Recipes betacf) --

def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c if abs(c) > 1e-300 else 1e300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c if abs(c) > 1e-300 else 1e300
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def t_two_sided_p(t: float, df: float) -> float:
    return _betai(df / 2, 0.5, df / (df + t * t))


def t_quantile(p: float, df: float) -> float:
    """Upper quantile: t such that P(T > t) = 1 - p, by bisection."""
    lo, hi = 0.0, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 1 - t_two_sided_p(mid, df) / 2 < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def welch(x: list[float], y: list[float]) -> tuple[float, float, float, float]:
    """(mean diff x-y, se, df, two-sided p)."""
    mx, my = statistics.mean(x), statistics.mean(y)
    vx, vy = statistics.variance(x), statistics.variance(y)
    se2 = vx / len(x) + vy / len(y)
    if se2 == 0:
        return mx - my, 0.0, float("inf"), 0.0 if mx != my else 1.0
    df = se2 ** 2 / ((vx / len(x)) ** 2 / (len(x) - 1) + (vy / len(y)) ** 2 / (len(y) - 1))
    t = (mx - my) / math.sqrt(se2)
    return mx - my, math.sqrt(se2), df, t_two_sided_p(abs(t), df)


def holm(ps: list[float]) -> list[float]:
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj = [0.0] * len(ps)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(ps) - rank) * ps[i]))
        adj[i] = running
    return adj


# -- helpers --

def read_csv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else None


def r4(v):
    return None if v is None else round(v, 4)


# -- sections --

def mechanism(adir: Path, ds: str) -> list[dict]:
    rows = []
    for cap_dir in sorted((adir / ds / "seeds").glob("*")):
        for order in ORDERS:
            f = cap_dir / f"{order}.csv"
            if not f.is_file():
                continue
            rr = read_csv(f)
            kb = sum(int(r["qlen"]) for r in rr) / 1000
            seeds = sum(int(r["n_seeds"]) for r in rr)
            anchors = [int(r["anchors"]) for r in rr]
            rows.append({
                "dataset": ds, "mid_occ_cap": cap_dir.name, "order": order, "n_reads": len(rr),
                "seeds_per_kb": round(seeds / kb, 2),
                "frac_seeds_filtered": round(sum(int(r["n_filtered"]) for r in rr) / seeds, 5),
                "anchors_per_read_mean": round(statistics.mean(anchors), 1),
                "anchors_per_read_median": statistics.median(anchors),
                "anchors_per_kb": round(sum(anchors) / kb, 1),
                "occ_before_filter_per_kb": round(sum(int(r["occ_all"]) for r in rr) / kb, 1),
                "km_cap_mb_mean": round(statistics.mean(int(r["km_cap_bytes"]) for r in rr) / 2**20, 2),
                "km_cap_mb_max": round(max(int(r["km_cap_bytes"]) for r in rr) / 2**20, 2),
            })
    base = {(r["mid_occ_cap"]): r for r in rows if r["order"] == "native"}
    for r in rows:
        b = base.get(r["mid_occ_cap"])
        r["anchors_ratio"] = round(r["anchors_per_kb"] / b["anchors_per_kb"], 4) if b else None
    return rows


def timing(adir: Path, ds: str, rng: random.Random) -> tuple[list[dict], list[dict]]:
    rows, per_read_rows = [], []
    for cap_dir in sorted((adir / ds / "timing").glob("*")):
        cap = cap_dir.name
        totals: dict[str, list[float]] = {}
        per_read: dict[str, dict[str, float]] = {}
        for order in ORDERS:
            reps = sorted(cap_dir.glob(f"{order}.r*.csv"))
            if not reps:
                continue
            acc: dict[str, list[float]] = {}
            for f in reps:
                rr = read_csv(f)
                totals.setdefault(order, []).append(sum(float(r["time_sec"]) for r in rr))
                for r in rr:
                    acc.setdefault(r["read"], []).append(float(r["time_sec"]))
            per_read[order] = {k: statistics.mean(v) for k, v in acc.items()}
        if "native" not in totals:
            continue
        nat = [math.log(v) for v in totals["native"]]
        block = []
        for order in ORDERS:
            if order not in totals:
                continue
            lx = [math.log(v) for v in totals[order]]
            row = {
                "dataset": ds, "mid_occ_cap": cap, "order": order, "n_reps": len(lx),
                "sample_time_sec_mean": round(statistics.mean(totals[order]), 2),
                "sample_time_sec_sd": round(statistics.stdev(totals[order]), 2) if len(lx) > 1 else None,
            }
            if order != "native" and len(lx) > 1 and len(nat) > 1:
                d, se, df, p = welch(lx, nat)
                tq = t_quantile(0.95, df) if math.isfinite(df) else 1.645
                lo, hi = math.exp(d - tq * se), math.exp(d + tq * se)
                row.update({
                    "ratio": round(math.exp(d), 4), "ci90_lo": round(lo, 4), "ci90_hi": round(hi, 4),
                    "welch_p": p,
                    "tost_equivalent_5pct": lo >= 1 - TIME_MARGIN and hi <= 1 / (1 - TIME_MARGIN),
                })
                # paired per read: log2(order / native), averaged over replicates per read
                common = [k for k in per_read[order] if k in per_read["native"]]
                lr = [math.log2(per_read[order][k] / per_read["native"][k]) for k in common
                      if per_read[order][k] > 0 and per_read["native"][k] > 0]
                med = statistics.median(lr)
                boots = sorted(statistics.median(rng.choices(lr, k=len(lr))) for _ in range(500))
                row.update({
                    "per_read_median_log2_ratio": round(med, 4),
                    "per_read_ci95_lo": round(boots[12], 4), "per_read_ci95_hi": round(boots[487], 4),
                })
            block.append(row)
        tested = [r for r in block if "welch_p" in r]
        for r, adj in zip(tested, holm([r["welch_p"] for r in tested])):
            r["welch_p_holm"] = round(adj, 6)
            r["welch_p"] = round(r["welch_p"], 6)
        rows.extend(block)
        # per-read mechanism: log anchors vs log time, within order
        for order in per_read:
            sf = adir / ds / "seeds" / cap / f"{order}.csv"
            if not sf.is_file():
                continue
            anchors = {r["read"]: int(r["anchors"]) for r in read_csv(sf)}
            xs, ys = [], []
            for k, t in per_read[order].items():
                if k in anchors and anchors[k] > 0 and t > 0:
                    xs.append(math.log(anchors[k])); ys.append(math.log(t))
            per_read_rows.append({"dataset": ds, "mid_occ_cap": cap, "order": order,
                                  "n_reads": len(xs), "pearson_log_anchors_log_time": r4(pearson(xs, ys))})
    return rows, per_read_rows


def accuracy(results: Path) -> list[dict]:
    f = results / "analysis_accuracy_paired.csv"
    if not f.is_file():
        return []
    rows = [r for r in read_csv(f) if r["order"] != "native"]
    out = []
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        groups.setdefault((r["dataset"], r["mid_occ_cap"], r["order"], r["subset"]), []).append(r)
    for (ds, cap, order, subset), rr in sorted(groups.items()):
        n = int(rr[0]["n_reads"])
        diffs = [float(r["diff_correct_rate"]) for r in rr]
        # pooled over seeds: sum discordant counts (each seed is its own paired comparison)
        b = sum(int(r["b_native_only"]) for r in rr)
        c = sum(int(r["c_order_only"]) for r in rr)
        N = n * len(rr)
        d = (c - b) / N
        se = math.sqrt(max(b + c - (b - c) ** 2 / N, 0)) / N
        out.append({
            "dataset": ds, "mid_occ_cap": cap, "order": order, "subset": subset, "n_seeds": len(rr),
            "n_reads": n, "diff_correct_rate_mean": round(statistics.mean(diffs), 6),
            "diff_ci90_lo": round(d - 1.645 * se, 6), "diff_ci90_hi": round(d + 1.645 * se, 6),
            "tost_equivalent_0.1pp": d - 1.645 * se >= -ACC_MARGIN and d + 1.645 * se <= ACC_MARGIN,
            "mcnemar_p_min": min(float(r["mcnemar_p"]) for r in rr),
        })
    for key in {(r["dataset"], r["mid_occ_cap"], r["subset"]) for r in out}:
        block = [r for r in out if (r["dataset"], r["mid_occ_cap"], r["subset"]) == key]
        for r, adj in zip(block, holm([r["mcnemar_p_min"] for r in block])):
            r["mcnemar_p_min_holm"] = round(adj, 6)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--datasets", default="chm13_chrX,grch38_wg,chm13_chrX_ont")
    args = ap.parse_args(argv)
    adir = args.results_dir / "analysis"
    rng = random.Random(11)
    mech, tim, prr = [], [], []
    for ds in args.datasets.split(","):
        if not (adir / ds).is_dir():
            continue
        mech += mechanism(adir, ds)
        t, p = timing(adir, ds, rng)
        tim += t
        prr += p
    write_csv(args.results_dir / "analysis_mechanism.csv", mech)
    write_csv(args.results_dir / "analysis_timing.csv", tim)
    write_csv(args.results_dir / "analysis_per_read_anchor_time.csv", prr)
    write_csv(args.results_dir / "analysis_accuracy_tests.csv", accuracy(args.results_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
