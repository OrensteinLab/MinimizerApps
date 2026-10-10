#!/usr/bin/env python3
"""Statistics for the scheme analysis (no scipy needed): reads the outputs of
run_seed_experiments.sh, accuracy_paired.py and lib/export.py under
--results-dir (csvs/analysis/<ds>/{seeds,timing}/ and csvs/analysis_*.csv) and
writes analysis_*.csv tables into <results-dir>/csvs.

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
- accuracy: per (dataset, cap, order), from accuracy_paired.py's
  analysis_accuracy_order.csv: the paired correct-rate difference vs native
  with reads as the independent unit -- per read, d_i = (mean over the
  order's seeds of correct_i) - native correct_i (all seeds share one native
  run and one read set, so seeds are not independent replicates and are never
  pooled as n * n_seeds pairs); mean(d), SE = sd(d) / sqrt(n), its 90% CI
  (t, df = n - 1) and TOST at +-0.1 percentage point (equivalent iff that CI
  lies inside [-0.001, 0.001]); p from the paired t-test on d over reads, or
  the exact McNemar test for a single-run (seed-independent) order; Holm over
  the orders of each (dataset, cap, subset) block. Datasets without truth
  (D3, real reads) have no accuracy rows.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from order_names import ID_FIELDS, parse_order_dir  # noqa: E402
from units import gb_from_bytes  # noqa: E402

# A "variant" is an order's run-directory name (order_names.order_dir_name):
# `native`, `xor`, `greedymini/k15_at1-cg0`, `mod-minimizer-xor/t9`, ... Its
# files under analysis/ use the name with "/" as "__"; every output row gets
# the separate identity columns (order, gm_k, projection, gm_table, mod_t).
CANONICAL = ["native", "xor", "xor-pcg-hash", "greedymini", "mod-minimizer-xor"]


def variant_key(variant: str) -> tuple:
    order = variant.split("/")[0]
    return (CANONICAL.index(order) if order in CANONICAL else len(CANONICAL), variant)


def id_cols(variant: str) -> dict:
    d = parse_order_dir(variant)
    return {f: d[f] for f in ID_FIELDS}


def variants_in(directory: Path, pattern: str) -> list[str]:
    """Variants having a file `<label>` matching `pattern` (a regex with one
    group for the label) in `directory`, in canonical order."""
    rx = re.compile(pattern)
    found = {m.group(1).replace("__", "/") for f in directory.iterdir() if (m := rx.fullmatch(f.name))}
    return sorted(found, key=variant_key)


def label(variant: str) -> str:
    return variant.replace("/", "__")
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
        for order in variants_in(cap_dir, r"(.+)\.csv"):
            f = cap_dir / f"{label(order)}.csv"
            rr = read_csv(f)
            kb = sum(int(r["qlen"]) for r in rr) / 1000
            seeds = sum(int(r["n_seeds"]) for r in rr)
            anchors = [int(r["anchors"]) for r in rr]
            rows.append({
                "dataset": ds, **id_cols(order), "mid_occ_cap": cap_dir.name, "n_reads": len(rr),
                "seeds_per_kb": round(seeds / kb, 2),
                "frac_seeds_filtered": round(sum(int(r["n_filtered"]) for r in rr) / seeds, 5),
                "anchors_per_read_mean": round(statistics.mean(anchors), 1),
                "anchors_per_read_median": statistics.median(anchors),
                "anchors_per_kb": round(sum(anchors) / kb, 1),
                "occ_before_filter_per_kb": round(sum(int(r["occ_all"]) for r in rr) / kb, 1),
                "km_cap_gb_mean": gb_from_bytes(statistics.mean(int(r["km_cap_bytes"]) for r in rr), 6),
                "km_cap_gb_max": gb_from_bytes(max(int(r["km_cap_bytes"]) for r in rr), 6),
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
        for order in variants_in(cap_dir, r"(.+)\.r\d+\.csv"):
            reps = sorted(cap_dir.glob(f"{label(order)}.r*.csv"))
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
        for order in sorted(totals, key=variant_key):
            if order not in totals:
                continue
            lx = [math.log(v) for v in totals[order]]
            row = {
                "dataset": ds, **id_cols(order), "mid_occ_cap": cap, "n_reps": len(lx),
                "sample_wall_s_mean": round(statistics.mean(totals[order]), 2),
                "sample_wall_s_std": round(statistics.stdev(totals[order]), 2) if len(lx) > 1 else None,
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
                if lr:  # no reads timed in both runs: leave the per-read columns blank
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
            sf = adir / ds / "seeds" / cap / f"{label(order)}.csv"
            if not sf.is_file():
                continue
            anchors = {r["read"]: int(r["anchors"]) for r in read_csv(sf)}
            xs, ys = [], []
            for k, t in per_read[order].items():
                if k in anchors and anchors[k] > 0 and t > 0:
                    xs.append(math.log(anchors[k])); ys.append(math.log(t))
            per_read_rows.append({"dataset": ds, **id_cols(order), "mid_occ_cap": cap,
                                  "n_reads": len(xs), "pearson_log_anchors_log_time": r4(pearson(xs, ys))})
    return rows, per_read_rows


def accuracy(results: Path) -> list[dict]:
    f = results / "csvs" / "analysis_accuracy_order.csv"
    if not f.is_file():
        if (results / "csvs" / "analysis_accuracy_paired.csv").is_file():
            print(f"warning: no {f.name}; rerun analysis/accuracy_paired.py (the per-seed "
                  "analysis_accuracy_paired.csv alone can't give a reads-as-unit test)")
        return []
    out = []
    for r in read_csv(f):
        if r["order"] == "native" or not r.get("diff_correct_rate") or not r.get("p"):
            continue  # native itself, or an empty subset (no in-repeat reads)
        d, se, df = float(r["diff_correct_rate"]), float(r["diff_se"]), int(r["df"])
        tq = t_quantile(0.95, df) if df > 0 else 1.645
        lo, hi = d - tq * se, d + tq * se
        out.append({
            "dataset": r["dataset"], **{f: r[f] for f in ID_FIELDS}, "mid_occ_cap": r["mid_occ_cap"], "subset": r["subset"],
            "n_seeds": int(r["n_seeds"]), "n_reads": int(r["n_reads"]),
            "diff_correct_rate_mean": round(d, 6), "diff_se": round(se, 8),
            "diff_ci90_lo": round(lo, 6), "diff_ci90_hi": round(hi, 6),
            "tost_equivalent_0.1pp": lo >= -ACC_MARGIN and hi <= ACC_MARGIN,
            "test": r["test"], "p": float(r["p"]),
        })
    for key in {(r["dataset"], r["mid_occ_cap"], r["subset"]) for r in out}:
        block = [r for r in out if (r["dataset"], r["mid_occ_cap"], r["subset"]) == key]
        for r, adj in zip(block, holm([r["p"] for r in block])):
            r["p_holm"] = round(adj, 6)
            r["p"] = round(r["p"], 6)
    return sorted(out, key=lambda r: (r["dataset"], r["mid_occ_cap"], r["subset"],
                                      CANONICAL.index(r["order"]) if r["order"] in CANONICAL else len(CANONICAL),
                                      r["gm_k"], r["projection"], r["mod_t"]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    # D1, D2, D3 (lib/ and README.md: Datasets); D3 (real reads, no truth) gets
    # mechanism/timing rows only -- accuracy_paired.py skips it.
    ap.add_argument("--datasets", default="chm13_chrX,grch38_wg,chm13_rel3_ont")
    args = ap.parse_args(argv)
    adir = args.results_dir / "csvs" / "analysis"
    (args.results_dir / "csvs").mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)
    mech, tim, prr = [], [], []
    for ds in args.datasets.split(","):
        if not (adir / ds).is_dir():
            continue
        mech += mechanism(adir, ds)
        t, p = timing(adir, ds, rng)
        tim += t
        prr += p
    write_csv(args.results_dir / "csvs" / "analysis_mechanism.csv", mech)
    write_csv(args.results_dir / "csvs" / "analysis_timing.csv", tim)
    write_csv(args.results_dir / "csvs" / "analysis_per_read_anchor_time.csv", prr)
    write_csv(args.results_dir / "csvs" / "analysis_accuracy_tests.csv", accuracy(args.results_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
