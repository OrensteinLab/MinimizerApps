#!/usr/bin/env python3
"""Significance tests of the wall-time and CPU-time differences between minimizer orders, per app.

  order_tests.py --app {minimap2,sshash,kraken2,kmc} --results-dir DIR [--exclude-cold-gb G] [--out FILE.csv]

For every context (dataset, repeat-filter condition, ...) and timed metric (wall seconds and CPU
seconds of each step) it compares the candidate orders -- every GreedyMini and mod-minimizer
variant, plus xor and xor-pcg-hash -- with
  * the baseline order (native for minimap2/sshash, default for kmc; kraken2 has no such run, so
    its baseline is xor), and
  * xor, the random-order control (so a gain is not just "any non-default order");
using Welch's t-test on log times over the seeds/repeats of each order (ratio = ratio of geometric
means, 90% CI, p, Holm-adjusted over the candidates of a block, TOST equivalence at +-5%: the order
is "equivalent" iff the 90% CI of the ratio lies in [1/1.05, 1.05]). Groups need >= 2 runs; blocks
with fewer are reported with n and no test. No scipy (common/stats_util.py).

Timing runs are noisy on a shared host and a baseline with few repeats has little power: read the
n columns, and prefer the CPU-time rows to settle whether a wall-time difference is the order's
work (both differ) or load/cache (only wall time differs). --exclude-cold-gb G drops runs whose
step read more than G GB from disk (all four apps record it per step), a sensitivity check for the
page-cache confound; only the affected step of a run is dropped. Zero or missing times are skipped.
TOST equivalence is a separate test per comparison and is not Holm-adjusted.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats_util import holm, t_quantile, welch  # noqa: E402

TOST_MARGIN = 1.05
CANDIDATE_PREFIXES = ("greedymini", "mod-minimizer")
CANDIDATES_EXTRA = ("xor", "xor-pcg-hash")


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def label_of(r):
    """order + variant: greedymini/at0-cg1[/binary], mod-minimizer-xor/9, ..."""
    parts = [r.get("order") or ""]
    for c in ("gm_k", "projection", "gm_lut", "mod_t"):
        v = r.get(c)
        if v not in (None, "", "nan"):
            parts.append(("k%d" % int(float(v))) if c == "gm_k" else str(int(float(v))) if c == "mod_t" else str(v))
    return "/".join(parts)


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# -- per-app loaders: rows of {ctx, label, rep, metrics{name: value}, read_gb} ------------------

def load_minimap2(results):
    rows = []
    for mp in sorted((results / "runs").glob("**/meta.json")):
        m = json.loads(mp.read_text())
        rel = mp.relative_to(results / "runs").parts
        if rel[0] == "smoke":
            continue
        met, rd = {}, {"build": 0.0, "map": 0.0}
        for step in ("build", "map"):
            s = m.get(step) or {}
            met[step + "_wall_s"] = num(s.get("wall_time_sec"))
            if s.get("user_s") is not None and s.get("sys_s") is not None:
                met[step + "_cpu_s"] = s["user_s"] + s["sys_s"]
            rd[step] += (s.get("fs_inputs") or 0) * 512 / 1e9
        if met.get("map_cpu_s") is None:     # old meta without CPU seconds: parse the kept time logs
            sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "minimap2" / "lib"))
            from write_meta import parse_time_v
            for step in ("build", "map"):
                t = parse_time_v(mp.parent / (step + "_time.txt"))
                if t.get("user_s") is not None:
                    met[step + "_cpu_s"] = t["user_s"] + t["sys_s"]
                if t.get("fs_inputs") is not None and not (m.get(step) or {}).get("fs_inputs"):
                    rd[step] = t["fs_inputs"] * 512 / 1e9
        rows.append({"ctx": (rel[0], rel[1]), "label": label_of({**m, "gm_lut": None}), "rep": m.get("seed_idx", 0),
                     "metrics": met, "read_gb": rd})
    return rows


def load_csv_app(path, ctx_cols, metric_cols, keep=None, read_cols=None):
    """read_cols: {step prefix of the metrics (build, query, ...): read_gb column of that step}."""
    rows = []
    for r in read_csv(path):
        if keep and not keep(r):
            continue
        rows.append({"ctx": tuple(r[c] for c in ctx_cols), "label": label_of(r), "rep": r.get("seed_idx", "0"),
                     "metrics": {c: num(r.get(c)) for c in metric_cols},
                     "read_gb": {st: (num(r.get(c)) or 0.0) for st, c in (read_cols or {}).items()}})
    return rows


def load(app, results):
    csvs = results / "csvs"
    if app == "minimap2":
        return load_minimap2(results), "native"
    if app == "sshash":
        return (load_csv_app(csvs / "sshash_runs.csv", ["dataset"],
                             ["build_wall_s", "build_cpu_s", "query_wall_s", "query_cpu_s"], None,
                             {"build": "build_read_gb", "query": "query_read_gb"}), "native")
    if app == "kmc":
        keep = lambda r: (r["n_bins"] == "1024" and r["threads"] == "24" and r["strict_memory"] == "False"
                          and r["mem_limit_gb"] in ("12", "32"))
        return (load_csv_app(csvs / "kmc3_runs.csv", ["dataset"],
                             ["count_wall_s", "count_cpu_s", "count_stage1_s", "count_stage2_s"], keep, {"count": "read_gb"}), "default")
    if app == "kraken2":
        a = load_csv_app(csvs / "kraken2_build_runs.csv", ["library"], ["build_wall_s", "build_cpu_s"],
                         None, {"build": "build_read_gb"})
        b = load_csv_app(csvs / "kraken2_runs.csv", ["dataset", "minimum_hit_groups"],
                         ["classify_run_s", "classify_wall_s", "classify_cpu_s"],
                         None, {"classify": "classify_read_gb"})
        return a + b, "xor"
    raise SystemExit("unknown app " + app)


# -- tests -----------------------------------------------------------------------------------

def compare(a, b):
    """Welch on log times of a vs b (>= 2 each): ratio of geometric means a/b, 90% CI, p."""
    la, lb = [math.log(x) for x in a], [math.log(x) for x in b]
    d, se, df, p = welch(la, lb)
    t = t_quantile(0.95, df) if math.isfinite(df) else 1.645
    lo, hi = math.exp(d - t * se), math.exp(d + t * se)
    return math.exp(d), lo, hi, p


def verdict(ratio, lo, hi, p_holm):
    if lo >= 1 / TOST_MARGIN and hi <= TOST_MARGIN:
        return "equivalent (+-5%)"
    if p_holm < 0.05:
        return "faster" if ratio < 1 else "slower"
    return "not different"


def run(app, results, exclude_cold_gb):
    rows, baseline = load(app, results)
    groups = defaultdict(list)       # (ctx, metric, label) -> values
    dropped = 0
    for r in rows:
        for k, v in r["metrics"].items():
            if v is None or v <= 0:
                continue
            step = k.split("_")[0]    # build_wall_s -> build; classify_processing_s -> classify; count_* -> count
            if exclude_cold_gb is not None and (r["read_gb"].get(step) or 0) > exclude_cold_gb:
                dropped += 1          # only the step that read from disk is dropped, not the run's other steps
                continue
            groups[(r["ctx"], k, r["label"])].append(v)
    if dropped:
        print("excluded %d step measurements that read > %g GB from disk" % (dropped, exclude_cold_gb))
    out = []
    keys = sorted({(c, m) for c, m, _ in groups})
    for ctx, metric in keys:
        labels = sorted(l for c, m, l in groups if (c, m) == (ctx, metric))
        cands = [l for l in labels if l.startswith(CANDIDATE_PREFIXES) or l in CANDIDATES_EXTRA]
        for ref in dict.fromkeys([baseline, "xor"]):
            if ref not in labels:
                continue
            block = [l for l in cands if l != ref]
            res = []
            for l in block:
                a, b = groups[(ctx, metric, l)], groups[(ctx, metric, ref)]
                if len(a) >= 2 and len(b) >= 2:
                    res.append((l,) + compare(a, b))
                else:
                    res.append((l, None, None, None, None))
            tested = [x for x in res if x[1] is not None]
            adj = dict(zip([x[0] for x in tested], holm([x[4] for x in tested]))) if tested else {}
            for l, ratio, lo, hi, p in res:
                row = {"app": app, "context": "/".join(ctx), "metric": metric, "order": l, "vs": ref,
                       "n_order": len(groups[(ctx, metric, l)]), "n_vs": len(groups[(ctx, metric, ref)])}
                if ratio is None:
                    row.update(ratio="", ci90_lo="", ci90_hi="", p="", p_holm="", verdict="n<2 (no test)")
                else:
                    ph = adj[l]
                    row.update(ratio=round(ratio, 4), ci90_lo=round(lo, 4), ci90_hi=round(hi, 4), p=round(p, 6),
                               p_holm=round(ph, 6), verdict=verdict(ratio, lo, hi, ph))
                out.append(row)
    return out, baseline


def summarise(out, baseline):
    cnt = defaultdict(lambda: defaultdict(int))
    for r in out:
        fam = "greedymini" if r["order"].startswith("greedymini") else (
            "mod-minimizer" if r["order"].startswith("mod-minimizer") else r["order"])
        cnt[(r["metric"], fam, r["vs"])][r["verdict"]] += 1
    lines = []
    for (metric, fam, vs), c in sorted(cnt.items()):
        lines.append("%-24s %-16s vs %-8s  %s" % (metric, fam, vs,
                     "  ".join("%s:%d" % (k, v) for k, v in sorted(c.items()))))
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--app", required=True, choices=["minimap2", "sshash", "kraken2", "kmc"])
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--exclude-cold-gb", type=float, default=None,
                    help="drop runs that read more than this many GB from disk (sensitivity check)")
    ap.add_argument("--out", type=Path, default=None, help="CSV of every test (default: <results-dir>/csvs/order_tests_<app>.csv)")
    a = ap.parse_args(argv)
    out, baseline = run(a.app, a.results_dir, a.exclude_cold_gb)
    dest = a.out or a.results_dir / "csvs" / ("order_tests_%s.csv" % a.app)
    dest.parent.mkdir(parents=True, exist_ok=True)
    cols = ["app", "context", "metric", "order", "vs", "n_order", "n_vs", "ratio", "ci90_lo", "ci90_hi", "p", "p_holm", "verdict"]
    with open(dest, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    print("wrote %s (%d tests; baseline %s)" % (dest, len(out), baseline))
    for line in summarise(out, baseline):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
