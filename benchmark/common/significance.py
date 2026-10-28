#!/usr/bin/env python3
"""Whether each scheme changes each metric of the paper's figures significantly, paired by seed.

    significance.py APP=RESULTS_DIR [APP=RESULTS_DIR ...]

Every run of a non-default scheme is paired with the run of the default scheme of the same setting
(dataset, ...) and seed -- the same round of the interleaved protocol (the mean of the default's runs
when it has several; minimap2: only its fixed-threshold block) -- and gives d = ln(metric / metric of
the default). Per (setting, metric, scheme), over the n paired seeds:
  * ratio = exp(mean d), the geometric mean of the ratios to the default, with its 95% CI;
  * p = two-sided paired t-test of mean d = 0 (df = n - 1; p = 0 when all d are equal and non-zero);
  * p_holm = Holm-adjusted over the schemes of the same setting and metric;
  * equivalent = TOST at +-5%: the 90% CI of the ratio lies in [1/1.05, 1.05];
  * verdict = "lower" / "higher" (p_holm < 0.05), else "equivalent" (TOST), else "inconclusive".
Size metrics barely vary between seeds, so their tests are decisive even for small changes; the test
matters for the times, which vary between runs on a shared server.

Each app's export adds the sheet "significance" to its workbook (sheets()) and writes
<RESULTS_DIR>/csvs/significance_<app>.csv; run by hand, the script writes the CSVs and prints the
verdicts. Standard library only.
"""

import argparse
import csv
import math
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from density_effect import ROWS, VARIANT, mean_of, num  # noqa: E402
from stats_util import holm, t_quantile, t_two_sided_p  # noqa: E402
from xlsx_util import write_csv  # noqa: E402

ALPHA = 0.05
TOST_MARGIN = 1.05

# app -> list of (runs csv, setting columns, default order, the metrics of the paper's figures)
SPEC = {
    "kraken2": [
        ("kraken2_build_runs.csv", ("library",), "xor", ["build_wall_s", "index_size_gb", "build_max_rss_gb"]),
        ("kraken2_runs.csv", ("dataset", "minimum_hit_groups"), "xor", ["classify_run_s", "classify_max_rss_gb"]),
    ],
    "minimap2": [
        ("minimap2_runs.csv", ("dataset", "mid_occ_threshold", "k"), "native",
         ["build_wall_s", "map_run_s", "index_size_gb", "map_max_rss_gb"]),
    ],
    "kmc": [
        ("kmc3_runs.csv", ("dataset", "n_bins", "threads", "mem_limit_gb"), "default",
         ["count_wall_s", "count_stage1_s", "count_stage2_s", "tmp_disk_gb"]),
    ],
    "sshash": [
        ("sshash_runs.csv", ("dataset",), "native",
         ["build_wall_s", "query_run_s", "bits_per_kmer", "build_max_rss_gb"]),
    ],
}

# rows not tested: settings and schemes the paper does not report
SKIP = {
    "kraken2": lambda r: r.get("minimum_hit_groups") not in (None, "", "2"),   # the paper: minimum hit groups 2
    "sshash": lambda r: r.get("gm_k") == "13",                                  # 13-mer GreedyMini tables: not used
}

COLS = ["setting", "metric", "scheme", "n", "ratio", "ci95_lo", "ci95_hi", "p", "p_holm", "equivalent", "verdict"]
DOC = ("paired by seed with the default scheme of the same setting: ratio = geometric mean of metric / default "
       "over the seeds, 95% CI, paired t-test on the log ratios (p), Holm-adjusted over the schemes of a setting "
       "and metric (p_holm), TOST at +-5% (equivalent); verdict lower / higher if p_holm < 0.05, else equivalent "
       "or inconclusive (benchmark/common/significance.py)")


def log_ratios(results, app):
    """{(setting, metric, scheme): [ln(metric / default), one per paired seed]} of one result set."""
    out = defaultdict(list)
    for fname, ctx, dflt, metrics in SPEC[app]:
        path = Path(results) / "csvs" / fname
        if not path.is_file():
            continue
        with open(path, newline="") as f:
            rows = ROWS.get(fname, list)(list(csv.DictReader(f)))
        rows = [r for r in rows if not SKIP.get(app, lambda _r: False)(r)]
        base = defaultdict(list)
        for r in rows:
            if r.get("order") == dflt:
                base[(tuple(r.get(c, "") for c in ctx), r.get("seed_idx"))].append(r)
        for r in rows:
            if r.get("order") == dflt:
                continue
            setting = tuple(r.get(c, "") for c in ctx)
            b = base.get((setting, r.get("seed_idx")))
            if not b:
                continue
            scheme = "/".join(r.get(c, "") for c in VARIANT if r.get(c, "") not in ("", "nan"))
            for m in metrics:
                v, vb = num(r.get(m)), mean_of(b, m)
                if v and vb and v > 0 and vb > 0:
                    out[("/".join(setting), m, scheme)].append(math.log(v / vb))
    return out


def paired_test(d):
    """(ratio, ci95_lo, ci95_hi, p, equivalent) of the log ratios d of n >= 2 paired seeds."""
    n, md = len(d), statistics.mean(d)
    sd = statistics.stdev(d)
    if sd < 1e-12:   # identical ratios in every seed (a size that does not depend on the seed)
        r = round(math.exp(md), 5)
        return r, r, r, (0.0 if abs(md) > 1e-12 else 1.0), (1 / TOST_MARGIN <= math.exp(md) <= TOST_MARGIN)
    se = sd / math.sqrt(n)
    p = t_two_sided_p(abs(md / se), n - 1)
    h95, h90 = t_quantile(0.975, n - 1) * se, t_quantile(0.95, n - 1) * se
    equivalent = math.exp(md - h90) >= 1 / TOST_MARGIN and math.exp(md + h90) <= TOST_MARGIN
    return round(math.exp(md), 5), round(math.exp(md - h95), 5), round(math.exp(md + h95), 5), p, equivalent


def rows_of(results, app):
    """The significance rows of one result set, Holm-adjusted per (setting, metric)."""
    groups = defaultdict(list)
    for (setting, metric, scheme), d in log_ratios(results, app).items():
        row = {"setting": setting, "metric": metric, "scheme": scheme, "n": len(d)}
        if len(d) >= 2:
            row["ratio"], row["ci95_lo"], row["ci95_hi"], row["p"], row["equivalent"] = paired_test(d)
        groups[(setting, metric)].append(row)
    order = [m for _, _, _, ms in SPEC[app] for m in ms]
    out = []
    for (setting, metric) in sorted(groups, key=lambda k: (order.index(k[1]), k[0])):
        rows = sorted(groups[(setting, metric)], key=lambda r: r["scheme"])
        tested = [r for r in rows if "p" in r]
        for r, adj in zip(tested, holm([r["p"] for r in tested])):
            r["p_holm"] = adj
            if adj < ALPHA:
                r["verdict"] = "lower" if r["ratio"] < 1 else "higher"
            else:
                r["verdict"] = "equivalent" if r["equivalent"] else "inconclusive"
        for r in tested:
            r["p"], r["p_holm"] = float("%.3g" % r["p"]), float("%.3g" % r["p_holm"])
        out += rows
    return out


def results_of(app, results):
    """The significance rows of one result set; also writes csvs/significance_<app>.csv."""
    rows = rows_of(results, app)
    write_csv(os.path.join(str(results), "csvs", "significance_%s.csv" % app), COLS, rows)
    return rows


def sheets(app, results):
    """The sheet an app's export adds to its workbook (empty list when there is nothing to test)."""
    rows = results_of(app, results)
    if not rows:
        return []
    return [{"name": "significance", "fields": COLS, "rows": rows, "banner": ("note", DOC),
             "block_fields": ["metric"], "merge_fields": ["metric"]}]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="+", metavar="APP=RESULTS_DIR")
    a = ap.parse_args(argv)
    for item in a.sets:
        app, _, results = item.partition("=")
        if app not in SPEC or not results:
            ap.error("expected APP=RESULTS_DIR with APP in %s, got %r" % (", ".join(SPEC), item))
        for r in results_of(app, results):
            print("%-8s %-22s %-36s %-34s n=%d ratio=%s p_holm=%s %s" % (
                app, r["metric"], r["setting"], r["scheme"], r["n"], r.get("ratio"), r.get("p_holm"),
                r.get("verdict", "")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
