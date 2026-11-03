#!/usr/bin/env python3
"""How closely each runtime and space metric follows the density, across schemes (paired by seed).

    density_effect.py APP=RESULTS_DIR [APP=RESULTS_DIR ...] [--out WORKBOOK.xlsx]
    e.g. density_effect.py kraken2=/localdata/shay_sirek/kraken2/results minimap2=... --out density_effect.xlsx

Every run of a non-default scheme is paired with the run of the default scheme of the same setting
(dataset, ...) and seed (the same round of the interleaved protocol; the mean of the default's runs when
it has several; minimap2: only its fixed-threshold block, see ROWS), and gives one point
    x = density factor / density factor of the default - 1,   y = metric / metric of the default - 1,
the relative changes of the density and of the metric. Per (app, metric) it reports n, the number of
schemes, and the Pearson r of y and x (r near 1: the metric follows the density; near 0: it does not),
over all the settings (the ratios remove the scale of each dataset) and over the runs of each setting.
The density differs between schemes by up to ~30% and hardly between seeds of a scheme, so r compares
schemes; they also differ in other properties (the cost of computing the order, the repetitiveness and
locality of the selected k-mers), so r describes the density together with those differences.

Each app's export adds the sheets "density_effect" (the correlations) and "density_effect_points" to its workbook
(sheets()) and writes <RESULTS_DIR>/csvs/density_effect_<app>.csv; run by hand, the script writes the CSVs and,
with --out, one workbook of all the given apps. Standard library only.
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scheme_tests import minimap2_fixed_block  # noqa: E402
from xlsx_book import write_book  # noqa: E402
from xlsx_util import write_csv  # noqa: E402

VARIANT = ("scheme", "gm_k", "projection", "gm_table", "gm_lut", "mod_t")

# app -> list of (runs csv, setting columns, default scheme, metrics). A setting is a dataset with its parameters;
# minimap2's includes the occurrence threshold, and only its fixed-threshold block (threshold fixed to the
# default scheme's value) is paired: a derived threshold differs between schemes, and one that happens to
# equal the default's must not join the fixed block (ROWS picks the block by position, not by value).
SPEC = {
    "kraken2": [
        ("kraken2_build_runs.csv", ("library",), "xor", ["build_wall_s", "build_max_rss_gb", "index_size_gb"]),
        ("kraken2_runs.csv", ("dataset", "minimum_hit_groups"), "xor",
         ["classify_run_s", "classify_wall_s", "classify_max_rss_gb"]),
    ],
    "minimap2": [
        ("minimap2_runs.csv", ("dataset", "mid_occ_threshold", "k"), "native",
         ["build_wall_s", "build_max_rss_gb", "index_size_gb", "map_run_s", "map_wall_s", "map_max_rss_gb"]),
    ],
    "kmc": [
        ("kmc3_runs.csv", ("dataset", "n_bins", "threads", "mem_limit_gb"), "default",
         ["count_wall_s", "count_stage1_s", "count_stage2_s", "tmp_disk_gb", "count_max_rss_gb"]),
    ],
    "sshash": [
        ("sshash_runs.csv", ("dataset",), "native",
         ["build_wall_s", "build_max_rss_gb", "bits_per_kmer", "query_run_s", "query_wall_s",
          "bench_positive_ns", "bench_negative_ns"]),
    ],
}


# runs csv -> the rows of it that are paired (all rows when absent)
ROWS = {"minimap2_runs.csv": lambda rows: [r for r, f in zip(rows, minimap2_fixed_block(rows)[0]) if f]}

# The runs the paper does not report, left out of the correlations here and of the tests of significance.py: the
# GreedyMini-based schemes with the polarity at1-cg0 (the paper reports at0-cg1, chosen arbitrarily; the two are compared
# in its supplementary table polarity_comparison), KMC3's GreedyMini with the DNA-alphabet table, SSHash's GreedyMini
# orders on 13-mers, and Kraken2's classification at minimum hit groups 0.
UNREPORTED_PROJECTION = "at1-cg0"
SKIP = {
    "kraken2": lambda r: r.get("minimum_hit_groups") not in (None, "", "2"),
    "kmc": lambda r: r.get("gm_lut") == "dna",
    "sshash": lambda r: r.get("gm_k") == "13",
}


def unreported(r, app):
    """Whether the paper does not report the run r of app (SKIP)."""
    return r.get("projection") == UNREPORTED_PROJECTION or SKIP.get(app, lambda _r: False)(r)


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def mean_of(rows, col):
    """Mean of a column over the default's runs of one setting and seed, or None."""
    v = [x for x in (num(r.get(col)) for r in rows) if x is not None]
    return sum(v) / len(v) if v else None


def pearson(a, b):
    n = len(a)
    if n < 3:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    saa = sum((x - ma) ** 2 for x in a)
    sbb = sum((y - mb) ** 2 for y in b)
    if saa <= 1e-18 or sbb <= 1e-18:
        return None
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / math.sqrt(saa * sbb)


def points(results, app):
    """[(metric, point dict)] of one app's result set."""
    out = defaultdict(list)
    for fname, ctx, dflt, metrics in SPEC[app]:
        path = Path(results) / "csvs" / fname
        if not path.is_file():
            continue
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f))
        rows = [r for r in ROWS.get(fname, list)(rows) if not unreported(r, app)]
        base = defaultdict(list)   # the default's runs of a setting and seed
        for r in rows:
            if r.get("scheme") == dflt:
                base[(tuple(r.get(c, "") for c in ctx), r.get("seed_idx"))].append(r)
        for r in rows:
            if r.get("scheme") == dflt:
                continue
            b = base.get((tuple(r.get(c, "") for c in ctx), r.get("seed_idx")))
            if not b:
                continue
            d, db = num(r.get("density_factor")), mean_of(b, "density_factor")
            if not (d and db):
                continue
            for m in metrics:
                v, vb = num(r.get(m)), mean_of(b, m)
                if v is None or not vb:
                    continue
                out[m].append({"setting": "/".join(r.get(c, "") for c in ctx),
                               "scheme": "/".join(r.get(c, "") for c in VARIANT if r.get(c, "") not in ("", "nan")),
                               "seed_idx": r.get("seed_idx"), "metric": m,
                               "metric_ratio": round(v / vb, 5), "density_ratio": round(d / db, 5)})
    order = [m for _, _, _, ms in SPEC[app] for m in ms]
    return [(m, out[m]) for m in order if out[m]]


def correlate(app, metric, setting, pts):
    x = [p["density_ratio"] - 1 for p in pts]
    y = [p["metric_ratio"] - 1 for p in pts]
    r = pearson(x, y)
    return {"app": app, "metric": metric, "setting": setting, "n": len(pts),
            "schemes": len({p["scheme"] for p in pts}), "r": None if r is None else round(r, 3)}


def rows_of(app, metric, pts):
    """The correlation over all the runs of a metric (setting "all"), then over the runs of each setting."""
    by = defaultdict(list)
    for p in pts:
        by[p["setting"]].append(p)
    return [correlate(app, metric, "all", pts)] + [correlate(app, metric, s, by[s]) for s in sorted(by)]


DOC = "r = Pearson correlation of (metric/default - 1) with (density/default - 1), over the runs paired with the default scheme's run of the same setting and seed, over all settings (all) and per setting (benchmark/common/density_effect.py)"

COLS = ["app", "metric", "setting", "n", "schemes", "r"]
PCOLS = ["setting", "scheme", "seed_idx", "metric", "metric_ratio", "density_ratio"]


def results_of(app, results):
    """(correlation rows, points) of one result set; also writes csvs/density_effect_<app>.csv."""
    rows, pts_all = [], []
    for metric, pts in points(results, app):
        rows += rows_of(app, metric, pts)
        pts_all += pts
    write_csv(os.path.join(str(results), "csvs", "density_effect_%s.csv" % app), COLS, rows)
    return rows, pts_all


def sheets(app, results):
    """The two sheets an app's export adds to its workbook (empty list when there is nothing to correlate)."""
    rows, pts = results_of(app, results)
    if not rows:
        return []
    return [{"name": "density_effect", "fields": COLS, "rows": rows, "banner": ("note", DOC),
             "block_fields": ["metric"], "merge_fields": ["metric"]},
            {"name": "density_effect_points", "fields": PCOLS, "rows": pts,
             "block_fields": ["metric"], "merge_fields": ["metric"]}]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sets", nargs="+", metavar="APP=RESULTS_DIR")
    ap.add_argument("--out", default=None, help="workbook with the summary and the points of every app")
    a = ap.parse_args(argv)
    summary, books = [], []
    for item in a.sets:
        app, _, results = item.partition("=")
        if app not in SPEC or not results:
            ap.error("expected APP=RESULTS_DIR with APP in %s, got %r" % (", ".join(SPEC), item))
        rows, pts_all = results_of(app, results)
        summary += rows
        books.append({"name": "%s points" % app, "fields": PCOLS, "rows": pts_all,
                      "block_fields": ["metric"], "merge_fields": ["metric"]})
        for r in rows:
            print("%-8s %-20s %-28s n=%4d  r=%6s" % (app, r["metric"], r["setting"], r["n"], r["r"]))
    if a.out:
        write_book(a.out, [{"name": "density effect", "fields": COLS, "rows": summary,
                            "block_fields": ["app"], "merge_fields": ["app"]}] + books)
        print("wrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
