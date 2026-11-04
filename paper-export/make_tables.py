#!/usr/bin/env python3
"""The supplementary tables of the paper (Tables S1, S2, ...) as CSV files, computed from the result sets of the benchmarks.

    make_tables.py --results kraken2=DIR minimap2=DIR kmc=DIR sshash=DIR [--out DIR] [--labels FILE]

Each DIR is a result set of benchmark/<method> (the paper's: the full runs together with the profiling runs of
benchmark/<method>/experiments/profile_run.sh). The tables are computed from the CSVs that the export of each method
writes to DIR/csvs (runs, summaries, schemes, profiles) and, for Minimap2's mapping profile, from the map_profile.log
of each run. The statistics are computed here, with benchmark/common:
  * significance tests (significance.py): every run of a scheme is paired with the run of the default scheme of the
    same setting and seed; per scheme, dataset and metric, the geometric mean of the ratios to the default with its
    95% CI, a two-sided paired t-test on the log ratios, Bonferroni-corrected over the schemes compared with the default
    on the same dataset and metric, and TOST at +-5% for equivalence;
  * the comparison of the two polarities of Phi: the same test, pairing each GreedyMini-based scheme's run with A,T
    mapped to 1 with its run with A,T mapped to 0 of the same seed;
  * correlations with density (density_effect.py): the Pearson correlation of the paired ratios of the density and
    of a metric to those of the default scheme.

Writes --out/S<nn>_<name>.csv (default: <repo>/supplementary_tables), numbered in the order of their first citation in
the paper (TABLES), the list of the tables with their captions between the markers of --out/README.md, and, with
--labels, their LaTeX labels for the paper (\\newlabel{tab:<name>}{{\\href{<CSV on GitHub>}{S<n>}}{}}, which needs
hyperref, loaded with implicit=false so that the labels keep the two fields of plain LaTeX). The rows are those of the paper: minimum
hit groups 2 for Kraken2, the occurrence threshold of the default scheme for Minimap2 (and also the derived one in
minimap2_build_mapping), and the GreedyMini-based schemes with the polarity A,T -> 0 (common.POLARITY).
Standard library only.
"""

import argparse
import collections
import csv
import json
import math
import re
import statistics
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (APPS, DATASETS, GM_SCHEMES, KRAKEN2_LIBRARY, LABELS, METHOD, POLARITIES, POLARITY, REPO,  # noqa: E402
                    ResultSets, dataset_name, default_thresholds, fl, key_label, kw, label, label_key, minimap2_blocks,
                    scheme_order)
import density_effect  # noqa: E402
import significance  # noqa: E402
from stats_util import bonferroni  # noqa: E402

# ---------------------------------------------------------------- cells and tables

class Stat:
    """A mean and its standard deviation, as text; written as the two columns '<name> mean' and '<name> std', or as one
    column '<name>' when no row of the table has a standard deviation."""

    def __init__(self, mean, std=""):
        self.mean, self.std = mean, std


def stat(mean, sd, fmt, scale=1.0):
    """mean +- sd, both in the format of the mean; no sd when it is 0 (one run, or a value that does not vary); an sd
    that rounds to 0 in that format is given to one significant digit (e.g. 0.0003)."""
    mean, sd = mean * scale, (sd or 0.0) * scale
    if abs(sd) < 1e-9:
        return Stat(fmt.format(mean))
    t = fmt.format(sd)
    if float(t) == 0:
        t = "%.*f" % (-math.floor(math.log10(abs(sd))), sd)
    return Stat(fmt.format(mean), t)


def stat_of(row, col, fmt, scale=1.0):
    """The Stat of column col of a summary row (col_mean and col_std, or col)."""
    v = fl(row.get(col + "_mean", row.get(col)))
    if v is None:
        raise ValueError("missing %s in %s %s" % (col, row.get("dataset"), row.get("scheme")))
    return stat(v, fl(row.get(col + "_std")), fmt, scale)


def sd1(sd):
    """A standard deviation (CSV text) to one significant digit, rounded half up ('' when 0 or missing)."""
    if fl(sd) in (None, 0.0):
        return ""
    v = Decimal(sd)
    d = -v.adjusted()
    r = v.quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    if r.adjusted() > v.adjusted():   # 0.0095 -> 0.01
        d -= 1
        r = v.quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    return "{:f}".format(r)


def mean_text(m, decimals):
    """A mean (CSV text) to at most `decimals` decimals, without digits beyond the CSV's own precision."""
    t = m.split(".")[1] if "." in m else ""
    return "%.*f" % (min(decimals, len(t)) if t else decimals, fl(m))


def signed(r):
    """A correlation to two decimals, rounded half up, with its sign ('' when undefined)."""
    if r in (None, ""):
        return ""
    return "%+.2f" % Decimal(str(r)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)


def p_text(p):
    return "<0.001" if p < 0.001 else "%.3f" % p


class Table:
    """A supplementary table: its caption (for the README), columns and rows (cells: text, Stat, or None where a value
    does not apply)."""

    def __init__(self, name, caption, columns, rows):
        self.name, self.caption, self.columns, self.rows = name, caption, columns, rows
        for r in rows:
            assert len(r) == len(columns), (name, r)

    def csv_rows(self):
        split = [any(isinstance(r[i], Stat) and r[i].std for r in self.rows) for i in range(len(self.columns))]
        head = []
        for c, s in zip(self.columns, split):
            head += [c + " mean", c + " std"] if s else [c]
        out = [head]
        for r in self.rows:
            cells = []
            for v, s in zip(r, split):
                if isinstance(v, Stat):
                    cells += [v.mean, v.std] if s else [v.mean]
                else:
                    cells += ["" if v is None else str(v)] + ([""] if s else [])
            out.append(cells)
        return out


def scheme_rows(rows, app, keep=lambda r: True, last=False):
    """[(label, row)] of the reported schemes among rows, in the paper's order; the first row of each label (the last
    with last=True)."""
    by = {}
    for r in rows:
        lab = keep(r) and label(r, app)
        if lab and (last or lab not in by):
            by[lab] = r
    return [(lab, by[lab]) for lab in sorted(by, key=label_key)]


def by_scheme(app, rows, cols, keep=lambda r: True, block=None, last=False):
    """Table rows [block?, scheme, order, Stat...] of the reported schemes; cols: [(column, format[, scale])]."""
    return [([block] if block is not None else []) + list(scheme_order(lab, app)) + [stat_of(r, *c) for c in cols]
            for lab, r in scheme_rows(rows, app, keep, last)]


SCHEME_COLS = ["scheme", "order"]


# ---------------------------------------------------------------- density and statistics (S1 to S4)

SCHEMES_CSV = {"kraken2": "kraken2_schemes.csv", "minimap2": "minimap2_schemes.csv", "kmc": "kmc3_schemes.csv",
               "sshash": "sshash_schemes.csv"}


def expected_density_polarities(R):
    rows = []
    for app in APPS:
        rs = [r for r in R.csv(app, SCHEMES_CSV[app]) if r["scheme"] in GM_SCHEMES
              and r.get("gm_table", "") != "w11_k13_canonical" and r.get("gm_lut", "binary") == "binary"]
        keys = sorted({(int(r["k"]), int(r["w"]), r["scheme"]) for r in rs}, key=lambda t: (-t[0], list(GM_SCHEMES).index(t[2])))
        for k, w, s in keys:
            both = {r["projection"]: r for r in rs if (int(r["k"]), int(r["w"]), r["scheme"]) == (k, w, s)}
            assert sorted(both) == list(POLARITIES), (app, k, s)
            for p in POLARITIES:
                r = both[p]
                rows.append([METHOD[app], k, w, *scheme_order(GM_SCHEMES[s], app), p,
                             Stat(mean_text(r["expected_density_factor_mean"], 5), sd1(r["expected_density_factor_std"])),
                             "yes" if p == POLARITY else None])
    return Table("expected_density_polarities",
                 "All methods: expected density factor of each polarity of Phi. GreedyMini-based schemes, counted by the "
                 "selection code of each method on 10^9 uniformly random bases, a different string for each of the five "
                 "seeds. Columns: the method, k and w; the scheme and its order; the polarity of Phi (at0-cg1: {A,T} mapped to 0; "
                 "at1-cg0: {A,T} mapped to 1); the expected density factor (mean and standard deviation over the seeds); and "
                 "whether the polarity is the one reported in the paper (at0-cg1, chosen arbitrarily).",
                 ["method", "k", "w"] + SCHEME_COLS + ["polarity", "expected density factor", "reported in the paper"], rows)


# the benchmark datasets of each method, as columns of <app>_schemes.csv: (column prefix, dataset id, k)
PARTICULAR = {"kraken2": [("strain_exclusion", None, 31)],
              "minimap2": [("chm13_chrX", "chm13_chrX", 19), ("grch38_wg", "grch38_wg", 19), ("chm13_rel3_ont", "chm13_rel3_ont", 15)],
              "kmc": [(d, d, 9) for d, _ in DATASETS["kmc"]], "sshash": [(d, d, 21) for d, _ in DATASETS["sshash"]]}


def particular_density(R):
    rows = []
    for app in APPS:
        rs = R.csv(app, SCHEMES_CSV[app])
        for pre, ds, k in PARTICULAR[app]:
            col = pre + "_particular_density_factor"
            for lab, r in scheme_rows(rs, app, lambda r: int(r["k"]) == k and fl(r.get(col + "_mean")) is not None):
                rows.append([METHOD[app], dataset_name(app, ds) if ds else KRAKEN2_LIBRARY, *kw(app, ds),
                             *scheme_order(lab, app), Stat("%.3f" % fl(r[col + "_mean"]), sd1(r.get(col + "_std", "")))])
    return Table("particular_density",
                 "All methods: particular density factor of each scheme on each dataset. Columns: the method, dataset, k and "
                 "w; the scheme and its order; and the particular density factor (mean and standard deviation over the 5 "
                 "seeds; the default schemes of Minimap2 and KMC3 have no seed). The GreedyMini orders are built as in Table 2 "
                 "of the paper, with {A,T} mapped to 0 by Phi. Minimap2 uses a GreedyMini order on t-mers in the "
                 "mod-minimizer only with k=19 (CHM13 chrX, GRCh38).",
                 ["method", "dataset", "k", "w"] + SCHEME_COLS + ["particular density factor"], rows)


# The metrics of the two polarities: app -> [(runs csv, dataset of a run, runs used, [(column, name, format)])]
POLARITY_METRICS = {
    "kraken2": [
        ("kraken2_build_runs.csv", lambda r: None, lambda r: True,
         [("density_factor", "particular density factor", "{:.3f}"), ("build_wall_s", "build wall time (s)", "{:.0f}"),
          ("build_max_rss_gb", "build peak memory (GB)", "{:.2f}"), ("index_size_gb", "database size (GB)", "{:.3f}")]),
        ("kraken2_runs.csv", lambda r: r["dataset"], lambda r: r["minimum_hit_groups"] == "2",
         [("classify_run_s", "classification processing time (s)", "{:.2f}"),
          ("classify_max_rss_gb", "classification peak memory (GB)", "{:.2f}"),
          ("genus_f1", "genus F1 (%)", "{:.2f}", 100), ("species_f1", "species F1 (%)", "{:.2f}", 100)]),
    ],
    "minimap2": [
        ("minimap2_runs.csv", lambda r: r["dataset"], None,   # None: the block at the occurrence threshold of the default scheme
         [("density_factor", "particular density factor", "{:.3f}"), ("build_wall_s", "index build wall time (s)", "{:.1f}"),
          ("build_max_rss_gb", "index build peak memory (GB)", "{:.2f}"), ("index_size_gb", "index size (GB)", "{:.3f}"),
          ("map_run_s", "mapping processing time (s)", "{:.1f}"), ("map_max_rss_gb", "mapping peak memory (GB)", "{:.2f}"),
          ("frac_mapped", "mapped reads (%)", "{:.2f}", 100), ("frac_correct", "correctly mapped reads (%)", "{:.2f}", 100)]),
    ],
    "kmc": [
        ("kmc3_runs.csv", lambda r: r["dataset"], lambda r: r["gm_lut"] == "binary",
         [("density_factor", "particular density factor", "{:.3f}"), ("count_wall_s", "total wall time (s)", "{:.1f}"),
          ("count_stage1_s", "stage 1 time (s)", "{:.1f}"), ("count_stage2_s", "stage 2 time (s)", "{:.1f}"),
          ("count_max_rss_gb", "peak memory (GB)", "{:.2f}"), ("tmp_disk_gb", "temporary disk (GB)", "{:.2f}")]),
    ],
    "sshash": [
        ("sshash_runs.csv", lambda r: r["dataset"], lambda r: r["gm_table"] != "w11_k13_canonical",
         [("density_factor", "particular density factor", "{:.3f}"), ("build_wall_s", "build wall time (s)", "{:.1f}"),
          ("build_max_rss_gb", "build peak memory (GB)", "{:.2f}"), ("bits_per_kmer", "index size (bits per L-mer)", "{:.3f}"),
          ("query_run_s", "streaming query processing time (s)", "{:.2f}"),
          ("query_max_rss_gb", "streaming query peak memory (GB)", "{:.2f}"),
          ("bench_positive_ns", "positive lookup (ns per L-mer)", "{:.1f}"),
          ("bench_negative_ns", "negative lookup (ns per L-mer)", "{:.1f}")]),
    ],
}


def polarity_comparison(R):
    """Every metric of the GreedyMini-based schemes with each polarity of Phi, and the paired test of at1-cg0 / at0-cg1."""
    rows = []
    for app in APPS:
        tests = []
        for fname, dataset, use, metrics in POLARITY_METRICS[app]:
            csv_rows = R.csv(app, fname)
            if use is None:
                csv_rows, use = minimap2_blocks(csv_rows)[0], lambda r: True
            runs = collections.defaultdict(lambda: collections.defaultdict(dict))   # (dataset, scheme) -> polarity -> seed -> run
            for r in csv_rows:
                if r["scheme"] in GM_SCHEMES and use(r):
                    seeds = runs[(dataset(r), r["scheme"])][r["projection"]]
                    assert r["seed_idx"] not in seeds, ("two runs of a seed", app, fname, r["scheme"], r["seed_idx"])
                    seeds[r["seed_idx"]] = r
            for (ds, s), pol in runs.items():
                a, b = pol[POLARITIES[0]], pol[POLARITIES[1]]
                seeds = sorted(set(a) & set(b))
                for i, (col, name, fmt, *scale) in enumerate(metrics):
                    scale = scale[0] if scale else 1.0
                    va, vb = [fl(a[x].get(col)) for x in seeds], [fl(b[x].get(col)) for x in seeds]
                    if all(v is None for v in va + vb):   # the metric does not apply (e.g. no truth for the reads)
                        continue
                    if len(seeds) < 2 or not all(v and v > 0 for v in va + vb):
                        print("polarity_comparison: no test of %s %s %s %s (fewer than 2 seeds, or a value missing or 0)"
                              % (app, ds, s, col), file=sys.stderr)
                        continue
                    logs = [math.log(y / x) for x, y in zip(va, vb)]
                    tests.append({"ds": ds, "scheme": s, "order": (fname, i), "name": name, "mean_log": statistics.mean(logs),
                                  "a": stat(statistics.mean(va), statistics.stdev(va), fmt, scale),
                                  "b": stat(statistics.mean(vb), statistics.stdev(vb), fmt, scale),
                                  "test": significance.paired_test(logs)})
        families = collections.defaultdict(list)   # Bonferroni over the schemes of a dataset and metric
        for t in tests:
            families[(t["ds"], t["name"])].append(t)
        for fam in families.values():
            for t, p in zip(fam, bonferroni([t["test"][3] for t in fam])):
                t["p"] = p
        dorder = [None] + [d for d, _ in DATASETS[app]]
        files = [f for f, *_ in POLARITY_METRICS[app]]
        for t in sorted(tests, key=lambda t: (dorder.index(t["ds"]), list(GM_SCHEMES).index(t["scheme"]),
                                               files.index(t["order"][0]), t["order"][1])):
            ratio, lo, hi, _, eq = t["test"]
            rows.append([METHOD[app], dataset_name(app, t["ds"]) if t["ds"] else KRAKEN2_LIBRARY, *kw(app, t["ds"]),
                         *scheme_order(GM_SCHEMES[t["scheme"]], app), t["name"], t["a"], t["b"], "%.3f" % ratio,
                         "[%.3f, %.3f]" % (lo, hi), p_text(t["p"]), significance.verdict(t["mean_log"], t["p"], eq)])
    return Table("polarity_comparison",
                 "All methods: the two polarities of Phi. Each metric of the GreedyMini-based schemes on each dataset, with "
                 "{A,T} mapped to 0 (at0-cg1, reported in the paper) and to 1 (at1-cg0). Columns: the method, dataset, k and w; "
                 "the scheme and its order; the metric, with its unit; its mean and standard deviation over the 5 seeds with "
                 "each polarity; the ratio at1-cg0 / at0-cg1 (geometric mean over the paired seeds) and its 95% confidence "
                 "interval; the Bonferroni-corrected p-value of the paired t-test on the log ratios (multiplied by the number of "
                 "GreedyMini-based schemes on the dataset and metric); and the verdict as in the significance tests (lower or "
                 "higher if significant, else equivalent within ±5% or inconclusive). Minimap2 with the occurrence threshold "
                 "of the default scheme; KMC3 with the binary table.",
                 ["method", "dataset", "k", "w"] + SCHEME_COLS + ["metric", POLARITIES[0], POLARITIES[1],
                                                                   "ratio at1-cg0 / at0-cg1", "95% CI of the ratio",
                                                                   "p-value (Bonferroni)", "verdict"], rows)


METRIC_NAME = {"build_wall_s": "build wall time", "index_size_gb": "index size", "build_max_rss_gb": "build peak memory",
               "classify_run_s": "classification processing time", "classify_wall_s": "classification total wall time",
               "classify_max_rss_gb": "classification peak memory", "map_run_s": "mapping processing time",
               "map_wall_s": "mapping total wall time", "map_max_rss_gb": "mapping peak memory",
               "count_wall_s": "total wall time", "count_stage1_s": "stage 1 time",
               "count_stage2_s": "stage 2 time", "tmp_disk_gb": "temporary disk", "count_max_rss_gb": "peak memory",
               "query_run_s": "streaming query processing time", "query_wall_s": "streaming query total wall time",
               "query_max_rss_gb": "streaming query peak memory", "bits_per_kmer": "index size (bits per L-mer)",
               "bench_positive_ns": "positive lookup time", "bench_negative_ns": "negative lookup time",
               "density_factor": "particular density factor"}


def metric_name(app, m):
    if app == "kraken2" and m == "index_size_gb":
        return "database size"
    return METRIC_NAME.get(m, m.replace("_", " "))


# the space metrics of the figures (memory and storage); the others are runtimes
SPACE_METRICS = {"index_size_gb", "build_max_rss_gb", "classify_max_rss_gb", "map_max_rss_gb", "tmp_disk_gb",
                 "bits_per_kmer"}


def significance_tests(R, space):
    """The tests of benchmark/common/significance.py, computed from each result set, of the runtime
    metrics of the figures or (space) of their memory and storage metrics."""
    rows = []
    for app in APPS:
        tested = [r for r in significance.rows_of(str(R[app]), app)
                  if r.get("verdict") and (r["metric"] in SPACE_METRICS) == space]
        dorder = [d for d, _ in DATASETS[app]]
        for m in dict.fromkeys(r["metric"] for r in tested):
            block = []
            for r in tested:
                lab = key_label(r["scheme"], app)
                if r["metric"] != m:
                    continue
                ds = r["setting"].split("/")[0]
                ds = ds if ds in dorder else None   # Kraken2's build: the library
                block.append(((dorder.index(ds) if ds else -1, label_key(lab)),
                              [METHOD[app], metric_name(app, m), dataset_name(app, ds) if ds else KRAKEN2_LIBRARY,
                               *kw(app, ds), *scheme_order(lab, app), "%.3f" % r["ratio"],
                               "[%.3f, %.3f]" % (r["ci95_lo"], r["ci95_hi"]), p_text(r["p_bonferroni"]), r["verdict"]]))
            rows += [r for _, r in sorted(block, key=lambda b: b[0])]
    return Table("significance_space" if space else "significance_runtime",
                 "All methods: significance tests of the %s. Change of each %s metric of the figures from the default scheme "
                 "(Methods). Columns: the method, metric, dataset, k and w; the scheme and its order; the ratio of the metric "
                 "to that of the default scheme (geometric mean over the 5 paired seeds) and its 95%% confidence interval; the "
                 "Bonferroni-corrected p-value of the paired t-test on the log ratios (multiplied by the number of schemes "
                 "compared with the default scheme on the dataset and metric); and the verdict: lower or higher "
                 "(significant, adjusted p<0.05), equivalent (90%% confidence interval within ±5%%), or inconclusive."
                 % (("memory and storage usage", "space") if space else ("runtime", "runtime")),
                 ["method", "metric", "dataset", "k", "w"] + SCHEME_COLS + ["ratio to the default", "95% CI of the ratio",
                                                                            "p-value (Bonferroni)", "verdict"], rows)


# ---------------------------------------------------------------- correlations with density

CORRELATION_NOTE = {
    "kraken2": " The database is built once, from the strain-exclusion library, so its metrics have no per-dataset values; "
               "the classification is given for each dataset, at minimum hit groups 2.",
    "minimap2": " Only the runs with the occurrence threshold fixed to the value of the default scheme are paired, since the "
                "threshold derived from each index differs between schemes.",
}


def density_correlations(R, app):
    """The correlations of benchmark/common/density_effect.py, over all runs and per dataset."""
    dorder = [d for d, _ in DATASETS[app]]
    corr = []
    for m, pts in density_effect.points(str(R[app]), app):
        corr += density_effect.rows_of(app, m, pts)
    per_ds = sorted({c["setting"] for c in corr if c["setting"].split("/")[0] in dorder},
                    key=lambda s: (dorder.index(s.split("/")[0]), s))
    rows = []
    for m in dict.fromkeys(c["metric"] for c in corr):
        by = {c["setting"]: c for c in corr if c["metric"] == m}
        rows.append([metric_name(app, m), signed(by["all"]["r"]), by["all"]["n"]] +
                    [signed(by[s]["r"]) if s in by else None for s in per_ds])
    return Table("%s_density_correlations" % ("kmc3" if app == "kmc" else app),
                 "%s: correlation with the particular density. How closely each metric follows the particular density factor "
                 "(Methods). Every run of a scheme is paired with the run of the default scheme of the same dataset, setting "
                 "and seed, and gives the ratios of its particular density factor and of its metric to those of the default "
                 "scheme. Columns: the metric; the Pearson correlation r of the two ratios over the runs of all datasets, and "
                 "the number of paired runs; and r over the runs of each dataset.%s"
                 % (METHOD[app], CORRELATION_NOTE.get(app, "")),
                 ["metric", "r, all datasets", "paired runs"] + ["r, " + dataset_name(app, s.split("/")[0]) for s in per_ds],
                 rows)


# ---------------------------------------------------------------- Kraken2

def kraken2_tables(R):
    build = R.csv("kraken2", "kraken2_build.csv")
    summary = R.csv("kraken2", "kraken2_summary.csv")
    yield Table("kraken2_build",
                "Kraken2: database build, from the strain-exclusion library, 16 threads. Columns: the scheme and its order; "
                "the particular density factor; the distinct selected k-mers (millions); the build wall time and its two "
                "steps: estimating the capacity of the hash table, and processing the library (computing the selected k-mers "
                "of every sequence and inserting them into the table), which with 1 to 2 s of the other steps of "
                "kraken2-build sum to the build wall time; the build CPU time; the build peak memory; and "
                "the database size. Mean and standard deviation over the 5 seeds.",
                SCHEME_COLS + ["particular density factor", "distinct selected k-mers (10^6)", "build wall time (s)",
                               "capacity estimation wall time (s)", "library processing wall time (s)", "build CPU time (s)",
                               "build peak memory (GB)", "database size (GB)"],
                by_scheme("kraken2", build, [("density_factor", "{:.3f}"), ("distinct_minimizers", "{:.0f}", 1e-6),
                                             ("build_wall_s", "{:.0f}"), ("build_estimate_s", "{:.0f}"), ("build_db_s", "{:.0f}"),
                                             ("build_cpu_s", "{:.0f}"), ("build_max_rss_gb", "{:.2f}"), ("index_size_gb", "{:.2f}")]))

    def mhg2(ds):
        return lambda r: r["dataset"] == ds and r["minimum_hit_groups"] == "2"

    # the database load (and with it the total time) depends on whether the database was in the page cache (< 5 s)
    load = sorted(v for v in (fl(r.get("classify_load_s")) for r in R.csv("kraken2", "kraken2_runs.csv")
                              if r["minimum_hit_groups"] == "2") if v)
    cached = sum(v < 5 for v in load)

    def load_range(vs):
        return "%.1f to %.1f s in %d runs" % (vs[0], vs[-1], len(vs)) if vs else "no runs"
    rows = []
    for ds, name in DATASETS["kraken2"]:
        rows += by_scheme("kraken2", summary, [("classify_run_s", "{:.2f}"), ("classify_wall_s", "{:.2f}"),
                                               ("classify_load_s", "{:.2f}"), ("classify_max_rss_gb", "{:.2f}"),
                                               ("classified_pct", "{:.2f}"), ("classify_mreads_per_min", "{:.2f}")],
                          mhg2(ds), name)
    yield Table("kraken2_classification",
                "Kraken2: classification, one thread, minimum hit groups 2 (the default of Kraken2). Columns: the dataset; "
                "the scheme and its order; the processing time (as timed by Kraken2: classifying the reads, without loading "
                "the database); the total wall time, which is the sum of the processing time and the database load time; the "
                "database load time; the classification peak memory; the classified reads (%% of the reads); and the speed "
                "(million reads per minute of processing time). Mean and standard deviation over the 5 seeds. The load time, "
                "and with it the total time, depends on whether the database was in the page cache: over the %d runs, it took "
                "%s, and, when read from disk, %s." % (len(load), load_range(load[:cached]), load_range(load[cached:])),
                ["dataset"] + SCHEME_COLS + ["processing time (s)", "total wall time (s)", "database load time (s)",
                                             "classification peak memory (GB)", "classified reads (%)",
                                             "speed (Mreads/min)"], rows)

    prof = R.csv("kraken2", "kraken2_profile.csv")
    rows = []
    for ds, name in DATASETS["kraken2"]:
        for lab, r in scheme_rows(prof, "kraken2", mhg2(ds)):
            rows.append([name, *scheme_order(lab, "kraken2")] +
                        [stat_of(r, c, f) for c, f in (("scan_pct_of_phases", "{:.1f}"), ("lookup_pct_of_phases", "{:.1f}"),
                                                        ("resolve_pct_of_phases", "{:.1f}"), ("scan_s", "{:.3f}"),
                                                        ("lookups_per_mbp", "{:.0f}"), ("lookup_miss_pct", "{:.1f}"),
                                                        ("ns_per_lookup", "{:.1f}"))] +
                        [stat_of(r, "non_uhs_pct", "{:.1f}") if "greedymini" in r["scheme"] else None])
    yield Table("kraken2_classification_profile",
                "Kraken2: classification profile, minimum hit groups 2 (Appendix (Profiling)). The classification of the "
                "reads has three timed steps: computing the selected k-mer of every L-window, looking the selected k-mers up "
                "in the hash table, and resolving the hits into a taxon; they take 57% to 76% of the processing time (Table "
                "S6), and the rest is reading the reads and writing the classification of each read. Columns: the dataset; "
                "the scheme and its order; the shares of the three steps (% of the time of the three steps, summing to 100); "
                "the time of computing the selected k-mers; the lookups per Mbp of reads; the lookups that miss the table (% "
                "of the lookups); the time per lookup; and, for the GreedyMini orders, the k-mers of the reads whose rank is "
                "infinite, i.e., outside the UHS (% of the k-mers). Mean and standard deviation over the 5 seeds.",
                ["dataset"] + SCHEME_COLS + ["computing selected k-mers (% of steps)", "lookups (% of steps)",
                                             "resolving hits (% of steps)", "computing selected k-mers (s)",
                                             "lookups per Mbp", "lookup misses (% of lookups)", "time per lookup (ns)",
                                             "k-mers outside the UHS (%)"], rows)


def kraken2_accuracy_tables(R):
    for name, fname, groups, caption in (
            ("kraken2_accuracy_strain_exclusion", "kraken2_accuracy_strain_exclusion.csv", ("all", "prokaryote", "virus"),
             "Kraken2: accuracy on the strain-exclusion reads. Simulated reads (Mason, held-out genomes), minimum hit groups "
             "2 (the default of Kraken2). Columns: the block (genome group and rank); the scheme and its order; and the "
             "sensitivity, the positive predictive value (PPV) and their harmonic mean F1, all in %, as the mean and standard "
             "deviation over the 5 seeds."),
            ("kraken2_accuracy_cami", "kraken2_accuracy_cami.csv", ("all", "prokaryote", "bacteria", "archaea", "virus"),
             "Kraken2: accuracy on the CAMI II reads. Simulated reads, minimum hit groups 2. Columns: the block (genome "
             "group and rank); the scheme and its order; and the sensitivity, the positive predictive value (PPV) and their "
             "harmonic mean F1, all in %, as the mean and standard deviation over the 5 seeds.")):
        acc = R.csv("kraken2", fname)
        rows = []
        for g in groups:
            for rank in ("genus", "species"):
                keep = (lambda g, rank: lambda r: r["group"] == g and r["rank"] == rank and r["minimum_hit_groups"] == "2")(g, rank)
                rows += by_scheme("kraken2", acc, [("sensitivity", "{:.2f}", 100), ("ppv", "{:.2f}", 100), ("f1", "{:.2f}", 100)], keep,
                                  "%s, %s level" % (g, rank))
        yield Table(name, caption, ["block"] + SCHEME_COLS + ["sensitivity (%)", "PPV (%)", "F1 (%)"], rows)

    # FDA-ARGOS: each of the 25 experiments (SRA run), scored against its SRA-assigned taxon
    fda = [r for r in R.csv("kraken2", "kraken2_accuracy_fda_argos.csv") if r["scope"] == "experiment"]
    experiments = sorted({(r["sra_run"], r["species"]) for r in fda})
    assert len(experiments) == 25, len(experiments)
    fractions = ("concordant", "discordant", "ancestor", "unclassified")
    rows = []
    for sra, species in experiments:
        rows += by_scheme("kraken2", fda, [(f + "_fraction", "{:.2f}", 100) for f in fractions],
                          lambda r: r["sra_run"] == sra and r["minimum_hit_groups"] == "2", "%s, %s" % (sra, species), last=True)
    yield Table("kraken2_accuracy_fda_argos",
                "Kraken2: accuracy on the FDA-ARGOS reads. Genus-level calls of each of the 25 FDA-ARGOS experiments, scored "
                "separately against the SRA-assigned taxon, at minimum hit groups 2, the default of Kraken2. Columns: the "
                "block (SRA run and species); the scheme and its order; and the reads (% of the reads of the experiment, "
                "summing to 100) called as the assigned genus or a descendant (concordant), off its lineage (discordant), as "
                "an ancestor of the genus (ancestor), or unclassified, as the mean and standard deviation over the 5 seeds.",
                ["block"] + SCHEME_COLS + [f + " (%)" for f in fractions], rows)


# ---------------------------------------------------------------- Minimap2

def aggregate(rows, app, metrics):
    """{(dataset, label): {metric: (mean, std)}} of runs; metrics: {name: f(run) -> value or None}."""
    vals = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        lab = label(r, app)
        if lab is None:
            continue
        for m, f in metrics.items():
            v = f(r)
            if v is not None:
                vals[(r["dataset"], lab)][m].append(v)
    return {k: {m: (statistics.mean(v), statistics.stdev(v) if len(v) > 1 else 0.0) for m, v in d.items()}
            for k, d in vals.items()}


def aggregate_rows(agg, app, cols, kw_cols=False):
    """Table rows [dataset, (k, w,) scheme, order, Stat...] of aggregate(); cols: [(metric, format)], a missing metric None."""
    rows = []
    for ds, name in DATASETS[app]:
        for lab in sorted({lab for d, lab in agg if d == ds}, key=label_key):
            v = agg[(ds, lab)]
            rows.append([name] + (list(kw(app, ds)) if kw_cols else []) + list(scheme_order(lab, app)) +
                        [stat(*v[m], fmt) if m in v else None for m, fmt in cols])
    return rows


def column(c):
    return lambda r: fl(r.get(c, ""))


def minimap2_tables(R, mmt):
    fixed, derived = minimap2_blocks(R.csv("minimap2", "minimap2_summary.csv"))
    fixed_note = "occurrence threshold %s (that of the default scheme)"

    # the index build (collecting the selected k-mers of the reference, from the build log), and the mapping profile
    # (-DMM_PHASE_PROF build, one thread, the first 3,000 reads, <run>/map_profile.log); runs at the default's threshold
    runs = minimap2_blocks(R.csv("minimap2", "minimap2_runs.csv"))[0]
    agg = aggregate(runs, "minimap2", {c: column(c) for c in ("build_collect_s", "build_collect_pct", "build_cpu_s")})
    phases = []
    for log in sorted(R["minimap2"].glob("runs/*/*_f*/**/map_profile.log")):
        meta = json.loads((log.parent / "meta.json").read_text())
        ph = re.search(r"sketch_s=(\S+) seed_s=(\S+) chain_s=(\S+) align_s=(\S+)", log.read_text())
        if ph and str(meta["mid_occ_cap"]) == mmt[meta["dataset"]]:
            phases.append({**{k: "" if meta.get(k) is None else str(meta[k])
                              for k in ("dataset", "scheme", "gm_k", "projection", "gm_table", "mod_t")},
                           **dict(zip(("sketch_s", "seed_s", "chain_s", "align_s"), ph.groups()))})
    prof = aggregate(phases, "minimap2", {c: column(c) for c in ("sketch_s", "seed_s", "chain_s", "align_s")})
    for r in R.csv("minimap2", "analysis_mechanism.csv"):   # seed 0 of each scheme
        lab = label(r, "minimap2")
        if lab and r["mid_occ_cap"] == mmt[r["dataset"]]:
            prof.setdefault((r["dataset"], lab), {}).update(
                {"seeds_per_kb": (fl(r["seeds_per_kb"]), 0.0), "anchors_per_read": (fl(r["anchors_per_read_mean"]) / 1e3, 0.0)})
    for k in agg:
        agg[k].update(prof.get(k, {}))
    yield Table("minimap2_build_mapping_profile",
                "Minimap2: index build and mapping profile (Appendix (Profiling)), with the occurrence threshold fixed to that "
                "of the default scheme. Columns: the dataset, k and w; the scheme and its order. Index build (the runs of "
                "Table S9, 16 threads): the time until Minimap2 has collected the selected k-mers of the reference (from its "
                "log; reading the reference and computing its selected k-mers, mostly in one thread) and its share of the "
                "build wall time, the rest being sorting the selected k-mers and building the hash tables; and the build CPU "
                "time. Mapping (a separate profiling run with one thread, on the first 3,000 reads of each dataset, all 1,000 "
                "for the ONT reads): the time of its four steps, which sum to the profiled mapping time: computing the "
                "selected k-mers of the reads, looking up the seeds in the index, chaining the anchors, and aligning the "
                "chains. Mean and standard deviation over the 5 seeds. And, for seed 0 of the scheme, the seeds per kb of "
                "reads and the anchors per read (thousands).",
                ["block", "k", "w"] + SCHEME_COLS + ["collecting selected k-mers (s)", "collecting (% of build wall time)",
                                                     "build CPU time (s)", "computing selected k-mers (s)",
                                                     "seed lookup (s)", "chaining (s)", "alignment (s)",
                                                     "seeds per kb of reads", "anchors per read (10^3)"],
                aggregate_rows(agg, "minimap2", [("build_collect_s", "{:.1f}"), ("build_collect_pct", "{:.1f}"),
                                                 ("build_cpu_s", "{:.0f}"), ("sketch_s", "{:.2f}"), ("seed_s", "{:.2f}"),
                                                 ("chain_s", "{:.2f}"), ("align_s", "{:.2f}"), ("seeds_per_kb", "{:.1f}"),
                                                 ("anchors_per_read", "{:.1f}")], kw_cols=True))

    cost = [("density_factor", "{:.3f}"), ("build_wall_s", "{:.2f}"), ("map_run_s", "{:.1f}"), ("build_max_rss_gb", "{:.2f}"),
            ("map_max_rss_gb", "{:.2f}"), ("index_size_gb", "{:.3f}")]
    rows = []
    for ds, name in DATASETS["minimap2"]:
        k, w = kw("minimap2", ds)
        block = "%s, occurrence threshold fixed to that of the default scheme" % name
        for lab, r in scheme_rows(fixed, "minimap2", lambda r: r["dataset"] == ds):
            rows.append([block, k, w, *scheme_order(lab, "minimap2"), mmt[ds]] + [stat_of(r, *c) for c in cost])
        block = "%s, default occurrence threshold (derived from each index)" % name
        for lab, r in scheme_rows(derived, "minimap2", lambda r: r["dataset"] == ds):
            rows.append([block, k, w, *scheme_order(lab, "minimap2"), r["mid_occ_threshold"]] + [stat_of(r, *c) for c in cost])
    yield Table("minimap2_build_mapping",
                "Minimap2: index build and mapping, 16 threads, with the occurrence threshold fixed to the value of the default "
                "scheme or at its default, which Minimap2 derives from each index. Columns: the block (dataset and threshold "
                "setting), k and w; the scheme and its order; the occurrence threshold; the particular density factor; the "
                "index build wall time; the mapping processing time (the mapping wall time minus the index load time that "
                "Minimap2 reports); the peak memory of the index build and of the mapping; and the index size. Mean and "
                "standard deviation over the 5 seeds (the default scheme has no seed).",
                ["block", "k", "w"] + SCHEME_COLS + ["occurrence threshold", "particular density factor",
                                                     "index build wall time (s)", "mapping processing time (s)",
                                                     "index build peak memory (GB)", "mapping peak memory (GB)",
                                                     "index size (GB)"], rows)

    def at_default_threshold(datasets, cols):
        rows = []
        for ds, name in datasets:
            rows += by_scheme("minimap2", fixed, cols, lambda r: r["dataset"] == ds, "%s, %s" % (name, fixed_note % mmt[ds]))
        return rows

    yield Table("minimap2_accuracy_simulated",
                "Minimap2: accuracy on the simulated reads, with the occurrence threshold fixed to that of the default scheme. "
                "Columns: the block (dataset and threshold); the scheme and its order; the reads (% of the reads) mapped, "
                "mapped correctly (the primary alignment on the true contig, overlapping the true interval by at least 10% of "
                "their union), and mapped wrongly (mapped but not correctly), so that mapped = correct + wrong; the reads from "
                "repeats mapped correctly (% of the reads from repeats); the reads mapped correctly with an overlap of at least "
                "50% (% of the reads); and the mean identity of the primary alignments (%, 1 minus the per-base divergence). "
                "Mean and standard deviation over the 5 seeds (the default scheme has no seed).",
                ["block"] + SCHEME_COLS + ["mapped (%)", "correct (%)", "wrong (%)", "correct in repeats (% of repeat reads)",
                                           "correct with >=50% overlap (%)", "mean identity (%)"],
                at_default_threshold(DATASETS["minimap2"][:2],
                                     [("frac_mapped", "{:.2f}", 100), ("frac_correct", "{:.2f}", 100),
                                      ("frac_incorrectly_mapped", "{:.3f}", 100), ("frac_correct_in_repeat", "{:.1f}", 100),
                                      ("frac_correct_strict50", "{:.2f}", 100), ("mean_identity", "{:.2f}", 100)]))
    yield Table("minimap2_accuracy_ont",
                "Minimap2: accuracy on the ONT reads: 1,000 real reads of CHM13 rel3 mapped to GRCh38, without truth, with the "
                "occurrence threshold fixed to that of the default scheme. Columns: the block (dataset and threshold); the "
                "scheme and its order; the mapped reads (% of the reads); and the mean identity of the primary alignments (%, "
                "1 minus the per-base divergence). Mean and standard deviation over the 5 seeds (the default scheme has no "
                "seed).",
                ["block"] + SCHEME_COLS + ["mapped (%)", "mean identity (%)"],
                at_default_threshold(DATASETS["minimap2"][2:], [("frac_mapped", "{:.2f}", 100), ("mean_identity", "{:.2f}", 100)]))


# ---------------------------------------------------------------- KMC3

def kmc3_tables(R):
    summary = R.csv("kmc", "kmc3_summary.csv")
    rows = []
    for ds, name in DATASETS["kmc"]:
        rows += by_scheme("kmc", summary, [("density_factor", "{:.3f}"), ("n_super_kmers", "{:.1f}", 1e-6),
                                           ("avg_super_kmer_len", "{:.2f}"), ("count_stage1_s", "{:.1f}"),
                                           ("count_stage2_s", "{:.1f}"), ("count_wall_s", "{:.1f}"), ("count_cpu_s", "{:.0f}"),
                                           ("tmp_disk_gb", "{:.1f}")], lambda r: r["dataset"] == ds, name)
    yield Table("kmc3_counting",
                "KMC3: counting, L=28, 1024 bins, 24 threads. Columns: the dataset; the scheme and its order; the particular "
                "density factor; the superwindows (millions); the mean superwindow length (L-mers per superwindow); the time "
                "of stage 1 (computing the selected k-mers and distributing the superwindows into bins) and of stage 2 "
                "(counting the L-mers of each bin), as reported by KMC3, which sum to the total wall time; the total wall time; "
                "the CPU time; and the temporary disk used by the bins. Mean and standard deviation over the 5 seeds (the "
                "default scheme has no seed). The peak memory is not listed, because the memory of KMC3 is limited by design "
                "(12 GB, 32 GB for H. sapiens 3) and is the same for all schemes.",
                ["block"] + SCHEME_COLS + ["particular density factor", "superwindows (10^6)",
                                           "mean superwindow length (L-mers)", "stage 1 time (s)", "stage 2 time (s)",
                                           "total wall time (s)", "CPU time (s)", "temporary disk (GB)"], rows)

    # the cost of the order (tools/bench_signature_impls), the GreedyMini order with its binary table
    cost = R.csv("kmc", "kmc3_signature_cost.csv")
    rows = []
    for ds, name in DATASETS["kmc"]:
        variants = {r["variant"].replace(" (binary table)", ""): r for r in cost if r["dataset"] == ds and "DNA" not in r["variant"]}
        for lab in sorted(variants, key=label_key):
            r = variants[lab]
            rows.append([name, *scheme_order(lab, "kmc"),
                         stat(fl(r["signature_ns_per_base_mean"]), fl(r["signature_ns_per_base_std"]), "{:.2f}"),
                         stat(fl(r["window_ns_per_base_mean"]), fl(r["window_ns_per_base_std"]), "{:.2f}"),
                         "%.3f" % fl(r["window_ns_per_base_ratio"]), "%.1f" % fl(r["superwindows_per_kb"])])
    yield Table("kmc3_order_cost_profile",
                "KMC3: order cost profile (Appendix (Profiling)), on one core, on 200 Mbp of the reads of each dataset. "
                "Columns: the dataset; the scheme and its order; the time of ranking every k-mer (ns per base); the time of "
                "ranking every k-mer and selecting the k-mer of each L-window (ns per base), which includes the ranking; its "
                "ratio to that of the default order; and the superwindows per kb of reads. Mean and standard deviation over the "
                "5 seeds, each the median of 7 repetitions.",
                ["block"] + SCHEME_COLS + ["ranking (ns per base)", "ranking and selection (ns per base)",
                                           "ranking and selection, ratio to the default", "superwindows per kb"], rows)


# ---------------------------------------------------------------- SSHash

def sshash_tables(R):
    summary = R.csv("sshash", "sshash_summary.csv")
    rows = []
    for ds, name in DATASETS["sshash"]:
        for lab, r in scheme_rows(summary, "sshash", lambda r: r["dataset"] == ds):
            rows.append([name, *scheme_order(lab, "sshash")] +
                        [stat_of(r, *c) for c in (("density_factor", "{:.3f}"), ("distinct_minimizers", "{:.1f}", 1e-6),
                                                  ("build_wall_s", "{:.1f}"), ("build_max_rss_gb", "{:.2f}"),
                                                  ("bits_per_kmer", "{:.3f}"))] +
                        [stat_of(r, "bits_per_kmer", "{:.2f}", fl(r["n_kmers"]) / 8e9)])
    yield Table("sshash_build",
                "SSHash: index build, k=21, w=11, 16 threads. Columns: the collection; the scheme and its order (the GreedyMini "
                "order built for 15-mers, center-extended to 21-mers); the particular density factor; the distinct selected "
                "k-mers (millions); the build wall time; the build peak memory; and the index size, in bits per L-mer and in "
                "GB. Mean and standard deviation over the 5 seeds.",
                ["block"] + SCHEME_COLS + ["particular density factor", "distinct selected k-mers (10^6)",
                                           "build wall time (s)", "build peak memory (GB)", "index size (bits per L-mer)",
                                           "index size (GB)"], rows)

    rows = []
    for ds, name in DATASETS["sshash"]:
        rows += by_scheme("sshash", summary, [("bench_positive_ns", "{:.0f}"), ("bench_negative_ns", "{:.0f}"),
                                              ("bench_access_ns", "{:.0f}"), ("query_run_s", "{:.2f}")],
                          lambda r: r["dataset"] == ds, name)
    yield Table("sshash_lookups_queries",
                "SSHash: isolated lookups and streaming queries, one thread. Columns: the collection; the scheme and its order; "
                "the time of an isolated lookup of an L-mer present in the collection (positive) and of one absent from it "
                "(negative), and of accessing an L-mer by its identifier, in ns per L-mer; and the processing time of the "
                "streaming query of the real reads (as timed by SSHash, without loading the index). Mean and standard "
                "deviation over the 5 seeds.",
                ["block"] + SCHEME_COLS + ["positive lookup (ns per L-mer)", "negative lookup (ns per L-mer)",
                                           "access (ns per L-mer)", "streaming query processing time (s)"], rows)

    # the streaming queries (-DSSHASH_PROFILE build): shares of the L-mers of the reads
    def share(c):
        return lambda r: 100 * fl(r[c]) / fl(r["query_n_kmers"]) if fl(r.get(c)) is not None and fl(r.get("query_n_kmers")) else None

    def ns_per_kmer(r):
        t, n = fl(r.get("query_run_s")), fl(r.get("query_n_kmers"))
        return 1e9 * t / n if t and n else None

    agg = aggregate(R.csv("sshash", "sshash_runs.csv"), "sshash",
                    {"neg": share("query_n_negative_kmers"), "search": share("query_n_searches"),
                     "ext": share("query_n_extensions"),
                     "full": share("query_prof_full_lookups"), "cache": share("query_n_bucket_cache_hits"),
                     "pos": column("query_prof_positions_per_full_lookup"), "ns": ns_per_kmer})
    yield Table("sshash_streaming_query_profile",
                "SSHash: streaming query profile (Appendix (Profiling)), one thread. Each L-mer of the reads is either "
                "absent from the collection, or present and found by a lookup, or present and found by extending the previous "
                "L-mer along the stored string, without a lookup; these three shares sum to 100. A lookup is either a full "
                "lookup of the bucket of the selected k-mer or, when the selected k-mer is that of the previous lookup, a check "
                "of the positions of that bucket, which are cached. Columns: the collection; the scheme and its order; the "
                "L-mers absent, found by a lookup, and found by extension (% of the L-mers); the full lookups and the cached "
                "checks (% of the L-mers); the positions of the selected k-mer scanned per full lookup; and the processing time "
                "per L-mer. Mean and standard deviation over the 5 seeds.",
                ["block"] + SCHEME_COLS + ["absent (% of L-mers)", "found by lookup (% of L-mers)",
                                           "found by extension (% of L-mers)", "full lookups (% of L-mers)",
                                           "cached checks (% of L-mers)", "positions per full lookup",
                                           "processing time (ns per L-mer)"],
                aggregate_rows(agg, "sshash", [("neg", "{:.1f}"), ("search", "{:.2f}"), ("ext", "{:.1f}"), ("full", "{:.2f}"),
                                               ("cache", "{:.2f}"), ("pos", "{:.2f}"), ("ns", "{:.1f}")]))

    # buckets: the k-mers in singleton, light and heavy buckets (per seed), and the positions per selected k-mer
    kinds = ("singleton", "light", "heavy")
    hist = R.csv("sshash", "sshash_bucket_histogram.csv")
    buckets = R.csv("sshash", "sshash_buckets.csv")
    rows = []
    for ds, name in DATASETS["sshash"]:
        per_seed = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(float)))
        for r in hist:
            lab = r["dataset"] == ds and label(r, "sshash")
            if lab:
                per_seed[lab][r["seed_idx"]][(r["bucket_type"], "kmers")] += float(r["n_kmers_in_buckets"])
                per_seed[lab][r["seed_idx"]][(r["bucket_type"], "buckets")] += float(r["n_buckets"])
        occurrences = dict(scheme_rows(buckets, "sshash", lambda r: r["dataset"] == ds))
        for lab in sorted(per_seed, key=label_key):
            seeds = list(per_seed[lab].values())

            def mean_sd(v):
                return statistics.mean(v), statistics.stdev(v) if len(v) > 1 else 0.0

            n = [sum(s[(t, "buckets")] for t in kinds) / 1e6 for s in seeds]
            shares = [mean_sd([100 * s[(t, "kmers")] / sum(s[(u, "kmers")] for u in kinds) for s in seeds]) for t in kinds]
            rows.append([name, *scheme_order(lab, "sshash"), stat(*mean_sd(n), "{:.2f}"),
                         stat_of(occurrences[lab], "minimizer_positions_per_minimizer", "{:.3f}")] +
                        [stat(m, sd, "{:.1f}") for m, sd in shares])
    yield Table("sshash_buckets",
                "SSHash: buckets. A bucket holds the positions of one selected k-mer. Columns: the collection; the scheme and "
                "its order; the buckets (millions); the occurrences per selected k-mer (positions per bucket); and the L-mers "
                "in singleton, light and heavy buckets (% of the L-mers, summing to 100). Mean and standard deviation over the "
                "5 seeds.",
                ["block"] + SCHEME_COLS + ["buckets (10^6)", "occurrences per selected k-mer", "singleton buckets (% of L-mers)",
                                           "light buckets (% of L-mers)", "heavy buckets (% of L-mers)"], rows)

    # where the space is saved: bits per L-mer of each structure, saved relative to the default scheme
    breakdown = R.csv("sshash", "sshash_space_breakdown.csv")
    parts = (("tags",), ("mphf",), ("strings", "strings_offsets"), ("skew_fv", "skew_H"), ("light_L", "other"))
    rows = []
    for ds, name in DATASETS["sshash"]:
        for lab, r in scheme_rows(breakdown, "sshash", lambda r: r["dataset"] == ds):
            if lab != "default":
                shares = ["%.1f" % (sum(fl(r["share_%s_pct" % c]) for c in cs) + 0.0) for cs in parts]
                rows.append([name, *scheme_order(lab, "sshash"), "%.3f" % fl(r["bits_total"]), "%.3f" % fl(r["saved_total"])] +
                            [s.replace("-0.0", "0.0") for s in shares])
    yield Table("sshash_index_space_breakdown",
                "SSHash: index space breakdown, mean over the 5 seeds. Columns: the collection; the scheme and its order; the "
                "index size (bits per L-mer); the bits per L-mer saved relative to the default scheme (negative: the index "
                "grew); and the share of the saved bits of each structure (% of the saved bits, summing to 100): the tags of "
                "the selected k-mers, the minimal perfect hash function (MPHF), the strings and their offsets, the skew index, "
                "and the rest; a negative share is a structure that grew.",
                ["block"] + SCHEME_COLS + ["index size (bits per L-mer)", "saved (bits per L-mer)", "tags (% of saved)",
                                           "MPHF (% of saved)", "strings (% of saved)", "skew index (% of saved)",
                                           "other (% of saved)"], rows)


# ---------------------------------------------------------------- Table 2 of the paper

# its rows: (app, k), and the row spacing after each (a group of rows per method)
SEEDLESS_DEFAULT = {"minimap2", "kmc"}
DENSITY_ROWS = [("kraken2", 31, "\\\\[1.5mm]"), ("minimap2", 19, "\\\\"), ("minimap2", 15, "\\\\[1.5mm]"), ("kmc", 9, "\\\\[1.5mm]"),
                ("sshash", 21, "\\\\")]


def density_table(R):
    """The rows of Table 2 (LaTeX): the expected density factor (mean, and the standard deviation over the seeds to one
    significant digit) of each scheme; in bold the lowest of a row."""
    lines = []
    for app, k, end in DENSITY_ROWS:
        rows = R.csv(app, SCHEMES_CSV[app])
        got = dict(scheme_rows(rows, app, lambda r: int(r["k"]) == k))
        if app == "kraken2":   # its default order is xor: in the xor column
            got["xor"] = got.pop("default")
        means = {lab: fl(r["expected_density_factor_mean"]) for lab, r in got.items()}
        best = min(means, key=means.get)
        cells = []
        for lab in LABELS:
            if app == "kraken2" and lab == "default":
                cells.append("\\multicolumn{3}{c}{= \\texttt{xor}}")
                continue
            if lab not in got:
                cells.append("\\multicolumn{3}{c}{--}")
                continue
            whole, frac = ("%.3f" % means[lab]).split(".")
            if lab == best:
                whole, frac = "\\mathbf{%s}" % whole, "\\mathbf{%s}" % frac
            # the default schemes of Minimap2 and KMC3 have no seed: their spread over the seeds' strings is not reported
            sd = "" if lab == "default" and app in SEEDLESS_DEFAULT else sd1(got[lab]["expected_density_factor_std"])
            cells.append("$%s$ & $%s$ & %s" % (whole, frac, "${}\\pm%s$" % sd if sd else ""))
        lines.append("%s & %d & %d & %s %s" % (METHOD[app], k, kw(app)[1], " & ".join(cells), end))
    return lines


# ---------------------------------------------------------------- all tables

# the supplementary tables, numbered S1, S2, ... in the order of their first citation in the paper
TABLES = ["expected_density_polarities", "particular_density", "polarity_comparison", "significance_runtime",
          "kraken2_build", "kraken2_classification", "kraken2_classification_profile",
          "minimap2_build_mapping_profile", "minimap2_build_mapping",
          "kmc3_density_correlations", "kmc3_counting", "kmc3_order_cost_profile",
          "sshash_build", "sshash_streaming_query_profile", "sshash_lookups_queries", "significance_space",
          "kraken2_density_correlations", "minimap2_density_correlations", "sshash_density_correlations", "sshash_buckets",
          "sshash_index_space_breakdown",
          "kraken2_accuracy_strain_exclusion", "kraken2_accuracy_cami", "kraken2_accuracy_fda_argos",
          "minimap2_accuracy_simulated", "minimap2_accuracy_ont"]


def all_tables(R):
    mmt = default_thresholds(R)
    tables = [expected_density_polarities(R), particular_density(R), polarity_comparison(R),
              significance_tests(R, space=False), significance_tests(R, space=True)]
    tables += list(kraken2_tables(R)) + list(kraken2_accuracy_tables(R)) + list(minimap2_tables(R, mmt))
    tables += list(kmc3_tables(R)) + list(sshash_tables(R)) + [density_correlations(R, app) for app in APPS]
    by_name = {t.name: t for t in tables}
    assert sorted(by_name) == sorted(TABLES), set(by_name) ^ set(TABLES)
    return [by_name[n] for n in TABLES]


README_BEGIN, README_END = "<!-- tables:begin -->", "<!-- tables:end -->"


# the folder of the supplementary tables on GitHub, which each label of the paper links to
CSV_URL = "https://github.com/OrensteinLab/MinimizerApps/blob/main/supplementary_tables/"


def write(tables, out, labels=None):
    """The CSVs (replacing the old ones), the list of the tables in out/README.md, and the LaTeX labels."""
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("S[0-9][0-9]_*.csv"):
        f.unlink()
    for i, t in enumerate(tables, 1):
        with open(out / ("S%02d_%s.csv" % (i, t.name)), "w", newline="") as f:
            csv.writer(f, lineterminator="\n").writerows(t.csv_rows())
    readme = out / "README.md"
    if readme.is_file() and README_BEGIN in readme.read_text():
        text = readme.read_text()
        items = ["- **S%d** (`S%02d_%s.csv`): %s" % (i, i, t.name, t.caption) for i, t in enumerate(tables, 1)]
        a, b = text.index(README_BEGIN) + len(README_BEGIN), text.index(README_END)
        readme.write_text(text[:a] + "\n" + "\n".join(items) + "\n" + text[b:])
    else:
        print("no %s in %s: the list of the tables is not written" % (README_BEGIN, readme))
    if labels:
        Path(labels).write_text("".join("\\newlabel{tab:%s}{{\\href{%sS%02d_%s.csv}{S%d}}{}}\n" % (t.name, CSV_URL, i, t.name, i)
                                        for i, t in enumerate(tables, 1)))
    print("%d tables in %s" % (len(tables), out))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", nargs=4, required=True, metavar="APP=DIR", help="the result set of each of " + ", ".join(APPS))
    ap.add_argument("--out", type=Path, default=REPO / "supplementary_tables", help="folder of the CSVs and their README.md")
    ap.add_argument("--labels", help="also write the LaTeX labels of the tables to this file")
    ap.add_argument("--density-table", help="also write the rows of Table 2 of the paper (LaTeX) to this file")
    a = ap.parse_args(argv)
    try:
        R = ResultSets(a.results)
    except ValueError as e:
        ap.error(str(e))
    write(all_tables(R), a.out, a.labels)
    if a.density_table:
        Path(a.density_table).write_text("\n".join(density_table(R)) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
