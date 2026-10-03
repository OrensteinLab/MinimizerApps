#!/usr/bin/env python3
"""Kraken 2 benchmark export layer: sheets, CSVs and the results workbook.

Same layout as the other apps (benchmark/sshash/scripts/export.py,
benchmark/kmc/scripts/export.py) -- one workbook, `kraken2_results.xlsx`, and
one CSV per sheet, `kraken2_<sheet>.csv`:

  overview   column reference / legend (the role `index` has in the other apps)
  libraries  one row per Kraken library (db size, taxonomy, role, query sets run)
  datasets   one row per query set (reads, read length, truth, provenance)
  orders     one row per minimizer scheme/order (CLI, l/k/w, spaced seed, density)
  summary    library x dataset x order, mean/std over seeds + ratio to Random
  runs       one row per seed x library x dataset x order (raw per-run numbers)
  ...then the detail sheets: vs_reference, timing_correlation, classify_repeats,
  strain_exclusion_{summary,runs}, fda_argos_{summary,runs}, cami_{summary,runs},
  mock_{summary,runs}.

compare_benchmark.py collects the raw numbers of one library (one
`run_benchmark.sh` invocation) and calls `export_library()`, which writes that
library's tables under <run_root>/tables/ and its workbook
<run_root>/benchmark_comparison.xlsx; `export_results()` then merges every
library under <base>/runs/*/tables/ into <base>/results/kraken2_results.xlsx
(+ kraken2_*.csv). The summary-level sheets are always *derived* from the base
tables (`derive_tables`), so a per-library export and the merged export use
the same code.

Python 3.6 compatible (run_benchmark.sh calls `python3`).
"""

from __future__ import print_function

import csv
import math
import os
import re
import struct
import sys
from collections import OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
_COMMON = os.path.join(os.path.dirname(HERE), "common")
for _p in (HERE, _COMMON):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from bench_const import REF_SCHEME, SCHEME_ORDER  # noqa: E402
from xlsx_book import write_book  # noqa: E402
from xlsx_util import write_csv  # noqa: E402

PREFIX = "kraken2_"
WORKBOOK = "kraken2_results.xlsx"
CORE_SHEETS = ("overview", "libraries", "datasets", "orders", "summary", "runs")
DETAIL_SHEETS = (
    "vs_reference",
    "timing_correlation",
    "classify_repeats",
    "strain_exclusion_summary",
    "strain_exclusion_runs",
    "fda_argos_summary",
    "fda_argos_runs",
    "cami_summary",
    "cami_runs",
    "mock_summary",
    "mock_runs",
)
ALL_SHEETS = CORE_SHEETS + DETAIL_SHEETS
# Base tables: what a per-library export persists and a merge reads back.
BASE_SHEETS = (
    "runs",
    "classify_repeats",
    "strain_exclusion_runs",
    "fda_argos_runs",
    "cami_runs",
    "mock_runs",
    "orders",
    "datasets",
)

# Pre-export-layout per-library CSV names (and binning leftovers) that this
# layer replaces; removed from <run_root>/tables/ when the new tables are written.
LEGACY_TABLES = (
    "contents.csv", "query_sets.csv", "scheme_summary.csv", "per_seed.csv",
    "vs_ref.csv", "vs_random.csv", "classify_repeats.csv",
    "accuracy_strain_exclusion.csv", "accuracy_strain_exclusion_summary.csv",
    "accuracy_fda_argos.csv", "accuracy_fda_argos_summary.csv",
    "accuracy_cami.csv", "accuracy_cami_summary.csv",
    "accuracy_mock.csv", "accuracy_mock_summary.csv",
)

LIBRARY_ORDER = ("standard", "standard_p10", "standard_p0.5", "strain_exclusion")
DATASET_ORDER = {
    "strain-exclusion": 0, "simulated": 0, "labeled": 0,
    "fda-argos": 1, "argos": 1, "cami": 2, "mock": 3, "speed": 4,
}
DOMAIN_ORDER = {"all": 0, "prokaryote": 1, "virus": 2, "viral": 2}
RANK_ORDER = {"genus": 0, "species": 1}

MOCK_HEADLINE_MIN_FRAC = 1e-4
TEXT_FIELDS = frozenset([
    "library", "dataset", "order", "scheme", "run_id", "seed", "context", "domain",
    "rank", "scope", "group", "sra_run", "species", "sources", "path", "table",
    "cli", "description", "sketching", "rank_order", "spaced_seed_mask", "role",
    "taxonomy", "db_path", "run_root", "query_sets_run", "orders_run", "libraries",
    "status", "truth", "provenance", "paired", "read_length_bp", "fragments_per_source",
    "metric", "sheet", "column", "meaning", "category", "top_off_target",
    "library_fraction", "flag", "fasta_dir", "n_sequences_note",
])


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def mean_std(values):
    vals = [float(v) for v in values if v is not None and v != ""]
    if not vals:
        return None, None
    n = len(vals)
    mu = sum(vals) / float(n)
    if n < 2:
        return mu, 0.0
    var = sum((x - mu) ** 2 for x in vals) / float(n - 1)
    return mu, math.sqrt(var)


def round_sig(x, sig=6):
    """Round a mean to `sig` significant digits (keeps tiny values readable)."""
    if x is None:
        return None
    if x == 0 or isinstance(x, int):
        return x
    return float("%.*g" % (sig, x))


def round_std(x):
    """Std rounded to >= 4 decimals and >= 4 significant digits (a 1e-5 spread
    over seeds must not collapse to 0.0) -- same rule as sshash's export."""
    if x is None:
        return None
    if x == 0:
        return 0.0
    return round(x, max(4, 3 - int(math.floor(math.log10(abs(x))))))


def clean_value(v):
    """Strip float noise (10 significant digits); NaN/inf -> None."""
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return None
        if v == 0:
            return 0.0
        return float("%.10g" % v)
    return v


def typed(v, field=None):
    """CSV text -> None / int / float / str (text-only fields stay strings)."""
    if v is None or v == "":
        return None
    if field in TEXT_FIELDS:
        return v
    try:
        return int(v)
    except ValueError:
        pass
    try:
        f = float(v)
        if math.isnan(f) or math.isinf(f):
            return v
        return f
    except ValueError:
        return v


def read_csv(path):
    """-> (fields, rows) with typed cells; ([], []) when the file is missing."""
    if not os.path.isfile(path):
        return [], []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = list(reader.fieldnames or [])
        rows = [OrderedDict((k, typed(r.get(k), k)) for k in fields) for r in reader]
    return fields, rows


def family(order):
    """Scheme family of an order id: 'GreedyMini_w5_k15_canonical' -> 'GreedyMini'."""
    for fam in sorted(SCHEME_ORDER, key=len, reverse=True):
        if order == fam or order.startswith(fam + "_"):
            return fam
    return order.split("_")[0]


def order_key(order):
    return (SCHEME_ORDER.get(family(order), 99), order)


def library_key(lib):
    try:
        return (LIBRARY_ORDER.index(lib), lib)
    except ValueError:
        return (len(LIBRARY_ORDER), lib)


def dataset_key(ds):
    return (DATASET_ORDER.get(ds, 9), ds)


def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def split_context(context):
    """'standard_p10/GreedyMini_w5_k15_canonical/3 [fda-argos]' ->
    (library, order, seed_idx-as-str)."""
    base = (context or "").split(" [")[0]
    parts = [p for p in base.split("/") if p]
    lib = parts[0] if parts else ""
    order = parts[1] if len(parts) > 1 else ""
    return lib, order, (parts[2] if len(parts) > 2 else "")


# ---------------------------------------------------------------------------
# metric catalogue: one place defines columns, categories, meanings, ratio sense
# ---------------------------------------------------------------------------
# (name, category, meaning, ratio direction: lower/higher better, or None)
# Categories double as the row-1 band labels of the runs/summary sheets.

CAT_ID = "id / params"
CAT_BUILD = "build"
CAT_CLASSIFY = "classify"
CAT_DENSITY = "density"
CAT_NEXTMIN = "nextminimizer"
CAT_VSREF = "vs reference"
CAT_ACC = "accuracy"

METRICS = [
    ("build_s", CAT_BUILD, "k2 build wall time (s)", "lower"),
    ("build_gb", CAT_BUILD, "k2 build peak RSS (GB)", "lower"),
    ("class_s", CAT_CLASSIFY, "classify wall time (s). runs: mean over the timed repeats of this seed (the cold round is excluded); summary: mean over seeds", "lower"),
    ("class_gb", CAT_CLASSIFY, "classify peak RSS (GB), same averaging as class_s", "lower"),
    ("mrd_min", CAT_CLASSIFY, "classification speed in million reads per minute", "higher"),
    ("cls_pct", CAT_CLASSIFY, "percent of reads classified (any taxon)", None),
    ("fac_expected", CAT_DENSITY, "minimizer density factor = density x (w+1) on random i.i.d. DNA (Monte Carlo; 2 = random minimizer theory). Order-level: measured by benchmark_minimizers.sh per seed, joined to the run by seed_idx where available", "lower"),
    ("fac_particular", CAT_DENSITY, "density factor on this library's FASTAs (benchmark_minimizers.sh -d)", "lower"),
    ("ns_per_minimizer", CAT_DENSITY, "Rank() cost in ns per minimizer (benchmark_minimizers.sh microbenchmark)", "lower"),
    ("nextmin_s", CAT_NEXTMIN, "time in NextMinimizer / Rank (s), needs a classify built with ENABLE_MINIMIZER_PROF=1", "lower"),
    ("nextmin_ns", CAT_NEXTMIN, "NextMinimizer cost in ns per call", "lower"),
    ("nextmin_frac", CAT_NEXTMIN, "NextMinimizer time as % of classify wall", "lower"),
    ("lookups_pm", CAT_NEXTMIN, "hash-table lookups per Mbp of reads", "lower"),
    ("lookup_frac", CAT_NEXTMIN, "lookup time as % of classify time (per-thread sums)", "lower"),
    ("mape_g", CAT_VSREF, "paired MAPE of genus abundance vs Random of the same seed_idx (%, abundance floor --min-abund)", None),
    ("mape_s", CAT_VSREF, "same at species rank", None),
    ("l1_g", CAT_VSREF, "L1 distance of genus abundance vs Random (%)", None),
    ("concordance", CAT_VSREF, "ClassConc: reads with identical taxid vs Random, % of all reads", None),
    ("classified_concordance", CAT_VSREF, "identical taxid among reads classified by both (%)", None),
    ("jaccard", CAT_VSREF, "Jaccard of the sets of classified read ids vs Random (%)", None),
    ("newly_classified_pct", CAT_VSREF, "reads unclassified by Random but classified by this order (% of all reads)", None),
    ("newly_unclassified_pct", CAT_VSREF, "reads classified by Random but unclassified by this order (% of all reads)", None),
]
# runs-only extras (not aggregated into summary): within-seed std, counts, phases
RUN_ONLY = [
    ("class_s_std", CAT_CLASSIFY, "std of classify wall over the timed repeats of this seed (s)"),
    ("class_gb_std", CAT_CLASSIFY, "std of classify RSS over the timed repeats (GB)"),
    ("mrd_min_std", CAT_CLASSIFY, "std of speed over repeats (Mread/min)"),
    ("cls_pct_std", CAT_CLASSIFY, "std of classified % over repeats"),
    ("n_classify", CAT_CLASSIFY, "number of timed classify repeats (cold round excluded)"),
    ("n_reads", CAT_CLASSIFY, "reads (read pairs) in the query set"),
    ("nextmin_s_std", CAT_NEXTMIN, "std of nextmin_s over repeats"),
    ("nextmin_cpm", CAT_NEXTMIN, "NextMinimizer calls per Mbp of reads"),
    ("scan_s", CAT_NEXTMIN, "scan phase (minimizer + tokens), sum of per-thread seconds"),
    ("lookup_s", CAT_NEXTMIN, "lookup phase (GetBatch), sum of per-thread seconds"),
    ("resolve_s", CAT_NEXTMIN, "resolve phase (token replay + ResolveTree), sum of per-thread seconds"),
    ("repeats_pm", CAT_NEXTMIN, "repeated-key hits per Mbp of reads"),
    ("uhs_frac", CAT_NEXTMIN, "% of Rank() calls on the finite-rank (UHS) path; GreedyMini orders only"),
    ("non_uhs_frac", CAT_NEXTMIN, "% of Rank() calls on the infinity early-out path"),
]
# headline accuracy columns pivoted into runs/summary: (name, meaning, ratio dir,
# dataset, selector of the detail row, source field)
ACC_COLS = []
for _rank in ("genus", "species"):
    for _m, _lbl in (("sensitivity", "sens"), ("ppv", "ppv"), ("f1", "f1")):
        ACC_COLS.append(("se_%s_%s" % (_lbl, _rank), "strain-exclusion",
                         "strain-exclusion %s at %s rank, domain all (Mason held-out genomes)" % (_m, _rank),
                         "higher", ("labeled", {"domain": "all", "rank": _rank}, _m)))
for _k, _lbl in (("concordant", "concordant"), ("discordant", "discordant"),
                 ("ancestor", "ancestor"), ("unclassified", "unclassified")):
    ACC_COLS.append(("argos_%s" % _lbl, "fda-argos",
                     "FDA-ARGOS genus (%s) fraction of fragments, unweighted mean over the experiments" % _lbl,
                     "higher" if _lbl == "concordant" else ("lower" if _lbl == "discordant" else None),
                     ("argos", {"rank": "genus", "scope": "unweighted_mean"}, _k + "_frac")))
for _rank in ("genus", "species"):
    for _m, _lbl in (("sensitivity", "sens"), ("ppv", "ppv"), ("f1", "f1"),
                     ("excluded_frac", "excluded")):
        ACC_COLS.append(("cami_%s_%s" % (_lbl, _rank), "cami",
                         "CAMI %s at %s rank, group all, reads whose true taxon is in the library%s" % (
                             _m, _rank, " (reads with no reference: reported, not scored)" if _m == "excluded_frac" else ""),
                         "higher" if _m != "excluded_frac" else None,
                         ("cami", {"group": "all", "rank": _rank}, _m)))
for _rank in ("genus", "species"):
    for _m, _lbl, _dir in (("f1", "f1", "higher"), ("bray_curtis", "bray_curtis", "lower"),
                           ("off_target_frac", "off_target", "lower")):
        ACC_COLS.append(("mock_%s_%s" % (_lbl, _rank), "mock",
                         "mock community (ZymoBIOMICS D6300) %s at %s rank, scope all, detection threshold min_frac=1e-4" % (_m, _rank),
                         _dir, ("mock", {"scope": "all", "rank": _rank, "min_frac": MOCK_HEADLINE_MIN_FRAC}, _m)))

RUN_ID_FIELDS = ["run_id", "library", "dataset", "order", "seed_idx", "seed"]
SUMMARY_ID_FIELDS = ["library", "dataset", "order", "n_seeds"]
METRIC_DIR = dict((m[0], m[3]) for m in METRICS)
for _c in ACC_COLS:
    METRIC_DIR[_c[0]] = _c[3]
METRIC_CAT = dict((m[0], m[1]) for m in METRICS)
for _m in RUN_ONLY:
    METRIC_CAT[_m[0]] = _m[1]
for _c in ACC_COLS:
    METRIC_CAT[_c[0]] = CAT_ACC
# which accuracy columns exist for which dataset (prune the others per dataset)
SUMMARY_METRICS = [m[0] for m in METRICS] + [c[0] for c in ACC_COLS]
# the paired-vs-reference metrics are already relative to Random: no ratio
NO_RATIO = frozenset(m[0] for m in METRICS if m[1] == CAT_VSREF)


def runs_field_order():
    return (RUN_ID_FIELDS + [m[0] for m in METRICS if m[1] in (CAT_BUILD,)]
            + [m[0] for m in METRICS if m[1] == CAT_CLASSIFY]
            + [m[0] for m in RUN_ONLY if m[1] == CAT_CLASSIFY]
            + [m[0] for m in METRICS if m[1] == CAT_DENSITY]
            + [m[0] for m in METRICS if m[1] == CAT_NEXTMIN]
            + [m[0] for m in RUN_ONLY if m[1] == CAT_NEXTMIN]
            + [m[0] for m in METRICS if m[1] == CAT_VSREF]
            + [c[0] for c in ACC_COLS])


def _dedupe(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def summary_field_order():
    fields = list(SUMMARY_ID_FIELDS)
    by_cat = OrderedDict()
    for m in SUMMARY_METRICS:
        by_cat.setdefault(METRIC_CAT[m], []).append(m)
    order = (CAT_BUILD, CAT_CLASSIFY, CAT_DENSITY, CAT_NEXTMIN, CAT_VSREF, CAT_ACC)
    for cat in order:
        for m in by_cat.get(cat, []):
            fields += [m + "_mean", m + "_std"] + ([] if m in NO_RATIO else [m + "_ratio"])
    return fields


def field_category(field):
    if field in RUN_ID_FIELDS or field in SUMMARY_ID_FIELDS:
        return CAT_ID
    base = field
    for suf in ("_mean", "_std", "_ratio"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return METRIC_CAT.get(base) or METRIC_CAT.get(field) or ""


# ---------------------------------------------------------------------------
# runs (+ classify_repeats) from the collected per-seed metrics
# ---------------------------------------------------------------------------

def _detail_pick(rows, selector, ds_key):
    """First detail row of (seed) matching `selector` -> value of the source field."""
    kind, match, src = selector
    for r in rows.get(kind, ()):
        ok = True
        for k, v in match.items():
            have = r.get(k)
            if isinstance(v, float):
                try:
                    if abs(float(have) - v) > 1e-12:
                        ok = False
                        break
                except (TypeError, ValueError):
                    ok = False
                    break
            elif (have if have is not None else "") != v:
                ok = False
                break
        if ok:
            val = r.get(src)
            return val if val not in ("", None) else None
    return None


def build_runs(library, per_seed_rows, acc_rows, dens_rows):
    """per_seed_rows: dicts from compare_benchmark (scheme, dataset, seed_idx,
    context, seed + collect_run_metrics keys incl. class_walls).
    acc_rows: {'labeled'|'argos'|'cami'|'mock': [rows with context/dataset/seed_idx]}.
    dens_rows: density_per_seed rows ({scheme, seed_idx, fac_*...})."""
    dens = {}
    for d in dens_rows:
        dens[(d.get("scheme"), _int(d.get("seed_idx")))] = d
    # index accuracy rows by (order, dataset, seed_idx)
    idx = defaultdict(lambda: defaultdict(list))
    for kind, rows in acc_rows.items():
        for r in rows:
            _lib, order, _ = split_context(r.get("context"))
            key = (order, r.get("dataset") or _kind_dataset(kind), _int(r.get("seed_idx")))
            idx[key][kind].append(r)
    fields = runs_field_order()
    runs, repeats = [], []
    for row in per_seed_rows:
        _lib, order, _ = split_context(row.get("context"))
        ds = row.get("dataset") or "speed"
        seed_idx = _int(row.get("seed_idx"))
        out = OrderedDict((f, None) for f in fields)
        out["run_id"] = "%s/%s/%s" % (library, order, seed_idx)
        out["library"] = library
        out["dataset"] = ds
        out["order"] = order
        out["seed_idx"] = seed_idx
        out["seed"] = _seed_hex(row.get("seed"))
        for f in fields:
            if f in RUN_ID_FIELDS:
                continue
            if f in row:
                out[f] = clean_value(row[f])
        # per-seed density row, else the scheme-level one (summary text: no seed_idx)
        d = (dens.get((family(order), seed_idx)) or dens.get((order, seed_idx))
             or dens.get((family(order), None)))
        if d:
            for f in ("fac_expected", "fac_particular", "ns_per_minimizer"):
                out[f] = clean_value(d.get(f))
        picked = idx.get((order, ds, seed_idx), {})
        for name, cds, _mean, _dir, selector in ACC_COLS:
            if cds != ds:
                continue
            out[name] = clean_value(_detail_pick(picked, selector, ds))
        runs.append(out)
        for i, w in enumerate(row.get("class_walls") or [], 1):
            repeats.append(OrderedDict([
                ("library", library), ("dataset", ds), ("order", order),
                ("seed_idx", seed_idx), ("repeat", i), ("class_s", w)]))
    runs.sort(key=_run_sort_key)
    repeats.sort(key=lambda r: (library_key(r["library"]), dataset_key(r["dataset"]),
                                order_key(r["order"]), r["seed_idx"] or 0, r["repeat"]))
    return runs, repeats


CLASSIFY_REPEAT_FIELDS = ["library", "dataset", "order", "seed_idx", "repeat", "class_s"]


def _kind_dataset(kind):
    return {"labeled": "strain-exclusion", "argos": "fda-argos", "cami": "cami", "mock": "mock"}[kind]


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _seed_hex(seed):
    s = str(seed or "").strip()
    if s and set(s) <= set("01") and len(s) > 8:
        return "0x%x" % int(s, 2)
    return s or None


def _run_sort_key(r):
    return (library_key(r.get("library") or ""), dataset_key(r.get("dataset") or ""),
            order_key(r.get("order") or ""), r.get("seed_idx") or 0)


# ---------------------------------------------------------------------------
# derived sheets: summary, vs_reference, timing_correlation, detail summaries
# ---------------------------------------------------------------------------

def summary_rows(runs):
    groups = OrderedDict()
    for r in runs:
        groups.setdefault((r.get("library"), r.get("dataset"), r.get("order")), []).append(r)
    out = []
    for (lib, ds, order), parts in groups.items():
        row = OrderedDict([("library", lib), ("dataset", ds), ("order", order),
                           ("n_seeds", len(set(p.get("seed_idx") for p in parts)))])
        for m in SUMMARY_METRICS:
            vals = [p.get(m) for p in parts if _num(p.get(m))]
            mu, sd = mean_std(vals)
            row[m + "_mean"] = round_sig(mu)
            row[m + "_std"] = round_std(sd)
        out.append(row)
    out.sort(key=lambda r: (library_key(r["library"]), dataset_key(r["dataset"]), order_key(r["order"])))
    attach_ratios(out)
    return out


def attach_ratios(rows):
    """{metric}_ratio = this order's mean / the reference (Random) mean of the
    same library and dataset (reference = 1)."""
    base = {}
    for r in rows:
        if family(r["order"]) == REF_SCHEME:
            base[(r["library"], r["dataset"])] = r
    for r in rows:
        b = base.get((r["library"], r["dataset"]))
        for m in SUMMARY_METRICS:
            if m in NO_RATIO:
                continue
            mu = r.get(m + "_mean")
            bm = b.get(m + "_mean") if b else None
            r[m + "_ratio"] = round(mu / bm, 4) if (_num(mu) and _num(bm) and bm) else None


def prune_empty(fields, rows, keep):
    """Drop columns without a single value (keep the id columns)."""
    out = []
    for f in fields:
        if f in keep or any(r.get(f) not in (None, "") for r in rows):
            out.append(f)
    return out


VS_REF_METRICS = [m[0] for m in METRICS if m[1] == CAT_VSREF]


def vs_reference_rows(runs):
    groups = OrderedDict()
    for r in runs:
        if family(r["order"]) == REF_SCHEME:
            continue
        groups.setdefault((r["library"], r["dataset"], r["order"]), []).append(r)
    out = []
    for (lib, ds, order), parts in groups.items():
        if not any(_num(p.get(m)) for p in parts for m in VS_REF_METRICS):
            continue
        row = OrderedDict([("library", lib), ("dataset", ds), ("order", order),
                           ("ref_order", REF_SCHEME), ("n_seeds", len(parts))])
        for m in VS_REF_METRICS:
            mu, sd = mean_std([p.get(m) for p in parts if _num(p.get(m))])
            row[m + "_mean"] = round_sig(mu)
            row[m + "_std"] = round_std(sd)
        out.append(row)
    out.sort(key=lambda r: (library_key(r["library"]), dataset_key(r["dataset"]), order_key(r["order"])))
    return out


def vs_reference_fields():
    return (["library", "dataset", "order", "ref_order", "n_seeds"]
            + [x for m in VS_REF_METRICS for x in (m + "_mean", m + "_std")])


def _pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / float(n), sum(ys) / float(n)
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (dx * dy)


CORREL_FIELDS = ["library", "dataset", "metric", "n_orders", "r_fac_expected", "r_fac_particular"]


def timing_correlation_rows(summary):
    """Pearson r (over orders) of the density factor vs build / classify time and
    RSS: does a lower density factor (fewer minimizers) mean a faster / smaller run?"""
    out = []
    groups = OrderedDict()
    for r in summary:
        groups.setdefault((r["library"], r["dataset"]), []).append(r)
    for (lib, ds), rows in groups.items():
        for metric in ("build_s", "build_gb", "class_s", "class_gb"):
            xs = {}
            for tag, col in (("e", "fac_expected_mean"), ("p", "fac_particular_mean")):
                pts = [(r[col], r[metric + "_mean"]) for r in rows
                       if _num(r.get(col)) and _num(r.get(metric + "_mean"))]
                xs[tag] = (_pearson([p[0] for p in pts], [p[1] for p in pts]), len(pts))
            n = max(xs["e"][1], xs["p"][1])
            if n < 3:
                continue
            out.append(OrderedDict([
                ("library", lib), ("dataset", ds), ("metric", metric), ("n_orders", n),
                ("r_fac_expected", round(xs["e"][0], 4) if xs["e"][0] is not None else None),
                ("r_fac_particular", round(xs["p"][0], 4) if xs["p"][0] is not None else None)]))
    return out


def _group_summary(rows, key_fields, metrics, order_field="order"):
    groups = OrderedDict()
    for r in rows:
        groups.setdefault(tuple(r.get(k) if r.get(k) is not None else "" for k in key_fields), []).append(r)
    fields = list(key_fields) + ["n_seeds"] + [x for m in metrics for x in (m + "_mean", m + "_std")]
    out = []
    for key, parts in groups.items():
        rec = OrderedDict(zip(key_fields, key))
        rec["n_seeds"] = len(set(p.get("seed_idx") for p in parts))
        for m in metrics:
            mu, sd = mean_std([p.get(m) for p in parts if _num(p.get(m))])
            rec[m + "_mean"] = round_sig(mu)
            rec[m + "_std"] = round_std(sd)
        out.append(rec)

    def sk(r):
        return (library_key(r.get("library") or ""), dataset_key(r.get("dataset") or ""),
                order_key(r.get("order") or ""), DOMAIN_ORDER.get(r.get("domain"), 9),
                RANK_ORDER.get(r.get("rank"), 9), str(r.get("scope") or ""),
                str(r.get("group") or ""), str(r.get("min_frac") if r.get("min_frac") is not None else ""))
    out.sort(key=sk)
    return fields, out


def strain_exclusion_summary(rows):
    return _group_summary(rows, ("library", "dataset", "order", "domain", "rank"),
                          ("sensitivity", "ppv", "f1"))


def _argos_headline(r):
    if r.get("rank") not in (None, "", "genus"):
        return False
    scope = r.get("scope") or ""
    if scope:
        return scope == "unweighted_mean"
    return not (r.get("sra_run") or "")


def fda_argos_summary(rows):
    head = [r for r in rows if _argos_headline(r)]
    fields, out = _group_summary(
        head, ("library", "dataset", "order"),
        ("concordant_frac", "discordant_frac", "ancestor_frac", "unclassified_frac"))
    n_exp = defaultdict(set)
    for r in rows:
        if (r.get("scope") or "") == "experiment" and r.get("sra_run"):
            n_exp[(r.get("library"), r.get("order"))].add(r["sra_run"])
    for rec in out:
        rec["n_experiments"] = len(n_exp.get((rec["library"], rec["order"]), ()))
    fields = ["library", "dataset", "order", "n_seeds", "n_experiments"] + [
        x for m in ("concordant", "discordant", "ancestor", "unclassified")
        for x in (m + "_mean", m + "_std")]
    for rec in out:
        for m in ("concordant", "discordant", "ancestor", "unclassified"):
            rec[m + "_mean"] = rec.pop(m + "_frac_mean")
            rec[m + "_std"] = rec.pop(m + "_frac_std")
    return fields, out


CAMI_METRICS = ("sensitivity", "ppv", "f1", "excluded_frac", "n_scored")
MOCK_METRICS = ("precision", "recall", "f1", "l1", "bray_curtis", "classified_frac",
                "on_target_frac", "off_target_frac")


def cami_summary(rows):
    return _group_summary(rows, ("library", "dataset", "order", "group", "rank"), CAMI_METRICS)


def mock_summary(rows):
    return _group_summary(rows, ("library", "dataset", "order", "scope", "rank", "min_frac"), MOCK_METRICS)


DETAIL_KEYS = {
    "strain_exclusion_runs": "labeled",
    "fda_argos_runs": "argos",
    "cami_runs": "cami",
    "mock_runs": "mock",
}


def detail_rows(library, kind_rows):
    """Prefix accuracy rows (from compare_benchmark.load_accuracy_sheets) with
    library / order columns. Returns (fields, rows)."""
    if not kind_rows:
        return [], []
    out = []
    for r in kind_rows:
        _lib, order, _ = split_context(r.get("context"))
        row = OrderedDict()
        row["library"] = library
        row["dataset"] = r.get("dataset") or None
        row["order"] = order
        for k, v in r.items():
            if k in ("scheme", "context", "dataset"):
                continue
            row[k] = clean_value(v)
        out.append(row)
    out.sort(key=lambda r: (dataset_key(r.get("dataset") or ""), order_key(r["order"]),
                            _int(r.get("seed_idx")) or 0, DOMAIN_ORDER.get(r.get("domain"), 9),
                            RANK_ORDER.get(r.get("rank"), 9), str(r.get("scope") or ""),
                            str(r.get("sra_run") or ""), str(r.get("group") or ""),
                            str(r.get("min_frac") if r.get("min_frac") is not None else "")))
    fields = []
    for r in out:
        for k in r:
            if k not in fields:
                fields.append(k)
    lead = ["library", "dataset", "order", "seed_idx", "seed"]
    fields = [f for f in lead if f in fields] + [f for f in fields if f not in lead]
    return fields, out


# ---------------------------------------------------------------------------
# orders
# ---------------------------------------------------------------------------

ORDER_INFO = {
    "Random": ("Minimizer", "RandomXor", "--sketching-scheme Minimizer --order RandomXor",
               "random order: the l-mer XORed with the seed (toggle mask), ranked by the hash of that value"),
    "RandomPcg": ("Minimizer", "RandomPcg", "--sketching-scheme Minimizer --order RandomPcg",
                  "random order with a PCG scramble of the seed-XORed l-mer"),
    "GreedyMini": ("Minimizer", "GreedyMini", "--sketching-scheme Minimizer --order GreedyMini --greedy-mini-file <table>.gm --gm-alphabet binary",
                   "GreedyMini table (binary alphabet, centre extension), ties broken by the seed"),
    "GreedyMiniDNA": ("Minimizer", "GreedyMini", "--sketching-scheme Minimizer --order GreedyMini --greedy-mini-file <table>.gm --gm-alphabet dna",
                      "GreedyMini table over the DNA alphabet"),
    "ModRandom": ("Mod", "RandomXor", "--sketching-scheme Mod --order RandomXor",
                  "mod-minimizer (Groot Koerkamp & Pibiri) on a random order"),
    "ModGreedyMini": ("Mod", "GreedyMini", "--sketching-scheme Mod --order GreedyMini --greedy-mini-file <table>.gm --gm-alphabet binary",
                      "mod-minimizer on a GreedyMini table"),
}


def order_params_from_db(run_root, order):
    """(k, l, spaced_seed_mask) of an order from its db/opts.k2d (any seed), or None."""
    odir = os.path.join(run_root, order)
    if not os.path.isdir(odir):
        return None
    for seed in sorted(os.listdir(odir)):
        p = os.path.join(odir, seed, "db", "opts.k2d")
        if os.path.isfile(p):
            try:
                with open(p, "rb") as f:
                    k, l, mask = struct.unpack("<QQQ", f.read(24))
                return k, l, mask
            except (OSError, struct.error):
                return None
    return None


def build_orders(library, run_root, runs, dens_rows):
    """One row per order seen in `runs`: CLI, l/k/w, spaced seed, density."""
    order_ids = _dedupe(r["order"] for r in runs)
    rows = []
    for order in sorted(order_ids, key=order_key):
        fam = family(order)
        sketching, rank_order, cli, desc = ORDER_INFO.get(fam, ("", "", "", ""))
        table = order[len(fam) + 1:] if order.startswith(fam + "_") else None
        if table:
            cli = cli.replace("<table>", "<orders>/GreedyMini/" + table)
        d_rows = [d for d in dens_rows if d.get("scheme") in (order, fam)]
        row = OrderedDict([
            ("order", order), ("scheme", fam), ("sketching", sketching),
            ("rank_order", rank_order), ("table", table), ("cli", cli),
            ("description", desc)])
        l = k = None
        for d in d_rows:
            if _num(d.get("l")) and _num(d.get("k")):
                l, k = int(d["l"]), int(d["k"])
                break
        mask = None
        opts = order_params_from_db(run_root, order) if run_root else None
        if opts:
            k, l, mask = opts[0], opts[1], opts[2]
        row["l"] = l
        row["k"] = k
        row["w"] = (k - l + 1) if (k is not None and l is not None) else None
        spaces = None
        if mask is not None and l:
            kept = bin(mask).count("1") // 2
            spaces = l - kept if mask else 0
        row["minimizer_spaces"] = spaces
        row["spaced_seed_mask"] = ("0x%x" % mask) if mask else ("none" if mask == 0 else None)
        e_mu, e_sd = mean_std([d.get("fac_expected") for d in d_rows if _num(d.get("fac_expected"))])
        p_mu, p_sd = mean_std([d.get("fac_particular") for d in d_rows if _num(d.get("fac_particular"))])
        row["expected_density_factor_mean"] = round_sig(e_mu)
        row["expected_density_factor_std"] = round_std(e_sd)
        row[library + "_particular_density_factor_mean"] = round_sig(p_mu)
        row[library + "_particular_density_factor_std"] = round_std(p_sd)
        rows.append(row)
    return rows


ORDERS_LEAD = ["order", "scheme", "sketching", "rank_order", "table", "cli", "description",
               "l", "k", "w", "minimizer_spaces", "spaced_seed_mask",
               "expected_density_factor_mean", "expected_density_factor_std"]


def orders_fields(rows):
    libs = []
    for r in rows:
        for k in r:
            m = re.match(r"^(.*)_particular_density_factor_mean$", k)
            if m and m.group(1) not in libs:
                libs.append(m.group(1))
    libs.sort(key=library_key)
    f = list(ORDERS_LEAD)
    for lib in libs:
        f += [lib + "_particular_density_factor_mean", lib + "_particular_density_factor_std"]
    return f


def merge_orders(lists):
    by = OrderedDict()
    for rows in lists:
        for r in rows:
            cur = by.setdefault(r["order"], OrderedDict())
            for k, v in r.items():
                if v is not None and cur.get(k) is None:
                    cur[k] = v
                elif k not in cur:
                    cur[k] = v
    return sorted(by.values(), key=lambda r: order_key(r["order"]))


def orders_category(field):
    if field in ("order", "scheme", "sketching", "rank_order", "table", "cli", "description"):
        return "order"
    if field in ("l", "k", "w", "minimizer_spaces", "spaced_seed_mask"):
        return "minimizer params"
    if field.startswith("expected_density_factor"):
        return "expected (theoretical)"
    marker = "_particular_density_factor_"
    if marker in field:
        return field.split(marker, 1)[0]
    return ""


# ---------------------------------------------------------------------------
# datasets (query sets) and libraries
# ---------------------------------------------------------------------------

DATASET_INFO = {
    "strain-exclusion": (
        "simulated reads (Mason) from the 50 Table S3 held-out genomes",
        "per-read taxid: the source genome (truth.tsv); genus/species sensitivity / PPV / F1",
        "Wood et al. 2019 Table S1; Mason 2 reads from Zenodo (download_data.sh / prepare_strex_data.sh)"),
    "fda-argos": (
        "real Illumina SRA runs (FDA-ARGOS, Table S5)",
        "per-experiment species (Table S5); genus concordance a-d, no sensitivity/PPV",
        "Wood et al. 2019 Table S5; ENA download of the sra_run accessions (prepare_fda_argos_reads.py)"),
    "cami": (
        "simulated reads (ART) of CAMI II marine genomes",
        "per-read true taxid (truth.tsv); library-aware sens / PPV / F1, excluded reads reported",
        "CAMI II marine sample, subsampled (prepare_cami_reads.py, SOURCE.txt)"),
    "mock": (
        "real reads of the ZymoBIOMICS D6300 mock community",
        "expected community composition only (composition.csv); profile-level precision/recall/L1/Bray-Curtis",
        "SRR12324251, subsampled (prepare_mock_reads.py, SOURCE.txt)"),
    "speed": (
        "unlabeled NCBI-library speed corpus",
        "none (timing only)",
        "data/speed/ or legacy data/benchmark_R*.fastq"),
}

DATASET_FIELDS = ["dataset", "sources", "n_sources", "fragments_per_source", "n_fragments",
                  "read_length_bp", "paired", "total_bases_r1_r2", "truth", "provenance",
                  "libraries", "path"]


def decorate_datasets(query_rows, runs):
    """query_rows: compare_benchmark.query_sets_sheet rows. Adds truth type,
    provenance and the libraries each set was run on; adds a row for a set that
    ran but has no data dir (e.g. `speed`)."""
    libs_by_ds = defaultdict(list)
    for r in runs:
        if r["library"] not in libs_by_ds[r["dataset"]]:
            libs_by_ds[r["dataset"]].append(r["library"])
    out = OrderedDict()
    for q in query_rows or []:
        out[q["dataset"]] = OrderedDict((f, q.get(f)) for f in DATASET_FIELDS)
    for ds in libs_by_ds:
        out.setdefault(ds, OrderedDict((f, None) for f in DATASET_FIELDS))
        out[ds]["dataset"] = ds
    for ds, row in out.items():
        info = DATASET_INFO.get(ds)
        if info:
            if not row.get("sources"):
                row["sources"] = info[0]
            row["truth"] = info[1]
            row["provenance"] = info[2]
        libs = sorted(libs_by_ds.get(ds, []), key=library_key)
        row["libraries"] = ", ".join(libs) if libs else None
    return sorted(out.values(), key=lambda r: dataset_key(r["dataset"]))


def merge_datasets(lists):
    by = OrderedDict()
    for rows in lists:
        for r in rows:
            cur = by.setdefault(r["dataset"], OrderedDict())
            for k, v in r.items():
                if k == "libraries":
                    have = [x.strip() for x in (cur.get(k) or "").split(",") if x.strip()]
                    for x in (v or "").split(","):
                        if x.strip() and x.strip() not in have:
                            have.append(x.strip())
                    cur[k] = ", ".join(sorted(have, key=library_key)) or None
                elif cur.get(k) is None:
                    cur[k] = v
    return sorted(by.values(), key=lambda r: dataset_key(r["dataset"]))


LIBRARY_INFO = {
    "standard": ("100% of the NCBI standard library (archaea, bacteria, human, viral, UniVec_Core)",
                 "-N (default)", "timing at full scale (speed corpus)"),
    "standard_p10": ("~10% random subsample of the standard library (SAMPLE_FRAC=0.1)",
                     "-m", "timing / density / speed experiments on a tractable library"),
    "standard_p0.5": ("~0.5% subsample of the standard library (SAMPLE_SMOKE_FRAC=0.005)",
                      "-z", "smoke tests only"),
    "strain_exclusion": ("Wood 2019 library with the 50 Table S3 taxa held out (strain_excluded.fna)",
                         "-E", "paper-style accuracy (strain-exclusion, FDA-ARGOS, CAMI, mock) and timing"),
}
LIBRARY_FIELDS = ["library", "role", "flag", "fraction", "taxonomy", "n_sequences", "n_taxa",
                  "n_fasta_files", "fasta_gb", "n_orders", "n_runs", "query_sets_run",
                  "orders_run", "status", "db_path", "run_root"]


def _scan_prelim(paths):
    seqs, taxa = 0, set()
    for p in paths:
        try:
            with open(p, errors="replace") as f:
                for line in f:
                    if not line.startswith("TAXID\t"):
                        continue
                    parts = line.split("\t")
                    seqs += 1
                    if len(parts) > 2:
                        taxa.add(parts[2])
        except OSError:
            continue
    return seqs, len(taxa)


def _library_files(db):
    """(fasta files, prelim_map files) of a Kraken library dir (recursive,
    one level of subdirectories as in library/<kind>/ and library/added/)."""
    fasta, prelim = [], []
    lib = os.path.join(db, "library")
    for root, _dirs, files in os.walk(lib):
        for f in files:
            full = os.path.join(root, f)
            if re.search(r"\.(fna|fa|fasta)$", f) and not os.path.islink(full):
                if os.path.getsize(full) > 0:
                    fasta.append(full)
            elif f.startswith("prelim_map") and f.endswith(".txt"):
                prelim.append(full)
    return fasta, prelim


def library_stats(db, cache):
    """Cheap library stats from prelim_map*.txt (one row per sequence) and the
    FASTA file sizes (no FASTA is read); cached by (path, size, mtime)."""
    key = os.path.abspath(db)
    fasta, prelim = _library_files(db)
    sig = [[p, os.path.getsize(p), int(os.path.getmtime(p))] for p in sorted(fasta + prelim)]
    hit = cache.get(key)
    if hit and hit.get("sig") == sig:
        return hit["stats"]
    seqs, ntax = _scan_prelim(prelim)
    stats = {
        "n_sequences": seqs or None,
        "n_taxa": ntax or None,
        "n_fasta_files": len(fasta),
        "fasta_gb": round(sum(os.path.getsize(p) for p in fasta) / 1e9, 3) if fasta else None,
    }
    cache[key] = {"sig": sig, "stats": stats}
    return stats


def _load_cache(path):
    import json
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _save_cache(path, cache):
    import json
    try:
        d = os.path.dirname(path)
        if not os.path.isdir(d):
            os.makedirs(d)
        with open(path, "w") as f:
            json.dump(cache, f)
    except OSError:
        pass


def build_libraries(base, runs, cache_path=None):
    """One row per library: every id in LIBRARY_ORDER with a db/<id> or runs/<id>
    dir, plus any other library present in `runs`."""
    cache = _load_cache(cache_path) if cache_path else {}
    ids = [i for i in LIBRARY_ORDER
           if os.path.isdir(os.path.join(base, "db", i)) or os.path.isdir(os.path.join(base, "runs", i))]
    for r in runs:
        if r["library"] not in ids:
            ids.append(r["library"])
    rows = []
    for lib in ids:
        db = os.path.join(base, "db", lib)
        info = LIBRARY_INFO.get(lib, ("", "", ""))
        lr = [r for r in runs if r["library"] == lib]
        st = {"n_sequences": None, "n_taxa": None, "n_fasta_files": None, "fasta_gb": None}
        if os.path.isdir(os.path.join(db, "library")):
            st = library_stats(db, cache)
        tax = None
        tdir = os.path.join(db, "taxonomy")
        if os.path.isdir(tdir):
            tax = "NCBI taxonomy (nodes.dmp / names.dmp)"
            if os.path.islink(tdir):
                tax += ", shared with " + os.path.basename(os.path.dirname(os.path.realpath(tdir)))
        orders = sorted(set(r["order"] for r in lr), key=order_key)
        sets = sorted(set(r["dataset"] for r in lr), key=dataset_key)
        rows.append(OrderedDict([
            ("library", lib), ("role", info[2] or None), ("flag", info[1] or None),
            ("fraction", info[0] or None), ("taxonomy", tax),
            ("n_sequences", st["n_sequences"]), ("n_taxa", st["n_taxa"]),
            ("n_fasta_files", st["n_fasta_files"]), ("fasta_gb", st["fasta_gb"]),
            ("n_orders", len(orders) or None), ("n_runs", len(set(r["run_id"] for r in lr)) or None),
            ("query_sets_run", ", ".join(sets) or None), ("orders_run", ", ".join(orders) or None),
            ("status", "runs + db" if lr and os.path.isdir(db) else "runs" if lr else
             "db only (no runs yet)" if os.path.isdir(db) else "runs dir only"),
            ("db_path", db if os.path.isdir(db) else None),
            ("run_root", os.path.join(base, "runs", lib) if os.path.isdir(os.path.join(base, "runs", lib)) else None)]))
    if cache_path:
        _save_cache(cache_path, cache)
    return rows


# ---------------------------------------------------------------------------
# overview (column reference)
# ---------------------------------------------------------------------------

SHEET_DOC = OrderedDict([
    ("overview", "this sheet: what every sheet and column means, grouped by category"),
    ("libraries", "one row per Kraken library: size, taxonomy, role in the benchmark, which query sets / orders ran on it"),
    ("datasets", "one row per query set: reads, read length, truth type, provenance"),
    ("orders", "one row per minimizer scheme/order: CLI flags, l/k/w, spaced seed, expected and per-library particular density factor"),
    ("summary", "library x dataset x order: mean / std over seeds, and `_ratio` = mean / the Random mean of the same library and dataset"),
    ("runs", "one row per seed x library x dataset x order: the raw per-run numbers"),
    ("vs_reference", "paired agreement with Random of the same seed (MAPE, concordance, Jaccard, new classified/unclassified), mean / std over seeds"),
    ("timing_correlation", "Pearson r over orders of density factor vs build / classify time and RSS (does fewer minimizers mean faster?)"),
    ("classify_repeats", "one row per timed classify repeat (wall seconds): the distribution behind class_s"),
    ("strain_exclusion_summary", "strain-exclusion sensitivity / PPV / F1 mean / std over seeds per domain x rank"),
    ("strain_exclusion_runs", "strain-exclusion per seed: TP / VP / FP / FN and sens / PPV / F1 per domain x rank"),
    ("fda_argos_summary", "FDA-ARGOS genus concordance (concordant / discordant / ancestor / unclassified) mean / std over seeds"),
    ("fda_argos_runs", "FDA-ARGOS per seed: per SRR experiment, the unweighted mean and the pooled fragment row"),
    ("cami_summary", "CAMI II marine library-aware sens / PPV / F1 and excluded fraction, mean / std over seeds per group x rank"),
    ("cami_runs", "CAMI II per seed: TP / VP / FP / FN, excluded reads and what was called for them"),
    ("mock_summary", "ZymoBIOMICS D6300 mock precision / recall / F1, L1, Bray-Curtis and off-target fraction mean / std per scope x rank x threshold"),
    ("mock_runs", "mock community profile metrics per seed"),
])

COMMON_DOC = [
    ("all", CAT_ID, "library", "Kraken library id: standard | standard_p10 | standard_p0.5 | strain_exclusion (see `libraries`)"),
    ("all", CAT_ID, "dataset", "query set: strain-exclusion | fda-argos | cami | mock | speed (see `datasets`)"),
    ("all", CAT_ID, "order", "minimizer scheme/order id = the run directory name: Random, RandomPcg, GreedyMini_<table>, ModRandom, ModGreedyMini_<table> (see `orders`)"),
    ("runs", CAT_ID, "run_id", "<library>/<order>/<seed_idx>"),
    ("runs", CAT_ID, "seed_idx / seed", "seed index and its 64-bit toggle mask in hex (`k2 --toggle-mask`): XOR tie-break / random-order seed of every scheme"),
    ("summary", CAT_ID, "n_seeds", "number of seeds aggregated into the row"),
    ("summary", "all metrics", "{metric}_mean / {metric}_std", "mean / sample std over seeds of the per-seed value in `runs` (std rounded to >= 4 significant digits)"),
    ("summary", "all metrics", "{metric}_ratio", "this order's mean / the Random mean of the same library and dataset (Random = 1); colour scale: green = better, red = worse. Not given for the paired vs-reference metrics (already relative to Random)"),
    ("summary", CAT_DENSITY, "fac_* / ns_per_minimizer (summary)", "mean / std over the seeds that have a density measurement (tables/density_per_seed.csv; an order measured only once, e.g. from minimizer_density_summary.txt, repeats that value, std 0)"),
    ("runs", CAT_ACC, "se_* / argos_* / cami_* / mock_*", "headline accuracy of the dataset the row belongs to (blank for the other datasets); the full domain x rank / per-SRR / per-group detail is in the *_runs / *_summary sheets"),
    ("libraries", "libraries", "n_sequences / n_taxa / fasta_gb", "sequences and distinct taxids in the library's prelim_map*.txt, FASTA size on disk (GB, library.fna files; not re-read)"),
    ("orders", "order", "cli", "k2 build flags that select the order (<orders> = the repo's orders/ dir)"),
    ("orders", "minimizer params", "l / k / w / minimizer_spaces / spaced_seed_mask", "minimizer length l, k-mer length k, window w = k - l + 1, spaced-seed spaces s (pattern 1*(l-2s) + (01)*s, so l - s kept positions) and its mask, read from the order's db/opts.k2d"),
    ("orders", "expected (theoretical)", "expected_density_factor_mean / std", "density factor x (w+1) on random i.i.d. DNA, mean / std over the density seeds (independent of the library)"),
    ("orders", "<library>", "<library>_particular_density_factor_mean / std", "the same factor measured on that library's FASTAs"),
    ("orders", "order", "scheme / sketching / rank_order / table / description", "scheme family (Random, RandomPcg, GreedyMini, ModRandom, ModGreedyMini), k2 --sketching-scheme (Minimizer | Mod), the --order it ranks by, the GreedyMini table (w<w>_k<gm_k>_<variant>.gm, blank for random orders) and a one-line description"),
    ("libraries", "libraries", "role / flag / fraction / taxonomy", "what the library is for in the benchmark, the run_benchmark.sh flag that selects it, how much of NCBI it holds, and the taxonomy it uses"),
    ("libraries", "libraries", "n_fasta_files / n_orders / n_runs", "library.fna files found, orders and distinct runs (order x seed) present under runs/<library>"),
    ("libraries", "libraries", "query_sets_run / orders_run / status", "datasets and orders classified on this library; whether it has runs and/or a db"),
    ("libraries", "libraries", "db_path / run_root", "shared library FASTA dir (db/<id>) and the run tree (runs/<id>)"),
    ("datasets", "dataset", "sources / n_sources / fragments_per_source", "what the reads come from, how many source genomes / taxa / experiments, fragments per source"),
    ("datasets", "dataset", "n_fragments / read_length_bp / paired / total_bases_r1_r2 / path", "read pairs, read length (bp, or min-max), paired-end flag, bases in R1+R2, location of the FASTQs"),
    ("datasets", "dataset", "truth / provenance / libraries", "truth type used for scoring, where the reads come from, libraries the set was classified against"),
    ("vs_reference", CAT_ID, "ref_order / n_seeds", "the reference order of the comparison (Random) and the number of paired seeds"),
    ("timing_correlation", CAT_ID, "metric / n_orders", "time or RSS metric correlated with the density factor, and the number of orders in the correlation"),
    ("classify_repeats", CAT_ID, "seed_idx", "seed the repeat belongs to"),
]
SHEET_FIELD_DOC = [
    ("vs_reference", CAT_VSREF, "{metric}_mean / {metric}_std", "mean / std over seeds of the paired per-seed metric; one row per non-reference order"),
    ("timing_correlation", "correlation", "r_fac_expected / r_fac_particular", "Pearson r over the orders of that library and dataset (n_orders >= 3)"),
    ("classify_repeats", "classify", "repeat / class_s", "timed repeat number (the cold round is excluded) and its classify wall time (s)"),
]


def overview_rows():
    rows = []

    def add(sheet, cat, column, meaning):
        rows.append(OrderedDict([("sheet", sheet), ("category", cat), ("column", column), ("meaning", meaning)]))

    add("overview", "sheets", "(sheet)", SHEET_DOC["overview"])
    for name, doc in SHEET_DOC.items():
        if name != "overview":
            add("sheets", "sheets", name, doc)
    for sheet, cat, col, meaning in COMMON_DOC:
        add(sheet, cat, col, meaning)
    for name, cat, meaning, dir_ in METRICS:
        add("summary/runs", cat, name, meaning)
    for name, cat, meaning in RUN_ONLY:
        add("runs", cat, name, meaning)
    for name, ds, meaning, dir_, _sel in ACC_COLS:
        add("summary/runs", CAT_ACC, name, meaning)
    for sheet, cat, col, meaning in SHEET_FIELD_DOC:
        add(sheet, cat, col, meaning)
    add("strain_exclusion_*", "strain-exclusion", "domain / rank / n_truth / n_calls / tp / vp / fp / fn / sensitivity / ppv / f1",
        "TP = true rank taxon or descendant, VP = ancestor, FP = off-lineage, FN = unclassified; sensitivity = TP/(TP+VP+FN+FP), PPV = TP/(TP+FP) (Wood et al. 2019)")
    add("fda_argos_*", "fda-argos", "scope / sra_run / concordant / discordant / ancestor / unclassified",
        "scope: experiment (one SRR) | unweighted_mean (headline, mean of the experiments) | pooled (all fragments); fractions of fragments")
    add("cami_*", "cami", "group / rank / excluded_frac / n_scored / sensitivity / ppv / f1",
        "group: all | bacteria | archaea ...; excluded = reads whose true taxon has no reference in the library (reported, not scored)")
    add("mock_*", "mock", "scope / rank / min_frac / precision / recall / f1 / l1 / bray_curtis / off_target_frac",
        "profile-level agreement with the expected composition at detection threshold min_frac; L1 and Bray-Curtis distances, off-target fraction")
    return rows


OVERVIEW_FIELDS = ["sheet", "category", "column", "meaning"]


# ---------------------------------------------------------------------------
# sheets for one set of base tables
# ---------------------------------------------------------------------------

def derive_tables(base, base_dir=None, libraries=None):
    """base: {runs, classify_repeats, strain_exclusion_runs, fda_argos_runs,
    cami_runs, mock_runs, orders, datasets} -> OrderedDict name -> (fields, rows)
    for every sheet (empty detail sheets are omitted)."""
    runs = base.get("runs") or []
    t = OrderedDict()
    t["overview"] = (OVERVIEW_FIELDS, overview_rows())
    lib_rows = libraries if libraries is not None else build_libraries(base_dir, runs) if base_dir else []
    t["libraries"] = (LIBRARY_FIELDS, lib_rows)
    t["datasets"] = (DATASET_FIELDS, base.get("datasets") or [])
    orders = base.get("orders") or []
    t["orders"] = (orders_fields(orders), orders)
    summ = summary_rows(runs)
    sfields = prune_empty(summary_field_order(), summ, set(SUMMARY_ID_FIELDS))
    t["summary"] = (sfields, summ)
    rfields = prune_empty(runs_field_order(), runs, set(RUN_ID_FIELDS))
    t["runs"] = (rfields, runs)
    vr = vs_reference_rows(runs)
    if vr:
        t["vs_reference"] = (vs_reference_fields(), vr)
    corr = timing_correlation_rows(summ)
    if corr:
        t["timing_correlation"] = (CORREL_FIELDS, corr)
    if base.get("classify_repeats"):
        t["classify_repeats"] = (CLASSIFY_REPEAT_FIELDS, base["classify_repeats"])
    for name, summarize, sname in (
            ("strain_exclusion_runs", strain_exclusion_summary, "strain_exclusion_summary"),
            ("fda_argos_runs", fda_argos_summary, "fda_argos_summary"),
            ("cami_runs", cami_summary, "cami_summary"),
            ("mock_runs", mock_summary, "mock_summary")):
        rows = base.get(name) or []
        if not rows:
            continue
        sf, sr = summarize(rows)
        if sr:
            t[sname] = (sf, sr)
        fields = []
        for r in rows:
            for k in r:
                if k not in fields:
                    fields.append(k)
        t[name] = (fields, rows)
    return t


# ---------------------------------------------------------------------------
# workbook assembly (+ charts)
# ---------------------------------------------------------------------------

def _categories(sheet, fields):
    if sheet in ("runs", "summary"):
        return [field_category(f) for f in fields]
    if sheet == "orders":
        return [orders_category(f) for f in fields]
    if sheet == "overview":
        return ["reference"] * len(fields)
    if sheet == "libraries":
        return ["library"] * len(fields)
    if sheet == "datasets":
        return ["dataset"] * len(fields)
    if sheet == "vs_reference":
        return [CAT_ID if f in ("library", "dataset", "order", "ref_order", "n_seeds") else CAT_VSREF for f in fields]
    if sheet == "timing_correlation":
        return [CAT_ID if f in ("library", "dataset", "metric", "n_orders") else "correlation" for f in fields]
    if sheet == "classify_repeats":
        return [CAT_ID if f != "class_s" else CAT_CLASSIFY for f in fields]
    # detail accuracy sheets
    ids = ("library", "dataset", "order", "seed_idx", "seed", "n_seeds", "context")
    return [CAT_ID if f in ids else CAT_ACC for f in fields]


def sheet_spec(name, fields, rows, charts=None):
    spec = {"name": name, "fields": list(fields), "rows": rows,
            "categories": _categories(name, fields), "charts": charts or []}
    block = ()
    merge = ()
    if name in ("runs", "summary", "vs_reference", "timing_correlation", "classify_repeats") \
            or name.endswith("_runs") or name.endswith("_summary"):
        block = tuple(f for f in ("library", "dataset") if f in fields)
        merge = tuple(f for f in ("library", "dataset", "order", "metric") if f in fields)
        if name == "runs":
            merge = tuple(f for f in ("library", "dataset", "order") if f in fields)
    if name == "overview":
        block = ("sheet",)
        merge = ("sheet",)
        spec["wrap_fields"] = ["meaning"]
        spec["widths"] = {"sheet": 22, "category": 18, "column": 44, "meaning": 110}
        spec["freeze_cols"] = 0
    if name in ("libraries", "datasets", "orders"):
        spec["wrap_fields"] = [f for f in ("role", "fraction", "taxonomy", "sources", "truth", "provenance",
                                           "query_sets_run", "orders_run", "description", "cli") if f in fields]
        spec["widths"] = dict((f, w) for f, w in (
            ("role", 40), ("fraction", 40), ("taxonomy", 30), ("sources", 38), ("truth", 48),
            ("provenance", 52), ("query_sets_run", 28), ("orders_run", 36), ("description", 56),
            ("cli", 52), ("path", 40), ("db_path", 40), ("run_root", 40)) if f in fields)
    spec["block_fields"] = block
    spec["merge_fields"] = merge
    if name == "summary":
        spec["ratio_fields"] = dict(
            (f, METRIC_DIR[f[:-6]]) for f in fields if f.endswith("_ratio") and METRIC_DIR.get(f[:-6]))
    if name in ("runs", "summary"):
        spec["freeze_cols"] = 3
    return spec


# ---- charts (column charts with mean +/- std error bars; scatter vs classify time)

def _chart_num(v):
    if v in ("", None):
        return 0.0
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def _orders_in(rows):
    seen = []
    for r in rows:
        o = r.get("order") or ""
        if o and o not in seen:
            seen.append(o)
    return sorted(seen, key=order_key)


def _series(rows, orders, name, mean_key, std_key, match):
    means, stds = [], []
    for o in orders:
        hit = None
        for r in rows:
            if (r.get("order") or "") != o:
                continue
            if all((r.get(k) if r.get(k) is not None else "") == v for k, v in match.items()):
                hit = r
                break
        means.append(_chart_num(hit.get(mean_key) if hit else None))
        stds.append(_chart_num(hit.get(std_key) if hit and std_key else None))
    if not any(means) and not any(stds):
        return None
    return {"name": name, "means": means, "stds": stds}


def _by_dataset(rows, orders, metrics, tag):
    datasets = sorted(set(r.get("dataset") or "" for r in rows), key=dataset_key)
    charts = []
    for title, mk, sk in metrics:
        series = [s for s in (_series(rows, orders, ds, mk, sk, {"dataset": ds}) for ds in datasets if ds) if s]
        if series:
            charts.append({"title": "%s %s" % (title, tag), "y_title": title, "cats": orders, "series": series})
    return charts


def _one_series(rows, orders, title, mk, sk, tag):
    first_ds = None
    for r in rows:
        if r.get("dataset"):
            first_ds = r["dataset"]
            break
    s = _series(rows, orders, title, mk, sk, {"dataset": first_ds} if first_ds else {})
    if not s:
        return None
    return {"title": "%s %s" % (title, tag), "y_title": title, "cats": orders, "series": [s]}


def _class_times(summary, datasets):
    out = {}
    for r in summary:
        if r.get("dataset") in datasets and r["order"] not in out:
            out[r["order"]] = (_chart_num(r.get("class_s_mean")), _chart_num(r.get("class_s_std")))
    return out


def _scatter(rows, orders, title, y_title, mk, sk, class_times, match):
    series = []
    for o in orders:
        hit = None
        for r in rows:
            if (r.get("order") or "") == o and all((r.get(k) if r.get(k) is not None else "") == v for k, v in match.items()):
                hit = r
                break
        xs = class_times.get(o)
        if not hit or not xs:
            continue
        series.append({"name": o, "x": [xs[0]], "y": [_chart_num(hit.get(mk))], "stds": [_chart_num(hit.get(sk))]})
    if not series:
        return None
    return {"type": "scatter", "title": title, "x_title": "Classify wall (s)", "y_title": y_title, "series": series}


def summary_charts(summary):
    charts = []
    for lib in _dedupe(r["library"] for r in summary):
        rows = [r for r in summary if r["library"] == lib]
        orders = _orders_in(rows)
        tag = "[%s]" % lib
        for spec in (
            _one_series(rows, orders, "Build wall (s)", "build_s_mean", "build_s_std", tag),
            _one_series(rows, orders, "Build RSS (GB)", "build_gb_mean", "build_gb_std", tag),
        ):
            if spec:
                charts.append(spec)
        charts += _by_dataset(rows, orders, (("Classify wall (s)", "class_s_mean", "class_s_std"),
                                             ("Classify RSS (GB)", "class_gb_mean", "class_gb_std")), tag)
        fe = _series(rows, orders, "expected (random DNA)", "fac_expected_mean", "fac_expected_std", {"dataset": rows[0]["dataset"]})
        fp = _series(rows, orders, "particular (library)", "fac_particular_mean", "fac_particular_std", {"dataset": rows[0]["dataset"]})
        fac = [s for s in (fe, fp) if s]
        if fac:
            charts.append({"title": "Minimizer density factor " + tag, "y_title": "Density factor",
                           "y_min": 1.0, "y_max": 2.1, "cats": orders, "series": fac})
        charts += _by_dataset(rows, orders, (("Lookups / Mbp", "lookups_pm_mean", "lookups_pm_std"),), tag)
        charts += _by_dataset(rows, orders, (("NextMinimizer (s)", "nextmin_s_mean", "nextmin_s_std"),), tag)
        rk = _one_series(rows, orders, "Rank (ns)", "ns_per_minimizer_mean", "ns_per_minimizer_std", tag)
        if rk:
            charts.append(rk)
    return charts


def accuracy_charts(sheet, rows, summary):
    charts = []
    for lib in _dedupe(r["library"] for r in rows):
        lrows = [r for r in rows if r["library"] == lib]
        srows = [r for r in summary if r["library"] == lib]
        orders = _orders_in(lrows)
        tag = "[%s]" % lib
        if sheet == "strain_exclusion_summary":
            ct = _class_times(srows, ("strain-exclusion", "simulated", "labeled"))
            slices = sorted(set((r.get("domain") or "", r.get("rank") or "") for r in lrows),
                            key=lambda s: (DOMAIN_ORDER.get(s[0], 9), RANK_ORDER.get(s[1], 9), s))
            for title, mk, sk in (("Sensitivity", "sensitivity_mean", "sensitivity_std"),
                                  ("PPV", "ppv_mean", "ppv_std"), ("F1", "f1_mean", "f1_std")):
                series = []
                for dom, rank in slices:
                    s = _series(lrows, orders, "%s / %s" % (dom, rank) if dom else rank, mk, sk,
                                {"domain": dom, "rank": rank})
                    if s:
                        series.append(s)
                if series:
                    charts.append({"title": "%s %s" % (title, tag), "y_title": title, "cats": orders, "series": series})
                sc = _scatter(lrows, orders, "%s vs classify time (all / genus) %s" % (title, tag), title,
                              mk, sk, ct, {"domain": "all", "rank": "genus"})
                if sc:
                    charts.append(sc)
        elif sheet == "fda_argos_summary":
            ct = _class_times(srows, ("fda-argos", "argos"))
            for title, mk, sk in (("Concordant", "concordant_mean", "concordant_std"),
                                  ("Discordant", "discordant_mean", "discordant_std"),
                                  ("Ancestor", "ancestor_mean", "ancestor_std"),
                                  ("Unclassified", "unclassified_mean", "unclassified_std")):
                s = _series(lrows, orders, title, mk, sk, {})
                if s:
                    charts.append({"title": "%s (fraction of fragments) %s" % (title, tag),
                                   "y_title": "Fraction", "cats": orders, "series": [s]})
                sc = _scatter(lrows, orders, "%s vs classify time %s" % (title, tag), title, mk, sk, ct, {})
                if sc:
                    charts.append(sc)
        elif sheet == "cami_summary":
            ct = _class_times(srows, ("cami",))
            for title, mk, sk in (("Sensitivity", "sensitivity_mean", "sensitivity_std"),
                                  ("PPV", "ppv_mean", "ppv_std"), ("F1", "f1_mean", "f1_std")):
                series = [s for s in (_series(lrows, orders, "all / " + rk, mk, sk, {"group": "all", "rank": rk})
                                      for rk in ("genus", "species")) if s]
                if series:
                    charts.append({"title": "CAMI %s (in-library reads) %s" % (title, tag),
                                   "y_title": title, "cats": orders, "series": series})
                sc = _scatter(lrows, orders, "CAMI %s vs classify time (all / genus) %s" % (title, tag), title,
                              mk, sk, ct, {"group": "all", "rank": "genus"})
                if sc:
                    charts.append(sc)
        elif sheet == "mock_summary":
            ct = _class_times(srows, ("mock",))
            head = [r for r in lrows if _chart_num(r.get("min_frac")) == MOCK_HEADLINE_MIN_FRAC]
            for title, mk, sk in (("F1 (taxa)", "f1_mean", "f1_std"),
                                  ("Bray-Curtis", "bray_curtis_mean", "bray_curtis_std"),
                                  ("Off-target fraction", "off_target_frac_mean", "off_target_frac_std")):
                series = []
                for scope in ("all", "bacteria"):
                    for rk in ("genus", "species"):
                        s = _series(head, orders, "%s / %s" % (scope, rk), mk, sk, {"scope": scope, "rank": rk})
                        if s:
                            series.append(s)
                if series:
                    charts.append({"title": "Mock %s %s" % (title, tag), "y_title": title,
                                   "cats": orders, "series": series})
                sc = _scatter(head, orders, "Mock %s vs classify time (all / species) %s" % (title, tag), title,
                              mk, sk, ct, {"scope": "all", "rank": "species"})
                if sc:
                    charts.append(sc)
    return charts


def build_sheet_specs(tables):
    summary = tables["summary"][1] if "summary" in tables else []
    specs = []
    for name in ALL_SHEETS:
        if name not in tables:
            continue
        fields, rows = tables[name]
        charts = []
        if name == "summary":
            charts = summary_charts(rows)
        elif name in ("strain_exclusion_summary", "fda_argos_summary", "cami_summary", "mock_summary"):
            charts = accuracy_charts(name, rows, summary)
        specs.append(sheet_spec(name, fields, rows, charts))
    return specs


def write_outputs(tables, csv_dir, xlsx_path):
    """kraken2_<sheet>.csv for every sheet in csv_dir + the workbook."""
    for name, (fields, rows) in tables.items():
        write_csv(os.path.join(csv_dir, PREFIX + name + ".csv"), fields, rows)
    names = write_book(xlsx_path, build_sheet_specs(tables))
    return names


# ---------------------------------------------------------------------------
# entry points
# ---------------------------------------------------------------------------

def resolve_base_dir(run_root):
    """<base> of <base>/runs/<library> (derived from the path; KRAKEN2_BENCHMARK_DIR
    only when the run root is not laid out as runs/<library>)."""
    root = os.path.abspath(run_root)
    if os.path.basename(os.path.dirname(root)) == "runs":
        return os.path.dirname(os.path.dirname(root))
    env = os.environ.get("KRAKEN2_BENCHMARK_DIR")
    if env and os.path.isdir(os.path.join(env, "runs")):
        return env
    return None


def _remove_legacy(tables_dir):
    removed = []
    for name in LEGACY_TABLES:
        p = os.path.join(tables_dir, name)
        if os.path.isfile(p):
            try:
                os.remove(p)
                removed.append(name)
            except OSError:
                pass
    return removed


def export_library(run_root, out_xlsx, library, per_seed_rows, acc_rows, dens_rows, query_rows):
    """Write <run_root>/tables/kraken2_*.csv and out_xlsx for one library.
    Returns (xlsx_path, tables_dir, tables)."""
    tables_dir = os.path.join(run_root, "tables")
    base_dir = resolve_base_dir(run_root)
    runs, repeats = build_runs(library, per_seed_rows, acc_rows, dens_rows)
    base = {"runs": runs, "classify_repeats": repeats}
    for sheet, kind in DETAIL_KEYS.items():
        base[sheet] = detail_rows(library, acc_rows.get(kind))[1]
    base["orders"] = build_orders(library, run_root, runs, dens_rows)
    base["datasets"] = decorate_datasets(query_rows, runs)
    cache = os.path.join(base_dir, "results", ".libraries_cache.json") if base_dir else None
    libs = build_libraries(base_dir, runs, cache) if base_dir else build_libraries_standalone(library, runs)
    libs = [r for r in libs if r["library"] == library]
    tables = derive_tables(base, libraries=libs)
    for old in _remove_legacy(tables_dir):
        print("removed legacy table {}".format(old))
    names = write_outputs(tables, tables_dir, out_xlsx)
    return out_xlsx, tables_dir, tables


def build_libraries_standalone(library, runs):
    lr = [r for r in runs if r["library"] == library]
    return [OrderedDict([("library", library), ("n_runs", len(lr) or None),
                         ("query_sets_run", ", ".join(sorted(set(r["dataset"] for r in lr), key=dataset_key)) or None),
                         ("orders_run", ", ".join(sorted(set(r["order"] for r in lr), key=order_key)) or None),
                         ("status", "runs")])]


def export_results(base_dir, quiet=False):
    """Merge every <base>/runs/*/tables/kraken2_runs.csv into
    <base>/results/kraken2_results.xlsx + kraken2_*.csv. Returns (xlsx, csv_dir,
    n_libraries) or None when no library has been exported yet."""
    runs_dir = os.path.join(base_dir, "runs")
    if not os.path.isdir(runs_dir):
        return None
    parts = defaultdict(list)
    libs_found = []
    for lib in sorted(os.listdir(runs_dir), key=library_key):
        tdir = os.path.join(runs_dir, lib, "tables")
        if not os.path.isfile(os.path.join(tdir, PREFIX + "runs.csv")):
            continue
        libs_found.append(lib)
        for sheet in BASE_SHEETS:
            _f, rows = read_csv(os.path.join(tdir, PREFIX + sheet + ".csv"))
            if rows:
                parts[sheet].append(rows)
    if not libs_found:
        return None
    base = {}
    for sheet in ("runs", "classify_repeats", "strain_exclusion_runs", "fda_argos_runs",
                  "cami_runs", "mock_runs"):
        rows = [r for part in parts[sheet] for r in part]
        if sheet == "runs":
            rows.sort(key=_run_sort_key)
        base[sheet] = rows
    base["orders"] = merge_orders(parts["orders"])
    base["datasets"] = merge_datasets(parts["datasets"])
    results = os.path.join(base_dir, "results")
    tables = derive_tables(base, base_dir=base_dir,
                           libraries=build_libraries(base_dir, base["runs"],
                                                     os.path.join(results, ".libraries_cache.json")))
    write_outputs(tables, results, os.path.join(results, WORKBOOK))
    return os.path.join(results, WORKBOOK), results, len(libs_found)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Merge the per-library kraken2 tables into the results workbook")
    ap.add_argument("--base-dir", required=True, help="benchmark dir holding runs/<library>/tables/")
    args = ap.parse_args(argv)
    res = export_results(args.base_dir)
    if not res:
        print("no runs/*/tables/kraken2_runs.csv under {}".format(args.base_dir), file=sys.stderr)
        return 1
    print("wrote {} ({} libraries) and kraken2_*.csv under {}".format(res[0], res[2], res[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
