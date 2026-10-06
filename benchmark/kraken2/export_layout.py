#!/usr/bin/env python3
"""Kraken 2 benchmark export layer: sheets, CSVs and the results workbook.

Same layout and naming as the other apps (benchmark/sshash/scripts/export.py,
benchmark/kmc/scripts/export.py): one workbook, `kraken2_results.xlsx`, and one
CSV per sheet, `kraken2_<sheet>.csv`. Column names follow sshash: the raw
value carries its unit (`build_wall_s`, `classify_max_rss_gb`,
`classified_pct`), `classify`-type sheets add `<metric>_mean`, `<metric>_std`
(over seeds) and (for metrics with a direction) `<metric>_ratio` = mean / the
`xor` mean of the same library, dataset and minimum_hit_groups.

Every (seed, dataset, order) is classified once per `--minimum-hit-groups`
value (0 = the Kraken 2 paper setting, 2 = the k2 default): `minimum_hit_groups`
is an identity column / slice key of all classify-side sheets; rows are sorted
dataset -> minimum_hit_groups -> order.

  overview        column reference / legend (the role `index` has elsewhere)
  libraries       one row per Kraken library
  datasets        one row per query set
  orders          one row per minimizer scheme/order (CLI, l/k/w, spaced seed)
  build           library x order: build time / RSS / index size and the
                  density microbenchmark, mean / std / ratio over seeds
  build_runs      the per-seed rows behind `build`
  classify        library x dataset x minimum_hit_groups x order: classify-phase
                  measures (wall = DB load + processing) + accuracy
  classify_runs   the per-seed rows behind `classify`
  profile         NextMinimizer / phase profile (only for ENABLE_MINIMIZER_PROF builds)
  vs_reference    paired agreement with xor, library x dataset x mhg x order
  accuracy_<dataset> (strain_exclusion, fda_argos, cami, mock) accuracy per
                  slice (mhg / domain / group / SRR / threshold), mean / std over seeds
  paper_check     xor at minimum_hit_groups 0 on the strain_exclusion library vs
                  the published Kraken 2 values (Wood et al. 2019, Tables S1 / S5)

Build values are a property of (library, order), not of the query set, so they
live once on `build`; `classify`/`classify_runs` hold only classify-phase measures.

compare_benchmark.py collects the raw numbers of one library (one
`run_benchmark.sh` invocation) and calls `export_library()`, which writes that
library's tables under <run_root>/tables/ and its workbook
<run_root>/benchmark_comparison.xlsx; `export_results()` then merges every
library under <base>/runs/*/tables/ into <base>/results/kraken2_results.xlsx
(+ kraken2_*.csv). `classify` and `build` are always *derived* from the base
tables (`derive_tables`), so a per-library and the merged export share code.

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
from xlsx_book import write_book, DATA_START_ROW  # noqa: E402
from xlsx_util import write_csv, col_letter  # noqa: E402

PREFIX = "kraken2_"
WORKBOOK = "kraken2_results.xlsx"
SHEETS = (
    "overview",
    "libraries",
    "datasets",
    "orders",
    "build",
    "build_runs",
    "classify",
    "classify_runs",
    "profile",
    "vs_reference",
    "accuracy_strain_exclusion",
    "accuracy_fda_argos",
    "accuracy_cami",
    "accuracy_mock",
    "paper_check",
)
# accuracy sheet of each dataset kind (compare_benchmark's acc_rows key)
ACCURACY_SHEETS = OrderedDict([
    ("labeled", "accuracy_strain_exclusion"),
    ("argos", "accuracy_fda_argos"),
    ("cami", "accuracy_cami"),
    ("mock", "accuracy_mock"),
])
# Base tables: what a per-library export persists and a merge reads back
# (build / classify are derived from build_runs / classify_runs).
BASE_SHEETS = ("build_runs", "classify_runs", "profile", "vs_reference") + tuple(ACCURACY_SHEETS.values()) + (
    "orders", "datasets")

# paper_check: tolerance (absolute, on 0-1 fractions) for `within_tolerance`
PAPER_TOLERANCE = 0.01
PAPER_PUBLISHED = os.path.join(HERE, "tables", "paper_kraken2_published.csv")
PAPER_LIBRARY = "strain_exclusion"
PAPER_ORDER = REF_SCHEME
PAPER_MHG = 0

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
    "rank", "scope", "group", "sra_run", "species", "sources", "path", "gm_table",
    "cli", "description", "sketching", "rank_order", "spaced_seed_mask", "role",
    "taxonomy", "db_path", "run_root", "query_sets_run", "orders_run", "libraries",
    "status", "truth", "scoring", "truth_detail", "provenance", "paired", "read_length_bp",
    "fragments_per_source", "metric", "sheet", "column", "meaning", "category",
    "top_off_target", "fraction", "flag", "table", "row", "source", "note", "within_tolerance",
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
    """Scheme family of an order id (long or short): 'greedymini_w5_k15_canonical' -> 'greedymini'.
    The separator between family and GreedyMini table is `_`, the hyphen belongs to
    the family names (`xor` is a string prefix of `xor-pcg`, never a family match)."""
    for fam in sorted(SCHEME_ORDER, key=len, reverse=True):
        if order == fam or order.startswith(fam + "_"):
            return fam
    return order.split("_")[0]


GM_FAMILIES = ("greedymini", "mod-minimizer-greedymini")


def short_order(order):
    """Display name of an order: the run directory name without the GreedyMini
    table ('greedymini_w5_k15_canonical' -> 'greedymini'). Short names pass through."""
    for fam in GM_FAMILIES:
        if order.startswith(fam + "_"):
            return fam
    return order


def gm_table(order):
    """Table part of a long order name ('w5_k15_canonical'), else None."""
    for fam in GM_FAMILIES:
        if order.startswith(fam + "_"):
            return order[len(fam) + 1:]
    return None


def display_orders(long_names):
    """{long: short} for the run directory names of one library; two different
    long names that would collapse to the same short name keep their long name."""
    names = _dedupe(long_names)
    count = defaultdict(set)
    for n in names:
        count[short_order(n)].add(n)
    return dict((n, short_order(n) if len(count[short_order(n)]) == 1 else n) for n in names)


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
    """'standard_p10/greedymini_w5_k15_canonical/3 [fda-argos, mhg0]' ->
    (library, order, seed_idx-as-str)."""
    base = (context or "").split(" [")[0]
    parts = [p for p in base.split("/") if p]
    lib = parts[0] if parts else ""
    order = parts[1] if len(parts) > 1 else ""
    return lib, order, (parts[2] if len(parts) > 2 else "")


def _dedupe(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


# ---------------------------------------------------------------------------
# metric catalogue: one place defines columns, categories, meanings, ratio sense
# ---------------------------------------------------------------------------
# Metric = (name, category, meaning, direction, source key in compare_benchmark's
# per-run dict). direction: "lower" / "higher" = better; None = neutral (no ratio).

CAT_ID = "id / params"
CAT_BUILD = "build"
CAT_DENSITY = "density"
CAT_CLASSIFY = "classify"
CAT_MINIMIZERS = "minimizers"
CAT_SPACE = "space"
CAT_PROFILE = "profile"
CAT_VSREF = "vs reference"

# Accuracy column groups. The accuracy kinds differ per dataset and are NOT
# comparable across datasets; a column name is shared only where the definition
# is identical (strain-exclusion and CAMI both score reads against read-level truth).
BAND_READ = "accuracy: simulated reads, read-level truth (strain-exclusion; cami: library-aware)"
BAND_CAMI = "accuracy: cami (library-aware)"
BAND_ARGOS = "accuracy: fda-argos (real, genus concordance)"
BAND_MOCK = "accuracy: mock (real, known composition)"
ACC_BANDS = (BAND_READ, BAND_CAMI, BAND_ARGOS, BAND_MOCK)

# Column groups in the SSHash order: id / params, density, minimizers, space, build.
BUILD_METRICS = [
    ("expected_density_factor", CAT_DENSITY, "minimizer density factor = density x (w+1) on random i.i.d. DNA (Monte Carlo; 2 = random-minimizer theory), from tables/minimizer_density.csv joined to the run by order and seed_idx", "lower", None),
    ("particular_density_factor", CAT_DENSITY, "density factor on this library's FASTAs (benchmark_minimizers.sh -d). This is the density column the correlation rows of `build` / `classify` use", "lower", None),
    ("distinct_minimizers", CAT_MINIMIZERS, "number of distinct minimizers (keys) in the built database = the `size` field of the 32-byte hash.k2d header (read without loading the table). A lower density gives fewer distinct minimizers", "lower", None),
    ("index_size_gb", CAT_SPACE, "size of the built hash table (db/hash.k2d) in GB; blank when the db is not on disk", "lower", None),
    ("build_wall_s", CAT_BUILD, "k2 build wall time (s)", "lower", "build_s"),
    ("build_max_rss_gb", CAT_BUILD, "k2 build peak RSS (GB)", "lower", "build_gb"),
]
# build properties repeated (per order and seed) on classify / classify_runs
ORDER_PROPS = ["expected_density_factor", "particular_density_factor", "distinct_minimizers"]
CLASSIFY_METRICS = [
    ("classify_wall_s", CAT_CLASSIFY, "classify wall time (s) of the single classify run of this seed (/usr/bin/time -v): database load + processing. classify: mean over seeds", "lower", "class_s"),
    ("db_load_s", CAT_CLASSIFY, "database load time (s) = classify_wall_s - classify_processing_s (reading the hash table / taxonomy before the first read). Depends on the page cache of the host: judge it over the seeds", "lower", "db_load_s"),
    ("classify_processing_s", CAT_CLASSIFY, "time classifying the reads (s) = the `processed in` seconds of the classify summary line, without the database load", "lower", "class_proc_s"),
    ("classify_max_rss_gb", CAT_CLASSIFY, "classify peak RSS (GB) of the run", "lower", "class_gb"),
    ("classify_mreads_per_min", CAT_CLASSIFY, "classification speed in million reads (read pairs) per minute of classify_processing_s (database load excluded; wall time if the processing time is missing)", "higher", "mrd_min"),
    ("classified_pct", CAT_CLASSIFY, "percent of reads classified (any taxon). Neutral: more classified is not better or worse by itself (it can be false positives)", None, "cls_pct"),
]
RUN_ONLY = [
    ("n_reads", CAT_CLASSIFY, "reads (read pairs) in the query set", "n_reads"),
]
# profile sheet: NextMinimizer / phase instrumentation (classify built with
# ENABLE_MINIMIZER_PROF=1)
PROFILE_METRICS = [
    ("nextminimizer_s", CAT_PROFILE, "time in NextMinimizer / Rank (s)", "lower", "nextmin_s"),
    ("nextminimizer_ns_per_call", CAT_PROFILE, "NextMinimizer cost in ns per call", "lower", "nextmin_ns"),
    ("nextminimizer_pct_of_classify", CAT_PROFILE, "NextMinimizer time as % of classify wall", "lower", "nextmin_frac"),
    ("lookups_per_mbp", CAT_PROFILE, "hash-table lookups per Mbp of reads", "lower", "lookups_pm"),
    ("lookup_pct_of_classify", CAT_PROFILE, "lookup time as % of classify time (per-thread sums)", "lower", "lookup_frac"),
    ("scan_s", CAT_PROFILE, "scan phase (minimizer + tokens), sum of per-thread seconds", "lower", "scan_s"),
    ("lookup_s", CAT_PROFILE, "lookup phase (GetBatch), sum of per-thread seconds", "lower", "lookup_s"),
    ("resolve_s", CAT_PROFILE, "resolve phase (token replay + ResolveTree), sum of per-thread seconds", "lower", "resolve_s"),
    ("nextminimizer_calls_per_mbp", CAT_PROFILE, "NextMinimizer calls per Mbp of reads. Neutral: a workload count", None, "nextmin_cpm"),
    ("dup_key_hits_per_mbp", CAT_PROFILE, "consecutive identical minimizers per Mbp of reads, resolved without a hash lookup (TOK_REPEAT). Neutral: a workload count", None, "repeats_pm"),
    ("uhs_pct", CAT_PROFILE, "% of Rank() calls on the finite-rank (UHS) path; GreedyMini orders only. Neutral: a path split, not a quality", None, "uhs_frac"),
    ("non_uhs_pct", CAT_PROFILE, "% of Rank() calls on the infinity early-out path. Neutral", None, "non_uhs_frac"),
]
# paired agreement with xor of the same seed (vs_reference sheet only; no ratio)
VSREF_METRICS = [
    ("genus_mape_pct", CAT_VSREF, "paired MAPE of genus abundance vs xor of the same seed_idx, dataset and minimum_hit_groups (%, abundance floor --min-abund)", None, "mape_g"),
    ("species_mape_pct", CAT_VSREF, "same at species rank", None, "mape_s"),
    ("genus_l1_pct", CAT_VSREF, "L1 distance of genus abundance vs xor (%)", None, "l1_g"),
    ("concordance_pct", CAT_VSREF, "ClassConc: reads with identical taxid vs xor, % of all reads", None, "concordance"),
    ("classified_concordance_pct", CAT_VSREF, "identical taxid among reads classified by both (%)", None, "classified_concordance"),
    ("jaccard_pct", CAT_VSREF, "Jaccard of the sets of classified read ids vs xor (%)", None, "jaccard"),
    ("newly_classified_pct", CAT_VSREF, "reads unclassified by xor but classified by this order (% of all reads)", None, "newly_classified_pct"),
    ("newly_unclassified_pct", CAT_VSREF, "reads classified by xor but unclassified by this order (% of all reads)", None, "newly_unclassified_pct"),
]

# Headline accuracy columns of classify_runs / classify. One dict per column:
# name, band, meaning, direction, selectors {dataset: (detail kind, row match, source field)}.
# Every row of classify_runs / classify is one dataset, so the names carry no dataset prefix.
ACC_COLS = []


def _acc(name, band, meaning, direction, selectors):
    ACC_COLS.append({"name": name, "band": band, "meaning": meaning,
                     "dir": direction, "sel": selectors})


for _rank in ("genus", "species"):
    for _m, _src in (("sensitivity", "sensitivity"), ("ppv", "ppv"), ("f1", "f1")):
        _acc("%s_%s" % (_rank, _m), BAND_READ,
             "%s %s against read-level truth. strain-exclusion: domain all, Mason held-out genomes; "
             "cami: group all, over the reads whose true taxon is in the library" % (_rank, _m),
             "higher",
             {"strain-exclusion": ("labeled", {"domain": "all", "rank": _rank}, _src),
              "cami": ("cami", {"group": "all", "rank": _rank}, _src)})
for _rank in ("genus", "species"):
    _acc("%s_excluded_fraction" % _rank, BAND_CAMI,
         "cami: fraction of reads whose true %s has no sequence in the library (reported, not scored). "
         "Neutral: a property of library x query set, not of the order" % _rank, None,
         {"cami": ("cami", {"group": "all", "rank": _rank}, "excluded_frac")})
for _k, _dir in (("concordant", "higher"), ("discordant", "lower"),
                 ("ancestor", None), ("unclassified", None)):
    _acc("genus_%s_fraction" % _k, BAND_ARGOS,
         "fda-argos: fraction of fragments whose genus call is %s with the SRA-assigned taxon, "
         "unweighted mean over the experiments (no read-level truth)%s" % (
             _k, "" if _dir else ". Neutral: neither right nor wrong"),
         _dir, {"fda-argos": ("argos", {"rank": "genus", "scope": "unweighted_mean"}, _k + "_frac")})
for _rank in ("genus", "species"):
    for _m, _src, _dir in (("taxa_f1", "f1", "higher"), ("bray_curtis", "bray_curtis", "lower"),
                           ("off_target_fraction", "off_target_frac", "lower")):
        _acc("%s_%s" % (_rank, _m), BAND_MOCK,
             "mock (ZymoBIOMICS D6300), scope all, detection threshold min_frac=1e-4: %s at %s rank. "
             "A taxa-detection / profile metric, not read-level (so not named like the read-level f1)" % (_m, _rank),
             _dir, {"mock": ("mock", {"scope": "all", "rank": _rank, "min_frac": MOCK_HEADLINE_MIN_FRAC}, _src)})

RUN_ID_FIELDS = ["run_id", "library", "dataset", "minimum_hit_groups", "order", "seed_idx", "seed"]
BUILD_ID_FIELDS = ["run_id", "library", "order", "seed_idx", "seed"]
AGG_ID = {"classify": ["library", "dataset", "minimum_hit_groups", "order", "n_seeds"],
          "build": ["library", "order", "n_seeds"]}

METRIC_DIR = {}
METRIC_CAT = {}
for _m in BUILD_METRICS + CLASSIFY_METRICS + PROFILE_METRICS + VSREF_METRICS:
    METRIC_DIR[_m[0]] = _m[3]
    METRIC_CAT[_m[0]] = _m[1]
for _m in RUN_ONLY:
    METRIC_CAT[_m[0]] = _m[1]
for _c in ACC_COLS:
    METRIC_DIR[_c["name"]] = _c["dir"]
    METRIC_CAT[_c["name"]] = _c["band"]
BUILD_NAMES = [m[0] for m in BUILD_METRICS]
SUMMARY_NAMES = ORDER_PROPS + [m[0] for m in CLASSIFY_METRICS] + [c["name"] for c in ACC_COLS]
PROFILE_NAMES = [m[0] for m in PROFILE_METRICS]
VSREF_NAMES = [m[0] for m in VSREF_METRICS]
SUMMARY_CATS = (CAT_DENSITY, CAT_MINIMIZERS, CAT_CLASSIFY) + ACC_BANDS


def runs_field_order():
    return (RUN_ID_FIELDS + ORDER_PROPS + [m[0] for m in CLASSIFY_METRICS] + [m[0] for m in RUN_ONLY]
            + [c["name"] for c in ACC_COLS])


def build_runs_field_order():
    return BUILD_ID_FIELDS + BUILD_NAMES


def agg_field_order(ids, metrics, cat_order):
    """ids, then per metric `_mean`, `_std` and -- for metrics with a direction --
    `_ratio`, grouped by category."""
    fields = list(ids)
    by_cat = OrderedDict()
    for m in metrics:
        by_cat.setdefault(METRIC_CAT[m], []).append(m)
    for cat in cat_order:
        for m in by_cat.get(cat, []):
            fields += [m + "_mean", m + "_std"] + ([m + "_ratio"] if METRIC_DIR.get(m) else [])
    return fields


def summary_field_order():
    return agg_field_order(AGG_ID["classify"], SUMMARY_NAMES, SUMMARY_CATS)


def build_field_order():
    return agg_field_order(AGG_ID["build"], BUILD_NAMES, (CAT_DENSITY, CAT_MINIMIZERS, CAT_SPACE, CAT_BUILD))


def profile_field_order():
    return agg_field_order(AGG_ID["classify"], PROFILE_NAMES, (CAT_PROFILE,))


def vs_reference_field_order():
    return agg_field_order(AGG_ID["classify"], VSREF_NAMES, (CAT_VSREF,))


def field_category(field):
    if field in RUN_ID_FIELDS or field in BUILD_ID_FIELDS or field in AGG_ID["classify"]:
        return CAT_ID
    base = field
    for suf in ("_mean", "_std", "_ratio"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    return METRIC_CAT.get(field) or METRIC_CAT.get(base) or ""


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




# ---------------------------------------------------------------------------
# per-run rows from the collected per-seed metrics
# ---------------------------------------------------------------------------

def _index_size_gb(run_root, order, seed_idx):
    if not run_root:
        return None
    p = os.path.join(run_root, order, str(seed_idx), "db", "hash.k2d")
    try:
        return round(os.path.getsize(p) / 1e9, 4)
    except OSError:
        return None


def _distinct_minimizers(run_root, order, seed_idx):
    """`size` of the hash.k2d header (4 little-endian uint64: capacity, size,
    key_bits, value_bits), read from the first 32 bytes only."""
    if not run_root:
        return None
    p = os.path.join(run_root, order, str(seed_idx), "db", "hash.k2d")
    try:
        with open(p, "rb") as f:
            head = f.read(32)
        if len(head) < 32:
            return None
        return struct.unpack("<4Q", head)[1]
    except (IOError, OSError, struct.error):
        return None


def order_map(per_seed_rows, acc_rows):
    """{run directory name: display name} for every order of one library."""
    longs = [split_context(r.get("context"))[1] for r in per_seed_rows]
    for rows in acc_rows.values():
        longs += [split_context(r.get("context"))[1] for r in rows]
    return display_orders([x for x in longs if x])


def build_tables(library, per_seed_rows, acc_rows, dens_rows, run_root=None, omap=None):
    """per_seed_rows: dicts from compare_benchmark (scheme, dataset, seed_idx,
    context, seed + collect_run_metrics keys). acc_rows: {'labeled'|'argos'|
    'cami'|'mock': raw metric rows with context/dataset/seed_idx}. dens_rows:
    density microbenchmark rows ({scheme, seed_idx, fac_*}); the scheme is the
    family (xor, greedymini, ...) or an order name, long or short.

    Returns an OrderedDict of raw (per-seed) tables: runs, build_runs,
    profile_raw, vs_reference_raw."""
    omap = omap if omap is not None else order_map(per_seed_rows, acc_rows)
    dens = {}
    for d in dens_rows:
        dens[(d.get("scheme"), _int(d.get("seed_idx")))] = d
    idx = defaultdict(lambda: defaultdict(list))
    for kind, rows in acc_rows.items():
        for r in rows:
            _lib, order, _ = split_context(r.get("context"))
            order = omap.get(order, short_order(order))
            key = (order, r.get("dataset") or _kind_dataset(kind), _int(r.get("minimum_hit_groups")),
                   _int(r.get("seed_idx")))
            idx[key][kind].append(r)
    rfields = runs_field_order()
    runs, builds, prof, vsr = [], [], [], []
    build_cache = {}
    for row in per_seed_rows:
        _lib, long_order, _ = split_context(row.get("context"))
        order = omap.get(long_order, short_order(long_order))   # display name; paths use long_order
        ds = row.get("dataset") or "speed"
        mhg = _int(row.get("minimum_hit_groups"))
        seed_idx = _int(row.get("seed_idx"))
        ident = OrderedDict([
            ("run_id", "%s/%s/%s" % (library, order, seed_idx)), ("library", library),
            ("dataset", ds), ("minimum_hit_groups", mhg), ("order", order), ("seed_idx", seed_idx),
            ("seed", _seed_hex(row.get("seed")))])
        # build row: one per (order, seed), whatever the number of datasets / minimum_hit_groups
        if (order, seed_idx) not in build_cache:
            b = OrderedDict((f, ident.get(f)) for f in build_runs_field_order())
            for m in BUILD_METRICS:
                if m[4]:
                    b[m[0]] = clean_value(row.get(m[4]))
            b["distinct_minimizers"] = _distinct_minimizers(run_root, long_order, seed_idx)
            b["index_size_gb"] = _index_size_gb(run_root, long_order, seed_idx)
            # per-seed density row, else the scheme-level one (density summary text: no seed_idx)
            d = (dens.get((family(order), seed_idx)) or dens.get((order, seed_idx))
                 or dens.get((long_order, seed_idx)) or dens.get((family(order), None)))
            if d:
                b["expected_density_factor"] = clean_value(d.get("fac_expected"))
                b["particular_density_factor"] = clean_value(d.get("fac_particular"))
            build_cache[(order, seed_idx)] = b
            builds.append(b)
        b = build_cache[(order, seed_idx)]
        out = OrderedDict((f, ident.get(f)) for f in rfields)
        for name in ORDER_PROPS:
            out[name] = b.get(name)
        for m in CLASSIFY_METRICS:
            out[m[0]] = clean_value(row.get(m[4]))
        for m in RUN_ONLY:
            out[m[0]] = clean_value(row.get(m[3]))
        picked = idx.get((order, ds, mhg, seed_idx), {})
        for c in ACC_COLS:
            sel = c["sel"].get(ds)
            if sel:
                out[c["name"]] = clean_value(_detail_pick(picked, sel, ds))
        runs.append(out)
        for names, dest in ((PROFILE_METRICS, prof), (VSREF_METRICS, vsr)):
            o2 = OrderedDict((f, ident.get(f)) for f in RUN_ID_FIELDS)
            for m in names:
                o2[m[0]] = clean_value(row.get(m[4]))
            if any(o2[m[0]] is not None for m in names):
                dest.append(o2)
    for rows in (runs, builds, prof, vsr):
        rows.sort(key=_run_sort_key)
    return OrderedDict([("classify_runs", runs), ("build_runs", builds), ("profile_raw", prof),
                        ("vs_reference_raw", vsr)])


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
            r.get("minimum_hit_groups") or 0, order_key(r.get("order") or ""), r.get("seed_idx") or 0)


# ---------------------------------------------------------------------------
# aggregation: classify / build / profile / vs_reference
# ---------------------------------------------------------------------------

def aggregate(rows, metrics, keys=("library", "dataset", "minimum_hit_groups", "order")):
    """One row per `keys` group: mean / std over seeds of every metric, then the
    ratio to the reference order (xor) of the same library, dataset and
    minimum_hit_groups."""
    groups = OrderedDict()
    for r in rows:
        groups.setdefault(tuple(r.get(k) for k in keys), []).append(r)
    out = []
    for key, parts in groups.items():
        row = OrderedDict(zip(keys, key))
        row["n_seeds"] = len(set(p.get("seed_idx") for p in parts))
        for m in metrics:
            mu, sd = mean_std([p.get(m) for p in parts if _num(p.get(m))])
            row[m + "_mean"] = round_sig(mu)
            row[m + "_std"] = round_std(sd)
        out.append(row)
    out.sort(key=lambda r: (library_key(r["library"]), dataset_key(r.get("dataset") or ""),
                            r.get("minimum_hit_groups") or 0, order_key(r["order"])))
    attach_ratios(out, metrics, keys)
    return out


def attach_ratios(rows, metrics, keys=("library", "dataset", "minimum_hit_groups", "order")):
    """{metric}_ratio = this order's mean / the reference (xor) mean of the
    same library, dataset and minimum_hit_groups; the reference = 1. Only for metrics with a
    direction; blank when either mean is missing or the reference mean is 0."""
    group = [k for k in keys if k != "order"]
    base = {}
    for r in rows:
        if family(r["order"]) == REF_SCHEME:
            base[tuple(r.get(k) for k in group)] = r
    for r in rows:
        b = base.get(tuple(r.get(k) for k in group))
        for m in metrics:
            if not METRIC_DIR.get(m):
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


def vs_reference_aggregate(raw):
    """Paired metrics per non-reference order (the reference compares to itself)."""
    return aggregate([r for r in raw if family(r["order"]) != REF_SCHEME], VSREF_NAMES)


# ---------------------------------------------------------------------------
# accuracy_detail: every per-seed accuracy row, shared column names
# ---------------------------------------------------------------------------

MOCK_RENAME = {"precision": "taxa_precision", "recall": "taxa_recall", "f1": "taxa_f1",
               "l1": "l1_distance"}
DETAIL_ID = ("library", "dataset", "order", "seed_idx", "seed")
# slice keys of the accuracy sheets; minimum_hit_groups is the outermost one
# (dataset -> minimum_hit_groups -> slice -> order)
DETAIL_SLICE = ("minimum_hit_groups", "scope", "domain", "group", "sra_run", "species", "min_frac", "rank")


def detail_name(kind, key):
    if kind == "mock" and key in MOCK_RENAME:
        return MOCK_RENAME[key]
    if key.endswith("_frac") and key != "min_frac":
        return key[:-5] + "_fraction"
    return key


def detail_rows(library, acc_rows, omap=None):
    """{accuracy sheet: rows} -- every per-seed accuracy row of each dataset kind
    (strain-exclusion, fda-argos, cami, mock): library / dataset / order / seed +
    the slice (rank, domain, group, scope, SRR, threshold) + counts and metrics."""
    omap = omap if omap is not None else order_map([], acc_rows)
    result = OrderedDict()
    for kind, sheet in ACCURACY_SHEETS.items():
        out = []
        for r in acc_rows.get(kind) or []:
            _lib, order, _ = split_context(r.get("context"))
            row = OrderedDict()
            row["library"] = library
            row["dataset"] = r.get("dataset") or _kind_dataset(kind)
            row["order"] = omap.get(order, short_order(order))
            for k, v in r.items():
                if k in ("scheme", "context", "dataset", "order", "library"):
                    continue
                row[detail_name(kind, k)] = _seed_hex(v) if k == "seed" else clean_value(v)
            out.append(row)
        out.sort(key=lambda r: (_int(r.get("minimum_hit_groups")) or 0, order_key(r["order"]), _int(r.get("seed_idx")) or 0,
                                DOMAIN_ORDER.get(r.get("domain"), 9), RANK_ORDER.get(r.get("rank"), 9),
                                str(r.get("scope") or ""), str(r.get("sra_run") or ""),
                                str(r.get("group") or ""),
                                str(r.get("min_frac") if r.get("min_frac") is not None else "")))
        result[sheet] = out
    return result


COUNT_FIELDS = frozenset([
    "n", "n_truth", "n_calls", "n_scored", "tp", "fp", "fn", "vp", "no", "excluded",
    "ex_unclassified", "ex_consistent", "ex_off", "n_reads", "n_expected", "n_detected",
    "tp_taxa", "fp_taxa", "fn_taxa"])
# direction of the accuracy metrics (ratio + colour); everything else is neutral
DETAIL_DIR = {
    "sensitivity": "higher", "ppv": "higher", "f1": "higher",
    "concordant_fraction": "higher", "discordant_fraction": "lower",
    "taxa_precision": "higher", "taxa_recall": "higher", "taxa_f1": "higher",
    "on_target_fraction": "higher", "bray_curtis": "lower", "l1_distance": "lower",
    "off_target_fraction": "lower",
}
SCOPE_ORDER = {"unweighted_mean": 0, "all": 0, "pooled": 1, "bacteria": 1, "experiment": 2}
DETAIL_TEXT_DROP = frozenset(["top_off_target"])


def _slice_sort(r):
    mf = r.get("min_frac")
    return (r.get("minimum_hit_groups") if _num(r.get("minimum_hit_groups")) else -1,
            SCOPE_ORDER.get(r.get("scope"), 9), DOMAIN_ORDER.get(r.get("domain"), 9),
            str(r.get("group") or ""), str(r.get("sra_run") or ""), str(r.get("species") or ""),
            mf if _num(mf) else -1.0, RANK_ORDER.get(r.get("rank"), 9))


def detail_aggregate(rows):
    """Per-seed accuracy rows (detail_rows) -> one row per library x dataset x
    slice x order: <metric>_mean / _std over seeds (counts: _mean only) and
    <metric>_ratio to xor of the same slice (incl. minimum_hit_groups) for metrics
    with a direction. Sorted dataset -> minimum_hit_groups -> slice -> order, so
    a slice's orders are consecutive."""
    skip = set(DETAIL_ID) | set(DETAIL_SLICE) | DETAIL_TEXT_DROP
    metrics = []
    for r in rows:
        for k, v in r.items():
            if k not in skip and k not in metrics and _num(v):
                metrics.append(k)
    groups = OrderedDict()
    for r in rows:
        key = (r["library"], r["dataset"]) + tuple(r.get(k) for k in DETAIL_SLICE) + (r["order"],)
        groups.setdefault(key, []).append(r)
    out = []
    for key, parts in groups.items():
        row = OrderedDict([("library", key[0]), ("dataset", key[1])])
        for k, v in zip(DETAIL_SLICE, key[2:2 + len(DETAIL_SLICE)]):
            row[k] = v
        row["order"] = key[-1]
        row["n_seeds"] = len(set(p.get("seed_idx") for p in parts))
        for m in metrics:
            mu, sd = mean_std([p.get(m) for p in parts if _num(p.get(m))])
            row[m + "_mean"] = round_sig(mu)
            if m not in COUNT_FIELDS:
                row[m + "_std"] = round_std(sd)
        out.append(row)
    out.sort(key=lambda r: (library_key(r["library"]), dataset_key(r["dataset"]), _slice_sort(r),
                            order_key(r["order"])))
    base = {}
    for r in out:
        if family(r["order"]) == REF_SCHEME:
            base[(r["library"], r["dataset"]) + tuple(r.get(k) for k in DETAIL_SLICE)] = r
    for r in out:
        b = base.get((r["library"], r["dataset"]) + tuple(r.get(k) for k in DETAIL_SLICE))
        for m in metrics:
            if m in DETAIL_DIR:
                mu = r.get(m + "_mean")
                bm = b.get(m + "_mean") if b else None
                r[m + "_ratio"] = round(mu / bm, 4) if (_num(mu) and _num(bm) and bm) else None
    return out


def detail_fields(rows):
    """Columns of an aggregated accuracy sheet: library, dataset, the slice keys
    this dataset uses, order, n_seeds, then the metrics (columns blank in every
    row are dropped)."""
    seen = []
    for r in rows:
        for k, v in r.items():
            if k not in seen and v not in (None, ""):
                seen.append(k)
    lead = (["library", "dataset"] + [k for k in DETAIL_SLICE if k in seen] + ["order", "n_seeds"])
    for k in lead:
        if k not in seen:
            seen.append(k)
    rest = []
    for f in seen:
        if f in lead:
            continue
        base = re.sub(r"_(mean|std|ratio)$", "", f)
        if base not in rest:
            rest.append(base)
    out = list(lead)
    for base in rest:
        out += [c for c in (base + "_mean", base + "_std", base + "_ratio") if c in seen]
    return out + [f for f in seen if f not in out]


# ---------------------------------------------------------------------------
# orders
# ---------------------------------------------------------------------------

ORDER_INFO = {
    "xor": ("minimizer", "xor", "--sketching minimizer --order xor",
            "random order: the l-mer XORed with the seed (toggle mask), ranked by the hash of that value"),
    "xor-pcg": ("minimizer", "xor-pcg", "--sketching minimizer --order xor-pcg",
                "random order with a PCG scramble of the seed-XORed l-mer"),
    "greedymini": ("minimizer", "greedymini", "--sketching minimizer --order greedymini --greedymini-file <table>.gm",
                   "GreedyMini table (centre extension), ties broken by the seed"),
    "mod-minimizer-xor": ("mod-minimizer", "xor", "--sketching mod-minimizer --order xor",
                          "mod-minimizer (Groot Koerkamp & Pibiri) on a random (xor) order"),
    "mod-minimizer-greedymini": ("mod-minimizer", "greedymini",
                                 "--sketching mod-minimizer --order greedymini --greedymini-file <table>.gm",
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


def build_orders(library, run_root, runs, dens_rows, omap=None):
    """One row per order seen in `runs`: CLI, l / k / w, spaced seed and the
    density factors (mean / std over the density microbenchmark rows of the
    order: expected = random DNA, particular = this library's FASTAs; one
    particular column pair per library). Density is also on `build`."""
    long_of = dict((short, long_) for long_, short in (omap or {}).items())
    order_ids = _dedupe(r["order"] for r in runs)
    rows = []
    for order in sorted(order_ids, key=order_key):
        long_order = long_of.get(order, order)        # run directory name
        fam = family(order)
        sketching, rank_order, cli, desc = ORDER_INFO.get(fam, ("", "", "", ""))
        table = gm_table(long_order)
        if table:
            cli = cli.replace("<table>", "<orders>/GreedyMini/" + table)
        row = OrderedDict([
            ("order", order), ("scheme", fam), ("sketching", sketching),
            ("rank_order", rank_order), ("gm_table", table), ("cli", cli),
            ("description", desc)])
        l = k = None
        for d in dens_rows:
            if d.get("scheme") in (order, long_order, fam) and _num(d.get("l")) and _num(d.get("k")):
                l, k = int(d["l"]), int(d["k"])
                break
        mask = None
        opts = order_params_from_db(run_root, long_order) if run_root else None
        if opts:
            k, l, mask = opts[0], opts[1], opts[2]
        row["l"] = l
        row["k"] = k
        row["w"] = (k - l + 1) if (k is not None and l is not None) else None
        spaces = None
        if mask is not None and l:
            spaces = l - bin(mask).count("1") // 2 if mask else 0
        row["minimizer_spaces"] = spaces
        row["spaced_seed_mask"] = ("0x%x" % mask) if mask else ("none" if mask == 0 else None)
        mine = [d for d in dens_rows if d.get("scheme") in (order, long_order, fam)]
        for col, src in (("expected_density_factor", "fac_expected"),
                         (library + "_particular_density_factor", "fac_particular")):
            mu, sd = mean_std([d.get(src) for d in mine if _num(d.get(src))])
            row[col + "_mean"] = round_sig(mu)
            row[col + "_std"] = round_std(sd)
        rows.append(row)
    return rows


ORDERS_FIELDS = ["order", "scheme", "sketching", "rank_order", "gm_table", "cli", "description",
                 "l", "k", "w", "minimizer_spaces", "spaced_seed_mask",
                 "expected_density_factor_mean", "expected_density_factor_std"]


def orders_fields(rows):
    """ORDERS_FIELDS plus one particular-density column pair per library present."""
    libs = []
    for r in rows:
        for k in r:
            m = re.match(r"^(.*)_particular_density_factor_mean$", k)
            if m and m.group(1) not in libs:
                libs.append(m.group(1))
    libs.sort(key=library_key)
    out = list(ORDERS_FIELDS)
    for lib in libs:
        out += [lib + "_particular_density_factor_mean", lib + "_particular_density_factor_std"]
    return out


def orders_category(field):
    if field in ("order", "scheme", "sketching", "rank_order", "gm_table", "cli", "description"):
        return "order"
    if field.startswith("expected_density_factor"):
        return "density (random DNA)"
    m = re.match(r"^(.*)_particular_density_factor_", field)
    if m:
        return "density (%s)" % m.group(1)
    return "minimizer params"


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


# ---------------------------------------------------------------------------
# datasets (query sets) and libraries
# ---------------------------------------------------------------------------

# dataset -> (sources, truth, scoring, truth_detail, provenance)
DATASET_INFO = {
    "strain-exclusion": (
        "simulated reads (Mason) from the 50 Table S3 held-out genomes",
        "simulated", "read-level truth",
        "per-read taxid: the source genome (truth.tsv); genus/species sensitivity / PPV / F1",
        "Wood et al. 2019 Table S1; Mason 2 reads from Zenodo (download_data.sh / prepare_strex_data.sh)"),
    "fda-argos": (
        "real Illumina SRA runs (FDA-ARGOS, Table S5)",
        "real", "concordance with SRA taxon",
        "per-experiment species (Table S5), no per-read truth; genus concordance a-d, no sensitivity/PPV",
        "Wood et al. 2019 Table S5; ENA download of the sra_run accessions (prepare_fda_argos_reads.py)"),
    "cami": (
        "simulated reads (ART) of CAMI II marine genomes",
        "simulated", "library-aware read-level truth",
        "per-read true taxid (truth.tsv); sens / PPV / F1 over reads whose taxon is in the library, excluded reads reported",
        "CAMI II marine sample, subsampled (prepare_cami_reads.py, SOURCE.txt)"),
    "mock": (
        "real reads of the ZymoBIOMICS D6300 mock community",
        "real", "profile vs expected composition",
        "expected community composition only (composition.csv); profile-level precision/recall/F1, L1, Bray-Curtis, off-target",
        "SRR12324251, subsampled (prepare_mock_reads.py, SOURCE.txt)"),
    "speed": (
        "unlabeled NCBI-library speed corpus",
        "real (unlabeled)", "none (timing only)",
        "no truth",
        "data/speed/ or legacy data/benchmark_R*.fastq"),
}

DATASET_FIELDS = ["dataset", "truth", "scoring", "sources", "n_sources", "fragments_per_source",
                  "n_fragments", "read_length_bp", "paired", "total_bases_r1_r2", "truth_detail",
                  "provenance", "libraries", "path"]


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
            row["scoring"] = info[2]
            row["truth_detail"] = info[3]
            row["provenance"] = info[4]
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
    """One row per library that has rows in `runs` (a db/<id> dir or an empty
    runs/<id> without results does not list: libraries without runs, e.g. the
    standard ones, stay out of the workbook)."""
    cache = _load_cache(cache_path) if cache_path else {}
    ids = sorted(set(r["library"] for r in runs), key=library_key)
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
        keep = set(os.path.join(base, "db", r["library"]) for r in rows)
        _save_cache(cache_path, dict((k, v) for k, v in cache.items() if k in keep))
    return rows


# ---------------------------------------------------------------------------
# overview (column reference)
# ---------------------------------------------------------------------------

SHEET_DOC = OrderedDict([
    ("overview", "this sheet: what every sheet and column means, grouped by category"),
    ("libraries", "one row per Kraken library: size, taxonomy, role in the benchmark, which query sets / orders ran on it"),
    ("datasets", "one row per query set: reads, read length, truth type, how it is scored, provenance. Every dataset is classified twice, with `--minimum-hit-groups` 0 (the Kraken 2 paper setting) and 2 (the k2 default); see `minimum_hit_groups`"),
    ("orders", "one row per minimizer scheme/order: CLI flags, l / k / w, spaced seed and the density factors (also on `build`: the duplication is deliberate, `orders` is the per-order reference card)"),
    ("build", "library x order, mean / std over seeds and `_ratio` to xor: build time, build RSS, index size and the density microbenchmark. Build values do not depend on the query set or on minimum_hit_groups, so they appear here once, not on every dataset row"),
    ("build_runs", "the per-seed rows behind `build` (same split as summary / runs in the sshash / kmc apps)"),
    ("classify", "library x dataset x minimum_hit_groups x order, classify phase only: mean / std over seeds and `_ratio` to xor (same dataset and minimum_hit_groups) of classify wall / DB-load / processing time, RSS, speed, classified % and the headline accuracy of the dataset"),
    ("classify_runs", "one row per seed x library x dataset x minimum_hit_groups x order: the raw numbers of the single classify run behind `classify` (the replicates are the seeds: one run per seed, dataset and minimum_hit_groups)"),
    ("profile", "NextMinimizer / phase profile, only present when the run used a classify built with ENABLE_MINIMIZER_PROF=1 (nextminimizer_timing_<dataset>_mhg<m>.txt exist): library x dataset x minimum_hit_groups x order, mean / std over seeds and `_ratio` to xor"),
    ("vs_reference", "paired agreement with xor of the same seed, dataset and minimum_hit_groups (MAPE, concordance, Jaccard, newly classified / unclassified): library x dataset x minimum_hit_groups x order, mean / std over seeds"),
    ("accuracy_strain_exclusion", "strain-exclusion accuracy per slice (minimum_hit_groups x domain x rank): one row per order, mean / std over seeds and `_ratio` to xor of the same slice (TP / VP / FP / FN counts as means)"),
    ("accuracy_fda_argos", "FDA-ARGOS accuracy per slice (minimum_hit_groups x scope: unweighted mean / pooled / each SRR experiment): genus concordance, one row per order, mean / std over seeds and `_ratio` to xor"),
    ("accuracy_cami", "CAMI accuracy per slice (minimum_hit_groups x group x rank): sens / PPV / F1 and the excluded fraction, one row per order, mean / std over seeds and `_ratio` to xor"),
    ("accuracy_mock", "mock-community accuracy per slice (minimum_hit_groups x scope x detection threshold x rank): taxa precision / recall / F1, distances, off-target, one row per order, mean / std over seeds and `_ratio` to xor"),
    ("paper_check", "paper equivalence: xor at minimum_hit_groups 0 on the strain_exclusion library vs the published Kraken 2 values (Wood et al. 2019, Table S1 genus sensitivity / PPV / F1 of the prokaryote and viral domains; Table S5 genus concordant / discordant / ancestor / unclassified, unweighted mean over experiments). Published values: benchmark/kraken2/tables/paper_kraken2_published.csv. `within_tolerance` = |ours - published| <= the stated absolute tolerance (0.01)"),
])

COMMON_DOC = [
    ("all", CAT_ID, "library", "Kraken library id (only libraries with runs are exported, here strain_exclusion; see `libraries`)"),
    ("all", CAT_ID, "dataset", "query set: strain-exclusion | fda-argos | cami | mock | speed (see `datasets`)"),
    ("classify/classify_runs/profile/vs_reference/accuracy_*", CAT_ID, "minimum_hit_groups", "the `--minimum-hit-groups` value of the classify run: 0 = Kraken 2 paper setting (no hit-group threshold; this is the setting that must reproduce the published accuracy, see `paper_check`), 2 = the k2 default. Every (seed, dataset, order) is classified once per value. Rows are sorted dataset -> minimum_hit_groups -> order, so for one dataset and value the orders are consecutive and every `_ratio` is against xor at the same dataset and minimum_hit_groups. `build` / `build_runs` have no such column: building does not depend on it"),
    ("all", CAT_ID, "order", "minimizer scheme/order: xor, xor-pcg, greedymini, mod-minimizer-xor, mod-minimizer-greedymini (the run directory name without the GreedyMini table, which is the `gm_table` column of `orders`)"),
    ("classify_runs/build_runs", CAT_ID, "run_id", "<library>/<order>/<seed_idx>"),
    ("classify_runs/build_runs", CAT_ID, "seed_idx / seed", "seed index and its 64-bit toggle mask in hex (`k2 --toggle-mask`): XOR tie-break / random-order seed of every scheme"),
    ("build/classify/profile/vs_reference", CAT_ID, "n_seeds", "number of seeds aggregated into the row (the replicates: each seed contributes one classify run per dataset and minimum_hit_groups)"),
    ("naming", "naming", "units", "same convention as the sshash export: the raw value carries its unit as a suffix (_s seconds, _gb gigabytes, _pct percent, _ns_per_call, _per_mbp, _fraction 0-1); `classify`-type sheets add <metric>_mean, <metric>_std (sample std over seeds, rounded to >= 4 significant digits) and <metric>_ratio"),
    ("build/classify/profile", "ratio", "<metric>_ratio", "this order's mean / the xor mean of the same library (and dataset and minimum_hit_groups); reference = 1; blank when the reference mean is 0 or missing. Only for metrics with a direction: lower is better (times, RSS, index size, density factors, lookups, NextMinimizer, discordant, Bray-Curtis, off-target) or higher is better (speed, sensitivity, PPV, F1, concordant, taxa F1). Colour scale: green = better than xor, white = 1, red = worse (colours reversed for higher-is-better columns)"),
    ("build/classify/profile", "ratio", "neutral columns (no ratio, no colour)", "classified_pct (more classified is not better by itself), genus_ancestor_fraction and genus_unclassified_fraction (neither right nor wrong), *_excluded_fraction (reads without a reference: a property of library x query set), profile workload counts nextminimizer_calls_per_mbp / dup_key_hits_per_mbp and the path split uhs_pct / non_uhs_pct"),
    ("classify/classify_runs", "classify", "classify_wall_s / db_load_s / classify_processing_s", "one classify run per cell: wall = database load + processing. db_load_s = wall - processing; processing = the `processed in` seconds of the classify summary. Compare processing time between orders; db_load_s reflects the host's page cache and the index size, and is noisy for a single run (judge it over the seeds)"),
    ("build/build_runs/classify/classify_runs", "column order", "column groups", "as in the SSHash summary sheet: id / params, then density (expected, particular density factor), minimizers (distinct_minimizers), space (index_size_gb, build / build_runs only), build (build_wall_s, build_max_rss_gb; build / build_runs only), then on classify / classify_runs the classify measures and the accuracy groups. The density and minimizer columns come first because they are what the correlation rows below the data use"),
    ("build/classify", "correlation", "correlation with density / correlation with distinct minimizers (rows under the data)", "two rows right under the data, as in the SSHash summary sheet: in every `*_ratio` column a live formula IFERROR(CORREL(<reference ratio column>, <this ratio column>), \"\") with absolute ranges over all data rows (cached value stored next to it). Reference columns: particular_density_factor_ratio (the density factor measured on the benchmark library; the SSHash sheet uses its density_factor_ratio) and distinct_minimizers_ratio. Ratios are against xor within the same block, so correlating over all rows is meaningful. `build` has one row per order, so its correlations rest on n = 5 points (the n is in the row label); in `classify` the density and distinct-minimizer values repeat for every dataset / minimum_hit_groups block of the same order, so the effective number of independent points is the number of orders (5), not the number of rows. The values (and the number of pairs used) are also in kraken2_correlations.csv (columns sheet, row_label, column, value, n); the sheet CSVs themselves hold data only"),
    ("orders", "order", "cli", "k2 build flags that select the order (<orders> = the repo's orders/ dir)"),
    ("orders", "order", "scheme / sketching / rank_order / gm_table / description", "scheme family, k2 --sketching (minimizer | mod-minimizer), the --order it ranks by (xor | xor-pcg | greedymini), the GreedyMini table (gm_table, e.g. w5_k15_canonical for greedymini and w5_k6_canonical for mod-minimizer-greedymini; blank for the xor orders) and a one-line description"),
    ("orders", "density", "expected_density_factor_mean / _std", "density factor x (w+1) on random i.i.d. DNA, mean / std over the density microbenchmark rows of the order (tables/minimizer_density.csv; independent of the library; blank for orders without density rows)"),
    ("orders", "density", "<library>_particular_density_factor_mean / _std", "the same factor measured on that library's FASTAs, one column pair per library (never averaged across libraries)"),
    ("orders", "minimizer params", "l / k / w / minimizer_spaces / spaced_seed_mask", "minimizer length l, k-mer length k, window w = k - l + 1, spaced-seed spaces s (pattern 1*(l-2s) + (01)*s, so l - s kept positions) and its mask, read from the order's db/opts.k2d"),
    ("libraries", "libraries", "n_sequences / n_taxa / n_fasta_files / fasta_gb", "sequences and distinct taxids in the library's prelim_map*.txt, library.fna files found and their size on disk (GB; not re-read)"),
    ("libraries", "libraries", "role / flag / fraction / taxonomy", "what the library is for in the benchmark, the run_benchmark.sh flag that selects it, how much of NCBI it holds, and the taxonomy it uses"),
    ("libraries", "libraries", "n_orders / n_runs / query_sets_run / orders_run / status", "orders and distinct runs (order x seed) under runs/<library>, datasets and orders classified on it, whether it has runs and / or a db"),
    ("libraries", "libraries", "db_path / run_root", "shared library FASTA dir (db/<id>) and the run tree (runs/<id>)"),
    ("datasets", "dataset", "sources / n_sources / fragments_per_source", "what the reads come from, how many source genomes / taxa / experiments, fragments per source"),
    ("datasets", "dataset", "n_fragments / read_length_bp / paired / total_bases_r1_r2 / path", "read pairs, read length (bp, or min-max), paired-end flag, bases in R1+R2, location of the FASTQs"),
    ("datasets", "dataset", "truth / scoring", "truth: simulated | real. scoring: read-level truth (strain-exclusion) | concordance with SRA taxon (fda-argos) | library-aware read-level truth (cami) | profile vs expected composition (mock)"),
    ("datasets", "dataset", "truth_detail / provenance / libraries", "what the truth file holds, where the reads come from, libraries the set was classified against"),
    ("classify/classify_runs", "accuracy", "accuracy column groups", "one group per kind of accuracy, because the kinds differ per dataset and are NOT comparable across datasets (accuracy is never averaged across datasets). simulated reads with read-level truth (strain-exclusion; cami restricted to reads whose taxon is in the library): genus_/species_ sensitivity, ppv, f1. cami only: genus_/species_excluded_fraction. fda-argos (real reads, no per-read truth): genus_concordant / discordant / ancestor / unclassified _fraction. mock (real reads, known composition): genus_/species_ taxa_f1, bray_curtis, off_target_fraction (a taxa-detection F1, so named taxa_f1, not f1). A row fills only its own dataset's group; the other cells are blank, columns that no row of the sheet uses are dropped"),
    ("accuracy_*", CAT_ID, "<metric>_mean / <metric>_std / <metric>_ratio / n_seeds", "the accuracy sheets are aggregated over seeds (per-seed rows are not exposed here; the raw per-seed numbers stay in each run's results/*_metrics_mhg<m>.tsv). Rows are sorted dataset -> minimum_hit_groups -> slice -> order: the orders of one slice are consecutive, slice cells are merged and a border separates slices. Counts (n, tp, fp, ...) have a mean only. <metric>_ratio = order mean / xor mean of the same slice and minimum_hit_groups (blank when xor's mean is 0 or missing), colour green = better, red = worse"),
    ("accuracy_*", "ratio", "ratio direction and neutral columns", "higher is better (ratio coloured, reversed scale): sensitivity, ppv, f1, concordant_fraction, taxa_precision, taxa_recall, taxa_f1, on_target_fraction. Lower is better: discordant_fraction, bray_curtis, l1_distance, off_target_fraction. Neutral (mean / std only, no ratio, no colour): ancestor_fraction, unclassified_fraction, excluded_fraction (a property of library x query set), classified_fraction, above_rank_fraction and the read counts. top_off_target (text) is not aggregated"),
    ("accuracy_*", "slice", "minimum_hit_groups / rank / domain / group / scope / sra_run / species / min_frac", "the slice of the row: minimum_hit_groups (0 | 2); rank (genus | species); strain-exclusion domain (all | prokaryote | virus); cami group (all | bacteria | archaea ...); fda-argos scope (experiment = one SRR | unweighted_mean = headline | pooled) and sra_run / species; mock scope (all | bacteria) and detection threshold min_frac"),
    ("accuracy_*", "accuracy detail", "n / n_truth / n_calls / n_scored / tp / vp / fp / fn / no", "read counts of the slice. TP = true rank taxon or descendant, VP = ancestor, FP = off-lineage, FN = unclassified (Wood et al. 2019); no = reads without a call in a truth-less slice"),
    ("accuracy_*", "accuracy detail", "sensitivity / ppv / f1", "read-level: sensitivity = TP/(TP+VP+FN+FP), PPV = TP/(TP+FP), F1 their harmonic mean (strain-exclusion, cami)"),
    ("accuracy_*", "accuracy detail", "excluded / excluded_fraction / ex_unclassified / ex_consistent / ex_off", "cami: reads whose true taxon has no sequence in the library, their fraction, and what was called for them"),
    ("accuracy_*", "accuracy detail", "concordant_fraction / discordant_fraction / ancestor_fraction / unclassified_fraction", "fda-argos genus concordance of the fragments with the SRA-assigned taxon"),
    ("accuracy_*", "accuracy detail", "taxa_precision / taxa_recall / taxa_f1 / l1_distance / bray_curtis / classified_fraction / on_target_fraction / off_target_fraction / above_rank_fraction / unclassified_fraction / n_reads / n_expected / n_detected / tp_taxa / fp_taxa / fn_taxa / top_off_target", "mock profile-level metrics against the expected composition at detection threshold min_frac (taxa-level, not read-level)"),
    ("paper_check", "paper check", "table / row / domain / rank / metric", "which published number: Wood et al. 2019 supplementary table (S1 or S5), its row (classifier), the domain (S1: prokaryote | virus), the rank (genus) and the metric (sensitivity | ppv | f1 for S1; concordant | discordant | ancestor | unclassified for S5)"),
    ("paper_check", "paper check", "ours / ours_std / n_seeds / n_experiments", "our value: xor, minimum_hit_groups 0, strain_exclusion library, mean and std over the seeds (strain-exclusion domain slice of accuracy_strain_exclusion; FDA-ARGOS unweighted mean over the n_experiments experiments, Table S5 definition). n_experiments is the number of SRA experiments actually classified (the paper has 25)"),
    ("paper_check", "paper check", "published / difference / tolerance / within_tolerance", "published Kraken 2 value (tables/paper_kraken2_published.csv; `not in repo` rows have no value and no verdict), ours - published, the absolute tolerance (0.01) and `yes` / `no`. Values are reported as they are: a `no` is a real difference, not rounded away"),
    ("paper_check", "paper check", "source / note", "where the published value comes from (table, row) and caveats (e.g. fewer experiments than the paper)"),
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
    for name, cat, meaning, _dir, _src in BUILD_METRICS:
        add("build/build_runs/classify/classify_runs" if name in ORDER_PROPS else "build/build_runs",
            cat, name, meaning + (" (an order property: the same value repeats on every dataset / minimum_hit_groups row of the order in `classify` / `classify_runs`)"
                                  if name in ORDER_PROPS else ""))
    for name, cat, meaning, _dir, _src in CLASSIFY_METRICS:
        add("classify/classify_runs", cat, name, meaning)
    for name, cat, meaning, _src in RUN_ONLY:
        add("classify_runs", cat, name, meaning)
    for c in ACC_COLS:
        add("classify/classify_runs", c["band"], c["name"], c["meaning"])
    for name, cat, meaning, _dir, _src in PROFILE_METRICS:
        add("profile", cat, name, meaning)
    for name, cat, meaning, _dir, _src in VSREF_METRICS:
        add("vs_reference", cat, name, meaning)
    return rows


OVERVIEW_FIELDS = ["sheet", "category", "column", "meaning"]


# ---------------------------------------------------------------------------
# paper_check: xor at minimum_hit_groups 0 vs the published Kraken 2 numbers
# ---------------------------------------------------------------------------

PAPER_FIELDS = ["table", "row", "domain", "rank", "metric", "ours", "ours_std", "n_seeds",
                "n_experiments", "published", "difference", "tolerance", "within_tolerance",
                "source", "note"]
S1_DOMAINS = {"prokaryote": ("prokaryote",), "virus": ("virus", "viral")}


def read_published(path=None):
    """Rows of tables/paper_kraken2_published.csv (table, row, domain, rank, metric,
    value, source); value None = not in repo."""
    _f, rows = read_csv(path or PAPER_PUBLISHED)
    return rows


def paper_check_rows(base, published=None, tolerance=PAPER_TOLERANCE):
    """Rows of the `paper_check` sheet: for the reference order (xor) at
    minimum_hit_groups 0 on the strain_exclusion library, our strain-exclusion genus
    sensitivity / PPV / F1 (domain slices) and FDA-ARGOS genus concordance (unweighted
    mean over the experiments) next to the published Kraken 2 values. [] when that
    cell was not run."""
    published = published if published is not None else read_published()

    def mine(sheet):
        return [r for r in base.get(sheet) or []
                if r.get("library") == PAPER_LIBRARY and r.get("order") == PAPER_ORDER
                and r.get("minimum_hit_groups") == PAPER_MHG]

    se = [r for r in mine("accuracy_strain_exclusion") if r.get("rank") == "genus"]
    ar = [r for r in mine("accuracy_fda_argos") if r.get("rank") == "genus"]
    if not se and not ar:
        return []
    out = []
    for p in published:
        tab, domain, metric = p.get("table"), p.get("domain") or None, p.get("metric")
        ours = std = n_seeds = n_exp = None
        note = []
        if tab == "S1":
            for r in se:
                if r.get("domain") in S1_DOMAINS.get(domain, (domain,)):
                    ours, std, n_seeds = r.get(metric + "_mean"), r.get(metric + "_std"), r.get("n_seeds")
        elif tab == "S5":
            for r in ar:
                if r.get("scope") == "unweighted_mean":
                    ours = r.get(metric + "_fraction_mean")
                    std = r.get(metric + "_fraction_std")
                    n_seeds = r.get("n_seeds")
            n_exp = len(set(r.get("sra_run") for r in ar
                            if r.get("scope") == "experiment" and r.get("sra_run")))
            pub_n = p.get("n_experiments")
            if n_exp and pub_n and n_exp != pub_n:
                note.append("ours over %d experiments, published over %d" % (n_exp, pub_n))
        pub = p.get("value")
        pub = pub if _num(pub) else None
        row = OrderedDict([
            ("table", tab), ("row", p.get("row")), ("domain", domain), ("rank", p.get("rank")),
            ("metric", metric), ("ours", ours), ("ours_std", std), ("n_seeds", n_seeds),
            ("n_experiments", n_exp), ("published", pub), ("difference", None),
            ("tolerance", tolerance), ("within_tolerance", None), ("source", p.get("source")),
            ("note", None)])
        if pub is None:
            row["within_tolerance"] = "not in repo"
        elif ours is None:
            note.append("no xor / minimum_hit_groups 0 result for this value")
        else:
            row["difference"] = round(ours - pub, 6)
            row["within_tolerance"] = "yes" if abs(ours - pub) <= tolerance + 1e-12 else "no"
        row["note"] = "; ".join(note) or None
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# sheets for one set of base tables
# ---------------------------------------------------------------------------

def derive_tables(base, base_dir=None, libraries=None):
    """base: {build_runs, classify_runs, profile, vs_reference, accuracy_<dataset>...,
    orders, datasets} -> OrderedDict name -> (fields, rows) for every sheet (sheets
    without rows -- no profile data, no accuracy -- are omitted)."""
    runs = base.get("classify_runs") or []
    build_runs = base.get("build_runs") or []
    t = OrderedDict()
    t["overview"] = (OVERVIEW_FIELDS, overview_rows())
    lib_rows = libraries if libraries is not None else build_libraries(base_dir, runs) if base_dir else []
    t["libraries"] = (LIBRARY_FIELDS, lib_rows)
    t["datasets"] = (DATASET_FIELDS, base.get("datasets") or [])
    orders = base.get("orders") or []
    t["orders"] = (orders_fields(orders), orders)
    build = aggregate(build_runs, BUILD_NAMES, ("library", "order"))
    t["build"] = (prune_empty(build_field_order(), build, set(AGG_ID["build"])), build)
    t["build_runs"] = (prune_empty(build_runs_field_order(), build_runs, set(BUILD_ID_FIELDS)), build_runs)
    summ = aggregate(runs, SUMMARY_NAMES)
    t["classify"] = (prune_empty(summary_field_order(), summ, set(AGG_ID["classify"])), summ)
    t["classify_runs"] = (prune_empty(runs_field_order(), runs, set(RUN_ID_FIELDS)), runs)
    for name, order_fn in (("profile", profile_field_order), ("vs_reference", vs_reference_field_order)):
        rows = base.get(name) or []
        if rows:
            t[name] = (prune_empty(order_fn(), rows, set(AGG_ID["classify"])), rows)
    for sheet in ACCURACY_SHEETS.values():
        detail = base.get(sheet) or []
        if detail:
            t[sheet] = (detail_fields(detail), detail)
    pc = paper_check_rows(base, base.get("paper_published"))
    if pc:
        t["paper_check"] = (PAPER_FIELDS, pc)
    return t


# ---------------------------------------------------------------------------
# workbook assembly
# ---------------------------------------------------------------------------

def _categories(sheet, fields):
    if sheet in ("build", "build_runs", "classify", "classify_runs", "profile", "vs_reference"):
        return [field_category(f) for f in fields]
    if sheet == "orders":
        return [orders_category(f) for f in fields]
    if sheet == "overview":
        return ["reference"] * len(fields)
    if sheet == "paper_check":
        return ["paper check" if f not in ("table", "row", "domain", "rank", "metric") else CAT_ID
                for f in fields]
    if sheet == "libraries":
        return ["library"] * len(fields)
    if sheet == "datasets":
        return ["dataset"] * len(fields)
    # accuracy_<dataset>
    band = {"accuracy_strain_exclusion": BAND_READ, "accuracy_fda_argos": BAND_ARGOS,
            "accuracy_cami": BAND_CAMI, "accuracy_mock": BAND_MOCK}.get(sheet, "accuracy detail")
    return [CAT_ID if f in ("library", "dataset", "order", "n_seeds") else
            "slice" if f in DETAIL_SLICE else band for f in fields]


def sheet_spec(name, fields, rows):
    spec = {"name": name, "fields": list(fields), "rows": rows,
            "categories": _categories(name, fields)}
    block = ()
    merge = ()
    if name in ("classify", "classify_runs", "profile", "vs_reference"):
        # border after each (library, dataset, minimum_hit_groups) block; the dataset
        # cell spans both minimum_hit_groups blocks of the dataset
        block = tuple(f for f in ("library", "dataset", "minimum_hit_groups") if f in fields)
        merge = tuple(f for f in ("library", "dataset", "minimum_hit_groups", "order") if f in fields)
        spec["merge_scope"] = {"library": ("library",), "dataset": ("library", "dataset")}
    if name.startswith("accuracy_"):
        # border between slices; library / dataset / minimum_hit_groups cells span their whole block
        slices = tuple(f for f in DETAIL_SLICE if f in fields)
        block = ("library", "dataset") + slices
        merge = ("library", "dataset") + slices
        spec["merge_scope"] = {"library": ("library",), "dataset": ("library", "dataset"),
                               "minimum_hit_groups": ("library", "dataset", "minimum_hit_groups")}
    if name == "paper_check":
        spec["wrap_fields"] = [f for f in ("source", "note") if f in fields]
        spec["widths"] = {"source": 60, "note": 40, "row": 12}
    if name in ("build", "build_runs"):
        block = ("library",)
        merge = ("library", "order") if name == "build_runs" else ("library",)
    if name == "overview":
        block = ("sheet",)
        merge = ("sheet",)
        spec["wrap_fields"] = ["meaning"]
        spec["widths"] = {"sheet": 22, "category": 18, "column": 44, "meaning": 110}
    if name in ("libraries", "datasets", "orders"):
        spec["wrap_fields"] = [f for f in ("role", "fraction", "taxonomy", "sources", "truth_detail", "scoring",
                                           "provenance", "query_sets_run", "orders_run", "description", "cli")
                               if f in fields]
        spec["widths"] = dict((f, w) for f, w in (
            ("role", 40), ("fraction", 40), ("taxonomy", 30), ("sources", 38), ("truth_detail", 48),
            ("scoring", 26), ("truth", 14), ("provenance", 52), ("query_sets_run", 28),
            ("orders_run", 36), ("description", 56), ("cli", 52), ("path", 40), ("db_path", 40),
            ("run_root", 40)) if f in fields)
    spec["block_fields"] = block
    spec["merge_fields"] = merge
    if name in ("build", "classify", "profile"):
        spec["ratio_fields"] = dict(
            (f, METRIC_DIR[f[:-6]]) for f in fields if f.endswith("_ratio") and METRIC_DIR.get(f[:-6]))
    if name.startswith("accuracy_"):
        spec["ratio_fields"] = dict(
            (f, DETAIL_DIR[f[:-6]]) for f in fields if f.endswith("_ratio") and DETAIL_DIR.get(f[:-6]))
    if name in ("build", "build_runs", "classify", "classify_runs", "profile", "vs_reference"):
        spec["freeze_cols"] = {"build": 2, "build_runs": 2, "classify_runs": 5}.get(name, 4)
    return spec


# ---- correlation rows under the data of `build` and `classify` (as SSHash's summary sheet)
CORREL_SHEETS = ("build", "classify")
CORREL_DENSITY = "particular_density_factor_ratio"
CORREL_MINIMIZERS = "distinct_minimizers_ratio"
CORREL_ROWS = (("correlation with density", CORREL_DENSITY),
               ("correlation with distinct minimizers", CORREL_MINIMIZERS))
CORREL_FIELDS = ["sheet", "row_label", "column", "value", "n"]


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


def correlation_rows(name, fields, rows):
    """(footer rows for the sheet, rows of kraken2_correlations.csv). One footer row
    per reference column (particular density ratio, distinct minimizers ratio);
    each `*_ratio` column of the sheet gets IFERROR(CORREL(ref, column), "") with
    absolute ranges over the data rows and the python value (pairs where both
    cells are numbers) as cached value."""
    if name not in CORREL_SHEETS or not rows:
        return [], []
    first, last = DATA_START_ROW, DATA_START_ROW + len(rows) - 1
    footer, flat = [], []
    for label, ref in CORREL_ROWS:
        if ref not in fields:
            continue
        rl = col_letter(fields.index(ref) + 1)
        cells = []
        shown = label if name != "build" else "%s (n = %d)" % (label, len(rows))
        for f in fields:
            if not f.endswith("_ratio"):
                continue
            cl = col_letter(fields.index(f) + 1)
            pairs = [(r.get(ref), r.get(f)) for r in rows if _num(r.get(ref)) and _num(r.get(f))]
            val = _pearson([a for a, _b in pairs], [b for _a, b in pairs])
            val = round(val, 4) if val is not None else None
            formula = 'IFERROR(CORREL($%s$%d:$%s$%d,%s$%d:%s$%d),"")' % (rl, first, rl, last, cl, first, cl, last)
            cells.append((f, formula, val))
            flat.append(OrderedDict([("sheet", name), ("row_label", label), ("column", f),
                                     ("value", val), ("n", len(pairs))]))
        footer.append({"label": shown, "cells": cells})
    return footer, flat


def build_sheet_specs(tables):
    specs = []
    for name in SHEETS:
        if name in tables:
            spec = sheet_spec(name, tables[name][0], tables[name][1])
            footer, _flat = correlation_rows(name, tables[name][0], tables[name][1])
            if footer:
                spec["footer_rows"] = footer
            specs.append(spec)
    return specs


def write_outputs(tables, csv_dir, xlsx_path):
    """kraken2_<sheet>.csv for every sheet in csv_dir + the workbook; kraken2_*.csv
    files of sheets that no longer exist (e.g. profile when no timing files) are removed."""
    for name, (fields, rows) in tables.items():
        write_csv(os.path.join(csv_dir, PREFIX + name + ".csv"), fields, rows)
    flat = []
    for name in CORREL_SHEETS:
        if name in tables:
            flat += correlation_rows(name, tables[name][0], tables[name][1])[1]
    if flat:
        write_csv(os.path.join(csv_dir, PREFIX + "correlations.csv"), CORREL_FIELDS, flat)
    keep = set(PREFIX + n + ".csv" for n in tables) | set([PREFIX + "correlations.csv"] if flat else [])
    for f in os.listdir(csv_dir):
        if f.startswith(PREFIX) and f.endswith(".csv") and f not in keep:
            try:
                os.remove(os.path.join(csv_dir, f))
                print("removed stale table {}".format(f))
            except OSError:
                pass
    return write_book(xlsx_path, build_sheet_specs(tables))


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


def export_library(run_root, out_xlsx, library, per_seed_rows, acc_rows, dens_rows, query_rows):
    """Write <run_root>/tables/kraken2_*.csv and out_xlsx for one library.
    Returns (xlsx_path, tables_dir, tables)."""
    tables_dir = os.path.join(run_root, "tables")
    base_dir = resolve_base_dir(run_root)
    omap = order_map(per_seed_rows, acc_rows)
    raw = build_tables(library, per_seed_rows, acc_rows, dens_rows, run_root, omap)
    runs = raw["classify_runs"]
    base = {"classify_runs": runs, "build_runs": raw["build_runs"],
            "profile": aggregate(raw["profile_raw"], PROFILE_NAMES),
            "vs_reference": vs_reference_aggregate(raw["vs_reference_raw"]),
            "orders": build_orders(library, run_root, runs, dens_rows, omap),
            "datasets": decorate_datasets(query_rows, runs)}
    for sheet, drows in detail_rows(library, acc_rows, omap).items():
        base[sheet] = detail_aggregate(drows)
    cache = os.path.join(base_dir, "results", ".libraries_cache.json") if base_dir else None
    libs = build_libraries(base_dir, runs, cache) if base_dir else build_libraries_standalone(library, runs)
    libs = [r for r in libs if r["library"] == library]
    tables = derive_tables(base, libraries=libs)
    write_outputs(tables, tables_dir, out_xlsx)
    return out_xlsx, tables_dir, tables


def build_libraries_standalone(library, runs):
    lr = [r for r in runs if r["library"] == library]
    return [OrderedDict([("library", library), ("n_runs", len(lr) or None),
                         ("query_sets_run", ", ".join(sorted(set(r["dataset"] for r in lr), key=dataset_key)) or None),
                         ("orders_run", ", ".join(sorted(set(r["order"] for r in lr), key=order_key)) or None),
                         ("status", "runs")])]


def export_results(base_dir, quiet=False):
    """Merge every <base>/runs/*/tables/kraken2_classify_runs.csv into
    <base>/results/kraken2_results.xlsx + kraken2_*.csv. Returns (xlsx, csv_dir,
    n_libraries) or None when no library has been exported yet. A library whose
    tables have no `minimum_hit_groups` column (not produced by the current
    compare_benchmark.py) is skipped with a message until it is exported again."""
    runs_dir = os.path.join(base_dir, "runs")
    if not os.path.isdir(runs_dir):
        return None
    parts = defaultdict(list)
    libs_found = []
    for lib in sorted(os.listdir(runs_dir), key=library_key):
        tdir = os.path.join(runs_dir, lib, "tables")
        fields, _rows = read_csv(os.path.join(tdir, PREFIX + "classify_runs.csv"))
        if not fields:
            continue
        if "minimum_hit_groups" not in fields or "classify_wall_s" not in fields:
            print("skipping {}: tables without a minimum_hit_groups column, re-run "
                  "compare_benchmark.py on this library".format(lib), file=sys.stderr)
            continue
        libs_found.append(lib)
        for sheet in BASE_SHEETS:
            _f, rows = read_csv(os.path.join(tdir, PREFIX + sheet + ".csv"))
            if rows:
                parts[sheet].append(rows)
    if not libs_found:
        return None
    base = {}
    for sheet in ("build_runs", "classify_runs", "profile", "vs_reference") + tuple(ACCURACY_SHEETS.values()):
        rows = [r for part in parts[sheet] for r in part]
        if sheet in ("build_runs", "classify_runs"):
            rows.sort(key=_run_sort_key)
        base[sheet] = rows
    base["orders"] = merge_orders(parts["orders"])
    base["datasets"] = merge_datasets(parts["datasets"])
    results = os.path.join(base_dir, "results")
    tables = derive_tables(base, base_dir=base_dir,
                           libraries=build_libraries(base_dir, base["classify_runs"],
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
        print("no runs/*/tables/kraken2_classify_runs.csv with a minimum_hit_groups column under {} "
              "(run compare_benchmark.py on the libraries first)".format(args.base_dir), file=sys.stderr)
        return 1
    print("wrote {} ({} libraries) and kraken2_*.csv under {}".format(res[0], res[2], res[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
