#!/usr/bin/env python3
"""Export sshash benchmark runs to CSV and XLSX (stdlib only).

Writes under ``--results-dir`` (default ``$SSHASH_RESULTS_DIR`` or
``$SSHASH_BENCH_ROOT/results``):

- ``sshash_runs.csv``       -- one row per run (metrics + telemetry)
- ``sshash_summary.csv``    -- one row per dataset x order (mean/std across seeds)
- ``sshash_vs_reference.csv`` -- informational: measured vs published k-mer counts
- ``sshash_density_expected_orders.csv`` -- cached expected-density Monte Carlo
  results (independent of any run; see scripts/expected_density.py)
- ``sshash_buckets.csv``    -- bucket/binning stats (mean/std across seeds) per
  dataset x order: bucket count, average/max bucket size, skew-index share,
  and the bucket-related space breakdown (see scripts/export.py:parse_bucket_stats)
- ``sshash_results.xlsx``   -- ``index``, ``orders``, ``datasets``, ``summary``,
  ``buckets``, ``runs``, ``vs_reference``. The ``orders`` sheet carries each
  order's description/CLI plus its density -- particular (measured, per
  dataset x order) and expected (theoretical, per order) -- side by side.

Density: ``density_factor = (num_minimizer_positions / num_kmers) * (w + 1)``,
where ``w = k - m + 1`` (the window size baked into the "w21" of
``orders/GreedyMini/w21_k11_canonical``) and num_minimizer_positions is the
number of super-k-mers (minimizer-value change points) sshash build
--verbose reports. The raw ratio isn't kept as its own column --
density_factor is the reported metric. A uniformly random order's
density_factor -> 2 as w grows; GreedyMini targets a lower value by
construction.

Buckets/binning: a lower-density order has fewer distinct minimizers, hence
fewer buckets (``num_buckets``) -- but that means a larger *average* bucket
(``avg_bucket_size = num_kmers / num_buckets``) and often a bigger share of
k-mers in sshash's skew index (``kmers_in_skew_index_pct``, its structure for
large/heavy buckets), which can offset a lower density's query-time benefit
even when it still saves per-bucket space (``control_codewords``/
``mid_load_bits_per_kmer``).

Exits with an error (hard-fail gate, mirroring kmc's validate_custom_vs_lex)
if:

- ``num_kmers`` differs across any order/seed run of the same dataset, or
- any run's ``query_pos_match`` / ``query_neg_match`` fixed-set check is
  false (see scripts/equivalence.py), or
- any run's ``check_ok`` (``sshash check`` -> "EVERYTHING OK!") is false.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import statistics
import sys
import zipfile
from collections.abc import Sequence
from pathlib import Path

from scripts.common import bench_paths, die, load_dotenv
from scripts.config import seeds as bench_seeds
from scripts.density import density_factor as _density_factor
from scripts.expected_density import ensure_expected_orders
from scripts.orders import ORDERS
from xlsx_util import col_letter, root_rels_xml, write_csv as _shared_write_csv, xml_escape


# ---------------------------------------------------------------------------
# time.txt / build+query JSON parsing
# ---------------------------------------------------------------------------


def parse_time_file(path: Path) -> dict:
    """Parse GNU /usr/bin/time -v into elapsed and peak-RSS fields."""
    out: dict = {"time_elapsed_s": None, "max_rss_kb": None, "user_s": None, "sys_s": None}
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")
    for line in text.splitlines():
        s = line.strip()
        if "Elapsed (wall clock) time" in s:
            m = re.search(r"(\d+):(\d+):(\d+(?:\.\d+)?)$", s)
            if m:
                out["time_elapsed_s"] = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
            else:
                m = re.search(r"(\d+):(\d+(?:\.\d+)?)$", s)
                if m:
                    out["time_elapsed_s"] = int(m.group(1)) * 60 + float(m.group(2))
        elif "Maximum resident set size" in s:
            m = re.search(r"(\d+)\s*$", s)
            if m:
                out["max_rss_kb"] = float(m.group(1))
        elif s.startswith("User time (seconds)"):
            m = re.search(r"([0-9.]+)\s*$", s)
            out["user_s"] = float(m.group(1)) if m else None
        elif s.startswith("System time (seconds)"):
            m = re.search(r"([0-9.]+)\s*$", s)
            out["sys_s"] = float(m.group(1)) if m else None
    return out


def parse_json_blob(path: Path) -> dict:
    """sshash --verbose tools print one trailing `{...}` JSON line (stdout or stderr;
    lib/run.sh captures both combined). All values in it are JSON strings."""
    if not path.is_file():
        return {}
    text = path.read_text(errors="replace")
    m = re.search(r"\{.*\}", text, re.DOTALL)
    if not m:
        return {}
    try:
        return json.loads(m.group(0))
    except json.JSONDecodeError:
        return {}


def parse_bucket_stats(path: Path) -> dict:
    """Parse sshash build --verbose's bucket-size / skew-index diagnostics
    (src/builder/build_sparse_and_skew_index.cpp): how many buckets there
    are, how large the biggest is, and what share of buckets/k-mers land in
    the "skew index" (the structure for large/heavy buckets) -- explains why
    a lower-density order isn't automatically faster to query: fewer buckets
    can mean each one is bigger on average, and lookups in big (skew-index)
    buckets cost more."""
    out: dict = {
        "num_buckets": None,
        "buckets_in_skew_index": None,
        "max_bucket_size": None,
        "log2_max_bucket_size": None,
        "kmers_in_skew_index": None,
        "kmers_in_skew_index_pct": None,
        "skew_index_bits_per_kmer": None,
        "mid_load_bits_per_kmer": None,
        "control_codewords_bits_per_kmer": None,
    }
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")

    m = re.search(r"^num_buckets_larger_than_1_not_in_skew_index \d+/(\d+)", text, re.MULTILINE)
    if m:
        out["num_buckets"] = int(m.group(1))
    m = re.search(r"^num_buckets_in_skew_index (\d+)/\d+", text, re.MULTILINE)
    if m:
        out["buckets_in_skew_index"] = int(m.group(1))
    m = re.search(r"^max_bucket_size (\d+)", text, re.MULTILINE)
    if m:
        out["max_bucket_size"] = int(m.group(1))
    m = re.search(r"^log2_max_bucket_size (\d+)", text, re.MULTILINE)
    if m:
        out["log2_max_bucket_size"] = int(m.group(1))
    m = re.search(r"^num kmers in skew index = (\d+) \(([\d.]+)%\)", text, re.MULTILINE)
    if m:
        out["kmers_in_skew_index"] = int(m.group(1))
        out["kmers_in_skew_index_pct"] = float(m.group(2))
    m = re.search(r"skew_index: ([\d.eE+-]+) \[bits/kmer\]", text)
    if m:
        out["skew_index_bits_per_kmer"] = float(m.group(1))
    m = re.search(r"mid_load_buckets: ([\d.eE+-]+) \[bits/kmer\]", text)
    if m:
        out["mid_load_bits_per_kmer"] = float(m.group(1))
    m = re.search(r"control_codewords: ([\d.eE+-]+) \[bits/kmer\]", text)
    if m:
        out["control_codewords_bits_per_kmer"] = float(m.group(1))
    return out


def parse_minimizer_counts(path: Path) -> dict:
    """Parse the plain-text 'num_minimizers = X' / 'num_minimizer_positions =
    Y' lines `sshash build --verbose` prints (dictionary_builder.hpp) --
    num_minimizer_positions is the number of super-k-mers (minimizer-value
    change points), the input to "particular" (measured) density."""
    out: dict = {"num_minimizers": None, "num_minimizer_positions": None}
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")
    m = re.search(r"^num_minimizers = (\d+)", text, re.MULTILINE)
    if m:
        out["num_minimizers"] = int(m.group(1))
    m = re.search(r"^num_minimizer_positions = (\d+)", text, re.MULTILINE)
    if m:
        out["num_minimizer_positions"] = int(m.group(1))
    return out


def _num(v, cast=float):
    if v is None or v == "":
        return None
    try:
        return cast(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Per-run row loading
# ---------------------------------------------------------------------------

RUN_FIELDS = (
    # id / params
    "run_id",
    "dataset",
    "k",
    "m",
    "w",
    "order",
    "order_seed",
    "seed_idx",
    "gm_file",
    "threads",
    # space
    "num_kmers",
    "index_size_in_bytes",
    "index_bytes_on_disk",
    "bits_per_kmer",
    # density (raw minimizer_density/num_minimizers/num_minimizer_positions
    # dropped -- density_factor derives from them and is enough on its own)
    "density_factor",
    # buckets
    "num_buckets",
    "avg_bucket_size",
    "max_bucket_size",
    "log2_max_bucket_size",
    "kmers_in_skew_index_pct",
    "skew_index_bits_per_kmer",
    "mid_load_bits_per_kmer",
    "control_codewords_bits_per_kmer",
    # build
    "build_wall_s",
    "build_max_rss_kb",
    "build_total_build_time_s",
    "build_step1_s",
    "build_step2_s",
    "build_step3_s",
    "build_step4_s",
    "build_step5_s",
    "build_step6_s",
    "build_step7_s",
    # query (real FASTQ) -- kept next to build, the other real wall-clock op
    "query_wall_s",
    "query_max_rss_kb",
    "query_elapsed_ms",
    "query_num_kmers",
    "query_num_positive_kmers",
    "query_num_negative_kmers",
    "query_pos_match",
    "query_neg_match",
    # bench (sshash's own internal micro-benchmark, see "what is bench" note
    # on _METRIC_CATEGORY below)
    "bench_positive_ns",
    "bench_negative_ns",
    "bench_access_ns",
    "bench_iterator_ns",
    # correctness
    "check_ok",
    # telemetry
    "sshash_version",
    "host",
    "build_exit_code",
    "check_exit_code",
    "bench_exit_code",
    "query_exit_code",
    "python_wall_s",
)


def iter_run_dirs(runs_root: Path) -> list[Path]:
    if not runs_root.is_dir():
        return []
    found = [p.parent for p in runs_root.rglob("meta.json") if (p.parent / "build.log").is_file()]
    return sorted(set(found))


def _bool_env(v) -> bool | None:
    if v in (True, False):
        return v
    if v in ("true", "True", "1"):
        return True
    if v in ("false", "False", "0"):
        return False
    return None


def load_run_row(run_dir: Path, runs_root: Path) -> dict | None:
    meta_path = run_dir / "meta.json"
    if not meta_path.is_file():
        return None
    meta = json.loads(meta_path.read_text())
    build_json = parse_json_blob(run_dir / "build.log")
    bench_json = parse_json_blob(run_dir / "bench.log")
    query_json = parse_json_blob(run_dir / "query.log")
    equiv_json = {}
    equiv_path = run_dir / "equivalence.json"
    if equiv_path.is_file():
        try:
            equiv_json = json.loads(equiv_path.read_text())
        except json.JSONDecodeError:
            equiv_json = {}
    build_time = parse_time_file(run_dir / "build_time.txt")
    query_time = parse_time_file(run_dir / "query_time.txt")

    num_kmers = _num(build_json.get("num_kmers"), int)
    index_bytes_json = _num(build_json.get("index_size_in_bytes"), int)
    index_bytes_disk = meta.get("index_bytes_on_disk")
    bits_per_kmer = None
    if index_bytes_disk and num_kmers:
        bits_per_kmer = round(index_bytes_disk * 8 / num_kmers, 4)

    mini = parse_minimizer_counts(run_dir / "build.log")
    buckets = parse_bucket_stats(run_dir / "build.log")
    avg_bucket_size = None
    if buckets["num_buckets"] and num_kmers:
        avg_bucket_size = round(num_kmers / buckets["num_buckets"], 2)
    k_val, m_val = meta.get("k"), meta.get("m")
    w = int(k_val) - int(m_val) + 1 if k_val is not None and m_val is not None else None
    # Raw minimizer_density (= num_minimizer_positions / num_kmers) isn't
    # kept on the row -- density_factor (that times (w+1)) is enough.
    density_fac = None
    if mini["num_minimizer_positions"] and num_kmers and w is not None:
        density_fac = _density_factor(mini["num_minimizer_positions"], num_kmers, w)

    def step(key: str) -> float | None:
        v = _num(build_json.get(key), float)
        return round(v / 1_000_000, 6) if v is not None else None

    row = {
        "run_id": meta.get("run_id") or str(run_dir.relative_to(runs_root)),
        "dataset": meta.get("dataset"),
        "k": meta.get("k"),
        "m": meta.get("m"),
        "order": meta.get("order"),
        "order_seed": meta.get("order_seed"),
        "seed_idx": meta.get("seed_idx"),
        "gm_file": meta.get("gm_file") or "",
        "threads": meta.get("threads"),
        "num_kmers": num_kmers,
        "index_size_in_bytes": index_bytes_json,
        "index_bytes_on_disk": index_bytes_disk,
        "bits_per_kmer": bits_per_kmer,
        "w": w,
        "density_factor": density_fac,
        "num_buckets": buckets["num_buckets"],
        "avg_bucket_size": avg_bucket_size,
        "max_bucket_size": buckets["max_bucket_size"],
        "log2_max_bucket_size": buckets["log2_max_bucket_size"],
        "kmers_in_skew_index_pct": buckets["kmers_in_skew_index_pct"],
        "skew_index_bits_per_kmer": buckets["skew_index_bits_per_kmer"],
        "mid_load_bits_per_kmer": buckets["mid_load_bits_per_kmer"],
        "control_codewords_bits_per_kmer": buckets["control_codewords_bits_per_kmer"],
        "build_wall_s": build_time["time_elapsed_s"],
        "build_max_rss_kb": build_time["max_rss_kb"],
        "build_total_build_time_s": (
            round(_num(build_json.get("total_build_time_in_microsec"), float) / 1_000_000, 6)
            if build_json.get("total_build_time_in_microsec") is not None
            else None
        ),
        "build_step1_s": step("step 1 (encode strings)"),
        "build_step2_s": step("step 2 (compute minimizer tuples)"),
        "build_step3_s": step("step 3 (merging minimizer tuples)"),
        "build_step4_s": step("step 4 (build mphf)"),
        "build_step5_s": step("step 5 (replacing minimizer values with MPHF hashes)"),
        "build_step6_s": step("step 6 (merging minimizers tuples)"),
        "build_step7_s": step("step 7 (build sparse and skew index)"),
        "check_ok": meta.get("check_ok"),
        "bench_positive_ns": _num(bench_json.get("positive lookup (avg_nanosec_per_kmer)")),
        "bench_negative_ns": _num(bench_json.get("negative lookup (avg_nanosec_per_kmer)")),
        "bench_access_ns": _num(bench_json.get("access (avg_nanosec_per_kmer)")),
        "bench_iterator_ns": _num(bench_json.get("iterator (avg_nanosec_per_kmer)")),
        "query_wall_s": query_time["time_elapsed_s"],
        "query_max_rss_kb": query_time["max_rss_kb"],
        "query_elapsed_ms": _num(query_json.get("elapsed_millisec")),
        "query_num_kmers": _num(query_json.get("num_kmers"), int),
        "query_num_positive_kmers": _num(query_json.get("num_positive_kmers"), int),
        "query_num_negative_kmers": _num(query_json.get("num_negative_kmers"), int),
        "query_pos_match": _bool_env(
            meta.get("query_pos_match") or (equiv_json.get("pos") or {}).get("match")
        ),
        "query_neg_match": _bool_env(
            meta.get("query_neg_match") or (equiv_json.get("neg") or {}).get("match")
        ),
        "sshash_version": meta.get("sshash_version"),
        "host": meta.get("host"),
        "build_exit_code": meta.get("build_exit_code"),
        "check_exit_code": meta.get("check_exit_code"),
        "bench_exit_code": meta.get("bench_exit_code"),
        "query_exit_code": meta.get("query_exit_code"),
        "python_wall_s": meta.get("python_wall_s"),
    }
    return row


# ---------------------------------------------------------------------------
# Equivalence hard-fail gate
# ---------------------------------------------------------------------------


def validate_orders_equivalence(rows: list[dict]) -> tuple[int, list[str]]:
    """Every order/seed run of the same dataset must report the same
    ``num_kmers`` and pass its check/query-answer gates. Returns
    (n_comparisons, errors); callers die() hard on any error (see
    task spec: this is the sshash analogue of kmc's validate_custom_vs_lex,
    generalized since sshash has no seed-less baseline -- every order,
    including native, is just another group member here)."""
    errors: list[str] = []
    n_cmp = 0
    by_dataset: dict[str, list[dict]] = {}
    for r in rows:
        by_dataset.setdefault(r.get("dataset") or "", []).append(r)

    for dataset, group in sorted(by_dataset.items()):
        counts = {r.get("num_kmers") for r in group if r.get("num_kmers") is not None}
        n_cmp += len(group)
        if len(counts) > 1:
            detail = ", ".join(
                f"{r.get('order')}/seed{r.get('seed_idx')}={r.get('num_kmers')}" for r in group
            )
            errors.append(f"{dataset}: num_kmers differs across orders/seeds: {detail}")
        for r in group:
            if r.get("check_ok") is False:
                errors.append(f"{dataset} {r.get('run_id')}: sshash check did not report EVERYTHING OK!")
            if r.get("query_pos_match") is False:
                errors.append(f"{dataset} {r.get('run_id')}: query_pos_match is false")
            if r.get("query_neg_match") is False:
                errors.append(f"{dataset} {r.get('run_id')}: query_neg_match is false")
    return n_cmp, errors


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

SUMMARY_METRICS = (
    # space
    "bits_per_kmer",
    # density (raw minimizer_density dropped -- density_factor is enough)
    "density_factor",
    # buckets
    "num_buckets",
    "avg_bucket_size",
    "max_bucket_size",
    "kmers_in_skew_index_pct",
    "skew_index_bits_per_kmer",
    "mid_load_bits_per_kmer",
    "control_codewords_bits_per_kmer",
    # build
    "build_wall_s",
    "build_max_rss_kb",
    # query (real FASTQ) -- kept next to build
    "query_wall_s",
    "query_elapsed_ms",
    # bench (sshash's own internal micro-benchmark)
    "bench_positive_ns",
    "bench_negative_ns",
    "bench_access_ns",
)


def _mean_std(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return statistics.mean(values), statistics.stdev(values)


def summary_group_key(r: dict) -> tuple:
    return (r.get("dataset") or "", int(r.get("k") or 0), int(r.get("m") or 0), r.get("order") or "")


def summary_sort_key(row: dict) -> tuple:
    try:
        rank = ORDERS.index(row.get("order") or "")
    except ValueError:
        rank = len(ORDERS)
    return (row.get("dataset") or "", rank)


def summary_rows(rows: list[dict]) -> list[dict]:
    buckets: dict[tuple, list[dict]] = {}
    for r in rows:
        buckets.setdefault(summary_group_key(r), []).append(r)

    out: list[dict] = []
    for key, parts in buckets.items():
        head = parts[0]
        row: dict = {
            "dataset": head["dataset"],
            "k": head["k"],
            "m": head["m"],
            "w": head.get("w"),
            "order": head["order"],
            "n_seeds": len({p.get("seed_idx") for p in parts}),
            "n_kmers": head.get("num_kmers"),
        }
        for metric in SUMMARY_METRICS:
            vals = [float(p[metric]) for p in parts if p.get(metric) is not None]
            mean, std = _mean_std(vals)
            row[f"{metric}_mean"] = round(mean, 4) if mean is not None else None
            row[f"{metric}_std"] = round(std, 4) if std is not None else None
        out.append(row)
    return sorted(out, key=summary_sort_key)


def attach_native_ratios(rows: list[dict]) -> None:
    """`{metric}_ratio` = this order's mean / native's mean, same dataset (native = 1)."""
    native_by_dataset = {r["dataset"]: r for r in rows if r.get("order") == "native"}
    for r in rows:
        base = native_by_dataset.get(r["dataset"])
        for metric in SUMMARY_METRICS:
            mean = r.get(f"{metric}_mean")
            base_mean = base.get(f"{metric}_mean") if base else None
            rname = f"{metric}_ratio"
            if mean is None or not base_mean:
                r[rname] = None
            else:
                r[rname] = round(mean / base_mean, 4)


REFERENCE_KMERS = {
    # Approximate distinct-31-mer counts from the Zenodo 17582116 record's own
    # published table (informational only -- not order-specific, just a sanity
    # cross-check that our eulertigs-derived num_kmers is in the right ballpark).
    "cod": 502_500_000,
    "kestrel": 1_150_000_000,
    "human": 2_510_000_000,
}


def vs_reference_rows(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    by_dataset: dict[str, dict] = {}
    for r in rows:
        by_dataset.setdefault(r["dataset"], r)  # first row's num_kmers stands for the dataset
    for dataset, ref_kmers in REFERENCE_KMERS.items():
        r = by_dataset.get(dataset)
        measured = r.get("num_kmers") if r else None
        ratio = round(measured / ref_kmers, 4) if measured and ref_kmers else None
        out.append(
            {
                "dataset": dataset,
                "published_num_kmers_approx": ref_kmers,
                "measured_num_kmers": measured,
                "measured_vs_published_ratio": ratio,
            }
        )
    return out


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------


def write_csv(path: Path, rows: list[dict], fields: Sequence[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = list(fields) if fields else list(rows[0].keys())
    _shared_write_csv(path, fieldnames, rows)
    print(f"wrote {path} ({len(rows)} rows)")


def summary_field_order() -> list[str]:
    fields = ["dataset", "k", "m", "w", "order", "n_seeds", "n_kmers"]
    for metric in SUMMARY_METRICS:
        fields.extend((f"{metric}_mean", f"{metric}_std", f"{metric}_ratio"))
    return fields


def vs_reference_field_order() -> list[str]:
    return ["dataset", "published_num_kmers_approx", "measured_num_kmers", "measured_vs_published_ratio"]


# ---------------------------------------------------------------------------
# XLSX (hand-rolled OOXML, stdlib only -- see benchmark/common/xlsx_util.py)
# ---------------------------------------------------------------------------

# Row 1 = category band (merged across contiguous same-category fields), row 2
# = column headers, data from row 3. Mirrors benchmark/kmc/scripts/export.py's
# two-row-header + vertical-merge + dataset-border approach (simplified: no
# banner row, no per-cell colored fills, no z-score dxfs -- sshash doesn't
# need those).
XLSX_HEADER_ROW = 2
XLSX_DATA_START_ROW = 3

# Repeated across every row of a dataset's block (k/m/w are fixed for the
# whole benchmark; dataset obviously repeats) -- merge vertically instead of
# printing the same value on every row.
XLSX_MERGE_FIELDS = frozenset({"dataset", "k", "m", "w"})
# Sheets with one row per (dataset, ...): get a divider (bottom border) after
# each dataset's last row, and the vertical merge above.
XLSX_DATASET_BORDER_SHEETS = frozenset({"runs", "summary", "buckets"})

_ID_FIELDS = frozenset({
    "run_id", "dataset", "k", "m", "w", "order", "order_seed", "seed_idx",
    "gm_file", "threads", "n_seeds", "n_kmers",
})
_METRIC_CATEGORY = {
    # space
    "num_kmers": "space", "index_size_in_bytes": "space", "index_bytes_on_disk": "space", "bits_per_kmer": "space",
    # density
    "density_factor": "density",
    # buckets
    "num_buckets": "buckets", "avg_bucket_size": "buckets", "max_bucket_size": "buckets",
    "log2_max_bucket_size": "buckets", "kmers_in_skew_index_pct": "buckets",
    "skew_index_bits_per_kmer": "buckets", "mid_load_bits_per_kmer": "buckets",
    "control_codewords_bits_per_kmer": "buckets",
    # build
    "build_wall_s": "build", "build_max_rss_kb": "build", "build_total_build_time_s": "build",
    "build_step1_s": "build", "build_step2_s": "build", "build_step3_s": "build", "build_step4_s": "build",
    "build_step5_s": "build", "build_step6_s": "build", "build_step7_s": "build",
    # correctness
    "check_ok": "correctness",
    # bench (sshash's own internal micro-benchmark)
    "bench_positive_ns": "bench", "bench_negative_ns": "bench", "bench_access_ns": "bench", "bench_iterator_ns": "bench",
    # query (real FASTQ)
    "query_wall_s": "query", "query_max_rss_kb": "query", "query_elapsed_ms": "query",
    "query_num_kmers": "query", "query_num_positive_kmers": "query", "query_num_negative_kmers": "query",
    "query_pos_match": "query", "query_neg_match": "query",
    # telemetry
    "sshash_version": "telemetry", "host": "telemetry", "build_exit_code": "telemetry",
    "check_exit_code": "telemetry", "bench_exit_code": "telemetry", "query_exit_code": "telemetry",
    "python_wall_s": "telemetry",
}


def _field_category(field: str) -> str:
    if field in _ID_FIELDS:
        return "id / params"
    base = field
    for suffix in ("_mean", "_std", "_ratio"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return _METRIC_CATEGORY.get(base, "")


def _orders_field_category(field: str) -> str:
    if field in ("order", "description", "cli", "w"):
        return "order"
    if field.startswith("expected_density_factor"):
        return "expected (theoretical)"
    marker = "_particular_density_factor_"
    if marker in field:
        return field.split(marker, 1)[0]
    return ""


def sheet_categories(raw_name: str, fields: Sequence[str]) -> list[str] | None:
    """One category label per field, for the row-1 band; None to skip it."""
    if raw_name in ("runs", "summary", "buckets"):
        return [_field_category(f) for f in fields]
    if raw_name == "orders":
        return [_orders_field_category(f) for f in fields]
    if raw_name == "index":
        return ["reference"] * len(fields)
    if raw_name == "datasets":
        return ["dataset reference"] * len(fields)
    if raw_name == "vs_reference":
        return ["reference"] * len(fields)
    return None


def _dataset_last_row_indices(rows: Sequence[dict]) -> set[int]:
    """0-indexed (relative to data start) row indices where a dataset's block
    ends -- where the divider border goes."""
    out: set[int] = set()
    for i, row in enumerate(rows):
        if i == len(rows) - 1 or row.get("dataset") != rows[i + 1].get("dataset"):
            out.add(i)
    return out


def _merge_rows_equal(rows: Sequence[dict], i: int, j: int, field: str) -> bool:
    if rows[i].get("dataset") != rows[j].get("dataset"):
        return False
    return rows[i].get(field) == rows[j].get(field)


def _vertical_merge_plan(
    rows: Sequence[dict], fields: Sequence[str], merge_fields: frozenset[str]
) -> tuple[list[str], set[tuple[int, str]], set[tuple[int, str]]]:
    """Return (mergeCell refs, skipped (row_idx, field) cells, merge-anchor
    (row_idx, field) cells). Runs never cross a dataset boundary even if the
    value happens to repeat (k/m/w are constant for the whole benchmark)."""
    refs: list[str] = []
    skip: set[tuple[int, str]] = set()
    anchors: set[tuple[int, str]] = set()
    if len(rows) < 2:
        return refs, skip, anchors
    col = {name: i for i, name in enumerate(fields)}

    def close_run(field: str, run_start: int, run_end: int) -> None:
        if run_end <= run_start:
            return
        letter = col_letter(col[field] + 1)
        r1 = XLSX_DATA_START_ROW + run_start
        r2 = XLSX_DATA_START_ROW + run_end
        refs.append(f"{letter}{r1}:{letter}{r2}")
        anchors.add((run_start, field))
        for r in range(run_start + 1, run_end + 1):
            skip.add((r, field))

    for field in fields:
        if field not in merge_fields or field not in col:
            continue
        run_start = 0
        for i in range(1, len(rows)):
            if _merge_rows_equal(rows, i - 1, i, field):
                continue
            close_run(field, run_start, i - 1)
            run_start = i
        close_run(field, run_start, len(rows) - 1)
    return refs, skip, anchors


def _xlsx_merges_xml(refs: Sequence[str]) -> str:
    if not refs:
        return ""
    cells = "".join('<mergeCell ref="{}"/>'.format(r) for r in refs)
    return f'<mergeCells count="{len(refs)}">{cells}</mergeCells>'


def _xlsx_cell(value, col: int, row: int, *, style: int = 0) -> str:
    ref = f"{col_letter(col + 1)}{row}"
    s_attr = f' s="{style}"' if style else ""
    if value is None or value == "":
        return f'<c r="{ref}"{s_attr}/>'
    if isinstance(value, bool):
        return f'<c r="{ref}"{s_attr} t="b"><v>{1 if value else 0}</v></c>'
    if isinstance(value, int):
        return f'<c r="{ref}"{s_attr}><v>{value}</v></c>'
    if isinstance(value, float):
        return f'<c r="{ref}"{s_attr}><v>{value}</v></c>'
    text = xml_escape(str(value))
    return f'<c r="{ref}"{s_attr} t="inlineStr"><is><t xml:space="preserve">{text}</t></is></c>'


def _ratio_bounds(rows: Sequence[dict], field: str) -> tuple[float, float]:
    vals = [v for r in rows if (v := r.get(field)) is not None]
    lo = min([1.0, *vals])
    hi = max([1.0, *vals])
    if lo == 1.0:
        lo = 0.99
    if hi == 1.0:
        hi = 1.01
    return lo, hi


def _color_scale_xml(fields: Sequence[str], rows: Sequence[dict]) -> str:
    """Green (< 1, better than native) -> white (1) -> red (> 1, worse) scale
    on every `{metric}_ratio` column -- not a curated subset."""
    if not rows:
        return ""
    data_start = XLSX_DATA_START_ROW
    last_row = data_start + len(rows) - 1
    parts: list[str] = []
    for priority, field in enumerate(fields, start=1):
        if not field.endswith("_ratio"):
            continue
        col_idx = fields.index(field)
        letter = col_letter(col_idx + 1)
        lo, hi = _ratio_bounds(rows, field)
        sqref = f"{letter}{data_start}:{letter}{last_row}"
        parts.append(
            f'<conditionalFormatting sqref="{sqref}">'
            f'<cfRule type="colorScale" priority="{priority}"><colorScale>'
            f'<cfvo type="num" val="{lo}"/><cfvo type="num" val="1"/><cfvo type="num" val="{hi}"/>'
            '<color rgb="FF63BE7B"/><color rgb="FFFFFFFF"/><color rgb="FFF8696B"/>'
            "</colorScale></cfRule></conditionalFormatting>"
        )
    return "".join(parts)


def _category_band_row(fields: Sequence[str], categories: Sequence[str]) -> tuple[str, list[str]]:
    """Row 1: category labels, merged across contiguous fields sharing one."""
    if not fields:
        return "", []
    cells: list[str] = []
    merges: list[str] = []
    i = 0
    while i < len(fields):
        cat = categories[i]
        start = i
        while i < len(fields) and categories[i] == cat:
            i += 1
        end = i - 1
        cells.append(_xlsx_cell(cat, start, 1, style=2))
        if end > start:
            merges.append(f"{col_letter(start + 1)}1:{col_letter(end + 1)}1")
    return f'<row r="1">{"".join(cells)}</row>', merges


def _sheet_xml(raw_name: str, fields: list[str], rows: list[dict]) -> str:
    xml = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    ]
    if fields:
        xml.append(
            "<sheetViews><sheetView workbookViewId=\"0\">"
            f'<pane ySplit="{XLSX_HEADER_ROW}" topLeftCell="A{XLSX_DATA_START_ROW}" '
            'activePane="bottomLeft" state="frozen"/>'
            "</sheetView></sheetViews>"
        )
    xml.append("<sheetData>")
    merge_refs: list[str] = []
    if fields:
        categories = sheet_categories(raw_name, fields)
        if categories:
            band_xml, band_merges = _category_band_row(fields, categories)
            xml.append(band_xml)
            merge_refs.extend(band_merges)

        header_cells = "".join(_xlsx_cell(h, i, XLSX_HEADER_ROW, style=2) for i, h in enumerate(fields))
        xml.append(f'<row r="{XLSX_HEADER_ROW}">{header_cells}</row>')

        merge_fields = XLSX_MERGE_FIELDS if raw_name in XLSX_DATASET_BORDER_SHEETS else frozenset()
        vertical_refs, merge_skip, merge_anchors = _vertical_merge_plan(rows, fields, merge_fields)
        merge_refs.extend(vertical_refs)
        dataset_last = _dataset_last_row_indices(rows) if raw_name in XLSX_DATASET_BORDER_SHEETS else set()

        for r_i, row in enumerate(rows, start=XLSX_DATA_START_ROW):
            row_idx = r_i - XLSX_DATA_START_ROW
            border = row_idx in dataset_last
            style = 1 if border else 0
            cells: list[str] = []
            for c, h in enumerate(fields):
                skipped = (row_idx, h) in merge_skip
                # A cell merged away (not the run's anchor) still needs to be
                # emitted -- blank, just for its style -- if this row is also
                # a dataset-divider row, so the border closes off visually
                # instead of vanishing under the merge.
                if skipped and not border:
                    continue
                value = None if skipped else row.get(h)
                cells.append(_xlsx_cell(value, c, r_i, style=style))
            xml.append(f'<row r="{r_i}">{"".join(cells)}</row>')
    xml.append("</sheetData>")
    if merge_refs:
        xml.append(_xlsx_merges_xml(merge_refs))
    if fields and rows:
        cs = _color_scale_xml(fields, rows)
        if cs:
            xml.append(cs)
    xml.append("</worksheet>")
    return "\n".join(xml)


def _xlsx_styles_xml() -> str:
    """Minimal styles: default, a bottom-border style (dataset divider), and
    a bold style (category band + column headers)."""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2">'
        '<font><sz val="11"/><color theme="1"/><name val="Calibri"/></font>'
        '<font><b/><sz val="11"/><color theme="1"/><name val="Calibri"/></font>'
        "</fonts>"
        '<fills count="2">'
        '<fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        "</fills>"
        '<borders count="2">'
        "<border><left/><right/><top/><bottom/><diagonal/></border>"
        '<border><left/><right/><top/><bottom style="medium"><color auto="1"/></bottom><diagonal/></border>'
        "</borders>"
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="3">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1"/>'
        '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
        "</cellXfs>"
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        "</styleSheet>"
    )


def write_xlsx(path: Path, sheets: Sequence[tuple[str, list[str], list[dict]]]) -> None:
    """Write a multi-sheet .xlsx without extra Python packages."""
    names: list[str] = []
    payloads: list[str] = []
    for raw_name, fields, rows in sheets:
        name = (raw_name or "Sheet")[:31]
        names.append(name)
        payloads.append(_sheet_xml(raw_name, fields, rows))

    sheet_refs, rels, overrides = [], [], [
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]
    rels_head = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
        '<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>',
    ]
    for i, name in enumerate(names, start=1):
        rid = f"rIdSheet{i}"
        sheet_refs.append(f'<sheet name="{xml_escape(name)}" sheetId="{i}" r:id="{rid}"/>')
        rels.append(
            f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{i}.xml"/>'
        )
        overrides.append(
            f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    rels_xml = "\n".join([*rels_head, *rels, "</Relationships>"])
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets>{"".join(sheet_refs)}</sheets></workbook>'
    )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        f'{"".join(overrides)}</Types>'
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", root_rels_xml())
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels", rels_xml)
        z.writestr("xl/styles.xml", _xlsx_styles_xml())
        for i, body in enumerate(payloads, start=1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", body)
    print(f"wrote {path} ({len(names)} sheets)")


def index_rows() -> list[dict]:
    def row(sheet: str, column: str, meaning: str) -> dict:
        return {"sheet": sheet, "column": column, "meaning": meaning}

    return [
        row("all", "dataset / k / m", "dataset id, k-mer length, minimizer length (m fixed at 11 for every order -- see README 'Choice of m')"),
        row("all", "order", "native | xor-hash | xor-pcg-hash | greedymini -- see orders/README.md"),
        row("all", "num_kmers", "distinct k-mers in the built index; must be identical across every order/seed for a dataset (equivalence gate)"),
        row("summary", "n_seeds", "number of seed repeats aggregated into this row's mean/std"),
        row("summary/runs", "bits_per_kmer", "index_bytes_on_disk * 8 / num_kmers"),
        row("summary/runs", "bench_*_ns", "sshash bench's internal avg_nanosec_per_kmer (the paper's own headline query-time metric)"),
        row("summary/runs", "query_*", "secondary: sshash query against the dataset's real paired FASTQ (wall time via /usr/bin/time -v, elapsed_millisec from sshash's own JSON)"),
        row("summary", "{metric}_ratio", "this order's mean / native's mean for the same dataset (native = 1); color scale green (better) -> red (worse)"),
        row("runs", "query_pos_match / query_neg_match", "fixed positive/negative k-mer set (generated once per dataset from a native reference build) answered identically by this order/seed's index -- see scripts/equivalence.py; export hard-fails if any is false"),
        row("vs_reference", "published_num_kmers_approx", "Zenodo 17582116's own published distinct-31-mer count (informational only, not order-specific)"),
        row("summary/runs/orders", "w", "window size k - m + 1 (=21 for k=31,m=11; the 'w21' in orders/GreedyMini/w21_k11_canonical)"),
        row("summary/runs/orders", "density_factor", "particular (measured): (num_minimizer_positions/num_kmers) * (w+1) -- see scripts/density.py"),
        row("orders", "expected_density_factor_mean/std", "theoretical (Monte Carlo on i.i.d. random DNA, mean/std over 5 seeds), independent of dataset -- see scripts/expected_density.py"),
        row("buckets/summary/runs", "num_buckets / avg_bucket_size", "num_kmers / num_buckets; a lower-density order has fewer buckets, hence a larger average bucket -- see scripts/export.py:parse_bucket_stats"),
        row("buckets/summary/runs", "kmers_in_skew_index_pct", "share of k-mers whose bucket is large enough to go in sshash's skew index (the structure for heavy buckets); higher can mean slower lookups despite a smaller/denser index"),
        row("buckets/summary/runs", "skew_index_bits_per_kmer / mid_load_bits_per_kmer / control_codewords_bits_per_kmer", "space breakdown for the bucket-count/size tradeoff: fewer buckets means less control_codewords/mid_load overhead but a bigger skew_index share"),
    ]


ORDER_DESCRIPTIONS = {
    "native": {
        "description": "sshash's own unmodified baseline order (mixer_64-based)",
        "cli": "(no -O)",
    },
    "xor-hash": {
        "description": "XOR the minimizer with the seed, then rank canonically",
        "cli": "-O xor-hash -s SEED",
    },
    "xor-pcg-hash": {
        "description": "xor-hash's XOR step, then PCG RXS-M-XS 64-bit scramble before ranking",
        "cli": "-O xor-pcg-hash -s SEED",
    },
    "greedymini": {
        "description": "Ranks minimizers by the precomputed GreedyMini w21_k11 rank table (gm_k=11, matches -m 11 exactly)",
        "cli": "-O greedymini --greedy-mini-file orders/GreedyMini/w21_k11_canonical.gm -s SEED",
    },
}


DATASET_ORDER_PREFERENCE = ("cod", "kestrel", "human")


def _dataset_sort_key(dataset: str) -> tuple:
    try:
        return (0, DATASET_ORDER_PREFERENCE.index(dataset))
    except ValueError:
        return (1, dataset)


def summary_datasets(summary: list[dict]) -> list[str]:
    return sorted({r["dataset"] for r in summary if r.get("dataset")}, key=_dataset_sort_key)


def orders_rows(summary: list[dict], expected: dict[str, dict] | None = None) -> list[dict]:
    """One row per order (not per dataset x order -- w, expected density, and
    description/CLI don't vary by dataset, so repeating them per dataset is
    pure duplication): the order's description/CLI/w/expected density_factor,
    plus each dataset's particular (measured, mean/std over seeds)
    density_factor pivoted into its own columns (see orders_field_order).
    Raw minimizer_density isn't kept anywhere -- density_factor is enough."""
    expected = expected or {}
    datasets = summary_datasets(summary)
    by_order: dict[str, dict[str, dict]] = {}
    for r in summary:
        order = r.get("order")
        if order:
            by_order.setdefault(order, {})[r.get("dataset")] = r

    out = []
    for order in ORDERS:
        per_dataset = by_order.get(order)
        if not per_dataset:
            continue
        info = ORDER_DESCRIPTIONS.get(order, {"description": "", "cli": ""})
        exp = expected.get(order, {})
        w_val = next((r.get("w") for r in per_dataset.values() if r.get("w") is not None), None)
        row = {
            "order": order,
            "description": info["description"],
            "cli": info["cli"],
            "w": w_val,
            "expected_density_factor_mean": _num(exp.get("density_factor_avg")),
            "expected_density_factor_std": _num(exp.get("density_factor_std")),
        }
        for ds in datasets:
            r = per_dataset.get(ds)
            row[f"{ds}_particular_density_factor_mean"] = r.get("density_factor_mean") if r else None
            row[f"{ds}_particular_density_factor_std"] = r.get("density_factor_std") if r else None
        out.append(row)
    return out


def orders_field_order(summary: list[dict]) -> list[str]:
    fields = ["order", "description", "cli", "w", "expected_density_factor_mean", "expected_density_factor_std"]
    for ds in summary_datasets(summary):
        fields.extend((f"{ds}_particular_density_factor_mean", f"{ds}_particular_density_factor_std"))
    return fields


def bucket_rows(summary: list[dict]) -> list[dict]:
    """Dedicated bucket/binning view (mean +/- std over seeds), per dataset x
    order -- see scripts/export.py:parse_bucket_stats. Explains why a
    lower-density order isn't automatically faster to query: fewer buckets
    (from lower density) means a larger *average* bucket, and a bigger share
    of k-mers land in sshash's skew index (its structure for heavy buckets)."""
    out = []
    for r in summary:
        out.append(
            {
                "dataset": r.get("dataset"),
                "order": r.get("order"),
                "n_seeds": r.get("n_seeds"),
                "n_kmers": r.get("n_kmers"),
                "num_buckets_mean": r.get("num_buckets_mean"),
                "num_buckets_std": r.get("num_buckets_std"),
                "avg_bucket_size_mean": r.get("avg_bucket_size_mean"),
                "avg_bucket_size_std": r.get("avg_bucket_size_std"),
                "avg_bucket_size_ratio": r.get("avg_bucket_size_ratio"),
                "max_bucket_size_mean": r.get("max_bucket_size_mean"),
                "max_bucket_size_std": r.get("max_bucket_size_std"),
                "kmers_in_skew_index_pct_mean": r.get("kmers_in_skew_index_pct_mean"),
                "kmers_in_skew_index_pct_std": r.get("kmers_in_skew_index_pct_std"),
                "skew_index_bits_per_kmer_mean": r.get("skew_index_bits_per_kmer_mean"),
                "skew_index_bits_per_kmer_std": r.get("skew_index_bits_per_kmer_std"),
                "mid_load_bits_per_kmer_mean": r.get("mid_load_bits_per_kmer_mean"),
                "control_codewords_bits_per_kmer_mean": r.get("control_codewords_bits_per_kmer_mean"),
            }
        )
    return out


def bucket_field_order() -> list[str]:
    return [
        "dataset",
        "order",
        "n_seeds",
        "n_kmers",
        "num_buckets_mean",
        "num_buckets_std",
        "avg_bucket_size_mean",
        "avg_bucket_size_std",
        "avg_bucket_size_ratio",
        "max_bucket_size_mean",
        "max_bucket_size_std",
        "kmers_in_skew_index_pct_mean",
        "kmers_in_skew_index_pct_std",
        "skew_index_bits_per_kmer_mean",
        "skew_index_bits_per_kmer_std",
        "mid_load_bits_per_kmer_mean",
        "control_codewords_bits_per_kmer_mean",
    ]


def datasets_rows(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    seen: set[str] = set()
    for r in rows:
        ds = r.get("dataset")
        if not ds or ds in seen:
            continue
        seen.add(ds)
        out.append({"dataset": ds, "k": r.get("k"), "m": r.get("m"), "num_kmers": r.get("num_kmers")})
    return sorted(out, key=lambda r: r["dataset"])


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def cmd_evaluate(args: argparse.Namespace) -> None:
    out_root = Path(args.results_dir).resolve()
    runs_root = out_root / "runs"
    rows: list[dict] = []
    for run_dir in iter_run_dirs(runs_root):
        row = load_run_row(run_dir, runs_root)
        if row:
            rows.append(row)
    if not rows:
        die(f"no completed runs in {runs_root}")

    n_cmp, errors = validate_orders_equivalence(rows)
    if errors:
        for msg in errors:
            print(f"error: {msg}", file=sys.stderr)
        die(f"{len(errors)} equivalence/check failure(s) out of {n_cmp} run(s); see above")

    rows.sort(key=lambda r: (r.get("dataset") or "", r.get("order") or "", int(r.get("seed_idx") or 0)))
    summary = summary_rows(rows)
    attach_native_ratios(summary)
    vs_ref = vs_reference_rows(rows)

    m_val = next((r.get("m") for r in summary if r.get("m") is not None), None)
    w_val = next((r.get("w") for r in summary if r.get("w") is not None), None)
    expected: dict[str, dict] = {}
    if m_val is not None and w_val is not None:
        from scripts.config import default_gm_file

        expected = ensure_expected_orders(
            out_root, int(m_val), int(w_val), list(bench_seeds())[:5], gm_file=default_gm_file()
        )
    orders = orders_rows(summary, expected)
    buckets = bucket_rows(summary)

    write_csv(out_root / "sshash_runs.csv", rows, list(RUN_FIELDS))
    write_csv(out_root / "sshash_summary.csv", summary, summary_field_order())
    write_csv(out_root / "sshash_vs_reference.csv", vs_ref, vs_reference_field_order())
    write_csv(out_root / "sshash_buckets.csv", buckets, bucket_field_order())

    sheets = [
        ("index", ["sheet", "column", "meaning"], index_rows()),
        ("orders", orders_field_order(summary), orders),
        ("datasets", ["dataset", "k", "m", "num_kmers"], datasets_rows(rows)),
        ("summary", summary_field_order(), summary),
        ("buckets", bucket_field_order(), buckets),
        ("runs", list(RUN_FIELDS), rows),
        ("vs_reference", vs_reference_field_order(), vs_ref),
    ]
    write_xlsx(out_root / "sshash_results.xlsx", sheets)

    print(f"equivalence OK: {n_cmp} run(s) checked across {len({r['dataset'] for r in rows})} dataset(s)")


def main(argv: list[str] | None = None) -> int:
    load_dotenv()
    _, default_results = bench_paths()
    p = argparse.ArgumentParser(prog="scripts export", description="Export sshash benchmark results to CSV and XLSX")
    p.add_argument("--results-dir", default=str(default_results))
    args = p.parse_args(argv)
    cmd_evaluate(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
