#!/usr/bin/env python3
"""Export sshash benchmark runs to CSV and XLSX (stdlib only).

Reads ``<results>/runs`` and writes the tables to ``<results>/csvs`` and the workbook to
``<results>`` (``--results-dir`` defaults to ``$SSHASH_RESULTS_DIR`` or
``$SSHASH_BENCH_ROOT/results``):

- ``sshash_runs.csv``       -- one row per run (metrics + telemetry)
- ``sshash_orders.csv``     -- the ``orders`` sheet (one row per order variant, k/w/L)
- ``sshash_summary.csv``    -- one row per dataset x order (mean/std across seeds)
- ``sshash_density_expected_orders.csv`` -- cached expected-density Monte Carlo
  results (independent of any run; see scripts/expected_density.py)
- ``sshash_buckets.csv``    -- bucket/binning stats (mean/std across seeds) per
  dataset x order: bucket count, average/max bucket size, skew-index share,
  and the bucket-related space breakdown (see scripts/export.py:parse_bucket_stats)
- ``sshash_results.xlsx``   -- ``overview``, ``orders``, ``datasets``, ``summary``,
  ``runs``, ``buckets``, ``space_breakdown``, ``bucket_histogram``,
  ``super_kmer_histogram`` (written by benchmark/common/xlsx_book.py). The
  ``orders`` sheet carries each order variant's description/CLI plus its density --
  particular (measured, per dataset x variant) and expected (theoretical, per
  variant) -- side by side.

An order variant is identified by the five fields ``order, gm_k, projection,
gm_table, mod_t`` (benchmark/common/order_names.py): ``order`` is the canonical name
(native | xor | xor-pcg-hash | greedymini | mod-minimizer-xor |
mod-minimizer-greedymini), the others say which GreedyMini table / projection /
mod-minimizer t-mer length it ran with. Every grouping key below includes all five.

Density: ``density_factor = (n_minimizer_positions / n_kmers) * (w + 1)``,
where ``w = L - k + 1`` (=11 for L=31, k=21; the "w11" in
``orders/GreedyMini/w11_k15_canonical``) and n_minimizer_positions is the
number of super-k-mers (minimizer-value change points) sshash build
--verbose reports. The raw ratio isn't kept as its own column --
density_factor is the reported metric. A uniformly random order's
density_factor -> 2 as w grows; GreedyMini targets a lower value by
construction.

Buckets/binning: a lower-density order has fewer distinct minimizers, hence
fewer buckets (``n_buckets``) -- but that means a larger *average* bucket
(``avg_bucket_size = n_kmers / n_buckets``) and often a bigger share of
k-mers in sshash's skew index (``kmers_in_skew_index_pct``, its structure for
large/heavy buckets), which can offset a lower density's query-time benefit
even when it still saves per-bucket space (``control_codewords``/
``mid_load_bits_per_kmer``).

Exits with an error (hard-fail gate, mirroring kmc's validate_custom_vs_lex)
if:

- ``n_kmers`` differs across any order/seed run of the same dataset, or
- any run's ``query_pos_match`` / ``query_neg_match`` fixed-set check is
  false (see scripts/equivalence.py), or
- any run's ``check_ok`` (``sshash check`` -> "EVERYTHING OK!") is false.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import sys
from collections.abc import Sequence
from pathlib import Path

from scripts.common import bench_paths, die, load_dotenv
from scripts.config import seeds as bench_seeds
from scripts.density import density_factor as _density_factor
from scripts.expected_density import ensure_expected_orders
from scripts.config import gm_table_file
from scripts.orders import BASELINES, ORDERS, order_dir_name, variant_key
import order_names as on
import orders_sheet as osh
from units import gb_from_bytes, gb_from_kib
from xlsx_book import DATA_START_ROW, write_book
from xlsx_util import col_letter, write_csv as _shared_write_csv, xml_escape


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
        "minimizer_positions_in_skew_index": None,
        "minimizer_positions_in_skew_index_pct": None,
        "num_bits_per_offset": None,
        "mphf_bits_per_kmer": None,
        "strings_offsets_bits_per_kmer": None,
        "strings_bits_per_kmer": None,
        "begin_buckets_bits_per_kmer": None,
    }
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")

    # Occurrences (super-k-mer start positions) of the heavy minimizers: the
    # entries of H (heavy_load_buckets), each num_bits_per_offset bits wide.
    m = re.search(
        r"^num_minimizer_positions_of_buckets_in_skew_index (\d+)/\d+ \(([\d.eE+-]+)%\)",
        text,
        re.MULTILINE,
    )
    if m:
        out["minimizer_positions_in_skew_index"] = int(m.group(1))
        out["minimizer_positions_in_skew_index_pct"] = float(m.group(2))
    m = re.search(r"^num_bits_per_offset = (\d+)", text, re.MULTILINE)
    if m:
        out["num_bits_per_offset"] = int(m.group(1))

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
    # rest of the SPACE BREAKDOWN (dictionary::print_info): the order-independent strings and
    # offsets, the MPHF over the minimizers, and the (negligible) bucket-size start array
    for key, label in (
        ("mphf_bits_per_kmer", "mphf"),
        ("strings_offsets_bits_per_kmer", "strings_offsets"),
        ("strings_bits_per_kmer", "strings"),
        ("begin_buckets_bits_per_kmer", "begin_buckets_of_size"),
    ):
        m = re.search(rf"^\s+{label}: ([\d.eE+-]+) \[bits/kmer\]", text, re.MULTILINE)
        if m:
            out[key] = float(m.group(1))
    return out


HISTOGRAM_ID_FIELDS = ("dataset", "k", "w", "L", *on.ID_FIELDS, "modality", "seed_idx")
BUCKET_HISTOGRAM_FIELDS = [
    *HISTOGRAM_ID_FIELDS, "bucket_type", "log2_bucket_size", "bucket_size_range",
    "n_buckets", "n_kmers_in_buckets",
]
SUPER_KMER_HISTOGRAM_FIELDS = [*HISTOGRAM_ID_FIELDS, "kmers_in_super_kmer", "n_super_kmers"]
HEAVY_KMER_FRACTION = 0.01  # top share of k-mers treated as "large bucket", all at one (capped) cost


def parse_bucket_histograms(path: Path) -> dict:
    """Parse the full histograms `sshash build --verbose --bucket-stats` prints via
    buckets_statistics::print_full: per bucket size (= number of super-k-mers
    sharing one sampled m-mer) the number of buckets and of k-mers, and per
    super-k-mer length the number of super-k-mers. Buckets larger than
    buckets_statistics::MAX_BUCKET_SIZE (4096) are not tallied by sshash."""
    out: dict = {"bucket_sizes": [], "super_kmer_lengths": [], "heavy_lower": 64}
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")
    # 2^l: buckets above it are heavy (skew-index partition 0 starts there)
    m = re.search(r"^\s*partition = 0: num kmers in buckets of size > (\d+)", text, re.MULTILINE)
    if m:
        out["heavy_lower"] = int(m.group(1))
    for m in re.finditer(
        r"^buckets with (\d+) minimizer positions = (\d+) \([^)]*\) \| num_kmers = (\d+)",
        text, re.MULTILINE,
    ):
        out["bucket_sizes"].append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    for m in re.finditer(r"^super_kmers with (\d+) kmer = (\d+) ", text, re.MULTILINE):
        out["super_kmer_lengths"].append((int(m.group(1)), int(m.group(2))))
    return out


BUCKET_HISTOGRAM_CAP = 4096  # sshash's buckets_statistics::MAX_BUCKET_SIZE: bigger buckets are not tallied per size


def histogram_overflow(
    bucket_sizes: list[tuple[int, int, int]], n_kmers: int | None, n_buckets: int | None
) -> tuple[int, int] | None:
    """(buckets, kmers) in the buckets larger than BUCKET_HISTOGRAM_CAP, which
    `build --bucket-stats` does not list per size. Exact: every k-mer sits in
    exactly one bucket and every sampled m-mer owns one bucket, so the overflow
    is the run's totals minus what the histogram lists. None if it cannot be
    derived (no histogram / totals) or the numbers are inconsistent."""
    if not bucket_sizes or n_kmers is None or n_buckets is None:
        return None
    kmers = n_kmers - sum(n for _, _, n in bucket_sizes)
    buckets = n_buckets - sum(b for _, b, _ in bucket_sizes)
    if kmers < 0 or buckets < 0 or (kmers == 0) != (buckets == 0):
        return None
    return (buckets, kmers) if kmers else None


def log2_bucket_bins(
    bucket_sizes: list[tuple[int, int, int]],
    overflow: tuple[int, int] | None = None,
    heavy_lower: int = 64,
) -> list[dict]:
    """Bucket-size histogram in the paper's logarithmic partitioning: bin b holds
    the sizes in (2^(b-1), 2^b] (b = ceil(log2 size); b = 0 is the singleton
    size 1), the same intervals as the skew index partitions K_i = (2^i, 2^(i+1)].
    `bucket_type` follows Section 6: singleton (size 1), light (2 .. 2^l) and
    heavy (more than 2^l, `heavy_lower` = 2^l), so the heavy bins are exactly the
    skew-index partitions. The buckets above BUCKET_HISTOGRAM_CAP (`overflow`,
    see histogram_overflow) form one last row, since sshash does not tally them
    per size (it covers the bins from log2 = 13 up)."""
    bins: dict[int, list[int]] = {}
    for size, buckets, kmers in bucket_sizes:
        cur = bins.setdefault((size - 1).bit_length(), [0, 0])
        cur[0] += buckets
        cur[1] += kmers

    def kind(upper: int) -> str:
        return "singleton" if upper == 1 else "light" if upper <= heavy_lower else "heavy"

    out = []
    for b in sorted(bins):
        lo, hi = (1, 1) if b == 0 else ((1 << (b - 1)) + 1, 1 << b)
        out.append({
            "bucket_type": kind(hi), "log2_bucket_size": b,
            "bucket_size_range": str(lo) if lo == hi else f"{lo}-{hi}",
            "n_buckets": bins[b][0], "n_kmers_in_buckets": bins[b][1],
        })
    if overflow:
        out.append({
            "bucket_type": "heavy", "log2_bucket_size": (BUCKET_HISTOGRAM_CAP - 1).bit_length() + 1,
            "bucket_size_range": f">{BUCKET_HISTOGRAM_CAP}",
            "n_buckets": overflow[0], "n_kmers_in_buckets": overflow[1],
        })
    return out


def bucket_distribution_metrics(
    bucket_sizes: list[tuple[int, int, int]], overflow: tuple[int, int] | None = None
) -> dict:
    """Skew of the bucket-size distribution from the (size, buckets, kmers)
    histogram. `expected_bucket_cost` models the query cost per k-mer as the
    number of super-k-mers to check (its bucket size), except that the
    HEAVY_KMER_FRACTION of k-mers in the largest buckets all cost the same,
    `heavy_bucket_threshold` (the smallest size s with at most that share of
    k-mers in buckets larger than s). `overflow` = (buckets, kmers) of the
    buckets above BUCKET_HISTOGRAM_CAP (see histogram_overflow): they are the
    largest class, so the shares and the heavy threshold use the true totals."""
    if overflow:
        bucket_sizes = [*bucket_sizes, (BUCKET_HISTOGRAM_CAP + 1, *overflow)]
    out: dict = {
        "singleton_bucket_pct": None,
        "kmers_in_singleton_bucket_pct": None,
        "heavy_bucket_threshold": None,
        "expected_bucket_cost": None,
    }
    total_buckets = sum(b for _, b, _ in bucket_sizes)
    total_kmers = sum(n for _, _, n in bucket_sizes)
    if not total_buckets or not total_kmers:
        return out
    by_size = {s: (b, n) for s, b, n in bucket_sizes}
    out["singleton_bucket_pct"] = round(100.0 * by_size.get(1, (0, 0))[0] / total_buckets, 4)
    out["kmers_in_singleton_bucket_pct"] = round(100.0 * by_size.get(1, (0, 0))[1] / total_kmers, 4)
    threshold, above = max(by_size), 0
    for s in sorted(by_size, reverse=True):
        if above + by_size[s][1] > HEAVY_KMER_FRACTION * total_kmers:
            threshold = s
            break
        above += by_size[s][1]
    out["heavy_bucket_threshold"] = threshold
    out["expected_bucket_cost"] = round(
        sum(n * min(s, threshold) for s, _, n in bucket_sizes) / total_kmers, 4
    )
    return out


def _extension_pct(query_json: dict) -> float | None:
    """Share of positive k-mers of the streamed reads answered by extension
    along the stored string (sshash query's num_extensions) rather than a
    seed lookup (num_searches)."""
    searches = _num(query_json.get("num_searches"), int)
    extensions = _num(query_json.get("num_extensions"), int)
    if searches is None or extensions is None or searches + extensions == 0:
        return None
    return round(100.0 * extensions / (searches + extensions), 4)


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
    "w",
    "L",
    "order",
    "gm_k",
    "projection",
    "gm_table",
    "mod_t",
    "modality",
    "threads",
    "seed",
    "seed_idx",
    # space
    "n_kmers",
    "index_size_gb",
    "index_disk_gb",
    "bits_per_kmer",
    # density (the raw minimizer_density stays dropped -- density_factor derives
    # from it; the occurrence count Z and the repetition factor Z/M are kept
    # because density fixes Z but M = Z / (Z/M) also depends on the data)
    "density_factor",
    "n_minimizer_positions",
    "minimizer_positions_per_minimizer",
    # buckets
    "n_buckets",
    "avg_bucket_size",
    "max_bucket_size",
    "log2_max_bucket_size",
    "kmers_in_skew_index_pct",
    "n_buckets_in_skew_index",
    "minimizer_positions_in_skew_index",
    "minimizer_positions_in_skew_index_pct",
    "num_bits_per_offset",
    "skew_index_bits_per_kmer",
    "skew_h_bits_per_kmer",
    "skew_mphf_v_bits_per_kmer",
    "mid_load_bits_per_kmer",
    "control_codewords_bits_per_kmer",
    "mphf_bits_per_kmer",
    "strings_offsets_bits_per_kmer",
    "strings_bits_per_kmer",
    "begin_buckets_bits_per_kmer",
    "singleton_bucket_pct",
    "kmers_in_singleton_bucket_pct",
    "heavy_bucket_threshold",
    "expected_bucket_cost",
    # build
    "build_wall_s",
    "build_max_rss_gb",
    "build_steps_total_s",
    "build_step1_s",
    "build_step2_s",
    "build_step3_s",
    "build_step4_s",
    "build_step5_s",
    "build_step6_s",
    "build_step7_s",
    # query (real FASTQ) -- kept next to build, the other real wall-clock op
    "query_wall_s",
    "query_max_rss_gb",
    "query_elapsed_ms",
    "query_n_kmers",
    "query_n_positive_kmers",
    "query_n_negative_kmers",
    "query_n_searches",
    "query_n_extensions",
    "query_extension_pct",
    "query_n_skipped_singleton_lookups",
    "query_n_bucket_cache_hits",
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


# Datasets that are real runs on disk (so --dataset-file/-d smoke sweeps work
# for local smoke testing, see README "Smoke testing") but never belong in
# the shipped CSV/XLSX reports -- excluded here, not by hand-editing the
# exported files, so every export (CLI, `run`'s automatic export, re-export)
# stays consistent.
EXCLUDED_DATASETS = frozenset({"smoke"})


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

    n_kmers = _num(build_json.get("num_kmers"), int)
    index_bytes_json = _num(build_json.get("index_size_in_bytes"), int)
    index_bytes_disk = meta.get("index_bytes_on_disk")
    bits_per_kmer = None
    if index_bytes_disk and n_kmers:
        bits_per_kmer = round(index_bytes_disk * 8 / n_kmers, 4)

    mini = parse_minimizer_counts(run_dir / "build.log")
    buckets = parse_bucket_stats(run_dir / "build.log")
    avg_bucket_size = None
    if buckets["num_buckets"] and n_kmers:
        avg_bucket_size = round(n_kmers / buckets["num_buckets"], 2)
    hist_sizes = parse_bucket_histograms(run_dir / "build.log")["bucket_sizes"]
    bucket_dist = bucket_distribution_metrics(
        hist_sizes, histogram_overflow(hist_sizes, n_kmers, buckets["num_buckets"])
    )
    # meta.json keeps sshash's own names: k = k-mer length (-> L), m = minimizer length (-> k)
    L_val, k_val = meta.get("k"), meta.get("m")
    w = int(L_val) - int(k_val) + 1 if L_val is not None and k_val is not None else None
    # Raw minimizer_density (= n_minimizer_positions / n_kmers) isn't
    # kept on the row -- density_factor (that times (w+1)) is enough.
    density_fac = None
    if mini["num_minimizer_positions"] and n_kmers and w is not None:
        density_fac = _density_factor(mini["num_minimizer_positions"], n_kmers, w)

    # Z / M: occurrences per distinct minimizer -- how much the data repeats a
    # sampled m-mer (num_minimizers == num_buckets, M).
    positions_per_minimizer = None
    if mini["num_minimizer_positions"] and mini["num_minimizers"]:
        positions_per_minimizer = round(
            mini["num_minimizer_positions"] / mini["num_minimizers"], 4
        )
    # H = heavy_load_buckets: one num_bits_per_offset-wide entry per occurrence of a
    # heavy minimizer (build_sparse_and_skew_index.cpp), so its size is exact; the
    # rest of the skew index is the per-partition MPHFs f_i and the V_i arrays.
    skew_h_bits = skew_mphf_v_bits = None
    if (
        buckets["minimizer_positions_in_skew_index"] is not None
        and buckets["num_bits_per_offset"]
        and n_kmers
    ):
        skew_h_bits = round(
            buckets["minimizer_positions_in_skew_index"] * buckets["num_bits_per_offset"] / n_kmers,
            4,
        )
        if buckets["skew_index_bits_per_kmer"] is not None:
            skew_mphf_v_bits = round(buckets["skew_index_bits_per_kmer"] - skew_h_bits, 4)

    def step(key: str) -> float | None:
        v = _num(build_json.get(key), float)
        return round(v / 1_000_000, 6) if v is not None else None

    modality = meta.get("modality") or "canonical"
    try:
        ident = on.parse_id({f: meta.get(f) for f in on.ID_FIELDS}, baselines=BASELINES)
    except ValueError as e:
        die(f"{meta_path}: {e}")

    row = {
        "run_id": meta.get("run_id") or str(run_dir.relative_to(runs_root)),
        "dataset": meta.get("dataset"),
        "k": k_val,
        "w": w,
        "L": L_val,
        **ident,
        "modality": modality,
        "threads": meta.get("threads"),
        "seed": meta.get("order_seed"),
        "seed_idx": meta.get("seed_idx"),
        "n_kmers": n_kmers,
        "index_size_gb": gb_from_bytes(index_bytes_json),
        "index_disk_gb": gb_from_bytes(index_bytes_disk),
        "bits_per_kmer": bits_per_kmer,
        "density_factor": density_fac,
        "n_minimizer_positions": mini["num_minimizer_positions"],
        "minimizer_positions_per_minimizer": positions_per_minimizer,
        "n_buckets": buckets["num_buckets"],
        "avg_bucket_size": avg_bucket_size,
        "max_bucket_size": buckets["max_bucket_size"],
        "log2_max_bucket_size": buckets["log2_max_bucket_size"],
        "kmers_in_skew_index_pct": buckets["kmers_in_skew_index_pct"],
        "n_buckets_in_skew_index": buckets["buckets_in_skew_index"],
        "minimizer_positions_in_skew_index": buckets["minimizer_positions_in_skew_index"],
        "minimizer_positions_in_skew_index_pct": buckets["minimizer_positions_in_skew_index_pct"],
        "num_bits_per_offset": buckets["num_bits_per_offset"],
        "skew_index_bits_per_kmer": buckets["skew_index_bits_per_kmer"],
        "skew_h_bits_per_kmer": skew_h_bits,
        "skew_mphf_v_bits_per_kmer": skew_mphf_v_bits,
        "mphf_bits_per_kmer": buckets["mphf_bits_per_kmer"],
        "strings_offsets_bits_per_kmer": buckets["strings_offsets_bits_per_kmer"],
        "strings_bits_per_kmer": buckets["strings_bits_per_kmer"],
        "begin_buckets_bits_per_kmer": buckets["begin_buckets_bits_per_kmer"],
        "mid_load_bits_per_kmer": buckets["mid_load_bits_per_kmer"],
        "control_codewords_bits_per_kmer": buckets["control_codewords_bits_per_kmer"],
        **bucket_dist,
        "build_wall_s": build_time["time_elapsed_s"],
        "build_max_rss_gb": gb_from_kib(build_time["max_rss_kb"]),
        "build_steps_total_s": (
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
        "query_max_rss_gb": gb_from_kib(query_time["max_rss_kb"]),
        "query_elapsed_ms": _num(query_json.get("elapsed_millisec")),
        "query_n_kmers": _num(query_json.get("num_kmers"), int),
        "query_n_positive_kmers": _num(query_json.get("num_positive_kmers"), int),
        "query_n_negative_kmers": _num(query_json.get("num_negative_kmers"), int),
        "query_n_searches": _num(query_json.get("num_searches"), int),
        "query_n_extensions": _num(query_json.get("num_extensions"), int),
        "query_extension_pct": _extension_pct(query_json),
        "query_n_skipped_singleton_lookups": _num(
            query_json.get("num_skipped_singleton_lookups"), int
        ),
        "query_n_bucket_cache_hits": _num(query_json.get("num_bucket_cache_hits"), int),
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


def order_label(r: dict) -> str:
    """Readable id of an order variant (its run-directory segment)."""
    return order_dir_name(r.get("order") or "", r.get("gm_k"), r.get("projection"), r.get("mod_t"))


def validate_orders_equivalence(rows: list[dict]) -> tuple[int, list[str]]:
    """Every order/seed run of the same dataset must report the same
    ``n_kmers`` and pass its check/query-answer gates. Returns
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
        counts = {r.get("n_kmers") for r in group if r.get("n_kmers") is not None}
        n_cmp += len(group)
        if len(counts) > 1:
            detail = ", ".join(
                f"{order_label(r)}"
                f"/seed{r.get('seed_idx')}={r.get('n_kmers')}"
                for r in group
            )
            errors.append(f"{dataset}: n_kmers differs across orders/seeds: {detail}")
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
    # buckets (n_buckets, M, is also shown in the summary: the correlation rows use it)
    "n_buckets",
    "minimizer_positions_per_minimizer",
    "avg_bucket_size",
    "max_bucket_size",
    "kmers_in_skew_index_pct",
    "n_buckets_in_skew_index",
    "minimizer_positions_in_skew_index",
    "minimizer_positions_in_skew_index_pct",
    "skew_index_bits_per_kmer",
    "skew_h_bits_per_kmer",
    "skew_mphf_v_bits_per_kmer",
    "mid_load_bits_per_kmer",
    "control_codewords_bits_per_kmer",
    "singleton_bucket_pct",
    "kmers_in_singleton_bucket_pct",
    # build
    "build_wall_s",
    "build_max_rss_gb",
    # query (real FASTQ) -- kept next to build
    "query_wall_s",
    "query_elapsed_ms",
    "query_extension_pct",
    # bench (sshash's own internal micro-benchmark)
    "bench_positive_ns",
    "bench_negative_ns",
    "bench_access_ns",
)


def _round_std(x: float | None) -> float | None:
    """Round a standard deviation to at least 4 decimals and at least 4
    significant digits, so a tiny std (a GreedyMini order's 1e-5 spread over
    seeds) is not flattened to 0.0."""
    if x is None:
        return None
    if x == 0:
        return 0.0
    return round(x, max(4, 3 - int(math.floor(math.log10(abs(x))))))


def _mean_std(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return statistics.mean(values), statistics.stdev(values)


def summary_group_key(r: dict) -> tuple:
    return (
        r.get("dataset") or "",
        int(r.get("k") or 0),
        int(r.get("w") or 0),
        int(r.get("L") or 0),
        r.get("modality") or "",
        *variant_key(r),
    )


def variant_sort_key(r: dict) -> tuple:
    """Canonical order rank, then the variant's gm_k / projection / table / t."""
    try:
        rank = ORDERS.index(r.get("order") or "")
    except ValueError:
        rank = len(ORDERS)
    return (rank, r.get("gm_k") or 0, r.get("projection") or "", r.get("gm_table") or "",
            r.get("mod_t") or 0)


def row_sort_key(r: dict) -> tuple:
    modality_rank = 0 if (r.get("modality") or "canonical") == "canonical" else 1
    return (r.get("dataset") or "", modality_rank, *variant_sort_key(r))


def summary_sort_key(row: dict) -> tuple:
    return row_sort_key(row)


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
            "w": head.get("w"),
            "L": head.get("L"),
            **{f: head.get(f) for f in on.ID_FIELDS},
            "modality": head["modality"],
            "n_seeds": len({p.get("seed_idx") for p in parts}),
            "n_kmers": head.get("n_kmers"),
        }
        for metric in SUMMARY_METRICS:
            vals = [float(p[metric]) for p in parts if p.get(metric) is not None]
            mean, std = _mean_std(vals)
            row[f"{metric}_mean"] = round(mean, 4) if mean is not None else None
            row[f"{metric}_std"] = _round_std(std)
        out.append(row)
    return sorted(out, key=summary_sort_key)


def attach_native_ratios(rows: list[dict], baseline: str = "native") -> None:
    """`{metric}_ratio` = this variant's mean / the same-modality `baseline` order's
    mean for the same dataset (baseline = 1). Isolates each order's own effect from
    the canonical-vs-regular parsing-mode difference. The baseline ("native") can run
    under either modality, so it is looked up per (dataset, modality); it falls back
    to the canonical baseline if no regular-modality one was run."""
    base_by_key: dict[tuple, dict] = {}
    canonical_native_by_dataset: dict[str, dict] = {}
    for r in rows:
        if r.get("order") == baseline:
            base_by_key[(r["dataset"], r.get("modality"))] = r
            if r.get("modality") == "canonical":
                canonical_native_by_dataset[r["dataset"]] = r
    for r in rows:
        modality = r.get("modality") or "canonical"
        base = base_by_key.get((r["dataset"], modality)) or canonical_native_by_dataset.get(
            r["dataset"]
        )
        for metric in SUMMARY_METRICS:
            mean = r.get(f"{metric}_mean")
            base_mean = base.get(f"{metric}_mean") if base else None
            rname = f"{metric}_ratio"
            if mean is None or not base_mean:
                r[rname] = None
            else:
                r[rname] = round(mean / base_mean, 4)


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


# Bucket metrics that the summary sheet shows as well as the buckets sheet:
# the number of distinct minimizers M is a driver of its correlation rows.
SUMMARY_KEEP_BUCKET_METRICS = frozenset({"n_buckets"})


def summary_field_order() -> list[str]:
    fields = ["dataset", "k", "w", "L", *on.ID_FIELDS, "modality", "n_seeds", "n_kmers"]
    for metric in SUMMARY_METRICS:
        if _METRIC_CATEGORY.get(metric) == "buckets" and metric not in SUMMARY_KEEP_BUCKET_METRICS:
            continue  # shown in the buckets sheet (bucket_field_order), not repeated here
        fields.extend((f"{metric}_mean", f"{metric}_std", f"{metric}_ratio"))
    return fields


# ---------------------------------------------------------------------------
# XLSX (written by benchmark/common/xlsx_book.write_book: category band in row 1,
# column headers in row 2, data from DATA_START_ROW, vertical merges of the identity
# columns inside a dataset block, a border after each block, colour scales on the
# `_ratio` columns)
# ---------------------------------------------------------------------------

# Repeated across every row of a dataset's block (k/w/L are fixed for the
# whole benchmark; dataset obviously repeats) -- merged vertically instead of
# printing the same value on every row.
XLSX_MERGE_FIELDS = ("dataset", "k", "w", "L", "modality")
# the columns that name an order variant (order, gm_k, projection, gm_table, mod_t)
XLSX_VARIANT_FIELDS = on.ID_FIELDS

_ID_FIELDS = frozenset({
    "run_id", "dataset", "k", "w", "L", *on.ID_FIELDS, "modality", "threads", "seed", "seed_idx",
    "n_seeds", "n_kmers",
})
_METRIC_CATEGORY = {
    # space
    "n_kmers": "space", "index_size_gb": "space", "index_disk_gb": "space", "bits_per_kmer": "space",
    # density
    "density_factor": "density", "n_minimizer_positions": "density",
    "minimizer_positions_per_minimizer": "buckets",
    # buckets
    "n_buckets_in_skew_index": "buckets", "minimizer_positions_in_skew_index": "buckets",
    "minimizer_positions_in_skew_index_pct": "buckets", "num_bits_per_offset": "buckets",
    "skew_h_bits_per_kmer": "buckets", "skew_mphf_v_bits_per_kmer": "buckets",
    "mphf_bits_per_kmer": "space", "strings_offsets_bits_per_kmer": "space",
    "strings_bits_per_kmer": "space", "begin_buckets_bits_per_kmer": "space",
    "n_buckets": "buckets", "avg_bucket_size": "buckets", "max_bucket_size": "buckets",
    "log2_max_bucket_size": "buckets", "kmers_in_skew_index_pct": "buckets",
    "skew_index_bits_per_kmer": "buckets", "mid_load_bits_per_kmer": "buckets",
    "control_codewords_bits_per_kmer": "buckets",
    "singleton_bucket_pct": "buckets", "kmers_in_singleton_bucket_pct": "buckets",
    "heavy_bucket_threshold": "buckets", "expected_bucket_cost": "buckets",
    # build
    "build_wall_s": "build", "build_max_rss_gb": "build", "build_steps_total_s": "build",
    "build_step1_s": "build", "build_step2_s": "build", "build_step3_s": "build", "build_step4_s": "build",
    "build_step5_s": "build", "build_step6_s": "build", "build_step7_s": "build",
    # correctness
    "check_ok": "correctness",
    # bench (sshash's own internal micro-benchmark)
    "bench_positive_ns": "bench", "bench_negative_ns": "bench", "bench_access_ns": "bench", "bench_iterator_ns": "bench",
    # query (real FASTQ)
    "query_wall_s": "query", "query_max_rss_gb": "query", "query_elapsed_ms": "query",
    "query_n_kmers": "query", "query_n_positive_kmers": "query", "query_n_negative_kmers": "query",
    "query_n_searches": "query", "query_n_extensions": "query", "query_extension_pct": "query",
    "query_n_skipped_singleton_lookups": "query", "query_n_bucket_cache_hits": "query",
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
    if field in (*on.ID_FIELDS, "description", "cli", "k", "w", "L"):
        return "order"
    if field.startswith("expected_density_factor"):
        return "expected (theoretical)"
    marker = "_particular_density_factor_"
    if marker in field:
        return field.split(marker, 1)[0]
    return ""


# Row-1 groups of the buckets sheet, in column order (see BUCKET_METRIC_GROUPS).
BUCKET_METRIC_GROUPS = (
    # every minimizer has one tag (control codeword); the three bucket kinds of the paper follow
    ("all buckets (one tag each)", (
        "n_buckets", "minimizer_positions_per_minimizer", "avg_bucket_size", "max_bucket_size",
        "control_codewords_bits_per_kmer",
    )),
    ("singleton buckets (position in the tag)", (
        "singleton_bucket_pct", "kmers_in_singleton_bucket_pct",
    )),
    ("light buckets (array L)", ("mid_load_bits_per_kmer",)),
    ("heavy buckets (skew index: f_i, V_i, H)", (
        "kmers_in_skew_index_pct", "n_buckets_in_skew_index",
        "minimizer_positions_in_skew_index", "minimizer_positions_in_skew_index_pct",
        "skew_index_bits_per_kmer", "skew_h_bits_per_kmer", "skew_mphf_v_bits_per_kmer",
    )),
)
BUCKET_METRICS = tuple(m for _, group in BUCKET_METRIC_GROUPS for m in group)
# Green/red scale only where lower is better (space, number of buckets); the other
# bucket ratios (average size, singleton share, ...) are descriptive, not a gain or a loss.
BUCKET_COLORED_RATIOS = frozenset({
    "n_buckets", "control_codewords_bits_per_kmer", "mid_load_bits_per_kmer",
    "skew_index_bits_per_kmer", "skew_h_bits_per_kmer", "skew_mphf_v_bits_per_kmer",
})


def _bucket_field_category(field: str) -> str:
    if field in ("dataset", *on.ID_FIELDS, "modality", "n_seeds", "n_kmers"):
        return "id / params"
    base = field
    for suffix in ("_mean", "_std", "_ratio"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    for label, group in BUCKET_METRIC_GROUPS:
        if base in group:
            return label
    return ""


def sheet_categories(raw_name: str, fields: Sequence[str]) -> list[str] | None:
    """One category label per field, for the row-1 band; None to skip it."""
    if raw_name == "buckets":
        return [_bucket_field_category(f) for f in fields]
    if raw_name == "space_breakdown":
        return [_space_field_category(f) for f in fields]
    if raw_name == "summary":
        # n_buckets (M) next to the density column, as its own group
        return ["minimizers" if f.startswith("n_buckets_") else _field_category(f) for f in fields]
    if raw_name == "runs":
        return [_field_category(f) for f in fields]
    if raw_name == "orders":
        return [_orders_field_category(f) for f in fields]
    if raw_name == "overview":
        return ["reference"] * len(fields)
    if raw_name == "datasets":
        return ["dataset reference"] * len(fields)
    return None


def _pearson(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Pearson correlation; None when undefined (fewer than 3 points or a constant series)."""
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 1e-18 or syy <= 1e-18:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return round(max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy))), 3)


# Two rows right under the summary table: for every `{metric}_ratio` column, its
# Pearson correlation (Excel's built-in CORREL, so it stays live and traceable)
# with the density ratio and with the distinct-minimizer ratio, over the dataset
# x order rows. The ratios (this order / native, same dataset) are
# dimensionless, so rows of different datasets are comparable.
CORRELATION_DRIVERS = (
    ("correlation with density", "density_factor_ratio"),
    ("correlation with distinct minimizers", "n_buckets_ratio"),
)
# The buckets sheet gets one such row, for its tag and light-array space ratios.
BUCKET_CORRELATION_DRIVERS = (("correlation with distinct minimizers", "n_buckets_ratio"),)
BUCKET_CORRELATION_COLUMNS = frozenset({
    "control_codewords_bits_per_kmer_ratio", "mid_load_bits_per_kmer_ratio",
})


def correlation_rows(
    fields: Sequence[str],
    summary: Sequence[dict],
    drivers: Sequence[tuple[str, str]] = CORRELATION_DRIVERS,
    columns: frozenset | None = None,
) -> list[dict]:
    """One xlsx_book footer row per driver: {"label", "cells": [(field, formula,
    cached value)]} with the CORREL formula per `_ratio` column and its cached value
    (so a viewer that does not recalculate still shows it)."""
    n = len(summary)
    if not n:
        return []
    first, last = DATA_START_ROW, DATA_START_ROW + n - 1
    out = []
    for label, driver in drivers:
        dl = col_letter(list(fields).index(driver) + 1)
        cells: list[tuple[str, str, float | None]] = []
        for c, f in enumerate(fields):
            if not f.endswith("_ratio") or (columns is not None and f not in columns):
                continue
            cl = col_letter(c + 1)
            pairs = [(r[driver], r[f]) for r in summary if r.get(driver) is not None and r.get(f) is not None]
            value = _pearson([a for a, _ in pairs], [b for _, b in pairs])
            # a constant column has no correlation: Excel's CORREL gives #DIV/0!, shown blank
            formula = f'IFERROR(CORREL(${dl}${first}:${dl}${last},{cl}${first}:{cl}${last}),"")'
            cells.append((f, formula, value))
        out.append({"label": label, "cells": cells})
    return out


# ---------------------------------------------------------------------------
# Charts (hand-written DrawingML, same package layout openpyxl emits)
# ---------------------------------------------------------------------------

CHART_GRID_COL0 = len(BUCKET_HISTOGRAM_FIELDS) + 1  # first chart column: right of the raw table
CHART_PITCH = (9, 16)                # grid pitch in (columns, rows) -- room for one chart
CHART_SIZE_EMU = (5_400_000, 2_880_000)
CHART_DATASET_ORDER = ("cod", "kestrel", "human", "ncbi-virus", "se", "hprc")
# bar colours per canonical order; variants of one order (tables / projections) cycle through its list
CHART_COLORS = {
    "native": ["7F7F7F"], "xor": ["ED7D31"], "xor-pcg-hash": ["FFC000"],
    "greedymini": ["70AD47", "4472C4", "A9D18E", "8FAADC"],
    "mod-minimizer-xor": ["C00000"], "mod-minimizer-greedymini": ["7030A0", "B084D6"],
}


def bucket_histogram_charts(bucket_hist: list[dict]) -> list[dict]:
    """One chart per (dataset, order variant): number of buckets per log2 bucket-size
    bin (columns on a log axis, the last bin is the aggregated '>4096') for the
    first seed, in a grid -- datasets down, variants across -- to the right of the
    raw table (rows of `bucket_hist` start at DATA_START_ROW)."""
    sizes_col = col_letter(BUCKET_HISTOGRAM_FIELDS.index("bucket_size_range") + 1)
    count_col = col_letter(BUCKET_HISTOGRAM_FIELDS.index("n_buckets") + 1)
    blocks: dict[tuple, list[int]] = {}
    first_seed: dict[tuple, int] = {}
    variants: dict[tuple, dict] = {}
    for idx, r in enumerate(bucket_hist):
        key = (r["dataset"], variant_key(r))
        variants[variant_key(r)] = r
        seed = int(r.get("seed_idx") or 0)
        if key not in first_seed or seed < first_seed[key]:
            first_seed[key] = seed
            blocks[key] = []
        if seed == first_seed[key]:
            blocks[key].append(idx)
    datasets = [d for d in CHART_DATASET_ORDER if any(k[0] == d for k in blocks)]
    datasets += sorted({k[0] for k in blocks} - set(datasets))
    ordered = sorted(variants.values(), key=variant_sort_key)
    seen_per_order: dict[str, int] = {}
    colors: dict[tuple, str] = {}
    for v in ordered:
        palette = CHART_COLORS.get(v["order"], ["4472C4"])
        n = seen_per_order.get(v["order"], 0)
        seen_per_order[v["order"]] = n + 1
        colors[variant_key(v)] = palette[n % len(palette)]
    charts = []
    for i, ds in enumerate(datasets):
        for j, v in enumerate(ordered):
            rows = blocks.get((ds, variant_key(v)))
            if not rows:
                continue
            r0, r1 = DATA_START_ROW + rows[0], DATA_START_ROW + rows[-1]
            charts.append({
                "title": f"{ds} / {order_label(v)} (seed {first_seed[(ds, variant_key(v))]})",
                "cat": f"${sizes_col}${r0}:${sizes_col}${r1}",
                "val": f"${count_col}${r0}:${count_col}${r1}",
                "col": CHART_GRID_COL0 + j * CHART_PITCH[0],
                "row": 2 + i * CHART_PITCH[1],
                "color": colors[variant_key(v)],
            })
    return charts


def _chart_text(text: str, size: int = 1000, bold: bool = False) -> str:
    b = ' b="1"' if bold else ""
    return (
        '<tx><rich><a:bodyPr xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"/>'
        '<a:p xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
        f'<a:pPr><a:defRPr sz="{size}"{b}/></a:pPr><a:r><a:rPr lang="en-US" sz="{size}"{b}/>'
        f"<a:t>{xml_escape(text)}</a:t></a:r></a:p></rich></tx>"
    )


def _chart_xml(sheet_name: str, chart: dict) -> str:
    """Column chart of buckets per log2 bucket-size bin (category = the bin's
    size range, value axis log10 so the head of the distribution and the heavy
    tail are both visible)."""
    sheet = "'" + sheet_name.replace("'", "''") + "'"
    color = chart["color"]
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<chartSpace xmlns="http://schemas.openxmlformats.org/drawingml/2006/chart" '
        'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<chart><title>{_chart_text(chart["title"], 1000, True)}<overlay val="0"/></title>'
        '<autoTitleDeleted val="0"/><plotArea><layout/>'
        '<barChart><barDir val="col"/><grouping val="clustered"/><varyColors val="0"/>'
        f'<ser><idx val="0"/><order val="0"/><spPr><a:solidFill><a:srgbClr val="{color}"/></a:solidFill></spPr>'
        '<invertIfNegative val="0"/>'
        f'<cat><strRef><f>{sheet}!{chart["cat"]}</f></strRef></cat>'
        f'<val><numRef><f>{sheet}!{chart["val"]}</f></numRef></val></ser>'
        '<gapWidth val="40"/><axId val="111"/><axId val="222"/></barChart>'
        '<catAx><axId val="111"/><scaling><orientation val="minMax"/></scaling><delete val="0"/>'
        '<axPos val="b"/>'
        f'<title>{_chart_text("bucket size (log2 bins; heavy = skew-index partitions)", 900)}<overlay val="0"/></title>'
        '<numFmt formatCode="General" sourceLinked="0"/><majorTickMark val="out"/><minorTickMark val="none"/>'
        '<tickLblPos val="nextTo"/><crossAx val="222"/><crosses val="autoZero"/><auto val="1"/>'
        '<lblAlgn val="ctr"/><lblOffset val="100"/><noMultiLvlLbl val="0"/></catAx>'
        '<valAx><axId val="222"/><scaling><logBase val="10"/><orientation val="minMax"/></scaling>'
        '<delete val="0"/><axPos val="l"/><majorGridlines/>'
        f'<title>{_chart_text("buckets, log", 900)}<overlay val="0"/></title>'
        '<numFmt formatCode="General" sourceLinked="0"/><majorTickMark val="out"/><minorTickMark val="none"/>'
        '<tickLblPos val="nextTo"/><crossAx val="111"/><crosses val="autoZero"/>'
        '<crossBetween val="between"/></valAx>'
        '</plotArea><plotVisOnly val="1"/><dispBlanksAs val="gap"/></chart></chartSpace>'
    )


def charts_for_sheet(charts: list[dict]):
    """xlsx_book `charts` callable: the chart dicts of bucket_histogram_charts as
    ready chartSpace parts (the writer needs the sanitised sheet name for the cell refs)."""
    cx, cy = CHART_SIZE_EMU

    def make(sheet_name: str) -> list[dict]:
        return [{"xml": _chart_xml(sheet_name, ch), "col": ch["col"], "row": ch["row"], "cx": cx, "cy": cy}
                for ch in charts]

    return make


def overview_rows() -> list[dict]:
    def row(sheet: str, column: str, meaning: str) -> dict:
        return {"sheet": sheet, "column": column, "meaning": meaning}

    return [
        row("all", "dataset / k / w / L", "dataset id and the unified parameter names shared by all apps: k = minimizer length (sshash's m, fixed at 21 for every order; see README 'Choice of m'), w = window size in minimizer-length k-mers = L - k + 1 (11), L = span of a window = k + w - 1 = sshash's k-mer length (31)"),
        row("all", "order", "canonical minimizer order name: native (sshash's own baseline) | xor | xor-pcg-hash | greedymini | mod-minimizer-xor | mod-minimizer-greedymini (strand-symmetric mod-minimizer). Never carries parameters: a variant is the five columns order, gm_k, projection, gm_table, mod_t"),
        row("all", "gm_k", "GreedyMini table k (the table's gm_k-mer is extended to the m-mer around its centre for greedymini, which needs m - gm_k even; for mod-minimizer-greedymini it is the t-mer table); blank for the other orders"),
        row("all", "projection", "GreedyMini DNA->binary projection polarity (sshash --gm-projection): at0-cg1 = A,T->0 and C,G->1 (default), at1-cg0 = flipped. GreedyMini orders only; blank otherwise. Choices: config.yaml's greedymini / mod_minimizer_greedymini projections (or sshash_bench.sh --gm-projections)"),
        row("all", "gm_table", "stem of the GreedyMini table file under orders/GreedyMini (e.g. w11_k15_canonical); GreedyMini orders only. Choices: config.yaml's tables (or --gm-tables)"),
        row("all", "mod_t", "t-mer length of the mod-minimizer orders (--mod-t): 4 + (m - 4) mod w for mod-minimizer-xor, the table's gm_k for mod-minimizer-greedymini; blank otherwise"),
        row("all", "modality", "canonical (single-probe, default) or regular (two-probe -- sshash_bench.sh -M/--modality, sshash build --regular). Separate axis from order; the app parameter, listed after the order columns (then threads in runs) and before the seed columns."),
        row("all", "n_kmers", "distinct k-mers in the built index; must be identical across every order/seed for a dataset (equivalence gate)"),
        row("summary", "n_seeds", "number of seed repeats aggregated into this row's mean/std"),
        row("summary/runs", "index_size_gb / index_disk_gb / build_max_rss_gb / query_max_rss_gb", "sizes and peak memory in decimal GB (1 GB = 10^9 bytes): index_size_gb = the index size build reports, index_disk_gb = the index file on disk, *_max_rss_gb = peak resident memory (GNU time, KiB converted). The per-k-mer metrics (bits_per_kmer, *_bits_per_kmer) are efficiency ratios and stay in bits"),
        row("summary/runs", "bits_per_kmer", "index_disk_gb * 1e9 * 8 / n_kmers (the on-disk index file size, os.path.getsize)"),
        row("summary/runs", "bench_*_ns", "sshash bench's internal avg_nanosec_per_kmer (the paper's own headline query-time metric)"),
        row("summary/runs", "query_*", "secondary: sshash query against the dataset's real paired FASTQ (wall time via /usr/bin/time -v, elapsed_millisec from sshash's own JSON)"),
        row("summary", "{metric}_ratio", "this variant's mean / the same-modality native's mean for the same dataset (either native baseline = 1), isolating the order effect from the canonical-vs-regular parsing-mode difference; color scale green (better) -> red (worse)"),
        row("runs", "query_pos_match / query_neg_match", "fixed positive/negative k-mer set (generated once per dataset from a native reference build) answered identically by this order/seed's index -- see scripts/equivalence.py; export hard-fails if any is false"),
        row("summary/runs/orders", "w", "window size L - k + 1 (=11 for L=31, k=21; the 'w11' in orders/GreedyMini/w11_k15_canonical)"),
        row("summary/runs/orders", "density_factor", "particular (measured): (n_minimizer_positions/n_kmers) * (w+1) -- see scripts/density.py"),
        row("orders", "k / w / L", "unified names shared by all apps: k = m (minimizer length, 21), w = window size in m-mers = sshash k - m + 1 (11), L = k + w - 1 = sshash's k-mer length (31). One row per (order variant, k, w, L)"),
        row("orders", "expected_density_factor_mean/std", "theoretical (Monte Carlo on i.i.d. random DNA, mean/std over 5 seeds), independent of dataset -- see scripts/expected_density.py"),
        row("buckets/summary/runs", "n_buckets / avg_bucket_size", "n_kmers / n_buckets; a lower-density order has fewer buckets, hence a larger average bucket -- see scripts/export.py:parse_bucket_stats"),
        row("buckets/summary/runs", "kmers_in_skew_index_pct", "share of k-mers whose bucket is large enough to go in sshash's skew index (the structure for heavy buckets); higher can mean slower lookups despite a smaller/denser index"),
        row("bucket_histogram", "bucket_type / log2_bucket_size / bucket_size_range / n_buckets / n_kmers_in_buckets", "one row per run x log2 bucket-size bin, the paper's logarithmic partitioning: bin b holds the sizes (2^(b-1), 2^b] (b = ceil(log2 size); 0 = singleton), the same intervals as the skew-index partitions. bucket_type: singleton (size 1), light (2 .. 2^l) or heavy (> 2^l, i.e. in the skew index; 2^l = 64 here). Bucket size = number of sampled-m-mer occurrences (super-k-mers). Only runs built with `--bucket-stats` (the harness passes it) have rows. sshash does not tally buckets larger than 4096 per size, so one aggregated row `>4096` (log2 >= 13) carries them (derived from the run's totals): the k-mer column then sums exactly to n_kmers and the bucket column to the number of sampled m-mers"),
        row("super_kmer_histogram", "kmers_in_super_kmer / n_super_kmers", "one row per run x super-k-mer length (k-mers per super-k-mer): how many super-k-mers have that length"),
        row("runs", "n_minimizer_positions", "Z: total occurrences of the sampled m-mers (super-k-mer start positions, build.log n_minimizer_positions); density_factor is Z / n_kmers * (w+1)"),
        row("buckets/runs", "minimizer_positions_per_minimizer", "Z / M (n_minimizer_positions / n_buckets): occurrences per distinct minimizer, i.e. how many positions a lookup scans in a light bucket and how much the data repeats a sampled m-mer. Density fixes Z, but M follows density only while this stays about the same"),
        row("buckets/runs", "n_buckets_in_skew_index / minimizer_positions_in_skew_index(_pct)", "heavy buckets (more than 2^l occurrences) and their occurrences, i.e. the entries of H (heavy_load_buckets); build.log num_buckets_in_skew_index / num_minimizer_positions_of_buckets_in_skew_index; the pct is the share of all Z"),
        row("runs", "num_bits_per_offset", "width in bits of one position in S (build.log), which is also the width of an H entry"),
        row("buckets/runs", "skew_h_bits_per_kmer / skew_mphf_v_bits_per_kmer", "split of skew_index_bits_per_kmer: H = minimizer_positions_in_skew_index * num_bits_per_offset / n_kmers (exact, H is a compact vector of that width); the rest is the per-partition MPHFs f_i and the V_i arrays"),
        row("summary", "correlation with density / distinct minimizers", "two rows right under the table: for every {metric}_ratio column, Excel CORREL (Pearson) between that column and density_factor_ratio / n_buckets_ratio over the dataset x order rows (ratios are order / native for the same dataset, so datasets are comparable); blank for a constant column. No colour scale"),
        row("buckets", "{metric}_ratio", "this order's mean / the same-modality native's mean for the same dataset, for every bucket metric (the gain in space for the bits/k-mer breakdown and n_buckets); green (lower) -> red (higher) only on the space and n_buckets ratios, where lower is better"),
        row("buckets", "correlation with distinct minimizers", "one row right under the table: Excel CORREL between n_buckets_ratio and control_codewords_bits_per_kmer_ratio / mid_load_bits_per_kmer_ratio (the tag array and the light array L) over the dataset x order rows"),
        row("space_breakdown", "bits_* / saved_* / share_*_pct", "where the bits of the index go (mean bits/k-mer over the seeds, from the build log's SPACE BREAKDOWN: strings, strings_offsets, mphf, tags = control_codewords, light_L = mid_load_buckets array L, skew_fv = the skew index's MPHFs f_i and V_i arrays, skew_H = H, other = begin_buckets_of_size), the saving of each component over the same-modality native (native - order: positive = saved) and its share of the total saving (%). The weights (3e-6 bits/k-mer) are left out"),
        row("space_breakdown", "saved_by_minimizers / saved_by_occurrences", "the saving split by what sizes the component: by_minimizers = mphf + tags (one entry per distinct minimizer M); by_occurrences = light_L + skew_fv + skew_H (sized by the minimizer occurrences / the heavy tail, so it can be negative); the strings and offsets do not depend on the order"),
        row("buckets/summary/runs", "n_buckets", "(also in the summary sheet: the correlation rows use its ratio) number of distinct sampled m-mers (num_minimizers, M): the MPHF and control-codeword array both have M entries"),
        row("buckets/summary/runs", "singleton_bucket_pct / kmers_in_singleton_bucket_pct", "share of buckets, and of k-mers, whose sampled m-mer occurs in exactly one super-k-mer (the fast singleton path, no OFFSETS access)"),
        row("runs", "heavy_bucket_threshold / expected_bucket_cost", "cost model (runs sheet only): the top 1% of k-mers (largest buckets) all cost one capped bucket size (= heavy_bucket_threshold); the rest cost their bucket size (number of super-k-mers to check). expected_bucket_cost = k-mer-weighted mean of that cost; lower is better"),
        row("buckets/summary/runs", "skew_index_bits_per_kmer / mid_load_bits_per_kmer / control_codewords_bits_per_kmer", "space breakdown for the bucket-count/size tradeoff: fewer buckets means less control_codewords/mid_load overhead but a bigger skew_index share"),
    ]


ORDER_DESCRIPTIONS = {
    "native": {
        "description": "sshash's own unmodified baseline order (mixer_64-based)",
        "cli": "(no -O)",
    },
    "xor": {
        "description": "XOR the minimizer with the seed, then rank canonically",
        "cli": "-O xor -s SEED",
    },
    "xor-pcg-hash": {
        "description": "xor's XOR step, then PCG RXS-M-XS 64-bit scramble before ranking",
        "cli": "-O xor-pcg-hash -s SEED",
    },
}


def order_info(v: dict) -> dict:
    """ORDER_DESCRIPTIONS entry for the order of variant `v` (an ID_FIELDS dict); a
    GreedyMini / mod-minimizer variant gets a description generated from its table,
    projection and t."""
    order = v.get("order") or ""
    if order in ORDER_DESCRIPTIONS:
        return ORDER_DESCRIPTIONS[order]
    flipped = v.get("projection") == "at1-cg0"
    table = f"{v.get('gm_table')}.gm" if v.get("gm_table") else None
    proj_cli = " --gm-projection at1-cg0" if flipped else ""
    flip_txt = (", with the DNA->binary projection polarity flipped (AT->1, CG->0 instead of AT->0, CG->1)"
                if flipped else "")
    if on.is_mod_minimizer(order):
        common = ("Strand-symmetric mod-minimizer (mod-sampling, Groot Koerkamp & Pibiri 2024, with a "
                  "per-k-mer preferred strand -- more G/T; this repo, not published): the k-mer's t-mers "
                  "are ranked, the m-mer at x_L mod w is sampled (x_L = leftmost minimal t-mer), or the "
                  "mirror rule on a reverse-preferred k-mer, so a k-mer and its reverse complement share "
                  "the bucket; ")
        if order == "mod-minimizer-xor":
            return {"description": common + "t-mers ranked by xor, t = 4 + (m - 4) mod w",
                    "cli": "-O mod-minimizer-xor --mod-t T -s SEED"}
        return {"description": common + f"t-mers ranked by the GreedyMini {table} table, t = its gm_k" + flip_txt,
                "cli": f"-O mod-minimizer-greedymini --greedy-mini-file orders/GreedyMini/{table} --mod-t T"
                + proj_cli + " -s SEED"}
    if order == "greedymini":
        return {
            "description": (
                f"Ranks minimizers by the GreedyMini {table} table extended to the m-mer around its "
                "centre gm_k-mer (ties by min(Sigma, reverse Sigma) then Gamma)" + flip_txt
            ),
            "cli": f"-O greedymini --greedy-mini-file orders/GreedyMini/{table}" + proj_cli + " -s SEED",
        }
    return {"description": "", "cli": ""}


DATASET_ORDER_PREFERENCE = ("cod", "kestrel", "human", "ncbi-virus", "se", "hprc")


def _dataset_sort_key(dataset: str) -> tuple:
    try:
        return (0, DATASET_ORDER_PREFERENCE.index(dataset))
    except ValueError:
        return (1, dataset)


def summary_datasets(summary: list[dict]) -> list[str]:
    return sorted({r["dataset"] for r in summary if r.get("dataset")}, key=_dataset_sort_key)


def orders_rows(summary: list[dict], expected: dict[tuple, dict] | None = None) -> list[dict]:
    """One row per order variant (`expected` is keyed by variant_key) -- not per
    dataset x variant -- w, expected density, and
    description/CLI don't vary by dataset, so repeating them per dataset is
    pure duplication: the order's description/CLI/w/expected density_factor,
    plus each dataset's particular (measured, mean/std over seeds)
    density_factor pivoted into its own columns (see orders_field_order). Raw
    minimizer_density isn't kept anywhere -- density_factor is enough."""
    expected = expected or {}
    datasets = summary_datasets(summary)
    by_variant: dict[tuple, dict[str, dict]] = {}
    heads: dict[tuple, dict] = {}
    for r in summary:
        if r.get("order"):
            key = variant_key(r)
            heads[key] = r
            per = by_variant.setdefault(key, {})
            # a regular-modality run never replaces the canonical one
            if r.get("dataset") not in per or (r.get("modality") or "canonical") == "canonical":
                per[r.get("dataset")] = r

    out = []
    for key in sorted(by_variant, key=lambda k: variant_sort_key(heads[k])):
        per_dataset = by_variant[key]
        info = order_info(heads[key])
        exp = expected.get(key, {})
        head = next((r for r in per_dataset.values() if r.get("k") is not None and r.get("w") is not None), {})
        deterministic = heads[key].get("order") in BASELINES
        exp_std = _num(exp.get("density_factor_std"))
        row = {
            **{f: heads[key].get(f) for f in on.ID_FIELDS},
            "description": info["description"],
            "cli": info["cli"],
            "expected_density_factor_mean": _num(exp.get("density_factor_mean")),
            "expected_density_factor_std": 0.0 if deterministic and exp_std in (None, "") else exp_std,
        }
        # unified names: k = m (minimizer length), w = k_sshash - m + 1, L = k_sshash
        if head:
            osh.with_klw(row, head["k"], head["w"])
        for ds in datasets:
            r = per_dataset.get(ds)
            row[f"{ds}_particular_density_factor_mean"] = r.get("density_factor_mean") if r else None
            std = r.get("density_factor_std") if r else None
            row[f"{ds}_particular_density_factor_std"] = 0.0 if r and deterministic and std in (None, "") else std
        out.append(row)
    return out


def orders_field_order(summary: list[dict]) -> list[str]:
    return list(osh.orders_fields(summary_datasets(summary)))




def bucket_rows(summary: list[dict]) -> list[dict]:
    """Dedicated bucket/binning view (mean +/- std over seeds), per dataset x
    order -- see scripts/export.py:parse_bucket_stats. Explains why a
    lower-density order isn't automatically faster to query: fewer buckets
    (from lower density) means a larger *average* bucket, and a bigger share
    of k-mers land in sshash's skew index (its structure for heavy buckets).
    These columns are the only place the bucket metrics are shown (the summary
    sheet leaves them out)."""
    out = []
    for r in summary:
        row = {k: r.get(k) for k in ("dataset", *on.ID_FIELDS, "modality", "n_seeds", "n_kmers")}
        for metric in BUCKET_METRICS:
            row[f"{metric}_mean"] = r.get(f"{metric}_mean")
            row[f"{metric}_std"] = r.get(f"{metric}_std")
            row[f"{metric}_ratio"] = r.get(f"{metric}_ratio")
        out.append(row)
    return out


def bucket_field_order() -> list[str]:
    fields = ["dataset", *on.ID_FIELDS, "modality", "n_seeds", "n_kmers"]
    for metric in BUCKET_METRICS:
        fields.extend((f"{metric}_mean", f"{metric}_std", f"{metric}_ratio"))
    return fields


# Space breakdown: where the bits of the index go, and which part of the saving of an order
# over native comes from what. Components are the dictionary::print_info entries of the
# build log (bits per k-mer); the weights (3e-6) are left out.
SPACE_COMPONENTS = (
    ("strings", "strings_bits_per_kmer"),
    ("strings_offsets", "strings_offsets_bits_per_kmer"),
    ("mphf", "mphf_bits_per_kmer"),
    ("tags", "control_codewords_bits_per_kmer"),
    ("light_L", "mid_load_bits_per_kmer"),
    ("skew_fv", "skew_mphf_v_bits_per_kmer"),
    ("skew_H", "skew_h_bits_per_kmer"),
    ("other", "begin_buckets_bits_per_kmer"),
)
SPACE_BY_MINIMIZERS = ("mphf", "tags")            # one entry per distinct minimizer (M)
SPACE_BY_OCCURRENCES = ("light_L", "skew_fv", "skew_H")  # sized by minimizer occurrences / heavy tail
SPACE_ID_FIELDS = ("dataset", *on.ID_FIELDS, "modality", "n_seeds", "n_kmers")


def space_rows(rows: list[dict], summary: list[dict], baseline: str = "native") -> list[dict]:
    """Per dataset x variant: mean bits/k-mer of every component, the saving over the
    same-modality native (native - order, so positive = saved), each component's share
    of the total saving, and the saving split by what sizes the component (see
    SPACE_BY_MINIMIZERS / SPACE_BY_OCCURRENCES). Values, not formulas: the means are
    over the seeds' runs."""
    runs_by_key: dict[tuple, list[dict]] = {}
    for r in rows:
        runs_by_key.setdefault(summary_group_key(r), []).append(r)
    means: dict[tuple, dict] = {}
    for srow in summary:
        parts = runs_by_key.get(summary_group_key(srow), [])
        comp = {}
        for name, field in SPACE_COMPONENTS:
            vals = [float(p[field]) for p in parts if p.get(field) is not None]
            comp[name] = (sum(vals) / len(vals)) if vals else None
        means[(srow["dataset"], srow.get("modality"), variant_key(srow))] = comp

    def total(comp: dict) -> float | None:
        vals = [v for v in comp.values() if v is not None]
        return sum(vals) if len(vals) == len(comp) else None

    out = []
    for srow in summary:
        key = (srow["dataset"], srow.get("modality"), variant_key(srow))
        comp = means[key]
        base_id = variant_key(on.make_id(baseline, baselines=BASELINES))
        base = means.get((srow["dataset"], srow.get("modality"), base_id)) or means.get(
            (srow["dataset"], "canonical", base_id)
        )
        row: dict = {k: srow.get(k) for k in SPACE_ID_FIELDS}
        for name, _ in SPACE_COMPONENTS:
            row[f"bits_{name}"] = round(comp[name], 4) if comp[name] is not None else None
        tot, btot = total(comp), (total(base) if base else None)
        row["bits_total"] = round(tot, 4) if tot is not None else None
        saved = {}
        for name, _ in SPACE_COMPONENTS:
            if base and base[name] is not None and comp[name] is not None:
                saved[name] = base[name] - comp[name]
            else:
                saved[name] = None
        saved_total = (btot - tot) if (tot is not None and btot is not None) else None
        for name, _ in SPACE_COMPONENTS:
            row[f"saved_{name}"] = round(saved[name], 4) if saved[name] is not None else None
        row["saved_total"] = round(saved_total, 4) if saved_total is not None else None
        has_saving = saved_total is not None and abs(saved_total) > 1e-9
        for name, _ in SPACE_COMPONENTS:
            row[f"share_{name}_pct"] = (
                round(100.0 * saved[name] / saved_total, 2)
                if has_saving and saved[name] is not None else None
            )
        by_m = sum(saved[n] for n in SPACE_BY_MINIMIZERS) if all(saved[n] is not None for n in SPACE_BY_MINIMIZERS) else None
        by_occ = sum(saved[n] for n in SPACE_BY_OCCURRENCES) if all(saved[n] is not None for n in SPACE_BY_OCCURRENCES) else None
        row["saved_by_minimizers"] = round(by_m, 4) if by_m is not None else None
        row["saved_by_occurrences"] = round(by_occ, 4) if by_occ is not None else None
        row["share_by_minimizers_pct"] = round(100.0 * by_m / saved_total, 2) if has_saving and by_m is not None else None
        row["share_by_occurrences_pct"] = round(100.0 * by_occ / saved_total, 2) if has_saving and by_occ is not None else None
        out.append(row)
    return out


def space_field_order() -> list[str]:
    names = [n for n, _ in SPACE_COMPONENTS]
    return [
        *SPACE_ID_FIELDS,
        *(f"bits_{n}" for n in names), "bits_total",
        *(f"saved_{n}" for n in names), "saved_total",
        *(f"share_{n}_pct" for n in names),
        "saved_by_minimizers", "saved_by_occurrences",
        "share_by_minimizers_pct", "share_by_occurrences_pct",
    ]


def _space_field_category(field: str) -> str:
    if field in SPACE_ID_FIELDS:
        return "id / params"
    if field.startswith(("saved_by_", "share_by_")):
        return "what sizes the saving (M vs occurrences)"
    if field.startswith("bits_"):
        return "bits per k-mer"
    if field.startswith("saved_"):
        return "saved vs native (bits per k-mer)"
    if field.startswith("share_"):
        return "share of the total saving (%)"
    return ""


def datasets_rows(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    seen: set[str] = set()
    for r in rows:
        ds = r.get("dataset")
        if not ds or ds in seen:
            continue
        seen.add(ds)
        out.append({"dataset": ds, "k": r.get("k"), "w": r.get("w"), "L": r.get("L"), "n_kmers": r.get("n_kmers")})
    return sorted(out, key=lambda r: r["dataset"])


# ---------------------------------------------------------------------------
# sheets (xlsx_book sheet dicts)
# ---------------------------------------------------------------------------

BASELINE_ORDER = "native"  # the ratio baseline: every `_ratio` is this variant / native

_WRAP_WIDTHS = {"description": 70, "cli": 60, "meaning": 110, "column": 40, "sheet": 22}


def _sheet(name: str, fields: Sequence[str], rows: list[dict], **extra) -> dict:
    fields = list(fields)
    sh = {"name": name, "fields": fields, "rows": rows,
          "categories": sheet_categories(name, fields) or [""] * len(fields)}
    sh.update(extra)
    return sh


def _dataset_block_sheet(name: str, fields: Sequence[str], rows: list[dict], **extra) -> dict:
    """A sheet with one block per dataset (border after it, dataset/k/m/w/modality merged
    inside it); the variant columns merge across the seeds of one variant (runs sheet)."""
    fields = list(fields)
    merge = [f for f in XLSX_MERGE_FIELDS if f in fields]
    scope = {}
    if name == "runs":
        scope = {f: ["dataset", "modality", *on.ID_FIELDS] for f in on.ID_FIELDS if f in fields}
        merge += list(scope)
    return _sheet(name, fields, rows, merge_fields=merge, merge_scope=scope,
                  block_fields=["dataset"], **extra)


def build_sheets(rows, summary, orders, buckets, space, bucket_hist, super_kmer_hist) -> list[dict]:
    """The workbook: overview, orders, datasets, summary, runs (+ buckets, space_breakdown,
    bucket_histogram, super_kmer_histogram)."""
    summary_fields = summary_field_order()
    bucket_fields = bucket_field_order()
    ratio = {f: "lower" for f in summary_fields if f.endswith("_ratio")}
    bucket_ratio = {f: "lower" for f in bucket_fields
                    if f.endswith("_ratio") and f[: -len("_ratio")] in BUCKET_COLORED_RATIOS}
    return [
        _sheet("overview", ["sheet", "column", "meaning"], overview_rows(),
               wrap_fields=["meaning"], widths=_WRAP_WIDTHS),
        _sheet("orders", orders_field_order(summary), orders, wrap_fields=["description", "cli"],
               widths=_WRAP_WIDTHS),
        _sheet("datasets", ["dataset", "k", "w", "L", "n_kmers"], datasets_rows(rows)),
        _dataset_block_sheet("summary", summary_fields, summary, ratio_fields=ratio,
                             ratio_baseline=1.0, footer_rows=correlation_rows(summary_fields, summary)),
        _dataset_block_sheet("runs", list(RUN_FIELDS), rows),
        _dataset_block_sheet("buckets", bucket_fields, buckets, ratio_fields=bucket_ratio, ratio_baseline=1.0,
                             footer_rows=correlation_rows(bucket_fields, buckets, BUCKET_CORRELATION_DRIVERS,
                                                          BUCKET_CORRELATION_COLUMNS)),
        _dataset_block_sheet("space_breakdown", space_field_order(), space),
        _sheet("bucket_histogram", BUCKET_HISTOGRAM_FIELDS, bucket_hist,
               charts=charts_for_sheet(bucket_histogram_charts(bucket_hist))),
        _sheet("super_kmer_histogram", SUPER_KMER_HISTOGRAM_FIELDS, super_kmer_hist),
    ]


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def cmd_evaluate(args: argparse.Namespace) -> None:
    out_root = Path(args.results_dir).resolve()
    runs_root = out_root / "runs"
    csvs_dir = out_root / "csvs"
    xlsx_dir = out_root  # workbook sits directly in results/
    csvs_dir.mkdir(parents=True, exist_ok=True)
    xlsx_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    bucket_hist: list[dict] = []
    super_kmer_hist: list[dict] = []
    for run_dir in iter_run_dirs(runs_root):
        row = load_run_row(run_dir, runs_root)
        if row and row.get("dataset") not in EXCLUDED_DATASETS:
            rows.append(row)
            hist = parse_bucket_histograms(run_dir / "build.log")
            ident = {f: row.get(f) for f in HISTOGRAM_ID_FIELDS}
            over = histogram_overflow(hist["bucket_sizes"], row.get("n_kmers"), row.get("n_buckets"))
            bucket_hist.extend(
                {**ident, **b}
                for b in log2_bucket_bins(hist["bucket_sizes"], over, hist["heavy_lower"])
            )
            super_kmer_hist.extend(
                {**ident, "kmers_in_super_kmer": l, "n_super_kmers": c}
                for l, c in hist["super_kmer_lengths"]
            )
    if not rows:
        die(f"no completed runs in {runs_root}")

    n_cmp, errors = validate_orders_equivalence(rows)
    if errors:
        for msg in errors:
            print(f"error: {msg}", file=sys.stderr)
        die(f"{len(errors)} equivalence/check failure(s) out of {n_cmp} run(s); see above")

    rows.sort(key=lambda r: (*row_sort_key(r), int(r.get("seed_idx") or 0)))
    summary = summary_rows(rows)
    attach_native_ratios(summary, BASELINE_ORDER)

    m_val = next((r.get("k") for r in summary if r.get("k") is not None), None)  # minimizer length
    w_val = next((r.get("w") for r in summary if r.get("w") is not None), None)
    expected: dict[tuple, dict] = {}
    if m_val is not None and w_val is not None:
        seen: dict[tuple, dict] = {}
        for r in summary:
            if r.get("order"):
                seen.setdefault(variant_key(r), {**{f: r.get(f) for f in on.ID_FIELDS},
                                                 "gm_file": r.get("gm_table") and gm_table_file(r["gm_table"])})
        expected = ensure_expected_orders(
            out_root, int(m_val), int(w_val), list(bench_seeds())[:5], variants=list(seen.values()),
        )
    orders = orders_rows(summary, expected)
    required = [f for f in orders_field_order(summary)
                if f in ("order", "k", "w", "L") or "density_factor" in f]
    holes = osh.missing_cells(orders, required)
    if holes:
        shown = ", ".join(f"{orders[i].get('order')}/{orders[i].get('gm_table') or ''}:{f}" for i, f in holes[:12])
        die(f"orders sheet has {len(holes)} empty density/identity cell(s): {shown}")
    buckets = bucket_rows(summary)
    space = space_rows(rows, summary, BASELINE_ORDER)

    write_csv(csvs_dir / "sshash_runs.csv", rows, list(RUN_FIELDS))
    write_csv(csvs_dir / "sshash_orders.csv", orders, orders_field_order(summary))
    write_csv(csvs_dir / "sshash_summary.csv", summary, summary_field_order())
    write_csv(csvs_dir / "sshash_buckets.csv", buckets, bucket_field_order())
    write_csv(csvs_dir / "sshash_space_breakdown.csv", space, space_field_order())
    write_csv(csvs_dir / "sshash_bucket_histogram.csv", bucket_hist, BUCKET_HISTOGRAM_FIELDS)
    write_csv(csvs_dir / "sshash_super_kmer_histogram.csv", super_kmer_hist,
              SUPER_KMER_HISTOGRAM_FIELDS)

    sheets = build_sheets(rows, summary, orders, buckets, space, bucket_hist, super_kmer_hist)
    names = write_book(xlsx_dir / "sshash_results.xlsx", sheets)
    print(f"wrote {xlsx_dir / 'sshash_results.xlsx'} ({len(names)} sheets)")

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
