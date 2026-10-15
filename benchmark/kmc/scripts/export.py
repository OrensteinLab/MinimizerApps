#!/usr/bin/env python3
"""Export KMC3 benchmark runs to CSV and XLSX (stdlib only).

Reads ``<results>/runs`` and writes CSVs to ``<results>/csvs`` and the workbook to
``<results>`` (``--results-dir`` default ``$KMC3_RESULTS_DIR`` or ``$KMC3_BENCH_ROOT/results``):

- ``kmc3_runs.csv`` — one row per run (metrics + telemetry)
- ``kmc3_summary.csv`` — one row per run (per seed; headline metrics)
- ``kmc3_vs_paper.csv`` — lex runs vs ``paper_kmc3.csv`` (when paper CSV exists)
- ``kmc3_orders.csv`` — the ``orders`` sheet (one row per order variant, k, w, L)
- ``kmc3_results.xlsx`` — ``overview`` (column reference), ``orders``, ``datasets``,
  ``summary`` (mean/std), ``readlen_check``, ``gm_lut``, ``runs``, ``vs_paper``
  (written with benchmark/common/xlsx_book.write_book)

Every sheet/CSV uses the unified parameter names k = minimizer length (kmc -p),
w = window, L = k + w - 1 = kmc's k-mer length (kmc -k); internally this module keeps
kmc's own k/p and renames on output (``unify_rows``).

Orders are the canonical names of benchmark/common/order_names.py (``xor``,
``xor-pcg-hash``, ``greedymini``) plus kmc's baselines ``default`` / ``lex``; the
GreedyMini identity is the separate columns ``gm_k`` (= minimizer length k of the
table), ``projection``, ``gm_table`` and ``gm_lut`` (dna | binary LUT backend).
- ``kmc3_density_particular.csv`` — per-dataset particular DNA vs binary factors
  (from ``scripts density``; GM-DNA rows get ``binary_density_factor`` for bin/DNA comparison)

Exits with an error if custom-order k-mer counts differ from the lex baseline.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import statistics
import sys
from collections.abc import Sequence
from pathlib import Path

import order_names
from paths import bench_paths
import orders_sheet
from units import gb_from_blocks, gb_from_bytes, gb_from_kib
from scripts.common import die, load_dotenv
from scripts.config import data_dir
from correlation import correlation_rows
from xlsx_book import DATA_START_ROW as XLSX_DATA_START_ROW, write_book
from xlsx_util import col_letter, write_csv as _shared_write_csv
from scripts.density import (
    attach_density_refs,
    load_particular_rows,
    particular_csv_path,
    particular_index,
    particular_variant_index,
)
from scripts.expected_density import expected_orders_csv_path, load_expected_orders
from scripts.metrics import enrich_super_kmer_metrics
from scripts.schemes import (
    BASELINES,
    GM_ID_FIELDS,
    ID_FIELDS,
    SCHEME_ORDER,
    order_from_json,
    resolve_scheme_label,
    variant_key,
    variant_sort_key,
)

# Ratios are taken against kmc's shipped lex order, the `default` rows.
BASELINE_ORDER = "default"


def results_dir(path: str) -> Path:
    return Path(path).resolve()


def parse_time_file(path: Path) -> dict:
    """Parse GNU /usr/bin/time -v into elapsed, RSS, CPU, and I/O fields."""
    out: dict = {
        "time_elapsed_s": None,
        "max_rss_kb": None,
        "user_s": None,
        "sys_s": None,
        "cpu_pct": None,
        "major_pagefaults": None,
        "minor_pagefaults": None,
        "fs_inputs": None,
        "fs_outputs": None,
        "ctx_voluntary": None,
        "ctx_involuntary": None,
    }
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
            out["user_s"] = _parse_trailing_float(s)
        elif s.startswith("System time (seconds)"):
            out["sys_s"] = _parse_trailing_float(s)
        elif "Percent of CPU this job got" in s:
            m = re.search(r"(\d+(?:\.\d+)?)%", s)
            if m:
                out["cpu_pct"] = float(m.group(1))
        elif "Major (requiring I/O) page faults" in s:
            out["major_pagefaults"] = _parse_trailing_float(s)
        elif "Minor (reclaiming a frame) page faults" in s:
            out["minor_pagefaults"] = _parse_trailing_float(s)
        elif s.startswith("File system inputs"):
            out["fs_inputs"] = _parse_trailing_float(s)
        elif s.startswith("File system outputs"):
            out["fs_outputs"] = _parse_trailing_float(s)
        elif s.startswith("Voluntary context switches"):
            out["ctx_voluntary"] = _parse_trailing_float(s)
        elif s.startswith("Involuntary context switches"):
            out["ctx_involuntary"] = _parse_trailing_float(s)
    return out


def _parse_trailing_float(line: str) -> float | None:
    m = re.search(r"([0-9.]+)\s*$", line)
    return float(m.group(1)) if m else None


def parse_kmc_json(path: Path) -> dict:
    raw = json.loads(path.read_text())
    total = raw.get("Total", "")
    m = re.match(r"([0-9.]+)s", str(total))
    total_s = float(m.group(1)) if m else None
    tmp = raw.get("Tmp_size", "")
    m = re.match(r"([0-9.]+)MB", str(tmp))
    # KMC writes tmpSize / 1e6 as "MB"
    tmp_bytes = float(m.group(1)) * 1_000_000 if m else None
    tmp_sm = raw.get("Tmp_size_strict_memory", "")
    m = re.match(r"([0-9.]+)MB", str(tmp_sm))
    tmp_sm_bytes = float(m.group(1)) * 1_000_000 if m else None
    stats = raw.get("Stats") or {}
    return {
        "kmc_total_s": total_s,
        "tmp_bytes": tmp_bytes,
        "tmp_strict_mem_bytes": tmp_sm_bytes,
        "n_below_ci": stats.get("#k-mers_below_min_threshold"),
        "n_above_cx": stats.get("#k-mers_above_max_threshold"),
        "n_unique": stats.get("#Unique_k-mers"),
        "n_unique_counted": stats.get("#Unique_counted_k-mers"),
        "n_total_kmers": stats.get("#Total no. of k-mers"),
        "n_reads": stats.get("#Total_reads"),
        "n_super_kmers": stats.get("#Total_super-k-mers"),
        "stage1_s": _parse_sec(raw.get("1st_stage")),
        "stage2_s": _parse_sec(raw.get("2nd_stage")),
        "stage3_s": _parse_sec(raw.get("3rd_stage")),
        "json_signature_order": raw.get("signature_order"),
        "json_signature_order_label": raw.get("signature_order_label"),
        "json_signature_order_mode": raw.get("signature_order_mode"),
        "json_signature_order_gm_lut": raw.get("signature_order_gm_lut"),
        "json_signature_order_seed": raw.get("signature_order_seed"),
    }


def _parse_sec(value: str | None) -> float | None:
    if not value:
        return None
    m = re.match(r"([0-9.]+)s", str(value))
    return float(m.group(1)) if m else None

COUNT_FIELDS = (
    "n_unique",
    "n_unique_counted",
    "n_total_kmers",
    "n_below_ci",
    "n_above_cx",
    "n_reads",
)

SUMMARY_XLSX_METRIC_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("super-k-mer", ("n_super_kmers", "avg_super_kmer_len", "density_factor")),
    ("tmp_disk_gb", ("tmp_disk_gb",)),
    ("runtime", ("count_wall_s", "count_cpu_s", "count_stage1_s", "count_stage2_s", "count_stage3_s")),
    ("memory", ("count_max_rss_gb",)),
)
DENSITY_REF_FIELDS = ("binary_density_factor",)
# Summary XLSX: derived from counts; Excel formulas overwrite on export for live recalc.
XLSX_FORMULA_MEAN_METRICS = frozenset({"avg_super_kmer_len", "density_factor"})
XLSX_FORMULA_FIELDS = frozenset(
    {"w", "avg_read_len", "avg_super_kmer_len_mean", "density_factor_mean"}
)
RATIO_COLOR_SCALE_METRICS = frozenset(
    {"n_super_kmers", "density_factor", "tmp_disk_gb", "count_wall_s", "count_cpu_s", "count_stage1_s", "count_stage2_s", "count_stage3_s"}
)
XLSX_RATIO_COLOR_SHEETS = frozenset({"summary", "readlen_check", "gm_lut", "runs"})
SUMMARY_ALL_AGG_METRICS = tuple(m for _c, ms in SUMMARY_XLSX_METRIC_GROUPS for m in ms)
RATIO_METRICS = SUMMARY_ALL_AGG_METRICS
def _with_ratio(fields: Sequence[str], metrics: Sequence[str] = RATIO_METRICS) -> tuple[str, ...]:
    want = set(metrics)
    out: list[str] = []
    for field in fields:
        out.append(field)
        if field in want:
            out.append(f"{field}_ratio")
    return tuple(out)


def _fields_from_groups(groups: Sequence[tuple[str, Sequence[str]]]) -> list[str]:
    return [field for _cat, fields in groups for field in fields]


SUMMARY_PARAM_FIELDS = (
    "dataset",
    "n_total_kmers",
    "avg_read_len",
    "n_reads",
    "input_format",
    "p",
    "w",
    "k",
    *ID_FIELDS,
    "n_bins",
    "threads",
    "strict_memory",
    "mem_limit_gb",
    "n_repeats",
)
RUNS_FIELD_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (
        # Matches SUMMARY_PARAM_FIELDS' order for every shared field (dataset ..
        # order, gm_k, projection, gm_table, gm_lut, then the app params n_bins, threads,
        # strict_memory, mem_limit_gb), so "runs" and "summary" line up
        # column-for-column as far as their differing row granularity allows;
        # runs-only extras (run_id, label, signature_order_label, seed_idx,
        # seed) are placed per the rule: dataset, order, params, seeds.
        "id / params",
        (
            "run_id",
            "dataset",
            "n_total_kmers",
            "avg_read_len",
            "n_reads",
            "input_format",
            "label",
            "p",
            "w",
            "k",
            *ID_FIELDS,
            "signature_order_label",
            "n_bins",
            "threads",
            "strict_memory",
            "mem_limit_gb",
            "seed_idx",
            "seed",
        ),
    ),
    (
        "super-k-mer",
        _with_ratio(
            (
                "n_super_kmers",
                "avg_super_kmer_len",
                "density_factor",
                "binary_density_factor",
                "n_unique_counted",
                "n_unique",
                "n_below_ci",
                "n_above_cx",
            )
        ),
    ),
    ("tmp_disk_gb", _with_ratio(("tmp_disk_gb", "tmp_disk_strict_mem_gb"))),
    ("runtime", _with_ratio(("count_wall_s", "count_cpu_s", "count_stage1_s", "count_stage2_s", "count_stage3_s"))),
    ("memory", _with_ratio(("count_max_rss_gb",))),
    (
        "order",
        (
            "signature_order",
            "signature_order_mode",
            "custom_order_file",
        ),
    ),
    (
        "telemetry",
        (
            "user",
            "sys",
            "cpu_pct",
            "major_pagefaults",
            "minor_pagefaults",
            "fs_inputs",
            "read_gb",
            "fs_outputs",
            "count_pinned_cpus",
            "ctx_voluntary",
            "ctx_involuntary",
            "harness_wall_s",
            "kmc_version",
            "host",
            "exit_code",
        ),
    ),
)
VS_PAPER_FIELD_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (
        "id / params",
        (
            "dataset",
            "k",
            "input_format",
            "n_bins",
            "threads",
            "strict_memory",
            "mem_limit_gb",
            "paper_platform",
        ),
    ),
    (
        "runtime",
        ("n_repeats", "count_wall_s", "count_wall_s_std", "paper_count_wall_s", "count_wall_s_ratio", "count_stage1_s", "count_stage2_s"),
    ),
    ("memory", ("count_max_rss_gb", "paper_count_max_rss_gb")),
    ("tmp_disk_gb", ("tmp_disk_gb", "paper_disk_gb")),
    ("super-k-mer", ("n_unique_counted",)),
)
_SUMMARY_CSV_METRICS = (
    "n_super_kmers",
    "avg_super_kmer_len",
    "density_factor",
    "tmp_disk_gb",
    "count_wall_s",
    "count_cpu_s",
    "count_stage1_s",
    "count_stage2_s",
    "count_stage3_s",
    "count_max_rss_gb",
)
CATEGORY_LABELS: dict[str, str] = {
    "id / params": "params",
    "super-k-mer": "super-k-mer",
    "tmp_disk_gb": "disk (GB)",
    "runtime": "runtime (s)",
    "memory": "memory (GB)",
    "binary / DNA": "bin/DNA",
    "layout": "layout",
    "runtime (s)": "runtime (s)",
    "k-mer counts": "k-mers",
    "disk (GB)": "disk (GB)",
    "order": "order",
    "telemetry": "telemetry",
}


def order_of(row: dict) -> str:
    """Order of a run: canonical name or baseline ("default" = KMC's shipped
    is_allowed-protected lex, "lex" = true canonical lex with no denylist --
    see the orders sheet)."""
    return resolve_scheme_label(row)


def scheme_sort_key(row: dict) -> tuple:
    """Same primary grouping as summary_sort_key (dataset -> bins -> threads ->
    mem -> strict_mem -> scheme), so "runs" rows appear in the same per-config
    blocks, in the same order, as their aggregated "summary" row -- k/input_format
    disambiguate within a config (summary_sort_key doesn't need them: one dataset
    is normally run at a single k/format per sweep), seed_idx orders repeats
    within a scheme (summary has no seed dimension, it's already aggregated)."""
    variant = variant_sort_key(row)
    seed = row.get("seed")
    return (
        row.get("dataset", ""),
        int(row.get("n_bins") or 0),
        int(row.get("threads") or 0),
        int(row.get("mem_limit_gb") or 0),
        int(bool(row.get("strict_memory"))),
        int(row.get("k") or 0),
        row.get("input_format", ""),
        variant,
        int(row.get("seed_idx") or 0),
        seed or 0,
    )


def decorate_export_row(row: dict) -> dict:
    out = dict(row)
    out["order"] = order_of(row)
    for field in GM_ID_FIELDS:
        out[field] = row.get(field) if row.get(field) not in ("",) else None
    out["p"] = row.get("signature_len")
    out["w"] = row.get("minimizer_window_w")
    return out


def sort_export_rows(rows: list[dict]) -> list[dict]:
    return sorted((decorate_export_row(r) for r in rows), key=scheme_sort_key)


def summary_sort_key(row: dict) -> tuple:
    return (
        row.get("dataset", ""),
        int(row.get("n_bins") or 0),
        int(row.get("threads") or 0),
        int(row.get("mem_limit_gb") or 0),
        int(row.get("strict_memory") or 0),
        variant_sort_key(row),
    )


def config_columns(r: dict) -> dict:
    return {
        "n_bins": r.get("n_bins"),
        "threads": r.get("threads"),
        "strict_memory": r.get("strict_memory"),
        "mem_limit_gb": r.get("mem_limit_gb"),
    }


def config_group_key(r: dict) -> tuple:
    return (r["dataset"], int(r["k"]), r["input_format"], r.get("config_tag") or "")


def agg_group_key(r: dict) -> tuple:
    """Config + the full order identity (order, gm_k, projection, gm_table, gm_lut),
    so e.g. the two projections or LUT backends of greedymini never pool."""
    return (*config_group_key(r), *(str(x) if x is not None else "" for x in variant_key(r)))


def _mean_std(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return statistics.mean(values), statistics.stdev(values)


def _derive_summary_count_metrics(row: dict) -> None:
    """Derive w, avg super-k-mer length, and density factor from aggregated counts."""
    try:
        k = int(row.get("k") or 0)
        p = int(row.get("p") or 0)
    except (TypeError, ValueError):
        return
    if k <= 0:
        return
    w = k - p + 1
    row["w"] = w
    ns = row.get("n_super_kmers_mean")
    nt = row.get("n_total_kmers")
    if ns is None or nt is None:
        return
    try:
        ns_f, nt_f = float(ns), float(nt)
    except (TypeError, ValueError):
        return
    if ns_f > 0:
        row["avg_super_kmer_len_mean"] = round(nt_f / ns_f, 4)
    if nt_f > 0:
        row["density_factor_mean"] = round(ns_f / nt_f * (w + 1), 4)


def summary_xlsx_rows(rows: list[dict]) -> list[dict]:
    """Aggregate per dataset/config/scheme: mean and std across seeds (xlsx summary sheet)."""
    buckets: dict[tuple, list[dict]] = {}
    for r in rows:
        buckets.setdefault(agg_group_key(r), []).append(r)

    out: list[dict] = []
    for key in sorted(buckets):
        parts = buckets[key]
        head = parts[0]
        row: dict = {
            "dataset": head["dataset"],
            "input_format": head["input_format"],
            "k": head["k"],
            "p": head.get("p") if head.get("p") is not None else head.get("signature_len"),
            "w": None,
            **config_columns(head),
            **{f: head.get(f) for f in ID_FIELDS},
            "n_repeats": len({p.get("seed_idx") for p in parts}),
            "n_reads": head.get("n_reads"),
            "avg_read_len": head.get("avg_read_len"),
            "binary_density_factor": head.get("binary_density_factor"),
        }
        nt_vals = [float(p["n_total_kmers"]) for p in parts if p.get("n_total_kmers") is not None]
        nt_mean, _ = _mean_std(nt_vals)
        row["n_total_kmers"] = round(nt_mean, 4) if nt_mean is not None else None
        for _cat, metrics in SUMMARY_XLSX_METRIC_GROUPS:
            for metric in metrics:
                vals = [float(p[metric]) for p in parts if p.get(metric) is not None]
                if metric in XLSX_FORMULA_MEAN_METRICS:
                    _, std = _mean_std(vals)
                    row[f"{metric}_std"] = round(std, 4) if std is not None else None
                else:
                    mean, std = _mean_std(vals)
                    row[f"{metric}_mean"] = round(mean, 4) if mean is not None else None
                    row[f"{metric}_std"] = round(std, 4) if std is not None else None
        _derive_summary_count_metrics(row)
        out.append(row)
    return sorted(out, key=summary_sort_key)


def _xlsx_match_key(row: dict, *, with_seed: bool) -> tuple:
    key = (
        row.get("dataset"),
        int(row.get("k") or 0),
        row.get("input_format"),
        int(row.get("n_bins") or 0),
        int(row.get("threads") or 0),
        int(row.get("mem_limit_gb") or 0),
        int(bool(row.get("strict_memory"))),
    )
    if with_seed:
        return (*key, int(row.get("seed_idx") or 0))
    return key


def _apply_lex_ratio_formulas(
    rows: list[dict],
    fields: Sequence[str],
    *,
    metric_field,
    with_seed: bool,
) -> None:
    """Excel `{metric}_ratio` = this / lex (lex row = 1) for the same config [+ seed]."""
    col = {name: _xlsx_col(i) for i, name in enumerate(fields)}
    lex_at: dict[tuple, int] = {}
    for i, row in enumerate(rows, start=XLSX_DATA_START_ROW):
        if row.get("order") == BASELINE_ORDER:
            lex_at[_xlsx_match_key(row, with_seed=with_seed)] = i
    for i, row in enumerate(rows, start=XLSX_DATA_START_ROW):
        lex_i = lex_at.get(_xlsx_match_key(row, with_seed=with_seed))
        if lex_i is None:
            continue
        idx = i - XLSX_DATA_START_ROW
        is_lex = row.get("order") == BASELINE_ORDER
        for metric in RATIO_METRICS:
            rname = f"{metric}_ratio"
            mcol = col.get(metric_field(metric))
            rcol = col.get(rname)
            if rname not in col or not mcol or not rcol or mcol == rcol:
                continue
            if is_lex:
                rows[idx][rname] = {"f": "=1", "v": 1.0}
                rows[idx][f"__{rname}_value"] = 1.0
            else:
                # Shadow plain-float value (cells hold formulas, not numbers) so the
                # writer's color-scale min/max can be computed in Python
                # (xlsx_book reads row["__<field>_value"]); also the cached value.
                this_val = _as_float(_plain(row.get(metric_field(metric))))
                lex_val = _as_float(_plain(rows[lex_i - XLSX_DATA_START_ROW].get(metric_field(metric))))
                cached = this_val / lex_val if this_val is not None and lex_val else None
                rows[idx][rname] = {"f": f"={mcol}{i}/{mcol}{lex_i}", "v": cached}
                if cached is not None:
                    rows[idx][f"__{rname}_value"] = cached


def _is_blank_key(v) -> bool:
    if isinstance(v, tuple):
        return False
    return v is None or v == ""


def _xlsx_formula_ref_row(
    rows: Sequence[dict],
    merge_fields: frozenset[str],
    field: str,
    row_idx: int,
) -> int:
    """Excel row to reference when a vertically merged column has no cell on this row
    (same merge group -- the setting block, or the dataset for dataset-level cells --
    and equal non-blank value: exactly xlsx_book's merge rule)."""
    if field not in merge_fields:
        return XLSX_DATA_START_ROW + row_idx
    anchor = row_idx
    while anchor > 0:
        group = XLSX_DATASET_MERGE_SCOPE.get(field, XLSX_BLOCK_FIELDS)
        if any(rows[anchor].get(g) != rows[anchor - 1].get(g) for g in group):
            break
        cur = _merge_field_value(rows[anchor], field)
        # xlsx_book never merges blank cells (and merges only equal values)
        if _is_blank_key(cur) or cur != _merge_field_value(rows[anchor - 1], field):
            break
        anchor -= 1
    return XLSX_DATA_START_ROW + anchor


def _plain(v):
    """The value itself, or a formula cell's cached value."""
    if isinstance(v, dict):
        return v.get("v")
    return v


def _apply_derived_metric_formulas(
    rows: list[dict],
    fields: Sequence[str],
    *,
    ns_field: str,
    avg_field: str,
    df_field: str,
    merge_fields: frozenset[str],
) -> None:
    """Excel formulas for w, avg lengths, and density factor."""
    col = {name: _xlsx_col(i) for i, name in enumerate(fields)}
    k_col, p_col, w_col = col.get("k"), col.get("p"), col.get("w")
    ns_col, nt_col = col.get(ns_field), col.get("n_total_kmers")
    nr_col, ar_col = col.get("n_reads"), col.get("avg_read_len")
    avg_col, df_col = col.get(avg_field), col.get(df_field)
    if not all((k_col, p_col, w_col, ns_col, nt_col, avg_col, df_col)):
        return
    for i, _row in enumerate(rows, start=XLSX_DATA_START_ROW):
        idx = i - XLSX_DATA_START_ROW
        kr = _xlsx_formula_ref_row(rows, merge_fields, "k", idx)
        pr = _xlsx_formula_ref_row(rows, merge_fields, "p", idx)
        wr = _xlsx_formula_ref_row(rows, merge_fields, "w", idx)
        ntr = _xlsx_formula_ref_row(rows, merge_fields, "n_total_kmers", idx)
        nrr = _xlsx_formula_ref_row(rows, merge_fields, "n_reads", idx)
        # cached values = the plain numbers the formulas replace
        cached = {f: _plain(rows[idx].get(f)) for f in ("w", avg_field, df_field, "avg_read_len")}
        rows[idx]["w"] = {"f": f"={k_col}{kr}-{p_col}{pr}+1", "v": cached["w"]}
        rows[idx][avg_field] = {"f": f"={nt_col}{ntr}/{ns_col}{i}", "v": cached[avg_field]}
        rows[idx][df_field] = {"f": f"={ns_col}{i}/{nt_col}{ntr}*({w_col}{wr}+1)", "v": cached[df_field]}
        if nr_col and ar_col:
            rows[idx]["avg_read_len"] = {
                "f": f"={nt_col}{ntr}/{nr_col}{nrr}+({k_col}{kr}-1)",
                "v": cached["avg_read_len"],
            }


def apply_summary_derived_formulas(rows: list[dict], fields: Sequence[str]) -> None:
    """Excel formulas for w, avg lengths, density factor, and /lex ratios on summary."""
    _apply_derived_metric_formulas(
        rows,
        fields,
        ns_field="n_super_kmers_mean",
        avg_field="avg_super_kmer_len_mean",
        df_field="density_factor_mean",
        merge_fields=XLSX_MERGE_FIELDS["summary"],
    )
    _apply_lex_ratio_formulas(rows, fields, metric_field=lambda m: f"{m}_mean", with_seed=False)


# Rows under the summary table: per `{metric}_ratio` column, CORREL with these
# driver ratios (shared helper, as the other apps' summary sheets).
SUMMARY_CORRELATION_DRIVERS = (
    ("correlation with density", "density_factor_ratio"),
)


def summary_correlation_rows(rows: Sequence[dict], fields: Sequence[str]) -> list[dict]:
    """xlsx_book footer rows; the cached values are the computed `__<ratio>_value`s."""
    return correlation_rows(fields, rows, SUMMARY_CORRELATION_DRIVERS,
                            value_of=lambda r, f: r.get(f"__{f}_value"))


def apply_runs_xlsx_formulas(rows: list[dict], fields: Sequence[str]) -> None:
    """Excel formulas for w, avg lengths, density factor, and /lex ratios on runs."""
    _apply_derived_metric_formulas(
        rows,
        fields,
        ns_field="n_super_kmers",
        avg_field="avg_super_kmer_len",
        df_field="density_factor",
        merge_fields=XLSX_MERGE_FIELDS["runs"],
    )
    _apply_lex_ratio_formulas(rows, fields, metric_field=lambda m: m, with_seed=True)


def summary_xlsx_field_order() -> list[str]:
    fields = list(SUMMARY_PARAM_FIELDS)
    for _cat, metrics in SUMMARY_XLSX_METRIC_GROUPS:
        for metric in metrics:
            fields.extend((f"{metric}_mean", f"{metric}_std", f"{metric}_ratio"))
            if metric == "density_factor":
                fields.append("binary_density_factor")
    return fields


def runs_field_order() -> list[str]:
    return _fields_from_groups(RUNS_FIELD_GROUPS)


def summary_field_order() -> list[str]:
    return [
        "dataset",
        "n_total_kmers",
        "avg_read_len",
        "input_format",
        "p",
        "w",
        "k",
        *ID_FIELDS,
        "signature_order_label",
        "n_bins",
        "threads",
        "strict_memory",
        "mem_limit_gb",
        *_with_ratio(("n_super_kmers", "avg_super_kmer_len", "density_factor")),
        "binary_density_factor",
        *_with_ratio(("tmp_disk_gb", "count_wall_s", "count_cpu_s", "count_stage1_s", "count_stage2_s", "count_stage3_s", "count_max_rss_gb")),
    ]


def vs_paper_field_order() -> list[str]:
    return _fields_from_groups(VS_PAPER_FIELD_GROUPS)


def summary_rows(rows: list[dict]) -> list[dict]:
    """One row per dataset / config / scheme with headline metrics."""
    out: list[dict] = []
    for r in rows:
        row: dict = {
            "dataset": r["dataset"],
            "n_total_kmers": r.get("n_total_kmers"),
            "avg_read_len": r.get("avg_read_len"),
            "n_reads": r.get("n_reads"),
            "input_format": r["input_format"],
            "k": r["k"],
            **config_columns(r),
            "p": r.get("p") if r.get("p") is not None else r.get("signature_len"),
            "w": r.get("w") if r.get("w") is not None else r.get("minimizer_window_w"),
            **{f: r.get(f) for f in ID_FIELDS},
            "signature_order_label": r.get("signature_order_label", ""),
            "binary_density_factor": r.get("binary_density_factor"),
            "n_unique_counted": r.get("n_unique_counted"),
        }
        for metric in _SUMMARY_CSV_METRICS:
            row[metric] = r.get(metric)
            row[f"{metric}_ratio"] = r.get(f"{metric}_ratio")
        out.append(row)
    return sorted(out, key=summary_sort_key)


def iter_run_dirs(runs_root: Path) -> list[Path]:
    """Find leaf run directories (nested layout under runs/)."""
    if not runs_root.is_dir():
        return []
    found = [
        meta_path.parent
        for meta_path in runs_root.rglob("meta.json")
        if (meta_path.parent / "kmc.json").is_file()
    ]
    return sorted(set(found))


def run_id_from_dir(runs_root: Path, run_dir: Path) -> str:
    try:
        return str(run_dir.relative_to(runs_root))
    except ValueError:
        return run_dir.name


def load_run_row(run_dir: Path, *, runs_root: Path | None = None) -> dict | None:
    meta_path = run_dir / "meta.json"
    json_path = run_dir / "kmc.json"
    if not meta_path.is_file() or not json_path.is_file():
        return None
    meta = json.loads(meta_path.read_text())
    stats = parse_kmc_json(json_path)
    tinfo = parse_time_file(run_dir / "time.txt")
    rss_kb = tinfo["max_rss_kb"]
    wall = stats["kmc_total_s"] if stats["kmc_total_s"] is not None else tinfo["time_elapsed_s"]
    sig_order = meta.get("signature_order") or stats.get("json_signature_order") or "default"
    order_label = meta.get("signature_order_label") or stats.get("json_signature_order_label") or ""
    order_mode = meta.get("signature_order_mode") or stats.get("json_signature_order_mode") or ""
    # canonical order + separate gm columns (meta.json written by scripts/meta.py);
    # a meta without `order` falls back to kmc.json's self-report.
    order = meta.get("order") or order_from_json(sig_order, order_mode, order_label)
    order_names.validate_order(order, BASELINES)
    is_gm = order_names.is_greedymini(order)
    gm_k = meta.get("gm_k") if is_gm else None
    projection = (meta.get("projection") or order_names.DEFAULT_PROJECTION) if is_gm else None
    gm_lut = (meta.get("gm_lut") or stats.get("json_signature_order_gm_lut") or "") if is_gm else ""
    order_seed = meta.get("signature_order_seed")
    if order_seed is None:
        order_seed = stats.get("json_signature_order_seed")
    ds = meta.get("dataset")
    k = meta.get("k")
    kind = meta.get("input_format")
    rid = meta.get("run_id") or (
        run_id_from_dir(runs_root, run_dir) if runs_root is not None else run_dir.name
    )
    row = {
        "run_id": rid,
        "dataset": ds,
        "label": meta.get("label"),
        "k": k,
        "input_format": kind,
        "threads": meta.get("threads"),
        "mem_limit_gb": meta.get("mem_limit_gb"),
        "n_bins": meta.get("n_bins"),
        "strict_memory": meta.get("strict_memory"),
        "config_tag": meta.get("config_tag", ""),
        "signature_len": meta.get("signature_len"),
        "order": order,
        "gm_k": gm_k,
        "projection": projection,
        "gm_table": (meta.get("gm_table") or None) if is_gm else None,
        "gm_lut": gm_lut or None,
        "signature_order": sig_order,
        "signature_order_mode": order_mode,
        "signature_order_label": order_label,
        "seed": order_seed,
        "seed_idx": meta.get("seed_idx"),
        "custom_order_file": meta.get("custom_order_file", ""),
        "kmc_version": meta.get("kmc_version"),
        "host": meta.get("host"),
        "exit_code": meta.get("exit_code"),
        "count_wall_s": wall,
        "kmc_total_s": stats["kmc_total_s"],
        "time_elapsed_s": tinfo["time_elapsed_s"],
        "harness_wall_s": meta.get("python_wall_s"),
        "user": tinfo["user_s"],
        "sys": tinfo["sys_s"],
        "cpu_pct": tinfo["cpu_pct"],
        "count_stage1_s": stats["stage1_s"],
        "count_stage2_s": stats["stage2_s"],
        "count_stage3_s": stats["stage3_s"],
        # CPU seconds (user + system, all threads; /usr/bin/time -v) of the timed kmc call: with wall
        # time it separates "ran longer" from "waited longer" (contention, input reads)
        "count_cpu_s": (round(tinfo["user_s"] + tinfo["sys_s"], 2)
                        if tinfo["user_s"] is not None and tinfo["sys_s"] is not None else None),
        "count_max_rss_gb": gb_from_kib(rss_kb),
        "tmp_disk_gb": gb_from_bytes(stats["tmp_bytes"]),
        "tmp_disk_strict_mem_gb": gb_from_bytes(stats["tmp_strict_mem_bytes"]),
        "major_pagefaults": tinfo["major_pagefaults"],
        "minor_pagefaults": tinfo["minor_pagefaults"],
        "fs_inputs": tinfo["fs_inputs"],
        # file-system reads in decimal GB ('File system inputs' counts 512-byte blocks); ~0 = input served from the page cache
        "read_gb": gb_from_blocks(tinfo["fs_inputs"]),
        "fs_outputs": tinfo["fs_outputs"],
        # CPUs the timed kmc call was taskset to (cpus.txt, written by lib/run.sh); absent in older runs
        "count_pinned_cpus": ((run_dir / "cpus.txt").read_text().strip() or None) if (run_dir / "cpus.txt").is_file() else None,
        "ctx_voluntary": tinfo["ctx_voluntary"],
        "ctx_involuntary": tinfo["ctx_involuntary"],
        "n_unique": stats["n_unique"],
        "n_unique_counted": stats["n_unique_counted"],
        "n_total_kmers": stats["n_total_kmers"],
        "n_reads": stats["n_reads"],
        "n_below_ci": stats["n_below_ci"],
        "n_above_cx": stats["n_above_cx"],
        "n_super_kmers": stats["n_super_kmers"],
    }
    enrich_super_kmer_metrics(row)
    return row


def _ratio(a, b) -> float | None:
    try:
        if a is None or b is None:
            return None
        bf = float(b)
        if bf == 0:
            return None
        return round(float(a) / bf, 4)
    except (TypeError, ValueError):
        return None


def _attach_lex_ratios(
    rows: list[dict],
    *,
    row_key,
    metric_value,
    ratio_name,
) -> None:
    lex_at: dict[tuple, dict] = {}
    for row in rows:
        if row.get("order") == BASELINE_ORDER:
            lex_at[row_key(row)] = row
    for row in rows:
        is_lex = row.get("order") == BASELINE_ORDER
        base = lex_at.get(row_key(row))
        for metric in RATIO_METRICS:
            rname = ratio_name(metric)
            val = metric_value(row, metric)
            if is_lex:
                row[rname] = 1.0 if val is not None else None
            elif base is None:
                row[rname] = None
            else:
                row[rname] = _ratio(val, metric_value(base, metric))


def attach_lex_ratios(rows: list[dict]) -> None:
    """Set `{metric}_ratio` = value / lex for the same config and seed (lex = 1)."""
    _attach_lex_ratios(
        rows,
        row_key=lambda r: (config_group_key(r), int(r.get("seed_idx") or 0)),
        metric_value=lambda r, m: r.get(m),
        ratio_name=lambda m: f"{m}_ratio",
    )


def is_gm_lut_run(row: dict) -> bool:
    """True for GM DNA-vs-binary LUT sweep runs (not the main GreedyMini scheme)."""
    rid = (row.get("run_id") or "").replace("\\", "/")
    if "/gm_lut_" in rid:
        return True
    tag = (row.get("config_tag") or "").strip()
    return tag.startswith("gm_lut_")


def rows_without_gm_lut(rows: Sequence[dict]) -> list[dict]:
    return [r for r in rows if not is_gm_lut_run(r)]


def is_readlen_derived_dataset(dataset: str) -> bool:
    """True for read-length-confound-check derived datasets (see README).

    `readlen-dataset`/`readlen-check` always name their output
    `{source}_shortsub` / `{source}_longsub` (lib/readlen.sh); this also
    covers manually built file-layout variants like `..._shortsub_mf11`.
    """
    return "_shortsub" in dataset or "_longsub" in dataset


def gm_lut_layout_tag(row: dict) -> str:
    tag = (row.get("config_tag") or "").strip()
    for prefix in ("gm_lut_dna_", "gm_lut_binary_"):
        if tag.startswith(prefix):
            return tag[len(prefix) :]
    rid = (row.get("run_id") or "").replace("\\", "/")
    for part in rid.split("/"):
        for prefix in ("gm_lut_dna_", "gm_lut_binary_"):
            if part.startswith(prefix):
                return part[len(prefix) :]
    return tag


def gm_lut_layout_key(row: dict) -> tuple:
    return (
        row.get("dataset"),
        int(row.get("k") or 0),
        row.get("input_format"),
        gm_lut_layout_tag(row),
    )


def attach_summary_lex_ratios(rows: list[dict]) -> None:
    """Set `{metric}_ratio` on summary-format rows = mean / lex mean (same hardware config)."""
    _attach_lex_ratios(
        rows,
        row_key=lambda r: _xlsx_match_key(r, with_seed=False),
        metric_value=lambda r, m: r.get(f"{m}_mean"),
        ratio_name=lambda m: f"{m}_ratio",
    )


def prepare_agg_xlsx_sheet(rows: list[dict]) -> list[str]:
    """Numeric /lex ratios, drop empty columns, then Excel formulas (summary + gm_lut)."""
    attach_summary_lex_ratios(rows)
    fields = prune_empty_fields(rows, summary_xlsx_field_order())
    apply_summary_derived_formulas(rows, fields)
    return fields


def prepare_runs_xlsx_sheet(rows: list[dict]) -> list[str]:
    """Drop empty columns, then Excel formulas for runs."""
    fields = prune_empty_fields(rows, runs_field_order())
    apply_runs_xlsx_formulas(rows, fields)
    return fields


def gm_lut_summary_rows(all_rows: Sequence[dict], main_rows: Sequence[dict]) -> list[dict]:
    """Summary-format sheet: lex + greedymini per gm_lut (dna | binary) and projection, per GM LUT layout (ratio vs lex)."""
    lut_runs = [r for r in all_rows if is_gm_lut_run(r)]
    if not lut_runs:
        return []

    lex_pool: dict[tuple, list[dict]] = {}
    for row in main_rows:
        if row.get("order") == BASELINE_ORDER:
            lex_pool.setdefault(_xlsx_match_key(row, with_seed=False), []).append(row)

    by_layout: dict[tuple, dict[tuple, list[dict]]] = {}
    for row in lut_runs:
        mode = ((row.get("gm_lut") or "dna").strip().lower(), row.get("projection") or "")
        by_layout.setdefault(gm_lut_layout_key(row), {}).setdefault(mode, []).append(row)

    bundle: list[dict] = []
    for layout_key in sorted(by_layout.keys()):
        modes = by_layout[layout_key]
        ref = next((rs[0] for rs in modes.values() if rs), None)
        if ref is None:
            continue
        hw = _xlsx_match_key(ref, with_seed=False)
        if hw in lex_pool:
            bundle.extend(lex_pool[hw])
        for mode in sorted(modes, key=lambda m: (m[0] != "dna", m)):
            bundle.extend(modes[mode])

    if not bundle:
        return []
    return summary_xlsx_rows(bundle)


# Unified parameter names (common/orders_sheet.py): k = minimizer length, w = window,
# L = k + w - 1 = span (kmc's k-mer length). Internally this module keeps kmc's own
# names (k = k-mer length, p = signature length, as in meta.json / the CLI); every
# sheet and CSV is renamed once on the way out. Simultaneous swap: old k -> L, old p -> k.
UNIFIED_PARAM_NAMES = {"k": "L", "p": "k"}


def unify_name(field: str) -> str:
    return UNIFIED_PARAM_NAMES.get(field, field)


def unify_rows(rows: Sequence[dict]) -> list[dict]:
    return [{unify_name(f): v for f, v in r.items()} for r in rows]


def write_csv(path: Path, rows: list[dict], fields: Sequence[str] | None = None, *, unified: bool = False) -> None:
    """`unified=True`: rows/fields already use the unified names (orders)."""
    if not rows:
        path.write_text("")
        return
    fieldnames = list(fields) if fields else list(rows[0].keys())
    if not unified:
        fieldnames = [unify_name(f) for f in fieldnames]
        rows = unify_rows(rows)
    _shared_write_csv(path, fieldnames, rows)
    print(f"wrote {path} ({len(rows)} rows)")



def _xlsx_col(n: int) -> str:
    return col_letter(n + 1)


XLSX_COLOR_SHEETS = frozenset({"summary", "readlen_check", "gm_lut", "runs"})
XLSX_MERGE_DATASET = (
    "dataset",
    "n_total_kmers",
    "n_reads",
    "input_format",
    "label",
)
XLSX_MERGE_MINIMIZER = ("k", "w", "L")  # unified names (make_sheet renames before use)
XLSX_MERGE_CONFIG = ("n_bins", "threads", "strict_memory", "mem_limit_gb")
XLSX_MERGE_FIELDS: dict[str, frozenset[str]] = {
    "summary": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, "n_repeats")
    ),
    "readlen_check": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, "n_repeats")
    ),
    "gm_lut": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, "n_repeats")
    ),
    # runs: the id columns (order, gm_k, projection, gm_table, gm_lut) merge over the
    # seed repeats of one variant -- see XLSX_ID_MERGE_SCOPE
    "runs": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, *ID_FIELDS)
    ),
}
# An id column merges across rows that agree on the whole config + the full
# order identity (so two variants never share a merged id cell).
XLSX_ID_MERGE_SCOPE = {
    f: ("dataset", "L", "input_format", "n_bins", "threads", "mem_limit_gb", "strict_memory", *ID_FIELDS)
    for f in ID_FIELDS
}
# Block = one setting = dataset + input format + config (n_bins, threads, memory): a medium
# bottom border after each; merges never cross it, except the dataset-level cells below.
XLSX_BLOCK_FIELDS = ("dataset", "input_format", "n_bins", "threads", "strict_memory", "mem_limit_gb")
# dataset-level cells span all the settings of one dataset
XLSX_DATASET_MERGE_SCOPE = {
    f: ("dataset",) for f in ("dataset", "n_total_kmers", "n_reads", "label", "avg_read_len")
}
# {metric}_ratio columns that get the green(<1) / white(1) / red(>1) scale: all "lower is better".
RATIO_COLOR_SCALE_FIELD_METRICS = RATIO_COLOR_SCALE_METRICS
# n_bins / threads cell fills (data cells only; headers and category bands are not filled)
_BINS_FILL = {512: "FFCCFFFF", 1024: "FF99FFFF", 1536: "FF66FFFF"}
_THREADS_FILL = {12: "FFFFCCFF", 24: "FFFF99FF"}


def _config_fill(field: str, _row: dict, value):
    try:
        v = int(value)
    except (TypeError, ValueError):
        return None
    if field == "n_bins":
        return _BINS_FILL.get(v)
    if field == "threads":
        return _THREADS_FILL.get(v)
    return None


def _is_empty_cell(v) -> bool:
    return v is None or v == ""


def _merge_value_key(v) -> str | int | float | bool:
    if isinstance(v, dict) and "f" in v:
        return str(v["f"])
    if v is None:
        return ""
    return v


def _merge_field_value(row: dict, field: str):
    """Comparable value for merge detection (handles formula-backed columns)."""
    if field == "avg_read_len":
        return (row.get("n_total_kmers"), row.get("n_reads"), row.get("L"))
    if field == "w":
        return (row.get("L"), row.get("k"))
    return _merge_value_key(row.get(field))





def _agg_metric_base(field: str) -> str | None:
    for suffix in ("_mean", "_std", "_ratio"):
        if field.endswith(suffix):
            return field[: -len(suffix)]
    return None


# Kept even when every row is empty (runs without cpus.txt): the column is part of the documented schema.
XLSX_KEEP_FIELDS = frozenset({"count_pinned_cpus"})


def prune_empty_fields(rows: Sequence[dict], fields: Sequence[str]) -> list[str]:
    """Drop columns whose values are all empty; drop ratio/std when mean/base is empty."""
    if not rows:
        return list(fields)

    def col_empty(field: str) -> bool:
        return all(_is_empty_cell(r.get(field)) for r in rows)

    drop: set[str] = set()
    fields = list(fields)
    for field in fields:
        if field in XLSX_KEEP_FIELDS:
            continue
        base = _agg_metric_base(field)
        if base:
            mean_f = f"{base}_mean"
            if mean_f in fields and col_empty(mean_f):
                if base in XLSX_FORMULA_MEAN_METRICS:
                    continue
                for ff in fields:
                    if ff == mean_f or ff.startswith(f"{base}_"):
                        drop.add(ff)
        elif col_empty(field):
            if field in XLSX_FORMULA_FIELDS:
                continue
            drop.add(field)
            ratio_f = f"{field}_ratio"
            if ratio_f in fields:
                drop.add(ratio_f)
    for field in fields:
        if field not in drop and col_empty(field):
            if field in XLSX_FORMULA_FIELDS or field in XLSX_KEEP_FIELDS:
                continue
            base = _agg_metric_base(field)
            if base in XLSX_FORMULA_MEAN_METRICS:
                continue
            if base and field.endswith(("_ratio", "_std")):
                mean_f = f"{base}_mean"
                if mean_f in fields and not col_empty(mean_f):
                    continue
            drop.add(field)
    return [f for f in fields if f not in drop]


def _summary_field_category(field: str) -> str:
    if field in SUMMARY_PARAM_FIELDS:
        return "id / params"
    if field in DENSITY_REF_FIELDS:
        return "super-k-mer"
    base = _agg_metric_base(field) or field
    for cat, metrics in SUMMARY_XLSX_METRIC_GROUPS:
        if base in metrics:
            return cat
    return ""


def _category_display(category: str) -> str:
    return CATEGORY_LABELS.get(category, category)


def _categories_from_groups(fields: Sequence[str], groups: Sequence[tuple[str, Sequence[str]]]) -> list[str]:
    lookup = {name: cat for cat, names in groups for name in names}
    return [lookup.get(field, "") for field in fields]

def _xlsx_field_categories(raw_name: str, fields: Sequence[str]) -> list[str] | None:
    """Raw category per column (shown through CATEGORY_LABELS in row 1 of the sheet)."""
    if raw_name in ("summary", "readlen_check"):
        return [_summary_field_category(f) for f in fields]
    if raw_name == "runs":
        return _categories_from_groups(fields, RUNS_FIELD_GROUPS)
    if raw_name == "gm_lut":
        return [_summary_field_category(f) for f in fields]
    if raw_name == "vs_paper":
        return _categories_from_groups(fields, VS_PAPER_FIELD_GROUPS)
    if raw_name == "overview":
        return ["metrics reference"] * len(fields)
    if raw_name == "orders":
        return ["order reference"] * len(fields)
    if raw_name == "datasets":
        return ["dataset reference"] * len(fields)
    return None




_WIDE_TEXT_FIELDS = {"meaning": 70, "notes": 80, "description": 80, "column": 34, "sheet": 20}


def make_sheet(
    name: str,
    rows: list[dict],
    fields: Sequence[str],
    *,
    footer_rows: list[dict] | None = None,
) -> dict:
    """xlsx_book sheet spec for one sheet of the workbook.

    The agg sheets (summary / readlen_check / gm_lut) and runs carry formulas
    authored for the data starting at DATA_START_ROW,
    {metric}_ratio colour scales (green < 1 < red: every ratio metric is "lower
    is better"), setting blocks
    (medium border after each dataset+config setting, merges never cross one except the dataset cells) and
    n_bins / threads cell fills. Category bands and headers are not filled."""
    fields = list(fields)
    categories = _xlsx_field_categories(name, fields)
    if name not in ("orders", "overview"):  # those two are authored with the unified names
        fields = [unify_name(f) for f in fields]
        rows = unify_rows(rows)
    sheet: dict = {
        "name": name,
        "fields": fields,
        "rows": rows,
    }
    if categories:
        sheet["categories"] = [_category_display(c) for c in categories]
    wide = {f: _WIDE_TEXT_FIELDS[f] for f in fields if f in _WIDE_TEXT_FIELDS}
    if name in ("overview", "orders", "datasets") and wide:
        sheet["widths"] = wide
        sheet["wrap_fields"] = [f for f in wide if f in ("meaning", "notes", "description")]
    if name in XLSX_COLOR_SHEETS:
        sheet["block_fields"] = list(XLSX_BLOCK_FIELDS)
        sheet["merge_fields"] = sorted(XLSX_MERGE_FIELDS[name])
        scope = {f: list(sc) for f, sc in XLSX_DATASET_MERGE_SCOPE.items()}
        if name == "runs":
            scope.update({f: list(sc) for f, sc in XLSX_ID_MERGE_SCOPE.items()})
        sheet["merge_scope"] = scope
        sheet["merge_key"] = lambda row, field: _merge_field_value(row, field)
        sheet["cell_fill"] = _config_fill
        sheet["ratio_fields"] = {
            f: "lower"
            for f in fields
            if f.endswith("_ratio") and f[: -len("_ratio")] in RATIO_COLOR_SCALE_METRICS
        }
        sheet["ratio_baseline"] = 1.0
        sheet["freeze_cols"] = 0
    if footer_rows:
        sheet["footer_rows"] = footer_rows
    return sheet


def validate_custom_vs_lex(rows: list[dict]) -> tuple[int, list[str]]:
    """Check custom schemes match lex k-mer counts; return (n_comparisons, errors)."""
    by_key: dict[tuple, dict[str, dict]] = {}
    for r in rows:
        label = r.get("signature_order_label") or ""
        slot = "lex" if r.get("order") == BASELINE_ORDER or not label else f"custom:{label}"
        by_key.setdefault(config_group_key(r), {})[slot] = r

    errors: list[str] = []
    n_cmp = 0
    for key in sorted(by_key):
        group = by_key[key]
        lex = group.get("lex")
        for slot in sorted(s for s in group if s.startswith("custom:")):
            n_cmp += 1
            custom = group[slot]
            ds, k, kind, cfg = key
            if lex is None:
                print(
                    f"warning: {ds} k={k} {kind} config={cfg} custom={slot}: no lex baseline",
                    file=sys.stderr,
                )
                continue
            for field in COUNT_FIELDS:
                a, b = lex.get(field), custom.get(field)
                if not (a is not None and b is not None and a == b):
                    errors.append(
                        f"{ds} k={k} {kind} config={cfg} custom={slot}: {field} lex={a} custom={b}"
                    )
    return n_cmp, errors


def vs_paper_rows(rows: list[dict], paper_rows: list[dict], platform: str) -> list[dict]:
    """Lex runs (averaged across repeats, if any) joined to published KMC3 numbers."""
    paper_index = {
        (r["dataset"], int(r["k"]), r["platform"], r["input_format"]): r for r in paper_rows
    }
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        if order_of(r) != BASELINE_ORDER:
            continue
        groups.setdefault(config_group_key(r), []).append(r)

    out: list[dict] = []
    for reps in groups.values():
        head = reps[0]
        key = (head["dataset"], int(head["k"]), platform, head["input_format"])
        p = paper_index.get(key)
        wt_vals = [float(r["count_wall_s"]) for r in reps if r.get("count_wall_s") is not None]
        wt_mean, wt_std = _mean_std(wt_vals)
        s1_vals = [float(r["count_stage1_s"]) for r in reps if r.get("count_stage1_s") is not None]
        s2_vals = [float(r["count_stage2_s"]) for r in reps if r.get("count_stage2_s") is not None]
        rss_vals = [float(r["count_max_rss_gb"]) for r in reps if r.get("count_max_rss_gb") is not None]
        disk_vals = [float(r["tmp_disk_gb"]) for r in reps if r.get("tmp_disk_gb") is not None]
        paper_t = float(p["paper_time_s"]) if p else None
        out.append(
            {
                "dataset": head["dataset"],
                "k": head["k"],
                "input_format": head["input_format"],
                **config_columns(head),
                "paper_platform": platform,
                "n_repeats": len({r.get("seed_idx") for r in reps}),
                "count_wall_s": round(wt_mean, 4) if wt_mean is not None else None,
                "count_wall_s_std": round(wt_std, 4) if wt_std is not None else None,
                "paper_count_wall_s": paper_t,
                "count_wall_s_ratio": _ratio(wt_mean, paper_t),
                "count_stage1_s": round(_mean_std(s1_vals)[0], 4) if s1_vals else None,
                "count_stage2_s": round(_mean_std(s2_vals)[0], 4) if s2_vals else None,
                "count_max_rss_gb": round(_mean_std(rss_vals)[0], 4) if rss_vals else None,
                "paper_count_max_rss_gb": float(p["paper_ram_gb"]) if p else None,
                "tmp_disk_gb": round(_mean_std(disk_vals)[0], 4) if disk_vals else None,
                "paper_disk_gb": float(p["paper_disk_gb"]) if p else None,
                "n_unique_counted": head.get("n_unique_counted"),
            }
        )
    return sorted(out, key=lambda r: (r["dataset"], int(r["k"]), r["input_format"]))


def load_density_context(out_root: Path) -> dict:
    """Load per-dataset particular DNA/binary density cache."""
    particular = particular_index(load_particular_rows(particular_csv_path(out_root)))
    if not particular:
        print(
            f"warning: no {particular_csv_path(out_root).name}; "
            "run `python -m scripts density particular` for per-dataset binary DF",
            file=sys.stderr,
        )
    return particular


def overview_field_order() -> list[str]:
    return ["sheet", "column", "meaning", "notes"]


PINNING_NOTE = (
    "Every timed step is taskset to one NUMA node, one thread per physical core plus 2 spare cores "
    "(benchmark/common/pin_cpus.sh); the CPU list is in count_pinned_cpus"
)

# column -> (meaning, notes). Metric columns are generated below.
_METRIC_MEANING: dict[str, tuple[str, str]] = {
    "n_super_kmers": ("number of super-k-mers KMC's splitter produced", ""),
    "avg_super_kmer_len": ("n_total_kmers / n_super_kmers -- longer is better (fewer, denser super-k-mers)", ""),
    "density_factor": (
        "n_super_kmers / n_total_kmers * (w + 1) -- lower is better (closer-to-theoretical minimizer density)", ""),
    "tmp_disk_gb": ("KMC temp directory size", "decimal GB (bytes / 10^9)"),
    "count_wall_s": ("wall time of the timed kmc call", "seconds"),
    "count_cpu_s": (
        "CPU time of the timed kmc call (user + system, summed over all threads; GNU time -v). Wall time = CPU "
        "time / cores actually running, so wall time rises with contention or input reads while CPU time mostly "
        "does not: compare orders on both -- a wall-time difference without a CPU-time one is load or I/O, "
        "not the order", "seconds"),
    "count_stage1_s": ("KMC internal stage 1 time (splitting)", "seconds"),
    "count_stage2_s": ("KMC internal stage 2 time (counting)", "seconds"),
    "count_stage3_s": ("KMC internal stage 3 time (only reported by some KMC versions/modes)", "seconds"),
    "count_max_rss_gb": ("peak RSS of the timed kmc call (GNU time -v)", "decimal GB (KiB * 1024 / 10^9; GNU time reports KiB)"),
}
_ORDER_SHEET_MEANING: dict[str, tuple[str, str]] = {
    "order": (
        "signature/minimizer ordering scheme (canonical names, orders/README.md)",
        "default (kmc's shipped lex, is_allowed-protected) | lex (pure canonical lex, no denylist) | xor (random "
        "xor seed) | xor-pcg-hash (xor seed then PCG hash) | greedymini (GreedyMini UHS table) -- see the orders "
        "sheet for what each one actually does"),
    "gm_k": ("greedymini only: table k (= minimizer length k)", ""),
    "projection": ("greedymini only: DNA->binary projection (at0-cg1 default | at1-cg0 = kmc -oat1)", ""),
    "gm_table": ("greedymini only: table file stem (w{w}_k{k}_canonical)", ""),
    "gm_lut": (
        "greedymini only: table LUT backend (dna = 4^k DNA table | binary = 2^k rank LUT + runtime projection)",
        "gm_k / projection / gm_table / gm_lut are separate columns, never packed into the order name; two rows "
        "that differ in any of them are different variants and are never pooled"),
}
_COLUMN_MEANING: dict[str, tuple[str, str]] = {
    "dataset": ("dataset id", ""),
    "k": ("minimizer (signature) length (kmc -p)", "unified parameter names shared with the other apps; the CLI flags and run-dir names (k28_p9) keep kmc's own"),
    "w": ("window in k-mers, w = L - k + 1", ""),
    "L": ("k-mer length L = k + w - 1 (kmc -k)", ""),
    "input_format": ("gzip or plain FASTQ input", ""),
    "label": ("dataset label as given in the run's meta (e.g. read-length variants)", ""),
    "n_total_kmers": ("total number of k-mers (length L) in the dataset", ""),
    "n_reads": ("number of reads", ""),
    "avg_read_len": ("average read length in bases", "Excel formula in summary/gm_lut/runs where derivable"),
    "gc_fraction": ("GC fraction of the sampled bases", ""),
    "mem_gb": ("kmc -m memory limit used for the dataset's canonical run", "decimal GB"),
    "n_bins": ("kmc -n number of bins", "config_tag in CSVs encodes the config as {input}_n{bins}_t{threads}_m{mem}[_sm]"),
    "threads": ("kmc -t number of threads", ""),
    "strict_memory": ("kmc -sm strict-memory mode on/off", ""),
    "mem_limit_gb": ("kmc -m memory limit", "decimal GB"),
    "n_repeats": (
        "number of runs (repeats) aggregated into the row's mean/std; for the default baseline the repeats are "
        "runs of the same configuration that only measure timing noise", ""),
    "run_id": ("run directory relative to runs/ (dataset/k_p/config/scheme/repeat)", ""),
    "seed_idx": ("0..n-1 repeat index", ""),
    "seed": ("actual hex seed used for that repeat (custom schemes only)", ""),
    "signature_order_label": ("the kmc -ol label of the signature order", ""),
    "signature_order": ("kmc -o signature order argument", ""),
    "signature_order_mode": ("how the signature order is supplied to kmc", ""),
    "custom_order_file": ("custom order (GreedyMini table) file used, if any", ""),
    "binary_density_factor": ("GreedyMini rows only: density factor of the binary-alphabet (2-symbol) projection of the table", ""),
    "n_unique_counted": ("number of distinct k-mers counted by kmc", ""),
    "n_unique": ("number of unique k-mers (kmc stats)", ""),
    "n_below_ci": ("k-mers below the counter lower bound -ci (excluded)", ""),
    "n_above_cx": ("k-mers above the counter upper bound -cx (excluded)", ""),
    "tmp_disk_strict_mem_gb": ("KMC temp directory size in strict-memory mode", "decimal GB; empty when strict-memory was not used"),
    "user": ("user CPU time of the timed call (GNU time -v)", "seconds"),
    "sys": ("system CPU time of the timed call (GNU time -v)", "seconds"),
    "cpu_pct": ("CPU utilisation of the timed call (GNU time -v)", "percent"),
    "major_pagefaults": ("major page faults (GNU time -v)", "count; needed disk I/O"),
    "minor_pagefaults": ("minor page faults (GNU time -v)", "count"),
    "fs_inputs": ("file system inputs (GNU time -v)", "512-byte blocks; ~0 = input served from the page cache"),
    "read_gb": ("file system reads of the timed call = fs_inputs * 512 bytes", "decimal GB"),
    "fs_outputs": ("file system outputs (GNU time -v)", "512-byte blocks"),
    "ctx_voluntary": ("voluntary context switches (GNU time -v)", "count"),
    "ctx_involuntary": ("involuntary context switches (GNU time -v)", "count"),
    "harness_wall_s": ("wall time of the whole harness step around kmc (includes startup and bookkeeping)", "seconds"),
    "kmc_version": ("KMC version string", ""),
    "host": ("host the run was executed on", ""),
    "exit_code": ("exit code of the timed kmc call", "0 = success"),
    "count_pinned_cpus": ("CPU list the timed step was pinned to (from the run's cpus.txt)", PINNING_NOTE + "; empty for runs without it"),
    "paper_platform": ("platform of the published KMC 3 numbers (xeon_hdd | xeon_ssd)", ""),
    "paper_count_wall_s": ("published KMC 3 paper wall time", "seconds"),
    "paper_count_max_rss_gb": ("published KMC 3 paper memory", "decimal GB"),
    "paper_disk_gb": ("published KMC 3 paper temp disk usage", "decimal GB"),
    "description": ("what the order variant does", ""),
    "cli": ("kmc command-line flags that select this variant", ""),
    "is_allowed_denylist": ("whether kmc's is_allowed denylist (forbidden signatures) applies to this order", ""),
    "expected_density_factor_mean": (
        "theoretical density factor on i.i.d. random DNA, mean over the benchmark's seeds for "
        "xor/xor-pcg-hash/greedymini (0 std for the deterministic default/lex)", ""),
    "expected_density_factor_std": ("sample std of expected_density_factor over the seeds", ""),
}
_RATIO_NOTE = "Excel formula, recalculated on open; color scale: green < 1 (better) < red"


def _legend_entry(sheet: str, col: str) -> tuple[str, str]:
    if sheet == "vs_paper":
        if col == "count_wall_s":
            return ("mean wall time of the default (lex) runs of this dataset/L/config across the n_repeats repeats", "seconds")
        if col == "count_wall_s_std":
            return ("sample std of count_wall_s across the n_repeats repeats", "seconds")
        if col == "count_wall_s_ratio":
            return (
                "this run's count_wall_s / paper_count_wall_s",
                "< 1 = faster than the 2017 paper's hardware; do not expect 1.0. NOT the summary sheet's "
                "count_wall_s_ratio, which compares to this run's own lex baseline")
        if col in _METRIC_MEANING:
            return (_METRIC_MEANING[col][0] + " (default/lex runs, mean across repeats)", _METRIC_MEANING[col][1])
    if col in _ORDER_SHEET_MEANING:
        return _ORDER_SHEET_MEANING[col]
    if col in _COLUMN_MEANING:
        return _COLUMN_MEANING[col]
    for suf in ("_mean", "_std", "_ratio"):
        if col.endswith(suf) and col[: -len(suf)] in _METRIC_MEANING:
            m = col[: -len(suf)]
            what, units = _METRIC_MEANING[m]
            if suf == "_mean":
                return (f"mean of {m} across the row's n_repeats repeats: {what}", units)
            if suf == "_std":
                return (f"sample std of {m} across the row's n_repeats repeats", units)
            extra = ""
            if m == "density_factor":
                extra = ("; the summary sheet shows, under the table, Excel =CORREL of density_factor_ratio with each "
                         "other {metric}_ratio column (ratios, not raw values, because raw values span very different "
                         "absolute scales per dataset)")
            base = "this run's" if sheet == "runs" else "this row's mean"
            return (f"{base} {m} / the matching default (lex) row's value, same dataset/L/input/config", _RATIO_NOTE + extra)
    if col in _METRIC_MEANING:
        return _METRIC_MEANING[col]
    if col.endswith("_particular_density_factor_mean"):
        ds = col[: -len("_particular_density_factor_mean")]
        return (f"density factor on real {ds} reads, mean over the seeds of the kmc runs (a density-only scan, std 0, "
                "where a variant has no run on the dataset)", "")
    if col.endswith("_particular_density_factor_std"):
        ds = col[: -len("_particular_density_factor_std")]
        return (f"sample std of the {ds} particular density factor over the seeds", "")
    raise KeyError(f"no legend entry for column {col!r} (sheet {sheet})")


def overview_rows(sheet_fields: dict[str, Sequence[str]]) -> list[dict]:
    """Legend: one row per (column, meaning) over every data sheet's actual (unified) headers."""
    order: list[tuple] = []
    sheets_of: dict[tuple, list[str]] = {}
    for sheet, fields in sheet_fields.items():
        for col in fields:
            meaning, notes = _legend_entry(sheet, col)
            key = (col, meaning, notes)
            if key not in sheets_of:
                order.append(key)
                sheets_of[key] = []
            sheets_of[key].append(sheet)
    return [
        {"sheet": " / ".join(sheets_of[key]), "column": key[0], "meaning": key[1], "notes": key[2]}
        for key in order
    ]


ORDER_DESCRIPTIONS: dict[str, tuple[str, str]] = {
    # order -> (description, has is_allowed denylist: yes/no)
    "default": (
        "KMC's shipped lexicographic order: canonical (reverse-complement-min) integer value "
        "of the p-mer, EXCEPT is_allowed() (kmc_api/mmer.h) maps poly-A/poly-T-type "
        "low-complexity motifs to a sentinel 'infinity' rank so they can never be the "
        "minimizer -- a denylist from the original KMC2 paper (Sec. 2.2) specifically to "
        "bound the largest bin's size.",
        "yes",
    ),
    "lex": (
        "True canonical lexicographic order, MIN(kmer, revcomp(kmer)), no is_allowed "
        "denylist. Implemented as -or xor -os0 (verified order-equivalent offline; no KMC "
        "core change needed). Always seed 0 -- deterministic, not a real random scheme; "
        "repeated -n times purely for wall-time noise stats, like default.",
        "no",
    ),
    "xor": (
        "Random order: rank(kmer) = kmer XOR seed, canonicalized. No denylist. Plain XOR "
        "has no mixing guarantee, so seed-to-seed density varies substantially -- a bad "
        "seed can accidentally rank a common repeat motif (e.g. poly-A) very low.",
        "no",
    ),
    "xor-pcg-hash": (
        "Random order: rank(kmer) = RXS-M-XS64(kmer XOR seed), canonicalized. No denylist. "
        "PCG's output-mixing avalanche makes it far more seed-stable than plain xor.",
        "no",
    ),
    "greedymini": (
        "GreedyMini universal-hitting-set table (gm_k = p), optimized for low expected density "
        "on a binary alphabet, extended to DNA via the AT/GC projection (`projection`) + "
        "canonical/seed tie-break. No denylist. The LUT backend (`gm_lut`) is a storage "
        "choice, not a different order: dna = precomputed DNA-sized (4^p) LUT; binary = "
        "runtime-projected binary-alphabet (2^p) rank LUT (kmc's default backend for a binary "
        ".gm source) with identical measured density. The flipped at1-cg0 projection "
        "(kmc -oat1) IS a different DNA order from at0-cg1.",
        "no",
    ),
}

def _as_float(v) -> float | None:
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


_DATASET_SORT_PRIORITY = ("f_vesca", "g_gallus", "m_balbisiana", "h_sapiens_2", "h_sapiens_3")


def _dataset_sort_key(dataset: str) -> tuple:
    try:
        return (0, _DATASET_SORT_PRIORITY.index(dataset))
    except ValueError:
        return (1, dataset)


ORDER_CLI: dict[str, str] = {
    "default": "kmc -k<L> -p<k> (shipped order, no -or/-co)",
    "lex": "kmc -k<L> -p<k> -orxor -os0",
    "xor": "kmc -k<L> -p<k> -orxor -os<seed>",
    "xor-pcg-hash": "kmc -k<L> -p<k> -orxor-pcg-hash -os<seed>",
    "greedymini": "kmc -k<L> -p<k> -co<gm_table>[_dna].gm -ol<label> -os<seed>",
}
ORDER_CLI_AT1 = " -oat1"  # appended for projection at1-cg0 (binary .gm only)

# Deterministic orders: no seed, so their expected-density std is 0 (the
# cache leaves it blank because there is nothing to average over).
DETERMINISTIC_ORDERS = frozenset({"default", "lex"})


def orders_field_order(datasets: Sequence[str]) -> list[str]:
    """Shared `orders` layout (common/orders_sheet.py) minus `mod_t` (a
    mod-minimizer parameter kmc has no use for), plus kmc's `is_allowed_denylist`."""
    fields = [f for f in orders_sheet.orders_fields(datasets, with_gm_lut=True) if f != "mod_t"]
    fields.insert(fields.index("cli") + 1, "is_allowed_denylist")
    return fields


def orders_density_fields(datasets: Sequence[str]) -> list[str]:
    """The density columns that must never be empty."""
    return [*orders_sheet.EXPECTED_FIELDS, *orders_sheet.particular_fields(datasets)]


def check_orders_complete(rows: Sequence[dict], datasets: Sequence[str]) -> None:
    """Refuse to write an orders sheet with an empty density cell."""
    missing = orders_sheet.missing_cells(rows, orders_density_fields(datasets))
    if not missing:
        return
    lines = []
    for i, field in missing:
        r = rows[i]
        ident = "/".join(str(r.get(f)) for f in ("order", "gm_lut", "projection") if r.get(f))
        lines.append(f"  {ident} (k={r.get('k')} w={r.get('w')} L={r.get('L')}): {field}")
    die(
        f"orders sheet has {len(missing)} empty density cell(s):\n" + "\n".join(lines) + "\n"
        "fill expected cells with `python -m scripts expected-density --force`, and particular "
        "cells of (order, dataset) pairs without a kmc run with `python -m scripts density "
        "particular-orders --order O [--projection P] --datasets D`"
    )


def _mean_std_or_zero(values: Sequence[float]) -> tuple[float, float]:
    """(mean, sample std); std is 0 for a single (deterministic) value."""
    mean = statistics.mean(values)
    return round(mean, 4), (round(statistics.stdev(values), 4) if len(values) > 1 else 0.0)


def orders_rows(main_rows: Sequence[dict], out_root: Path) -> tuple[list[dict], list[str]]:
    """One row per (order variant, k, w, L), k = p (minimizer length), w = window,
    L = k + w - 1 (kmc's k-mer length): description + cli, theoretical expected
    density_factor on i.i.d. random DNA (mean/std over the benchmark's predefined
    seeds for xor/xor-pcg-hash/greedymini; std 0 for the deterministic default/lex),
    and a particular (real-data) density_factor mean/std per dataset.

    Particular = measured density_factor of the kmc runs (mean/std over seeds;
    bins/threads/input format do not change it, so every config counts, with
    repeats of one seed averaged first). A (variant, dataset) with no kmc run falls
    back to the density-only scan rows of kmc3_density_particular.csv
    (`scripts density particular-orders`), std 0 (one scan)."""
    scans = particular_variant_index(load_particular_rows(particular_csv_path(out_root)))

    usable = [r for r in main_rows if _as_float(r.get("density_factor")) is not None]
    datasets = sorted({r["dataset"] for r in main_rows}, key=_dataset_sort_key)

    def klw(r: dict) -> tuple[int, int, int]:
        p = int(r.get("p") or r.get("signature_len") or 0)
        k = int(r.get("k") or 0)
        return p, int(r.get("w") or (k - p + 1)), k

    # (variant, p, w, L) -> dataset -> seed_idx -> [density_factor over configs]
    by_id: dict[tuple, dict[str, dict]] = {}
    id_rows: dict[tuple, dict] = {}
    for r in main_rows:
        ident = (variant_key(r), *klw(r))
        id_rows.setdefault(ident, r)
    for r in usable:
        ident = (variant_key(r), *klw(r))
        by_id.setdefault(ident, {}).setdefault(r["dataset"], {}).setdefault(str(r.get("seed_idx")), []).append(
            _as_float(r["density_factor"])
        )

    expected = load_expected_orders(out_root)
    if not expected:
        print(
            f"warning: no {expected_orders_csv_path(out_root).name}; "
            "run `python -m scripts expected-density` for per-order expected density",
            file=sys.stderr,
        )

    out: list[dict] = []
    for ident in sorted(id_rows, key=lambda i: (variant_sort_key(id_rows[i]), i[1:])):
        scheme, p, w, L = ident
        order, gm_k, projection, gm_table, gm_lut = scheme
        desc, denylist = ORDER_DESCRIPTIONS.get(order, ("", ""))
        cli = ORDER_CLI.get(order, "")
        if projection == "at1-cg0":
            cli += ORDER_CLI_AT1
        row: dict = {
            "order": order,
            "gm_k": gm_k,
            "projection": projection,
            "gm_table": gm_table,
            "gm_lut": gm_lut,
            "description": desc,
            "cli": cli,
            "is_allowed_denylist": denylist,
        }
        orders_sheet.with_klw(row, p, w)
        if row["L"] != L:
            die(f"orders: k={p} w={w} does not give L={L} for {order}")
        exp = expected.get((order, projection or ""), {})
        if int(exp.get("k") or 0) == p and int(exp.get("w") or 0) == w:
            row["expected_density_factor_mean"] = _as_float(exp.get("density_factor_mean"))
            std = _as_float(exp.get("density_factor_std"))
            row["expected_density_factor_std"] = 0.0 if std is None and order in DETERMINISTIC_ORDERS else std
        for ds in datasets:
            seeds = by_id.get(ident, {}).get(ds)
            if seeds:
                vals = [statistics.mean(v) for v in seeds.values()]
            else:
                scan = scans.get((order, projection or "", ds, L, p))
                v = _as_float(scan.get("density_factor_dna")) if scan else None
                vals = [v] if v is not None else []
            if vals:
                mean, std = _mean_std_or_zero(vals)
                row[f"{ds}_particular_density_factor_mean"] = mean
                row[f"{ds}_particular_density_factor_std"] = std
        out.append(row)

    return out, datasets


def datasets_field_order() -> list[str]:
    return [
        "dataset",
        "label",
        "p",
        "w",
        "k",
        "input_format",
        "n_total_kmers",
        "n_reads",
        "avg_read_len",
        "gc_fraction",
        "mem_gb",
        "paper_count_max_rss_gb",
        "paper_disk_gb",
        "paper_count_wall_s",
    ]


def datasets_rows(main_rows: Sequence[dict], out_root: Path, paper_rows: Sequence[dict]) -> list[dict]:
    """One row per dataset: size/read-length stats (from its own default run) + paper reference."""
    particular = particular_index(load_particular_rows(particular_csv_path(out_root)))
    paper_index: dict[tuple, dict] = {}
    for r in paper_rows:
        try:
            key = (r["dataset"], int(r["k"]), r.get("input_format"))
        except (KeyError, TypeError, ValueError):
            continue
        paper_index[key] = r

    canon = [
        r
        for r in main_rows
        if order_of(r) == BASELINE_ORDER
        and int(r.get("n_bins") or 0) == 512
        and int(r.get("threads") or 0) == 12
        and not r.get("strict_memory")
    ]
    by_dataset: dict[str, dict] = {}
    for r in canon:
        by_dataset.setdefault(r["dataset"], r)  # first (any repeat has the same stats)

    out: list[dict] = []
    for ds in sorted(by_dataset, key=_dataset_sort_key):
        r = by_dataset[ds]
        k = int(r["k"])
        part = particular.get((ds, k, int(r.get("p") or 0)))
        paper = paper_index.get((ds, k, r.get("input_format")))
        out.append(
            {
                "dataset": ds,
                "label": r.get("label"),
                "k": k,
                "p": r.get("p"),
                "w": r.get("w"),
                "input_format": r.get("input_format"),
                "n_total_kmers": r.get("n_total_kmers"),
                "n_reads": r.get("n_reads"),
                "avg_read_len": r.get("avg_read_len"),
                "gc_fraction": _as_float(part.get("gc_fraction")) if part else None,
                "mem_gb": r.get("mem_limit_gb"),
                "paper_count_max_rss_gb": _as_float(paper.get("paper_ram_gb")) if paper else None,
                "paper_disk_gb": _as_float(paper.get("paper_disk_gb")) if paper else None,
                "paper_count_wall_s": _as_float(paper.get("paper_time_s")) if paper else None,
            }
        )
    return out


def cmd_evaluate(args: argparse.Namespace) -> None:
    results_root = results_dir(args.results_dir)
    runs_root = results_root / "runs"
    out_root = results_root / "csvs"  # exported csv tables + density caches
    xlsx_root = results_root  # workbook sits directly in results/
    for d in (runs_root, out_root, xlsx_root):
        d.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for run_dir in iter_run_dirs(runs_root):
        row = load_run_row(run_dir, runs_root=runs_root)
        if not row:
            continue
        rows.append(row)
    if not rows:
        die(f"no completed runs in {runs_root}")

    particular = load_density_context(out_root)
    attach_density_refs(rows, particular=particular)

    rows = sort_export_rows(rows)
    main_rows = rows_without_gm_lut(rows)
    attach_lex_ratios(main_rows)
    summary = summary_rows(main_rows)
    # The read-length-confound-check derived datasets (h_sapiens_2_shortsub /
    # _longsub / _shortsub_mf11, ...) get their own xlsx sheet instead of
    # mixing into "summary" -- CSVs (kmc3_summary.csv, kmc3_runs.csv) still
    # carry every row, this only changes xlsx sheet assignment.
    paper_dataset_rows = [r for r in main_rows if not is_readlen_derived_dataset(r.get("dataset") or "")]
    readlen_rows = [r for r in main_rows if is_readlen_derived_dataset(r.get("dataset") or "")]
    summary_xlsx = summary_xlsx_rows(paper_dataset_rows)
    readlen_xlsx = summary_xlsx_rows(readlen_rows) if readlen_rows else []
    runs_xlsx = [dict(r) for r in main_rows]
    gm_lut_xlsx = gm_lut_summary_rows(rows, main_rows)
    orders_xlsx, orders_datasets = orders_rows(main_rows, out_root)
    check_orders_complete(orders_xlsx, orders_datasets)
    n_lex_cmp, count_errors = validate_custom_vs_lex(main_rows)

    paper_path = Path(args.paper_csv)
    vs_paper: list[dict] = []
    paper_kmc3_rows: list[dict] = []
    if not paper_path.is_file():
        print(f"warning: no paper CSV at {paper_path}; skip vs-paper export", file=sys.stderr)
    else:
        paper_kmc3_rows = list(csv.DictReader(paper_path.open()))
        vs_paper = vs_paper_rows(rows, paper_kmc3_rows, args.paper_platform)
    datasets_xlsx = datasets_rows(main_rows, out_root, paper_kmc3_rows)

    if count_errors:
        for msg in count_errors:
            print(f"error: {msg}", file=sys.stderr)
        die(
            f"{len(count_errors)} default vs custom count mismatch(es); "
            "unique/total k-mer counts must be identical (super-k-mers may differ)"
        )

    by: dict[tuple, dict[str, dict]] = {}
    for r in rows:
        by.setdefault((r["dataset"], r["k"], r.get("signature_order_label") or ""), {})[r["input_format"]] = r
    for key, variants in sorted(by.items()):
        if "gzip" in variants and "plain" in variants:
            a = variants["gzip"]["n_unique_counted"]
            b = variants["plain"]["n_unique_counted"]
            if a is not None and b is not None and a != b:
                print(
                    f"warning: unique counted k-mers differ gzip vs plain for {key}: {a} vs {b}",
                    file=sys.stderr,
                )

    write_csv(out_root / "kmc3_runs.csv", main_rows, prune_empty_fields(runs_xlsx, runs_field_order()))
    if gm_lut_xlsx:
        attach_summary_lex_ratios(gm_lut_xlsx)
        gm_lut_csv_fields = prune_empty_fields(gm_lut_xlsx, summary_xlsx_field_order())
        write_csv(out_root / "kmc3_gm_lut.csv", gm_lut_xlsx, gm_lut_csv_fields)
    attach_summary_lex_ratios(summary_xlsx)
    write_csv(out_root / "kmc3_summary.csv", summary_xlsx, prune_empty_fields(summary_xlsx, summary_xlsx_field_order()))
    write_csv(out_root / "kmc3_orders.csv", orders_xlsx, orders_field_order(orders_datasets), unified=True)
    if vs_paper:
        write_csv(out_root / "kmc3_vs_paper.csv", vs_paper, vs_paper_field_order())
    for stale in ("kmc3_stages.csv", "kmc3_order_compare.csv", "kmc3_detail.csv", "kmc3_vs_lex.csv"):
        stale_path = out_root / stale
        if stale_path.is_file():
            stale_path.unlink()
            print(f"removed deprecated {stale_path}")

    summary_fields = prepare_agg_xlsx_sheet(summary_xlsx)
    summary_footer = summary_correlation_rows(summary_xlsx, summary_fields)
    readlen_fields = prepare_agg_xlsx_sheet(readlen_xlsx) if readlen_xlsx else []
    runs_fields = prepare_runs_xlsx_sheet(runs_xlsx)
    gm_lut_fields = prepare_agg_xlsx_sheet(gm_lut_xlsx) if gm_lut_xlsx else []
    orders_fields = orders_field_order(orders_datasets)  # complete by check_orders_complete
    datasets_fields = prune_empty_fields(datasets_xlsx, datasets_field_order())
    sheets = [
        make_sheet("orders", orders_xlsx, orders_fields),
        make_sheet("datasets", datasets_xlsx, datasets_fields),
        make_sheet("summary", summary_xlsx, summary_fields, footer_rows=summary_footer),
    ]
    if readlen_xlsx:
        sheets.append(make_sheet("readlen_check", readlen_xlsx, readlen_fields))
    if gm_lut_xlsx:
        sheets.append(make_sheet("gm_lut", gm_lut_xlsx, gm_lut_fields))
    sheets.append(make_sheet("runs", runs_xlsx, runs_fields))
    sheets.append(make_sheet("vs_paper", vs_paper, vs_paper_field_order()))
    xlsx_path = xlsx_root / "kmc3_results.xlsx"
    legend = overview_rows({sh["name"]: sh["fields"] for sh in sheets})
    sheets.insert(0, make_sheet("overview", legend, overview_field_order()))
    names = write_book(xlsx_path, sheets)
    print(f"wrote {xlsx_path} ({len(names)} sheets)")

    if n_lex_cmp:
        print(f"default vs custom counts match for {n_lex_cmp} comparison(s)")


def main(argv: list[str] | None = None) -> int:
    load_dotenv()
    default_results = bench_paths("kmc")["results"]
    p = argparse.ArgumentParser(prog="scripts export", description="Export kmc3 benchmark results to CSV and XLSX")
    p.add_argument("--results-dir", default=str(default_results))
    p.add_argument("--paper-csv", default=str(data_dir() / "paper_kmc3.csv"))
    p.add_argument("--paper-platform", default="xeon_hdd", choices=["xeon_hdd", "xeon_ssd"])
    args = p.parse_args(argv)
    cmd_evaluate(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

