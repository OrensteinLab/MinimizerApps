#!/usr/bin/env python3
"""Export KMC3 benchmark runs to CSV and XLSX (stdlib only).

Writes under ``--results-dir`` (default ``$KMC3_RESULTS_DIR`` or ``$KMC3_BENCH_ROOT/results``):

- ``kmc3_runs.csv`` — one row per run (metrics + telemetry)
- ``kmc3_summary.csv`` — one row per run (per seed; headline metrics)
- ``kmc3_vs_paper.csv`` — lex runs vs ``paper_kmc3.csv`` (when paper CSV exists)
- ``kmc3_results.xlsx`` — ``index`` (column reference), ``orders``, ``datasets``,
  ``summary`` (mean/std), ``readlen_check``, ``runs``, ``gm_lut``, ``vs_paper``
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
import zipfile
from collections.abc import Sequence
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape

from scripts.common import die, load_dotenv
from scripts.config import data_dir
from xlsx_util import col_letter, write_csv as _shared_write_csv
from scripts.density import (
    attach_density_refs,
    load_particular_rows,
    particular_csv_path,
    particular_index,
)
from scripts.expected_density import expected_orders_csv_path, load_expected_orders
from scripts.metrics import enrich_super_kmer_metrics
from scripts.schemes import (
    PROJECTION_AT0,
    PROJECTION_AT1,
    SCHEME_ORDER,
    is_at1_projection,
    resolve_scheme_label,
)


def bench_paths() -> tuple[Path, Path]:
    root = os.environ.get("KMC3_BENCH_ROOT", "").strip()
    if not root:
        die("KMC3_BENCH_ROOT is not set (see benchmark/.env)")
    root_path = Path(root)
    data = os.environ.get("KMC3_DATA_DIR", "").strip() or str(root_path)
    results = os.environ.get("KMC3_RESULTS_DIR", "").strip() or str(root_path / "results")
    return Path(data), Path(results)


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
    ("disk", ("disk",)),
    ("runtime", ("wall_time", "stage1", "stage2", "stage3")),
    ("memory", ("ram",)),
)
DENSITY_REF_FIELDS = ("binary_density_factor",)
# Summary XLSX: derived from counts; Excel formulas overwrite on export for live recalc.
XLSX_FORMULA_MEAN_METRICS = frozenset({"avg_super_kmer_len", "density_factor"})
XLSX_FORMULA_FIELDS = frozenset(
    {"w", "avg_read_len", "avg_super_kmer_len_mean", "density_factor_mean"}
)
RATIO_COLOR_SCALE_METRICS = frozenset(
    {"n_super_kmers", "density_factor", "disk", "wall_time", "stage1", "stage2", "stage3"}
)
XLSX_RATIO_COLOR_SHEETS = frozenset({"summary", "readlen_check", "gm_lut", "runs"})
SUMMARY_ALL_AGG_METRICS = tuple(m for _c, ms in SUMMARY_XLSX_METRIC_GROUPS for m in ms)
RATIO_METRICS = SUMMARY_ALL_AGG_METRICS
# Runtime metrics get a {metric}_z column (Welch's t-stat vs lex, pooled over both
# sides' seed/repeat std) so a ratio far from 1 can be told apart from seed noise.
# |z| >= 2 is a reasonable "probably real" cutoff. See README: "Run order and the
# lex baseline".
Z_SCORE_METRICS = ("wall_time", "stage1", "stage2")


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
    "k",
    "p",
    "w",
    "n_bins",
    "threads",
    "strict_memory",
    "mem_limit",
    "order_label",
    "n_seeds",
)
RUNS_FIELD_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (
        # Matches SUMMARY_PARAM_FIELDS' order for every shared field (dataset ..
        # order_label), so "runs" and "summary" line up column-for-column as far as
        # their differing row granularity allows; runs-only extras (run_id, label,
        # seed_idx, signature_order_seed, gm_lut_mode, gm_projection) are prepended/appended around it.
        "id / params",
        (
            "run_id",
            "dataset",
            "n_total_kmers",
            "avg_read_len",
            "n_reads",
            "input_format",
            "k",
            "p",
            "w",
            "n_bins",
            "threads",
            "strict_memory",
            "mem_limit",
            "order_label",
            "label",
            "seed_idx",
            "signature_order_seed",
            "gm_lut_mode",
            "gm_projection",
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
    ("disk", _with_ratio(("disk", "disk_strict_mem"))),
    ("runtime", _with_ratio(("wall_time", "stage1", "stage2", "stage3"))),
    ("memory", _with_ratio(("ram", "max_rss_kb"))),
    (
        "order",
        (
            "signature_order",
            "signature_order_label",
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
            "fs_outputs",
            "ctx_voluntary",
            "ctx_involuntary",
            "python_wall",
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
            "mem_limit",
            "paper_platform",
        ),
    ),
    (
        "runtime",
        ("n_seeds", "wall_time", "wall_time_std", "paper_time", "wall_time_ratio", "stage1", "stage2"),
    ),
    ("memory", ("ram", "paper_ram")),
    ("disk", ("disk", "paper_disk")),
    ("super-k-mer", ("n_unique_counted",)),
)
_SUMMARY_CSV_METRICS = (
    "n_super_kmers",
    "avg_super_kmer_len",
    "density_factor",
    "disk",
    "wall_time",
    "stage1",
    "stage2",
    "stage3",
    "ram",
)
CATEGORY_LABELS: dict[str, str] = {
    "id / params": "params",
    "super-k-mer": "super-k-mer",
    "disk": "disk (GB)",
    "runtime": "runtime (s)",
    "memory": "memory (GB)",
    "GM-DNA": "GM-DNA",
    "GM-bin": "GM-bin",
    "binary / DNA": "bin/DNA",
    "layout": "layout",
    "runtime (s)": "runtime (s)",
    "k-mer counts": "k-mers",
    "disk (GB)": "disk (GB)",
    "order": "order",
    "telemetry": "telemetry",
}
XLSX_HEADER_ROW = 2
XLSX_DATA_START_ROW = 3


def order_label(row: dict) -> str:
    """Human-readable scheme label for export; see scripts.schemes.SCHEME_ORDER
    ("default" = KMC's shipped is_allowed-protected order, "lexicographic" =
    true canonical lex with no denylist -- see the orders sheet)."""
    return resolve_scheme_label(row)


def scheme_sort_key(row: dict) -> tuple:
    """Same primary grouping as summary_sort_key (dataset -> bins -> threads ->
    mem -> strict_mem -> scheme), so "runs" rows appear in the same per-config
    blocks, in the same order, as their aggregated "summary" row -- k/input_format
    disambiguate within a config (summary_sort_key doesn't need them: one dataset
    is normally run at a single k/format per sweep), seed_idx orders repeats
    within a scheme (summary has no seed dimension, it's already aggregated)."""
    label = order_label(row)
    try:
        rank = SCHEME_ORDER.index(label)
    except ValueError:
        rank = len(SCHEME_ORDER)
    seed = row.get("signature_order_seed")
    return (
        row.get("dataset", ""),
        int(row.get("n_bins") or 0),
        int(row.get("threads") or 0),
        int(row.get("mem_limit") or 0),
        int(bool(row.get("strict_memory"))),
        int(row.get("k") or 0),
        row.get("input_format", ""),
        rank,
        int(row.get("seed_idx") or 0),
        seed or 0,
    )


def decorate_export_row(row: dict) -> dict:
    out = dict(row)
    out["order_label"] = order_label(row)
    out["p"] = row.get("signature_len")
    out["w"] = row.get("minimizer_window_w")
    return out


def sort_export_rows(rows: list[dict]) -> list[dict]:
    return sorted((decorate_export_row(r) for r in rows), key=scheme_sort_key)


def summary_sort_key(row: dict) -> tuple:
    label = row.get("order_label", "")
    try:
        rank = SCHEME_ORDER.index(label)
    except ValueError:
        rank = len(SCHEME_ORDER)
    return (
        row.get("dataset", ""),
        int(row.get("n_bins") or 0),
        int(row.get("threads") or 0),
        int(row.get("mem_limit") or 0),
        int(row.get("strict_memory") or 0),
        rank,
    )


def config_columns(r: dict) -> dict:
    return {
        "n_bins": r.get("n_bins"),
        "threads": r.get("threads"),
        "strict_memory": r.get("strict_memory"),
        "mem_limit": r.get("mem_limit"),
    }


def config_group_key(r: dict) -> tuple:
    return (r["dataset"], int(r["k"]), r["input_format"], r.get("config_tag") or "")


def agg_group_key(r: dict) -> tuple:
    return (*config_group_key(r), r.get("order_label") or "")


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
            "order_label": head.get("order_label"),
            "n_seeds": len({p.get("seed_idx") for p in parts}),
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
        int(row.get("mem_limit") or 0),
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
        if (row.get("order_label") or "") == "default":
            lex_at[_xlsx_match_key(row, with_seed=with_seed)] = i
    for i, row in enumerate(rows, start=XLSX_DATA_START_ROW):
        lex_i = lex_at.get(_xlsx_match_key(row, with_seed=with_seed))
        if lex_i is None:
            continue
        idx = i - XLSX_DATA_START_ROW
        is_lex = (row.get("order_label") or "") == "default"
        for metric in RATIO_METRICS:
            rname = f"{metric}_ratio"
            mcol = col.get(metric_field(metric))
            rcol = col.get(rname)
            if rname not in col or not mcol or not rcol or mcol == rcol:
                continue
            if is_lex:
                rows[idx][rname] = {"f": "=1"}
                rows[idx][f"__{rname}_value"] = 1.0
            else:
                rows[idx][rname] = {"f": f"={mcol}{i}/{mcol}{lex_i}"}
                # Shadow plain-float value (cells hold formulas, not numbers) so the
                # color scale's min/max can be computed in Python -- see
                # _xlsx_ratio_color_scales_xml.
                this_val = _as_float(row.get(metric_field(metric)))
                lex_val = _as_float(rows[lex_i - XLSX_DATA_START_ROW].get(metric_field(metric)))
                if this_val is not None and lex_val:
                    rows[idx][f"__{rname}_value"] = this_val / lex_val


def _apply_lex_zscore_formulas(
    rows: list[dict],
    fields: Sequence[str],
    *,
    metrics: Sequence[str],
    merge_fields: frozenset[str],
) -> None:
    """Excel `{metric}_z` = Welch's t-stat of (this mean - lex mean) vs pooled seed std.

    Both sides need >= 2 seeds/repeats for a non-blank value; otherwise there is no
    seed-noise estimate to compare against and the cell is left blank via IFERROR.

    `n_seeds` is vertically merged (see XLSX_MERGE_FIELDS): only the first of a run of
    equal values actually has a cell written, so a non-anchor row's own `n_seeds` cell
    is blank in the XML. Formulas must resolve through `_xlsx_formula_ref_row` (like
    `_apply_derived_metric_formulas` already does for k/p/w/n_total_kmers) instead of
    referencing `n_col{i}` directly, or they silently divide by a blank cell and
    IFERROR blanks the result -- which is why z came out empty for merged-away rows.
    """
    col = {name: _xlsx_col(i) for i, name in enumerate(fields)}
    n_col = col.get("n_seeds")
    if not n_col:
        return
    lex_at: dict[tuple, int] = {}
    for i, row in enumerate(rows, start=XLSX_DATA_START_ROW):
        if (row.get("order_label") or "") == "default":
            lex_at[_xlsx_match_key(row, with_seed=False)] = i
    for i, row in enumerate(rows, start=XLSX_DATA_START_ROW):
        lex_i = lex_at.get(_xlsx_match_key(row, with_seed=False))
        if lex_i is None:
            continue
        idx = i - XLSX_DATA_START_ROW
        is_lex = (row.get("order_label") or "") == "default"
        for metric in metrics:
            zname = f"{metric}_z"
            mean_col, std_col = col.get(f"{metric}_mean"), col.get(f"{metric}_std")
            zcol = col.get(zname)
            if not (mean_col and std_col and zcol):
                continue
            if is_lex:
                rows[idx][zname] = {"f": "=0"}
                continue
            n_row = _xlsx_formula_ref_row(rows, merge_fields, "n_seeds", idx)
            n_lex_row = _xlsx_formula_ref_row(rows, merge_fields, "n_seeds", lex_i - XLSX_DATA_START_ROW)
            rows[idx][zname] = {
                "f": (
                    "=IFERROR(("
                    f"{mean_col}{i}-{mean_col}{lex_i}"
                    f")/SQRT(({std_col}{i}^2/{n_col}{n_row})+({std_col}{lex_i}^2/{n_col}{n_lex_row})),\"\")"
                )
            }


def _xlsx_formula_ref_row(
    rows: Sequence[dict],
    merge_fields: frozenset[str],
    field: str,
    row_idx: int,
) -> int:
    """Excel row to reference when a vertically merged column has no cell on this row."""
    if field not in merge_fields:
        return XLSX_DATA_START_ROW + row_idx
    anchor = row_idx
    while anchor > 0:
        if rows[anchor].get("dataset") != rows[anchor - 1].get("dataset"):
            break
        if _merge_field_value(rows[anchor], field) != _merge_field_value(rows[anchor - 1], field):
            break
        anchor -= 1
    return XLSX_DATA_START_ROW + anchor


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
        rows[idx]["w"] = {"f": f"={k_col}{kr}-{p_col}{pr}+1"}
        rows[idx][avg_field] = {"f": f"={nt_col}{ntr}/{ns_col}{i}"}
        rows[idx][df_field] = {"f": f"={ns_col}{i}/{nt_col}{ntr}*({w_col}{wr}+1)"}
        if nr_col and ar_col:
            rows[idx]["avg_read_len"] = {"f": f"={nt_col}{ntr}/{nr_col}{nrr}+({k_col}{kr}-1)"}


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
    _apply_lex_zscore_formulas(
        rows, fields, metrics=Z_SCORE_METRICS, merge_fields=XLSX_MERGE_FIELDS["summary"]
    )


SUMMARY_CORRELATION_BANNER_LABEL = "Correlation (density_factor_ratio vs disk_ratio):"


def summary_correlation_banner(rows: Sequence[dict], fields: Sequence[str]) -> dict | None:
    """One live `=CORREL(...)` Excel formula (absolute references) over this
    sheet's own density_factor_ratio and disk_ratio columns -- a single
    sheet-wide statistic, not a per-row metric, so it's rendered once as a
    banner row above the table (see write_xlsx's `banners` param) instead of
    a repeated column. Those two metrics (not raw density_factor, which spans
    wildly different absolute scales per dataset) are the scale-invariant pair
    found to track almost perfectly (r~1) across every scheme/dataset in this
    benchmark. Row refs here assume the unshifted (no-banner) row layout --
    write_xlsx shifts them to match the banner's actual placement."""
    if not rows or "density_factor_ratio" not in fields or "disk_ratio" not in fields:
        return None
    col = {name: _xlsx_col(i) for i, name in enumerate(fields)}
    d_col, k_col = col["density_factor_ratio"], col["disk_ratio"]
    last_row = XLSX_DATA_START_ROW + len(rows) - 1
    formula = f"=CORREL(${d_col}${XLSX_DATA_START_ROW}:${d_col}${last_row},${k_col}${XLSX_DATA_START_ROW}:${k_col}${last_row})"
    return {"f": formula}


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
            if metric in Z_SCORE_METRICS:
                fields.append(f"{metric}_z")
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
        "k",
        "p",
        "w",
        "n_bins",
        "threads",
        "strict_memory",
        "mem_limit",
        "order_label",
        "signature_order_label",
        *_with_ratio(("n_super_kmers", "avg_super_kmer_len", "density_factor")),
        "binary_density_factor",
        *_with_ratio(("disk", "wall_time", "stage1", "stage2", "stage3", "ram")),
    ]


def vs_paper_field_order() -> list[str]:
    return _fields_from_groups(VS_PAPER_FIELD_GROUPS)


def summary_rows(rows: list[dict]) -> list[dict]:
    """One row per dataset / config / scheme with headline metrics."""
    out: list[dict] = []
    for r in rows:
        label = r.get("order_label") or resolve_scheme_label(r)
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
            "order_label": label,
            "signature_order_label": r.get("signature_order_label", ""),
            "binary_density_factor": r.get("binary_density_factor"),
            "n_unique_counted": r.get("n_unique_counted"),
        }
        for metric in _SUMMARY_CSV_METRICS:
            row[metric] = r.get(metric)
            row[f"{metric}_ratio"] = r.get(f"{metric}_ratio")
        out.append(row)
    return sorted(out, key=summary_sort_key)


def _gb_from_kb(rss_kb: float | None) -> float | None:
    return None if rss_kb is None else round(rss_kb / (1024 * 1024), 3)


def _gb_from_bytes(n: float | None) -> float | None:
    return None if n is None else round(n / (1000 ** 3), 3)


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
    order = meta.get("signature_order") or stats.get("json_signature_order") or "default"
    order_label = meta.get("signature_order_label") or stats.get("json_signature_order_label") or ""
    order_mode = meta.get("signature_order_mode") or stats.get("json_signature_order_mode") or ""
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
        "mem_limit": meta.get("mem_limit_gb"),
        "n_bins": meta.get("n_bins"),
        "strict_memory": meta.get("strict_memory"),
        "config_tag": meta.get("config_tag", ""),
        "signature_len": meta.get("signature_len"),
        "signature_order": order,
        "signature_order_mode": order_mode,
        "signature_order_label": order_label,
        "signature_order_seed": order_seed,
        "seed_idx": meta.get("seed_idx"),
        "custom_order_file": meta.get("custom_order_file", ""),
        "kmc_version": meta.get("kmc_version"),
        "host": meta.get("host"),
        "exit_code": meta.get("exit_code"),
        "wall_time": wall,
        "kmc_total_s": stats["kmc_total_s"],
        "time_elapsed_s": tinfo["time_elapsed_s"],
        "python_wall": meta.get("python_wall_s"),
        "user": tinfo["user_s"],
        "sys": tinfo["sys_s"],
        "cpu_pct": tinfo["cpu_pct"],
        "stage1": stats["stage1_s"],
        "stage2": stats["stage2_s"],
        "stage3": stats["stage3_s"],
        "ram": _gb_from_kb(rss_kb),
        "max_rss_kb": rss_kb,
        "disk": _gb_from_bytes(stats["tmp_bytes"]),
        "disk_strict_mem": _gb_from_bytes(stats["tmp_strict_mem_bytes"]),
        "major_pagefaults": tinfo["major_pagefaults"],
        "minor_pagefaults": tinfo["minor_pagefaults"],
        "fs_inputs": tinfo["fs_inputs"],
        "fs_outputs": tinfo["fs_outputs"],
        "ctx_voluntary": tinfo["ctx_voluntary"],
        "ctx_involuntary": tinfo["ctx_involuntary"],
        "n_unique": stats["n_unique"],
        "n_unique_counted": stats["n_unique_counted"],
        "n_total_kmers": stats["n_total_kmers"],
        "n_reads": stats["n_reads"],
        "n_below_ci": stats["n_below_ci"],
        "n_above_cx": stats["n_above_cx"],
        "n_super_kmers": stats["n_super_kmers"],
        "gm_lut_mode": meta.get("gm_lut_mode", ""),
        "gm_projection": "",
    }
    row["gm_projection"] = gm_projection_from_meta(meta, row)
    enrich_super_kmer_metrics(row)
    return row


def gm_projection_from_meta(meta: dict, row: dict) -> str:
    """DNA->binary projection of a GreedyMini run ("" for non-GM runs).

    lib/run.sh passes ``-oat1`` to kmc for GM_PROJECTION=at1-cg0 and meta.json
    records the full command line; older runs (and at0 runs) never had it, so
    they stay at the default at0-cg1. The scheme dir / label suffix
    (``greedymini_binary_at1`` / ``..._at1_s0x..``) is the fallback for a meta
    without a command."""
    probe = dict(row)
    probe["gm_projection"] = ""
    if not resolve_scheme_label(probe).startswith("greedymini"):
        return ""
    command = meta.get("command")
    if isinstance(command, list) and any(str(a) == "-oat1" for a in command):
        return PROJECTION_AT1
    return PROJECTION_AT1 if is_at1_projection(probe) else PROJECTION_AT0


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
        if (row.get("order_label") or "") == "default":
            lex_at[row_key(row)] = row
    for row in rows:
        is_lex = (row.get("order_label") or "") == "default"
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


def is_gm_binary_run(row: dict) -> bool:
    if (row.get("gm_lut_mode") or "").strip().lower() == "binary":
        return True
    return order_label(row).startswith("greedymini_binary")


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
    """Summary-format sheet: lex + GM-DNA + GM-bin (+ GM-bin-at1) per GM LUT layout (ratio vs lex)."""
    lut_runs = [r for r in all_rows if is_gm_lut_run(r)]
    if not lut_runs:
        return []

    lex_pool: dict[tuple, list[dict]] = {}
    for row in main_rows:
        if (row.get("order_label") or "") == "default":
            lex_pool.setdefault(_xlsx_match_key(row, with_seed=False), []).append(row)

    by_layout: dict[tuple, dict[str, list[dict]]] = {}
    for row in lut_runs:
        mode = (row.get("gm_lut_mode") or "").strip().lower()
        if not mode:
            mode = "binary" if is_gm_binary_run(row) else "dna"
        if is_at1_projection(row):
            mode += "_at1"
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
        for mode in ("dna", "binary", *sorted(m for m in modes if m not in ("dna", "binary"))):
            if mode in modes:
                bundle.extend(modes[mode])

    if not bundle:
        return []
    return summary_xlsx_rows(bundle)


def write_csv(path: Path, rows: list[dict], fields: Sequence[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = list(fields) if fields else list(rows[0].keys())
    _shared_write_csv(path, fieldnames, rows)
    print(f"wrote {path} ({len(rows)} rows)")


def _xlsx_col(n: int) -> str:
    return col_letter(n + 1)


XLSX_STYLE_DEFAULT = 0
XLSX_STYLE_DATASET_BORDER = 1
# cellXfs pairs per fill: even=plain, odd=+dataset bottom border
XLSX_STYLE_BINS = {512: 2, 1024: 4, 1536: 6}
XLSX_STYLE_THREADS = {12: 8, 24: 10}
XLSX_STYLE_MERGE_OFFSET = 12
XLSX_STYLE_COUNT = 12
XLSX_COLOR_SHEETS = frozenset({"summary", "readlen_check", "gm_lut", "runs"})
XLSX_DATASET_BORDER_SHEETS = frozenset({"summary", "readlen_check", "gm_lut", "runs"})
XLSX_MERGE_DATASET = (
    "dataset",
    "n_total_kmers",
    "n_reads",
    "input_format",
    "label",
)
XLSX_MERGE_MINIMIZER = ("k", "p", "w")
XLSX_MERGE_CONFIG = ("n_bins", "threads", "strict_memory", "mem_limit")
XLSX_MERGE_FIELDS: dict[str, frozenset[str]] = {
    "summary": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, "n_seeds")
    ),
    "readlen_check": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, "n_seeds")
    ),
    "gm_lut": frozenset(
        (*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG, "n_seeds")
    ),
    "runs": frozenset((*XLSX_MERGE_DATASET, "avg_read_len", *XLSX_MERGE_MINIMIZER, *XLSX_MERGE_CONFIG)),
}


def _dataset_last_row_indices(rows: Sequence[dict]) -> set[int]:
    """Zero-based indices of the last row per dataset block."""
    if not rows:
        return set()
    out: set[int] = set()
    for i, row in enumerate(rows):
        if i == len(rows) - 1 or row.get("dataset") != rows[i + 1].get("dataset"):
            out.add(i)
    return out


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
        return (row.get("n_total_kmers"), row.get("n_reads"), row.get("k"))
    if field == "w":
        k = row.get("k")
        p = row.get("p") if row.get("p") is not None else row.get("signature_len")
        return (k, p)
    return _merge_value_key(row.get(field))


def _merge_rows_equal(rows: Sequence[dict], i: int, j: int, field: str) -> bool:
    if rows[i].get("dataset") != rows[j].get("dataset"):
        return False
    return _merge_field_value(rows[i], field) == _merge_field_value(rows[j], field)


def _vertical_merge_plan(
    rows: Sequence[dict],
    fields: Sequence[str],
    merge_fields: frozenset[str],
    row_offset: int = 0,
) -> tuple[list[str], set[tuple[int, str]], set[tuple[int, str]]]:
    """Return (mergeCell refs, skipped cells, merge anchor cells to center)."""
    refs: list[str] = []
    skip: set[tuple[int, str]] = set()
    anchors: set[tuple[int, str]] = set()
    if len(rows) < 2:
        return refs, skip, anchors
    col = {name: _xlsx_col(i) for i, name in enumerate(fields)}
    data_start = XLSX_DATA_START_ROW + row_offset

    def _close_run(field: str, letter: str, run_start: int, run_end: int) -> None:
        if run_end <= run_start:
            return
        r1 = data_start + run_start
        r2 = data_start + run_end
        refs.append(f"{letter}{r1}:{letter}{r2}")
        anchors.add((run_start, field))
        for r in range(run_start + 1, run_end + 1):
            skip.add((r, field))

    for field in fields:
        if field not in merge_fields:
            continue
        letter = col.get(field)
        if not letter:
            continue
        run_start = 0
        for i in range(1, len(rows)):
            if _merge_rows_equal(rows, i - 1, i, field):
                continue
            _close_run(field, letter, run_start, i - 1)
            run_start = i
        _close_run(field, letter, run_start, len(rows) - 1)
    return refs, skip, anchors


def _ratio_color_scale_fields(fields: Sequence[str]) -> list[str]:
    """Ratio columns that get red (high) → green (low) color scales."""
    return [
        f
        for f in fields
        if f.endswith("_ratio") and f.removesuffix("_ratio") in RATIO_COLOR_SCALE_METRICS
    ]


def _ratio_color_scale_bounds(rows: Sequence[dict], field: str) -> tuple[float, float]:
    """Low/high stops for `field`'s color scale, clamped to always include 1 as a
    genuine interior point. A data-driven max/min (Excel `cfvo type="min"/"max"`)
    degenerates when every value sits on one side of 1 -- e.g. gm_lut's disk_ratio
    never exceeds 1, so a dynamic max would coincide exactly with the val=1 midpoint
    and the white/red half of the scale would never render. Computed in Python from
    shadow `__{field}_value` floats (see _apply_lex_ratio_formulas) since the cells
    themselves hold Excel formulas, not values."""
    vals = [v for r in rows if (v := r.get(f"__{field}_value")) is not None]
    lo = min([1.0, *vals])
    hi = max([1.0, *vals])
    if lo == 1.0:
        lo = 0.99
    if hi == 1.0:
        hi = 1.01
    return lo, hi


def _xlsx_ratio_color_scales_xml(fields: Sequence[str], rows: Sequence[dict], row_offset: int = 0) -> str:
    """3-color ratio scale: green below 1, white at 1, red above 1."""
    if not rows:
        return ""
    data_start = XLSX_DATA_START_ROW + row_offset
    last_row = data_start + len(rows) - 1
    col = {name: _xlsx_col(i) for i, name in enumerate(fields)}
    parts: list[str] = []
    for priority, field in enumerate(_ratio_color_scale_fields(fields), start=1):
        letter = col.get(field)
        if not letter:
            continue
        lo, hi = _ratio_color_scale_bounds(rows, field)
        sqref = f"{letter}{data_start}:{letter}{last_row}"
        parts.append(
            f'<conditionalFormatting sqref="{sqref}">'
            f'<cfRule type="colorScale" priority="{priority}">'
            "<colorScale>"
            f'<cfvo type="num" val="{lo}"/>'
            '<cfvo type="num" val="1"/>'
            f'<cfvo type="num" val="{hi}"/>'
            '<color rgb="FF63BE7B"/>'
            '<color rgb="FFFFFFFF"/>'
            '<color rgb="FFF8696B"/>'
            "</colorScale>"
            "</cfRule>"
            "</conditionalFormatting>"
        )
    n_ratio = len(parts)
    for offset, metric in enumerate(Z_SCORE_METRICS):
        field = f"{metric}_z"
        letter = col.get(field)
        if not letter:
            continue
        priority = n_ratio + offset + 1
        sqref = f"{letter}{data_start}:{letter}{last_row}"
        # |z| >= 2: bold + tint, flagging "probably a real difference from lex,
        # not seed noise" (Welch's t-stat computed in _apply_lex_zscore_formulas).
        parts.append(
            f'<conditionalFormatting sqref="{sqref}">'
            f'<cfRule type="cellIs" dxfId="0" priority="{priority}" operator="greaterThanOrEqual">'
            "<formula>2</formula></cfRule>"
            f'<cfRule type="cellIs" dxfId="1" priority="{priority + 1}" operator="lessThanOrEqual">'
            "<formula>-2</formula></cfRule>"
            "</conditionalFormatting>"
        )
    return "".join(parts)


def _xlsx_merges_xml(refs: Sequence[str]) -> str:
    if not refs:
        return ""
    cells = "".join(f'<mergeCell ref="{r}"/>' for r in refs)
    return f'<mergeCells count="{len(refs)}">{cells}</mergeCells>'


def _xlsx_data_cell_style(
    field: str,
    row: dict,
    *,
    dataset_border: bool,
    colorize: bool,
    merged: bool = False,
) -> int:
    """Pick cellXf for n_bins (cyan) / threads (magenta) fills, optional dataset bottom border."""
    if dataset_border and field not in ("n_bins", "threads"):
        style = XLSX_STYLE_DATASET_BORDER
    elif not colorize or field not in ("n_bins", "threads"):
        style = XLSX_STYLE_DATASET_BORDER if dataset_border else XLSX_STYLE_DEFAULT
    else:
        try:
            val = int(row.get(field) or 0)
        except (TypeError, ValueError):
            val = 0
        if field == "n_bins":
            base = XLSX_STYLE_BINS.get(val)
            if base is None:
                style = XLSX_STYLE_DATASET_BORDER if dataset_border else XLSX_STYLE_DEFAULT
            else:
                style = base + (1 if dataset_border else 0)
        else:
            base = XLSX_STYLE_THREADS.get(val)
            if base is None:
                style = XLSX_STYLE_DATASET_BORDER if dataset_border else XLSX_STYLE_DEFAULT
            else:
                style = base + (1 if dataset_border else 0)
    if merged:
        style += XLSX_STYLE_MERGE_OFFSET
    return style


def _xlsx_cell(value, col: int, row: int, *, style: int = XLSX_STYLE_DEFAULT) -> str:
    ref = f"{_xlsx_col(col)}{row}"
    s_attr = f' s="{style}"' if style else ""
    if isinstance(value, dict) and "f" in value:
        return f'<c r="{ref}"{s_attr}><f>{xml_escape(str(value["f"]))}</f></c>'
    if value is None or value == "":
        return f'<c r="{ref}"{s_attr}/>'
    if isinstance(value, bool):
        return f'<c r="{ref}"{s_attr} t="b"><v>{1 if value else 0}</v></c>'
    if isinstance(value, int) and not isinstance(value, bool):
        return f'<c r="{ref}"{s_attr}><v>{value}</v></c>'
    if isinstance(value, float):
        return f'<c r="{ref}"{s_attr}><v>{value}</v></c>'
    text = xml_escape(str(value), {"\"": "&quot;"})
    return f'<c r="{ref}"{s_attr} t="inlineStr"><is><t xml:space="preserve">{text}</t></is></c>'


def _xlsx_xf_xml(
    fill_id: int,
    border_id: int,
    *,
    apply_fill: bool = False,
    apply_border: bool = False,
    center: bool = False,
) -> str:
    attrs = [
        'numFmtId="0"',
        'fontId="0"',
        f'fillId="{fill_id}"',
        f'borderId="{border_id}"',
        'xfId="0"',
    ]
    if apply_fill:
        attrs.append('applyFill="1"')
    if apply_border:
        attrs.append('applyBorder="1"')
    if center:
        attrs.append('applyAlignment="1"')
        return f'<xf {" ".join(attrs)}><alignment horizontal="center" vertical="center"/></xf>'
    return f'<xf {" ".join(attrs)}/>'


def _xlsx_styles_xml() -> str:
    """Styles: default, dataset border, fills, and centered variants for merged cells."""
    fill_rgb = [
        None,
        None,
        "FFCCFFFF",
        "FF99FFFF",
        "FF66FFFF",
        "FFFFCCFF",
        "FFFF99FF",
    ]
    fill_parts = [
        '<fill><patternFill patternType="none"/></fill>',
        '<fill><patternFill patternType="gray125"/></fill>',
    ]
    for rgb in fill_rgb[2:]:
        fill_parts.append(
            f'<fill><patternFill patternType="solid">'
            f'<fgColor rgb="{rgb}"/><bgColor indexed="64"/></patternFill></fill>'
        )
    xfs = [
        _xlsx_xf_xml(0, 0),
        _xlsx_xf_xml(0, 1, apply_border=True),
    ]
    for fill_id in (2, 3, 4):
        xfs.append(_xlsx_xf_xml(fill_id, 0, apply_fill=True))
        xfs.append(_xlsx_xf_xml(fill_id, 1, apply_fill=True, apply_border=True))
    for fill_id in (5, 6):
        xfs.append(_xlsx_xf_xml(fill_id, 0, apply_fill=True))
        xfs.append(_xlsx_xf_xml(fill_id, 1, apply_fill=True, apply_border=True))
    xfs.extend(
        [
            _xlsx_xf_xml(0, 0, center=True),
            _xlsx_xf_xml(0, 1, apply_border=True, center=True),
        ]
    )
    for fill_id in (2, 3, 4):
        xfs.append(_xlsx_xf_xml(fill_id, 0, apply_fill=True, center=True))
        xfs.append(_xlsx_xf_xml(fill_id, 1, apply_fill=True, apply_border=True, center=True))
    for fill_id in (5, 6):
        xfs.append(_xlsx_xf_xml(fill_id, 0, apply_fill=True, center=True))
        xfs.append(_xlsx_xf_xml(fill_id, 1, apply_fill=True, apply_border=True, center=True))
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><color theme="1"/><name val="Calibri"/></font></fonts>'
        f'<fills count="{len(fill_parts)}">{"".join(fill_parts)}</fills>'
        '<borders count="2">'
        "<border><left/><right/><top/><bottom/><diagonal/></border>"
        '<border><left/><right/><top/>'
        '<bottom style="medium"><color auto="1"/></bottom>'
        "<diagonal/></border>"
        "</borders>"
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        f'<cellXfs count="{len(xfs)}">{"".join(xfs)}</cellXfs>'
        "<dxfs count=\"2\">"
        # dxf 0: |z| >= 2 "worse than lex" — bold + light red (matches the ratio color scale's red = worse)
        '<dxf><font><b/><color rgb="FF9C0006"/></font>'
        '<fill><patternFill><bgColor rgb="FFFFC7CE"/></patternFill></fill></dxf>'
        # dxf 1: |z| >= 2 "better than lex" — bold + light green
        '<dxf><font><b/><color rgb="FF006100"/></font>'
        '<fill><patternFill><bgColor rgb="FFC6EFCE"/></patternFill></fill></dxf>'
        "</dxfs>"
        "</styleSheet>"
    )


def _xlsx_freeze_header_rows(row_offset: int = 0) -> str:
    header_row = XLSX_HEADER_ROW + row_offset
    data_row = XLSX_DATA_START_ROW + row_offset
    return (
        "<sheetViews>"
        '<sheetView workbookViewId="0">'
        f'<pane ySplit="{header_row}" topLeftCell="A{data_row}" '
        'activePane="bottomLeft" state="frozen"/>'
        f'<selection pane="bottomLeft" activeCell="A{data_row}" '
        f'sqref="A{data_row}"/>'
        "</sheetView>"
        "</sheetViews>"
    )


def _agg_metric_base(field: str) -> str | None:
    for suffix in ("_mean", "_std", "_ratio", "_z"):
        if field.endswith(suffix):
            return field[: -len(suffix)]
    return None


def prune_empty_fields(rows: Sequence[dict], fields: Sequence[str]) -> list[str]:
    """Drop columns whose values are all empty; drop ratio/std when mean/base is empty."""
    if not rows:
        return list(fields)

    def col_empty(field: str) -> bool:
        return all(_is_empty_cell(r.get(field)) for r in rows)

    drop: set[str] = set()
    for field in fields:
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
            if field in XLSX_FORMULA_FIELDS:
                continue
            base = _agg_metric_base(field)
            if base in XLSX_FORMULA_MEAN_METRICS:
                continue
            if base and field.endswith(("_ratio", "_std", "_z")):
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


def _agg_header_label(metric: str, suffix: str, _group_metrics: Sequence[str]) -> str:
    tag = "avg" if suffix == "_mean" else "std"
    return f"{metric}_{tag}"


def _category_display(category: str) -> str:
    return CATEGORY_LABELS.get(category, category)


def _build_xlsx_header_labels() -> dict[str, str]:
    labels: dict[str, str] = {
        "run_id": "id",
        "dataset": "dataset",
        "label": "label",
        "input_format": "fmt",
        "k": "k",
        "p": "p",
        "w": "w",
        "order_label": "order",
        "n_bins": "bins",
        "strict_memory": "strict_memory",
        "threads": "threads",
        "mem_limit": "mem_limit",
        "seed_idx": "seed",
        "signature_order_seed": "seed_hex",
        "gm_lut_mode": "lut",
        "gm_projection": "projection",
        "n_super_kmers": "n_super_kmers",
        "n_total_kmers": "n_total_kmers",
        "avg_read_len": "avg_read_len",
        "avg_super_kmer_len": "avg_super_kmer_len",
        "density_factor": "density_factor",
        "binary_density_factor": "binary_density_factor",
        "n_unique_counted": "n_unique_counted",
        "paper_time": "paper_time",
        "paper_ram": "paper_ram",
        "paper_disk": "paper_disk",
        "wall_time_ratio": "wall_time_ratio",
        "n_unique": "unique",
        "n_reads": "reads",
        "n_below_ci": "below",
        "n_above_cx": "above",
        "disk": "disk",
        "disk_strict_mem": "disk_strict_mem",
        "wall_time": "wall_time",
        "stage1": "stage1",
        "stage2": "stage2",
        "stage3": "stage3",
        "ram": "ram",
        "max_rss_kb": "rss_kb",
        "signature_order": "sig",
        "signature_order_label": "sig_lbl",
        "signature_order_mode": "sig_mode",
        "custom_order_file": "gm_file",
        "user": "user",
        "sys": "sys",
        "cpu_pct": "cpu",
        "major_pagefaults": "maj_pf",
        "minor_pagefaults": "min_pf",
        "fs_inputs": "fs_in",
        "fs_outputs": "fs_out",
        "ctx_voluntary": "ctx_vol",
        "ctx_involuntary": "ctx_inv",
        "python_wall": "python_wall",
        "kmc_version": "kmc_ver",
        "host": "host",
        "exit_code": "exit",
        "n_seeds": "seeds",
        "stage1_ratio": "stage1_ratio",
        "stage2_ratio": "stage2_ratio",
        "disk_ratio": "disk_ratio",
    }
    for _cat, metrics in SUMMARY_XLSX_METRIC_GROUPS:
        for metric in metrics:
            for suffix in ("_mean", "_std"):
                labels[f"{metric}{suffix}"] = _agg_header_label(metric, suffix, metrics)
    return labels


XLSX_HEADER_LABELS = _build_xlsx_header_labels()


def _xlsx_header_label(field: str, sheet: str = "") -> str:
    if field.endswith("_z") and sheet in ("summary", "gm_lut"):
        base = field[: -len("_z")]
        return f"{base}_z(lex)"
    if field in XLSX_HEADER_LABELS:
        return XLSX_HEADER_LABELS[field]
    for prefix in ("dna_", "bin_"):
        if field.startswith(prefix):
            rest = field[len(prefix) :]
            return XLSX_HEADER_LABELS.get(rest, rest)
    return field


def _categories_from_groups(fields: Sequence[str], groups: Sequence[tuple[str, Sequence[str]]]) -> list[str]:
    lookup = {name: cat for cat, names in groups for name in names}
    return [lookup.get(field, "") for field in fields]


def _xlsx_category_band_row(
    fields: Sequence[str], categories: Sequence[str], row_offset: int = 0
) -> tuple[str, list[str]]:
    """Return (category header band row XML, mergeCell refs). Normally row 1;
    shifted down by `row_offset` on a sheet with a banner row above it."""
    if not fields:
        return "", []
    band_row = 1 + row_offset
    cells: list[str] = []
    merges: list[str] = []
    i = 0
    while i < len(fields):
        cat = categories[i]
        start = i
        while i < len(fields) and categories[i] == cat:
            i += 1
        end = i - 1
        cells.append(_xlsx_cell(_category_display(cat), start, band_row))
        if end > start:
            merges.append(f"{_xlsx_col(start)}{band_row}:{_xlsx_col(end)}{band_row}")
    return f'<row r="{band_row}">{"".join(cells)}</row>', merges


def _xlsx_field_categories(raw_name: str, fields: Sequence[str]) -> list[str] | None:
    if raw_name in ("summary", "readlen_check"):
        return [_summary_field_category(f) for f in fields]
    if raw_name == "runs":
        return _categories_from_groups(fields, RUNS_FIELD_GROUPS)
    if raw_name == "gm_lut":
        return [_summary_field_category(f) for f in fields]
    if raw_name == "vs_paper":
        return _categories_from_groups(fields, VS_PAPER_FIELD_GROUPS)
    if raw_name == "index":
        return ["metrics reference"] * len(fields)
    if raw_name == "orders":
        return ["order reference"] * len(fields)
    if raw_name == "datasets":
        return ["dataset reference"] * len(fields)
    return None


_XLSX_CELL_REF_RE = re.compile(r"(\$?)([A-Z]{1,3})(\$?)(\d+)")


def _xlsx_shift_formula_row_refs(value, delta: int):
    """Shift every cell-reference row number in a `{"f": "=..."}` formula cell
    by `delta`; leave plain values (and the ratio-color shadow `__x_value`
    floats) untouched. Used so a sheet's pre-baked formulas (built assuming no
    banner row above the table) still point at the right cells once a banner
    pushes the category band/header/data rows down."""
    if not (isinstance(value, dict) and "f" in value):
        return value
    shifted = _XLSX_CELL_REF_RE.sub(
        lambda m: f"{m.group(1)}{m.group(2)}{m.group(3)}{int(m.group(4)) + delta}", value["f"]
    )
    return {**value, "f": shifted}


def _xlsx_banner_row_xml(label: str, formula_cell: dict) -> str:
    """Row 1: a single sheet-wide statistic shown once, above the category
    band/header/table, instead of as a repeated per-row column."""
    cells = _xlsx_cell(label, 0, 1) + _xlsx_cell(formula_cell, 1, 1)
    return f'<row r="1">{cells}</row>'


def write_xlsx(
    path: Path,
    sheets: Sequence[tuple[str, list[dict]]],
    field_orders: dict[str, list[str]] | None = None,
    banners: dict[str, tuple[str, dict]] | None = None,
) -> None:
    """Write a multi-sheet .xlsx without extra Python packages."""
    field_orders = field_orders or {}
    banners = banners or {}
    used: dict[str, int] = {}
    names: list[str] = []
    payloads: list[str] = []

    for raw_name, rows in sheets:
        name = (raw_name or "Sheet")[:31]
        if name in used:
            used[name] += 1
            name = f"{name[:28]}_{used[name]}"
        else:
            used[name] = 1
        names.append(name)
        fields: list[str] = list(field_orders.get(raw_name, []))
        if rows and not fields:
            for row in rows:
                for k in row.keys():
                    if k not in fields:
                        fields.append(k)
        if rows and fields and raw_name not in field_orders:
            fields = prune_empty_fields(rows, fields)
        banner = banners.get(raw_name)
        row_offset = 1 if banner else 0
        if row_offset and rows:
            rows = [
                {k: _xlsx_shift_formula_row_refs(v, row_offset) for k, v in row.items()} for row in rows
            ]
        xml_rows = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>']
        xml_rows.append(
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        )
        if raw_name in ("index", "orders", "datasets", "summary", "readlen_check", "runs", "gm_lut", "vs_paper") and fields:
            xml_rows.append(_xlsx_freeze_header_rows(row_offset))
        xml_rows.append("<sheetData>")
        merge_refs: list[str] = []
        if banner:
            banner_label, banner_formula = banner
            banner_formula = _xlsx_shift_formula_row_refs(banner_formula, row_offset)
            xml_rows.append(_xlsx_banner_row_xml(banner_label, banner_formula))
        if fields:
            header_row = XLSX_HEADER_ROW + row_offset
            data_start = XLSX_DATA_START_ROW + row_offset
            categories = _xlsx_field_categories(raw_name, fields)
            if categories:
                cat_row, cat_merge = _xlsx_category_band_row(fields, categories, row_offset)
                xml_rows.append(cat_row)
                merge_refs.extend(cat_merge)
            merge_fields = XLSX_MERGE_FIELDS.get(raw_name, frozenset())
            vertical_refs, merge_skip, merge_anchors = _vertical_merge_plan(
                rows, fields, merge_fields, row_offset
            )
            merge_refs.extend(vertical_refs)
            cells = "".join(_xlsx_cell(_xlsx_header_label(h, raw_name), i, header_row) for i, h in enumerate(fields))
            xml_rows.append(f'<row r="{header_row}">{cells}</row>')
            dataset_last = _dataset_last_row_indices(rows) if raw_name in XLSX_DATASET_BORDER_SHEETS else set()
            colorize = raw_name in XLSX_COLOR_SHEETS
            for r_i, row in enumerate(rows, start=data_start):
                row_idx = r_i - data_start
                border = row_idx in dataset_last
                row_cells: list[str] = []
                for c, h in enumerate(fields):
                    if (row_idx, h) in merge_skip:
                        continue
                    row_cells.append(
                        _xlsx_cell(
                            row.get(h),
                            c,
                            r_i,
                            style=_xlsx_data_cell_style(
                                h,
                                row,
                                dataset_border=border,
                                colorize=colorize,
                                merged=(row_idx, h) in merge_anchors,
                            ),
                        )
                    )
                xml_rows.append(f'<row r="{r_i}">{"".join(row_cells)}</row>')
        xml_rows.append("</sheetData>")
        merge_xml = _xlsx_merges_xml(merge_refs)
        if merge_xml:
            xml_rows.append(merge_xml)
        if raw_name in XLSX_RATIO_COLOR_SHEETS and fields and rows:
            ratio_fmt = _xlsx_ratio_color_scales_xml(fields, rows, row_offset)
            if ratio_fmt:
                xml_rows.append(ratio_fmt)
        xml_rows.append("</worksheet>")
        payloads.append("\n".join(xml_rows))

    sheet_refs: list[str] = []
    overrides = [
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]
    rels = [
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
    rels.append("</Relationships>")

    styles = _xlsx_styles_xml()

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
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels", "\n".join(rels))
        z.writestr("xl/styles.xml", styles)
        for i, body in enumerate(payloads, start=1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", body)

    print(f"wrote {path} ({len(names)} sheets)")


def validate_custom_vs_lex(rows: list[dict]) -> tuple[int, list[str]]:
    """Check custom schemes match lex k-mer counts; return (n_comparisons, errors)."""
    by_key: dict[tuple, dict[str, dict]] = {}
    for r in rows:
        order = r.get("signature_order") or "default"
        label = r.get("signature_order_label") or ""
        slot = "lex" if order == "default" or not label else f"custom:{label}"
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
        if order_label(r) != "default":
            continue
        groups.setdefault(config_group_key(r), []).append(r)

    out: list[dict] = []
    for reps in groups.values():
        head = reps[0]
        key = (head["dataset"], int(head["k"]), platform, head["input_format"])
        p = paper_index.get(key)
        wt_vals = [float(r["wall_time"]) for r in reps if r.get("wall_time") is not None]
        wt_mean, wt_std = _mean_std(wt_vals)
        s1_vals = [float(r["stage1"]) for r in reps if r.get("stage1") is not None]
        s2_vals = [float(r["stage2"]) for r in reps if r.get("stage2") is not None]
        ram_vals = [float(r["ram"]) for r in reps if r.get("ram") is not None]
        disk_vals = [float(r["disk"]) for r in reps if r.get("disk") is not None]
        paper_t = float(p["paper_time_s"]) if p else None
        out.append(
            {
                "dataset": head["dataset"],
                "k": head["k"],
                "input_format": head["input_format"],
                **config_columns(head),
                "paper_platform": platform,
                "n_seeds": len({r.get("seed_idx") for r in reps}),
                "wall_time": round(wt_mean, 4) if wt_mean is not None else None,
                "wall_time_std": round(wt_std, 4) if wt_std is not None else None,
                "paper_time": paper_t,
                "wall_time_ratio": _ratio(wt_mean, paper_t),
                "stage1": round(_mean_std(s1_vals)[0], 4) if s1_vals else None,
                "stage2": round(_mean_std(s2_vals)[0], 4) if s2_vals else None,
                "ram": round(_mean_std(ram_vals)[0], 4) if ram_vals else None,
                "paper_ram": float(p["paper_ram_gb"]) if p else None,
                "disk": round(_mean_std(disk_vals)[0], 4) if disk_vals else None,
                "paper_disk": float(p["paper_disk_gb"]) if p else None,
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


def index_field_order() -> list[str]:
    return ["sheet", "column", "meaning", "notes"]


def index_rows() -> list[dict]:
    """Static reference for every non-obvious column in this workbook."""

    def row(sheet: str, column: str, meaning: str, notes: str = "") -> dict:
        return {"sheet": sheet, "column": column, "meaning": meaning, "notes": notes}

    return [
        row("all", "dataset / k / p / w", "dataset id, k-mer length, minimizer length, window (w = k - p + 1)"),
        row("all", "input_format", "gzip or plain FASTQ input"),
        row(
            "all", "n_bins / threads / strict_memory / mem_limit",
            "kmc -n / -t / -sm / -m for this run",
            "config_tag in CSVs encodes these as {input}_n{bins}_t{threads}_m{mem}[_sm]",
        ),
        row(
            "all", "order_label", "signature/minimizer ordering scheme",
            "default (kmc's shipped lex, is_allowed-protected) | lexicographic "
            "(pure canonical lex, no denylist) | hash (random xor) | pcg (xor+PCG) | "
            "greedymini_dna / greedymini_binary (GreedyMini UHS table) | "
            "greedymini_binary_at1 (same binary table, flipped at1-cg0 projection, kmc -oat1) "
            "-- see the orders sheet for what each one actually does",
        ),
        row(
            "summary / gm_lut", "n_seeds", "number of repeat runs aggregated into this row's mean/std",
            "lexicographic is repeated the same way custom schemes are (see README: "
            "'Run order and the lex baseline') so it also gets a real n_seeds/std, not always 1",
        ),
        row(
            "summary / gm_lut", "{metric}_mean / {metric}_std",
            "mean and sample std of the metric across this row's n_seeds repeats",
        ),
        row(
            "summary", "correlation banner (row 1)",
            "Excel =CORREL(density_factor_ratio range, disk_ratio range), shown once above "
            "the table (not a column)",
            "a sheet-wide statistic, not a per-row one, so it isn't repeated down every row; "
            "recalculates live on open. density_factor_ratio and disk_ratio (not raw "
            "density_factor, which spans very different absolute scales per dataset) are the "
            "scale-invariant pair found to track almost perfectly (r~1) across every "
            "scheme/dataset in this benchmark",
        ),
        row(
            "summary / gm_lut / runs", "{metric}_ratio",
            "this row's value / the matching default (lex) row's value, same dataset/k/input/config",
            "Excel formula, recalculated on open; color scale: green < 1 (better) < red",
        ),
        row(
            "summary / gm_lut", "wall_time_z / stage1_z / stage2_z",
            "Welch's t-statistic of (this mean - lex mean) against both sides' pooled seed/repeat std",
            "|z| >= 2 (highlighted red/green) is a reasonable 'probably a real difference, not "
            "seed noise' cutoff; a ratio far from 1 with small |z| is more likely noise, "
            "especially on short/fast configs where a handful of seeds gives a poor std estimate. "
            "Blank when either side has too few repeats (e.g. n_seeds=1) to estimate a std from.",
        ),
        row("summary / gm_lut / runs", "n_super_kmers", "number of super-k-mers KMC's splitter produced"),
        row(
            "summary / gm_lut / runs", "avg_super_kmer_len",
            "n_total_kmers / n_super_kmers -- longer is better (fewer, denser super-k-mers)",
        ),
        row(
            "summary / gm_lut / runs", "density_factor",
            "n_super_kmers / n_total_kmers * (w + 1) -- lower is better (closer-to-theoretical minimizer density)",
        ),
        row(
            "summary / gm_lut / runs", "binary_density_factor",
            "GreedyMini rows only: density factor of the binary-alphabet (2-symbol) projection of the table",
        ),
        row("summary / gm_lut / runs", "disk", "KMC temp directory size", "decimal GB (bytes / 1000^3)"),
        row(
            "summary / gm_lut / runs", "wall_time / stage1 / stage2 / stage3",
            "total wall time and KMC's internal stage timings (stage1=splitting, stage2=counting)", "seconds",
        ),
        row("summary / gm_lut / runs", "ram", "peak RSS (GNU time -v)", "GiB (KiB / 1024^2)"),
        row(
            "runs", "seed_idx / signature_order_seed",
            "0..n-1 repeat index, and the actual hex seed used for that repeat (custom schemes only)",
        ),
        row(
            "vs_paper", "wall_time / wall_time_std / n_seeds",
            "default (lex) runs only, averaged across repeats for this dataset/k/config",
            "NOT the same as the summary sheet's wall_time_ratio: vs_paper compares to the "
            "*published 2017 KMC3 paper* numbers, summary compares to *this run's own lex baseline*",
        ),
        row(
            "vs_paper", "wall_time_ratio", "this run's wall_time / the paper's published wall_time_s",
            "< 1 = faster than the 2017 paper's hardware; do not expect 1.0, hardware has moved on",
        ),
        row(
            "orders", "(sheet)", "one row per scheme: what it is, whether it has the is_allowed "
            "denylist, a small-sample offline particular-density diagnostic (GreedyMini rows only), "
            "theoretical expected density_factor on i.i.d. random DNA (mean/std across this "
            "benchmark's predefined seeds for hash/pcg/greedymini -- deterministic, std blank, for "
            "default/lexicographic), and real measured density_factor per dataset (mean across "
            "seeds, canonical bins=512/t12 config)",
        ),
        row(
            "datasets", "(sheet)", "one row per dataset: size, read-length, and config stats "
            "(from its own default/lex run), plus published KMC3-paper reference numbers where available",
        ),
    ]


ORDER_DESCRIPTIONS: dict[str, tuple[str, str]] = {
    # order_label -> (description, has is_allowed denylist: yes/no)
    "default": (
        "KMC's shipped lexicographic order: canonical (reverse-complement-min) integer value "
        "of the p-mer, EXCEPT is_allowed() (kmc_api/mmer.h) maps poly-A/poly-T-type "
        "low-complexity motifs to a sentinel 'infinity' rank so they can never be the "
        "minimizer -- a denylist from the original KMC2 paper (Sec. 2.2) specifically to "
        "bound the largest bin's size.",
        "yes",
    ),
    "lexicographic": (
        "True canonical lexicographic order, MIN(kmer, revcomp(kmer)), no is_allowed "
        "denylist. Implemented as -or hash -os0 (verified order-equivalent offline; no KMC "
        "core change needed). Always seed 0 -- deterministic, not a real random scheme; "
        "repeated -n times purely for wall-time noise stats, like default.",
        "no",
    ),
    "hash": (
        "Random order: rank(kmer) = kmer XOR seed, canonicalized. No denylist. Plain XOR "
        "has no mixing guarantee, so seed-to-seed density varies substantially -- a bad "
        "seed can accidentally rank a common repeat motif (e.g. poly-A) very low.",
        "no",
    ),
    "pcg": (
        "Random order: rank(kmer) = RXS-M-XS64(kmer XOR seed), canonicalized. No denylist. "
        "PCG's output-mixing avalanche makes it far more seed-stable than plain hash.",
        "no",
    ),
    "greedymini_dna": (
        "GreedyMini universal-hitting-set table, optimized for low expected density on a "
        "binary alphabet, extended to DNA via AT/GC projection + canonical/seed tie-break. "
        "No denylist. Stored as a precomputed DNA-sized (4^p) LUT.",
        "no",
    ),
    "greedymini_binary": (
        "Same GreedyMini order as greedymini_dna -- identical measured density -- stored "
        "instead as a runtime-projected binary-alphabet (2^p) LUT (kmc's default backend "
        "for a binary .gm source); a storage-format optimization, not a different order.",
        "no",
    ),
    "greedymini_binary_at1": (
        "The binary GreedyMini table under the flipped DNA->binary projection (A/T->1, "
        "C/G->0; kmc -oat1, GM_PROJECTION=at1-cg0) instead of the default at0-cg1 -- a "
        "different DNA order from greedymini_binary. Expected density left blank: the "
        "expected-density Monte Carlo reads the at0 _dna.gm sort key and has no at1 model.",
        "no",
    ),
}

# Schemes the expected-density cache (scripts/expected_density.py) has no model
# for: their expected_density_* cells are left blank on purpose, never filled
# from another scheme's entry.
NO_EXPECTED_DENSITY = frozenset({"greedymini_binary_at1"})


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


def orders_field_order(datasets: Sequence[str]) -> list[str]:
    return [
        "order_label",
        "description",
        "is_allowed_denylist",
        "particular_density_binary_diagnostic_avg",
        "expected_density_avg",
        "expected_density_std",
        *[f"density_{d}" for d in datasets],
    ]


def orders_rows(main_rows: Sequence[dict], out_root: Path) -> tuple[list[dict], list[str]]:
    """One row per scheme: description + a small-sample particular-density diagnostic
    (GreedyMini rows only) + theoretical expected density_factor on i.i.d. random DNA
    (mean/std over this benchmark's predefined seeds for hash/pcg/greedymini; a single
    deterministic value, std=None, for default/lexicographic) + real measured
    density_factor per dataset (mean across seeds, canonical bins=512/t12 config)."""
    particular = particular_index(load_particular_rows(particular_csv_path(out_root)))
    bin_diag_vals = [
        v for r in particular.values() if (v := _as_float(r.get("density_factor_binary"))) is not None
    ]
    bin_diag_avg = round(sum(bin_diag_vals) / len(bin_diag_vals), 4) if bin_diag_vals else None

    canon = [
        r
        for r in main_rows
        if int(r.get("n_bins") or 0) == 512
        and int(r.get("threads") or 0) == 12
        and not r.get("strict_memory")
    ]
    datasets = sorted({r["dataset"] for r in canon}, key=_dataset_sort_key)

    by_scheme_ds_density: dict[tuple, list[float]] = {}
    for r in canon:
        label = r.get("order_label") or order_label(r)
        d = _as_float(r.get("density_factor"))
        if d is not None:
            by_scheme_ds_density.setdefault((label, r["dataset"]), []).append(d)

    present = [s for s in SCHEME_ORDER if any(k[0] == s for k in by_scheme_ds_density)]

    # Read-only: like kmc3_density_particular.csv, this cache is built by a
    # separate explicit step (`python -m scripts expected-density`), not
    # auto-computed here -- it's a real (if one-time) Monte Carlo cost, and
    # this function needs to stay cheap/pure.
    expected = load_expected_orders(out_root)
    if not expected:
        print(
            f"warning: no {expected_orders_csv_path(out_root).name}; "
            "run `python -m scripts expected-density` for per-order expected density",
            file=sys.stderr,
        )

    out: list[dict] = []
    for scheme in present:
        desc, denylist = ORDER_DESCRIPTIONS.get(scheme, ("", ""))
        row: dict = {"order_label": scheme, "description": desc, "is_allowed_denylist": denylist}
        if scheme.startswith("greedymini"):
            row["particular_density_binary_diagnostic_avg"] = bin_diag_avg
        exp = {} if scheme in NO_EXPECTED_DENSITY else expected.get(scheme, {})
        row["expected_density_avg"] = _as_float(exp.get("density_factor_avg"))
        row["expected_density_std"] = _as_float(exp.get("density_factor_std"))
        for ds in datasets:
            dvals = by_scheme_ds_density.get((scheme, ds))
            if dvals:
                row[f"density_{ds}"] = round(sum(dvals) / len(dvals), 4)
        out.append(row)

    return out, datasets


def datasets_field_order() -> list[str]:
    return [
        "dataset",
        "label",
        "k",
        "p",
        "w",
        "input_format",
        "n_total_kmers",
        "n_reads",
        "avg_read_len",
        "gc_fraction",
        "mem_gb",
        "paper_ram_gb",
        "paper_disk_gb",
        "paper_time_s",
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
        if (r.get("order_label") or order_label(r)) == "default"
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
                "mem_gb": r.get("mem_limit"),
                "paper_ram_gb": _as_float(paper.get("paper_ram_gb")) if paper else None,
                "paper_disk_gb": _as_float(paper.get("paper_disk_gb")) if paper else None,
                "paper_time_s": _as_float(paper.get("paper_time_s")) if paper else None,
            }
        )
    return out


def cmd_evaluate(args: argparse.Namespace) -> None:
    out_root = results_dir(args.results_dir)
    runs_root = out_root / "runs"
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

    write_csv(out_root / "kmc3_runs.csv", main_rows, runs_field_order())
    if gm_lut_xlsx:
        attach_summary_lex_ratios(gm_lut_xlsx)
        gm_lut_csv_fields = prune_empty_fields(gm_lut_xlsx, summary_xlsx_field_order())
        write_csv(out_root / "kmc3_gm_lut.csv", gm_lut_xlsx, gm_lut_csv_fields)
    write_csv(out_root / "kmc3_summary.csv", summary, summary_field_order())
    if vs_paper:
        write_csv(out_root / "kmc3_vs_paper.csv", vs_paper, vs_paper_field_order())
    for stale in ("kmc3_stages.csv", "kmc3_order_compare.csv", "kmc3_detail.csv", "kmc3_vs_lex.csv"):
        stale_path = out_root / stale
        if stale_path.is_file():
            stale_path.unlink()
            print(f"removed deprecated {stale_path}")

    summary_fields = prepare_agg_xlsx_sheet(summary_xlsx)
    summary_banner = summary_correlation_banner(summary_xlsx, summary_fields)
    readlen_fields = prepare_agg_xlsx_sheet(readlen_xlsx) if readlen_xlsx else []
    runs_fields = prepare_runs_xlsx_sheet(runs_xlsx)
    gm_lut_fields = prepare_agg_xlsx_sheet(gm_lut_xlsx) if gm_lut_xlsx else []
    orders_fields = prune_empty_fields(orders_xlsx, orders_field_order(orders_datasets))
    datasets_fields = prune_empty_fields(datasets_xlsx, datasets_field_order())
    sheets: list[tuple[str, list[dict]]] = [("index", index_rows())]
    sheets.append(("orders", orders_xlsx))
    sheets.append(("datasets", datasets_xlsx))
    sheets.append(("summary", summary_xlsx))
    if readlen_xlsx:
        sheets.append(("readlen_check", readlen_xlsx))
    if gm_lut_xlsx:
        sheets.append(("gm_lut", gm_lut_xlsx))
    sheets.append(("runs", runs_xlsx))
    sheets.append(("vs_paper", vs_paper))
    field_orders: dict[str, list[str]] = {
        "index": index_field_order(),
        "orders": orders_fields,
        "datasets": datasets_fields,
        "summary": summary_fields,
        "runs": runs_fields,
        "vs_paper": vs_paper_field_order(),
    }
    if readlen_xlsx:
        field_orders["readlen_check"] = readlen_fields
    if gm_lut_xlsx:
        field_orders["gm_lut"] = gm_lut_fields

    banners: dict[str, tuple[str, dict]] = {}
    if summary_banner is not None:
        banners["summary"] = (SUMMARY_CORRELATION_BANNER_LABEL, summary_banner)

    write_xlsx(out_root / "kmc3_results.xlsx", sheets, field_orders=field_orders, banners=banners)

    if n_lex_cmp:
        print(f"default vs custom counts match for {n_lex_cmp} comparison(s)")


def main(argv: list[str] | None = None) -> int:
    load_dotenv()
    _, default_results = bench_paths()
    p = argparse.ArgumentParser(prog="scripts export", description="Export kmc3 benchmark results to CSV and XLSX")
    p.add_argument("--results-dir", default=str(default_results))
    p.add_argument("--paper-csv", default=str(data_dir() / "paper_kmc3.csv"))
    p.add_argument("--paper-platform", default="xeon_hdd", choices=["xeon_hdd", "xeon_ssd"])
    args = p.parse_args(argv)
    cmd_evaluate(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
