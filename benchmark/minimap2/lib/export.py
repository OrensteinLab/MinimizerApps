#!/usr/bin/env python3
"""Aggregate every runs/**/meta.json into per-sheet CSVs and one XLSX,
mirroring sshash's/kmc's own export.py layout (a lighter version: no
buckets/space-breakdown sheets, no equivalence hard-fail gate -- minimap2's
benchmark doesn't have those concepts). Writes, under --results-dir:

- minimap2_overview.csv, minimap2_datasets.csv, minimap2_summary.csv,
  minimap2_orders.csv, minimap2_runs.csv -- one CSV per sheet, same rows as
  the XLSX below
- minimap2_results.xlsx -- overview (per-column legend, like sshash's own
  'index' sheet), datasets, summary, orders, runs, with sshash's
  own layout conventions: a category band (row 1) grouping columns into
  sections, vertical merges on columns that repeat within a dataset's
  block (dataset/k/w), a bottom-border divider after each dataset's last
  row (summary/runs only -- datasets/orders are one row per dataset/order
  already), and the 'smoke' dataset excluded throughout (fast-iteration/CI
  dataset, not a source of real accuracy/density conclusions -- see
  README.md's "Datasets" section).

summary's {metric}_ratio columns (this order's mean / native's mean for the
same dataset -- native itself is always 1) isolate each order's own effect
from run-to-run noise on the same machine, and get a green(better than
native)->white(1)->red(worse) color scale. Polarity is metric-specific, not
uniformly "lower is better" the way sshash's own metrics are: density_factor/
build+map time+memory/frac_incorrectly_mapped are lower-is-better, but
frac_mapped/frac_correct*/mean_identity are higher-is-better -- see
LOWER_IS_BETTER/HIGHER_IS_BETTER below, both consulted by attach_native_ratios
and _color_scale_xml so "green" always actually means "better than native"
and never gets flipped for an accuracy metric.

datasets: one row per dataset, describing what it actually is (genome size,
read count/source, repeat region), read straight off disk (genome.fa.fai via
samtools faidx, truth.tsv, repeats.bed) rather than hand-maintained, so it
can't go stale.

summary: one row per (dataset, k, w, order) -- mean/std of every numeric
metric across seeds (n_seeds=1 for native/greedymini, which don't vary with
seed; see lib/orders.sh's order_is_seed_dependent).

orders: one row per order (not per dataset x order) -- w/t/description/cli
and expected (theoretical) density_factor don't vary by dataset, plus each
dataset's particular (measured) density_factor pivoted into its own column
pair. expected_density.py computes the theoretical value via Monte Carlo
simulation on i.i.d. random DNA (mirroring sshash's own scripts/
expected_density.py), reusing benchmark/common/density.py's shared
primitives -- the same ones order.cpp's own C++ is built on -- rather than
re-deriving minimap2's order formulas a third time; cached to
minimap2_density_expected_orders.csv since it's a pure function of k/w/t
(independent of any actual run).

density_factor is on minimap2's own sketching scale: lib/run.sh always passes
`-x map-pb`, i.e. HPC, so minimap2 picks minimizers over homopolymer-
compressed k-mers. density_factor = total occurrences * (w+1) / n_hpc_kmers,
where n_hpc_kmers counts the compressed k-mers minimap2 actually sketches
(lib/hpc_kmers.py) -- same scale as the Monte Carlo expected value, which
also compresses first. write_meta.py's meta.json density_factor instead
divides by the *raw* length (minimap2 -d's "average spacing"), understating
it by the genome's compression ratio (~0.70 on human) and any N gaps; that's
kept as density_factor_raw_bp in the runs sheet only.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
import subprocess
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from xlsx_util import col_letter, root_rels_xml, write_csv, xml_escape  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hpc_kmers import hpc_kmers  # noqa: E402

# Grouped so each row-1 category band is one contiguous span:
# id/params -> density -> total -> build -> map -> accuracy.
RUN_FIELDS = [
    "dataset", "k", "w", "mid_occ_cap", "order", "mod_t", "seed", "seed_idx", "n_reads", "n_reads_in_repeat",
    "density_factor", "density_factor_raw_bp", "distinct_minimizers", "avg_occurrences",
    "total_wall_time_sec",
    "build_wall_time_sec", "build_peak_rss_gb",
    "mid_occ", "map_wall_time_sec", "map_peak_rss_gb",
    "frac_mapped", "frac_correct", "frac_incorrectly_mapped", "frac_correct_in_repeat",
    "frac_correct_strict50", "frac_correct_strict50_in_repeat",
    "mean_identity",
]

# mid_occ_cap: the map-time `-f` condition (lib/run.sh's MID_OCC) -- this
# label for minimap2's own per-index default, else the fixed absolute cap.
# Every per-block computation (native ratio, divider, merge, correlation) is
# keyed on (dataset, mid_occ_cap), so a fixed-cap control never gets compared
# against the other condition's native.
DEFAULT_MID_OCC_CAP = "default"

KB_PER_GB = 1024 * 1024  # binary GB, matching how peak_rss_kb itself is binary (ru_maxrss-derived)

SUMMARY_METRICS = (
    "density_factor", "avg_occurrences",
    "total_wall_time_sec",
    "build_wall_time_sec", "build_peak_rss_gb",
    "map_wall_time_sec", "map_peak_rss_gb",
    "frac_mapped", "frac_correct", "frac_incorrectly_mapped", "frac_correct_in_repeat",
    "frac_correct_strict50", "frac_correct_strict50_in_repeat",
    "mean_identity",
)

# summary's correlation grid (below the table): Pearson r across each
# dataset's orders (one point per order = its mean over seeds), via Excel's
# own CORREL so the sheet stays live if values are edited.
CORR_X = ("density_factor", "avg_occurrences")
CORR_HIGHLIGHT_ABOVE = 0.8  # grid cells with r > this get a green fill
CORR_Y = ("total_wall_time_sec", "build_wall_time_sec", "build_peak_rss_gb", "map_wall_time_sec", "map_peak_rss_gb")

# A ratio < 1 (this order's mean below native's) means "better than native"
# for these -- fewer/faster/smaller/less-mis-mapped.
LOWER_IS_BETTER = frozenset({
    "density_factor", "avg_occurrences",
    "frac_incorrectly_mapped",
    "total_wall_time_sec",
    "build_wall_time_sec", "build_peak_rss_gb",
    "map_wall_time_sec", "map_peak_rss_gb",
})
# A ratio > 1 means "better than native" for these -- more correctly mapped,
# higher identity. Together with LOWER_IS_BETTER this must cover every
# SUMMARY_METRICS entry exactly once (asserted below), since the color scale
# needs an explicit direction for each -- unlike sshash, whose own metrics
# are all lower-is-better, this benchmark has real accuracy metrics too.
HIGHER_IS_BETTER = frozenset({
    "frac_mapped", "frac_correct", "frac_correct_in_repeat",
    "frac_correct_strict50", "frac_correct_strict50_in_repeat",
    "mean_identity",
})
assert LOWER_IS_BETTER | HIGHER_IS_BETTER == set(SUMMARY_METRICS)
assert not (LOWER_IS_BETTER & HIGHER_IS_BETTER)

ORDER_SORT = [
    "native", "xor-hash", "xor-pcg-hash", "greedymini",
    "mod-minimizer-xor-hash", "mod-minimizer-pcg-xor-hash",
]

# t-mer length for mod-minimizer's inner order O_t; fixed at 5 throughout
# this benchmark (see orders/README.md's r=4 formula: t=r+(k-r)%w=5 for
# k=15,w=10) -- not tracked per-row in summary_rows, so kept here instead of
# threading mod_t through that function just for the 'orders' sheet's sake.
_ORDER_T = {"mod-minimizer-xor-hash": 5, "mod-minimizer-pcg-xor-hash": 5}

ORDER_DESCRIPTIONS = {
    "native": {
        "description": "minimap2's own unmodified baseline order (mm_hash64-based)",
        "cli": "(no --minimizer-order)",
    },
    "xor-hash": {
        "description": "XOR the k-mer with the seed, then rank canonically",
        "cli": "--minimizer-order xor-hash --order-seed SEED",
    },
    "xor-pcg-hash": {
        "description": "xor-hash's XOR step, then PCG RXS-M-XS 64-bit scramble before ranking",
        "cli": "--minimizer-order xor-pcg-hash --order-seed SEED",
    },
    "greedymini": {
        "description": "Ranks k-mers by the precomputed GreedyMini w10_k15 rank table (binary alphabet, dual-orientation projection, no Extend-k)",
        "cli": "--minimizer-order greedymini --greedymini-file orders/GreedyMini/w10_k15_canonical.gm",
    },
    "mod-minimizer-xor-hash": {
        "description": "Mod-minimizer (WABI 2024): ranks t-mers by xor-hash, samples the k-mer at the winning t-mer's position mod w",
        "cli": "--minimizer-order mod-minimizer-xor-hash --order-seed SEED --mod-t 5",
    },
    "mod-minimizer-pcg-xor-hash": {
        "description": "Mod-minimizer: ranks t-mers by xor-pcg-hash, samples the k-mer at the winning t-mer's position mod w",
        "cli": "--minimizer-order mod-minimizer-pcg-xor-hash --order-seed SEED --mod-t 5",
    },
}

DATASET_ORDER_PREFERENCE = ("chm13_chrX", "grch38_wg")


def _dataset_sort_key(dataset: str) -> tuple:
    try:
        return (0, DATASET_ORDER_PREFERENCE.index(dataset))
    except ValueError:
        return (1, dataset)


def summary_datasets(summary: list[dict]) -> list[str]:
    return sorted({r["dataset"] for r in summary if r.get("dataset")}, key=_dataset_sort_key)


def _num(v, cast=float):
    if v is None or v == "":
        return None
    try:
        return cast(v)
    except (TypeError, ValueError):
        return None


DATASET_DESCRIPTIONS = {
    "chm13_chrX": "Real T2T CHM13 chrX v0.7 (Winnowmap's own D1: simulated PacBio reads, centromeric-repeat comparison)",
    "grch38_wg": "Real GRCh38 no-alt-analysis-set whole genome (Winnowmap's own D2: simulated PacBio reads); repeat region = a verified UCSC genomicSuperDups chr9 segmental duplication pair",
}

# 'smoke' is a fast-iteration/CI dataset (lib/make_dataset.py), not a source
# of real accuracy/density conclusions -- excluded from every exported sheet,
# same as sshash's own EXCLUDED_DATASETS.
EXCLUDED_DATASETS = frozenset({"smoke"})


def _kb_to_gb(kb) -> float | None:
    return round(kb / KB_PER_GB, 3) if kb is not None else None


def _total_time(build_sec, map_sec) -> float | None:
    """Index build + mapping wall time: the end-to-end cost of one run."""
    if build_sec is None or map_sec is None:
        return None
    return round(build_sec + map_sec, 2)


def flatten(meta: dict) -> dict:
    acc = meta.get("accuracy") or {}
    build = meta.get("build") or {}
    mapp = meta.get("map") or {}
    return {
        "dataset": meta.get("dataset"),
        "mid_occ_cap": meta.get("mid_occ_cap") or DEFAULT_MID_OCC_CAP,
        "mid_occ": meta.get("mid_occ"),
        "order": meta.get("order"),
        "seed": meta.get("seed"),
        "seed_idx": meta.get("seed_idx"),
        "k": meta.get("k"),
        "w": meta.get("w"),
        # run.sh passes --mod-t to every run, but it only means anything for
        # the mod-minimizer orders -- blank elsewhere rather than a misleading 5.
        "mod_t": meta.get("mod_t") if str(meta.get("order") or "").startswith("mod-minimizer-") else None,
        "distinct_minimizers": meta.get("distinct_minimizers"),
        "avg_occurrences": meta.get("avg_occurrences"),
        # raw-bp scale as written by write_meta.py; attach_hpc_density()
        # overwrites density_factor with the HPC-k-mer scale
        "density_factor": meta.get("density_factor"),
        "density_factor_raw_bp": meta.get("density_factor"),
        "total_wall_time_sec": _total_time(build.get("wall_time_sec"), mapp.get("wall_time_sec")),
        "build_wall_time_sec": build.get("wall_time_sec"),
        "build_peak_rss_gb": _kb_to_gb(build.get("peak_rss_kb")),
        "map_wall_time_sec": mapp.get("wall_time_sec"),
        "map_peak_rss_gb": _kb_to_gb(mapp.get("peak_rss_kb")),
        "n_reads": acc.get("n_reads"),
        "n_reads_in_repeat": acc.get("n_reads_in_repeat"),
        "frac_mapped": acc.get("frac_mapped"),
        "frac_correct": acc.get("frac_correct"),
        "frac_correct_in_repeat": acc.get("frac_correct_in_repeat"),
        "frac_incorrectly_mapped": acc.get("frac_incorrectly_mapped"),
        "frac_correct_strict50": acc.get("frac_correct_strict50"),
        "frac_correct_strict50_in_repeat": acc.get("frac_correct_strict50_in_repeat"),
        "mean_identity": acc.get("mean_identity"),
    }


def load_runs(results_dir: Path) -> list[dict]:
    rows = []
    for meta_path in sorted((results_dir / "runs").glob("**/meta.json")):
        meta = json.loads(meta_path.read_text())
        row = flatten(meta)
        if row.get("dataset") in EXCLUDED_DATASETS:
            continue
        if row["mid_occ"] is None:
            row["mid_occ"] = _map_log_mid_occ(meta_path.parent / "map.log")
        rows.append(row)
    return sorted(rows, key=lambda r: (*_block_sort_key(r), _order_rank(r["order"] or ""), int(r.get("seed_idx") or 0)))


def attach_hpc_density(runs: list[dict], hpc: dict[str, int]) -> None:
    """density_factor on the HPC-k-mer scale (see the module docstring):
    total occurrences = distinct_minimizers * avg_occurrences, both straight
    from minimap2 -d's own mm_idx_stat line. Datasets whose n_hpc_kmers is
    unknown keep the raw-bp value."""
    for r in runs:
        n_hpc = hpc.get(r.get("dataset"))
        distinct, avg_occ, w = _num(r.get("distinct_minimizers")), _num(r.get("avg_occurrences")), _num(r.get("w"), int)
        if n_hpc and distinct and avg_occ and w:
            r["density_factor"] = round(distinct * avg_occ * (w + 1) / n_hpc, 6)


def _map_log_mid_occ(map_log: Path) -> int | None:
    """Runs from before write_meta.py recorded mid_occ: read minimap2's own
    '[M::mm_mapopt_update::...] mid_occ = N' line instead."""
    if not map_log.is_file():
        return None
    m = re.search(r"mid_occ = (\d+)", map_log.read_text(errors="replace"))
    return int(m.group(1)) if m else None


def _block(row: dict) -> tuple:
    return (row.get("dataset"), row.get("mid_occ_cap"))


def _block_sort_key(row: dict) -> tuple:
    cap = row.get("mid_occ_cap")
    return (row.get("dataset") or "", cap != DEFAULT_MID_OCC_CAP, str(cap))


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------


def _mean_std(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return statistics.mean(values), statistics.stdev(values)


def _order_rank(order: str) -> int:
    try:
        return ORDER_SORT.index(order)
    except ValueError:
        return len(ORDER_SORT)


def summary_rows(runs: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for r in runs:
        key = (r.get("dataset"), r.get("k"), r.get("w"), r.get("mid_occ_cap"), r.get("order"))
        groups.setdefault(key, []).append(r)

    out = []
    for (dataset, k, w, cap, order), parts in groups.items():
        mids = [p["mid_occ"] for p in parts if p.get("mid_occ") is not None]
        row = {
            "dataset": dataset,
            "k": k,
            "w": w,
            "mid_occ_cap": cap,
            "order": order,
            # effective cap (mean over seeds); informational, no std/ratio
            "mid_occ_mean": round(statistics.mean(mids), 1) if mids else None,
            "n_seeds": len({p.get("seed_idx") for p in parts}),
            "n_reads": parts[0].get("n_reads"),
        }
        for metric in SUMMARY_METRICS:
            vals = [float(p[metric]) for p in parts if p.get(metric) is not None]
            mean, std = _mean_std(vals)
            row[f"{metric}_mean"] = round(mean, 6) if mean is not None else None
            row[f"{metric}_std"] = round(std, 6) if std is not None else None
        out.append(row)
    return sorted(out, key=lambda r: (*_block_sort_key(r), _order_rank(r["order"] or "")))


def attach_native_ratios(summary: list[dict]) -> None:
    """{metric}_ratio = this order's mean / native's mean for the same
    dataset (native = 1 exactly). native has no seed variance of its own, so
    it's an unambiguous, always-present baseline (unlike sshash, which needs
    a modality fallback -- this benchmark has no modality axis)."""
    native_by_block = {_block(r): r for r in summary if r.get("order") == "native"}
    for r in summary:
        base = native_by_block.get(_block(r))
        for metric in SUMMARY_METRICS:
            mean = r.get(f"{metric}_mean")
            base_mean = base.get(f"{metric}_mean") if base else None
            rname = f"{metric}_ratio"
            r[rname] = round(mean / base_mean, 4) if mean is not None and base_mean else None


def summary_field_order() -> list[str]:
    fields = ["dataset", "k", "w", "mid_occ_cap", "order", "n_seeds", "n_reads"]
    for metric in SUMMARY_METRICS:
        if metric == "map_wall_time_sec":
            fields.append("mid_occ_mean")
        fields.extend((f"{metric}_mean", f"{metric}_std", f"{metric}_ratio"))
    return fields


# ---------------------------------------------------------------------------
# datasets (read straight off disk, not hand-maintained)
# ---------------------------------------------------------------------------


def _faidx(samtools: str | None, genome: Path) -> Path | None:
    fai = genome.with_suffix(genome.suffix + ".fai")
    if fai.is_file():
        return fai
    if not samtools:
        return None
    try:
        subprocess.run([samtools, "faidx", str(genome)], check=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return fai if fai.is_file() else None


def _count_lines(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open("rb") as f:
        return sum(1 for _ in f)


def _repeat_bp(path: Path) -> tuple[int, int] | None:
    """(n_regions, total_bp) from a BED file."""
    if not path.is_file():
        return None
    n, total = 0, 0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        total += int(parts[2]) - int(parts[1])
        n += 1
    return n, total


def dataset_rows(data_dir: Path, samtools: str | None, dataset_names: set[str], k: int) -> list[dict]:
    out = []
    for name in sorted(dataset_names):
        d = data_dir / name
        genome = d / "genome.fa"
        row = {
            "dataset": name,
            "description": DATASET_DESCRIPTIONS.get(name, ""),
            "n_contigs": None,
            "genome_len_bp": None,
            "hpc_len_bp": None,
            "n_hpc_kmers": None,
            "n_reads": None,
            "n_repeat_regions": None,
            "repeat_bp": None,
        }
        fai = _faidx(samtools, genome) if genome.is_file() else None
        if fai:
            lines = fai.read_text().splitlines()
            row["n_contigs"] = len(lines)
            row["genome_len_bp"] = sum(int(line.split("\t")[1]) for line in lines if line.strip())
        if genome.is_file():
            hpc = hpc_kmers(genome, k)
            row["hpc_len_bp"], row["n_hpc_kmers"] = hpc["hpc_bp"], hpc["n_hpc_kmers"]
        n_reads = _count_lines(d / "truth.tsv")
        if n_reads is not None:
            row["n_reads"] = max(0, n_reads - 1)  # minus header
        rep = _repeat_bp(d / "repeats.bed")
        if rep:
            row["n_repeat_regions"], row["repeat_bp"] = rep
        out.append(row)
    return out


def dataset_field_order() -> list[str]:
    return [
        "dataset", "description", "n_contigs", "genome_len_bp", "hpc_len_bp", "n_hpc_kmers",
        "n_reads", "n_repeat_regions", "repeat_bp",
    ]


# ---------------------------------------------------------------------------
# orders -- one row per order (not per dataset x order): w/t/description/cli
# and expected (theoretical) density don't vary by dataset, so repeating
# them per dataset would be pure duplication (see expected_density.py's own
# docstring). Each dataset's particular (measured, mean/std over seeds)
# density_factor is pivoted into its own pair of columns instead -- same
# shape as sshash's own scripts/export.py:orders_rows.
# ---------------------------------------------------------------------------


def orders_rows(summary: list[dict], expected: dict[str, dict] | None = None) -> list[dict]:
    expected = expected or {}
    datasets = summary_datasets(summary)
    by_order: dict[str, dict[str, dict]] = {}
    for r in summary:
        order = r.get("order")
        if order and r.get("mid_occ_cap") == DEFAULT_MID_OCC_CAP:
            by_order.setdefault(order, {})[r.get("dataset")] = r

    out = []
    for order in ORDER_SORT:
        per_dataset = by_order.get(order)
        if not per_dataset:
            continue
        info = ORDER_DESCRIPTIONS.get(order, {"description": "", "cli": ""})
        exp = expected.get(order, {})
        w_val = next((r.get("w") for r in per_dataset.values() if r.get("w") is not None), None)
        k_val = next((r.get("k") for r in per_dataset.values() if r.get("k") is not None), None)
        row = {
            "order": order,
            "description": info["description"],
            "cli": info["cli"],
            "k": k_val,
            "w": w_val,
            "t": _ORDER_T.get(order),
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
    fields = ["order", "description", "cli", "k", "w", "t", "expected_density_factor_mean", "expected_density_factor_std"]
    for ds in summary_datasets(summary):
        fields.extend((f"{ds}_particular_density_factor_mean", f"{ds}_particular_density_factor_std"))
    return fields


# ---------------------------------------------------------------------------
# Column -> category, for the row-1 band (runs/summary only; datasets gets a
# single uniform category, matching sshash's own treatment of its simple
# reference sheets).
# ---------------------------------------------------------------------------

_ID_FIELDS = frozenset({
    "dataset", "order", "seed", "seed_idx", "k", "w", "mod_t", "mid_occ_cap", "n_seeds", "n_reads", "n_reads_in_repeat",
})
_METRIC_CATEGORY = {
    "distinct_minimizers": "density", "avg_occurrences": "density", "density_factor": "density",
    "density_factor_raw_bp": "density",
    "total_wall_time_sec": "total",
    "build_wall_time_sec": "build", "build_peak_rss_gb": "build",
    "map_wall_time_sec": "map", "map_peak_rss_gb": "map", "mid_occ": "map",
    "frac_mapped": "accuracy", "frac_correct": "accuracy",
    "frac_correct_in_repeat": "accuracy", "frac_incorrectly_mapped": "accuracy",
    "frac_correct_strict50": "accuracy", "frac_correct_strict50_in_repeat": "accuracy",
    "mean_identity": "accuracy",
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
    if field in ("order", "description", "cli", "k", "w", "t"):
        return "order"
    if field.startswith("expected_density_factor"):
        return "expected (theoretical)"
    marker = "_particular_density_factor_"
    if marker in field:
        return field.split(marker, 1)[0]
    return ""


# ---------------------------------------------------------------------------
# overview -- one row per column (or column family), sheet/column/meaning,
# same shape as sshash's own scripts/export.py:index_rows.
# ---------------------------------------------------------------------------

OVERVIEW_FIELDS = ["sheet", "column", "meaning"]


def overview_rows() -> list[dict]:
    def row(sheet: str, column: str, meaning: str) -> dict:
        return {"sheet": sheet, "column": column, "meaning": meaning}

    return [
        # --- id / params
        row("all", "dataset", "chm13_chrX = Winnowmap's D1 (CHM13 chrX v0.7), grch38_wg = Winnowmap's D2 (GRCh38 whole genome); both with PBSIM2-simulated PacBio reads (see datasets sheet)"),
        row("all", "k / w", "minimizer k-mer length and window size (k=15, w=10 -- minimap2's own map-pb/map-ont default window); HPC (homopolymer-compressed) k-mers on via -x map-pb"),
        row("summary/runs", "mid_occ_cap", "map-time repeat-filter condition: 'default' = minimap2 picks mid_occ itself per index; a number = fixed with -f <N> for every order (set to native's own default value, so all orders filter repeats by the same rule). Each (dataset, mid_occ_cap) is its own block: own native baseline, divider, correlation rows"),
        row("all", "order", "minimizer order: native (minimap2's own hash) | xor-hash | xor-pcg-hash | greedymini | mod-minimizer-xor-hash | mod-minimizer-pcg-xor-hash -- see orders sheet"),
        row("runs", "mod_t", "t-mer length of mod-minimizer's inner order (5); blank for non-mod orders, merged across each order's seeds"),
        row("runs", "seed / seed_idx", "order seed (XOR key / PCG tiebreak) and its index in config.yaml's shared seed list; native and greedymini are seed-independent, so run once (seed_idx 0)"),
        row("summary", "n_seeds", "number of seed runs aggregated into this row's mean/std (1 for native/greedymini, 5 otherwise)"),
        row("summary/runs", "n_reads / n_reads_in_repeat", "simulated reads in the dataset / of those, reads whose true origin overlaps repeats.bed (D1: centromeric DXZ1 array; D2: a chr9 segmental duplication pair)"),
        # --- density
        row("summary/runs/orders", "density_factor", "particular (measured) density on minimap2's own sketching scale: total minimizer occurrences (distinct_minimizers * avg_occurrences, from minimap2 -d) * (w+1) / n_hpc_kmers. -x map-pb turns on HPC, so minimap2 picks minimizers over homopolymer-compressed k-mers -- same scale as expected_density_factor; 1 = theoretical minimum, ~2 = random order"),
        row("runs", "density_factor_raw_bp", "same occurrences over the raw genome length instead ((w+1) / minimap2 -d's 'average spacing'): understates density by the genome's HPC compression ratio (hpc_len_bp / genome_len_bp, ~0.70 on human) and N gaps -- informational only"),
        row("runs", "distinct_minimizers", "distinct minimizer keys in the index (informational; undercounts selected positions whenever a minimizer repeats, so not used for density)"),
        row("summary/runs", "avg_occurrences", "mean reference positions per distinct minimizer (minimap2 -d's 'average occurrences'); higher = fewer, more repeated minimizers -- lower-density orders are higher here"),
        # --- total / build / map
        row("summary/runs", "total_wall_time_sec", "build_wall_time_sec + map_wall_time_sec: end-to-end cost of one run"),
        row("summary/runs", "build_wall_time_sec / build_peak_rss_gb", "index build (minimap2 -d -x map-pb -k 15 -w 10 -t 16): wall time and peak memory (GB = 1024^3 bytes), from /usr/bin/time -v"),
        row("summary/runs", "mid_occ (summary: mid_occ_mean)", "effective repeat-filter threshold minimap2 used when mapping (its own '[M::mm_mapopt_update] mid_occ = N' log line): any minimizer occurring > mid_occ times in the reference is skipped as a seed -- see 'mid_occ (algorithm)' below"),
        row("summary/runs", "map_wall_time_sec / map_peak_rss_gb", "mapping (minimap2 -x map-pb -c -t 16 [-f cap]): wall time and peak memory incl. the loaded index, from /usr/bin/time -v"),
        # --- accuracy
        row("summary/runs", "frac_mapped", "reads with any primary alignment / n_reads"),
        row("summary/runs", "frac_correct", "reads whose primary alignment is on the true contig and overlaps the true interval by >= 10% of their union / n_reads -- the Winnowmap / minimap2 papers' own criterion"),
        row("summary/runs", "frac_incorrectly_mapped", "mapped but not correct (frac_correct criterion) / n_reads"),
        row("summary/runs", "frac_correct_in_repeat", "frac_correct restricted to the n_reads_in_repeat reads (the Winnowmap paper's headline comparison)"),
        row("summary/runs", "frac_correct_strict50[_in_repeat]", "stricter variant: overlap >= 50% of the true span (overall / in-repeat)"),
        row("summary/runs", "mean_identity", "mean of 1 - de (gap-compressed per-base divergence from PAF's de:f: tag, dv:f: fallback) over primary alignments"),
        # --- summary aggregates
        row("summary", "{metric}_mean / {metric}_std", "mean / sample standard deviation over the row's n_seeds runs (std = 0 when n_seeds = 1)"),
        row("summary", "{metric}_ratio", "this order's mean / native's mean in the same (dataset, mid_occ_cap) block (native = 1, white). Colored green = better than native, red = worse; direction is per metric -- lower is better for density/avg_occurrences/time/memory/frac_incorrectly_mapped, higher is better for frac_mapped/frac_correct*/mean_identity"),
        row("summary", "correlation grid (below the table)", f"per (dataset, mid_occ_cap) block: Pearson r via Excel CORREL over that block's 6 order means, rows = density_factor_mean / avg_occurrences_mean, columns = total / build / map time and memory means; cells with r > {CORR_HIGHLIGHT_ABOVE} are green. Only 6 points per r -- a trend indicator, not a significance test"),
        # --- orders / datasets
        row("orders", "description / cli", "what the order does and the minimap2 flags that select it"),
        row("orders", "t", "mod-minimizer t-mer length (blank for non-mod orders)"),
        row("orders", "expected_density_factor_mean/std", "theoretical density_factor: Monte Carlo on i.i.d. random DNA (HPC-modeled), mean/std over the 5 seeds -- see lib/expected_density.py"),
        row("orders", "{dataset}_particular_density_factor_mean/std", "measured density_factor on that dataset (default-cap runs; density doesn't depend on the map-time cap)"),
        row("datasets", "n_contigs / genome_len_bp / n_reads / n_repeat_regions / repeat_bp", "read off disk: genome.fa.fai, truth.tsv, repeats.bed"),
        row("datasets", "hpc_len_bp / n_hpc_kmers", "genome length after homopolymer compression (N gaps excluded) / compressed k-mers minimap2 sketches (per N-free segment: hpc length - k + 1) -- lib/hpc_kmers.py, cached next to genome.fa"),
        # --- mid_occ algorithm
        row("mid_occ (algorithm)", "what it filters", "at map time each read minimizer is looked up in the index; each hit becomes a chaining anchor. A minimizer with more than mid_occ reference occurrences is dropped (seed.c mm_collect_matches), and one repeated > mid_occ times within the read itself and > 1% of its minimizers is dropped too (mm_seed_mz_flt). Bases covered by dropped minimizers count as rep_len, which lowers MAPQ"),
        row("mid_occ (algorithm)", "rescue", "if mid_occ < max_max_occ (4095, -U), inside each run of dropped minimizers the lowest-occurrence ones are kept, ~1 per 500 bp of read (--occ-dist, mm_seed_select); never above 4095. Active on D1 (mid_occ ~1244-1440); off on D2 (mid_occ ~9465-15000 > 4095: a hard cutoff)"),
        row("mid_occ (algorithm)", "default value", "without -f, mid_occ = occurrence count at the top 0.02% (-f 2e-4) of the index's distinct minimizers, clamped to [10, 1e6] (index.c mm_idx_cal_max_occ) -- so it's computed per index, per order"),
        row("mid_occ (algorithm)", "why it drifts with the order", "a lower-density order has fewer distinct minimizers, each occurring more often (higher avg_occurrences), so the same 0.02% quantile lands higher (D1: native 1244 vs greedymini 1417; D2: 9465 vs 13305). More repetitive hits then pass, adding chaining work without changing where reads map; fixing -f to native's value (mid_occ_cap) removes this confound"),
    ]


def sheet_categories(raw_name: str, fields: list[str]) -> list[str] | None:
    """One category label per field, for the row-1 band; None to skip it."""
    if raw_name == "overview":
        return ["reference"] * len(fields)
    if raw_name in ("runs", "summary"):
        return [_field_category(f) for f in fields]
    if raw_name == "orders":
        return [_orders_field_category(f) for f in fields]
    if raw_name == "datasets":
        return ["dataset reference"] * len(fields)
    return None


# ---------------------------------------------------------------------------
# XLSX (hand-rolled OOXML via benchmark/common/xlsx_util.py, no extra deps).
# Row 1 = category band (merged across contiguous same-category fields), row
# 2 = column headers, data from row 3 -- mirrors sshash's own
# two-row-header + vertical-merge + dataset-border + ratio-color-scale
# approach.
# ---------------------------------------------------------------------------

BAND_ROW = 1
HEADER_ROW = 2
DATA_START_ROW = 3

# Repeated across every row of a dataset's block (k/w are fixed for the whole
# benchmark; dataset obviously repeats) -- merge vertically instead of
# printing the same value on every row. 'order'/'seed'/'seed_idx' stay
# unmerged so each row's own identity is always visible.
XLSX_MERGE_FIELDS = frozenset({"dataset", "k", "w", "mid_occ_cap", "n_reads", "n_reads_in_repeat", "mod_t"})
# Merged only within one order's block (not just one dataset's): mod_t is
# per-order (blank for non-mod orders), so two adjacent mod orders sharing
# t=5 still each get their own merged cell.
XLSX_PER_ORDER_MERGE_FIELDS = frozenset({"mod_t"})
# Sheets with one row per (dataset, order[, seed]): get a divider (bottom
# border) after each dataset's last row, and the vertical merge above.
# 'datasets' is one row per dataset already, so it gets neither.
XLSX_DATASET_BORDER_SHEETS = frozenset({"summary", "runs"})
# Explicit column widths, for text-heavy sheets only.
XLSX_COL_WIDTHS = {"overview": (20, 42, 160)}


def _cell(value, col: int, row: int, *, style: int = 0) -> str:
    ref = f"{col_letter(col + 1)}{row}"
    s_attr = f' s="{style}"' if style else ""
    if value is None or value == "":
        return f"<c r=\"{ref}\"{s_attr}/>"
    if isinstance(value, bool):
        return f"<c r=\"{ref}\"{s_attr} t=\"b\"><v>{1 if value else 0}</v></c>"
    if isinstance(value, (int, float)):
        return f"<c r=\"{ref}\"{s_attr}><v>{value}</v></c>"
    text = xml_escape(str(value))
    return f"<c r=\"{ref}\"{s_attr} t=\"inlineStr\"><is><t xml:space=\"preserve\">{text}</t></is></c>"


def _dataset_last_row_indices(rows: list[dict]) -> set[int]:
    """0-indexed (relative to data start) row indices where a dataset's block
    ends -- where the divider border goes."""
    out: set[int] = set()
    for i, row in enumerate(rows):
        if i == len(rows) - 1 or _block(row) != _block(rows[i + 1]):
            out.add(i)
    return out


def _merge_rows_equal(rows: list[dict], i: int, j: int, field: str) -> bool:
    if rows[i].get("dataset") != rows[j].get("dataset"):
        return False
    if field != "dataset" and _block(rows[i]) != _block(rows[j]):
        return False
    if field in XLSX_PER_ORDER_MERGE_FIELDS:
        if rows[i].get("order") != rows[j].get("order") or rows[i].get(field) is None:
            return False
    return rows[i].get(field) == rows[j].get(field)


def _vertical_merge_plan(rows: list[dict], fields: list[str]) -> tuple[list[str], set[tuple[int, str]]]:
    """Return (mergeCell refs, skipped (row_idx, field) cells). Runs never
    cross a dataset boundary even if the value happens to repeat."""
    refs: list[str] = []
    skip: set[tuple[int, str]] = set()
    if len(rows) < 2:
        return refs, skip
    col = {name: i for i, name in enumerate(fields)}

    def close_run(field: str, run_start: int, run_end: int) -> None:
        if run_end <= run_start:
            return
        letter = col_letter(col[field] + 1)
        r1 = DATA_START_ROW + run_start
        r2 = DATA_START_ROW + run_end
        refs.append(f"{letter}{r1}:{letter}{r2}")
        for r in range(run_start + 1, run_end + 1):
            skip.add((r, field))

    for field in fields:
        if field not in XLSX_MERGE_FIELDS or field not in col:
            continue
        run_start = 0
        for i in range(1, len(rows)):
            if _merge_rows_equal(rows, i - 1, i, field):
                continue
            close_run(field, run_start, i - 1)
            run_start = i
        close_run(field, run_start, len(rows) - 1)
    return refs, skip


def _xlsx_merges_xml(refs: list[str]) -> str:
    if not refs:
        return ""
    cells = "".join(f'<mergeCell ref="{r}"/>' for r in refs)
    return f'<mergeCells count="{len(refs)}">{cells}</mergeCells>'


def _category_band_row(fields: list[str], categories: list[str]) -> tuple[str, list[str]]:
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
        cells.append(_cell(cat, start, BAND_ROW, style=2))
        if end > start:
            merges.append(f"{col_letter(start + 1)}{BAND_ROW}:{col_letter(end + 1)}{BAND_ROW}")
    return f'<row r="{BAND_ROW}">{"".join(cells)}</row>', merges


def _ratio_bounds(rows: list[dict], field: str) -> tuple[float, float]:
    vals = [v for r in rows if (v := r.get(field)) is not None]
    lo = min([1.0, *vals])
    hi = max([1.0, *vals])
    if lo == 1.0:
        lo = 0.99
    if hi == 1.0:
        hi = 1.01
    return lo, hi


def _color_scale_xml(fields: list[str], rows: list[dict]) -> str:
    """Color scale on every `{metric}_ratio` column: green = better than
    native (ratio 1, white), red = worse -- which side of 1 is green depends
    on the metric's own polarity (LOWER_IS_BETTER vs HIGHER_IS_BETTER), so
    an accuracy metric's color never gets flipped backwards from a
    time/density metric's."""
    if not rows:
        return ""
    data_start = DATA_START_ROW
    last_row = data_start + len(rows) - 1
    parts: list[str] = []
    green, red, white = "FF63BE7B", "FFF8696B", "FFFFFFFF"
    for priority, field in enumerate(fields, start=1):
        if not field.endswith("_ratio"):
            continue
        base = field[: -len("_ratio")]
        col_idx = fields.index(field)
        letter = col_letter(col_idx + 1)
        lo, hi = _ratio_bounds(rows, field)
        sqref = f"{letter}{data_start}:{letter}{last_row}"
        lo_color, hi_color = (red, green) if base in HIGHER_IS_BETTER else (green, red)
        parts.append(
            f'<conditionalFormatting sqref="{sqref}">'
            f'<cfRule type="colorScale" priority="{priority}"><colorScale>'
            f'<cfvo type="num" val="{lo}"/><cfvo type="num" val="1"/><cfvo type="num" val="{hi}"/>'
            f'<color rgb="{lo_color}"/><color rgb="{white}"/><color rgb="{hi_color}"/>'
            "</colorScale></cfRule></conditionalFormatting>"
        )
    return "".join(parts)


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    return sxy / (sxx * syy) ** 0.5


def _correlation_block(fields: list[str], rows: list[dict], start_row: int) -> tuple[list[str], list[str], str]:
    """summary's correlation grid, placed below the data table: per dataset,
    one row per CORR_X metric, one column per CORR_Y metric, each cell an
    Excel CORREL over that dataset's own row range (one point per order).
    The Python-computed Pearson r is stored as the formula's cached value so
    readers that don't recalculate still show it. Returns (row xml, merge
    refs, conditional-formatting xml)."""
    col = {name: i for i, name in enumerate(fields)}
    xs = [f"{m}_mean" for m in CORR_X if f"{m}_mean" in col]
    ys = [f"{m}_mean" for m in CORR_Y if f"{m}_mean" in col]
    if not xs or not ys or not rows:
        return [], [], ""

    blocks: list[tuple[tuple, int, int]] = []  # ((dataset, cap), first idx, last idx)
    for i, row in enumerate(rows):
        if blocks and blocks[-1][0] == _block(row):
            blocks[-1] = (blocks[-1][0], blocks[-1][1], i)
        else:
            blocks.append((_block(row), i, i))

    out: list[str] = []
    merges: list[str] = []
    r = start_row
    title = "correlation: Pearson r (CORREL) across each dataset's orders, one point per order's mean"
    out.append(f'<row r="{r}">{_cell(title, 0, r, style=2)}</row>')
    r += 1
    header = [_cell("dataset", 0, r, style=2), _cell("mid_occ_cap", 1, r, style=2), _cell("x \\ y", 2, r, style=2)]
    header += [_cell(y, 3 + j, r, style=2) for j, y in enumerate(ys)]
    out.append(f'<row r="{r}">{"".join(header)}</row>')
    r += 1
    grid_first = r
    for (dataset, cap), i1, i2 in blocks:
        r1, r2 = DATA_START_ROW + i1, DATA_START_ROW + i2
        for xi, x in enumerate(xs):
            last = xi == len(xs) - 1
            base, num = (1, 4) if last else (0, 3)
            cells = [
                _cell(dataset if xi == 0 else None, 0, r, style=base),
                _cell(cap if xi == 0 else None, 1, r, style=base),
                _cell(x, 2, r, style=base),
            ]
            xl = col_letter(col[x] + 1)
            xv = [_num(rows[i].get(x)) for i in range(i1, i2 + 1)]
            for j, y in enumerate(ys):
                yl = col_letter(col[y] + 1)
                yv = [_num(rows[i].get(y)) for i in range(i1, i2 + 1)]
                pairs = [(a, b) for a, b in zip(xv, yv) if a is not None and b is not None]
                rho = _pearson([a for a, _ in pairs], [b for _, b in pairs])
                ref = f"{col_letter(4 + j)}{r}"
                cached = f"<v>{round(rho, 6)}</v>" if rho is not None else ""
                cells.append(f'<c r="{ref}" s="{num}"><f>CORREL({xl}{r1}:{xl}{r2},{yl}{r1}:{yl}{r2})</f>{cached}</c>')
            out.append(f'<row r="{r}">{"".join(cells)}</row>')
            r += 1
        if len(xs) > 1:
            merges.append(f"A{r - len(xs)}:A{r - 1}")
            merges.append(f"B{r - len(xs)}:B{r - 1}")
    sqref = f"D{grid_first}:{col_letter(3 + len(ys))}{r - 1}"
    # No color scale here -- just flag strong positive correlation (dxfId 0 =
    # the green fill in _xlsx_styles_xml's <dxfs>).
    cf = (
        f'<conditionalFormatting sqref="{sqref}">'
        f'<cfRule type="cellIs" dxfId="0" priority="{len(fields) + 1}" operator="greaterThan">'
        f"<formula>{CORR_HIGHLIGHT_ABOVE}</formula></cfRule></conditionalFormatting>"
    )
    return out, merges, cf


def _sheet_xml(raw_name: str, fields: list[str], rows: list[dict]) -> str:
    xml = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    ]
    if fields:
        xml.append(
            '<sheetViews><sheetView workbookViewId="0">'
            f'<pane ySplit="{HEADER_ROW}" topLeftCell="A{DATA_START_ROW}" '
            'activePane="bottomLeft" state="frozen"/>'
            "</sheetView></sheetViews>"
        )
    if raw_name in XLSX_COL_WIDTHS:
        cols = "".join(
            f'<col min="{i}" max="{i}" width="{w}" customWidth="1"/>'
            for i, w in enumerate(XLSX_COL_WIDTHS[raw_name], start=1)
        )
        xml.append(f"<cols>{cols}</cols>")
    xml.append("<sheetData>")
    merge_refs: list[str] = []
    if fields:
        categories = sheet_categories(raw_name, fields)
        if categories:
            band_xml, band_merges = _category_band_row(fields, categories)
            xml.append(band_xml)
            merge_refs.extend(band_merges)

        header_cells = "".join(_cell(h, i, HEADER_ROW, style=2) for i, h in enumerate(fields))
        xml.append(f'<row r="{HEADER_ROW}">{header_cells}</row>')

        border_sheet = raw_name in XLSX_DATASET_BORDER_SHEETS
        vertical_refs, merge_skip = _vertical_merge_plan(rows, fields) if border_sheet else ([], set())
        merge_refs.extend(vertical_refs)
        dataset_last = _dataset_last_row_indices(rows) if border_sheet else set()

        for r_i, row in enumerate(rows, start=DATA_START_ROW):
            row_idx = r_i - DATA_START_ROW
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
                cells.append(_cell(value, c, r_i, style=style))
            xml.append(f'<row r="{r_i}">{"".join(cells)}</row>')
    corr_cf = ""
    if raw_name == "summary" and fields and rows:
        corr_rows, corr_merges, corr_cf = _correlation_block(fields, rows, DATA_START_ROW + len(rows) + 1)
        xml.extend(corr_rows)
        merge_refs.extend(corr_merges)
    xml.append("</sheetData>")
    if merge_refs:
        xml.append(_xlsx_merges_xml(merge_refs))
    if fields and rows:
        cs = _color_scale_xml(fields, rows)
        if cs:
            xml.append(cs)
    if corr_cf:
        xml.append(corr_cf)
    xml.append("</worksheet>")
    return "\n".join(xml)


def _xlsx_styles_xml() -> str:
    """Minimal styles: default, a bottom-border style (dataset divider), and
    a bold style (column headers)."""
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
        '<cellXfs count="5">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1"/>'
        '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
        # 3/4: correlation grid cells, "0.00" (built-in numFmt 2), plain / divider
        '<xf numFmtId="2" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
        '<xf numFmtId="2" fontId="0" fillId="0" borderId="1" xfId="0" applyNumberFormat="1" applyBorder="1"/>'
        "</cellXfs>"
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        # dxf 0: green fill, for the correlation grid's r > CORR_HIGHLIGHT_ABOVE rule
        '<dxfs count="1"><dxf><font><color rgb="FF006100"/></font>'
        '<fill><patternFill><bgColor rgb="FFC6EFCE"/></patternFill></fill></dxf></dxfs>'
        "</styleSheet>"
    )


def write_xlsx(path: Path, sheets: list[tuple[str, list[str], list[dict]]]) -> None:
    names, payloads = [], []
    for raw_name, fields, rows in sheets:
        names.append(raw_name[:31])
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


def _load_seeds(minimap2_bench_dir: Path) -> list[int]:
    """First 5 '  - 0x...' lines under config.yaml's 'seeds:' key -- same
    file lib/constants.sh's own load_bench_config() reads, kept in sync by
    reading that one file rather than duplicating the literal list."""
    cfg = minimap2_bench_dir / "config.yaml"
    fallback = [0x2545F4914F6CDD1D, 0x9E3779B97F4A7C15, 0xDEADBEEFCAFEBABE, 0x1234567890ABCDEF, 0xA5A5A5A5A5A5A5A5]
    if not cfg.is_file():
        return fallback
    seeds: list[int] = []
    in_seeds = False
    for line in cfg.read_text().splitlines():
        if line.startswith("seeds:"):
            in_seeds = True
            continue
        if in_seeds:
            stripped = line.strip()
            if line[:1] in (" ", "\t") and stripped.startswith("-"):
                val = stripped[1:].strip()
                seeds.append(int(val, 16) if val.lower().startswith("0x") else int(val))
            elif stripped:
                break
    return seeds or fallback


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--data-dir", type=Path, default=None, help="dataset root, for the 'datasets' sheet (default: $MINIMAP2_DATA_DIR or --results-dir/../data)")
    ap.add_argument("--samtools", default=None, help="samtools binary, for genome.fa.fai (falls back to $SAMTOOLS or 'samtools' on PATH)")
    ap.add_argument("--gm-file", type=Path, default=None, help="GreedyMini table, for the 'orders' sheet's expected density (default: orders/GreedyMini/w10_k15_canonical.gm)")
    args = ap.parse_args(argv)

    runs = load_runs(args.results_dir)
    if not runs:
        print(f"error: no completed runs in {args.results_dir / 'runs'}", file=sys.stderr)
        return 1

    import os
    data_dir = args.data_dir or Path(os.environ.get("MINIMAP2_DATA_DIR", args.results_dir.parent / "data"))
    samtools = args.samtools or os.environ.get("SAMTOOLS") or "samtools"
    dataset_names = {r["dataset"] for r in runs if r.get("dataset")}
    k_val = next((r.get("k") for r in runs if r.get("k") is not None), 15)
    datasets = dataset_rows(data_dir, samtools, dataset_names, int(k_val))
    attach_hpc_density(runs, {d["dataset"]: d["n_hpc_kmers"] for d in datasets})
    summary = summary_rows(runs)
    attach_native_ratios(summary)

    minimap2_bench_dir = Path(__file__).resolve().parents[1]
    repo_root = Path(__file__).resolve().parents[3]
    gm_file = args.gm_file or (repo_root / "orders" / "GreedyMini" / "w10_k15_canonical.gm")
    w_val = next((r.get("w") for r in runs if r.get("w") is not None), 10)
    mod_t_val = next((r.get("mod_t") for r in runs if r.get("mod_t")), 5)
    seeds = _load_seeds(minimap2_bench_dir)

    import expected_density
    expected = expected_density.ensure_expected_orders(
        args.results_dir, int(k_val), int(w_val), seeds,
        gm_file=gm_file, mod_t=int(mod_t_val),
    )
    orders = orders_rows(summary, expected)

    write_csv(args.results_dir / "minimap2_overview.csv", OVERVIEW_FIELDS, overview_rows())
    write_csv(args.results_dir / "minimap2_datasets.csv", dataset_field_order(), datasets)
    write_csv(args.results_dir / "minimap2_summary.csv", summary_field_order(), summary)
    write_csv(args.results_dir / "minimap2_runs.csv", RUN_FIELDS, runs)
    write_csv(args.results_dir / "minimap2_orders.csv", orders_field_order(summary), orders)
    print(f"wrote {len(datasets)} dataset row(s), {len(summary)} summary row(s), {len(runs)} run row(s), {len(orders)} order row(s) under {args.results_dir}")

    write_xlsx(args.results_dir / "minimap2_results.xlsx", [
        ("overview", OVERVIEW_FIELDS, overview_rows()),
        ("datasets", dataset_field_order(), datasets),
        ("summary", summary_field_order(), summary),
        ("orders", orders_field_order(summary), orders),
        ("runs", RUN_FIELDS, runs),
    ])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
