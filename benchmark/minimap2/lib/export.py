#!/usr/bin/env python3
"""Aggregate every runs/**/meta.json into per-sheet CSVs and one XLSX,
mirroring sshash's/kmc's own export.py layout (a lighter version: no
buckets/space-breakdown sheets, no equivalence hard-fail gate -- minimap2's
benchmark doesn't have those concepts). Writes, under --results-dir (csvs/ subdir created here; workbook directly in it):

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
LOWER_IS_BETTER/HIGHER_IS_BETTER below; summary_sheet passes each ratio
column's polarity to benchmark/common/xlsx_book.write_book (ratio_fields) so
"green" always actually means "better than native" and never gets flipped for
an accuracy metric. The workbook layout (category band, merges, dividers,
colour scales, footer rows) is entirely xlsx_book's.

Order identity is five separate columns (benchmark/common/order_names.py):
order, gm_k, projection, gm_table, mod_t. greedymini runs with each AT/CG
projection are separate rows (and run directories).

datasets: one row per dataset, describing what it actually is (genome size,
read count/source, repeat region), read straight off disk (genome.fa.fai via
samtools faidx, truth.tsv, repeats.bed) rather than hand-maintained, so it
can't go stale.

summary: one row per (dataset, k, w, order, gm_k, projection, gm_table, mod_t) -- mean/std of every numeric
metric across seeds (n_seeds=1 for native, which doesn't vary with seed; see
lib/orders.sh's order_is_seed_dependent).

orders: one row per (order variant, k, w, L) -- greedymini once per projection, each order once per (k, w, L) the
datasets were run with (L = k + w - 1; layout: benchmark/common/orders_sheet.py) -- with the expected (theoretical,
for that row's own k, w, mod_t) density_factor, then per dataset its particular (measured) density_factor on the rows
of its own (k, w, mod_t); export fails if an applicable cell is empty. expected_density.py
computes the theoretical value via Monte Carlo simulation on i.i.d. random
DNA for that dataset's k, w, mod-minimizer t and HPC setting (k and t differ
by preset: map-pb k=19 t=9, map-ont k=15 t=5), reusing
benchmark/common/density.py's shared primitives -- the same ones order.cpp's
own C++ is built on; cached per (k, w, t, HPC) plus the GreedyMini table,
seeds and Monte Carlo sizes (expected_density.py cache_key), since it is independent of any
actual run.

density_factor is on minimap2's own sketching scale: density_factor = total
occurrences * (w+1) / n_sketch_kmers, where n_sketch_kmers counts the k-mers
minimap2 actually sketches (lib/sketch_kmers.py) -- homopolymer-compressed
ones for an HPC preset (map-pb, the default; see each dataset's preset in
lib/constants.sh:dataset_preset), raw ones otherwise (map-ont). That's the
same scale as the dataset's Monte Carlo expected value (the orders sheet's expected_*). write_meta.py's meta.json density_factor
instead divides by the *raw* length (minimap2 -d's "average spacing"),
understating HPC density by the genome's compression ratio (~0.70 on human)
and any N gaps; that's kept as density_factor_raw_bp in the runs sheet only.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from correlation import correlation_rows  # noqa: E402
from order_names import BASELINES, ID_FIELDS, make_id  # noqa: E402
from orders_sheet import EXPECTED_FIELDS, missing_cells, orders_fields, particular_fields, span_length, with_klw  # noqa: E402
from units import gb_from_bytes, gb_from_kib  # noqa: E402
from xlsx_book import DATA_START_ROW, write_book  # noqa: E402
from xlsx_util import col_letter, write_csv  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sketch_kmers import sketch_kmers  # noqa: E402

# minimap2 presets that turn on HPC (-H) k-mers (minimap2/options.c).
HPC_PRESETS = frozenset({"map-pb", "map10k", "ava-pb"})


def dataset_preset(dataset_dir: Path, runs: list[dict]) -> str:
    """The preset the dataset's runs recorded in meta.json (lib/run.sh), else
    lib/constants.sh:dataset_preset's rule: the data/<dataset>/preset file,
    else map-pb."""
    recorded = {r["preset"] for r in runs if r.get("dataset") == dataset_dir.name and r.get("preset")}
    if len(recorded) > 1:
        raise SystemExit(f"error: {dataset_dir.name} has runs with different presets: {sorted(recorded)}")
    if recorded:
        return recorded.pop()
    f = dataset_dir / "preset"
    return f.read_text().strip() if f.is_file() else "map-pb"

# Grouped so each row-1 category band is one contiguous span:
# id (dataset, order, app params, seeds) -> density -> total -> build -> map -> accuracy.
RUN_FIELDS = [
    "dataset", "n_reads", "n_reads_in_repeat", "k", "w", "L", "order", "gm_k", "projection", "gm_table", "mod_t", "mid_occ_threshold", "threads", "seed", "seed_idx",
    "density_factor", "density_factor_raw_bp", "distinct_minimizers", "avg_occurrences",
    "total_wall_s",
    "build_wall_s", "build_max_rss_gb", "index_size_gb",
    "map_wall_s", "map_max_rss_gb",
    "frac_mapped", "frac_correct", "frac_incorrectly_mapped", "frac_correct_in_repeat",
    "frac_correct_strict50", "frac_correct_strict50_in_repeat",
    "mean_identity",
]

# mid_occ_threshold (shown after the order columns): the effective repeat-filter threshold minimap2 mapped
# with. Internally every row also carries mid_occ_cap, the map-time `-f` condition
# (lib/run.sh's MID_OCC): this label for minimap2's own per-index default, else the fixed
# absolute cap. Every per-block computation (native ratio, divider, merge, correlation) is
# keyed on (dataset, mid_occ_cap) -- not displayed, since the threshold shows it: identical
# for every order in a fixed-cap block, per-index in a default block -- so a fixed-cap
# control never gets compared against the other condition's native.
DEFAULT_MID_OCC_CAP = "default"
MOD_ORDERS = frozenset({"mod-minimizer-xor"})


SUMMARY_METRICS = (
    "density_factor", "distinct_minimizers",
    "total_wall_s",
    "build_wall_s", "build_max_rss_gb", "index_size_gb",
    "map_wall_s", "map_max_rss_gb",
    "frac_mapped", "frac_correct", "frac_incorrectly_mapped", "frac_correct_in_repeat",
    "frac_correct_strict50", "frac_correct_strict50_in_repeat",
    "mean_identity",
)

# Rows right under the summary table (sshash's layout): for every
# `{metric}_ratio` column, its Pearson correlation (Excel's built-in CORREL, so
# it stays live and traceable) with these driver ratios, over all summary rows.
CORRELATION_DRIVERS = (
    ("correlation with density", "density_factor_ratio"),
    ("correlation with distinct minimizers", "distinct_minimizers_ratio"),
)

# A ratio < 1 (this order's mean below native's) means "better than native"
# for these -- fewer/faster/smaller/less-mis-mapped.
LOWER_IS_BETTER = frozenset({
    "density_factor", "distinct_minimizers",
    "frac_incorrectly_mapped",
    "total_wall_s",
    "build_wall_s", "build_max_rss_gb", "index_size_gb",
    "map_wall_s", "map_max_rss_gb",
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
    "native", "xor", "xor-pcg-hash", "greedymini", "mod-minimizer-xor",
]
PROJECTION_SORT = ["at0-cg1", "at1-cg0"]

ORDER_DESCRIPTIONS = {
    "native": {
        "description": "minimap2's own unmodified baseline order (mm_hash64-based)",
        "cli": "(no --minimizer-order)",
    },
    "xor": {
        "description": "XOR the k-mer with the seed, then rank canonically",
        "cli": "--minimizer-order xor --order-seed SEED",
    },
    "xor-pcg-hash": {
        "description": "xor's XOR step, then PCG RXS-M-XS 64-bit scramble before ranking",
        "cli": "--minimizer-order xor-pcg-hash --order-seed SEED",
    },
    "greedymini": {
        "description": "Ranks k-mers by the precomputed GreedyMini table gm_table (gm_k = its k) under the binary AT/CG projection (at0-cg1: A,T -> 0, C,G -> 1; at1-cg0: flipped, every table key complemented, as sshash's --gm-projection), centre-extended to k (k = gm_k + 2i: the table ranks the middle gm_k bases, then Sigma and Gamma of the whole k-mer break ties; greedymini_centre.h, as in kraken2/sshash)",
        "cli": "--minimizer-order greedymini --greedymini-file orders/GreedyMini/<gm_table>.gm --gm-projection <projection> --order-seed SEED",
    },
    "mod-minimizer-xor": {
        "description": "Strand-symmetric mod-minimizer (mod-sampling, WABI 2024, with the centre tie rule of Cologni & Pibiri as in SSHash; this repo): ranks t-mers by xor, t == k mod w; the minimal t-mer at x samples the k-mer at offset x mod w, a tie the offset closest to the window centre (two equally close ones: the left if the window <= its reverse complement), so a sequence and its reverse complement select the same loci and the scheme is forward",
        "cli": "--minimizer-order mod-minimizer-xor --order-seed SEED --mod-t T (T = 4 + (k-4) % w)",
    },
}

DATASET_ORDER_PREFERENCE = ("chm13_chrX", "grch38_wg", "chm13_rel3_ont")


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
    "chm13_chrX": "Winnowmap D1: T2T CHM13 chrX v0.7 with PBSIM 1.0.3 PacBio CLR reads (10% error, 15 kb mean, 100-25,000 bp, 5x), -x map-pb; repeat set = the DXZ1 centromeric array (the paper's Table 3)",
    "grch38_wg": "Winnowmap D2: GRCh38 with ALT contigs (UCSC hg38, 455 sequences) with the same PBSIM 1.0.3 reads at 0.5x, reads containing 'NN' dropped (paftools pbsim2fq), -x map-pb; repeat set = a UCSC genomicSuperDups chr9 segmental duplication pair (our addition)",
    "chm13_rel3_ont": "Winnowmap D3: 1,000 real CHM13 rel3 ultra-long ONT reads (uniform random sample, >= 1 kbp) mapped to GRCh38 with ALT contigs (D2's genome) with -x map-ont; no ground truth: unmapped %, time, memory and density only, as in the paper",
}

# 'smoke' is a fast-iteration/CI dataset (lib/make_dataset.py), not a source
# of real accuracy/density conclusions -- excluded from every exported sheet,
# same as sshash's own EXCLUDED_DATASETS.
EXCLUDED_DATASETS = frozenset({"smoke"})


def _kb_to_gb(kb) -> float | None:
    """peak_rss_kb (KiB, /usr/bin/time) -> decimal GB (1 GB = 10^9 bytes)."""
    return gb_from_kib(kb)


def _total_time(build_sec, map_sec) -> float | None:
    """Index build + mapping wall time: the end-to-end cost of one run."""
    if build_sec is None or map_sec is None:
        return None
    return round(build_sec + map_sec, 2)


def flatten(meta: dict) -> dict:
    try:
        idd = make_id(meta.get("order"), meta.get("gm_k"), meta.get("projection"), meta.get("gm_table"),
                      meta.get("mod_t"), baselines=BASELINES, default_projection=False)
    except ValueError as e:
        raise SystemExit(f"error: {e}")
    acc = meta.get("accuracy") or {}
    build = meta.get("build") or {}
    mapp = meta.get("map") or {}
    return {
        "dataset": meta.get("dataset"),
        "mid_occ_cap": meta.get("mid_occ_cap") or DEFAULT_MID_OCC_CAP,
        "mid_occ_threshold": meta.get("mid_occ"),
        # order, gm_k, projection, gm_table, mod_t (order_names.ID_FIELDS; only the
        # fields that apply to the order are set, e.g. mod_t only for mod-minimizer-*)
        **idd,
        "seed": meta.get("seed"),
        "seed_idx": meta.get("seed_idx"),
        "k": meta.get("k"),
        "w": meta.get("w"),
        "L": span_length(meta["k"], meta["w"]) if meta.get("k") is not None and meta.get("w") is not None else None,
        "preset": meta.get("preset"),  # absent in runs from before presets were recorded
        # minimap2 -t; absent in runs from before it was recorded (load_runs
        # then reads it off the build log's own CMD line)
        "threads": meta.get("threads"),
        "distinct_minimizers": meta.get("distinct_minimizers"),
        "avg_occurrences": meta.get("avg_occurrences"),
        # raw-bp scale as written by write_meta.py; attach_sketch_density()
        # overwrites density_factor with the sketched-k-mer scale
        "density_factor": meta.get("density_factor"),
        "density_factor_raw_bp": meta.get("density_factor"),
        "total_wall_s": _total_time(build.get("wall_time_sec"), mapp.get("wall_time_sec")),
        "build_wall_s": build.get("wall_time_sec"),
        "build_max_rss_gb": _kb_to_gb(build.get("peak_rss_kb")),
        "map_wall_s": mapp.get("wall_time_sec"),
        "map_max_rss_gb": _kb_to_gb(mapp.get("peak_rss_kb")),
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


def load_runs(results_dir: Path, datasets: set[str] | None = None) -> list[dict]:
    rows = []
    for meta_path in sorted((results_dir / "runs").glob("**/meta.json")):
        meta = json.loads(meta_path.read_text())
        row = flatten(meta)
        if row.get("dataset") in EXCLUDED_DATASETS or (datasets and row.get("dataset") not in datasets):
            continue
        if row["mid_occ_threshold"] is None:
            row["mid_occ_threshold"] = _map_log_mid_occ(meta_path.parent / "map.log")
        if row["threads"] is None:
            row["threads"] = _log_threads(meta_path.parent / "build.log")
        mmi = meta_path.parent / "index.mmi"
        row["index_size_gb"] = gb_from_bytes(mmi.stat().st_size) if mmi.is_file() else None
        rows.append(row)
    return sorted(rows, key=lambda r: (*_block_sort_key(r), _order_rank(r), int(r.get("seed_idx") or 0)))


def attach_sketch_density(runs: list[dict], n_kmers: dict[str, int]) -> None:
    """density_factor on the sketched-k-mer scale (see the module docstring):
    total occurrences = distinct_minimizers * avg_occurrences, both straight
    from minimap2 -d's own mm_idx_stat line. Datasets whose n_sketch_kmers is
    unknown keep the raw-bp value."""
    for r in runs:
        n_sk = n_kmers.get(r.get("dataset"))
        distinct, avg_occ, w = _num(r.get("distinct_minimizers")), _num(r.get("avg_occurrences")), _num(r.get("w"), int)
        if n_sk and distinct and avg_occ and w:
            r["density_factor"] = round(distinct * avg_occ * (w + 1) / n_sk, 6)


def _map_log_mid_occ(map_log: Path) -> int | None:
    """Runs from before write_meta.py recorded mid_occ: read minimap2's own
    '[M::mm_mapopt_update::...] mid_occ = N' line instead."""
    if not map_log.is_file():
        return None
    m = re.search(r"mid_occ = (\d+)", map_log.read_text(errors="replace"))
    return int(m.group(1)) if m else None


def _log_threads(log: Path) -> int | None:
    """Runs from before write_meta.py recorded threads: the -t on minimap2's
    own '[M::main] CMD: ...' log line."""
    if not log.is_file():
        return None
    m = re.search(r"\[M::main\] CMD:.*?\s-t\s*(\d+)", log.read_text(errors="replace"))
    return int(m.group(1)) if m else None


def _threads_text(runs: list[dict]) -> str:
    """The -t value(s) the exported runs used, for the overview text."""
    vals = sorted({int(r["threads"]) for r in runs if r.get("threads") is not None})
    return "/".join(str(v) for v in vals) if vals else "<threads>"


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


def _variant(row: dict) -> tuple:
    """The order identity (order_names.ID_FIELDS) of a run/summary row."""
    return tuple(row.get(f) for f in ID_FIELDS)


def _order_rank(row: dict) -> tuple:
    """Sort key: canonical order position, then gm_k, projection, gm_table, mod_t."""
    order = row.get("order") or ""
    try:
        o = ORDER_SORT.index(order)
    except ValueError:
        o = len(ORDER_SORT)
    proj = row.get("projection")
    p = PROJECTION_SORT.index(proj) if proj in PROJECTION_SORT else len(PROJECTION_SORT)
    return (o, row.get("gm_k") or 0, p, row.get("gm_table") or "", row.get("mod_t") or 0)


def summary_rows(runs: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for r in runs:
        key = (r.get("dataset"), r.get("k"), r.get("w"), r.get("mid_occ_cap"), *_variant(r))
        groups.setdefault(key, []).append(r)

    out = []
    for (dataset, k, w, cap, *ident), parts in groups.items():
        mids = [p["mid_occ_threshold"] for p in parts if p.get("mid_occ_threshold") is not None]
        row = {
            "dataset": dataset,
            "k": k,
            "w": w,
            "L": span_length(k, w) if k is not None and w is not None else None,
            "mid_occ_cap": cap,  # block key only (not a column, see RUN_FIELDS' note)
            **dict(zip(ID_FIELDS, ident)),
            # effective repeat-filter threshold (mean over seeds); no std/ratio
            "mid_occ_threshold": round(statistics.mean(mids), 1) if mids else None,
            "n_seeds": len({p.get("seed_idx") for p in parts}),
            "n_reads": parts[0].get("n_reads"),
        }
        for metric in SUMMARY_METRICS:
            vals = [float(p[metric]) for p in parts if p.get(metric) is not None]
            mean, std = _mean_std(vals)
            row[f"{metric}_mean"] = round(mean, 6) if mean is not None else None
            row[f"{metric}_std"] = round(std, 6) if std is not None else None
        out.append(row)
    return sorted(out, key=lambda r: (*_block_sort_key(r), _order_rank(r)))


def attach_native_ratios(summary: list[dict]) -> None:
    """{metric}_ratio = this order's mean / native's mean for the same
    (dataset, mid_occ_cap) block (native = 1 exactly). native has no seed variance of its own, so
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
    fields = ["dataset", "n_reads", "k", "w", "L", *ID_FIELDS, "mid_occ_threshold", "n_seeds"]
    for metric in SUMMARY_METRICS:
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


def dataset_params(runs: list[dict]) -> dict[str, tuple[int, int, int | None]]:
    """{dataset: (k, w, t)} from the runs (k/w differ by preset; t is the
    mod-minimizer t recorded by its runs)."""
    out: dict[str, tuple[int, int, int | None]] = {}
    for r in runs:
        ds, k, w = r.get("dataset"), _num(r.get("k"), int), _num(r.get("w"), int)
        if not ds or k is None or w is None:
            continue
        t = _num(r.get("mod_t"), int)
        prev = out.get(ds)
        if prev and (prev[0], prev[1]) != (k, w):
            raise SystemExit(f"error: {ds} has runs with different k/w: {prev[:2]} vs {(k, w)}")
        out[ds] = (k, w, t if t is not None else (prev[2] if prev else None))
    return out


def _count_fasta_reads(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open("rb") as f:
        return sum(1 for line in f if line.startswith(b">"))


def dataset_rows(data_dir: Path, samtools: str | None, runs: list[dict]) -> list[dict]:
    params = dataset_params(runs)
    out = []
    for name in sorted(params, key=_dataset_sort_key):
        k, w, t = params[name]
        d = data_dir / name
        genome = d / "genome.fa"
        row = {
            "dataset": name,
            "description": DATASET_DESCRIPTIONS.get(name, ""),
            "n_contigs": None,
            "genome_len_bp": None,
            "preset": dataset_preset(d, runs),
            "k": k,
            "w": w,
            "L": span_length(k, w),
            "mod_t": t,
            "sketch_len_bp": None,
            "n_sketch_kmers": None,
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
            sk = sketch_kmers(genome, k, hpc=row["preset"] in HPC_PRESETS)
            row["sketch_len_bp"], row["n_sketch_kmers"] = sk["sketch_bp"], sk["n_sketch_kmers"]
        n_reads = _count_lines(d / "truth.tsv")
        if n_reads is not None:
            row["n_reads"] = max(0, n_reads - 1)  # minus header
        else:  # real reads, no truth (the paper's D3)
            row["n_reads"] = _count_fasta_reads(d / "reads.fa")
        rep = _repeat_bp(d / "repeats.bed")
        if rep:
            row["n_repeat_regions"], row["repeat_bp"] = rep
        out.append(row)
    return out


def dataset_field_order() -> list[str]:
    return [
        "dataset", "description", "preset", "k", "w", "L", "mod_t", "n_contigs", "genome_len_bp", "sketch_len_bp", "n_sketch_kmers",
        "n_reads", "n_repeat_regions", "repeat_bp",
    ]


# ---------------------------------------------------------------------------
# orders -- one row per order; per dataset, the expected (theoretical, for that
# dataset's k, t and HPC setting -- see expected_density.py) and particular
# (measured, mean/std over seeds) density_factor side by side. Same shape as
# sshash's own scripts/export.py:orders_rows.
# ---------------------------------------------------------------------------


def _variant_cli(row: dict) -> str:
    """The minimap2 flags selecting one order variant (greedymini's table and
    projection are its own parameters, not part of the order name)."""
    info = ORDER_DESCRIPTIONS.get(row["order"], {"cli": ""})
    cli = info["cli"]
    if row["order"] == "greedymini" and row.get("gm_table"):
        cli = cli.replace("<gm_table>", row["gm_table"]).replace("<projection>", row.get("projection") or "at0-cg1")
    return cli


def dataset_combos(datasets: list[dict]) -> list[dict]:
    """Distinct (k, w, L, mod_t, hpc) combinations the datasets were run with,
    in dataset order, each with its datasets. mod_t = the datasets' recorded
    t, else 4 + (k-4) % w."""
    combos: dict[tuple, dict] = {}
    for d in datasets:
        k, w = int(d["k"]), int(d["w"])
        t = int(d["mod_t"]) if d.get("mod_t") is not None else 4 + (k - 4) % w
        hpc = d["preset"] in HPC_PRESETS
        c = combos.setdefault((k, w, t), {"k": k, "w": w, "L": span_length(k, w), "mod_t": t, "hpc": hpc, "datasets": []})
        if c["hpc"] != hpc:
            raise SystemExit(f"error: datasets of k={k} w={w} t={t} mix HPC and non-HPC presets")
        c["datasets"].append(d["dataset"])
    return list(combos.values())


def _plain_variant(row: dict) -> tuple:
    """(order, gm_k, projection, gm_table): an order variant without mod_t,
    which belongs to the (k, w, L) combination, not to the order."""
    return tuple(row.get(f) for f in ID_FIELDS if f != "mod_t")


def orders_rows(summary: list[dict], datasets: list[dict],
                expected: dict[tuple, dict[tuple, dict]] | None = None) -> list[dict]:
    """One row per (order variant, k, w, L) -- see benchmark/common/orders_sheet.py.
    datasets: dataset_rows(); expected: {(k, w, mod_t): {(order, projection or ''):
    expected_density row}}. A dataset's particular cells are filled on the rows of
    its own (k, w, mod_t); on other combinations they stay empty."""
    expected = expected or {}
    combos = dataset_combos(datasets)
    names = [d["dataset"] for d in datasets]
    measured: dict[tuple, dict] = {}   # (dataset, plain variant, mod_t or None) -> summary row
    variants: set[tuple] = set()
    for r in summary:
        if r.get("order") and r.get("mid_occ_cap") == DEFAULT_MID_OCC_CAP:
            pv = _plain_variant(r)
            variants.add(pv)
            t = r.get("mod_t") if r["order"] in MOD_ORDERS else None
            measured[(r.get("dataset"), pv, t)] = r

    def rank(pv):
        return _order_rank(dict(zip((f for f in ID_FIELDS if f != "mod_t"), pv)))

    out = []
    for pv in sorted(variants, key=rank):
        ident = dict(zip((f for f in ID_FIELDS if f != "mod_t"), pv))
        info = ORDER_DESCRIPTIONS.get(ident["order"], {"description": "", "cli": ""})
        seen_klw = set()
        for c in combos:
            if ident["order"] not in MOD_ORDERS:  # no mod_t: one row per (k, w, L)
                if (c["k"], c["w"]) in seen_klw:
                    continue
                seen_klw.add((c["k"], c["w"]))
            t = c["mod_t"]
            row = {"order": ident["order"], **{f: ident[f] for f in ("gm_k", "projection", "gm_table")}}
            with_klw(row, c["k"], c["w"])
            row["mod_t"] = t if ident["order"] in MOD_ORDERS else None
            row["description"], row["cli"] = info["description"], _variant_cli(ident)
            if ident["order"] in MOD_ORDERS:
                row["cli"] = row["cli"].replace("T (T = 4 + (k-4) % w)", f"{t}")
            exp = expected.get((c["k"], c["w"], t), {}).get((ident["order"], ident["projection"] or ""), {})
            row["expected_density_factor_mean"] = _num(exp.get("density_factor_mean"))
            row["expected_density_factor_std"] = _num(exp.get("density_factor_std"))
            for ds in names:
                r = measured.get((ds, pv, t if ident["order"] in MOD_ORDERS else None)) if ds in c["datasets"] else None
                mean = r.get("density_factor_mean") if r else None
                std = r.get("density_factor_std") if r else None
                if mean not in (None, "") and std in (None, ""):
                    std = 0.0
                row[f"{ds}_particular_density_factor_mean"] = mean
                row[f"{ds}_particular_density_factor_std"] = std
            out.append(row)
    return out


def check_orders_complete(rows: list[dict], datasets: list[dict]) -> None:
    """Refuse an orders sheet with holes in applicable cells: expected on every
    row, a dataset's particular cells on the rows of its own (k, w, mod_t)."""
    own = {}
    for c in dataset_combos(datasets):
        for ds in c["datasets"]:
            own[ds] = (c["k"], c["w"], c["mod_t"])
    bad = []
    for r in rows:
        fields = [*EXPECTED_FIELDS, *particular_fields([ds for ds, key in own.items()
                                                   if key[:2] == (r["k"], r["w"]) and (r["mod_t"] is None or key[2] == r["mod_t"])])]
        bad += [(r, f) for _, f in missing_cells([r], fields)]
    if bad:
        detail = "; ".join(f"{r['order']} {r.get('projection') or ''} k={r['k']} w={r['w']} t={r['mod_t']}: {f}" for r, f in bad[:12])
        raise SystemExit(f"error: orders sheet has {len(bad)} empty applicable cell(s): {detail}")


def orders_field_order(summary: list[dict]) -> list[str]:
    return list(orders_fields(summary_datasets(summary)))


# ---------------------------------------------------------------------------
# Column -> category, for the row-1 band (runs/summary only; datasets gets a
# single uniform category, matching sshash's own treatment of its simple
# reference sheets).
# ---------------------------------------------------------------------------

_ID_FIELDS = frozenset({
    "dataset", *ID_FIELDS, "seed", "seed_idx", "threads", "k", "w", "L", "mid_occ_threshold", "n_seeds", "n_reads", "n_reads_in_repeat",
})
_METRIC_CATEGORY = {
    "distinct_minimizers": "density", "avg_occurrences": "density", "density_factor": "density",
    "density_factor_raw_bp": "density",
    "total_wall_s": "total",
    "build_wall_s": "build", "build_max_rss_gb": "build", "index_size_gb": "build",
    "map_wall_s": "map", "map_max_rss_gb": "map",
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
    if field in (*ID_FIELDS, "k", "w", "L", "description", "cli", *EXPECTED_FIELDS):
        return "order"
    for marker in ("_expected_density_factor_", "_particular_density_factor_"):
        if marker in field:
            return field.split(marker, 1)[0]
    return ""


# ---------------------------------------------------------------------------
# overview -- one row per column (or column family), sheet/column/meaning,
# same shape as sshash's own scripts/export.py:index_rows.
# ---------------------------------------------------------------------------

OVERVIEW_FIELDS = ["sheet", "column", "meaning"]


def overview_rows(runs: list[dict] | None = None) -> list[dict]:
    threads = _threads_text(runs or [])

    def row(sheet: str, column: str, meaning: str) -> dict:
        return {"sheet": sheet, "column": column, "meaning": meaning}

    return [
        # --- id / params
        row("all", "dataset", "the Winnowmap paper's datasets (Jain et al. 2020): chm13_chrX = D1 (CHM13 chrX v0.7) and grch38_wg = D2 (GRCh38 whole genome), PacBio CLR reads simulated with PBSIM 1.0.3 as in the paper (10% error, 15 kb mean, lengths 100-25,000), mapped with -x map-pb; chm13_rel3_ont = D3, 1,000 real CHM13 rel3 ultra-long ONT reads (random sample, >= 1 kbp) mapped to GRCh38 with -x map-ont, no ground truth (see datasets sheet)"),
        row("summary/runs", "n_reads / n_reads_in_repeat", "reads in the dataset / of those, simulated reads whose true origin overlaps repeats.bed (chm13_chrX: the DXZ1 centromeric array, the paper's Table 3; grch38_wg: a chr9 segmental duplication pair, our addition; chm13_rel3_ont: none)"),
        row("all", "k / w / L", "minimizer length k, window size w and L = k + w - 1 (the span of a window; derived, right after w in every sheet that has k and w -- datasets, summary, runs, orders): each preset's own default, as the paper ran minimap2 (map-pb k=19, map-ont k=15; w=10); HPC (homopolymer-compressed) k-mers on for map-pb, off for map-ont (datasets sheet: preset, k, w, L)"),
        row("all", "order", "minimizer order (canonical names): native (minimap2's own hash, the baseline) | xor | xor-pcg-hash | greedymini | mod-minimizer-xor -- see orders sheet. The next four columns are the order's parameters, separate from its name; blank when they do not apply"),
        row("all", "gm_k / gm_table", "greedymini only: the GreedyMini table stem (orders/GreedyMini/<gm_table>.gm, e.g. w10_k15_canonical) and its k (gm_k = 15 here; at -k 19 the table is centre-extended to k = gm_k + 2i)"),
        row("all", "projection", "greedymini only: DNA -> binary projection, --gm-projection at0-cg1 (A,T -> 0, C,G -> 1, the default) or at1-cg0 (flipped: every table key complemented, as sshash). Each projection is a separate set of rows"),
        row("all", "mod_t", "mod-minimizer only: t-mer length of its inner order, t = 4 + (k-4) % w (9 at k=19, 5 at k=15); merged across each order's seeds"),
        row("summary/runs", "mid_occ_threshold", "effective repeat-filter threshold minimap2 mapped with (its own '[M::mm_mapopt_update] mid_occ = N' log line; summary: mean over seeds): any minimizer occurring more than this many times in the reference is skipped as a seed -- see 'mid_occ (algorithm)' below. In a default block (minimap2 picks it per index) it differs per order; in a fixed-cap block (-f <N> for every order, N = native's own default threshold) it is the same for all orders. Each such block is its own group of rows: own native baseline, divider; the ratios are computed within it"),
        row("runs", "threads", "minimap2 -t for build and map (meta.json; older runs: the build log's CMD line)"),
        row("runs", "seed / seed_idx", "order seed (XOR key; for greedymini the Gamma tie-break toggle) and its index in config.yaml's shared seed list; native is seed-independent, so run once (seed_idx 0)"),
        row("summary", "n_seeds", "number of seed runs aggregated into this row's mean/std (1 for native, 5 otherwise)"),
        # --- density
        row("summary/runs/orders", "density_factor", "particular (measured) density on minimap2's own sketching scale: total minimizer occurrences (distinct_minimizers * avg_occurrences, from minimap2 -d) * (w+1) / n_sketch_kmers (datasets sheet): the k-mers minimap2 actually sketches, homopolymer-compressed under -x map-pb, raw under -x map-ont -- same scale as the orders sheet's expected value of the matching kind; 1 = theoretical minimum, ~2 = random order"),
        row("runs", "density_factor_raw_bp", "same occurrences over the raw genome length instead ((w+1) / minimap2 -d's 'average spacing'): understates density by the genome's HPC compression ratio (sketch_len_bp / genome_len_bp under map-pb, ~0.70 on human) and N gaps -- informational only"),
        row("runs", "distinct_minimizers", "distinct minimizer keys in the index (informational; undercounts selected positions whenever a minimizer repeats, so not used for density)"),
        row("runs", "avg_occurrences", "mean reference positions per distinct minimizer (minimap2 -d's 'average occurrences'); higher = fewer, more repeated minimizers -- lower-density orders are higher here. Raw run value only (not aggregated in the summary)"),
        # --- total / build / map
        row("summary/runs", "total_wall_s", "build_wall_s + map_wall_s: end-to-end cost of one run"),
        row("summary/runs", "build_wall_s / build_max_rss_gb", f"index build (minimap2 -d -x <preset> -k <k> -w 10 -t {threads}; runs sheet: threads): wall time and peak memory (decimal GB, 1 GB = 10^9 bytes; /usr/bin/time's KiB converted via bytes), from /usr/bin/time -v"),
        row("summary/runs", "map_wall_s / map_max_rss_gb", f"mapping (minimap2 -x <preset> -c -t {threads} [-f cap], as the paper's -cx map-pb / -cx map-ont): wall time and peak memory incl. the loaded index, from /usr/bin/time -v"),
        # --- accuracy
        row("summary/runs", "frac_mapped", "reads with any primary alignment / n_reads (1 - frac_mapped = the paper's 'unmapped reads'); the only accuracy metric for chm13_rel3_ont (real reads, no truth), as in the paper"),
        row("summary/runs", "frac_correct", "reads whose primary alignment is on the true contig and overlaps the true interval by >= 10% of their union / n_reads -- the Winnowmap / minimap2 papers' own criterion"),
        row("summary/runs", "frac_incorrectly_mapped", "mapped but not correct (frac_correct criterion) / n_reads -- the paper's 'incorrectly mapped reads'; simulated datasets only"),
        row("summary/runs", "frac_correct_in_repeat", "frac_correct restricted to the n_reads_in_repeat reads (the Winnowmap paper's headline comparison)"),
        row("summary/runs", "frac_correct_strict50[_in_repeat]", "stricter variant: overlap >= 50% of the true span (overall / in-repeat)"),
        row("summary/runs", "mean_identity", "mean of 1 - de (gap-compressed per-base divergence from PAF's de:f: tag, dv:f: fallback) over primary alignments"),
        # --- summary aggregates
        row("summary", "{metric}_mean / {metric}_std", "mean / sample standard deviation over the row's n_seeds runs (std = 0 when n_seeds = 1)"),
        row("summary", "{metric}_ratio", "this order's mean / native's mean in the same block (same dataset and repeat-filter condition; native = 1, white). Colored green = better than native, red = worse; direction is per metric -- lower is better for density/distinct minimizers/time/memory/frac_incorrectly_mapped, higher is better for frac_mapped/frac_correct*/mean_identity"),
        row("summary", "correlation with density / distinct minimizers", "two rows right under the table (sshash's layout): for every {metric}_ratio column, Excel CORREL (Pearson) between that column and density_factor_ratio / distinct_minimizers_ratio over all summary rows (ratios are order / native within one block, so rows are comparable); blank for a constant column. No colour scale"),
        row("summary/runs", "index_size_gb", "size of the built index (.mmi, decimal GB = 10^9 bytes): 16 B per distinct key (hash table) + 8 B per occurrence of a repeated key (position array) + the 4-bit packed genome"),
        # --- orders / datasets
        row("orders", "description / cli", "what the order does and the minimap2 flags that select it; one row per (order variant, k, w, L) (greedymini once per projection)"),
        row("orders", "k / w / L / mod_t", "minimizer length k, window w (k-mers) and L = k + w - 1 (the span of a window), and the mod-minimizer t used at that (k, w, L), blank for the other orders (t = 4 + (k-4) % w; 9 at k=19, 5 at k=15). Each order appears once per (k, w, L) the datasets were run with"),
        row("orders", "expected_density_factor_mean/std", "theoretical density_factor for that row's own k, w, mod_t (HPC for map-pb datasets, raw for map-ont): Monte Carlo on i.i.d. random DNA, mean/std over the 5 seeds (std 0 for native) -- see lib/expected_density.py"),
        row("orders", "{dataset}_particular_density_factor_mean/std", "measured density_factor on that dataset (default-cap runs; density doesn't depend on the map-time cap); filled on the rows of the (k, w, L, mod_t) the dataset was run with (native std 0), empty on rows of another (k, w, L)"),
        row("datasets", "n_contigs / genome_len_bp / n_reads / n_repeat_regions / repeat_bp", "read off disk: genome.fa.fai, truth.tsv (reads.fa when there is no truth), repeats.bed"),
        row("datasets", "k / w / L / mod_t", "the dataset's k (its preset's default), w, L = k + w - 1 and mod-minimizer t, as used by every run on it"),
        row("datasets", "preset", "minimap2 -x preset for build and map (data/<dataset>/preset; default map-pb): map-pb = PacBio CLR with HPC k-mers, map-ont = ONT without"),
        row("datasets", "sketch_len_bp / n_sketch_kmers", "genome length minimap2 sketches (N gaps excluded; homopolymer-compressed under an HPC preset) / k-mers it sketches (per N-free segment: length - k + 1) -- lib/sketch_kmers.py, cached next to genome.fa"),
        # --- mid_occ algorithm
        row("mid_occ (algorithm)", "what it filters", "at map time each read minimizer is looked up in the index; each hit becomes a chaining anchor. A minimizer with more than mid_occ reference occurrences is dropped (seed.c mm_collect_matches), and one repeated > mid_occ times within the read itself and > 1% of its minimizers is dropped too (mm_seed_mz_flt). Bases covered by dropped minimizers count as rep_len, which lowers MAPQ"),
        row("mid_occ (algorithm)", "rescue", "if mid_occ < max_max_occ (4095, -U), inside each run of dropped minimizers the lowest-occurrence ones are kept, ~1 per 500 bp of read (--occ-dist, mm_seed_select); never above 4095. So the rescue is active when a dataset's mid_occ is below 4095 and the filter is a hard cutoff above it (see the mid_occ column)"),
        row("mid_occ (algorithm)", "default value", "without -f, mid_occ = occurrence count at the top 0.02% (-f 2e-4) of the index's distinct minimizers, clamped to [10, 1e6] (index.c mm_idx_cal_max_occ) -- so it's computed per index, per order"),
        row("mid_occ (algorithm)", "why it drifts with the order", "a lower-density order has fewer distinct minimizers, each occurring more often (higher avg_occurrences), so the same 0.02% quantile lands higher (compare mid_occ across orders in the default-cap block). More repetitive hits then pass, adding chaining work without changing where reads map; fixing -f to native's value (mid_occ_cap) removes this confound"),
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
# XLSX: sheets are plain dicts for benchmark/common/xlsx_book.write_book (row 1
# category band, row 2 headers, data from row 3, vertical merges, block
# dividers, ratio colour scales, footer rows).
# ---------------------------------------------------------------------------

# Repeated across every row of a (dataset, mid_occ_cap) block -- merged
# vertically instead of printing the same value on every row. 'order' and
# 'seed'/'seed_idx' stay unmerged so each row's own identity is always visible.
XLSX_MERGE_FIELDS = ["dataset", "n_reads", "n_reads_in_repeat", "k", "w", "L"]
# The order parameters are merged only within one order variant's rows (same
# order, gm_k, projection, gm_table, mod_t), so two adjacent variants that
# happen to share a value still each get their own merged cell.
XLSX_PARAM_FIELDS = ["gm_k", "projection", "gm_table", "mod_t"]
# summary/runs: a divider after each (dataset, mid_occ_cap) block; 'datasets'
# and 'orders' are one row per dataset/order already, so no blocks there.
XLSX_BLOCK_FIELDS = ["dataset", "mid_occ_cap"]


def _ratio_polarity(fields: list[str]) -> dict[str, str]:
    """{`{metric}_ratio`: 'lower' | 'higher'}: which side of native (1) is green."""
    out = {}
    for f in fields:
        if f.endswith("_ratio"):
            out[f] = "higher" if f[: -len("_ratio")] in HIGHER_IS_BETTER else "lower"
    return out


def _sheet(name: str, fields: list[str], rows: list[dict], **opts) -> dict:
    cats = sheet_categories(name, fields)
    sheet = {"name": name, "fields": fields, "rows": rows}
    if cats:
        sheet["categories"] = cats
    sheet.update(opts)
    return sheet


def _block_opts(fields: list[str]) -> dict:
    """Block/merge options of the summary and runs sheets: blocks by
    (dataset, mid_occ_cap); the dataset cell spans its caps; the order
    parameters merge within one order variant's rows."""
    variant = ["dataset", "mid_occ_cap", *ID_FIELDS]
    scope = {"dataset": ["dataset"]}
    scope.update({f: variant for f in XLSX_PARAM_FIELDS})
    return {
        "block_fields": XLSX_BLOCK_FIELDS,
        "merge_fields": [f for f in (*XLSX_MERGE_FIELDS, *XLSX_PARAM_FIELDS) if f in fields],
        "merge_scope": scope,
    }


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
    ap.add_argument("--datasets", default=None, help="comma-separated datasets to export (default: all)")
    ap.add_argument("--gm-file", type=Path, default=None, help="GreedyMini table, for the 'orders' sheet's expected density (default: orders/GreedyMini/<gm_table>.gm of the dataset's greedymini runs, else w10_k15_canonical)")
    args = ap.parse_args(argv)

    runs = load_runs(args.results_dir, set(args.datasets.split(",")) if args.datasets else None)
    if not runs:
        print(f"error: no completed runs in {args.results_dir / 'runs'}", file=sys.stderr)
        return 1

    import os
    data_dir = args.data_dir or Path(os.environ.get("MINIMAP2_DATA_DIR", args.results_dir.parent / "data"))
    samtools = args.samtools or os.environ.get("SAMTOOLS") or "samtools"
    datasets = dataset_rows(data_dir, samtools, runs)
    attach_sketch_density(runs, {d["dataset"]: d["n_sketch_kmers"] for d in datasets})
    summary = summary_rows(runs)
    attach_native_ratios(summary)

    minimap2_bench_dir = Path(__file__).resolve().parents[1]
    repo_root = Path(__file__).resolve().parents[3]
    gm_dir = repo_root / "orders" / "GreedyMini"
    seeds = _load_seeds(minimap2_bench_dir)

    import expected_density
    expected: dict[tuple, dict[tuple, dict]] = {}
    for c in dataset_combos(datasets):
        tables = sorted({r["gm_table"] for r in runs if r["dataset"] in c["datasets"] and r.get("gm_table")})
        gm_file = args.gm_file or gm_dir / f"{tables[0] if tables else 'w10_k15_canonical'}.gm"
        expected[(c["k"], c["w"], c["mod_t"])] = expected_density.ensure_expected_orders(
            args.results_dir, c["k"], c["w"], seeds, gm_file=gm_file, mod_t=c["mod_t"], is_hpc=c["hpc"],
        )
    orders = orders_rows(summary, datasets, expected)
    check_orders_complete(orders, datasets)

    csvs_dir, xlsx_dir = args.results_dir / "csvs", args.results_dir
    csvs_dir.mkdir(parents=True, exist_ok=True)
    xlsx_dir.mkdir(parents=True, exist_ok=True)
    write_csv(csvs_dir / "minimap2_overview.csv", OVERVIEW_FIELDS, overview_rows(runs))
    write_csv(csvs_dir / "minimap2_datasets.csv", dataset_field_order(), datasets)
    write_csv(csvs_dir / "minimap2_summary.csv", summary_field_order(), summary)
    write_csv(csvs_dir / "minimap2_runs.csv", RUN_FIELDS, runs)
    write_csv(csvs_dir / "minimap2_orders.csv", orders_field_order(summary), orders)
    print(f"wrote {len(datasets)} dataset row(s), {len(summary)} summary row(s), {len(runs)} run row(s), {len(orders)} order row(s) under {csvs_dir}")

    # sheet order as sshash's export: legend, orders, datasets, summary, runs
    summary_fields = summary_field_order()
    write_book(xlsx_dir / "minimap2_results.xlsx", [
        _sheet("overview", OVERVIEW_FIELDS, overview_rows(runs),
               widths={"sheet": 20, "column": 42, "meaning": 160}),
        _sheet("orders", orders_field_order(summary), orders),
        _sheet("datasets", dataset_field_order(), datasets),
        _sheet("summary", summary_fields, summary, **_block_opts(summary_fields),
               ratio_fields=_ratio_polarity(summary_fields), ratio_baseline=1.0,
               footer_rows=correlation_rows(summary_fields, summary, CORRELATION_DRIVERS)),
        _sheet("runs", RUN_FIELDS, runs, **_block_opts(RUN_FIELDS)),
    ])
    print(f"wrote {xlsx_dir / 'minimap2_results.xlsx'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
