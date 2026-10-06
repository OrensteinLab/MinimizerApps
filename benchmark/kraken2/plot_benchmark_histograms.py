#!/usr/bin/env python3
"""Per-seed summary figures for the Kraken 2 scheme benchmark.

Every (scheme, seed, dataset, minimum-hit-groups) cell is classified exactly
once, so there are no timing histograms over repeats: the replicates are the
seeds. Reads the same run manifest as compare_benchmark.py and writes PNGs, each
one a grid (rows = --minimum-hit-groups, columns = dataset) with one dot per
seed and a bar at the mean over seeds, per scheme:

  classify_processing_s.png   classification time without the database load
  db_load_s.png               database load time (wall - processing)
  classify_wall_s.png         stacked mean: database load + processing
  classify_max_rss_gb.png     classify peak RSS
  fda-argos_genus_concordance.png   Table S5 genus concordant / discordant /
                              ancestor / unclassified (mean over seeds), one panel per
                              minimum-hit-groups
and, only when nextminimizer_timing_<dataset>_mhg<m>.txt files exist (classify
built with ENABLE_MINIMIZER_PROF=1):
  nextmin_s.png, nextmin_ns_per_call.png, lookup_s.png, lookups_per_Mbp.png

Invoked automatically by run_benchmark.sh after compare. matplotlib only.

Usage:
  ./plot_benchmark_histograms.py --base-dir /path/to/benchmark/runs/<library>
  ./plot_benchmark_histograms.py --base-dir ... --out-dir .../figures
"""

from __future__ import print_function

import argparse
import os
import sys

# Reuse parsers from compare_benchmark (same directory).
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from compare_benchmark import (  # noqa: E402
    DATASET_ORDER,
    SCHEME_ORDER,
    discover_classify_cells,
    load_accuracy_sheets,
    load_manifest,
    mean_std,
    nextmin_path,
    parse_nextmin_timing_file,
    parse_processing_s,
    parse_time_file,
    perf_path,
    summarize_argos_concordance,
)


SCHEME_COLORS = {
    "xor": "#4C78A8",
    "xor-pcg": "#B279A2",
    "greedymini": "#F58518",
    "mod-minimizer-xor": "#E45756",
    "mod-minimizer-greedymini": "#72B7B2",
}


def collect_cells(runs):
    """{(dataset, mhg): {scheme: {metric: [value per seed]}}} from the single
    classify run of every seed; metrics missing for a run are skipped."""
    cells = {}
    for run in runs:
        out = run["outdir"]
        for ds, mhg in discover_classify_cells(out):
            slot = cells.setdefault((ds, mhg), {}).setdefault(run["scheme"], {})
            cf = perf_path(out, ds, mhg)
            wall, rss = parse_time_file(cf)
            proc = parse_processing_s(cf)
            vals = {
                "wall": wall,
                "proc": proc,
                "db_load": max(0.0, wall - proc) if (wall is not None and proc is not None) else None,
                "rss": rss,
            }
            nf = nextmin_path(out, ds, mhg)
            nm = parse_nextmin_timing_file(nf) if nf else None
            if nm:
                vals["nextmin_s"] = nm.get("nextmin_s")
                vals["nextmin_ns"] = nm.get("ns_per_call")
                vals["lookup_s"] = nm.get("lookup_s")
                vals["lookups_pm"] = nm.get("lookups_per_Mbp")
            for k, v in vals.items():
                if v is not None:
                    slot.setdefault(k, []).append(v)
    return cells


def _scheme_order_key(name):
    return SCHEME_ORDER.get(name, 99), name


def _color(scheme):
    return SCHEME_COLORS.get(scheme, "#54A24B")


def _grid(cells):
    datasets = sorted(set(d for d, _m in cells), key=lambda d: (DATASET_ORDER.get(d, 9), d))
    mhgs = sorted(set(m for _d, m in cells))
    return datasets, mhgs


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_seed_grid(cells, key, ylabel, title, out_path):
    """Per-seed dots + mean bar per scheme; rows = minimum-hit-groups, cols = dataset."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not any(s.get(key) for cell in cells.values() for s in cell.values()):
        print("skip {}: no data".format(key), file=sys.stderr)
        return False
    datasets, mhgs = _grid(cells)
    schemes = sorted(set(s for cell in cells.values() for s in cell), key=_scheme_order_key)
    fig, axes = plt.subplots(
        len(mhgs), len(datasets), figsize=(3.6 * len(datasets) + 1.2, 3.0 * len(mhgs) + 1.4),
        squeeze=False, sharex=True,
    )
    for r, mhg in enumerate(mhgs):
        for c, ds in enumerate(datasets):
            ax = axes[r][c]
            cell = cells.get((ds, mhg), {})
            for i, sch in enumerate(schemes):
                vals = cell.get(sch, {}).get(key) or []
                if not vals:
                    continue
                n = len(vals)
                xs = [i + (j - (n - 1) / 2.0) * 0.08 for j in range(n)]
                ax.scatter(xs, vals, s=14, color=_color(sch), alpha=0.8, zorder=3)
                mu, _sd = mean_std(vals)
                ax.plot([i - 0.3, i + 0.3], [mu, mu], color=_color(sch), linewidth=2, zorder=2)
            ax.set_title("{}  (mhg {})".format(ds, mhg), fontsize=9)
            ax.set_xticks(list(range(len(schemes))))
            ax.set_xticklabels(schemes, rotation=40, ha="right", fontsize=7)
            if c == 0:
                ax.set_ylabel(ylabel)
            ax.set_ylim(bottom=0)
            _style(ax)
    fig.suptitle(title + "  (dots = seeds, bar = mean over seeds)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote {}".format(out_path))
    return True


def plot_wall_breakdown(cells, out_path):
    """Stacked bars of the mean DB-load and processing time per scheme."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not any(s.get("proc") for cell in cells.values() for s in cell.values()):
        print("skip classify_wall_s: no data", file=sys.stderr)
        return False
    datasets, mhgs = _grid(cells)
    schemes = sorted(set(s for cell in cells.values() for s in cell), key=_scheme_order_key)
    fig, axes = plt.subplots(
        len(mhgs), len(datasets), figsize=(3.6 * len(datasets) + 1.2, 3.0 * len(mhgs) + 1.4),
        squeeze=False, sharex=True,
    )
    for r, mhg in enumerate(mhgs):
        for c, ds in enumerate(datasets):
            ax = axes[r][c]
            cell = cells.get((ds, mhg), {})
            for i, sch in enumerate(schemes):
                d = mean_std(cell.get(sch, {}).get("db_load") or [])[0]
                p = mean_std(cell.get(sch, {}).get("proc") or [])[0]
                if p is None:
                    continue
                ax.bar(i, d or 0.0, color=_color(sch), alpha=0.45, edgecolor="white",
                       label="database load" if (i == 0 and r == 0 and c == 0) else None)
                ax.bar(i, p, bottom=d or 0.0, color=_color(sch), edgecolor="white",
                       label="processing" if (i == 0 and r == 0 and c == 0) else None)
            ax.set_title("{}  (mhg {})".format(ds, mhg), fontsize=9)
            ax.set_xticks(list(range(len(schemes))))
            ax.set_xticklabels(schemes, rotation=40, ha="right", fontsize=7)
            if c == 0:
                ax.set_ylabel("classify wall (s)")
            _style(ax)
    fig.suptitle("Classify wall = database load (light) + processing (dark); mean over seeds",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote {}".format(out_path))
    return True


def plot_argos_concordance(runs, out_path):
    """Grouped bars: Table S5 genus (a-d) per scheme (mean over seeds), one panel per
    minimum-hit-groups."""
    sheets = load_accuracy_sheets(runs)
    argos_rows = []
    for name, _fields, rows in sheets:
        if name == "accuracy_fda_argos":
            argos_rows = rows
            break
    if not argos_rows:
        return False
    _fields, summary = summarize_argos_concordance(argos_rows)
    if not summary:
        print("skip fda-argos concordance: no genus rows", file=sys.stderr)
        return False

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["concordant", "discordant", "ancestor", "unclassified"]
    keys = ["concordant_mean", "discordant_mean", "ancestor_mean", "unclassified_mean"]
    mhgs = sorted(set(r["minimum_hit_groups"] for r in summary))
    fig, axes = plt.subplots(1, len(mhgs), figsize=(7.2 * len(mhgs), 4.6), squeeze=False)
    for ax, mhg in zip(axes[0], mhgs):
        rows = [r for r in summary if r["minimum_hit_groups"] == mhg]
        n_sch = len(rows)
        x = list(range(len(labels)))
        width = min(0.8 / max(n_sch, 1), 0.22)
        for i, rec in enumerate(rows):
            vals = [100.0 * float(rec[k] or 0.0) for k in keys]
            offset = (i - (n_sch - 1) / 2.0) * width
            ax.bar([xi + offset for xi in x], vals, width, label=rec["scheme"] or "?",
                   color=_color(rec["scheme"] or ""), edgecolor="white", linewidth=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Fragments (%)")
        ax.set_title("FDA-ARGOS genus concordance (Table S5), minimum-hit-groups {}".format(mhg),
                     fontsize=9)
        ax.set_ylim(0, 105)
        ax.legend(frameon=False, fontsize=8)
        _style(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote {}".format(out_path))
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", required=True, help="Benchmark run root (runs/<library>)")
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--out-dir", default=None,
                    help="Figure output directory (default: BASE_DIR/figures)")
    args = ap.parse_args()

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("Error: matplotlib is required (pip install matplotlib / conda install matplotlib).",
              file=sys.stderr)
        return 1

    base = args.base_dir
    manifest = args.manifest or os.path.join(base, ".benchmark_run_manifest")
    out_dir = args.out_dir or os.path.join(base, "figures")
    if not os.path.isfile(manifest):
        print("Error: manifest not found: {}".format(manifest), file=sys.stderr)
        return 1
    runs = load_manifest(manifest)
    if not runs:
        print("Error: empty manifest", file=sys.stderr)
        return 1

    cells = collect_cells(runs)
    if not cells:
        print("Error: no classification_performance_<dataset>_mhg<N>.txt found", file=sys.stderr)
        return 1
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    plot_seed_grid(cells, "proc", "classify processing (s)",
                   "Classify processing time (database load excluded)",
                   os.path.join(out_dir, "classify_processing_s.png"))
    plot_seed_grid(cells, "db_load", "database load (s)", "Database load time",
                   os.path.join(out_dir, "db_load_s.png"))
    plot_wall_breakdown(cells, os.path.join(out_dir, "classify_wall_s.png"))
    plot_seed_grid(cells, "rss", "peak RSS (GB)", "Classify peak RSS",
                   os.path.join(out_dir, "classify_max_rss_gb.png"))
    # Phase figures only exist for ENABLE_MINIMIZER_PROF builds.
    for key, ylabel, title, fname in (
        ("nextmin_s", "NextMinimizer time (s)", "NextMinimizer (rank calculation) time", "nextmin_s.png"),
        ("nextmin_ns", "ns / NextMinimizer call", "Ranking cost per window", "nextmin_ns_per_call.png"),
        ("lookup_s", "GetBatch lookup time (s)", "Hash-lookup phase time", "lookup_s.png"),
        ("lookups_pm", "distinct lookups / Mbp", "Minimizer density (TOK_LOOKUP / Mbp)", "lookups_per_Mbp.png"),
    ):
        plot_seed_grid(cells, key, ylabel, title, os.path.join(out_dir, fname))
    plot_argos_concordance(runs, os.path.join(out_dir, "fda-argos_genus_concordance.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
