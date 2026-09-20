#!/usr/bin/env python3
"""Plot histogram figures for Kraken 2 scheme benchmark classify timings.

Reads the same run manifest as compare_benchmark.py and writes PNGs:
  classify_wall_s.png
  nextmin_s.png
  nextmin_ns_per_call.png
  lookup_s.png
  lookups_per_Mbp.png
  classify_vs_nextmin.png   (multi-panel summary)

Cold runs (classification_performance_cold.txt / nextminimizer_timing_cold.txt)
are excluded — only timed *_NN.txt repeats are plotted.

Invoked automatically by run_benchmark.sh after compare.

Usage:
  ./plot_benchmark_histograms.py --base-dir /path/to/benchmark
  ./plot_benchmark_histograms.py --base-dir ... --out-dir .../figures --bins 20
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
    SCHEME_ORDER,
    discover_classify_datasets,
    list_classify_time_files,
    list_nextmin_timing_files,
    load_accuracy_sheets,
    load_manifest,
    mean_std,
    parse_nextmin_timing_file,
    parse_time_file,
    summarize_argos_concordance,
)


SCHEME_COLORS = {
    "Random": "#4C78A8",
    "RandomXor": "#4C78A8",  # alias
    "GreedyMini": "#F58518",
    "GreedyMiniDNA": "#54A24B",
    "ModRandom": "#E45756",
    "ModRandomXor": "#E45756",  # alias
    "ModGreedyMini": "#72B7B2",
    "RandomPcg": "#B279A2",
}


def collect_repeat_series(runs, dataset="speed"):
    """Per scheme: lists of classify wall / nextmin metrics across all seeds×repeats."""
    metric_keys = (
        "class_s",
        "nextmin_s",
        "ns_per_call",
        "calls_per_Mbp",
        "scan_s",
        "lookup_s",
        "resolve_s",
        "lookups_per_Mbp",
        "repeats_per_Mbp",
        "lookup_frac_pct",
    )
    by_scheme = {}
    for run in runs:
        scheme = run["scheme"]
        if scheme not in by_scheme:
            slot = {k: [] for k in metric_keys}
            slot["contexts"] = []
            by_scheme[scheme] = slot
        slot = by_scheme[scheme]
        out = run["outdir"]
        for cf in list_classify_time_files(out, dataset):
            wall, _ = parse_time_file(cf)
            if wall is not None:
                slot["class_s"].append(wall)
                slot["contexts"].append(run["context"])
        for nf in list_nextmin_timing_files(out, dataset):
            nm = parse_nextmin_timing_file(nf)
            if not nm:
                continue
            for key in metric_keys:
                if key == "class_s":
                    continue
                if key in nm:
                    slot[key].append(nm[key])
    return by_scheme


def _scheme_order_key(name):
    return SCHEME_ORDER.get(name, 99), name


def _color(scheme):
    return SCHEME_COLORS.get(scheme, "#54A24B")


def _annotate_stats(ax, scheme, values, color, y_frac):
    mu, sd = mean_std(values)
    if mu is None:
        return
    label = "{}  n={}  mean={:.2f}".format(scheme, len(values), mu)
    if sd is not None and sd > 0:
        label += "±{:.2f}".format(sd)
    ax.axvline(mu, color=color, linestyle="--", linewidth=1.2, alpha=0.9)
    ax.text(
        0.98,
        y_frac,
        label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color=color,
    )


def plot_overlay_hist(by_scheme, key, xlabel, title, out_path, bins):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = []
    for scheme in sorted(by_scheme.keys(), key=_scheme_order_key):
        vals = by_scheme[scheme].get(key) or []
        if vals:
            series.append((scheme, vals))
    if not series:
        print("skip {}: no data".format(key), file=sys.stderr)
        return False

    all_vals = [v for _, vals in series for v in vals]
    # Shared bin edges; clip extreme outliers for readability (keep points, widen slightly)
    lo, hi = min(all_vals), max(all_vals)
    if hi <= lo:
        hi = lo + 1e-6

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for i, (scheme, vals) in enumerate(series):
        color = _color(scheme)
        ax.hist(
            vals,
            bins=bins,
            range=(lo, hi),
            alpha=0.45,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            label=scheme,
        )
        _annotate_stats(ax, scheme, vals, color, 0.96 - 0.08 * i)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count (classify repeats)")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote {}".format(out_path))
    return True


def plot_summary_panels(by_scheme, out_path, bins):
    """Multi-panel figure: wall, nextmin, lookup phase, density."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [
        ("class_s", "Classify wall time (s)", "Query wall time"),
        ("nextmin_s", "NextMinimizer time (s)", "Rank (NextMinimizer)"),
        ("lookup_s", "GetBatch lookup time (s)", "Hash lookup phase"),
        ("lookups_per_Mbp", "Distinct lookups / Mbp", "Minimizer density"),
    ]
    have = [
        (k, xl, t)
        for k, xl, t in panels
        if any(by_scheme[s].get(k) for s in by_scheme)
    ]
    if not have:
        return False

    ncols = len(have)
    fig, axes = plt.subplots(1, ncols, figsize=(4.8 * ncols, 4.2), squeeze=False)
    for ax, (key, xlabel, title) in zip(axes[0], have):
        series = []
        for scheme in sorted(by_scheme.keys(), key=_scheme_order_key):
            vals = by_scheme[scheme].get(key) or []
            if vals:
                series.append((scheme, vals))
        all_vals = [v for _, vals in series for v in vals]
        lo, hi = min(all_vals), max(all_vals)
        if hi <= lo:
            hi = lo + 1e-6
        for i, (scheme, vals) in enumerate(series):
            color = _color(scheme)
            ax.hist(
                vals,
                bins=bins,
                range=(lo, hi),
                alpha=0.45,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                label=scheme,
            )
            _annotate_stats(ax, scheme, vals, color, 0.96 - 0.08 * i)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Count")
        ax.set_title(title)
        ax.legend(frameon=False, loc="upper left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote {}".format(out_path))
    return True


def plot_argos_concordance(runs, out_path):
    """Grouped bars: Table S5 genus (a–d) per scheme (mean over seeds)."""
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

    labels = [
        "concordant",
        "discordant",
        "ancestor",
        "unclassified",
    ]
    keys = [
        "concordant_mean",
        "discordant_mean",
        "ancestor_mean",
        "unclassified_mean",
    ]
    schemes = [r["scheme"] or "?" for r in summary]
    n_sch = len(schemes)
    x = list(range(len(labels)))
    width = min(0.8 / max(n_sch, 1), 0.22)
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for i, rec in enumerate(summary):
        vals = [100.0 * float(rec[k] or 0.0) for k in keys]
        offset = (i - (n_sch - 1) / 2.0) * width
        ax.bar(
            [xi + offset for xi in x],
            vals,
            width,
            label=rec["scheme"] or "?",
            color=_color(rec["scheme"] or ""),
            edgecolor="white",
            linewidth=0.4,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fragments (%)")
    ax.set_title("FDA-ARGOS genus concordance (Table S5; no sensitivity/PPV)")
    ax.set_ylim(0, 105)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote {}".format(out_path))
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-dir", required=True, help="Benchmark base directory")
    ap.add_argument("--manifest", default=None)
    ap.add_argument(
        "--out-dir",
        default=None,
        help="Figure output directory (default: BASE_DIR/figures)",
    )
    ap.add_argument("--bins", type=int, default=20)
    args = ap.parse_args()

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print(
            "Error: matplotlib is required. Install with:\n"
            "  pip install matplotlib\n"
            "  # or: conda install matplotlib",
            file=sys.stderr,
        )
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

    datasets = []
    for r in runs:
        for d in discover_classify_datasets(r["outdir"]):
            if d not in datasets:
                datasets.append(d)
    if not datasets:
        datasets = ["speed"]

    os.makedirs(out_dir, exist_ok=True)
    for ds in datasets:
        by_scheme = collect_repeat_series(runs, dataset=ds)
        prefix = "" if ds == "speed" else ds + "_"
        tag = "" if ds == "speed" else " [{}]".format(ds)
        plot_overlay_hist(
            by_scheme,
            "class_s",
            "Classify wall time (s)",
            "Query-phase wall time distribution{}".format(tag),
            os.path.join(out_dir, "{}classify_wall_s.png".format(prefix)),
            args.bins,
        )
        plot_overlay_hist(
            by_scheme,
            "nextmin_s",
            "NextMinimizer time (s)",
            "NextMinimizer (rank calculation) time distribution{}".format(tag),
            os.path.join(out_dir, "{}nextmin_s.png".format(prefix)),
            args.bins,
        )
        plot_overlay_hist(
            by_scheme,
            "ns_per_call",
            "ns / NextMinimizer call",
            "Ranking cost per window{}".format(tag),
            os.path.join(out_dir, "{}nextmin_ns_per_call.png".format(prefix)),
            args.bins,
        )
        plot_overlay_hist(
            by_scheme,
            "lookup_s",
            "GetBatch lookup time (s)",
            "Hash-lookup phase time distribution{}".format(tag),
            os.path.join(out_dir, "{}lookup_s.png".format(prefix)),
            args.bins,
        )
        plot_overlay_hist(
            by_scheme,
            "lookups_per_Mbp",
            "Distinct lookups / Mbp",
            "Minimizer density (TOK_LOOKUP / Mbp){}".format(tag),
            os.path.join(out_dir, "{}lookups_per_Mbp.png".format(prefix)),
            args.bins,
        )
        plot_summary_panels(
            by_scheme,
            os.path.join(out_dir, "{}classify_vs_nextmin.png".format(prefix)),
            args.bins,
        )
    plot_argos_concordance(
        runs, os.path.join(out_dir, "fda-argos_genus_concordance.png")
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
