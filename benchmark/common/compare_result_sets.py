#!/usr/bin/env python3
"""Run-to-run stability of two result sets of one app (e.g. one made under an earlier measurement
protocol and one under the current one), to check whether a protocol change makes timings more stable.

  compare_result_sets.py --app {kmc,kraken2,minimap2,sshash} --set A=DIR --set B=DIR [--out-dir DIR]

DIR is a result set with csvs/<app>_runs.csv (an archived summary is enough; minimap2 is read
from minimap2_runs.csv in both sets so the contexts match). Per group of repeats (context,
metric, scheme variant; the repeats are the seeds) it reports n, median, CV = sd/mean, robust CV
= 1.4826 MAD/median and max/min. Seeds change the xor/pcg orders themselves, so their spread
includes real seed-to-seed variation; the baseline scheme (native, kmc default) is the same run
repeated, so its spread is pure measurement noise (column pure_noise).

Per metric it then pairs the groups present in both sets: median CV per set, share of groups with
CV > 5%, and a paired Wilcoxon signed-rank test of log CV (second set vs first). Because the
benchmark reports schemes as ratios to a baseline, the deciding figure is the precision of those
ratios: scheme_tests.py's 90% confidence interval of every scheme/baseline ratio is computed on
both sets, and its width (hi/lo) is compared the same way, together with how far the ratios of
the two sets agree and how many verdicts differ.

Writes <out-dir>/stability_groups_<app>.csv, stability_summary_<app>.csv and
stability_ratio_ci_<app>.csv (default out-dir: the second set's csvs/).
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "kraken2"))
import scheme_tests as ot  # noqa: E402
from compare_benchmark import wilcoxon_signed_rank  # noqa: E402  (exact for n <= 20)

CV_HIGH = 0.05
BASELINES = {"minimap2": "native", "sshash": "native", "kmc": "default"}


def load_set(app: str, results: Path):
    """(rows, baseline) of scheme_tests, minimap2 always from its runs CSV."""
    if app == "minimap2":
        return ot.load_minimap2_csv(results / "csvs" / "minimap2_runs.csv"), "native"
    return ot.load(app, results)


def groups_of(rows) -> dict[tuple, list[float]]:
    g = defaultdict(list)
    for r in rows:
        for metric, v in r["metrics"].items():
            if v is not None and v > 0:
                g[("/".join(r["ctx"]), metric, r["label"])].append(v)
    return g


def spread(vals: list[float]) -> dict:
    med = statistics.median(vals)
    out = {"n": len(vals), "median": med, "cv": None, "robust_cv": None, "max_min": max(vals) / min(vals)}
    if len(vals) >= 2:
        mean = statistics.mean(vals)
        out["cv"] = statistics.stdev(vals) / mean
        out["robust_cv"] = 1.4826 * statistics.median(abs(v - med) for v in vals) / med
    return out


def paired(a: list[float], b: list[float]) -> dict:
    """Second set (b) vs first (a) over paired positive values: Wilcoxon on log b - log a."""
    la, lb = [math.log(x) for x in a], [math.log(x) for x in b]
    w = wilcoxon_signed_rank(lb, la)
    return {"n_pairs": len(a), "n_second_lower": sum(1 for x, y in zip(la, lb) if y < x),
            "p_wilcoxon": None if w is None else w["p_two"]}


def fmt(x, nd=4):
    return "" if x is None else round(x, nd)


def compare(app: str, sets: list[tuple[str, Path]]):
    (n1, d1), (n2, d2) = sets
    loaded = {name: load_set(app, d) for name, d in sets}
    g1, g2 = groups_of(loaded[n1][0]), groups_of(loaded[n2][0])
    baseline = loaded[n1][1]

    group_rows = []
    for key in sorted(set(g1) | set(g2)):
        ctx, metric, label = key
        row = {"app": app, "context": ctx, "metric": metric, "scheme": label,
               "pure_noise": label == BASELINES.get(app)}
        for name, g in ((n1, g1), (n2, g2)):
            s = spread(g[key]) if key in g else {}
            for f in ("n", "median", "cv", "robust_cv", "max_min"):
                row[f"{f}_{name}"] = fmt(s.get(f), 6) if f != "n" else s.get(f, 0)
        group_rows.append(row)

    summary = []
    by_metric = defaultdict(list)
    for r in group_rows:
        if r[f"cv_{n1}"] not in ("", None) and r[f"cv_{n2}"] not in ("", None):
            by_metric[(r["metric"], False)].append(r)
            if r["pure_noise"]:
                by_metric[(r["metric"], True)].append(r)
    for (metric, noise_only), rs in sorted(by_metric.items()):
        a = [max(r[f"cv_{n1}"], 1e-6) for r in rs]
        b = [max(r[f"cv_{n2}"], 1e-6) for r in rs]
        p = paired(a, b)
        summary.append({"app": app, "metric": metric, "groups": "baseline only (pure noise)" if noise_only else "all schemes",
                        "n_groups": len(rs),
                        f"median_cv_{n1}": fmt(statistics.median(a)), f"median_cv_{n2}": fmt(statistics.median(b)),
                        f"share_cv_gt5_{n1}": fmt(sum(x > CV_HIGH for x in a) / len(a), 3),
                        f"share_cv_gt5_{n2}": fmt(sum(x > CV_HIGH for x in b) / len(b), 3),
                        f"n_{n2}_lower": p["n_second_lower"], "p_wilcoxon": fmt(p["p_wilcoxon"], 6)})

    # precision of the scheme/baseline ratios the benchmark reports
    t1 = {(r["context"], r["metric"], r["scheme"], r["vs"]): r
          for r in ot.run(app, d1, None, loaded=loaded[n1])[0] if r["ratio"] != ""}
    t2 = {(r["context"], r["metric"], r["scheme"], r["vs"]): r
          for r in ot.run(app, d2, None, loaded=loaded[n2])[0] if r["ratio"] != ""}
    ci_rows = []
    for key in sorted(set(t1) & set(t2)):
        a, b = t1[key], t2[key]
        ci_rows.append({"app": app, "context": key[0], "metric": key[1], "scheme": key[2], "vs": key[3],
                        f"ratio_{n1}": a["ratio"], f"ratio_{n2}": b["ratio"],
                        f"ci_width_{n1}": fmt(a["ci90_hi"] / a["ci90_lo"]),
                        f"ci_width_{n2}": fmt(b["ci90_hi"] / b["ci90_lo"]),
                        f"verdict_{n1}": a["verdict"], f"verdict_{n2}": b["verdict"]})
    by_metric = defaultdict(list)
    for r in ci_rows:
        by_metric[r["metric"]].append(r)
    for metric, rs in sorted(by_metric.items()):
        a = [r[f"ci_width_{n1}"] for r in rs]
        b = [r[f"ci_width_{n2}"] for r in rs]
        p = paired(a, b)
        gap = [abs(math.log(r[f"ratio_{n2}"] / r[f"ratio_{n1}"])) for r in rs]
        summary.append({"app": app, "metric": metric, "groups": "scheme/baseline ratio: 90% CI width (hi/lo)",
                        "n_groups": len(rs),
                        f"median_cv_{n1}": fmt(statistics.median(a)), f"median_cv_{n2}": fmt(statistics.median(b)),
                        f"share_cv_gt5_{n1}": fmt(sum(x > 1.10 for x in a) / len(a), 3),
                        f"share_cv_gt5_{n2}": fmt(sum(x > 1.10 for x in b) / len(b), 3),
                        f"n_{n2}_lower": p["n_second_lower"], "p_wilcoxon": fmt(p["p_wilcoxon"], 6),
                        "median_ratio_gap": fmt(statistics.median(gap)),
                        "verdicts_differ": sum(r[f"verdict_{n1}"] != r[f"verdict_{n2}"] for r in rs)})
    return group_rows, summary, ci_rows, baseline


def write(path: Path, rows: list[dict]):
    if not rows:
        return
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--app", required=True, choices=["kmc", "kraken2", "minimap2", "sshash"])
    ap.add_argument("--set", action="append", required=True, metavar="NAME=DIR",
                    help="a result set; give exactly two, the reference first")
    ap.add_argument("--out-dir", type=Path, default=None)
    a = ap.parse_args(argv)
    sets = [(s.split("=", 1)[0], Path(s.split("=", 1)[1])) for s in a.set]
    if len(sets) != 2:
        ap.error("give exactly two --set NAME=DIR")
    groups, summary, ci_rows, _ = compare(a.app, sets)
    out = a.out_dir or sets[1][1] / "csvs"
    out.mkdir(parents=True, exist_ok=True)
    write(out / f"stability_groups_{a.app}.csv", groups)
    write(out / f"stability_summary_{a.app}.csv", summary)
    write(out / f"stability_ratio_ci_{a.app}.csv", ci_rows)
    (n1, _), (n2, _) = sets
    print(f"{a.app}: {len(groups)} groups, {len(ci_rows)} scheme/baseline ratios in both sets -> {out}")
    print(f"  {'metric':<22} {'groups':<44} {'n':>4} {n1:>8} {n2:>8}  {'>thr ' + n1:>9} {'>thr ' + n2:>9} "
          f"{n2 + ' lower':>9} {'p':>8}")
    for r in summary:
        print(f"  {r['metric']:<22} {r['groups']:<44} {r['n_groups']:>4} {r[f'median_cv_{n1}']:>8} "
              f"{r[f'median_cv_{n2}']:>8}  {r[f'share_cv_gt5_{n1}']:>9} {r[f'share_cv_gt5_{n2}']:>9} "
              f"{r[f'n_{n2}_lower']:>9} {r['p_wilcoxon']:>8}"
              + (f"  ratio gap {r['median_ratio_gap']}, verdicts differ {r['verdicts_differ']}"
                 if 'verdicts_differ' in r else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
