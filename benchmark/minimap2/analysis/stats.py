#!/usr/bin/env python3
"""Statistics for the scheme analysis (no scipy needed): reads the outputs of
run_seed_experiments.sh, accuracy_paired.py and lib/export.py under
--results-dir (csvs/analysis/<ds>/{seeds,timing}/ and csvs/analysis_*.csv) and
writes analysis_*.csv tables into <results-dir>/csvs.

- mechanism: per (dataset, cap, scheme), from the --dbg-seed-occ pass: read
  minimizers in the index per kb, fraction dropped by mid_occ, anchors per
  read (chaining input), reference occurrences before filtering, thread-local
  memory per read.
- timing: per (dataset, cap, scheme), total sample time per replicate; ratio to
  native; Welch t-test on log totals (Holm over the schemes of a block); 90% CI
  of the ratio and TOST equivalence at +-5% (equivalent iff the 90% CI lies in
  [0.95, 1/0.95]); paired per-read median log2 ratio with a bootstrap CI.
- per-read mechanism: within each scheme, Pearson r between log anchors and log
  per-read time across reads.
- accuracy: per (dataset, cap, scheme), from accuracy_paired.py's
  analysis_accuracy_scheme.csv: the paired correct-rate difference vs native
  with reads as the independent unit -- per read, d_i = (mean over the
  scheme's seeds of correct_i) - native correct_i (all seeds share one native
  run and one read set, so seeds are not independent replicates and are never
  pooled as n * n_seeds pairs); mean(d), SE = sd(d) / sqrt(n), its 90% CI
  (t, df = n - 1) and TOST at +-0.1 percentage point (equivalent iff that CI
  lies inside [-0.001, 0.001]); p from the paired t-test on d over reads, or
  the exact McNemar test for a single-run (seed-independent) scheme; Holm over
  the schemes of each (dataset, cap, subset) block. Datasets without truth
  (D3, real reads) have no accuracy rows.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from scheme_names import ID_FIELDS, parse_scheme_dir  # noqa: E402
from units import gb_from_bytes  # noqa: E402

# A "variant" is a scheme's run-directory name (scheme_names.scheme_dir_name):
# `native`, `xor`, `greedymini/k15_at1-cg0`, `mod-minimizer-xor/t9`, ... Its
# files under analysis/ use the name with "/" as "__"; every output row gets
# the separate identity columns (scheme, gm_k, projection, gm_table, mod_t).
CANONICAL = ["native", "xor", "xor-pcg-hash", "greedymini", "mod-minimizer-xor"]


def variant_key(variant: str) -> tuple:
    scheme = variant.split("/")[0]
    return (CANONICAL.index(scheme) if scheme in CANONICAL else len(CANONICAL), variant)


def id_cols(variant: str) -> dict:
    d = parse_scheme_dir(variant)
    return {f: d[f] for f in ID_FIELDS}


def variants_in(directory: Path, pattern: str) -> list[str]:
    """Variants having a file `<label>` matching `pattern` (a regex with one
    group for the label) in `directory`, in canonical order."""
    rx = re.compile(pattern)
    found = {m.group(1).replace("__", "/") for f in directory.iterdir() if (m := rx.fullmatch(f.name))}
    return sorted(found, key=variant_key)


def label(variant: str) -> str:
    return variant.replace("/", "__")
TIME_MARGIN = 0.05
ACC_MARGIN = 0.001


from stats_util import holm, t_quantile, t_two_sided_p, welch  # noqa: E402,F401  (shared, scipy-free)


# -- helpers --

def read_csv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else None


def r4(v):
    return None if v is None else round(v, 4)


# -- sections --

def mechanism(adir: Path, ds: str) -> list[dict]:
    rows = []
    for cap_dir in sorted((adir / ds / "seeds").glob("*")):
        for scheme in variants_in(cap_dir, r"(.+)\.csv"):
            f = cap_dir / f"{label(scheme)}.csv"
            rr = read_csv(f)
            kb = sum(int(r["qlen"]) for r in rr) / 1000
            seeds = sum(int(r["n_seeds"]) for r in rr)
            anchors = [int(r["anchors"]) for r in rr]
            rows.append({
                "dataset": ds, **id_cols(scheme), "mid_occ_cap": cap_dir.name, "n_reads": len(rr),
                "seeds_per_kb": round(seeds / kb, 2),
                "frac_seeds_filtered": round(sum(int(r["n_filtered"]) for r in rr) / seeds, 5),
                "anchors_per_read_mean": round(statistics.mean(anchors), 1),
                "anchors_per_read_median": statistics.median(anchors),
                "anchors_per_kb": round(sum(anchors) / kb, 1),
                "occ_before_filter_per_kb": round(sum(int(r["occ_all"]) for r in rr) / kb, 1),
                "km_cap_gb_mean": gb_from_bytes(statistics.mean(int(r["km_cap_bytes"]) for r in rr), 6),
                "km_cap_gb_max": gb_from_bytes(max(int(r["km_cap_bytes"]) for r in rr), 6),
            })
    base = {(r["mid_occ_cap"]): r for r in rows if r["scheme"] == "native"}
    for r in rows:
        b = base.get(r["mid_occ_cap"])
        r["anchors_ratio"] = round(r["anchors_per_kb"] / b["anchors_per_kb"], 4) if b else None
    return rows


def timing(adir: Path, ds: str, rng: random.Random) -> tuple[list[dict], list[dict]]:
    rows, per_read_rows = [], []
    for cap_dir in sorted((adir / ds / "timing").glob("*")):
        cap = cap_dir.name
        totals: dict[str, list[float]] = {}
        per_read: dict[str, dict[str, float]] = {}
        for scheme in variants_in(cap_dir, r"(.+)\.r\d+\.csv"):
            reps = sorted(cap_dir.glob(f"{label(scheme)}.r*.csv"))
            acc: dict[str, list[float]] = {}
            for f in reps:
                rr = read_csv(f)
                totals.setdefault(scheme, []).append(sum(float(r["time_sec"]) for r in rr))
                for r in rr:
                    acc.setdefault(r["read"], []).append(float(r["time_sec"]))
            per_read[scheme] = {k: statistics.mean(v) for k, v in acc.items()}
        if "native" not in totals:
            continue
        nat = [math.log(v) for v in totals["native"]]
        block = []
        for scheme in sorted(totals, key=variant_key):
            if scheme not in totals:
                continue
            lx = [math.log(v) for v in totals[scheme]]
            row = {
                "dataset": ds, **id_cols(scheme), "mid_occ_cap": cap, "n_reps": len(lx),
                "sample_wall_s_mean": round(statistics.mean(totals[scheme]), 2),
                "sample_wall_s_std": round(statistics.stdev(totals[scheme]), 2) if len(lx) > 1 else None,
            }
            if scheme != "native" and len(lx) > 1 and len(nat) > 1:
                d, se, df, p = welch(lx, nat)
                tq = t_quantile(0.95, df) if math.isfinite(df) else 1.645
                lo, hi = math.exp(d - tq * se), math.exp(d + tq * se)
                row.update({
                    "ratio": round(math.exp(d), 4), "ci90_lo": round(lo, 4), "ci90_hi": round(hi, 4),
                    "welch_p": p,
                    "tost_equivalent_5pct": lo >= 1 - TIME_MARGIN and hi <= 1 / (1 - TIME_MARGIN),
                })
                # paired per read: log2(scheme / native), averaged over replicates per read
                common = [k for k in per_read[scheme] if k in per_read["native"]]
                lr = [math.log2(per_read[scheme][k] / per_read["native"][k]) for k in common
                      if per_read[scheme][k] > 0 and per_read["native"][k] > 0]
                if lr:  # no reads timed in both runs: leave the per-read columns blank
                    med = statistics.median(lr)
                    boots = sorted(statistics.median(rng.choices(lr, k=len(lr))) for _ in range(500))
                    row.update({
                        "per_read_median_log2_ratio": round(med, 4),
                        "per_read_ci95_lo": round(boots[12], 4), "per_read_ci95_hi": round(boots[487], 4),
                    })
            block.append(row)
        tested = [r for r in block if "welch_p" in r]
        for r, adj in zip(tested, holm([r["welch_p"] for r in tested])):
            r["welch_p_holm"] = round(adj, 6)
            r["welch_p"] = round(r["welch_p"], 6)
        rows.extend(block)
        # per-read mechanism: log anchors vs log time, within scheme
        for scheme in per_read:
            sf = adir / ds / "seeds" / cap / f"{label(scheme)}.csv"
            if not sf.is_file():
                continue
            anchors = {r["read"]: int(r["anchors"]) for r in read_csv(sf)}
            xs, ys = [], []
            for k, t in per_read[scheme].items():
                if k in anchors and anchors[k] > 0 and t > 0:
                    xs.append(math.log(anchors[k])); ys.append(math.log(t))
            per_read_rows.append({"dataset": ds, **id_cols(scheme), "mid_occ_cap": cap,
                                  "n_reads": len(xs), "pearson_log_anchors_log_time": r4(pearson(xs, ys))})
    return rows, per_read_rows


def accuracy(results: Path) -> list[dict]:
    f = results / "csvs" / "analysis_accuracy_scheme.csv"
    if not f.is_file():
        if (results / "csvs" / "analysis_accuracy_paired.csv").is_file():
            print(f"warning: no {f.name}; rerun analysis/accuracy_paired.py (the per-seed "
                  "analysis_accuracy_paired.csv alone can't give a reads-as-unit test)")
        return []
    out = []
    for r in read_csv(f):
        if r["scheme"] == "native" or not r.get("diff_correct_rate") or not r.get("p"):
            continue  # native itself, or an empty subset (no in-repeat reads)
        d, se, df = float(r["diff_correct_rate"]), float(r["diff_se"]), int(r["df"])
        tq = t_quantile(0.95, df) if df > 0 else 1.645
        lo, hi = d - tq * se, d + tq * se
        out.append({
            "dataset": r["dataset"], **{f: r[f] for f in ID_FIELDS}, "mid_occ_cap": r["mid_occ_cap"], "subset": r["subset"],
            "n_seeds": int(r["n_seeds"]), "n_reads": int(r["n_reads"]),
            "diff_correct_rate_mean": round(d, 6), "diff_se": round(se, 8),
            "diff_ci90_lo": round(lo, 6), "diff_ci90_hi": round(hi, 6),
            "tost_equivalent_0.1pp": lo >= -ACC_MARGIN and hi <= ACC_MARGIN,
            "test": r["test"], "p": float(r["p"]),
        })
    for key in {(r["dataset"], r["mid_occ_cap"], r["subset"]) for r in out}:
        block = [r for r in out if (r["dataset"], r["mid_occ_cap"], r["subset"]) == key]
        for r, adj in zip(block, holm([r["p"] for r in block])):
            r["p_holm"] = round(adj, 6)
            r["p"] = round(r["p"], 6)
    return sorted(out, key=lambda r: (r["dataset"], r["mid_occ_cap"], r["subset"],
                                      CANONICAL.index(r["scheme"]) if r["scheme"] in CANONICAL else len(CANONICAL),
                                      r["gm_k"], r["projection"], r["mod_t"]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    # D1, D2, D3 (lib/ and README.md: Datasets); D3 (real reads, no truth) gets
    # mechanism/timing rows only -- accuracy_paired.py skips it.
    ap.add_argument("--datasets", default="chm13_chrX,grch38_wg,chm13_rel3_ont")
    args = ap.parse_args(argv)
    adir = args.results_dir / "csvs" / "analysis"
    (args.results_dir / "csvs").mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)
    mech, tim, prr = [], [], []
    for ds in args.datasets.split(","):
        if not (adir / ds).is_dir():
            continue
        mech += mechanism(adir, ds)
        t, p = timing(adir, ds, rng)
        tim += t
        prr += p
    write_csv(args.results_dir / "csvs" / "analysis_mechanism.csv", mech)
    write_csv(args.results_dir / "csvs" / "analysis_timing.csv", tim)
    write_csv(args.results_dir / "csvs" / "analysis_per_read_anchor_time.csv", prr)
    write_csv(args.results_dir / "csvs" / "analysis_accuracy_tests.csv", accuracy(args.results_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
