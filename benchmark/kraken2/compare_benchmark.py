#!/usr/bin/env python3
"""Compare Kraken 2 minimizer-scheme benchmarks across random-order seeds.

Manifest columns (tab-separated):
  scheme  context  outdir  [seed]  [seed_idx]  [n_classify]

Classification may be repeated per seed; timing files are
  classification_performance_01.txt … _NN.txt
A preceding cold run may exist as classification_performance_cold.txt
(and nextminimizer_timing_cold.txt); those are excluded from stats.
When multiple repeats are present, Class(s)/Mrd/min are mean±std within the seed.
When multiple seeds are present, scheme summary is mean±std over seeds
(using each seed's classify mean).

NextMinimizer / phase instrumentation (when classify built with
ENABLE_MINIMIZER_PROF=1 → -DKRAKEN2_MINIMIZER_PROF) writes:
  nextminimizer_timing_01.txt … _NN.txt
  (key=value: nextmin_s, scan_s, lookup_s, resolve_s,
   lookups, repeats, lookups_per_Mbp, lookup_frac_pct,
   uhs, non_uhs, uhs_frac_pct, non_uhs_frac_pct, …)
  Path is set via env KRAKEN2_MINIMIZER_PROF.
  uhs = GM Rank finite pack (toggle/tb used); non_uhs = ∞ early-out.

Paired metrics (MAPE, concordance, ClassConc, Jaccard, NewClass/NewUnclass)
compare each seed to Random with the same seed_idx when available
(using last classify outputs).

Section 4 reports exact paired Wilcoxon signed-rank tests (by seed_idx) on
Class wall, Rank ns/call, and lookups/Mbp — no scipy required.
Seed = toggle mask for all schemes (Random order; GreedyMini Γ / Mod tb XOR).

Metrics follow Wood et al. Genome Biol 2019 (PMC6883579).
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from bench_const import REF_SCHEME, SCHEME_ORDER  # noqa: E402
from xlsxutil import load_table, write_tables, _col_letter  # noqa: E402


DATASET_ORDER = {
    "strain-exclusion": 0,
    "simulated": 0,
    "labeled": 0,
    "fda-argos": 1,
    "argos": 1,
    "speed": 2,
}
DOMAIN_ORDER = {"all": 0, "prokaryote": 1, "virus": 2, "viral": 2}
RANK_ORDER = {"genus": 0, "species": 1}
SCOPE_ORDER = {"experiment": 0, "unweighted_mean": 1, "pooled": 2}


def parse_time_file(path):
    if not os.path.isfile(path):
        return None, None
    wall_s = None
    rss_gb = None
    with open(path) as f:
        text = f.read()
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)", text)
    if m:
        parts = m.group(1).strip().split(":")
        try:
            if len(parts) == 2:
                wall_s = float(parts[0]) * 60 + float(parts[1])
            elif len(parts) == 3:
                wall_s = (
                    float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
                )
        except ValueError:
            wall_s = None
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if m:
        rss_gb = int(m.group(1)) / (1024.0 * 1024.0)
    return wall_s, rss_gb


def parse_classify_counts(path):
    if not os.path.isfile(path):
        return None, None
    n_reads = None
    n_class = None
    with open(path) as f:
        for line in f:
            m = re.search(r"(\d+)\s+sequences\s+\(", line)
            if m and "processed" in line:
                n_reads = int(m.group(1))
            m = re.search(r"(\d+)\s+sequences classified", line)
            if m:
                n_class = int(m.group(1))
    return n_reads, n_class


def parse_nextmin_timing_file(path):
    """Parse nextminimizer_timing_XX.txt (key=value) or stderr timing lines."""
    if not os.path.isfile(path):
        return None
    data = {}
    with open(path) as f:
        text = f.read()
    keys = (
        "nextmin_s",
        "calls",
        "ns_per_call",
        "calls_per_Mbp",
        "scan_s",
        "lookup_s",
        "resolve_s",
        "lookups",
        "repeats",
        "skips",
        "lookups_per_Mbp",
        "repeats_per_Mbp",
        "lookup_frac_pct",
        "uhs",
        "non_uhs",
        "uhs_frac_pct",
        "non_uhs_frac_pct",
        "bases",
        "sequences",
    )
    for key in keys:
        m = re.search(r"(?:^|\n)" + key + r"=([0-9.eE+-]+)", text)
        if m:
            try:
                data[key] = float(m.group(1))
            except ValueError:
                pass
    # Fallback: stderr lines embedded in classification_performance_*.txt
    if "nextmin_s" not in data:
        m = re.search(
            r"NEXTMIN_TIMING\s+nextmin_s=([0-9.eE+-]+)\s+calls=(\d+)\s+"
            r"ns_per_call=([0-9.eE+-]+)\s+calls_per_Mbp=([0-9.eE+-]+)"
            r"(?:\s+bases=\d+\s+sequences=\d+)?"
            r"(?:\s+uhs=(\d+)\s+non_uhs=(\d+)\s+"
            r"uhs_frac_pct=([0-9.eE+-]+)\s+non_uhs_frac_pct=([0-9.eE+-]+))?",
            text,
        )
        if m:
            data["nextmin_s"] = float(m.group(1))
            data["calls"] = float(m.group(2))
            data["ns_per_call"] = float(m.group(3))
            data["calls_per_Mbp"] = float(m.group(4))
            if m.group(5) is not None:
                data["uhs"] = float(m.group(5))
                data["non_uhs"] = float(m.group(6))
                data["uhs_frac_pct"] = float(m.group(7))
                data["non_uhs_frac_pct"] = float(m.group(8))
    if "scan_s" not in data:
        m = re.search(
            r"PHASE_TIMING\s+scan_s=([0-9.eE+-]+)\s+lookup_s=([0-9.eE+-]+)\s+"
            r"resolve_s=([0-9.eE+-]+)\s+lookups=(\d+)\s+repeats=(\d+)\s+"
            r"skips=(\d+)\s+lookups_per_Mbp=([0-9.eE+-]+)\s+"
            r"repeats_per_Mbp=([0-9.eE+-]+)\s+lookup_frac_pct=([0-9.eE+-]+)",
            text,
        )
        if m:
            data["scan_s"] = float(m.group(1))
            data["lookup_s"] = float(m.group(2))
            data["resolve_s"] = float(m.group(3))
            data["lookups"] = float(m.group(4))
            data["repeats"] = float(m.group(5))
            data["skips"] = float(m.group(6))
            data["lookups_per_Mbp"] = float(m.group(7))
            data["repeats_per_Mbp"] = float(m.group(8))
            data["lookup_frac_pct"] = float(m.group(9))
    if "nextmin_s" not in data and "scan_s" not in data:
        return None
    return data


def parse_kraken_report(path):
    rows = {}
    if not os.path.isfile(path):
        return 0, rows
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                pct = float(parts[0].strip())
                clade = int(parts[1])
                taxon = int(parts[2])
            except ValueError:
                continue
            rows[parts[4].strip()] = {
                "pct": pct,
                "clade": clade,
                "taxon": taxon,
                "rank": parts[3].strip(),
                "name": parts[5].strip(),
            }
    total = 0
    if "0" in rows:
        total += rows["0"]["clade"]
    if "1" in rows:
        total += rows["1"]["clade"]
    if total <= 0 and rows:
        total = max(r["clade"] for r in rows.values())
    return total, rows


def abundance_at_rank(rows, total, rank_letter):
    abund = {}
    if total <= 0:
        return abund
    for taxid, row in rows.items():
        if row["rank"] == rank_letter:
            abund[taxid] = float(row["clade"]) / float(total)
    if "0" in rows:
        abund["0"] = float(rows["0"]["clade"]) / float(total)
    return abund


def mape(reference, estimate, min_frac):
    keys = [k for k, v in reference.items() if v >= min_frac]
    if not keys:
        return None, 0
    err = 0.0
    for k in keys:
        t = reference[k]
        err += abs(t - estimate.get(k, 0.0)) / t
    return 100.0 * err / len(keys), len(keys)


def l1_percent(reference, estimate):
    keys = set(reference) | set(estimate)
    return 50.0 * sum(
        abs(reference.get(k, 0.0) - estimate.get(k, 0.0)) for k in keys
    )


def concordance_vs_ref(ref_path, other_path):
    n = same = ref_c = other_c = 0
    both_c = both_same = new_c = new_u = 0
    ref_cls_ids = set()
    other_cls_ids = set()
    with open(ref_path) as fr, open(other_path) as fo:
        for lr, lo in zip(fr, fo):
            pr = lr.rstrip("\n").split("\t")
            po = lo.rstrip("\n").split("\t")
            if len(pr) < 3 or len(po) < 3:
                continue
            n += 1
            ref_cls = pr[0] == "C"
            other_cls = po[0] == "C"
            if ref_cls:
                ref_c += 1
                ref_cls_ids.add(pr[1])
            if other_cls:
                other_c += 1
                other_cls_ids.add(po[1])
            if ref_cls and other_cls:
                both_c += 1
                if pr[2] == po[2]:
                    both_same += 1
            elif other_cls:
                new_c += 1
            elif ref_cls:
                new_u += 1
            if pr[0] == po[0] and pr[2] == po[2]:
                same += 1
    inter = len(ref_cls_ids & other_cls_ids)
    union = len(ref_cls_ids | other_cls_ids)
    return {
        "n": n,
        "concordance": 100.0 * same / n if n else 0.0,
        "classified_concordance": 100.0 * both_same / both_c if both_c else 0.0,
        "jaccard": 100.0 * inter / union if union else 0.0,
        "newly_classified_pct": 100.0 * new_c / n if n else 0.0,
        "newly_unclassified_pct": 100.0 * new_u / n if n else 0.0,
        "classified_pct_ref": 100.0 * ref_c / n if n else 0.0,
        "classified_pct": 100.0 * other_c / n if n else 0.0,
    }


def mean_std(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return None, None
    n = len(vals)
    mu = sum(vals) / float(n)
    if n < 2:
        return mu, 0.0
    var = sum((x - mu) ** 2 for x in vals) / float(n - 1)
    return mu, math.sqrt(var)


def wilcoxon_signed_rank(a, b):
    """Exact paired Wilcoxon signed-rank test on paired samples a[i] vs b[i].

    Zero differences are dropped (Wilcoxon). Average ranks used for tied |diff|.
    For n≤20 non-zero pairs, p-values are exact (enumerate 2^n sign assignments).
    For n>20, a normal approximation with continuity correction is used.

    Returns dict with n, mean_diff (a-b), n_lt (count a<b), W_plus, W_minus,
    p_two, p_less (P(W+ ≤ obs) supporting a<b), p_greater; or None if n<1.
    """
    if len(a) != len(b) or not a:
        return None
    diffs = [x - y for x, y in zip(a, b)]
    nz = [(abs(d), 1 if d > 0 else -1) for d in diffs if d != 0]
    n = len(nz)
    if n < 1:
        return None

    # Average ranks of |diff|
    order = sorted(range(n), key=lambda i: nz[i][0])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and nz[order[j]][0] == nz[order[i]][0]:
            j += 1
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j

    w_plus = sum(r for r, (_, s) in zip(ranks, nz) if s > 0)
    w_minus = sum(r for r, (_, s) in zip(ranks, nz) if s < 0)
    mean_diff = sum(diffs) / float(len(diffs))
    n_lt = sum(1 for d in diffs if d < 0)

    if n <= 20:
        all_w = [
            sum(ranks[bit] for bit in range(n) if mask & (1 << bit))
            for mask in range(1 << n)
        ]
        n_null = float(len(all_w))
        p_less = sum(1 for w in all_w if w <= w_plus + 1e-12) / n_null
        p_greater = sum(1 for w in all_w if w >= w_plus - 1e-12) / n_null
        p_two = min(1.0, 2.0 * min(p_less, p_greater))
    else:
        mu = n * (n + 1) / 4.0
        var = n * (n + 1) * (2 * n + 1) / 24.0
        # continuity correction
        z_less = (w_plus + 0.5 - mu) / math.sqrt(var)
        z_greater = (w_plus - 0.5 - mu) / math.sqrt(var)
        # Φ(z) via erf
        def _phi(z):
            return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

        p_less = _phi(z_less)
        p_greater = 1.0 - _phi(z_greater)
        p_two = min(1.0, 2.0 * min(p_less, p_greater))

    return {
        "n": n,
        "mean_diff": mean_diff,
        "n_lt": n_lt,
        "n_pairs": len(diffs),
        "W_plus": w_plus,
        "W_minus": w_minus,
        "p_two": p_two,
        "p_less": p_less,
        "p_greater": p_greater,
    }


def series_by_seed(runs, metrics, key):
    """Map seed_idx -> metric[key] for runs/metrics zipped in parallel."""
    out = {}
    for r, m in zip(runs, metrics):
        v = m.get(key)
        if v is not None:
            out[r["seed_idx"]] = v
    return out


def paired_series(series_a, series_b):
    """Return aligned (a_vals, b_vals) for shared seed_idx keys."""
    keys = sorted(set(series_a.keys()) & set(series_b.keys()))
    return [series_a[k] for k in keys], [series_b[k] for k in keys]


def fmt(x, nd=2, na="—"):
    if x is None:
        return na
    return "{:.{nd}f}".format(x, nd=nd)


def fmt_ms(mean, std, nd=2):
    if mean is None:
        return "—"
    if std is None or std == 0.0:
        return fmt(mean, nd)
    return "{:.{nd}f}±{:.{nd}f}".format(mean, std, nd=nd)


def relocate_outdir(outdir, run_root):
    """Map pre-reorg result dirs (…/p10/Scheme/1/results) onto the run root."""
    if os.path.isdir(outdir):
        return outdir
    parts = [p for p in os.path.normpath(outdir).split(os.sep) if p]
    for i in range(len(parts)):
        cand = os.path.join(run_root, *parts[i:])
        if os.path.isdir(cand):
            return cand
    return outdir


def load_manifest(path):
    run_root = os.path.dirname(os.path.abspath(path))
    run_name = os.path.basename(run_root)
    runs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            outdir = relocate_outdir(parts[2], run_root)
            context = parts[1]
            if not os.path.isdir(parts[2]) and os.path.isdir(outdir):
                slash = context.find("/")
                if slash != -1:
                    context = run_name + context[slash:]
            run = {
                "scheme": parts[0],
                "context": context,
                "outdir": outdir,
                "seed": parts[3] if len(parts) > 3 else "",
                "seed_idx": int(parts[4]) if len(parts) > 4 and parts[4] else 0,
            }
            runs.append(run)
    runs.sort(
        key=lambda r: (
            SCHEME_ORDER.get(r["scheme"], 99),
            r["seed_idx"],
            r["context"],
        )
    )
    return runs


def discover_classify_datasets(outdir):
    """Dataset ids with timed classify files: 'speed' (legacy names) or labeled/argos."""
    found = []
    try:
        names = os.listdir(outdir)
    except OSError:
        return ["speed"]
    for name in names:
        m = re.match(r"classification_performance_([A-Za-z][\w-]*)_\d+\.txt$", name)
        if m:
            found.append(m.group(1))
        elif re.match(r"classification_performance_\d+\.txt$", name):
            found.append("speed")
    if not found:
        legacy = os.path.join(outdir, "classification_performance.txt")
        if os.path.isfile(legacy):
            return ["speed"]
        return ["speed"]
    order = {
        "strain-exclusion": 0,
        "simulated": 0,
        "labeled": 0,
        "fda-argos": 1,
        "argos": 1,
        "speed": 2,
    }
    uniq = []
    for d in found:
        if d not in uniq:
            uniq.append(d)
    return sorted(uniq, key=lambda d: (order.get(d, 9), d))


def list_classify_time_files(outdir, dataset="speed"):
    """Return sorted timed classify files (repeats). Excludes *_cold.txt."""
    files = []
    if dataset == "speed":
        pat = r"classification_performance_\d+\.txt$"
    else:
        pat = r"classification_performance_" + re.escape(dataset) + r"_\d+\.txt$"
    try:
        for name in sorted(os.listdir(outdir)):
            if re.match(pat, name):
                files.append(os.path.join(outdir, name))
    except OSError:
        pass
    if files:
        return files
    if dataset != "speed":
        return []
    legacy = os.path.join(outdir, "classification_performance.txt")
    if os.path.isfile(legacy) and not os.path.islink(legacy):
        # Ignore cold-only legacy name
        if os.path.basename(os.path.realpath(legacy)).endswith("_cold.txt"):
            return []
        return [legacy]
    # Symlink to last repeat still counts as one file if no numbered files found
    if os.path.isfile(legacy):
        target = os.path.realpath(legacy)
        if os.path.isfile(target) and not target.endswith("_cold.txt"):
            return [target]
    return []


def list_nextmin_timing_files(outdir, dataset="speed"):
    """Return sorted timed NextMinimizer files. Excludes *_cold.txt."""
    files = []
    if dataset == "speed":
        pat = r"nextminimizer_timing_\d+\.txt$"
    else:
        pat = r"nextminimizer_timing_" + re.escape(dataset) + r"_\d+\.txt$"
    try:
        for name in sorted(os.listdir(outdir)):
            if re.match(pat, name):
                files.append(os.path.join(outdir, name))
    except OSError:
        pass
    if files:
        return files
    if dataset != "speed":
        return list_classify_time_files(outdir, dataset)
    legacy = os.path.join(outdir, "nextminimizer_timing.txt")
    if os.path.isfile(legacy):
        target = os.path.realpath(legacy)
        if os.path.isfile(target) and not target.endswith("_cold.txt"):
            return [target]
    # Fall back to parsing NEXTMIN_TIMING lines from timed classification_performance_*.txt
    return list_classify_time_files(outdir, dataset)


def classification_paths(outdir, dataset="speed"):
    if dataset == "speed":
        return (
            os.path.join(outdir, "benchmark_classification.txt"),
            os.path.join(outdir, "benchmark_report.txt"),
        )
    return (
        os.path.join(outdir, "{}_classification.txt".format(dataset)),
        os.path.join(outdir, "{}_report.txt".format(dataset)),
    )


def resolve_build_performance_path(outdir):
    """build_performance.txt in outdir (legacy: ../fixed/results/ fallback)."""
    direct = os.path.join(outdir, "build_performance.txt")
    if os.path.isfile(direct):
        return direct
    scheme_dir = os.path.dirname(os.path.dirname(os.path.abspath(outdir)))
    shared = os.path.join(scheme_dir, "fixed", "results", "build_performance.txt")
    if os.path.isfile(shared):
        return shared
    return direct


def collect_run_metrics(run, min_abund, ref_by_seed, dataset="speed"):
    out = run["outdir"]
    bw, br = parse_time_file(resolve_build_performance_path(out))

    class_files = list_classify_time_files(out, dataset)
    class_walls = []
    class_rss = []
    speeds = []
    cls_pcts = []
    n_reads = n_class = None
    for cf in class_files:
        cw, cr = parse_time_file(cf)
        nr, nc = parse_classify_counts(cf)
        if nr is not None:
            n_reads = nr
        if nc is not None:
            n_class = nc
        if cw is not None:
            class_walls.append(cw)
            if cw > 0 and nr:
                speeds.append((nr / 1e6) / (cw / 60.0))
        if cr is not None:
            class_rss.append(cr)
        if nr and nc is not None:
            cls_pcts.append(100.0 * nc / nr)

    class_s_mean, class_s_std = mean_std(class_walls)
    class_gb_mean, class_gb_std = mean_std(class_rss)
    speed_mean, speed_std = mean_std(speeds)
    cls_pct_mean, cls_pct_std = mean_std(cls_pcts)
    n_classify = len(class_walls)

    nextmin_files = list_nextmin_timing_files(out, dataset)
    nextmin_s_vals = []
    nextmin_ns_vals = []
    nextmin_calls_vals = []
    nextmin_cpm_vals = []
    nextmin_frac_vals = []
    scan_s_vals = []
    lookup_s_vals = []
    resolve_s_vals = []
    lookups_vals = []
    repeats_vals = []
    lookups_pm_vals = []
    repeats_pm_vals = []
    lookup_frac_vals = []
    uhs_vals = []
    non_uhs_vals = []
    uhs_frac_vals = []
    non_uhs_frac_vals = []
    for i, nf in enumerate(nextmin_files):
        nm = parse_nextmin_timing_file(nf)
        if not nm:
            continue
        ns = nm.get("nextmin_s")
        npc = nm.get("ns_per_call")
        calls = nm.get("calls")
        cpm = nm.get("calls_per_Mbp")
        if ns is not None:
            nextmin_s_vals.append(ns)
            if i < len(class_walls) and class_walls[i] and class_walls[i] > 0:
                nextmin_frac_vals.append(100.0 * ns / class_walls[i])
        if npc is not None:
            nextmin_ns_vals.append(npc)
        if calls is not None:
            nextmin_calls_vals.append(calls)
        if cpm is not None:
            nextmin_cpm_vals.append(cpm)
        for key, bucket in (
            ("scan_s", scan_s_vals),
            ("lookup_s", lookup_s_vals),
            ("resolve_s", resolve_s_vals),
            ("lookups", lookups_vals),
            ("repeats", repeats_vals),
            ("lookups_per_Mbp", lookups_pm_vals),
            ("repeats_per_Mbp", repeats_pm_vals),
            ("lookup_frac_pct", lookup_frac_vals),
            ("uhs", uhs_vals),
            ("non_uhs", non_uhs_vals),
            ("uhs_frac_pct", uhs_frac_vals),
            ("non_uhs_frac_pct", non_uhs_frac_vals),
        ):
            if key in nm:
                bucket.append(nm[key])

    nextmin_s_mean, nextmin_s_std = mean_std(nextmin_s_vals)
    nextmin_ns_mean, nextmin_ns_std = mean_std(nextmin_ns_vals)
    nextmin_calls_mean, nextmin_calls_std = mean_std(nextmin_calls_vals)
    nextmin_cpm_mean, nextmin_cpm_std = mean_std(nextmin_cpm_vals)
    nextmin_frac_mean, nextmin_frac_std = mean_std(nextmin_frac_vals)
    scan_s_mean, scan_s_std = mean_std(scan_s_vals)
    lookup_s_mean, lookup_s_std = mean_std(lookup_s_vals)
    resolve_s_mean, resolve_s_std = mean_std(resolve_s_vals)
    lookups_mean, lookups_std = mean_std(lookups_vals)
    repeats_mean, repeats_std = mean_std(repeats_vals)
    lookups_pm_mean, lookups_pm_std = mean_std(lookups_pm_vals)
    repeats_pm_mean, repeats_pm_std = mean_std(repeats_pm_vals)
    lookup_frac_mean, lookup_frac_std = mean_std(lookup_frac_vals)
    uhs_mean, uhs_std = mean_std(uhs_vals)
    non_uhs_mean, non_uhs_std = mean_std(non_uhs_vals)
    uhs_frac_mean, uhs_frac_std = mean_std(uhs_frac_vals)
    non_uhs_frac_mean, non_uhs_frac_std = mean_std(non_uhs_frac_vals)

    mape_g = mape_s = l1_g = conc = None
    class_conc = jaccard = newly_cls = newly_uncls = None
    class_path, report_path = classification_paths(out, dataset)
    ref = ref_by_seed.get(run["seed_idx"]) or ref_by_seed.get(0)
    if ref and run["context"] != ref["context"]:
        _, ref_report = classification_paths(ref["outdir"], dataset)
        ref_class, _ = classification_paths(ref["outdir"], dataset)
        ref_total, ref_rows = parse_kraken_report(ref_report)
        ref_g = abundance_at_rank(ref_rows, ref_total, "G")
        ref_s = abundance_at_rank(ref_rows, ref_total, "S")
        total, rows = parse_kraken_report(report_path)
        g = abundance_at_rank(rows, total, "G")
        s = abundance_at_rank(rows, total, "S")
        mape_g, _ = mape(ref_g, g, min_abund)
        mape_s, _ = mape(ref_s, s, min_abund)
        l1_g = l1_percent(ref_g, g)
        if os.path.isfile(ref_class) and os.path.isfile(class_path):
            cv = concordance_vs_ref(ref_class, class_path)
            conc = cv["concordance"]
            class_conc = cv["classified_concordance"]
            jaccard = cv["jaccard"]
            newly_cls = cv["newly_classified_pct"]
            newly_uncls = cv["newly_unclassified_pct"]

    return {
        "dataset": dataset,
        "build_s": bw,
        "build_gb": br,
        "class_s": class_s_mean,
        "class_s_std": class_s_std,
        "class_gb": class_gb_mean,
        "class_gb_std": class_gb_std,
        "mrd_min": speed_mean,
        "mrd_min_std": speed_std,
        "cls_pct": cls_pct_mean,
        "cls_pct_std": cls_pct_std,
        "n_classify": n_classify,
        "class_walls": class_walls,
        "nextmin_s": nextmin_s_mean,
        "nextmin_s_std": nextmin_s_std,
        "nextmin_ns": nextmin_ns_mean,
        "nextmin_ns_std": nextmin_ns_std,
        "nextmin_calls": nextmin_calls_mean,
        "nextmin_calls_std": nextmin_calls_std,
        "nextmin_cpm": nextmin_cpm_mean,
        "nextmin_cpm_std": nextmin_cpm_std,
        "nextmin_frac": nextmin_frac_mean,
        "nextmin_frac_std": nextmin_frac_std,
        "uhs": uhs_mean,
        "uhs_std": uhs_std,
        "non_uhs": non_uhs_mean,
        "non_uhs_std": non_uhs_std,
        "uhs_frac": uhs_frac_mean,
        "uhs_frac_std": uhs_frac_std,
        "non_uhs_frac": non_uhs_frac_mean,
        "non_uhs_frac_std": non_uhs_frac_std,
        "scan_s": scan_s_mean,
        "scan_s_std": scan_s_std,
        "lookup_s": lookup_s_mean,
        "lookup_s_std": lookup_s_std,
        "resolve_s": resolve_s_mean,
        "resolve_s_std": resolve_s_std,
        "lookups": lookups_mean,
        "lookups_std": lookups_std,
        "repeats": repeats_mean,
        "repeats_std": repeats_std,
        "lookups_pm": lookups_pm_mean,
        "lookups_pm_std": lookups_pm_std,
        "repeats_pm": repeats_pm_mean,
        "repeats_pm_std": repeats_pm_std,
        "lookup_frac": lookup_frac_mean,
        "lookup_frac_std": lookup_frac_std,
        "mape_g": mape_g,
        "mape_s": mape_s,
        "l1_g": l1_g,
        "concordance": conc,
        "classified_concordance": class_conc,
        "jaccard": jaccard,
        "newly_classified_pct": newly_cls,
        "newly_unclassified_pct": newly_uncls,
        "n_reads": n_reads,
    }


SKIP_EXPORT_KEYS = frozenset({"class_walls"})

PER_SEED_FIELDS = [
    "scheme",
    "dataset",
    "seed_idx",
    "context",
    "seed",
    "build_s",
    "build_gb",
    "class_s",
    "class_s_std",
    "class_gb",
    "class_gb_std",
    "mrd_min",
    "mrd_min_std",
    "cls_pct",
    "cls_pct_std",
    "n_classify",
    "n_reads",
    "nextmin_s",
    "nextmin_s_std",
    "nextmin_ns",
    "nextmin_cpm",
    "nextmin_frac",
    "scan_s",
    "lookup_s",
    "resolve_s",
    "lookups_pm",
    "repeats_pm",
    "lookup_frac",
    "uhs_frac",
    "non_uhs_frac",
    "mape_g",
    "mape_s",
    "l1_g",
    "concordance",
    "classified_concordance",
    "jaccard",
    "newly_classified_pct",
    "newly_unclassified_pct",
]


def _scheme_base(name):
    return name.split(" [")[0]


def parse_density_summary(path):
    """Parse density summary or CSV: factor (expected/particular) and ns/minimizer.

    CSV may have one row per scheme×seed; values are averaged per scheme.
    """
    out = {}
    if not path or not os.path.isfile(path):
        return out
    if path.endswith(".csv"):
        by = {}
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = _scheme_base((row.get("scheme") or "").strip())
                if not name:
                    continue
                try:
                    rec = {
                        "fac_expected": float(row["fac_expected"]),
                        "fac_particular": float(row["fac_particular"]),
                        "ns_per_minimizer": float(row.get("ns_per_minimizer") or 0),
                    }
                except (TypeError, ValueError, KeyError):
                    continue
                by.setdefault(name, []).append(rec)
        for name, recs in by.items():
            n = float(len(recs))
            out[name] = {
                "fac_expected": sum(r["fac_expected"] for r in recs) / n,
                "fac_particular": sum(r["fac_particular"] for r in recs) / n,
                "ns_per_minimizer": sum(r["ns_per_minimizer"] for r in recs) / n,
            }
        return out
    in_table = False
    with open(path) as f:
        for line in f:
            if line.startswith("Scheme") and "dens_expected" in line:
                in_table = True
                continue
            if not in_table:
                continue
            s = line.strip()
            if not s or s.startswith("---") or s.startswith("Notes"):
                if s.startswith("Notes"):
                    break
                continue
            parts = s.split()
            if len(parts) < 8:
                continue
            try:
                rec = {
                    "fac_expected": float(parts[5]),
                    "fac_particular": float(parts[6]),
                    "ns_per_minimizer": float(parts[7]),
                }
            except (TypeError, ValueError):
                continue
            out[_scheme_base(parts[0])] = rec
    return out


def load_density_by_scheme(run_dir):
    for name in (
        os.path.join("tables", "density_per_seed.csv"),
        "density_per_seed.csv",
        "minimizer_density.csv",
        "minimizer_density_benchmark.txt",
        "minimizer_density_summary.txt",
    ):
        parsed = parse_density_summary(os.path.join(run_dir, name))
        if parsed:
            return parsed
    return {}


def pearson_r(xs, ys):
    n = min(len(xs), len(ys))
    if n < 2:
        return None
    xs = [float(x) for x in xs[:n]]
    ys = [float(y) for y in ys[:n]]
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0.0 or dy == 0.0:
        return None
    return num / (dx * dy)


def _a1_range(col, r0, r1):
    letter = _col_letter(col)
    return "%s%d:%s%d" % (letter, r0, letter, r1)


def timing_correl_cells(fields, rows):
    """Excel CORREL of scheme-level density factor vs build/classify time and RSS.

    Uses one row per scheme (first dataset for build; each dataset for classify).
    """
    if not rows or "fac_expected" not in fields:
        return {}
    if all(r.get("fac_expected") in ("", None) for r in rows):
        return {}
    col = {name: i for i, name in enumerate(fields, 1)}
    need = (
        "fac_expected",
        "fac_particular",
        "build_s_mean",
        "build_gb_mean",
        "class_s_mean",
        "class_gb_mean",
    )
    if any(k not in col for k in need):
        return {}
    by_ds = []
    seen = set()
    for i, row in enumerate(rows, 2):
        ds = row.get("dataset") or ""
        if ds not in seen:
            seen.add(ds)
            by_ds.append((ds, i, i))
        else:
            by_ds[-1] = (by_ds[-1][0], by_ds[-1][1], i)
    if not by_ds:
        return {}
    start_c = len(fields) + 2
    cells = {
        (start_c, 1): "Pearson r (n=schemes)",
        (start_c + 1, 1): "fac_expected",
        (start_c + 2, 1): "fac_particular",
    }
    first_ds, br0, br1 = by_ds[0]
    pairs = [
        ("build wall (s)", "build_s_mean", br0, br1),
        ("build RSS (GB)", "build_gb_mean", br0, br1),
    ]
    for ds, r0, r1 in by_ds:
        tag = " [%s]" % ds if ds else ""
        pairs.append(("classify wall (s)%s" % tag, "class_s_mean", r0, r1))
        pairs.append(("classify RSS (GB)%s" % tag, "class_gb_mean", r0, r1))

    def vals(field, r0, r1):
        return [
            _chart_num(rows[r - 2].get(field))
            for r in range(r0, r1 + 1)
        ]

    for i, (label, field, r0, r1) in enumerate(pairs, 2):
        cells[(start_c, i)] = label
        for j, dens in enumerate(("fac_expected", "fac_particular"), 1):
            formula = "CORREL(%s,%s)" % (
                _a1_range(col[dens], r0, r1),
                _a1_range(col[field], r0, r1),
            )
            cached = pearson_r(vals(dens, r0, r1), vals(field, r0, r1))
            cells[(start_c + j, i)] = {
                "formula": formula,
                "value": cached if cached is not None else "",
            }
    return cells


def _scheme_dataset(name):
    if " [" in name and name.endswith("]"):
        return name[name.index("[") + 1 : -1]
    return "speed"


def _scheme_sort_key(name):
    return (
        DATASET_ORDER.get(_scheme_dataset(name), 9),
        SCHEME_ORDER.get(_scheme_base(name), 99),
        name,
    )


def _row_sort_key(row):
    try:
        seed_idx = int(row.get("seed_idx") or 0)
    except (TypeError, ValueError):
        seed_idx = 0
    return (
        DATASET_ORDER.get(row.get("dataset") or "", 9),
        SCHEME_ORDER.get(row.get("scheme") or "", 99),
        seed_idx,
        DOMAIN_ORDER.get(row.get("domain") or "", 9),
        RANK_ORDER.get(row.get("rank") or "", 9),
        SCOPE_ORDER.get(row.get("scope") or "", 9),
        row.get("sra_run") or "",
    )


def build_export_sheets(by_scheme, metrics_by_scheme, ref_name, dens_by_scheme=None):
    per_seed = []
    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        runs = by_scheme.get(scheme, [])
        ms = metrics_by_scheme[scheme]
        for run, m in zip(runs, ms):
            row = {
                "scheme": _scheme_base(scheme),
                "dataset": m.get("dataset") or _scheme_dataset(scheme),
                "seed_idx": run.get("seed_idx"),
                "context": run.get("context"),
                "seed": run.get("seed"),
            }
            for k, v in m.items():
                if k not in SKIP_EXPORT_KEYS:
                    row[k] = v
            per_seed.append(row)

    summary_fields = [
        "scheme",
        "dataset",
        "n_seeds",
        "fac_expected",
        "fac_particular",
        "ns_per_minimizer",
        "build_s_mean",
        "build_s_std",
        "build_gb_mean",
        "build_gb_std",
        "class_s_mean",
        "class_s_std",
        "class_gb_mean",
        "class_gb_std",
        "mrd_min_mean",
        "mrd_min_std",
        "cls_pct_mean",
        "cls_pct_std",
        "nextmin_s_mean",
        "nextmin_s_std",
        "nextmin_ns_mean",
        "nextmin_ns_std",
        "lookups_pm_mean",
        "lookups_pm_std",
        "lookup_frac_mean",
        "lookup_frac_std",
    ]
    summary = []
    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        ms = metrics_by_scheme[scheme]
        def ms_of(key):
            return mean_std([m[key] for m in ms])

        b, bs = ms_of("build_s")
        bg, bgs = ms_of("build_gb")
        cs, css = ms_of("class_s")
        cg, cgs = ms_of("class_gb")
        sp, sps = ms_of("mrd_min")
        cp, cps = ms_of("cls_pct")
        nms, nmss = ms_of("nextmin_s")
        nns, nnss = ms_of("nextmin_ns")
        lp, lps = ms_of("lookups_pm")
        lf, lfs = ms_of("lookup_frac")
        dens = (dens_by_scheme or {}).get(_scheme_base(scheme), {})
        summary.append(
            {
                "scheme": _scheme_base(scheme),
                "dataset": (ms[0].get("dataset") if ms else None) or _scheme_dataset(scheme),
                "n_seeds": len(ms),
                "fac_expected": dens.get("fac_expected"),
                "fac_particular": dens.get("fac_particular"),
                "ns_per_minimizer": dens.get("ns_per_minimizer"),
                "build_s_mean": b,
                "build_s_std": bs,
                "build_gb_mean": bg,
                "build_gb_std": bgs,
                "class_s_mean": cs,
                "class_s_std": css,
                "class_gb_mean": cg,
                "class_gb_std": cgs,
                "mrd_min_mean": sp,
                "mrd_min_std": sps,
                "cls_pct_mean": cp,
                "cls_pct_std": cps,
                "nextmin_s_mean": nms,
                "nextmin_s_std": nmss,
                "nextmin_ns_mean": nns,
                "nextmin_ns_std": nnss,
                "lookups_pm_mean": lp,
                "lookups_pm_std": lps,
                "lookup_frac_mean": lf,
                "lookup_frac_std": lfs,
            }
        )

    sheets = [
        ("scheme_summary", summary_fields, summary),
        ("per_seed", PER_SEED_FIELDS, per_seed),
    ]
    return sheets


S1_METRICS = ("sensitivity", "ppv", "f1")

XLSX_TAB = {
    "contents": "0 contents",
    "query_sets": "1 query sets",
    "scheme_summary": "2 timing summary",
    "per_seed": "3 timing per seed",
    "accuracy_strain_exclusion_summary": "4 strain-exclusion summary",
    "accuracy_strain_exclusion": "5 strain-exclusion per seed",
    "accuracy_fda_argos_summary": "6 FDA-ARGOS summary",
    "accuracy_fda_argos": "7 FDA-ARGOS per experiment",
}

SHEET_BLURB = {
    "query_sets": (
        "Query FASTQs: sources (genomes or SRA runs), fragments per source, "
        "paired-end fragment count, and read length"
    ),
    "scheme_summary": (
        "Mean±std classify wall / RSS / NextMinimizer over seeds; "
        "density factor and Rank ns/minimizer are scheme-level, not per dataset"
    ),
    "per_seed": "One row per scheme × dataset × seed (timing)",
    "accuracy_strain_exclusion_summary": (
        "Mason held-out genomes: mean±std sensitivity / PPV / F1 "
        "(scheme × domain × rank)"
    ),
    "accuracy_strain_exclusion": (
        "Mason held-out genomes, per seed (TP/VP/FP/FN + sens/PPV/F1)"
    ),
    "accuracy_fda_argos_summary": (
        "FDA-ARGOS genus concordance mean±std over seeds: concordant, "
        "discordant, ancestor, unclassified"
    ),
    "accuracy_fda_argos": (
        "FDA-ARGOS per SRR, plus unweighted mean and pooled fragments, per seed"
    ),
}

SHEET_ORDER = (
    "query_sets",
    "scheme_summary",
    "per_seed",
    "accuracy_strain_exclusion_summary",
    "accuracy_strain_exclusion",
    "accuracy_fda_argos_summary",
    "accuracy_fda_argos",
)

ARGOS_FOUR = (
    "concordant_frac",
    "discordant_frac",
    "ancestor_frac",
    "unclassified_frac",
)


def summarize_strain_exclusion(labeled_rows):
    """Mean±std over seeds of Mason strain-exclusion sensitivity / PPV / F1."""
    groups = defaultdict(list)
    for row in labeled_rows:
        key = (
            row.get("scheme") or "",
            row.get("domain") or "all",
            row.get("rank") or "",
        )
        groups[key].append(row)
    fields = [
        "scheme",
        "domain",
        "rank",
        "n_seeds",
        "sensitivity_mean",
        "sensitivity_std",
        "ppv_mean",
        "ppv_std",
        "f1_mean",
        "f1_std",
    ]
    out = []
    for key in sorted(
        groups,
        key=lambda k: (SCHEME_ORDER.get(k[0], 99), k[0], k[1], k[2]),
    ):
        ms = groups[key]
        rec = {
            "scheme": key[0],
            "domain": key[1],
            "rank": key[2],
            "n_seeds": len(ms),
        }
        for k in S1_METRICS:
            mu, sd = mean_std(
                [float(r[k]) for r in ms if r.get(k) not in ("", None)]
            )
            rec[k + "_mean"] = mu
            rec[k + "_std"] = sd
        out.append(rec)
    return fields, out


def _argos_genus_headline(row):
    """Unweighted-mean genus row, or legacy pooled genus row."""
    if row.get("rank") not in (None, "", "genus"):
        return False
    scope = row.get("scope") or ""
    if scope:
        return scope == "unweighted_mean"
    return not (row.get("sra_run") or "")


def summarize_argos_concordance(argos_rows):
    """Mean±std over seeds of FDA-ARGOS genus concordance (a–d) per scheme."""
    by_scheme = defaultdict(list)
    n_exp_by_scheme = defaultdict(set)
    for row in argos_rows:
        scheme = row.get("scheme") or ""
        if (row.get("scope") or "") == "experiment" and row.get("sra_run"):
            n_exp_by_scheme[scheme].add(row.get("sra_run"))
        if not _argos_genus_headline(row):
            continue
        by_scheme[scheme].append(row)
    fields = ["scheme", "n_seeds", "n_experiments"] + [
        x
        for k in ARGOS_FOUR
        for x in (k.replace("_frac", "_mean"), k.replace("_frac", "_std"))
    ]
    rows = []
    for scheme in sorted(by_scheme, key=lambda s: (SCHEME_ORDER.get(s, 99), s)):
        ms = by_scheme[scheme]
        rec = {
            "scheme": scheme,
            "n_seeds": len(ms),
            "n_experiments": len(n_exp_by_scheme.get(scheme) or ()),
        }
        for k in ARGOS_FOUR:
            mu, sd = mean_std([float(r[k]) for r in ms if r.get(k) not in ("", None)])
            rec[k.replace("_frac", "_mean")] = mu
            rec[k.replace("_frac", "_std")] = sd
        rows.append(rec)
    return fields, rows


def load_accuracy_sheets(runs):
    """Concat per-seed strain_exclusion_metrics.tsv / fda_argos_metrics.tsv."""
    labeled = []
    argos = []
    labeled_fields = []
    argos_fields = []
    for r in runs:
        out = r["outdir"]
        for fname, bucket, kind in (
            ("strain_exclusion_metrics.tsv", labeled, "sim"),
            ("simulated_metrics.tsv", labeled, "sim"),
            ("labeled_metrics.tsv", labeled, "sim"),
            ("fda_argos_metrics.tsv", argos, "argos"),
            ("argos_metrics.tsv", argos, "argos"),
        ):
            path = os.path.join(out, fname)
            if not os.path.isfile(path):
                continue
            fields, rows = load_table(path)
            extra = ["seed_idx", "seed", "context"]
            for row in rows:
                row["seed_idx"] = r.get("seed_idx")
                row["seed"] = r.get("seed")
                row["context"] = r.get("context")
                bucket.append(row)
            merged = extra + [f for f in fields if f not in extra]
            if kind == "sim":
                labeled_fields = merged
            else:
                argos_fields = merged
    labeled.sort(key=_row_sort_key)
    argos.sort(key=_row_sort_key)
    sheets = []
    if labeled:
        sfields, srows = summarize_strain_exclusion(labeled)
        if srows:
            sheets.append(
                ("accuracy_strain_exclusion_summary", sfields, srows)
            )
        sheets.append(("accuracy_strain_exclusion", labeled_fields, labeled))
    if argos:
        sfields, srows = summarize_argos_concordance(argos)
        if srows:
            sheets.append(("accuracy_fda_argos_summary", sfields, srows))
        sheets.append(("accuracy_fda_argos", argos_fields, argos))
    return sheets


def resolve_data_dir(run_root):
    env = os.environ.get("KRAKEN2_BENCHMARK_DIR")
    if env:
        d = os.path.join(env, "data")
        if os.path.isdir(d):
            return d
    parent = os.path.dirname(os.path.dirname(os.path.abspath(run_root or ".")))
    d = os.path.join(parent, "data")
    if os.path.isdir(d):
        return d
    return None


def _fastq_n_and_len(path):
    n = 0
    lengths = Counter()
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                L = len(line.rstrip("\n\r"))
                n += 1
                lengths[L] += 1
    if not n:
        return 0, None, None
    uniq = sorted(lengths)
    if len(uniq) == 1:
        return n, uniq[0], uniq[0]
    return n, uniq[0], uniq[-1]


def _truth_source_stats(path):
    src = Counter()
    if not path or not os.path.isfile(path):
        return 0, None
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            seqid = line.split("\t", 1)[0]
            if seqid.startswith("ARGOS."):
                key = seqid.split(".")[1]
            else:
                key = seqid.rsplit(".", 1)[0]
            src[key] += 1
    if not src:
        return 0, None
    counts = sorted(src.values())
    n_src = len(src)
    if counts[0] == counts[-1]:
        return n_src, counts[0]
    return n_src, "%d–%d" % (counts[0], counts[-1])


def query_sets_sheet(data_dir):
    """One row per query FASTQ directory under data/."""
    if not data_dir or not os.path.isdir(data_dir):
        return None
    meta = (
        (
            "strain-exclusion",
            "held-out genomes (Mason)",
            os.path.join(data_dir, "strain-exclusion"),
        ),
        (
            "fda-argos",
            "SRA experiments (Table S5)",
            os.path.join(data_dir, "fda-argos"),
        ),
    )
    fields = [
        "dataset",
        "sources",
        "n_sources",
        "fragments_per_source",
        "n_fragments",
        "read_length_bp",
        "paired",
        "total_bases_r1_r2",
        "path",
    ]
    rows = []
    for ds, source_kind, d in meta:
        r1 = os.path.join(d, "R1.fastq")
        r2 = os.path.join(d, "R2.fastq")
        if not os.path.isfile(r1):
            continue
        n, len_min, len_max = _fastq_n_and_len(r1)
        n_src, per_src = _truth_source_stats(os.path.join(d, "truth.tsv"))
        if len_min == len_max:
            read_len = len_min
            total_bp = n * read_len * (2 if os.path.isfile(r2) else 1)
        else:
            read_len = "%d–%d" % (len_min, len_max)
            total_bp = None
        rows.append(
            {
                "dataset": ds,
                "sources": source_kind,
                "n_sources": n_src or None,
                "fragments_per_source": per_src,
                "n_fragments": n,
                "read_length_bp": read_len,
                "paired": "yes" if os.path.isfile(r2) else "no",
                "total_bases_r1_r2": total_bp,
                "path": d,
            }
        )
    if not rows:
        return None
    return ("query_sets", fields, rows)


def _contents_sheet(sheets):
    fields = ["tab", "csv", "section", "description"]
    rows = []
    for name, _fields, _rows in sheets:
        if name == "contents":
            continue
        rows.append(
            {
                "tab": XLSX_TAB.get(name, name)[:31],
                "csv": name + ".csv",
                "section": (
                    "query sets"
                    if name == "query_sets"
                    else "strain-exclusion"
                    if "strain_exclusion" in name
                    else "FDA-ARGOS"
                    if "fda_argos" in name
                    else "timing"
                ),
                "description": SHEET_BLURB.get(name, ""),
            }
        )
    return ("contents", fields, rows)


def organize_sheets(timing_sheets, extra_sheets):
    """Timing first, then strain-exclusion accuracy, then FDA-ARGOS concordance."""
    by_name = {}
    for trip in list(timing_sheets) + list(extra_sheets or []):
        by_name[trip[0]] = trip
    ordered = []
    for name in SHEET_ORDER:
        if name in by_name:
            ordered.append(by_name.pop(name))
    ordered.extend(by_name[k] for k in sorted(by_name))
    return [_contents_sheet(ordered)] + ordered


def _chart_num(v):
    if v in ("", None):
        return 0.0
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def _ordered_schemes(rows):
    seen = []
    have = set()
    for row in rows:
        s = row.get("scheme") or ""
        if s and s not in have:
            have.add(s)
            seen.append(s)
    return sorted(seen, key=lambda s: (SCHEME_ORDER.get(s, 99), s))


def _mean_std_series(rows, schemes, series_name, mean_key, std_key, match):
    means, stds = [], []
    for sch in schemes:
        hit = None
        for row in rows:
            if (row.get("scheme") or "") != sch:
                continue
            ok = True
            for k, v in match.items():
                if (row.get(k) or "") != v:
                    ok = False
                    break
            if ok:
                hit = row
                break
        means.append(_chart_num(hit.get(mean_key) if hit else None))
        stds.append(_chart_num(hit.get(std_key) if hit and std_key else None))
    if not any(means) and not any(stds):
        return None
    return {"name": series_name, "means": means, "stds": stds}


def _charts_by_dataset(rows, schemes, metrics):
    """Same layout as timing: clustered columns, one series per dataset."""
    datasets = []
    have = set()
    for row in rows:
        ds = row.get("dataset") or ""
        if ds and ds not in have:
            have.add(ds)
            datasets.append(ds)
    datasets.sort(key=lambda d: DATASET_ORDER.get(d, 9))
    charts = []
    for title, mean_k, std_k in metrics:
        series = []
        for ds in datasets:
            ser = _mean_std_series(
                rows, schemes, ds, mean_k, std_k, {"dataset": ds}
            )
            if ser:
                series.append(ser)
        if series:
            charts.append(
                {"title": title, "y_title": title, "cats": schemes, "series": series}
            )
    return charts


def _chart_one_series(rows, schemes, title, mean_k, std_k):
    """Library-level metric (same across query datasets): one series, no breakdown."""
    ser = _mean_std_series(rows, schemes, title, mean_k, std_k, {})
    if not ser:
        return None
    return {"title": title, "y_title": title, "cats": schemes, "series": [ser]}


def _class_time_by_scheme(timing_rows, dataset_names):
    """scheme -> (class_s_mean, class_s_std) for the first matching dataset."""
    want = set(dataset_names)
    out = {}
    for row in timing_rows or []:
        ds = row.get("dataset") or ""
        if ds not in want:
            continue
        sch = row.get("scheme") or ""
        if sch in out:
            continue
        out[sch] = (
            _chart_num(row.get("class_s_mean")),
            _chart_num(row.get("class_s_std")),
        )
    return out


def _scatter_metric_vs_classify(
    rows, schemes, title, y_title, mean_k, std_k, class_times, match
):
    """One point per scheme: x = classify wall, y = metric. Legend = scheme."""
    series = []
    for sch in schemes:
        hit = None
        for row in rows:
            if (row.get("scheme") or "") != sch:
                continue
            ok = True
            for k, v in match.items():
                if (row.get(k) or "") != v:
                    ok = False
                    break
            if ok:
                hit = row
                break
        if not hit:
            continue
        xs = class_times.get(sch)
        if not xs:
            continue
        series.append(
            {
                "name": sch,
                "x": [xs[0]],
                "y": [_chart_num(hit.get(mean_k))],
                "stds": [_chart_num(hit.get(std_k))],
            }
        )
    if not series:
        return None
    return {
        "type": "scatter",
        "title": title,
        "x_title": "Classify wall (s)",
        "y_title": y_title,
        "series": series,
    }


def charts_for_sheet(csv_name, rows, timing_rows=None):
    """Column charts with mean±std error bars; x-axis = scheme.

    Accuracy summary sheets also get scatter charts of each metric vs classify time.
    """
    if not rows:
        return []
    schemes = _ordered_schemes(rows)
    if not schemes:
        return []
    charts = []
    if csv_name == "scheme_summary":
        build_wall = _chart_one_series(
            rows, schemes, "Build wall (s)", "build_s_mean", "build_s_std"
        )
        build_mem = _chart_one_series(
            rows, schemes, "Build RSS (GB)", "build_gb_mean", "build_gb_std"
        )
        classify = _charts_by_dataset(
            rows,
            schemes,
            (
                ("Classify wall (s)", "class_s_mean", "class_s_std"),
                ("Classify RSS (GB)", "class_gb_mean", "class_gb_std"),
            ),
        )
        fac_exp = _mean_std_series(
            rows, schemes, "expected (random DNA)", "fac_expected", None, {}
        )
        fac_part = _mean_std_series(
            rows, schemes, "particular (DB library)", "fac_particular", None, {}
        )
        fac_series = [s for s in (fac_exp, fac_part) if s]
        density = []
        if fac_series:
            density.append(
                {
                    "title": "Minimizer density factor",
                    "y_title": "Density factor",
                    "y_min": 1.0,
                    "y_max": 2.1,
                    "cats": schemes,
                    "series": fac_series,
                }
            )
        lookups = _charts_by_dataset(
            rows,
            schemes,
            (("Lookups / Mbp", "lookups_pm_mean", "lookups_pm_std"),),
        )
        nextmin = _charts_by_dataset(
            rows,
            schemes,
            (("NextMinimizer (s)", "nextmin_s_mean", "nextmin_s_std"),),
        )
        rank = _chart_one_series(
            rows, schemes, "Rank (ns)", "ns_per_minimizer", None
        )
        charts = (
            ([build_wall] if build_wall else [])
            + ([build_mem] if build_mem else [])
            + classify
            + density
            + lookups
            + nextmin
            + ([rank] if rank else [])
        )
    elif csv_name == "accuracy_strain_exclusion_summary":
        slices = []
        for row in rows:
            sl = (row.get("domain") or "", row.get("rank") or "")
            if sl not in slices:
                slices.append(sl)
        slices.sort(
            key=lambda sl: (
                DOMAIN_ORDER.get(sl[0], 9),
                RANK_ORDER.get(sl[1], 9),
                sl,
            )
        )
        class_times = _class_time_by_scheme(
            timing_rows, ("strain-exclusion", "simulated", "labeled")
        )
        for title, mean_k, std_k in (
            ("Sensitivity", "sensitivity_mean", "sensitivity_std"),
            ("PPV", "ppv_mean", "ppv_std"),
            ("F1", "f1_mean", "f1_std"),
        ):
            series = []
            for domain, rank in slices:
                label = "{} / {}".format(domain, rank) if domain else rank
                ser = _mean_std_series(
                    rows,
                    schemes,
                    label,
                    mean_k,
                    std_k,
                    {"domain": domain, "rank": rank},
                )
                if ser:
                    series.append(ser)
            if series:
                charts.append(
                    {
                        "title": title,
                        "y_title": title,
                        "cats": schemes,
                        "series": series,
                    }
                )
            spec = _scatter_metric_vs_classify(
                rows,
                schemes,
                "{} vs classify time (all / genus)".format(title),
                title,
                mean_k,
                std_k,
                class_times,
                {"domain": "all", "rank": "genus"},
            )
            if spec:
                charts.append(spec)
    elif csv_name == "accuracy_fda_argos_summary":
        class_times = _class_time_by_scheme(
            timing_rows, ("fda-argos", "argos")
        )
        for title, mean_k, std_k in (
            ("Concordant", "concordant_mean", "concordant_std"),
            ("Discordant", "discordant_mean", "discordant_std"),
            ("Ancestor", "ancestor_mean", "ancestor_std"),
            ("Unclassified", "unclassified_mean", "unclassified_std"),
        ):
            ser = _mean_std_series(rows, schemes, title, mean_k, std_k, {})
            if ser:
                charts.append(
                    {
                        "title": title + " (fraction of fragments)",
                        "y_title": "Fraction",
                        "cats": schemes,
                        "series": [ser],
                    }
                )
            spec = _scatter_metric_vs_classify(
                rows,
                schemes,
                "{} vs classify time".format(title),
                title,
                mean_k,
                std_k,
                class_times,
                {},
            )
            if spec:
                charts.append(spec)
    return charts


def export_tables(
    out_txt, by_scheme, metrics_by_scheme, ref_name, extra_sheets=None, dens_by_scheme=None
):
    out_dir = os.path.dirname(os.path.abspath(out_txt)) or "."
    csv_dir = os.path.join(out_dir, "tables")
    xlsx = os.path.splitext(out_txt)[0] + ".xlsx"
    sheets = organize_sheets(
        build_export_sheets(
            by_scheme, metrics_by_scheme, ref_name, dens_by_scheme=dens_by_scheme
        ),
        extra_sheets,
    )
    xlsx_sheets = []
    charts_by_tab = {}
    extra_cells_by_tab = {}
    timing_rows = []
    timing_fields = []
    for name, fields, rows in sheets:
        if name == "scheme_summary":
            timing_rows = rows
            timing_fields = fields
            break
    correl = timing_correl_cells(timing_fields, timing_rows)
    for name, fields, rows in sheets:
        tab = XLSX_TAB.get(name, name)
        xlsx_sheets.append((tab, fields, rows))
        charts = charts_for_sheet(name, rows, timing_rows=timing_rows)
        if charts:
            charts_by_tab[tab] = charts
        if name == "scheme_summary" and correl:
            extra_cells_by_tab[tab] = correl
    write_tables(
        xlsx,
        csv_dir,
        sheets,
        xlsx_sheets=xlsx_sheets,
        charts_by_name=charts_by_tab,
        extra_cells_by_name=extra_cells_by_tab,
    )
    return xlsx, csv_dir

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-dir", required=True)
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-abund", type=float, default=1e-6)
    args = ap.parse_args()

    base = args.base_dir
    manifest = args.manifest or os.path.join(base, ".benchmark_run_manifest")
    out_path = args.out or os.path.join(base, "benchmark_comparison.txt")

    if not os.path.isfile(manifest):
        print("Error: manifest not found: {}".format(manifest), file=sys.stderr)
        return 1

    runs = load_manifest(manifest)
    if not runs:
        print("Error: empty manifest", file=sys.stderr)
        return 1

    by_scheme = defaultdict(list)
    for r in runs:
        by_scheme[r["scheme"]].append(r)

    orig_runs = runs
    datasets = []
    for r in runs:
        for d in discover_classify_datasets(r["outdir"]):
            if d not in datasets:
                datasets.append(d)
    if not datasets:
        datasets = ["speed"]
    use_ds_keys = not (len(datasets) == 1 and datasets[0] == "speed")

    ref_name = REF_SCHEME if REF_SCHEME in by_scheme else None
    if ref_name is None and "RandomXor" in by_scheme:
        ref_name = "RandomXor"  # older RandomXor-labeled manifests
    ref_runs = by_scheme.get(ref_name, []) if ref_name else []
    ref_by_seed = {}
    for r in ref_runs:
        ref_by_seed[r["seed_idx"] or len(ref_by_seed) + 1] = r
    if not ref_by_seed and runs:
        ref_by_seed[0] = runs[0]

    if use_ds_keys:
        flat = defaultdict(list)
        for ds in datasets:
            for scheme, rs in list(by_scheme.items()):
                tagged = []
                for r in rs:
                    rr = dict(r)
                    rr["context"] = "{} [{}]".format(r["context"], ds)
                    tagged.append(rr)
                flat["{} [{}]".format(scheme, ds)] = tagged
        by_scheme = flat

    # Column widths: never trim run/scheme names (GM DNA contexts are >30 chars).
    ctx_lens = [len(r["context"]) for r in orig_runs]
    if use_ds_keys:
        extra = max(len(d) for d in datasets) + 3
        ctx_lens = [n + extra for n in ctx_lens]
    run_w = max([len("Run")] + ctx_lens)
    scheme_w = max([len("Scheme")] + [len(s) for s in by_scheme.keys()])

    n_seeds = max(len(v) for v in by_scheme.values())
    n_classify_reps = 1
    for r in runs:
        nfile = os.path.join(r["outdir"], "n_classify_repeats.txt")
        if os.path.isfile(nfile):
            try:
                with open(nfile) as f:
                    n_classify_reps = max(n_classify_reps, int(f.read().strip()))
            except (ValueError, OSError):
                pass
        else:
            n_classify_reps = max(
                n_classify_reps,
                max(
                    len(list_classify_time_files(r["outdir"], d))
                    for d in discover_classify_datasets(r["outdir"])
                )
                or 0,
            )

    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    emit("=" * 78)
    emit("Kraken 2 minimizer-scheme benchmark — statistical summary")
    emit("Metrics: Wood et al. Genome Biol 2019 (PMC6883579)")
    emit(
        "Seeds: n={}  |  classify repeats/seed: n={}  |  datasets: {}  |  mean±sample-std".format(
            n_seeds, n_classify_reps, ",".join(datasets)
        )
    )
    emit("=" * 78)

    # Per-seed detail (classify columns are mean±std over repeats)
    emit("")
    emit("1. PER-SEED PERFORMANCE  (Class* = mean±std over {} classify repeats)".format(
        n_classify_reps
    ))
    emit(
        "{:<{rw}} {:>8} {:>7} {:>14} {:>14} {:>14} {:>10}".format(
            "Run", "Build(s)", "BldGB", "Class(s)", "ClsGB", "Mrd/min", "Class%",
            rw=run_w,
        )
    )
    emit("-" * (run_w + 68))

    metrics_by_scheme = defaultdict(list)
    for scheme in sorted(by_scheme.keys(), key=_scheme_sort_key):
        ds = _scheme_dataset(scheme) if use_ds_keys else "speed"
        for r in by_scheme[scheme]:
            m = collect_run_metrics(r, args.min_abund, ref_by_seed, dataset=ds)
            metrics_by_scheme[scheme].append(m)
            emit(
                "{:<{rw}} {:>8} {:>7} {:>14} {:>14} {:>14} {:>10}".format(
                    r["context"],
                    fmt(m["build_s"]),
                    fmt(m["build_gb"]),
                    fmt_ms(m["class_s"], m["class_s_std"]),
                    fmt_ms(m["class_gb"], m["class_gb_std"]),
                    fmt_ms(m["mrd_min"], m["mrd_min_std"]),
                    fmt_ms(m["cls_pct"], m["cls_pct_std"]),
                    rw=run_w,
                )
            )

    # Classify repeat detail (optional compact percentiles if many repeats)
    emit("")
    emit("1b. CLASSIFY REPEAT DISTRIBUTION (per seed, wall seconds)")
    emit(
        "{:<{rw}} {:>5} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
            "Run", "n", "mean", "std", "min", "p50", "p90", "max",
            rw=run_w,
        )
    )
    emit("-" * (run_w + 55))
    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        for r, m in zip(by_scheme[scheme], metrics_by_scheme[scheme]):
            walls = sorted(m.get("class_walls") or [])
            if not walls:
                emit(
                    "{:<{rw}} {:>5} {:>8}".format(
                        r["context"], 0, "—", rw=run_w
                    )
                )
                continue
            p50 = walls[(len(walls) - 1) // 2]
            p90 = walls[int(round(0.9 * (len(walls) - 1)))]
            emit(
                "{:<{rw}} {:>5} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
                    r["context"],
                    len(walls),
                    fmt(m["class_s"]),
                    fmt(m["class_s_std"]),
                    fmt(walls[0]),
                    fmt(p50),
                    fmt(p90),
                    fmt(walls[-1]),
                    rw=run_w,
                )
            )

    # NextMinimizer timing (rank calculation), when instrumentation files present
    has_nextmin = any(m.get("nextmin_s") is not None for ms in metrics_by_scheme.values() for m in ms)
    has_phase = any(m.get("scan_s") is not None for ms in metrics_by_scheme.values() for m in ms)
    has_uhs = any(m.get("uhs") is not None for ms in metrics_by_scheme.values() for m in ms)
    if has_nextmin:
        emit("")
        emit(
            "1c. NextMinimizer TIMING  (rank calc; mean±std over classify repeats; "
            "uhs = toggle/tb path used, non_uhs = ∞ early-out)"
        )
        if has_uhs:
            emit(
                "{:<{rw}} {:>12} {:>12} {:>14} {:>12} {:>12} {:>12} {:>12} {:>8}".format(
                    "Run", "NextMin(s)", "ns/call", "calls", "calls/Mbp",
                    "% of wall", "UHS", "non-UHS", "UHS%",
                    rw=run_w,
                )
            )
            emit("-" * (run_w + 98))
            for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
                for r, m in zip(by_scheme[scheme], metrics_by_scheme[scheme]):
                    emit(
                        "{:<{rw}} {:>12} {:>12} {:>14} {:>12} {:>12} {:>12} {:>12} {:>8}".format(
                            r["context"],
                            fmt_ms(m["nextmin_s"], m["nextmin_s_std"], nd=3),
                            fmt_ms(m["nextmin_ns"], m["nextmin_ns_std"], nd=2),
                            fmt_ms(m["nextmin_calls"], m["nextmin_calls_std"], nd=0),
                            fmt_ms(m["nextmin_cpm"], m["nextmin_cpm_std"], nd=1),
                            fmt_ms(m["nextmin_frac"], m["nextmin_frac_std"], nd=1),
                            fmt_ms(m["uhs"], m["uhs_std"], nd=0),
                            fmt_ms(m["non_uhs"], m["non_uhs_std"], nd=0),
                            fmt_ms(m["uhs_frac"], m["uhs_frac_std"], nd=2),
                            rw=run_w,
                        )
                    )
        else:
            emit(
                "{:<{rw}} {:>12} {:>12} {:>14} {:>12} {:>12}".format(
                    "Run", "NextMin(s)", "ns/call", "calls", "calls/Mbp", "% of wall",
                    rw=run_w,
                )
            )
            emit("-" * (run_w + 64))
            for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
                for r, m in zip(by_scheme[scheme], metrics_by_scheme[scheme]):
                    emit(
                        "{:<{rw}} {:>12} {:>12} {:>14} {:>12} {:>12}".format(
                            r["context"],
                            fmt_ms(m["nextmin_s"], m["nextmin_s_std"], nd=3),
                            fmt_ms(m["nextmin_ns"], m["nextmin_ns_std"], nd=2),
                            fmt_ms(m["nextmin_calls"], m["nextmin_calls_std"], nd=0),
                            fmt_ms(m["nextmin_cpm"], m["nextmin_cpm_std"], nd=1),
                            fmt_ms(m["nextmin_frac"], m["nextmin_frac_std"], nd=1),
                            rw=run_w,
                        )
                    )

    if has_phase:
        emit("")
        emit(
            "1d. CLASSIFY PHASES + MINIMIZER DENSITY  "
            "(scan/lookup/resolve are sum of per-thread times)"
        )
        emit(
            "{:<{rw}} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10}".format(
                "Run", "Scan(s)", "Lookup(s)", "Resolve(s)",
                "lookups/Mbp", "repeats/Mbp", "lookup%",
                rw=run_w,
            )
        )
        emit("-" * (run_w + 66))
        for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
            for r, m in zip(by_scheme[scheme], metrics_by_scheme[scheme]):
                emit(
                    "{:<{rw}} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10}".format(
                        r["context"],
                        fmt_ms(m["scan_s"], m["scan_s_std"], nd=3),
                        fmt_ms(m["lookup_s"], m["lookup_s_std"], nd=3),
                        fmt_ms(m["resolve_s"], m["resolve_s_std"], nd=3),
                        fmt_ms(m["lookups_pm"], m["lookups_pm_std"], nd=1),
                        fmt_ms(m["repeats_pm"], m["repeats_pm_std"], nd=1),
                        fmt_ms(m["lookup_frac"], m["lookup_frac_std"], nd=2),
                        rw=run_w,
                    )
                )

    # Statistical summary across seeds (of per-seed classify means)
    emit("")
    emit("2. SCHEME SUMMARY  (mean±std over seeds; Class uses per-seed means of repeats)")
    emit(
        "{:<{sw}} {:>12} {:>12} {:>12} {:>12} {:>12} {:>10}".format(
            "Scheme", "Build(s)", "BldGB", "Class(s)", "ClsGB", "Mrd/min", "Class%",
            sw=scheme_w,
        )
    )
    emit("-" * (scheme_w + 72))

    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        ms = metrics_by_scheme[scheme]
        mark = " *" if _scheme_base(scheme) == _scheme_base(ref_name or "") else ""

        def ms_of(key):
            return mean_std([m[key] for m in ms])

        b_s, b_ss = ms_of("build_s")
        b_g, b_gs = ms_of("build_gb")
        c_s, c_ss = ms_of("class_s")
        c_g, c_gs = ms_of("class_gb")
        sp, sps = ms_of("mrd_min")
        cp, cps = ms_of("cls_pct")
        emit(
            "{:<{sw}} {:>12} {:>12} {:>12} {:>12} {:>12} {:>10}{}".format(
                scheme,
                fmt_ms(b_s, b_ss),
                fmt_ms(b_g, b_gs),
                fmt_ms(c_s, c_ss),
                fmt_ms(c_g, c_gs),
                fmt_ms(sp, sps),
                fmt_ms(cp, cps),
                mark,
                sw=scheme_w,
            )
        )
    emit("  * = reference scheme ({})".format(ref_name or "none"))

    if has_nextmin:
        emit("")
        emit("2b. NextMinimizer SCHEME SUMMARY  (mean±std over seeds)")
        if has_uhs:
            emit(
                "{:<{sw}} {:>12} {:>12} {:>14} {:>12} {:>12} {:>12} {:>12} {:>8}".format(
                    "Scheme", "NextMin(s)", "ns/call", "calls", "calls/Mbp",
                    "% of wall", "UHS", "non-UHS", "UHS%",
                    sw=scheme_w,
                )
            )
            emit("-" * (scheme_w + 98))
        else:
            emit(
                "{:<{sw}} {:>12} {:>12} {:>14} {:>12} {:>12}".format(
                    "Scheme", "NextMin(s)", "ns/call", "calls", "calls/Mbp", "% of wall",
                    sw=scheme_w,
                )
            )
            emit("-" * (scheme_w + 64))
        for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
            ms = metrics_by_scheme[scheme]
            mark = " *" if _scheme_base(scheme) == _scheme_base(ref_name or "") else ""

            def ms_of(key):
                return mean_std([m[key] for m in ms])

            nm_s, nm_ss = ms_of("nextmin_s")
            nm_n, nm_ns = ms_of("nextmin_ns")
            nm_c, nm_cs = ms_of("nextmin_calls")
            nm_p, nm_ps = ms_of("nextmin_cpm")
            nm_f, nm_fs = ms_of("nextmin_frac")
            if has_uhs:
                uh, uhs = ms_of("uhs")
                nu, nus = ms_of("non_uhs")
                uf, ufs = ms_of("uhs_frac")
                emit(
                    "{:<{sw}} {:>12} {:>12} {:>14} {:>12} {:>12} {:>12} {:>12} {:>8}{}".format(
                        scheme,
                        fmt_ms(nm_s, nm_ss, nd=3),
                        fmt_ms(nm_n, nm_ns, nd=2),
                        fmt_ms(nm_c, nm_cs, nd=0),
                        fmt_ms(nm_p, nm_ps, nd=1),
                        fmt_ms(nm_f, nm_fs, nd=1),
                        fmt_ms(uh, uhs, nd=0),
                        fmt_ms(nu, nus, nd=0),
                        fmt_ms(uf, ufs, nd=2),
                        mark,
                        sw=scheme_w,
                    )
                )
            else:
                emit(
                    "{:<{sw}} {:>12} {:>12} {:>14} {:>12} {:>12}{}".format(
                        scheme,
                        fmt_ms(nm_s, nm_ss, nd=3),
                        fmt_ms(nm_n, nm_ns, nd=2),
                        fmt_ms(nm_c, nm_cs, nd=0),
                        fmt_ms(nm_p, nm_ps, nd=1),
                        fmt_ms(nm_f, nm_fs, nd=1),
                        mark,
                        sw=scheme_w,
                    )
                )
        emit(
            "  Note: nextmin_s is wall NextMinimizer() time with classify --threads 1."
        )
        emit(
            "  ns/call = ranking cost per window; calls/Mbp ≈ windows/Mbp (similar across schemes)."
        )
        if has_uhs:
            emit(
                "  UHS = Rank() used toggle/tb pack; non-UHS = ∞ early-out; "
                "UHS% = 100*uhs/(uhs+non_uhs). (GM schemes only.)"
            )

    if has_phase:
        emit("")
        emit("2c. PHASE + DENSITY SCHEME SUMMARY  (mean±std over seeds)")
        emit(
            "{:<{sw}} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10}".format(
                "Scheme", "Scan(s)", "Lookup(s)", "Resolve(s)",
                "lookups/Mbp", "repeats/Mbp", "lookup%",
                sw=scheme_w,
            )
        )
        emit("-" * (scheme_w + 66))
        for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
            ms = metrics_by_scheme[scheme]
            mark = " *" if _scheme_base(scheme) == _scheme_base(ref_name or "") else ""

            def ms_of(key):
                return mean_std([m[key] for m in ms])

            sc, scs = ms_of("scan_s")
            lk, lks = ms_of("lookup_s")
            rs, rss = ms_of("resolve_s")
            lp, lps = ms_of("lookups_pm")
            rp, rps = ms_of("repeats_pm")
            lf, lfs = ms_of("lookup_frac")
            emit(
                "{:<{sw}} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10}{}".format(
                    scheme,
                    fmt_ms(sc, scs, nd=3),
                    fmt_ms(lk, lks, nd=3),
                    fmt_ms(rs, rss, nd=3),
                    fmt_ms(lp, lps, nd=1),
                    fmt_ms(rp, rps, nd=1),
                    fmt_ms(lf, lfs, nd=2),
                    mark,
                    sw=scheme_w,
                )
            )
        emit(
            "  lookups = TOK_LOOKUP (distinct consecutive minimizers queried);"
        )
        emit(
            "  repeats = TOK_REPEAT; lookup% = 100*lookups/(lookups+repeats)."
        )
        emit(
            "  Lower lookups/Mbp ⇒ lower minimizer density at the hash-lookup stage."
        )
        emit(
            "  If GM has higher scan/nextmin but lower lookup_s, density is helping query."
        )

    emit("")
    emit("3. vs {} (paired by seed_idx) — mean±std".format(ref_name or "reference"))
    emit(
        "{:<{sw}} {:>14} {:>14} {:>12} {:>14} {:>12} {:>10} {:>10} {:>12}".format(
            "Scheme",
            "MAPE_genus%",
            "MAPE_species%",
            "L1_genus%",
            "Concordance%",
            "ClassConc%",
            "Jaccard%",
            "NewClass%",
            "NewUnclass%",
            sw=scheme_w,
        )
    )
    emit("-" * (scheme_w + 100))

    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        if _scheme_base(scheme) == _scheme_base(ref_name or REF_SCHEME):
            emit(
                "{:<{sw}} {:>14} {:>14} {:>12} {:>14} {:>12} {:>10} {:>10} {:>12}".format(
                    scheme,
                    "0.00±0.00",
                    "0.00±0.00",
                    "0.00±0.00",
                    "100.00±0.00",
                    "100.00±0.00",
                    "100.00±0.00",
                    "0.00±0.00",
                    "0.00±0.00",
                    sw=scheme_w,
                )
            )
            continue
        ms = metrics_by_scheme[scheme]
        mg, mgs = mean_std([m["mape_g"] for m in ms])
        mp, mps = mean_std([m["mape_s"] for m in ms])
        l1, l1s = mean_std([m["l1_g"] for m in ms])
        cc, ccs = mean_std([m["concordance"] for m in ms])
        clc, clcs = mean_std([m["classified_concordance"] for m in ms])
        jac, jacs = mean_std([m["jaccard"] for m in ms])
        nc, ncs = mean_std([m["newly_classified_pct"] for m in ms])
        nu, nus = mean_std([m["newly_unclassified_pct"] for m in ms])
        emit(
            "{:<{sw}} {:>14} {:>14} {:>12} {:>14} {:>12} {:>10} {:>10} {:>12}".format(
                scheme,
                fmt_ms(mg, mgs),
                fmt_ms(mp, mps),
                fmt_ms(l1, l1s),
                fmt_ms(cc, ccs),
                fmt_ms(clc, clcs),
                fmt_ms(jac, jacs),
                fmt_ms(nc, ncs),
                fmt_ms(nu, nus),
                sw=scheme_w,
            )
        )

    # 4. Paired Wilcoxon signed-rank (by seed_idx)
    schemes_sorted = sorted(
        metrics_by_scheme.keys(), key=_scheme_sort_key
    )
    wilcoxon_metrics = [("Class(s)", "class_s")]
    if has_nextmin:
        wilcoxon_metrics.append(("ns/call", "nextmin_ns"))
    if has_phase:
        wilcoxon_metrics.append(("lookups/Mbp", "lookups_pm"))

    # Build comparisons: each scheme vs Random on the same dataset, then
    # pairwise non-ref schemes (same dataset only).
    ref_base = _scheme_base(ref_name or REF_SCHEME)
    wilcoxon_pairs = []
    for scheme in schemes_sorted:
        if _scheme_base(scheme) == ref_base:
            continue
        ds = _scheme_dataset(scheme)
        ref_k = None
        for k in schemes_sorted:
            if _scheme_base(k) == ref_base and _scheme_dataset(k) == ds:
                ref_k = k
                break
        if ref_k:
            wilcoxon_pairs.append((scheme, ref_k))
    for i, sa in enumerate(schemes_sorted):
        if _scheme_base(sa) == ref_base:
            continue
        for sb in schemes_sorted[i + 1 :]:
            if _scheme_base(sb) == ref_base:
                continue
            if _scheme_dataset(sa) != _scheme_dataset(sb):
                continue
            wilcoxon_pairs.append((sa, sb))

    if wilcoxon_pairs and wilcoxon_metrics and n_seeds >= 2:
        emit("")
        emit(
            "4. PAIRED WILCOXON SIGNED-RANK  (by seed_idx; exact for n≤20 non-zero pairs)"
        )
        emit(
            "    H1 one-sided A<B: A is smaller (faster Class / cheaper Rank / denser)"
        )
        emit(
            "    Pairing: same seed_idx (matched toggle mask) across schemes."
        )
        cmp_w = max(
            len("Comparison"),
            max(len("{} vs {}".format(a, b)) for a, b in wilcoxon_pairs),
        )
        met_w = max(len("Metric"), max(len(m[0]) for m in wilcoxon_metrics))
        emit(
            "{:<{mw}}  {:<{cw}}  {:>3}  {:>8}  {:>10}  {:>9}  {:>8}  {:>8}  {:>8}".format(
                "Metric",
                "Comparison",
                "n",
                "A<B",
                "MeanΔ(A-B)",
                "W+/W-",
                "p_two",
                "p(A<B)",
                "p(A>B)",
                mw=met_w,
                cw=cmp_w,
            )
        )
        emit("-" * (met_w + cmp_w + 68))

        for met_label, met_key in wilcoxon_metrics:
            for sa, sb in wilcoxon_pairs:
                ser_a = series_by_seed(
                    by_scheme[sa], metrics_by_scheme[sa], met_key
                )
                ser_b = series_by_seed(
                    by_scheme[sb], metrics_by_scheme[sb], met_key
                )
                va, vb = paired_series(ser_a, ser_b)
                res = wilcoxon_signed_rank(va, vb)
                cmp_s = "{} vs {}".format(sa, sb)
                if res is None:
                    emit(
                        "{:<{mw}}  {:<{cw}}  {:>3}  {:>8}  {:>10}  {:>9}  {:>8}  {:>8}  {:>8}".format(
                            met_label,
                            cmp_s,
                            len(va),
                            "—",
                            "—",
                            "—",
                            "—",
                            "—",
                            "—",
                            mw=met_w,
                            cw=cmp_w,
                        )
                    )
                    continue
                emit(
                    "{:<{mw}}  {:<{cw}}  {:>3}  {:>4}/{:<3}  {:>+10.3f}  {:>4.0f}/{:<4.0f}  {:>8.4f}  {:>8.4f}  {:>8.4f}".format(
                        met_label,
                        cmp_s,
                        res["n"],
                        res["n_lt"],
                        res["n_pairs"],
                        res["mean_diff"],
                        res["W_plus"],
                        res["W_minus"],
                        res["p_two"],
                        res["p_less"],
                        res["p_greater"],
                        mw=met_w,
                        cw=cmp_w,
                    )
                )
        emit(
            "  Note: n = non-zero paired diffs; A<B = count where A < B (incl. zeros in denom)."
        )
        emit(
            "  Two-sided exact floor at n=5 is p=0.0625 (cannot reach α=0.05 even if unanimous)."
        )
        emit(
            "  One-sided α=0.05 needs n≥5 unanimous (p=0.03125) or larger n with mixed signs."
        )

    emit("")
    emit("Notes:")
    emit(
        "  • Seed = toggle mask (k2 --toggle-mask → build_db -T) for all schemes "
        "(Random order; GreedyMini Γ / Mod tb XOR tie-breaks)."
    )
    emit(
        "  • n={} seeds/scheme; each seed: 1 cold classify (excluded) + {} timed "
        "(Class mean±std within seed).".format(n_seeds, n_classify_reps)
    )
    emit("  • MAPE / concordance use the last classification outputs per seed.")
    emit("  • Lower MAPE / higher concordance ⇒ closer to {} under that seed.".format(
        ref_name or "reference"))
    emit(
        "  • ClassConc% = identical taxid among reads classified by both; "
        "Jaccard% = |ref∩other| / |ref∪other| over sets of classified read IDs; "
        "NewClass%/NewUnclass% = ref U→C / C→U as % of all reads."
    )
    emit(
        "  • Wilcoxon (section 4) pairs schemes by seed_idx on per-seed classify means."
    )
    if has_nextmin:
        emit(
            "  • NextMinimizer / phase timing requires classify built with "
            "ENABLE_MINIMIZER_PROF=1 (-DKRAKEN2_MINIMIZER_PROF); files: "
            "nextminimizer_timing_XX.txt via env KRAKEN2_MINIMIZER_PROF "
            "(cold: *_cold.txt, excluded)."
        )
        emit(
            "  • Phases: scan (minimizer+tokens), lookup (GetBatch), "
            "resolve (token replay+ResolveTree)."
        )
        emit(
            "  • UHS / non-UHS: GM Rank() toggle/tb pack vs ∞ early-out; "
            "UHS% = 100*uhs/(uhs+non_uhs)."
        )
    emit("")
    acc_sheets = load_accuracy_sheets(orig_runs)
    qs = query_sets_sheet(resolve_data_dir(base))
    extra_sheets = ([qs] if qs else []) + list(acc_sheets or [])
    if acc_sheets:
        emit("")
        emit("5. ACCURACY vs LABELS  (Wood et al. 2019; last timed classify per seed)")
        emit(
            "  Strain-exclusion: genus/species sensitivity, PPV, F1 "
            "(simulated truth)."
        )
        emit(
            "  FDA-ARGOS: genus concordance only — (a) concordant (b) discordant "
            "(c) ancestor (d) unclassified; no sensitivity/PPV."
        )
        for name, fields, rows in acc_sheets:
            emit("")
            emit("  {} (n={})".format(name, len(rows)))
            if name == "accuracy_strain_exclusion":
                emit("    (per-seed detail in tables/accuracy_strain_exclusion.csv)")
                continue
            if name == "accuracy_fda_argos":
                emit("    (per-SRR detail in tables/accuracy_fda_argos.csv)")
                continue
            if name == "accuracy_strain_exclusion_summary":
                show = [c for c in fields]
                display = rows
            elif name == "accuracy_fda_argos_summary":
                show = [c for c in fields]
                display = rows
            else:
                show = [
                    c
                    for c in (
                        "scheme",
                        "dataset",
                        "domain",
                        "rank",
                        "seed_idx",
                        "sensitivity",
                        "ppv",
                        "f1",
                    )
                    if c in fields
                ]
                display = rows
            if not show:
                continue

            def _cell(v):
                if isinstance(v, float):
                    v = "{:.6g}".format(v)
                return "{:<14}".format("" if v is None else v)[:14]

            emit("    " + "  ".join(_cell(c) for c in show))
            for row in display:
                emit("    " + "  ".join(_cell(row.get(c, "")) for c in show))

    emit("Output directories:")
    for r in runs:
        seed = r["seed"] or "?"
        emit("  {} (seed {}) -> {}".format(r["context"], seed, r["outdir"]))

    xlsx_path = csv_dir = None
    try:
        dens_by_scheme = load_density_by_scheme(base)
        xlsx_path, csv_dir = export_tables(
            out_path,
            by_scheme,
            metrics_by_scheme,
            ref_name,
            extra_sheets=extra_sheets,
            dens_by_scheme=dens_by_scheme,
        )
    except Exception as e:
        emit("")
        emit("Could not write CSV/XLSX tables: {}".format(e))

    try:
        with open(out_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        emit("")
        emit("Wrote {}".format(out_path))
        if xlsx_path:
            emit("Wrote {}".format(xlsx_path))
            emit("Wrote CSVs under {}".format(csv_dir))
    except OSError as e:
        fallback = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "benchmark_comparison.txt",
        )
        fallback = os.path.abspath(fallback)
        try:
            with open(fallback, "w") as f:
                f.write("\n".join(lines) + "\n")
            emit("")
            emit("Could not write {}: {}".format(out_path, e))
            emit("Wrote fallback {}".format(fallback))
        except OSError as e2:
            emit("")
            emit("Could not write summary file: {} / {}".format(e, e2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
