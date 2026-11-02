#!/usr/bin/env python3
"""Compare Kraken 2 minimizer-scheme benchmarks across toggle-mask seeds.

Manifest columns (tab-separated, written by run_benchmark.sh):
  scheme  context  outdir  seed  seed_idx  gm_k  projection  gm_table  mod_t
`scheme` is the canonical scheme (xor, xor-pcg-hash, greedymini,
mod-minimizer-xor, mod-minimizer-greedymini); the last four are the separate
identity fields of benchmark/common/scheme_names.py (blank when not applicable).

Every (scheme, seed, dataset) is classified exactly once per
--minimum-hit-groups value m in {0, 2}. Per cell the results dir holds
  classification_performance_<dataset>_mhg<m>.txt
      classify stderr (summary lines) followed by the GNU time -v block;
      wall = DB load + processing, processing = the "processed in" seconds
  <dataset>_classification_mhg<m>.txt, <dataset>_report_mhg<m>.txt
  <dataset>_metrics_mhg<m>.tsv (strain_exclusion, fda_argos, cami)
There are no repeats: the replicates are the seeds. Per-seed classify rows
are summarised as mean±std over seeds (n = number of seeds).

Phase instrumentation (classify built with ENABLE_MINIMIZER_PROF=1 ->
-DKRAKEN2_MINIMIZER_PROF; experiments/profile_run.sh) writes, via env
KRAKEN2_MINIMIZER_PROF, classify_profile_<dataset>_mhg<m>.txt
(key=value: scan_s, lookup_s, resolve_s, lookups, repeats, lookups_per_Mbp,
lookup_frac_pct, lookup_misses, lookup_miss_pct, uhs, non_uhs, ...). Absent for non-PROF
builds: the profile output is then skipped.

Paired metrics (MAPE, concordance, ClassConc, Jaccard, NewClass/NewUnclass)
compare each seed to `xor` with the same seed_idx, dataset and
minimum-hit-groups.

Section 4 reports exact paired Wilcoxon signed-rank tests (by seed_idx) on
Class wall, Rank ns/call, and lookups/Mbp -- no scipy required.
Seed = toggle mask for all schemes (xor order; greedymini Γ / mod-minimizer tb XOR).

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
from bench_const import REF_SCHEME, scheme_rank  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import step_record  # noqa: E402
import scheme_names  # noqa: E402
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))
import units  # noqa: E402
from xlsx_util import load_table  # noqa: E402
import export_layout  # noqa: E402


DATASET_ORDER = {
    "strain-exclusion": 0,
    "simulated": 0,
    "labeled": 0,
    "fda-argos": 1,
    "argos": 1,
    "cami": 2,
    "speed": 3,
}
DOMAIN_ORDER = {"all": 0, "prokaryote": 1, "virus": 2, "viral": 2}
RANK_ORDER = {"genus": 0, "species": 1}
SCOPE_ORDER = {"experiment": 0, "unweighted_mean": 1, "pooled": 2}


_LOG_TS = re.compile(r"^\[(?:INFO|DEBUG) - (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d),(\d{3})\]: (.*)$")


def parse_build_phases(path):
    """Phases of the k2 build from the time stamps of its own log (build_performance.txt): estimating the
    hash-table capacity (estimate_capacity), then build_db (computing the selected k-mers of every library
    sequence, inserting them into the hash table, writing the database). Seconds; a phase whose stamps are
    missing is left out."""
    import datetime
    stamps = {}
    try:
        with open(path, errors="replace") as f:
            text = f.read()
    except OSError:
        return {}
    for line in text.replace("[DEBUG", "\n[DEBUG").replace("[INFO", "\n[INFO").splitlines():
        m = _LOG_TS.match(line.strip())
        if not m:
            continue
        t = datetime.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S").timestamp() + int(m.group(2)) / 1000
        msg = m.group(3)
        for key, prefix in (("estimate", "Running: "), ("build", "Starting database build"),
                            ("done", "Finished building database")):
            if msg.startswith(prefix) and (key != "estimate" or "estimate_capacity" in msg):
                stamps.setdefault(key, t)
    out = {}
    if "estimate" in stamps and "build" in stamps:
        out["build_estimate_s"] = round(stamps["build"] - stamps["estimate"], 3)
    if "build" in stamps and "done" in stamps:
        out["build_db_s"] = round(stamps["done"] - stamps["build"], 3)
    return out


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
        rss_gb = units.gb_from_kib(int(m.group(1)))   # GNU time KiB -> decimal GB
    return wall_s, rss_gb


def parse_time_io(path):
    """(major page faults, file-system reads in decimal GB) of a GNU `time -v` log: the disk I/O the
    step needed (~0 = served from the page cache). 'File system inputs' counts 512-byte blocks."""
    if not os.path.isfile(path):
        return None, None
    with open(path) as f:
        text = f.read()
    m = re.search(r"Major \(requiring I/O\) page faults:\s*(\d+)", text)
    major = int(m.group(1)) if m else None
    m = re.search(r"File system inputs:\s*(\d+)", text)
    read_gb = units.gb_from_blocks(int(m.group(1))) if m else None
    return major, read_gb


def read_cpu_list(path):
    """First line of a pinned-CPU list file (cpus_build.txt / *_cpus.txt), or None when absent/empty."""
    if not path or not os.path.isfile(path):
        return None
    with open(path) as f:
        line = f.readline().strip()
    return line or None


def parse_time_cpu(path):
    """User + system CPU seconds (all threads) of a GNU `time -v` log, or None."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        text = f.read()
    u = re.search(r"User time \(seconds\):\s*([0-9.]+)", text)
    sy = re.search(r"System time \(seconds\):\s*([0-9.]+)", text)
    return round(float(u.group(1)) + float(sy.group(1)), 2) if u and sy else None


def parse_processing_s(path):
    """Seconds of the classify "N sequences (X Mbp) processed in S s" line (the
    classification itself, without loading the database), or None."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        for line in f:
            m = re.search(r"sequences\s+\([^)]*\)\s+processed in\s+([0-9.]+)s", line)
            if m:
                return float(m.group(1))
    return None


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


def parse_profile_file(path):
    """Parse classify_profile_XX.txt (key=value) or the PHASE_TIMING stderr line."""
    if not os.path.isfile(path):
        return None
    data = {}
    with open(path) as f:
        text = f.read()
    keys = (
        "scan_s",
        "lookup_s",
        "resolve_s",
        "lookups",
        "repeats",
        "skips",
        "lookups_per_Mbp",
        "repeats_per_Mbp",
        "lookup_frac_pct",
        "lookup_misses",
        "lookup_miss_pct",
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
    # Fallback: the stderr line embedded in classification_performance_*.txt
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
    if "scan_s" not in data:
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
            extra = dict(zip(("gm_k", "projection", "gm_table", "mod_t"),
                             (parts[5:9] + [""] * 4)[:4]))
            idd = scheme_names.parse_id(dict(extra, scheme=parts[0]))
            run = {
                "label": parts[0],
                "scheme": parts[0],
                "context": context,
                "outdir": outdir,
                "seed": parts[3] if len(parts) > 3 else "",
                "seed_idx": int(parts[4]) if len(parts) > 4 and parts[4] else 0,
            }
            run.update((f, idd[f]) for f in scheme_names.ID_FIELDS if f != "scheme")
            runs.append(run)
    # Display label: the scheme, plus " (<params>)" only when a library holds
    # several parameterisations of the same scheme (e.g. both projections).
    variants = {}
    for r in runs:
        r["params"] = scheme_names.scheme_dir_name(
            r["scheme"], r["gm_k"], r["projection"], r["mod_t"]).partition("/")[2]
        variants.setdefault(r["scheme"], set()).add(r["params"])
    for r in runs:
        if len(variants[r["scheme"]]) > 1 and r["params"]:
            r["label"] = "{} ({})".format(r["scheme"], r["params"])
    runs.sort(
        key=lambda r: (
            scheme_rank(r["label"]),
            r["label"],
            r["seed_idx"],
            r["context"],
        )
    )
    return runs


PERF_RE = re.compile(r"^classification_performance_(.+)_mhg(\d+)\.txt$")


def discover_classify_cells(outdir):
    """[(dataset, minimum_hit_groups)] of the classification_performance_<ds>_mhg<m>.txt
    files of one results dir, sorted dataset -> mhg."""
    found = []
    try:
        names = os.listdir(outdir)
    except OSError:
        return found
    for name in names:
        m = PERF_RE.match(name)
        if m and (m.group(1), int(m.group(2))) not in found:
            found.append((m.group(1), int(m.group(2))))
    return sorted(found, key=lambda c: (DATASET_ORDER.get(c[0], 9), c[0], c[1]))


def perf_path(outdir, dataset, mhg):
    return os.path.join(outdir, "classification_performance_{}_mhg{}.txt".format(dataset, mhg))


def profile_path(outdir, dataset, mhg):
    """classify_profile_<ds>_mhg<m>.txt (ENABLE_MINIMIZER_PROF=1 builds only), else None."""
    p = os.path.join(outdir, "classify_profile_{}_mhg{}.txt".format(dataset, mhg))
    return p if os.path.isfile(p) else None


def classification_paths(outdir, dataset, mhg):
    return (
        os.path.join(outdir, "{}_classification_mhg{}.txt".format(dataset, mhg)),
        os.path.join(outdir, "{}_report_mhg{}.txt".format(dataset, mhg)),
    )


def resolve_build_performance_path(outdir):
    return os.path.join(outdir, "build_performance.txt")


def collect_run_metrics(run, min_abund, ref_by_seed, dataset, mhg):
    """Metrics of ONE classify run (seed x dataset x minimum-hit-groups): wall,
    DB-load and processing time, RSS, classified %, optional phase
    profile and the paired agreement with the reference scheme."""
    out = run["outdir"]
    bw, br = parse_time_file(resolve_build_performance_path(out))

    cf = perf_path(out, dataset, mhg)
    class_s, class_gb = parse_time_file(cf)
    b_cpu = parse_time_cpu(resolve_build_performance_path(out))
    c_cpu = parse_time_cpu(cf)
    b_major, b_read_gb = parse_time_io(resolve_build_performance_path(out))
    c_major, c_read_gb = parse_time_io(cf)
    proc_s = parse_processing_s(cf)
    build_pinned = read_cpu_list(os.path.join(out, "cpus_build.txt"))
    class_pinned = read_cpu_list(os.path.splitext(cf)[0] + "_cpus.txt")
    # node / page-size / foreign-load record of each step (common/bench_step.py)
    records = dict(step_record.columns(os.path.join(out, "step_build.json"), "build"))
    records.update(step_record.columns(os.path.splitext(cf)[0] + "_step.json", "classify"))
    n_reads, n_class = parse_classify_counts(cf)
    db_load_s = None
    if class_s is not None and proc_s is not None:
        # /usr/bin/time prints the wall clock to 10 ms; processing is in ms
        db_load_s = max(0.0, class_s - proc_s)
    speed = None
    dur = proc_s if proc_s else class_s
    if dur and n_reads:
        speed = (n_reads / 1e6) / (dur / 60.0)
    cls_pct = 100.0 * n_class / n_reads if (n_reads and n_class is not None) else None

    nm = parse_profile_file(profile_path(out, dataset, mhg) or cf) or {}
    scan_pct = lookup_pct = resolve_pct = ns_per_lookup = None
    phases = [nm.get(k) for k in ("scan_s", "lookup_s", "resolve_s")]
    if all(v is not None for v in phases) and sum(phases) > 0:
        scan_pct, lookup_pct, resolve_pct = (100.0 * v / sum(phases) for v in phases)
    if nm.get("lookup_s") is not None and nm.get("lookups"):
        ns_per_lookup = 1e9 * nm["lookup_s"] / nm["lookups"]

    mape_g = mape_s = l1_g = conc = None
    class_conc = jaccard = newly_cls = newly_uncls = None
    class_path, report_path = classification_paths(out, dataset, mhg)
    ref = ref_by_seed.get(run["seed_idx"]) or ref_by_seed.get(0)
    if ref and ref["outdir"] != out:
        ref_class, ref_report = classification_paths(ref["outdir"], dataset, mhg)
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
        "minimum_hit_groups": mhg,
        "build_s": bw,
        "build_gb": br,
        "build_cpu_s": b_cpu,
        **parse_build_phases(resolve_build_performance_path(out)),
        "class_cpu_s": c_cpu,
        "build_major_pf": b_major,
        "build_read_gb": b_read_gb,
        "class_major_pf": c_major,
        "class_read_gb": c_read_gb,
        "class_s": class_s,
        "db_load_s": db_load_s,
        "build_pinned_cpus": build_pinned,
        "classify_pinned_cpus": class_pinned,
        **records,
        "class_proc_s": proc_s,
        "class_gb": class_gb,
        "mrd_min": speed,
        "cls_pct": cls_pct,
        "uhs": nm.get("uhs"),
        "non_uhs": nm.get("non_uhs"),
        "uhs_frac": nm.get("uhs_frac_pct"),
        "non_uhs_frac": nm.get("non_uhs_frac_pct"),
        "scan_s": nm.get("scan_s"),
        "lookup_s": nm.get("lookup_s"),
        "resolve_s": nm.get("resolve_s"),
        "lookups": nm.get("lookups"),
        "repeats": nm.get("repeats"),
        "lookups_pm": nm.get("lookups_per_Mbp"),
        "repeats_pm": nm.get("repeats_per_Mbp"),
        "lookup_frac": nm.get("lookup_frac_pct"),
        "miss_pct": nm.get("lookup_miss_pct"),
        "scan_pct": scan_pct,
        "lookup_pct": lookup_pct,
        "resolve_pct": resolve_pct,
        "ns_per_lookup": ns_per_lookup,
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


DENSITY_FILES = ("minimizer_density.csv",)     # results/runs/<library>/ (input written by the density microbenchmark)


def load_density_rows(run_dir):
    """Per-(scheme, seed) density rows ({scheme, seed_idx, k, w, L, fac_expected,
    fac_particular, ns_per_minimizer, ...}) of the run root, from the density
    microbenchmark CSV. Schemes missing from the CSV fall back to the
    scheme-level summary text (one row per scheme, no seed_idx)."""
    rows = []
    for name in DENSITY_FILES:
        path = os.path.join(run_dir, name)
        if not os.path.isfile(path):
            continue
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                scheme = _scheme_base((row.get("scheme") or "").strip())
                if not scheme:
                    continue
                rec = {"scheme": scheme}
                if "L" not in row and "l" in row:      # legacy header: l = minimizer length, k = k-mer length
                    row = {({"l": "k", "k": "L"}.get(k0, k0)): v0 for k0, v0 in row.items()}
                for k, v in row.items():
                    if k in ("scheme", "toggle_hex", "variant", "context") or v in (None, ""):
                        continue
                    if k in scheme_names.ID_FIELDS:   # gm_k, projection, gm_table, mod_t
                        rec[k] = int(v) if k in ("gm_k", "mod_t") else v
                        continue
                    try:
                        rec[k] = int(v) if k == "seed_idx" else float(v)
                    except ValueError:
                        pass
                if "fac_expected" in rec or "fac_particular" in rec:
                    rows.append(rec)
        if rows:
            break
    have = set(r["scheme"] for r in rows)
    for name in ("minimizer_density_benchmark.txt", "minimizer_density_summary.txt"):
        parsed = parse_density_summary(os.path.join(run_dir, name))
        if parsed:
            for scheme, rec in parsed.items():
                if scheme not in have:
                    row = {"scheme": scheme}
                    row.update(rec)
                    rows.append(row)
            break
    return rows


ACCURACY_KIND = {
    "accuracy_strain_exclusion": "labeled",
    "accuracy_fda_argos": "argos",
    "accuracy_cami": "cami",
}


def _cell_key(scheme, dataset, mhg):
    """by_scheme key of one (scheme, dataset, minimum-hit-groups) cell."""
    return "{} [{}, mhg{}]".format(scheme, dataset, mhg)


def _scheme_cell(name):
    """(dataset, mhg) of a _cell_key() name; ('speed', 0) for a bare scheme."""
    if " [" in name and name.endswith("]"):
        inner = name[name.index("[") + 1 : -1]
        ds, _, m = inner.rpartition(", mhg")
        try:
            return ds, int(m)
        except ValueError:
            pass
    return "speed", 0


def _scheme_dataset(name):
    return _scheme_cell(name)[0]


def _scheme_mhg(name):
    return _scheme_cell(name)[1]


def _scheme_sort_key(name):
    return (
        DATASET_ORDER.get(_scheme_dataset(name), 9),
        _scheme_mhg(name),
        scheme_rank(_scheme_base(name)),
        name,
    )


def _int0(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _row_sort_key(row):
    seed_idx = _int0(row.get("seed_idx"))
    return (
        DATASET_ORDER.get(row.get("dataset") or "", 9),
        _int0(row.get("minimum_hit_groups")),
        scheme_rank(row.get("scheme") or ""),
        seed_idx,
        DOMAIN_ORDER.get(row.get("domain") or "", 9),
        RANK_ORDER.get(row.get("rank") or "", 9),
        SCOPE_ORDER.get(row.get("scope") or "", 9),
        row.get("sra_run") or "",
        row.get("group") or "",
        str(row.get("min_frac") if row.get("min_frac") is not None else ""),
    )


def collect_per_seed_rows(by_scheme, metrics_by_scheme):
    """One dict per (scheme x dataset x minimum-hit-groups x seed): identity +
    collect_run_metrics(); the export layer pivots these into the
    `runs` sheet (export_layout.build_tables)."""
    per_seed = []
    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        runs = by_scheme.get(scheme, [])
        ms = metrics_by_scheme[scheme]
        for run, m in zip(runs, ms):
            row = {
                "dataset": m.get("dataset") or _scheme_dataset(scheme),
                "minimum_hit_groups": m.get("minimum_hit_groups", _scheme_mhg(scheme)),
                "seed_idx": run.get("seed_idx"),
                "context": run.get("context"),
                "seed": run.get("seed"),
            }
            row.update((f, run.get(f)) for f in scheme_names.ID_FIELDS)
            row.update(m)
            per_seed.append(row)
    return per_seed


S1_METRICS = ("sensitivity", "ppv", "f1")

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
            _int0(row.get("minimum_hit_groups")),
            row.get("domain") or "all",
            row.get("rank") or "",
        )
        groups[key].append(row)
    fields = [
        "scheme",
        "minimum_hit_groups",
        "domain",
        "rank",
        "n_repeats",
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
        key=lambda k: (k[1], scheme_rank(k[0])) + k,
    ):
        ms = groups[key]
        rec = {
            "scheme": key[0],
            "minimum_hit_groups": key[1],
            "domain": key[2],
            "rank": key[3],
            "n_repeats": len(ms),
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
        scheme = (row.get("scheme") or "", _int0(row.get("minimum_hit_groups")))
        if (row.get("scope") or "") == "experiment" and row.get("sra_run"):
            n_exp_by_scheme[scheme].add(row.get("sra_run"))
        if not _argos_genus_headline(row):
            continue
        by_scheme[scheme].append(row)
    fields = ["scheme", "minimum_hit_groups", "n_repeats", "n_experiments"] + [
        x
        for k in ARGOS_FOUR
        for x in (k.replace("_frac", "_mean"), k.replace("_frac", "_std"))
    ]
    rows = []
    for scheme in sorted(by_scheme, key=lambda s: (s[1], scheme_rank(s[0]), s[0])):
        ms = by_scheme[scheme]
        rec = {
            "scheme": scheme[0],
            "minimum_hit_groups": scheme[1],
            "n_repeats": len(ms),
            "n_experiments": len(n_exp_by_scheme.get(scheme) or ()),
        }
        for k in ARGOS_FOUR:
            mu, sd = mean_std([float(r[k]) for r in ms if r.get(k) not in ("", None)])
            rec[k.replace("_frac", "_mean")] = mu
            rec[k.replace("_frac", "_std")] = sd
        rows.append(rec)
    return fields, rows


CAMI_METRICS = ("sensitivity", "ppv", "f1", "excluded_frac", "n_scored")


def _summarize_by(rows, key_fields, metrics):
    """Mean±std over seeds of `metrics`, grouped by key_fields (scheme first)."""
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row.get(k) if row.get(k) is not None else "" for k in key_fields)].append(row)
    fields = list(key_fields) + ["n_repeats"] + [
        x for m in metrics for x in (m + "_mean", m + "_std")
    ]
    out = []
    for key in sorted(
        groups,
        key=lambda k: (_int0(k[1]), scheme_rank(k[0])) + tuple(str(x) for x in k),
    ):
        ms = groups[key]
        rec = dict(zip(key_fields, key))
        rec["n_repeats"] = len(ms)
        for m in metrics:
            mu, sd = mean_std([float(r[m]) for r in ms if r.get(m) not in ("", None)])
            rec[m + "_mean"] = mu
            rec[m + "_std"] = sd
        out.append(rec)
    return fields, out


def summarize_cami(rows):
    """CAMI (library-aware): mean±std sens / PPV / F1 / excluded fraction per scheme × group × rank."""
    return _summarize_by(rows, ("scheme", "minimum_hit_groups", "group", "rank"), CAMI_METRICS)


METRICS_RE = re.compile(r"^(strain_exclusion|fda_argos|cami)_metrics_mhg(\d+)\.tsv$")


def load_accuracy_sheets(runs):
    """Concat the per-seed <dataset>_metrics_mhg<m>.tsv files (both hit-group
    settings; the rows carry a minimum_hit_groups column)."""
    labeled = []
    argos = []
    cami = []
    labeled_fields = []
    argos_fields = []
    cami_fields = []
    kinds = {
        "strain_exclusion": (labeled, "sim"),
        "fda_argos": (argos, "argos"),
        "cami": (cami, "cami"),
    }
    for r in runs:
        out = r["outdir"]
        try:
            names = sorted(os.listdir(out))
        except OSError:
            continue
        for fname in names:
            m = METRICS_RE.match(fname)
            if not m:
                continue
            bucket, kind = kinds[m.group(1)]
            file_mhg = int(m.group(2))
            fields, rows = load_table(os.path.join(out, fname))
            extra = ["seed_idx", "seed", "context"]
            for row in rows:
                if row.get("minimum_hit_groups") in (None, ""):
                    row["minimum_hit_groups"] = file_mhg
                row["seed_idx"] = r.get("seed_idx")
                row["seed"] = r.get("seed")
                row["context"] = r.get("context")
                bucket.append(row)
            if "minimum_hit_groups" not in fields:
                fields = list(fields)
                fields.insert(fields.index("scheme") + 1 if "scheme" in fields else 0, "minimum_hit_groups")
            merged = extra + [f for f in fields if f not in extra]
            if kind == "sim":
                labeled_fields = merged
            elif kind == "cami":
                cami_fields = merged
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
    cami.sort(key=_row_sort_key)
    if cami:
        sfields, srows = summarize_cami(cami)
        if srows:
            sheets.append(("accuracy_cami_summary", sfields, srows))
        sheets.append(("accuracy_cami", cami_fields, cami))
    return sheets


def resolve_data_dir(run_root):
    """<ROOT>/data for <ROOT>/results/runs/<library> (KRAKEN2_DATA_DIR, else
    KRAKEN2_BENCH_ROOT/data, else derived from the run root)."""
    env = os.environ.get("KRAKEN2_DATA_DIR")
    if env and os.path.isdir(env):
        return env
    root = os.environ.get("KRAKEN2_BENCH_ROOT")
    if root and os.path.isdir(os.path.join(root, "data")):
        return os.path.join(root, "data")
    layout = export_layout.resolve_layout(run_root)
    if layout and os.path.isdir(layout["data"]):
        return layout["data"]
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


def _truth_taxon_stats(path):
    """(distinct true taxids, None) for read-level truth without per-source IDs (CAMI)."""
    taxa = set()
    if not path or not os.path.isfile(path):
        return 0, None
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > 1:
                taxa.add(parts[1])
    return len(taxa), None


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
        (
            "cami",
            "CAMI II marine genomes (distinct true taxids)",
            os.path.join(data_dir, "cami"),
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
        if ds == "cami":
            n_src, per_src = _truth_taxon_stats(os.path.join(d, "truth.tsv"))
        else:
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
        by_scheme[r["label"]].append(r)

    orig_runs = runs
    cells = []
    for r in runs:
        for c in discover_classify_cells(r["outdir"]):
            if c not in cells:
                cells.append(c)
    cells.sort(key=lambda c: (DATASET_ORDER.get(c[0], 9), c[0], c[1]))
    if not cells:
        print(
            "Error: no classification_performance_<dataset>_mhg<N>.txt under the results "
            "dirs of {} (this tree was not produced by the current run_benchmark.sh)".format(
                manifest
            ),
            file=sys.stderr,
        )
        return 1
    datasets = []
    for d, _m in cells:
        if d not in datasets:
            datasets.append(d)
    mhgs = sorted(set(m for _d, m in cells))

    ref_name = REF_SCHEME if REF_SCHEME in by_scheme else None
    ref_runs = by_scheme.get(ref_name, []) if ref_name else []
    ref_by_seed = {}
    for r in ref_runs:
        ref_by_seed[r["seed_idx"] or len(ref_by_seed) + 1] = r
    if not ref_by_seed and runs:
        ref_by_seed[0] = runs[0]

    # one cell per scheme x dataset x minimum-hit-groups
    flat = defaultdict(list)
    for ds, mhg in cells:
        for scheme, rs in list(by_scheme.items()):
            tagged = []
            for r in rs:
                if (ds, mhg) not in discover_classify_cells(r["outdir"]):
                    continue
                rr = dict(r)
                rr["context"] = "{} [{}, mhg{}]".format(r["context"], ds, mhg)
                tagged.append(rr)
            if tagged:
                flat[_cell_key(scheme, ds, mhg)] = tagged
    by_scheme = flat

    # Column widths: never trim run/scheme names.
    ctx_lens = [len(r["context"]) + max(len(d) for d in datasets) + 10 for r in orig_runs]
    run_w = max([len("Run")] + ctx_lens)
    scheme_w = max([len("Scheme")] + [len(s) for s in by_scheme.keys()])

    n_seeds = max(len(v) for v in by_scheme.values())

    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    emit("=" * 78)
    emit("Kraken 2 minimizer-scheme benchmark — statistical summary")
    emit("Metrics: Wood et al. Genome Biol 2019 (PMC6883579)")
    emit(
        "Seeds: n={}  |  one classify per seed x dataset x minimum-hit-groups  |  "
        "datasets: {}  |  minimum-hit-groups: {}  |  mean±sample-std over seeds".format(
            n_seeds, ",".join(datasets), ",".join(str(m) for m in mhgs)
        )
    )
    emit("=" * 78)

    # Per-seed detail (one classify run per row)
    emit("")
    emit("1. PER-SEED PERFORMANCE  (one classify run per row; Wall = DBLoad + Proc)")
    emit(
        "{:<{rw}} {:>8} {:>7} {:>9} {:>9} {:>9} {:>7} {:>9} {:>8}".format(
            "Run", "Build(s)", "BldGB", "Wall(s)", "DBLoad(s)", "Proc(s)", "ClsGB",
            "Mrd/min", "Class%",
            rw=run_w,
        )
    )
    emit("-" * (run_w + 76))

    metrics_by_scheme = defaultdict(list)
    for scheme in sorted(by_scheme.keys(), key=_scheme_sort_key):
        ds, mhg = _scheme_cell(scheme)
        for r in by_scheme[scheme]:
            m = collect_run_metrics(r, args.min_abund, ref_by_seed, ds, mhg)
            metrics_by_scheme[scheme].append(m)
            emit(
                "{:<{rw}} {:>8} {:>7} {:>9} {:>9} {:>9} {:>7} {:>9} {:>8}".format(
                    r["context"],
                    fmt(m["build_s"]),
                    fmt(m["build_gb"]),
                    fmt(m["class_s"]),
                    fmt(m["db_load_s"]),
                    fmt(m["class_proc_s"]),
                    fmt(m["class_gb"]),
                    fmt(m["mrd_min"]),
                    fmt(m["cls_pct"]),
                    rw=run_w,
                )
            )

    has_phase = any(m.get("scan_s") is not None for ms in metrics_by_scheme.values() for m in ms)
    # GreedyMini rank path: the ranker loads the table rank first and computes the full value (tie-break
    # fields) only when two finite ranks are equal; infinite ranks (k-mers outside the UHS) skip it
    has_uhs = any(m.get("uhs") is not None for ms in metrics_by_scheme.values() for m in ms)
    if has_uhs:
        emit("")
        emit(
            "1c. GREEDYMINI RANK PATH  (one classify run per row; UHS = finite table rank, "
            "non-UHS = infinite rank, no tie-break value computed)"
        )
        emit("{:<{rw}} {:>14} {:>14} {:>8}".format("Run", "UHS", "non-UHS", "UHS%", rw=run_w))
        emit("-" * (run_w + 40))
        for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
            for r, m in zip(by_scheme[scheme], metrics_by_scheme[scheme]):
                if m.get("uhs") is None:
                    continue
                emit(
                    "{:<{rw}} {:>14} {:>14} {:>8}".format(
                        r["context"], fmt(m["uhs"], 0), fmt(m["non_uhs"], 0), fmt(m["uhs_frac"], 2), rw=run_w,
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
                        fmt(m["scan_s"], 3),
                        fmt(m["lookup_s"], 3),
                        fmt(m["resolve_s"], 3),
                        fmt(m["lookups_pm"], 1),
                        fmt(m["repeats_pm"], 1),
                        fmt(m["lookup_frac"], 2),
                        rw=run_w,
                    )
                )

    # Statistical summary across seeds
    emit("")
    emit("2. SCHEME SUMMARY  (mean±std over seeds, per dataset and minimum-hit-groups)")
    emit(
        "{:<{sw}} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>10}".format(
            "Scheme", "Build(s)", "BldGB", "Wall(s)", "DBLoad(s)", "Proc(s)", "ClsGB",
            "Mrd/min", "Class%",
            sw=scheme_w,
        )
    )
    emit("-" * (scheme_w + 98))

    for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
        ms = metrics_by_scheme[scheme]
        mark = " *" if _scheme_base(scheme) == _scheme_base(ref_name or "") else ""

        def ms_of(key):
            return mean_std([m[key] for m in ms])

        b_s, b_ss = ms_of("build_s")
        b_g, b_gs = ms_of("build_gb")
        c_s, c_ss = ms_of("class_s")
        d_s, d_ss = ms_of("db_load_s")
        p_s, p_ss = ms_of("class_proc_s")
        c_g, c_gs = ms_of("class_gb")
        sp, sps = ms_of("mrd_min")
        cp, cps = ms_of("cls_pct")
        emit(
            "{:<{sw}} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>10}{}".format(
                scheme,
                fmt_ms(b_s, b_ss),
                fmt_ms(b_g, b_gs),
                fmt_ms(c_s, c_ss),
                fmt_ms(d_s, d_ss),
                fmt_ms(p_s, p_ss),
                fmt_ms(c_g, c_gs),
                fmt_ms(sp, sps),
                fmt_ms(cp, cps),
                mark,
                sw=scheme_w,
            )
        )
    emit("  * = reference scheme ({})".format(ref_name or "none"))

    if has_uhs:
        emit("")
        emit("2b. GREEDYMINI RANK PATH SCHEME SUMMARY  (mean±std over seeds)")
        emit("{:<{sw}} {:>14} {:>14} {:>8}".format("Scheme", "UHS", "non-UHS", "UHS%", sw=scheme_w))
        emit("-" * (scheme_w + 40))
        for scheme in sorted(metrics_by_scheme.keys(), key=_scheme_sort_key):
            ms = [m for m in metrics_by_scheme[scheme] if m.get("uhs") is not None]
            if not ms:
                continue
            (uh, uhs), (nu, nus), (uf, ufs) = (mean_std([m[k] for m in ms]) for k in ("uhs", "non_uhs", "uhs_frac"))
            emit(
                "{:<{sw}} {:>14} {:>14} {:>8}".format(
                    scheme, fmt_ms(uh, uhs, nd=0), fmt_ms(nu, nus, nd=0), fmt_ms(uf, ufs, nd=2), sw=scheme_w,
                )
            )
        emit("  UHS% = 100*uhs/(uhs+non_uhs): share of ranks that take the finite (table) path. GreedyMini schemes only.")

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
            "  If GM has a higher scan_s but a lower lookup_s, density is helping query."
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
    wilcoxon_metrics = [("Wall(s)", "class_s"), ("Proc(s)", "class_proc_s")]
    if has_phase:
        wilcoxon_metrics.append(("lookups/Mbp", "lookups_pm"))

    # Build comparisons: each scheme vs the reference on the same dataset and
    # minimum-hit-groups, then pairwise non-ref schemes (same cell only).
    ref_base = _scheme_base(ref_name or REF_SCHEME)
    wilcoxon_pairs = []
    for scheme in schemes_sorted:
        if _scheme_base(scheme) == ref_base:
            continue
        cell = _scheme_cell(scheme)
        ref_k = None
        for k in schemes_sorted:
            if _scheme_base(k) == ref_base and _scheme_cell(k) == cell:
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
            if _scheme_cell(sa) != _scheme_cell(sb):
                continue
            wilcoxon_pairs.append((sa, sb))

    if wilcoxon_pairs and wilcoxon_metrics and n_seeds >= 2:
        emit("")
        emit(
            "4. PAIRED WILCOXON SIGNED-RANK  (by seed_idx, within one dataset and minimum-hit-groups; exact for n≤20 non-zero pairs)"
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
        "(xor order; greedymini Γ / mod-minimizer tb XOR tie-breaks)."
    )
    emit(
        "  • n={} seeds/scheme; each (seed, dataset, minimum-hit-groups) cell is classified "
        "once: the replicates are the seeds (±std is across seeds). Wall = DB load + "
        "processing; Proc = the classify \"processed in\" time.".format(n_seeds)
    )
    emit(
        "  • minimum-hit-groups 0 reproduces the Kraken 2 paper setting; 2 is the k2 default."
    )
    emit("  • MAPE / concordance compare to the reference scheme of the same seed, dataset and minimum-hit-groups.")
    emit("  • Lower MAPE / higher concordance ⇒ closer to {} under that seed.".format(
        ref_name or "reference"))
    emit(
        "  • ClassConc% = identical taxid among reads classified by both; "
        "Jaccard% = |ref∩other| / |ref∪other| over sets of classified read IDs; "
        "NewClass%/NewUnclass% = ref U→C / C→U as % of all reads."
    )
    emit(
        "  • Wilcoxon (section 4) pairs schemes by seed_idx on the per-seed classify wall."
    )
    if has_phase:
        emit(
            "  • Phase timing requires classify built with "
            "ENABLE_MINIMIZER_PROF=1 (-DKRAKEN2_MINIMIZER_PROF); files: "
            "classify_profile_<dataset>_mhg<N>.txt via env KRAKEN2_MINIMIZER_PROF."
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
    if acc_sheets:
        emit("")
        emit("5. ACCURACY vs LABELS  (Wood et al. 2019; both minimum-hit-groups settings)")
        emit(
            "  Strain-exclusion: genus/species sensitivity, PPV, F1 "
            "(simulated truth)."
        )
        emit(
            "  FDA-ARGOS: genus concordance only — (a) concordant (b) discordant "
            "(c) ancestor (d) unclassified; no sensitivity/PPV."
        )
        emit(
            "  CAMI (supporting set): sens/PPV/F1 over reads whose true taxon is in the library "
            "(excluded = no reference; reported, not scored)."
        )
        for name, fields, rows in acc_sheets:
            emit("")
            emit("  {} (n={})".format(name, len(rows)))
            if name == "accuracy_strain_exclusion":
                emit("    (mean/std over seeds in results/csvs/kraken2_accuracy_strain_exclusion.csv)")
                continue
            if name == "accuracy_fda_argos":
                emit("    (per-SRR scope rows in results/csvs/kraken2_accuracy_fda_argos.csv)")
                continue
            if name == "accuracy_cami":
                emit("    (mean/std over seeds in results/csvs/kraken2_{}.csv)".format(name))
                continue
            if name == "accuracy_strain_exclusion_summary":
                show = [c for c in fields]
                display = rows
            elif name == "accuracy_fda_argos_summary":
                show = [c for c in fields]
                display = rows
            elif name == "accuracy_cami_summary":
                show = [c for c in fields if not c.endswith("_std")]
                display = [r for r in rows if r.get("group") in ("all", "bacteria", "archaea")]
            else:
                show = [
                    c
                    for c in (
                        "scheme",
                        "minimum_hit_groups",
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

            def _cell(v, w=14):
                if isinstance(v, float):
                    v = "{:.6g}".format(v)
                return "{:<{w}}".format("" if v is None else v, w=w)[:w]

            def _w(c):
                return 26 if c == "scheme" else 14

            emit("    " + "  ".join(_cell(c, _w(c)) for c in show))
            for row in display:
                emit("    " + "  ".join(_cell(row.get(c, ""), _w(c)) for c in show))

    emit("Output directories:")
    for r in runs:
        seed = r["seed"] or "?"
        emit("  {} (seed {}) -> {}".format(r["context"], seed, r["outdir"]))

    xlsx_path = csv_dir = None
    try:
        library = os.path.basename(os.path.normpath(base))
        acc_rows = {}
        for name, _fields, rows in acc_sheets or []:
            kind = ACCURACY_KIND.get(name)
            if kind:
                acc_rows[kind] = rows
        qs_sheet = qs[2] if qs else []
        xlsx_path, csv_dir, _tables = export_layout.export_library(
            base,
            library,
            collect_per_seed_rows(by_scheme, metrics_by_scheme),
            acc_rows,
            load_density_rows(base),
            qs_sheet,
        )
        pc = _tables.get("paper_check") if library == export_layout.PAPER_LIBRARY else None
        if pc:
            emit("")
            emit("6. PAPER CHECK  (xor, minimum-hit-groups 0, strain_exclusion vs Wood et al. 2019; "
                 "tolerance {} absolute)".format(export_layout.PAPER_TOLERANCE))
            emit("  {:<10} {:<11} {:<8} {:<26} {:>9} {:>9} {:>9}  {}".format(
                "table", "domain", "rank", "metric", "ours", "published", "diff", "within"))
            for r in pc[1]:
                emit("  {:<10} {:<11} {:<8} {:<26} {:>9} {:>9} {:>9}  {}".format(
                    r["table"], r.get("domain") or "-", r["rank"], r["metric"],
                    fmt(r.get("ours"), 4), fmt(r.get("published"), 4),
                    fmt(r.get("difference"), 4), r.get("within_tolerance") or "-"))
    except Exception as e:
        import traceback

        traceback.print_exc()
        emit("")
        emit("Could not write CSV/XLSX tables: {}".format(e))

    try:
        with open(out_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        emit("")
        emit("Wrote {}".format(out_path))
        if xlsx_path:
            emit("Wrote {} (all libraries)".format(xlsx_path))
            emit("Wrote CSVs (kraken2_<sheet>.csv) under {}".format(csv_dir))
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
