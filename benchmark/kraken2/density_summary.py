#!/usr/bin/env python3
"""Summarise density_sweep.sh CSVs: mean +- std (sample, over seed_idx) of the
density factor per (variant, spaces, scheme, l, w).

  density_summary.py a.csv b.csv ... [--col fac_expected|fac_particular]
Prints a markdown table, one row per l, one column per (variant, scheme)."""
import csv, sys, statistics, collections

args = sys.argv[1:]
col = "fac_expected"
if "--col" in args:
    i = args.index("--col"); col = args[i + 1]; del args[i:i + 2]
groups = collections.defaultdict(list)
for path in args:
    for r in csv.DictReader(open(path)):
        groups[(r["spaces"], int(r["w"]), r["variant"], r["scheme"], int(r["l"]))].append(float(r[col]))
for sp, w in sorted({(g[0], g[1]) for g in groups}):
    cols = sorted({(g[2], g[3]) for g in groups if g[0] == sp and g[1] == w})
    ls = sorted({g[4] for g in groups if g[0] == sp and g[1] == w})
    print(f"\n{col}, spaces={sp}, w={w}\n")
    print("| l | " + " | ".join(f"{v} {s}" for v, s in cols) + " |")
    print("|---|" + "---|" * len(cols))
    for l in ls:
        cells = []
        for v, s in cols:
            x = groups.get((sp, w, v, s, l))
            if not x: cells.append("-")
            elif len(x) == 1: cells.append(f"{x[0]:.4f}")
            else: cells.append(f"{statistics.mean(x):.4f} +- {statistics.stdev(x):.4f} (n={len(x)})")
        print(f"| {l} | " + " | ".join(cells) + " |")
