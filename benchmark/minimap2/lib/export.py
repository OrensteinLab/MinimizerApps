#!/usr/bin/env python3
"""Aggregate every runs/**/meta.json into minimap2_runs.csv, mirroring the
per-run-then-aggregate shape of sshash's/kmc's own export.py (a lighter
version: CSV only, no xlsx, since this harness has 5-7 orders x 1 seed each
by default rather than a full seed sweep)."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FIELDS = [
    "dataset", "order", "seed", "k", "w", "mod_t",
    "distinct_minimizers", "density_factor",
    "build_wall_time_sec", "build_peak_rss_kb",
    "map_wall_time_sec", "map_peak_rss_kb",
    "n_reads", "n_reads_in_repeat", "frac_mapped",
    "frac_correct", "frac_correct_in_repeat", "frac_incorrectly_mapped",
    "frac_correct_strict50", "frac_correct_strict50_in_repeat",
    "mean_identity",
]


def flatten(meta: dict) -> dict:
    acc = meta.get("accuracy") or {}
    build = meta.get("build") or {}
    mapp = meta.get("map") or {}
    return {
        "dataset": meta.get("dataset"),
        "order": meta.get("order"),
        "seed": meta.get("seed"),
        "k": meta.get("k"),
        "w": meta.get("w"),
        "mod_t": meta.get("mod_t"),
        "distinct_minimizers": meta.get("distinct_minimizers"),
        "density_factor": meta.get("density_factor"),
        "build_wall_time_sec": build.get("wall_time_sec"),
        "build_peak_rss_kb": build.get("peak_rss_kb"),
        "map_wall_time_sec": mapp.get("wall_time_sec"),
        "map_peak_rss_kb": mapp.get("peak_rss_kb"),
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    rows = []
    for meta_path in sorted((args.results_dir / "runs").glob("**/meta.json")):
        meta = json.loads(meta_path.read_text())
        rows.append(flatten(meta))

    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"wrote {args.out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
