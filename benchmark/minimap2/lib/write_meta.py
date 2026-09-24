#!/usr/bin/env python3
"""Assemble one run's meta.json from its /usr/bin/time -v logs, the parsed
minimizer-density count, and score_accuracy.py's JSON output. Mirrors the
shape of sshash's/kmc's own per-run meta.json (env-var-driven there via
scripts/meta.py; a plain CLI here since this harness has no Python-side
per-run driver yet)."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_time_v(path: Path) -> dict:
    if not path.is_file():
        return {}
    text = path.read_text()
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time.*?: (\S+)", text)
    if m:
        out["wall_time_str"] = m.group(1)
        out["wall_time_sec"] = _hms_to_sec(m.group(1))
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if m:
        out["peak_rss_kb"] = int(m.group(1))
    m = re.search(r"Percent of CPU this job got: (\d+)%", text)
    if m:
        out["cpu_percent"] = int(m.group(1))
    return out


def _hms_to_sec(s: str) -> float:
    parts = s.split(":")
    parts = [float(p) for p in parts]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, sec = parts[-3:]
    return h * 3600 + m * 60 + sec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", required=True)
    ap.add_argument("--seed")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--w", type=int, required=True)
    ap.add_argument("--mod-t", type=int, default=0)
    ap.add_argument("--dataset", default="smoke")
    ap.add_argument("--build-time", type=Path, required=True)
    ap.add_argument("--map-time", type=Path, required=True)
    ap.add_argument("--num-kmers", type=int, default=None, help="denominator for density_factor: genome length")
    ap.add_argument("--density", type=int, default=None, help="distinct minimizers reported by minimap2 -d")
    ap.add_argument("--accuracy", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    meta = {
        "dataset": args.dataset,
        "order": args.order,
        "seed": args.seed,
        "k": args.k,
        "w": args.w,
        "mod_t": args.mod_t or None,
        "build": parse_time_v(args.build_time),
        "map": parse_time_v(args.map_time),
        "distinct_minimizers": args.density,
    }
    if args.num_kmers and args.density:
        meta["density_factor"] = round(args.density / args.num_kmers * (args.w + 1), 6)
    if args.accuracy.is_file():
        meta["accuracy"] = json.loads(args.accuracy.read_text())

    args.out.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
