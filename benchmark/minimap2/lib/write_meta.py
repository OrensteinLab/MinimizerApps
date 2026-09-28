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
    ap.add_argument("--seed-idx", type=int, default=0, help="index into config.yaml's shared seed list (0 for seed-independent orders: native, greedymini)")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--w", type=int, required=True)
    ap.add_argument("--mod-t", type=int, default=0)
    ap.add_argument("--dataset", default="smoke")
    ap.add_argument("--mid-occ-cap", default="", help="-f value passed at map time (empty = minimap2's own per-index default)")
    ap.add_argument("--mid-occ", default="", help="effective mid_occ minimap2 reported in its map log")
    ap.add_argument("--build-time", type=Path, required=True)
    ap.add_argument("--map-time", type=Path, required=True)
    ap.add_argument("--density", type=int, default=None, help="distinct minimizers reported by minimap2 -d (informational only -- NOT the density_factor numerator, see below)")
    ap.add_argument("--avg-occurrences", type=float, default=None, help="'average occurrences' from minimap2 -d's stat line (distinct minimizers undercount total selected positions whenever a minimizer repeats -- see mm_idx_stat in index.c)")
    ap.add_argument("--avg-spacing", type=float, default=None, help="'average spacing' from minimap2 -d's stat line: genome_len / total_occurrences -- the correct density_factor input")
    ap.add_argument("--accuracy", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    meta = {
        "dataset": args.dataset,
        "order": args.order,
        "seed": args.seed,
        "seed_idx": args.seed_idx,
        "k": args.k,
        "w": args.w,
        "mod_t": args.mod_t or None,
        "mid_occ_cap": int(args.mid_occ_cap) if args.mid_occ_cap else None,
        "mid_occ": int(args.mid_occ) if args.mid_occ else None,
        "build": parse_time_v(args.build_time),
        "map": parse_time_v(args.map_time),
        "distinct_minimizers": args.density,
        "avg_occurrences": args.avg_occurrences,
    }
    if args.avg_spacing:
        # density = total_occurrences / genome_len = 1 / avg_spacing (mm_idx_stat
        # prints avg_spacing = genome_len / total_occurrences directly);
        # density_factor = density * (w+1). NOT distinct_minimizers/genome_len:
        # "distinct minimizers" collapses every repeated occurrence of the same
        # minimizer key to one count, so it always undercounts total selected
        # positions whenever the genome has any repeat content -- and it
        # undercounts *more* for an order that's *better* at consistently
        # re-selecting the same minimizer across a repeat's copies (exactly
        # the property these custom orders are supposed to have), which can
        # push a distinct-minimizer-based ratio below the proven theoretical
        # minimum of density_factor=1 (every scheme must select at least
        # 1 per w+1-window). Using avg_spacing avoids this entirely.
        #
        # Same formula as benchmark/common/density.py's density_factor(n_selected,
        # n_kmers, w) (= n_selected/n_kmers*(w+1) = (1/avg_spacing)*(w+1)), just
        # applied to a different input shape -- minimap2's own build log exposes
        # the ratio directly as avg_spacing, not separate n_selected/n_kmers
        # counts, and re-deriving n_kmers from genome length here would be less
        # accurate (minimap2 skips N-runs when sketching, so genome length isn't
        # the true k-mer window count) for no benefit. Not a duplicate
        # implementation, so it isn't routed through that shared function.
        meta["density_factor"] = round((args.w + 1) / args.avg_spacing, 6)
    if args.accuracy.is_file():
        meta["accuracy"] = json.loads(args.accuracy.read_text())

    args.out.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
