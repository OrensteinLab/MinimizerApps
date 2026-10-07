#!/usr/bin/env python3
"""Per-read seeding/chaining work and time for one (index, read sample) pair,
from minimap2's own debug output (single thread, so lines of one read are
contiguous):

  --print-qname   QR <read> <tid> <len>        before mapping a read
                  QM <read> <len> cap=..,..    thread-local (kalloc) memory after it
                  QT <read> <tid> <seconds>    wall time to map it
  --dbg-seed-occ  SF <qpos> <n> <flt>          one line per read minimizer found in
                                                the index: n = reference occurrences,
                                                flt = 1 if dropped by mid_occ (seed.c)

Columns written per read: read, qlen, time_sec, km_cap_bytes, and (with
--seeds) n_seeds (read minimizers found in the index), n_filtered, anchors
(sum of n over kept seeds = chaining input size, map.c collect_seed_hits),
occ_all (sum of n over all seeds, before filtering).

Usage: seed_stats.py --minimap2 BIN --index IDX --reads FA --out OUT.csv [--seeds] [--cap N]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seeds", action="store_true", help="also collect per-seed occurrence lines (slower; don't use for timing)")
    ap.add_argument("--cap", default="", help="fixed mid_occ (-f N); empty = minimap2's default")
    ap.add_argument("--preset", default="map-pb", help="minimap2 -x preset (the dataset's, as in the full runs)")
    args = ap.parse_args(argv)

    cmd = [args.minimap2, "-x", args.preset, "-c", "-t", "1", "--print-qname"]
    if args.seeds:
        cmd.append("--dbg-seed-occ")
    if args.cap:
        cmd += ["-f", args.cap]
    cmd += [args.index, args.reads]

    rows: list[dict] = []
    cur: dict | None = None
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stderr is not None
    for line in proc.stderr:
        tag = line[:3]
        if tag == "SF\t" and cur is not None:
            _, _qpos, n, flt = line.split("\t")
            n = int(n)
            cur["n_seeds"] += 1
            cur["occ_all"] += n
            if flt.strip() == "1":
                cur["n_filtered"] += 1
            else:
                cur["anchors"] += n
        elif tag == "QR\t":
            _, name, _tid, qlen = line.rstrip("\n").split("\t")
            cur = {"read": name, "qlen": int(qlen), "time_sec": None, "km_cap_bytes": None,
                   "n_seeds": 0, "n_filtered": 0, "anchors": 0, "occ_all": 0}
        elif tag == "QM\t" and cur is not None:
            cap = line.split("cap=", 1)[1].split(",", 1)[0]
            cur["km_cap_bytes"] = int(cap)
        elif tag == "QT\t" and cur is not None:
            cur["time_sec"] = float(line.rstrip("\n").split("\t")[3])
            rows.append(cur)
            cur = None
    if proc.wait() != 0:
        print(f"error: {' '.join(cmd)} exited {proc.returncode}", file=sys.stderr)
        return 1

    fields = ["read", "qlen", "time_sec", "km_cap_bytes"]
    if args.seeds:
        fields += ["n_seeds", "n_filtered", "anchors", "occ_all"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} reads)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
