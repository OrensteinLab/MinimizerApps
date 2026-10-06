#!/usr/bin/env python3
"""Helpers for the per-seed density table (tables/minimizer_density.csv).

  density_csv.py missing CSV SEED_IDX TOGGLE_HEX L K SPACES SCHEME...
      print the schemes that have no row for (scheme, seed_idx) with the same
      toggle / l / k / spaces in CSV (one per line)
  density_csv.py merge CSV NEW_CSV SPACES
      replace rows of the same (scheme, seed_idx) in CSV by the rows of
      NEW_CSV (written by kraken2/src/benchmark_minimizers -o X.csv); keeps
      all other rows; adds a `spaces` column; sorted by seed_idx, scheme
"""

from __future__ import print_function

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_const import scheme_sort_key  # noqa: E402


def read_rows(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path) as f:
        r = csv.DictReader(f)
        return list(r), list(r.fieldnames or [])


def same_num(a, b):
    try:
        return int(a, 0) == int(b, 0)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


def cmd_missing(argv):
    path, seed_idx, toggle, l, k, spaces = argv[:6]
    schemes = argv[6:]
    rows, _ = read_rows(path)
    for s in schemes:
        ok = False
        for r in rows:
            if (r.get("scheme") == s and same_num(r.get("seed_idx"), seed_idx)
                    and same_num(r.get("toggle_hex"), toggle)
                    and same_num(r.get("l"), l) and same_num(r.get("k"), k)
                    and same_num(r.get("spaces", spaces), spaces)
                    and r.get("fac_particular") not in (None, "")):
                ok = True
                break
        if not ok:
            print(s)


def cmd_merge(argv):
    path, new_path, spaces = argv[:3]
    old, old_fields = read_rows(path)
    new, new_fields = read_rows(new_path)
    if not new:
        sys.exit("density_csv.py: no rows in %s" % new_path)
    fields = list(new_fields)
    if "spaces" not in fields:
        fields.insert(fields.index("w") + 1 if "w" in fields else len(fields),
                      "spaces")
    for f in old_fields:
        if f not in fields:
            fields.append(f)
    for r in new:
        r.setdefault("spaces", spaces)
    keys = set((r["scheme"], int(r["seed_idx"])) for r in new)
    kept = [r for r in old
            if (r.get("scheme"), int(r.get("seed_idx") or 0)) not in keys]
    allrows = kept + new
    allrows.sort(key=lambda r: (int(r.get("seed_idx") or 0),
                                scheme_sort_key(r["scheme"])))
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n",
                           restval="")
        w.writeheader()
        for r in allrows:
            w.writerow(r)
    os.replace(tmp, path)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("missing", "merge"):
        sys.exit(__doc__)
    (cmd_missing if sys.argv[1] == "missing" else cmd_merge)(sys.argv[2:])


if __name__ == "__main__":
    main()
