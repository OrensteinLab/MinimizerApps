#!/usr/bin/env python3
"""Helpers for the per-seed density table (results/runs/<lib>/minimizer_density.csv).

  density_csv.py missing CSV SEED_IDX TOGGLE_HEX L K SPACES PROJECTION GM_FILE MOD_GM_FILE SCHEME...
      print the schemes (canonical scheme names) that have no row for the same
      (scheme, gm_k, projection, gm_table, mod_t, seed_idx) with the same
      toggle / l / k / spaces in CSV (one per line); a mod-minimizer row
      counts only with the current tie rule (mod_tie_rule = MOD_TIE_RULE)
  density_csv.py merge CSV NEW_CSV SPACES L K PROJECTION GM_FILE MOD_GM_FILE
      replace rows of the same (scheme, gm fields, seed_idx) in CSV by the rows
      of NEW_CSV (written by kraken2/src/benchmark_minimizers -o X.csv); keeps
      all other rows; adds the `spaces` column and the identity columns
      gm_k, projection, gm_table, mod_t (benchmark/common/scheme_names.py)
      and, on mod-minimizer rows, mod_tie_rule; sorted by seed_idx, scheme
  Columns use the unified names k = minimizer length (the tool's -l), w,
  L = k + w - 1 (the tool's -k); files with the old l / k header are mapped on
  read and rewritten with the new names.
"""

from __future__ import print_function

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_const import scheme_id, scheme_sort_key  # noqa: E402
import scheme_names  # noqa: E402

ID_COLS = ("gm_k", "projection", "gm_table", "mod_t")
# Tie rule of the mod-minimizer the rows were computed with; rows without it predate the center rule
# (G/T-strand rule) and are recomputed. Same stamp as benchmark/sshash/scripts/expected_density.py.
MOD_TIE_RULE = "center"


def is_mod(scheme):
    return (scheme or "").startswith("mod-minimizer")


def ident(scheme, projection, gm_file, mod_gm_file, l, k):
    """Identity columns (text) of a density row of `scheme`."""
    gm = {"greedymini": gm_file, "mod-minimizer-greedymini": mod_gm_file}.get(scheme)
    idd = scheme_id(scheme, gm or None, projection or None, int(l), int(k))
    return scheme_names.format_id(idd)


def row_key(r):
    return (r.get("scheme"), r.get("gm_k") or "", r.get("projection") or "",
            r.get("gm_table") or "", r.get("mod_t") or "", int(r.get("seed_idx") or 0))


def unify_fields(fields):
    """Legacy / tool header (l = minimizer length, k = k-mer length, w) -> the
    unified names k = minimizer length, w, L = k + w - 1 (k-mer length). A header
    that already has `L` is left alone."""
    if "L" in fields or "l" not in fields:
        return list(fields)
    return [{"l": "k", "k": "L"}.get(f, f) for f in fields]


def read_rows(path):
    """(rows, fields) with the unified names (see unify_fields)."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path) as f:
        r = csv.DictReader(f)
        old = list(r.fieldnames or [])
        new = unify_fields(old)
        rows = [dict(zip(new, [row.get(o) for o in old])) for row in r]
        return rows, new


def same_num(a, b):
    try:
        return int(a, 0) == int(b, 0)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


def cmd_missing(argv):
    # argv l / k are the tool's -l / -k (minimizer length, k-mer length)
    path, seed_idx, toggle, l, k, spaces, projection, gm_file, mod_gm_file = argv[:9]
    schemes = argv[9:]
    rows, _ = read_rows(path)
    for s in schemes:
        ok = False
        want = ident(s, projection, gm_file, mod_gm_file, l, k)
        for r in rows:
            if (r.get("scheme") == s and same_num(r.get("seed_idx"), seed_idx)
                    and all((r.get(c) or "") == want[c] for c in ID_COLS)
                    and same_num(r.get("toggle_hex"), toggle)
                    and same_num(r.get("k"), l) and same_num(r.get("L"), k)
                    and same_num(r.get("spaces", spaces), spaces)
                    and (not is_mod(s) or r.get("mod_tie_rule") == MOD_TIE_RULE)
                    and r.get("fac_particular") not in (None, "")):
                ok = True
                break
        if not ok:
            print(s)


def cmd_merge(argv):
    path, new_path, spaces, l, k, projection, gm_file, mod_gm_file = argv[:8]
    old, old_fields = read_rows(path)
    new, new_fields = read_rows(new_path)
    if not new:
        sys.exit("density_csv.py: no rows in %s" % new_path)
    fields = list(new_fields)
    if "spaces" not in fields:
        fields.insert(fields.index("L") + 1 if "L" in fields else len(fields),
                      "spaces")
    for c in ID_COLS:
        if c not in fields:
            fields.insert(fields.index("scheme") + 1 + ID_COLS.index(c), c)
    for f in old_fields:
        if f not in fields:
            fields.append(f)
    if "mod_tie_rule" not in fields:
        fields.append("mod_tie_rule")
    if all(c in fields for c in ("k", "w", "L")):      # k, w, L together
        fields = [f for f in fields if f not in ("w", "L")]
        i = fields.index("k") + 1
        fields[i:i] = ["w", "L"]
    for r in new:
        r.setdefault("spaces", spaces)
        idd = ident(r["scheme"], projection, gm_file, mod_gm_file, l, k)
        r.update((c, idd[c]) for c in ID_COLS)
        r["mod_tie_rule"] = MOD_TIE_RULE if is_mod(r["scheme"]) else ""
    keys = set(row_key(r) for r in new)
    kept = [r for r in old if row_key(r) not in keys]
    allrows = kept + new
    allrows.sort(key=lambda r: (int(r.get("seed_idx") or 0),
                                scheme_sort_key(r["scheme"]),
                                r.get("projection") or "", r.get("gm_k") or ""))
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
