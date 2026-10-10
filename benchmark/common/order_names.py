"""Canonical minimizer-order names and identity fields shared by kraken2, kmc,
sshash and minimap2 (stdlib only, Python 3.6 compatible).

An order is identified by FIVE separate fields (never packed into one string):

  order       one of CANONICAL (plus an app baseline: native / default / lex)
  gm_k        GreedyMini table k (int) or None        (greedymini orders only)
  projection  'at0-cg1' (A,T->0 C,G->1; default) | 'at1-cg0' (flipped) | None
  gm_table    table file stem, e.g. 'w11_k15_canonical' | None
  mod_t       mod-minimizer t-mer length (int) | None (mod-minimizer orders)

ID_FIELDS lists them in that column order. Canonical names:

  xor, xor-pcg-hash, greedymini, mod-minimizer-xor, mod-minimizer-greedymini

App baselines are not canonical orders; accept them by passing
`baselines=BASELINES` (or a subset) to the validating functions.

Public API
----------
  is_greedymini(order) -> bool        greedymini or mod-minimizer-greedymini
  is_mod_minimizer(order) -> bool
  validate_order(order, baselines=()) -> order   (ValueError if unknown)
  projection_from_flag(flag) -> 'at0-cg1' | 'at1-cg0'
      flag: True/False (True = at1), 'at1'/'at0', 'at1-cg0'/'at0-cg1', 1/0.
  parse_table_name(name) -> {'gm_table': stem, 'gm_k': int, 'gm_w': int|None}
      from 'w11_k15_canonical.gm' / a path / a bare stem.
  gm_fields(table=None, projection=None, gm_k=None) -> dict(gm_k, projection,
      gm_table): derive the gm columns from a table filename and/or a
      projection flag (projection defaults to 'at0-cg1' when a table or gm_k
      is given).
  make_id(order, gm_k=None, projection=None, gm_table=None, mod_t=None,
          baselines=(), default_projection=True) -> dict over ID_FIELDS
      validated + normalised (non-greedymini orders get gm_k/projection/
      gm_table None; non-mod orders get mod_t None).
  format_id(idd) -> {field: str}   CSV text ('' for None)
  parse_id(strs, baselines=()) -> dict   inverse of format_id (typed)
  order_dir_name(order, gm_k=None, projection=None, mod_t=None) -> str
  parse_order_dir(path) -> dict over ID_FIELDS (inverse, gm_table None)

Run-directory scheme (order_dir_name): a stable relative POSIX path
  <order>[/<params>]   params = '_'.join of  k<gm_k>, <projection>, t<mod_t>
  (those that apply, in that order). Examples:
    xor                                    native
    greedymini/k15_at0-cg1                 greedymini/k15_at1-cg0
    mod-minimizer-xor/t10                  mod-minimizer-greedymini/k6_at0-cg1_t6
Apps append their own parameter segments (k, w, seed, ...) after this.
"""

from __future__ import print_function

import os
import re

CANONICAL = ("xor", "xor-pcg-hash", "greedymini", "mod-minimizer-xor",
             "mod-minimizer-greedymini")
BASELINES = ("native", "default", "lex")
GREEDYMINI_ORDERS = ("greedymini", "mod-minimizer-greedymini")
MOD_ORDERS = ("mod-minimizer-xor", "mod-minimizer-greedymini")
PROJECTIONS = ("at0-cg1", "at1-cg0")
DEFAULT_PROJECTION = "at0-cg1"
ID_FIELDS = ("order", "gm_k", "projection", "gm_table", "mod_t")


def is_greedymini(order):
    return order in GREEDYMINI_ORDERS


def is_mod_minimizer(order):
    return order in MOD_ORDERS


def validate_order(order, baselines=()):
    if order in CANONICAL or order in tuple(baselines):
        return order
    raise ValueError("unknown order %r (canonical: %s; baselines allowed: %s)"
                     % (order, ", ".join(CANONICAL), ", ".join(baselines) or "none"))


def projection_from_flag(flag):
    """Normalise a projection flag to 'at0-cg1' / 'at1-cg0'."""
    if flag is None:
        return DEFAULT_PROJECTION
    if isinstance(flag, bool) or isinstance(flag, int):
        return "at1-cg0" if flag else "at0-cg1"
    f = str(flag).strip().lower()
    if f in ("at1", "at1-cg0", "1", "true"):
        return "at1-cg0"
    if f in ("at0", "at0-cg1", "0", "false", ""):
        return "at0-cg1"
    raise ValueError("unknown projection %r (use at0-cg1 or at1-cg0)" % (flag,))


_TABLE_RE = re.compile(r"^(?:w(\d+)_)?k(\d+)(?:_[A-Za-z0-9]+)*$")


def parse_table_name(name):
    """'w11_k15_canonical.gm' (or a path / bare stem) -> gm_table/gm_k/gm_w."""
    stem = os.path.basename(str(name))
    for ext in (".gm", ".txt"):
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
    m = _TABLE_RE.match(stem)
    if not m:
        raise ValueError("cannot parse a GreedyMini table name from %r" % (name,))
    return {"gm_table": stem, "gm_k": int(m.group(2)),
            "gm_w": int(m.group(1)) if m.group(1) else None}


def gm_fields(table=None, projection=None, gm_k=None):
    """gm_k / projection / gm_table from a table filename and/or projection flag."""
    out = {"gm_k": gm_k, "projection": None, "gm_table": None}
    if table:
        p = parse_table_name(table)
        out["gm_table"] = p["gm_table"]
        if gm_k is not None and gm_k != p["gm_k"]:
            raise ValueError("gm_k %s disagrees with table %s" % (gm_k, table))
        out["gm_k"] = p["gm_k"]
    if projection is not None or table or gm_k is not None:
        out["projection"] = projection_from_flag(projection)
    return out


def make_id(order, gm_k=None, projection=None, gm_table=None, mod_t=None,
            baselines=(), default_projection=True):
    """Validated, normalised identity dict over ID_FIELDS."""
    validate_order(order, baselines)
    if is_greedymini(order):
        if gm_table and gm_k is None:
            gm_k = parse_table_name(gm_table)["gm_k"]
        if gm_table:
            gm_table = parse_table_name(gm_table)["gm_table"]
        if projection is not None:
            projection = projection_from_flag(projection)
        elif default_projection:
            projection = DEFAULT_PROJECTION
        gm_k = None if gm_k is None else int(gm_k)
    else:
        gm_k = projection = gm_table = None
    mod_t = int(mod_t) if (is_mod_minimizer(order) and mod_t not in (None, "")) else None
    return {"order": order, "gm_k": gm_k, "projection": projection,
            "gm_table": gm_table or None, "mod_t": mod_t}


def format_id(idd):
    return dict((f, "" if idd.get(f) is None else str(idd[f])) for f in ID_FIELDS)


def parse_id(strs, baselines=()):
    def num(v):
        return None if v in (None, "") else int(v)

    def txt(v):
        return None if v in (None, "") else str(v)

    return make_id(strs["order"], num(strs.get("gm_k")), txt(strs.get("projection")),
                   txt(strs.get("gm_table")), num(strs.get("mod_t")), baselines,
                   default_projection=False)


def order_dir_name(order, gm_k=None, projection=None, mod_t=None):
    """Stable run-dir path '<order>[/k<gm_k>_<projection>_t<mod_t>]'."""
    parts = []
    if is_greedymini(order):
        if gm_k is not None:
            parts.append("k%d" % int(gm_k))
        parts.append(projection_from_flag(projection))
    if is_mod_minimizer(order) and mod_t is not None:
        parts.append("t%d" % int(mod_t))
    return order + "/" + "_".join(parts) if parts else order


def parse_order_dir(path):
    """Inverse of order_dir_name (gm_table is not encoded -> None)."""
    segs = str(path).strip("/").split("/")
    order = segs[0]
    gm_k = projection = mod_t = None
    if len(segs) > 1:
        for tok in segs[1].split("_"):
            if re.match(r"^k\d+$", tok):
                gm_k = int(tok[1:])
            elif re.match(r"^t\d+$", tok):
                mod_t = int(tok[1:])
            elif tok in PROJECTIONS:
                projection = tok
            else:
                raise ValueError("bad order dir segment %r" % (segs[1],))
    return {"order": order, "gm_k": gm_k, "projection": projection,
            "gm_table": None, "mod_t": mod_t}
