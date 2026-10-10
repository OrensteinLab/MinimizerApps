"""Order name -> sshash CLI flags (see ../../../orders/README.md).

Orders carry only their canonical name (benchmark/common/order_names.py):
native (sshash's own baseline), xor, xor-pcg-hash, greedymini,
mod-minimizer-xor, mod-minimizer-greedymini. A GreedyMini order is run once per
variant -- a (table, projection) pair from config.yaml's `greedymini` /
`mod_minimizer_greedymini` sections -- and each variant carries its own identity
fields gm_k / projection / gm_table / mod_t, never a suffix of the order name.

Every order (including ``native``) takes a seed with sshash -- unlike kmc,
where ``lex`` has no seed and is special-cased. So there is no "no seed"
branch here: every variant is run once per seed the harness sweeps.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import order_names as on
from scripts.config import gm_choices, gm_table_file

BASELINES = ("native",)
ORDERS = ("native", *on.CANONICAL)
# orders whose only parameter is the seed (flags: -O ORDER -s SEED)
_SEED_ONLY = ("xor", "xor-pcg-hash")


def variant_key(v: dict) -> tuple:
    """Hashable identity of a variant (its ID_FIELDS values)."""
    return tuple(v.get(f) for f in on.ID_FIELDS)


def table_gm_k(path: Path) -> int:
    """gm_k of a table file: its k from the file name, checked against the
    table size (2^gm_k entries)."""
    from scripts.density import read_gm_ranks, _gm_k_from_table

    k = on.parse_table_name(path.name)["gm_k"]
    if path.is_file() and _gm_k_from_table(read_gm_ranks(path)) != k:
        raise ValueError(f"{path.name}: the file name says gm_k={k} but the table is not 2^{k} entries")
    return k


def xor_mod_t(k: int, m: int) -> int:
    """t-mer length of mod-minimizer-xor at (k, m): 4 + (m - 4) mod w."""
    return 4 + (m - 4) % (k - m + 1)


def expand_variants(orders, k: int, m: int, tables=None, projections=None) -> list[dict]:
    """One dict per run variant of `orders` (canonical names or native), in
    ORDERS order: the ID_FIELDS plus "gm_file" (path of the table or None).
    `tables` / `projections` override config.yaml's choices for the GreedyMini
    orders (table file names or stems; at0-cg1 / at1-cg0)."""
    wanted = [o.strip() for o in orders if o.strip()]
    for o in wanted:
        on.validate_order(o, BASELINES)
    out: list[dict] = []
    for order in ORDERS:
        if order not in wanted:
            continue
        if not on.is_greedymini(order):
            mod_t = xor_mod_t(k, m) if order == "mod-minimizer-xor" else None
            v = on.make_id(order, mod_t=mod_t, baselines=BASELINES)
            v["gm_file"] = None
            out.append(v)
            continue
        cfg_tables, cfg_projs = gm_choices(order)
        for table in (tables if tables else cfg_tables):
            path = gm_table_file(table)
            gm_k = table_gm_k(path)
            for proj in (projections if projections else cfg_projs):
                v = on.make_id(order, gm_k, proj, path.name, gm_k if order == "mod-minimizer-greedymini" else None,
                               baselines=BASELINES)
                v["gm_file"] = path
                out.append(v)
    return out


def check_gm_order(order: str, m: int, gm_file: Path) -> str | None:
    """Error message if the plain greedymini order cannot run at minimizer length m
    with the table `gm_file` (None if fine or not such an order): needs m - gm_k even
    and >= 0, and a mirror-symmetric table."""
    if order != "greedymini":
        return None
    from scripts.density import gm_table_mirror_symmetric, read_gm_ranks, _gm_k_from_table

    ranks = read_gm_ranks(gm_file)
    gm_k = _gm_k_from_table(ranks)
    if m < gm_k or (m - gm_k) % 2:
        return (
            f"{order}: needs m - gm_k even and >= 0, but m={m} and the table "
            f"{gm_file.name} has gm_k={gm_k} (m - gm_k = {m - gm_k})"
        )
    if not gm_table_mirror_symmetric(ranks):
        return f"{order}: table {gm_file.name} is not mirror-symmetric (rank[key] != rank[reverse(key)])"
    return None


def mod_t(order: str, k: int, m: int, gm_file: Path | str | None = None) -> int:
    """The t-mer length of a mod-minimizer order at (k, m): xor -> 4 + (m - 4) mod w,
    greedymini -> the table's gm_k."""
    if order == "mod-minimizer-xor":
        return xor_mod_t(k, m)
    if order != "mod-minimizer-greedymini":
        raise ValueError(f"{order} is not a mod-minimizer order")
    path = Path(gm_file) if gm_file else gm_table_file(gm_choices(order)[0][0])
    return table_gm_k(path)


def order_dir_name(order: str, gm_k: int | None = None, projection: str | None = None,
                   mod_t_: int | None = None) -> str:
    """Run-directory segment of a variant (see order_names.order_dir_name)."""
    return on.order_dir_name(order, gm_k, projection, mod_t_)


def build_flags(order: str, seed: int, gm_file: str | Path | None = None, k: int = 31,
                m: int = 21, projection: str | None = None) -> list[str]:
    """sshash `build` flags for `order` with `seed` (every order takes -s) at (k, m).
    GreedyMini orders use `gm_file` (default: the first table config.yaml registers)
    and `projection` (default at0-cg1; only a non-default one is passed)."""
    on.validate_order(order, BASELINES)
    flags: list[str] = ["-s", str(seed)]
    if order == "native":
        return flags
    if order in _SEED_ONLY:
        return ["-O", order, *flags]
    if order == "mod-minimizer-xor":
        return ["-O", order, "--mod-t", str(xor_mod_t(k, m)), *flags]
    path = Path(gm_file) if gm_file else gm_table_file(gm_choices(order)[0][0])
    if not path.is_file():
        raise FileNotFoundError(f"GreedyMini table not found: {path}")
    out = ["-O", order, "--greedy-mini-file", str(path)]
    if order == "mod-minimizer-greedymini":
        out += ["--mod-t", str(table_gm_k(path))]
    proj = on.projection_from_flag(projection)
    if proj != on.DEFAULT_PROJECTION:
        out += ["--gm-projection", proj]
    return [*out, *flags]


def _csv(s: str | None) -> list[str] | None:
    return [x.strip() for x in s.split(",") if x.strip()] if s else None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts orders", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_flags = sub.add_parser("flags", help="print sshash build flags for an order/seed")
    p_flags.add_argument("order", choices=ORDERS)
    p_flags.add_argument("seed", type=lambda s: int(s, 0))
    p_flags.add_argument("--gm-file", default=None)
    p_flags.add_argument("--projection", default=None)
    p_flags.add_argument("--k", type=int, default=31)
    p_flags.add_argument("--m", type=int, default=21)

    p_t = sub.add_parser("mod-t", help="print the t-mer length of a mod-minimizer order at (k, m)")
    p_t.add_argument("order")
    p_t.add_argument("k", type=int)
    p_t.add_argument("m", type=int)
    p_t.add_argument("--gm-file", default=None)

    p_dir = sub.add_parser("dir-name", help="print the run-directory segment of a variant")
    p_dir.add_argument("order")
    p_dir.add_argument("--gm-k", type=int, default=None)
    p_dir.add_argument("--projection", default=None)
    p_dir.add_argument("--mod-t", type=int, default=None)

    sub.add_parser("list", help="list canonical order names (plus the native baseline)")

    p_var = sub.add_parser(
        "variants",
        help="tab-separated run variants of an order list: order, gm_file, projection, gm_k, "
        "gm_table, mod_t ('-' = empty)")
    p_var.add_argument("orders", help="comma-separated order names")
    p_var.add_argument("--k", type=int, required=True)
    p_var.add_argument("--m", type=int, required=True)
    p_var.add_argument("--gm-tables", default=None, help="comma-separated table names (override config.yaml)")
    p_var.add_argument("--gm-projections", default=None, help="comma-separated projections (override config.yaml)")

    p_chk = sub.add_parser("check-gm", help="validate a greedymini table for minimizer length m (exit 1 with a message if not)")
    p_chk.add_argument("order")
    p_chk.add_argument("m", type=int)
    p_chk.add_argument("gm_file")

    args = ap.parse_args(argv)
    try:
        return _dispatch(args)
    except (ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


def _dispatch(args: argparse.Namespace) -> int:
    if args.cmd == "flags":
        print(" ".join(build_flags(args.order, args.seed, args.gm_file, args.k, args.m, args.projection)))
        return 0
    if args.cmd == "mod-t":
        if not on.is_mod_minimizer(args.order):
            print(f"error: {args.order} is not a mod-minimizer order", file=sys.stderr)
            return 1
        print(mod_t(args.order, args.k, args.m, args.gm_file))
        return 0
    if args.cmd == "dir-name":
        print(order_dir_name(args.order, args.gm_k, args.projection, args.mod_t))
        return 0
    if args.cmd == "check-gm":
        msg = check_gm_order(args.order, args.m, Path(args.gm_file))
        if msg:
            print(f"error: {msg}", file=sys.stderr)
            return 1
        return 0
    if args.cmd == "variants":
        orders = _csv(args.orders) or []
        for v in expand_variants(orders, args.k, args.m, _csv(args.gm_tables), _csv(args.gm_projections)):
            cells = [v["order"], v["gm_file"], v["projection"], v["gm_k"], v["gm_table"], v["mod_t"]]
            print("\t".join("-" if c is None else str(c) for c in cells))
        return 0
    if args.cmd == "list":
        for o in ORDERS:
            print(o)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
