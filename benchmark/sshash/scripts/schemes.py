"""Scheme name -> sshash CLI flags (see ../../../orders/README.md).

Schemes carry only their canonical name (benchmark/common/scheme_names.py):
native (sshash's own baseline), xor, xor-pcg-hash, greedymini,
mod-minimizer-xor, mod-minimizer-greedymini. A GreedyMini scheme is run once per
variant -- a (table, projection) pair from config.yaml's `greedymini` /
`mod_minimizer_greedymini` sections -- and each variant carries its own identity
fields gm_k / projection / gm_table / mod_t, never a suffix of the scheme name.

Every scheme (including ``native``) takes a seed with sshash -- unlike kmc,
where ``lex`` has no seed and is special-cased. So there is no "no seed"
branch here: every variant is run once per seed the harness sweeps.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import scheme_names as sn
from scripts.config import gm_choices, gm_table_file

BASELINES = ("native",)
SCHEMES = ("native", *sn.CANONICAL)
# schemes whose only parameter is the seed (flags: --scheme SCHEME -s SEED)
_SEED_ONLY = ("xor", "xor-pcg-hash")


def variant_key(v: dict) -> tuple:
    """Hashable identity of a variant (its ID_FIELDS values)."""
    return tuple(v.get(f) for f in sn.ID_FIELDS)


def table_gm_k(path: Path) -> int:
    """gm_k of a table file: its k from the file name, checked against the
    table size (2^gm_k entries)."""
    from scripts.density import read_gm_ranks, _gm_k_from_table

    k = sn.parse_table_name(path.name)["gm_k"]
    if path.is_file() and _gm_k_from_table(read_gm_ranks(path)) != k:
        raise ValueError(f"{path.name}: the file name says gm_k={k} but the table is not 2^{k} entries")
    return k


def xor_mod_t(k: int, m: int) -> int:
    """t-mer length of mod-minimizer-xor at (k, m): 4 + (m - 4) mod w."""
    return 4 + (m - 4) % (k - m + 1)


def expand_variants(schemes, k: int, m: int, tables=None, projections=None) -> list[dict]:
    """One dict per run variant of `schemes` (canonical names or native), in
    SCHEMES order: the ID_FIELDS plus "gm_file" (path of the table or None).
    `tables` / `projections` override config.yaml's choices for the GreedyMini
    schemes (table file names or stems; at0-cg1 / at1-cg0)."""
    wanted = [o.strip() for o in schemes if o.strip()]
    for s in wanted:
        sn.validate_scheme(s, BASELINES)
    out: list[dict] = []
    for scheme in SCHEMES:
        if scheme not in wanted:
            continue
        if not sn.is_greedymini(scheme):
            mod_t = xor_mod_t(k, m) if scheme == "mod-minimizer-xor" else None
            v = sn.make_id(scheme, mod_t=mod_t, baselines=BASELINES)
            v["gm_file"] = None
            out.append(v)
            continue
        cfg_tables, cfg_projs = gm_choices(scheme)
        for table in (tables if tables else cfg_tables):
            path = gm_table_file(table)
            gm_k = table_gm_k(path)
            for proj in (projections if projections else cfg_projs):
                v = sn.make_id(scheme, gm_k, proj, path.name, gm_k if scheme == "mod-minimizer-greedymini" else None,
                               baselines=BASELINES)
                v["gm_file"] = path
                out.append(v)
    return out


def check_gm_scheme(scheme: str, m: int, gm_file: Path) -> str | None:
    """Error message if the plain greedymini scheme cannot run at minimizer length m
    with the table `gm_file` (None if fine or not such a scheme): needs m - gm_k even
    and >= 0, and a mirror-symmetric table."""
    if scheme != "greedymini":
        return None
    from scripts.density import gm_table_mirror_symmetric, read_gm_ranks, _gm_k_from_table

    ranks = read_gm_ranks(gm_file)
    gm_k = _gm_k_from_table(ranks)
    if not sn.gm_centre_fits(m, gm_k):
        return (
            f"{scheme}: needs m - gm_k even and >= 0, but m={m} and the table "
            f"{gm_file.name} has gm_k={gm_k} (m - gm_k = {m - gm_k})"
        )
    if not gm_table_mirror_symmetric(ranks):
        return f"{scheme}: table {gm_file.name} is not mirror-symmetric (rank[key] != rank[reverse(key)])"
    return None


def mod_t(scheme: str, k: int, m: int, gm_file: Path | str | None = None) -> int:
    """The t-mer length of a mod-minimizer scheme at (k, m): xor -> 4 + (m - 4) mod w,
    greedymini -> the table's gm_k."""
    if scheme == "mod-minimizer-xor":
        return xor_mod_t(k, m)
    if scheme != "mod-minimizer-greedymini":
        raise ValueError(f"{scheme} is not a mod-minimizer scheme")
    path = Path(gm_file) if gm_file else gm_table_file(gm_choices(scheme)[0][0])
    return table_gm_k(path)


def scheme_dir_name(scheme: str, gm_k: int | None = None, projection: str | None = None,
                    mod_t_: int | None = None) -> str:
    """Run-directory segment of a variant (see scheme_names.scheme_dir_name)."""
    return sn.scheme_dir_name(scheme, gm_k, projection, mod_t_)


def build_flags(scheme: str, seed: int, gm_file: str | Path | None = None, k: int = 31,
                m: int = 21, projection: str | None = None) -> list[str]:
    """sshash `build` flags for `scheme` with `seed` (every scheme takes -s) at (k, m).
    GreedyMini schemes use `gm_file` (default: the first table config.yaml registers)
    and `projection` (default at0-cg1; only a non-default one is passed)."""
    sn.validate_scheme(scheme, BASELINES)
    flags: list[str] = ["-s", str(seed)]
    if scheme == "native":
        return flags
    if scheme in _SEED_ONLY:
        return ["--scheme", scheme, *flags]
    if scheme == "mod-minimizer-xor":
        return ["--scheme", scheme, "--mod-t", str(xor_mod_t(k, m)), *flags]
    path = Path(gm_file) if gm_file else gm_table_file(gm_choices(scheme)[0][0])
    if not path.is_file():
        raise FileNotFoundError(f"GreedyMini table not found: {path}")
    out = ["--scheme", scheme, "--greedy-mini-file", str(path)]
    if scheme == "mod-minimizer-greedymini":
        out += ["--mod-t", str(table_gm_k(path))]
    proj = sn.projection_from_flag(projection)
    if proj != sn.DEFAULT_PROJECTION:
        out += ["--gm-projection", proj]
    return [*out, *flags]


def _csv(s: str | None) -> list[str] | None:
    return [x.strip() for x in s.split(",") if x.strip()] if s else None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts schemes", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_flags = sub.add_parser("flags", help="print sshash build flags for a scheme/seed")
    p_flags.add_argument("scheme", choices=SCHEMES)
    p_flags.add_argument("seed", type=lambda s: int(s, 0))
    p_flags.add_argument("--gm-file", default=None)
    p_flags.add_argument("--projection", default=None)
    p_flags.add_argument("--k", type=int, default=31)
    p_flags.add_argument("--m", type=int, default=21)

    p_t = sub.add_parser("mod-t", help="print the t-mer length of a mod-minimizer scheme at (k, m)")
    p_t.add_argument("scheme")
    p_t.add_argument("k", type=int)
    p_t.add_argument("m", type=int)
    p_t.add_argument("--gm-file", default=None)

    p_dir = sub.add_parser("dir-name", help="print the run-directory segment of a variant")
    p_dir.add_argument("scheme")
    p_dir.add_argument("--gm-k", type=int, default=None)
    p_dir.add_argument("--projection", default=None)
    p_dir.add_argument("--mod-t", type=int, default=None)

    sub.add_parser("list", help="list canonical scheme names (plus the native baseline)")

    p_var = sub.add_parser(
        "variants",
        help="tab-separated run variants of a scheme list: scheme, gm_file, projection, gm_k, "
        "gm_table, mod_t ('-' = empty)")
    p_var.add_argument("schemes", help="comma-separated scheme names")
    p_var.add_argument("--k", type=int, required=True)
    p_var.add_argument("--m", type=int, required=True)
    p_var.add_argument("--gm-tables", default=None, help="comma-separated table names (override config.yaml)")
    p_var.add_argument("--gm-projections", default=None, help="comma-separated projections (override config.yaml)")

    p_chk = sub.add_parser("check-gm", help="validate a greedymini table for minimizer length m (exit 1 with a message if not)")
    p_chk.add_argument("scheme")
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
        print(" ".join(build_flags(args.scheme, args.seed, args.gm_file, args.k, args.m, args.projection)))
        return 0
    if args.cmd == "mod-t":
        if not sn.is_mod_minimizer(args.scheme):
            print(f"error: {args.scheme} is not a mod-minimizer scheme", file=sys.stderr)
            return 1
        print(mod_t(args.scheme, args.k, args.m, args.gm_file))
        return 0
    if args.cmd == "dir-name":
        print(scheme_dir_name(args.scheme, args.gm_k, args.projection, args.mod_t))
        return 0
    if args.cmd == "check-gm":
        msg = check_gm_scheme(args.scheme, args.m, Path(args.gm_file))
        if msg:
            print(f"error: {msg}", file=sys.stderr)
            return 1
        return 0
    if args.cmd == "variants":
        schemes = _csv(args.schemes) or []
        for v in expand_variants(schemes, args.k, args.m, _csv(args.gm_tables), _csv(args.gm_projections)):
            cells = [v["scheme"], v["gm_file"], v["projection"], v["gm_k"], v["gm_table"], v["mod_t"]]
            print("\t".join("-" if c is None else str(c) for c in cells))
        return 0
    if args.cmd == "list":
        for o in SCHEMES:
            print(o)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
