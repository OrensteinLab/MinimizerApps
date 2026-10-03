"""Order name -> sshash CLI flags (see ../../../orders/README.md).

Every order (including ``native``) takes a seed with sshash -- unlike kmc,
where ``lex`` has no seed and is special-cased. So there is no "no seed"
branch here: every order in ORDERS is run once per seed the harness sweeps.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from scripts.config import greedymini_label_file, greedymini_labels

# Plain greedymini orders: "greedymini-<label>[-at1]" -- <label> picks the table from
# config.yaml's greedymini_labels (k13, k15), "-at1" the flipped DNA->binary polarity
# (sshash --gm-projection at1-cg0). Registering a table in config.yaml is enough to make a
# real, argparse-accepted order name here, matching lib/orders.sh. sshash extends the gm_k table to the m-mer by the
# centre extension, which needs m - gm_k even (see check_gm_order).
_AT1_SUFFIX = "-at1"


def _plain_greedymini_variants() -> tuple[str, ...]:
    bases = [f"greedymini-{label}" for label in greedymini_labels()]
    return (*bases, *(b + _AT1_SUFFIX for b in bases))


def parse_greedymini_order(order: str) -> tuple[str, str] | None:
    """(label, sshash --gm-projection value) for a plain greedymini order name
    greedymini-<label>[-at1]; None for any other order."""
    if not order.startswith("greedymini-"):
        return None
    proj = "at0-cg1"
    if order.endswith(_AT1_SUFFIX):
        proj = "at1-cg0"
        order = order[: -len(_AT1_SUFFIX)]
    return order[len("greedymini-"):], proj


def check_gm_order(order: str, m: int) -> str | None:
    """Error message if the plain greedymini order `order` cannot run at minimizer
    length m (None if fine or not such an order): needs m - gm_k even and >= 0, and a
    mirror-symmetric table."""
    plain = parse_greedymini_order(order)
    if plain is None:
        return None
    from scripts.density import gm_table_mirror_symmetric, read_gm_ranks, _gm_k_from_table

    path = greedymini_label_file(plain[0])
    ranks = read_gm_ranks(path)
    gm_k = _gm_k_from_table(ranks)
    if m < gm_k or (m - gm_k) % 2:
        return (
            f"{order}: needs m - gm_k even and >= 0, but m={m} and the table "
            f"{path.name} has gm_k={gm_k} (m - gm_k = {m - gm_k})"
        )
    if not gm_table_mirror_symmetric(ranks):
        return f"{order}: table {path.name} is not mirror-symmetric (rank[key] != rank[reverse(key)])"
    return None


ORDERS = (
    "native",
    "xor-hash",
    "xor-pcg-hash",
    *_plain_greedymini_variants(),
)


def order_dir_name(order: str) -> str:
    """Filesystem-safe directory segment for an order (see run_rel_path)."""
    return {
        "native": "native",
        "xor-hash": "xor_hash",
        "xor-pcg-hash": "xor_pcg_hash",
    }.get(order, order.replace("-", "_"))


def build_flags(order: str, seed: int, gm_file: str | None = None) -> list[str]:
    """sshash `build` flags for `order` with `seed` (every order takes -s)."""
    if order not in ORDERS:
        raise ValueError(f"unknown order: {order} (use one of {', '.join(ORDERS)})")
    flags: list[str] = ["-s", str(seed)]
    if order == "native":
        return flags
    if order in ("xor-hash", "xor-pcg-hash"):
        return ["-O", order, *flags]
    label, proj = parse_greedymini_order(order)  # type: ignore[misc]
    gm_path = Path(gm_file) if gm_file else greedymini_label_file(label)
    if not gm_path.is_file():
        raise FileNotFoundError(f"GreedyMini table not found: {gm_path}")
    out = ["-O", "greedymini", "--greedy-mini-file", str(gm_path)]
    if proj != "at0-cg1":
        out += ["--gm-projection", proj]
    return [*out, *flags]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts orders", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_flags = sub.add_parser("flags", help="print sshash build flags for an order/seed")
    p_flags.add_argument("order", choices=ORDERS)
    p_flags.add_argument("seed", type=lambda s: int(s, 0))
    p_flags.add_argument("--gm-file", default=None)

    p_dir = sub.add_parser("dir-name", help="print the run-directory segment for an order")
    p_dir.add_argument("order")

    sub.add_parser("list", help="list known order names")

    p_chk = sub.add_parser("check-gm", help="validate a greedymini-<label> order for minimizer length m (exit 1 with a message if not)")
    p_chk.add_argument("order")
    p_chk.add_argument("m", type=int)

    args = ap.parse_args(argv)
    if args.cmd == "flags":
        print(" ".join(build_flags(args.order, args.seed, args.gm_file)))
        return 0
    if args.cmd == "dir-name":
        print(order_dir_name(args.order))
        return 0
    if args.cmd == "check-gm":
        msg = check_gm_order(args.order, args.m)
        if msg:
            print(f"error: {msg}", file=sys.stderr)
            return 1
        return 0
    if args.cmd == "list":
        for o in ORDERS:
            print(o)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
