"""Order name -> sshash CLI flags (see ../../../orders/README.md).

Every order (including ``native``) takes a seed with sshash -- unlike kmc,
where ``lex`` has no seed and is special-cased. So there is no "no seed"
branch here: every order in ORDERS is run once per seed the harness sweeps.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from scripts.config import default_gm_file

ORDERS = ("native", "xor-hash", "xor-pcg-hash", "greedymini")


def order_dir_name(order: str) -> str:
    """Filesystem-safe directory segment for an order (see run_rel_path)."""
    return {
        "native": "native",
        "xor-hash": "xor_hash",
        "xor-pcg-hash": "xor_pcg_hash",
        "greedymini": "greedymini",
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
    # greedymini
    gm_path = Path(gm_file) if gm_file else default_gm_file()
    if not gm_path.is_file():
        raise FileNotFoundError(f"GreedyMini table not found: {gm_path}")
    return ["-O", "greedymini", "--greedy-mini-file", str(gm_path), *flags]


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

    args = ap.parse_args(argv)
    if args.cmd == "flags":
        print(" ".join(build_flags(args.order, args.seed, args.gm_file)))
        return 0
    if args.cmd == "dir-name":
        print(order_dir_name(args.order))
        return 0
    if args.cmd == "list":
        for o in ORDERS:
            print(o)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
