"""Order name -> sshash CLI flags (see ../../../orders/README.md).

Every order (including ``native``) takes a seed with sshash -- unlike kmc,
where ``lex`` has no seed and is special-cased. So there is no "no seed"
branch here: every order in ORDERS is run once per seed the harness sweeps.

The 4 ``mod-minimizer-*`` orders (Groot Koerkamp & Pibiri, WABI 2024) reuse
the 4 base orders as O_t (the order on t-mers), plus a ``t`` length: sshash's
own ``-O`` flag accepts all 8 names directly and decomposes them internally
(see sshash/include/custom_order_hasher.hpp::parse_order_name), so this file
just needs to also thread ``--t-length`` through for them, and sweep ``t``
itself (see ``order_t_values``).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from scripts.config import default_gm_file, mod_t_gm_file

ORDERS = (
    "native",
    "xor-hash",
    "xor-pcg-hash",
    "greedymini",
    "mod-minimizer-native-hash",
    "mod-minimizer-xor-hash",
    "mod-minimizer-pcg-xor-hash",
    "mod-minimizer-greedymini",
)

MOD_ORDERS = frozenset(o for o in ORDERS if o.startswith("mod-minimizer-"))

# t values swept for the 3 hash-based mod-minimizer orders: t=10 reproduces
# the mod-minimizer paper's own published Human/GRCh38 numbers; t=11 matches
# the GreedyMini t-order table's own length, for a fair comparison against
# mod-minimizer-greedymini (which is fixed at t=11 -- see order_t_values).
MOD_T_VALUES = (10, 11)


def order_dir_name(order: str) -> str:
    """Filesystem-safe directory segment for an order (see run_rel_path)."""
    return {
        "native": "native",
        "xor-hash": "xor_hash",
        "xor-pcg-hash": "xor_pcg_hash",
        "greedymini": "greedymini",
    }.get(order, order.replace("-", "_"))


def order_t_values(order: str) -> tuple[int | None, ...]:
    """t values to sweep for `order`: (None,) for the 4 base orders, (11,)
    for mod-minimizer-greedymini (fixed, matching w21_k11_canonical.gm's own
    gm_k), MOD_T_VALUES for the 3 hash-based mod-minimizer orders."""
    if order not in MOD_ORDERS:
        return (None,)
    if order == "mod-minimizer-greedymini":
        return (11,)
    return MOD_T_VALUES


def build_flags(order: str, seed: int, gm_file: str | None = None, t: int | None = None) -> list[str]:
    """sshash `build` flags for `order` with `seed` (every order takes -s)."""
    if order not in ORDERS:
        raise ValueError(f"unknown order: {order} (use one of {', '.join(ORDERS)})")
    flags: list[str] = ["-s", str(seed)]
    if order == "native":
        return flags
    if order in ("xor-hash", "xor-pcg-hash"):
        return ["-O", order, *flags]
    if order in ("mod-minimizer-native-hash", "mod-minimizer-xor-hash", "mod-minimizer-pcg-xor-hash"):
        if t is None:
            raise ValueError(f"order {order} requires t (see order_t_values)")
        return ["-O", order, "--t-length", str(t), *flags]
    # greedymini / mod-minimizer-greedymini: both need --greedy-mini-file
    if order == "mod-minimizer-greedymini" and t is None:
        raise ValueError(f"order {order} requires t (see order_t_values)")
    gm_path = Path(gm_file) if gm_file else (mod_t_gm_file() if order == "mod-minimizer-greedymini" else default_gm_file())
    if not gm_path.is_file():
        raise FileNotFoundError(f"GreedyMini table not found: {gm_path}")
    out = ["-O", order, "--greedy-mini-file", str(gm_path)]
    if t is not None:
        out += ["--t-length", str(t)]
    return [*out, *flags]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts orders", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_flags = sub.add_parser("flags", help="print sshash build flags for an order/seed")
    p_flags.add_argument("order", choices=ORDERS)
    p_flags.add_argument("seed", type=lambda s: int(s, 0))
    p_flags.add_argument("--gm-file", default=None)
    p_flags.add_argument("--t", type=int, default=None)

    p_dir = sub.add_parser("dir-name", help="print the run-directory segment for an order")
    p_dir.add_argument("order")

    sub.add_parser("list", help="list known order names")

    args = ap.parse_args(argv)
    if args.cmd == "flags":
        print(" ".join(build_flags(args.order, args.seed, args.gm_file, args.t)))
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
