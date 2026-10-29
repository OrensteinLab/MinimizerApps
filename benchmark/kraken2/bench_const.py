#!/usr/bin/env python3
"""Shared scheme display sequence, reference scheme, and toggle seeds for the benchmark harness.

Usage:
  python3 benchmark/kraken2/bench_const.py seeds
  python3 benchmark/kraken2/bench_const.py scheme-id SCHEME GM_FILE PROJECTION L K
      -> '|'-separated: run_dir  gm_k  projection  gm_table  mod_t
         (run_dir from scheme_names.scheme_dir_name; fields blank when n/a)

A scheme is one of the canonical names (benchmark/common/scheme_names.py: xor,
xor-pcg-hash, greedymini, mod-minimizer-xor, mod-minimizer-greedymini); the
GreedyMini parameters (gm_k, projection, gm_table) and the mod-minimizer t are
separate fields, never part of the scheme name.
"""

from __future__ import print_function

import os
import sys

_COMMON = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common")
if _COMMON not in sys.path:
    sys.path.insert(0, _COMMON)

import scheme_names  # noqa: E402

REF_SCHEME = "xor"

SCHEME_ORDER = dict((o, i) for i, o in enumerate(scheme_names.CANONICAL))

# Kraken 2 defaults of the benchmark (k2 build: l=31, k=35 -> w=5).
DEFAULT_L = 31
DEFAULT_K = 35
MOD_R = 4  # mod-minimizer t = r + ((l - r) % w), r = 4 (src/mod_minimizer.h)


def mod_t(l=DEFAULT_L, k=DEFAULT_K):
    return MOD_R + ((l - MOD_R) % (k - l + 1))


def scheme_rank(name):
    """Sort rank of a scheme label: the canonical rank of the scheme name, optionally followed by
    ' (params)' or ' [cell]' ('greedymini (k15_at1-cg0)')."""
    return SCHEME_ORDER.get(str(name).split(" ")[0], 99)


def scheme_id(scheme, gm_file=None, projection=None, l=DEFAULT_L, k=DEFAULT_K):
    """scheme_names.make_id() for one benchmark scheme (gm_file: the .gm table
    the scheme uses; mod_t derived from l / k for the mod-minimizer schemes)."""
    table = os.path.basename(gm_file) if gm_file else None
    t = mod_t(l, k) if scheme_names.is_mod_minimizer(scheme) else None
    return scheme_names.make_id(scheme, projection=projection, gm_table=table, mod_t=t)


def scheme_run_dir(idd):
    return scheme_names.scheme_dir_name(idd["scheme"], idd["gm_k"], idd["projection"], idd["mod_t"])

# Same 64-bit binary strings as historical run_benchmark.sh DEFAULT_SEEDS.
DEFAULT_SEEDS = (
    "1110001101111110001010001100010000100111000110110101101000101101",
    "0101101000101101000111001000101100111111100111100100101001110000",
    "1001111100111100011011100001101000101011100001000111010111010000",
    "0001101101111010100111000100111001011101011000101000001111110001",
    "0111110101001110001010101001000111001000101100111111000001010110",
    "0010110100011100100010110011111110011110010010100111000011100011",
    "0110111000011010001010111000010001110101110100001001111100111100",
    "1001110001001110010111010110001010000011111100010001101101111010",
    "1000101100111111100111100100101001110000111000110010110100011100",
    "0111010111010000100111110011110001101110000110100010101110000100",
    "0101110101100010100000111111000100011011011110101001110001001110",
    "0011111110011110010010100111000011100011001011010001110010001011",
)

DENSITY_METRICS = (
    "dens_expected",
    "dens_particular",
    "fac_expected",
    "fac_particular",
)

DENSITY_CHART_METRICS = (
    ("fac_expected", "fac_expected"),
    ("fac_particular", "fac_particular"),
)


def seed_hex(bits):
    return "0x%x" % int(bits, 2)


def scheme_sort_key(name):
    return (scheme_rank(name), name)


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if argv and argv[0] == "seeds":
        for s in DEFAULT_SEEDS:
            print(s)
        return 0
    if argv and argv[0] == "hex":
        n = int(argv[1]) if len(argv) > 1 else len(DEFAULT_SEEDS)
        for bits in DEFAULT_SEEDS[:n]:
            print(seed_hex(bits))
        return 0
    if argv and argv[0] == "scheme-id" and len(argv) == 6:
        scheme, gm_file, projection, l, k = argv[1:]
        idd = scheme_id(scheme, gm_file or None, projection or None, int(l), int(k))
        f = scheme_names.format_id(idd)
        print("|".join([scheme_run_dir(idd), f["gm_k"], f["projection"],
                        f["gm_table"], f["mod_t"]]))
        return 0
    sys.stderr.write("usage: bench_const.py seeds|hex [N]|scheme-id SCHEME GM_FILE PROJECTION L K\n")
    return 2


if __name__ == "__main__":
    sys.exit(main())
