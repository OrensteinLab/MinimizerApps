#!/usr/bin/env python3
"""Shared scheme order, reference scheme, and toggle seeds for the benchmark harness.

Usage:
  python3 benchmark/kraken2/bench_const.py seeds
"""

from __future__ import print_function

import sys

REF_SCHEME = "xor"

SCHEME_ORDER = {
    "xor": 0,
    "xor-pcg": 1,
    "greedymini": 2,
    "mod-minimizer-xor": 3,
    "mod-minimizer-greedymini": 4,
    "greedymini-at1": 5,
    "mod-minimizer-greedymini-at1": 6,
}

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
    return (SCHEME_ORDER.get(name, 99), name)


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
    sys.stderr.write("usage: bench_const.py seeds|hex [N]\n")
    return 2


if __name__ == "__main__":
    sys.exit(main())
