#!/usr/bin/env python3
"""Cross-validation of the greedymini order (scheme.cpp RankState<MM_SCHEME_GREEDYMINI>: the
shared GreedyMini centre extension, injective packing) against two
independent Python implementations:

- reference_value(): written per base from the definition -- for an l-mer
  x y z (|x| = |z| = i, |y| = gm_k), Sigma = 1 for C/G, Gamma = 1 for G/T,
  first base = most significant bit; value = rank(table[Sigma(y)]) << 2l |
  min over both strands of (Sigma << l | Gamma), the reverse complement
  reversing both planes and complementing Gamma; --gm-projection at1-cg0
  complements Sigma.
- benchmark/common/density.py greedymini_centre_hash (the shared bit-plane
  port); checked equal to the reference on every k-mer.

The compiled sketch (minimap2/sketch_dump driver, `make sketch_dump` in minimap2/) is compared with
- on N-free sequences, a brute-force minimizer scan over the reference
  values: per window of w k-mers, every position holding the window's
  minimum is selected (minimap2's robust winnowing keeps all tied minima);
- on sequences with N runs, a Python port of minimap2's own mm_sketch loop
  (sketch.c) fed the reference values. minimap2 does not flush the pending
  minimum at an N (l resets, min is overwritten), so a segment's last
  minimizer can be dropped -- upstream behaviour that scheme.cpp's
  sketch_single_window copies on purpose; the port is itself checked against
  the compiled native order (mm_hash64) on the same sequences.
Covers k = 15 (no extension), 17, 19, 21, 23 and both projections, plus
strand invariance of the value, and --order-seed (the Gamma toggle).
"""
from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SKETCH_DUMP = REPO / "minimap2" / "sketch_dump"
GM_FILE = REPO / "orders" / "GreedyMini" / "w10_k15_canonical.gm"
sys.path.insert(0, str(REPO / "benchmark" / "common"))
from density import greedymini_centre_hash, gm_k_from_table, read_gm_ranks  # noqa: E402

COMP = {0: 3, 1: 2, 2: 1, 3: 0}


def mm_hash64(key: int, mask: int) -> int:
    """Port of sketch.c's mm_hash64 (Thomas Wang's 64-bit integer mix), the
    native order; `(~key) & mask` gives C's fixed-width `~key` for any mask."""
    key = (~key + (key << 21)) & mask
    key = key ^ (key >> 24)
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ^ (key >> 14)
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ^ (key >> 28)
    key = (key + (key << 31)) & mask
    return key


def bits_to_int(bits: list[int]) -> int:
    v = 0
    for b in bits:
        v = (v << 1) | b
    return v


def reference_value(codes: list[int], gm_ranks: list[int], gm_k: int, flip: bool, toggle: int = 0) -> int:
    l = len(codes)
    toggle &= (1 << l) - 1
    i = (l - gm_k) // 2

    def planes(cs: list[int]) -> tuple[int, int]:
        sigma = [(1 if c in (1, 2) else 0) ^ (1 if flip else 0) for c in cs]
        gamma = [1 if c >= 2 else 0 for c in cs]
        return bits_to_int(sigma), bits_to_int(gamma)

    s, g = planes(codes)
    s_rc, g_rc = planes([COMP[c] for c in reversed(codes)])
    rank = gm_ranks[(s >> i) & ((1 << gm_k) - 1)]
    return (rank << (2 * l)) | min((s << l) | (g ^ toggle), (s_rc << l) | (g_rc ^ toggle))


def pack(codes: list[int]) -> int:
    v = 0
    for c in codes:
        v = (v << 2) | c
    return v


def brute_force_positions(seq: str, w: int, k: int, gm_ranks: list[int], gm_k: int, flip: bool,
                          toggle: int = 0) -> list[int]:
    """Sorted end positions (0-based raw index of a k-mer's last base)."""
    out: set[int] = set()
    start = 0
    for end in range(len(seq) + 1):
        if end < len(seq) and seq[end] != "N":
            continue
        seg = seq[start:end]
        codes = ["ACGT".index(c) for c in seg]
        vals = [reference_value(codes[p:p + k], gm_ranks, gm_k, flip, toggle) for p in range(len(codes) - k + 1)]
        for wst in range(len(vals) - w + 1):
            win = vals[wst:wst + w]
            m = min(win)
            out.update(start + wst + j + k - 1 for j, v in enumerate(win) if v == m)
        start = end + 1
    return sorted(out)


def mm_sketch_port(seq: str, w: int, k: int, value) -> list[int]:
    """Port of sketch.c mm_sketch (HPC off) with the per-k-mer value given by
    value(canonical_kmer_codes); returns sorted unique end positions."""
    INF = None
    shift1, mask = 2 * (k - 1), (1 << 2 * k) - 1
    kmer = [0, 0]
    buf = [(INF, -1)] * w  # (value, end position)
    mn, min_pos = (INF, -1), 0
    out: list[int] = []
    l = buf_pos = 0
    less_eq = lambda a, b: b is INF or (a is not INF and a <= b)  # noqa: E731
    for i, ch in enumerate(seq):
        info = (INF, -1)
        if ch != "N":
            c = "ACGT".index(ch)
            kmer[0] = ((kmer[0] << 2) | c) & mask
            kmer[1] = (kmer[1] >> 2) | ((3 ^ c) << shift1)
            if kmer[0] == kmer[1]:
                continue
            z = 0 if kmer[0] < kmer[1] else 1
            l += 1
            if l >= k:
                x = kmer[z]
                codes = [(x >> (2 * (k - 1 - p))) & 3 for p in range(k)]
                info = (value(codes), i)
        else:
            l = 0
        buf[buf_pos] = info
        if l == w + k - 1 and mn[0] is not INF:
            for j in list(range(buf_pos + 1, w)) + list(range(0, buf_pos)):
                if buf[j][0] == mn[0] and buf[j][1] != mn[1]:
                    out.append(buf[j][1])
        if less_eq(info[0], mn[0]):
            if l >= w + k and mn[0] is not INF:
                out.append(mn[1])
            mn, min_pos = info, buf_pos
        elif buf_pos == min_pos:
            if l >= w + k - 1 and mn[0] is not INF:
                out.append(mn[1])
            mn = (INF, -1)
            for j in list(range(buf_pos + 1, w)) + list(range(0, buf_pos + 1)):
                if less_eq(buf[j][0], mn[0]):
                    mn, min_pos = buf[j], j
            if l >= w + k - 1 and mn[0] is not INF:
                for j in list(range(buf_pos + 1, w)) + list(range(0, buf_pos + 1)):
                    if buf[j][0] == mn[0] and buf[j][1] != mn[1]:
                        out.append(buf[j][1])
        buf_pos = 0 if buf_pos + 1 == w else buf_pos + 1
    if mn[0] is not INF:
        out.append(mn[1])
    return sorted(set(out))


def run_binary(seq: str, w: int, k: int, flip: bool, scheme: str = "greedymini", seed: int = 0) -> list[int]:
    args = [str(SKETCH_DUMP), scheme, str(w), str(k), hex(seed), "0"]
    if scheme == "greedymini":
        args += [str(GM_FILE), "at1-cg0" if flip else "at0-cg1"]
    res = subprocess.run(args, input=seq, capture_output=True, text=True, check=True)
    return sorted(set(int(x) for x in res.stdout.split()))


def random_seq(rng: random.Random, n: int, with_n: bool) -> str:
    out = []
    while len(out) < n:
        if with_n and rng.random() < 0.002:
            out.extend("N" * rng.randint(1, 20))
        elif rng.random() < 0.05:  # homopolymer / low-complexity runs
            out.extend(rng.choice("ACGT") * rng.randint(2, 12))
        else:
            out.append("ACGT"[rng.randrange(4)])
    return "".join(out[:n])


def main() -> int:
    if not SKETCH_DUMP.is_file():
        print(f"error: {SKETCH_DUMP} not built (make sketch_dump)", file=sys.stderr)
        return 1
    gm_ranks = read_gm_ranks(GM_FILE)
    gm_k = gm_k_from_table(gm_ranks)
    rng = random.Random(20261007)
    failures = []
    n_cases = n_kmers = 0
    for k in (15, 17, 19, 21, 23):
        for flip in (False, True):
            # value: reference == shared port, and strand invariant
            for _ in range(3000):
                codes = [rng.randrange(4) for _ in range(k)]
                ref = reference_value(codes, gm_ranks, gm_k, flip)
                port = greedymini_centre_hash(pack(codes), k, gm_ranks, gm_k, flip)
                rc = reference_value([COMP[c] for c in reversed(codes)], gm_ranks, gm_k, flip)
                n_kmers += 1
                if ref != port or ref != rc:
                    failures.append(f"k={k} flip={flip} {codes}: ref={ref} port={port} rc={rc}")
                    break
            # compiled sketch vs brute force (N-free) / vs the mm_sketch port (N runs)
            gm_value = lambda codes, f=flip: reference_value(codes, gm_ranks, gm_k, f)  # noqa: E731
            for w, n, with_n in ((10, 3000, False), (5, 1500, False), (10, 3000, True), (5, 1500, True)):
                n_cases += 1
                seq = random_seq(rng, n, with_n)
                if with_n:
                    kmask = (1 << 2 * k) - 1
                    nat_exp = mm_sketch_port(seq, w, k, lambda codes: mm_hash64(pack(codes), kmask))
                    nat_got = run_binary(seq, w, k, flip, scheme="native")
                    if nat_exp != nat_got:
                        failures.append(f"k={k} w={w}: mm_sketch port != compiled native ({len(nat_exp)} vs {len(nat_got)})")
                    exp = mm_sketch_port(seq, w, k, gm_value)
                else:
                    exp = brute_force_positions(seq, w, k, gm_ranks, gm_k, flip)
                got = run_binary(seq, w, k, flip)
                if exp != got:
                    diff = sorted(set(exp) ^ set(got))[:5]
                    failures.append(f"k={k} w={w} flip={flip} N={with_n}: {len(exp)} expected vs {len(got)} got, first differing {diff}")
    # --order-seed: the Gamma toggle (low k bits of the seed), strand invariant.
    # It only reorders k-mers with equal rank and Sigma plane (rare on random
    # DNA, ~1% of the positions of a real genome at k = 19), so a change is
    # required somewhere, not in every case.
    n_changed = 0
    for k in (15, 19, 23):
        for seed in (0x9E3779B97F4A7C15, 0x2545F4914F6CDD1D):
            for _ in range(2000):
                codes = [rng.randrange(4) for _ in range(k)]
                ref = reference_value(codes, gm_ranks, gm_k, False, seed)
                rc = reference_value([COMP[c] for c in reversed(codes)], gm_ranks, gm_k, False, seed)
                port = greedymini_centre_hash(pack(codes), k, gm_ranks, gm_k, False, seed)
                n_kmers += 1
                if ref != rc or ref != port:
                    failures.append(f"k={k} seed={seed:#x} {codes}: ref={ref} rc={rc} port={port}")
                    break
            n_cases += 1
            seq = random_seq(rng, 3000, False)
            exp = brute_force_positions(seq, 10, k, gm_ranks, gm_k, False, seed)
            got = run_binary(seq, 10, k, False, seed=seed)
            if exp != got:
                failures.append(f"k={k} seed={seed:#x}: {len(exp)} expected vs {len(got)} got")
            n_changed += exp != brute_force_positions(seq, 10, k, gm_ranks, gm_k, False)
    if n_changed == 0:
        failures.append("--order-seed never changed the greedymini selection")
    print(f"checked {n_kmers} k-mer values (reference vs port vs reverse complement) and {n_cases} sketches")
    if failures:
        print(f"{len(failures)} FAILURE(S):", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print("all match")
    return 0


def test_greedymini_sketch():
    import pytest

    for f in (SKETCH_DUMP, GM_FILE):
        if not f.is_file():
            pytest.skip(f"{f} not built (make sketch_dump in minimap2/, the .gm tables in orders/GreedyMini/)")
    assert main() == 0


if __name__ == "__main__":
    raise SystemExit(main())
