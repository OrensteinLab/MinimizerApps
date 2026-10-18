#!/usr/bin/env python3
"""Brute-force cross-validation of sketch_mod_minimizer() (order.cpp), the
strand-symmetric mod-minimizer-xor sketch.

brute_force() applies the definition window by window (no deques, no
incremental state): for the window of L = w+k-1 bases starting at k-mer
ordinal `win`, the T = w+k-t t-mers are ranked as canonical t-mer XOR seed;
the preferred strand has more G/T bases (tie, even L: window <= its reverse
complement lexicographically); a forward-preferred window samples the k-mer at
offset x_L mod w (x_L = leftmost minimal t-mer, as in Groot Koerkamp & Pibiri,
WABI 2024), a reverse-preferred one at (w-1) - (((T-1) - x_R) mod w) (x_R =
rightmost minimal t-mer).

On many random (w, k, t, seed) configs the compiled code, via the sketch_dump
test driver (`make sketch_dump`), must emit exactly the brute force's loci,
including sequences with N runs (every N resets the scan) and even k
(palindromic k-mers are skipped, as in native mm_sketch()); and a sequence and
its reverse complement must select mirrored loci, except inside windows that
are their own reverse complement (no preferred strand).
"""
from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

MASK64 = (1 << 64) - 1
SKETCH_DUMP = Path(__file__).resolve().parent / "sketch_dump"


def xor(value: int, seed: int) -> int:
    return (value ^ seed) & MASK64


def revcomp(x: int, n: int) -> int:
    r = 0
    for _ in range(n):
        r = (r << 2) | (3 - (x & 3))
        x >>= 2
    return r


def canon(x: int, n: int) -> int:
    r = revcomp(x, n)
    return x if x <= r else r


def rolling_values(seq: list[int], length: int) -> list[int]:
    if len(seq) < length:
        return []
    mask = (1 << (2 * length)) - 1
    state = 0
    out = []
    for i, c in enumerate(seq):
        state = ((state << 2) | c) & mask
        if i + 1 >= length:
            out.append(state)
    return out


def brute_force(seq: list[int], w: int, k: int, t: int, seed: int) -> list[int]:
    """Sorted distinct k-mer ordinals selected by the definition above."""
    T = w + k - t
    L = w + k - 1
    t_raw = rolling_values(seq, t)
    out = set()
    for win in range(len(seq) - L + 1):
        bases = seq[win:win + L]
        gt = sum(1 for b in bases if b >= 2)
        if 2 * gt != L:
            fwd = 2 * gt > L
        else:
            fwd = bases <= [3 - b for b in reversed(bases)]
        ranks = [xor(canon(t_raw[win + x], t), seed) for x in range(T)]
        m = min(ranks)
        if fwd:
            off = ranks.index(m) % w
        else:
            x_r = T - 1 - ranks[::-1].index(m)
            off = (w - 1) - (((T - 1) - x_r) % w)
        out.add(win + off)
    return sorted(out)


def expected_end_positions(segments: list[tuple[int, list[int]]], w: int, k: int, t: int, seed: int) -> list[int]:
    """Brute force per N-free segment, as raw k-mer end positions; drops
    palindromic k-mers, which the C++ (like native) never emits."""
    out = []
    for start, seg in segments:
        k_raw = rolling_values(seg, k)
        for o in brute_force(seg, w, k, t, seed):
            if k_raw[o] != revcomp(k_raw[o], k):
                out.append(start + o + k - 1)
    return out


def make_sequence(rng: random.Random, seq_len: int, with_n: bool) -> tuple[str, list[tuple[int, list[int]]]]:
    """A random ACGT string (optionally with a few N runs) and its N-free
    segments as (raw start, 2-bit codes)."""
    chars = []
    i = 0
    while i < seq_len:
        if with_n and rng.random() < 0.002:
            n = rng.randint(1, 30)
            chars.extend("N" * n)
            i += n
        else:
            chars.append("ACGT"[rng.randrange(4)])
            i += 1
    s = "".join(chars)
    segments, cur, cur_start = [], [], 0
    for j, ch in enumerate(s):
        if ch == "N":
            if cur:
                segments.append((cur_start, cur))
            cur, cur_start = [], j + 1
        else:
            cur.append("ACGT".index(ch))
    if cur:
        segments.append((cur_start, cur))
    return s, segments


def run_binary(s: str, w: int, k: int, t: int, seed: int) -> list[int]:
    res = subprocess.run(
        [str(SKETCH_DUMP), "mod-minimizer-xor", str(w), str(k), hex(seed), str(t)],
        input=s, capture_output=True, text=True, check=True,
    )
    return [int(x) for x in res.stdout.split()]


def run_case(rng: random.Random, w: int, k: int, t: int, seed: int, seq_len: int) -> str | None:
    plain = make_sequence(rng, seq_len, with_n=False)
    with_n_seq = make_sequence(rng, seq_len, with_n=True)
    for with_n, (s_, segs) in ((False, plain), (True, with_n_seq)):
        exp = expected_end_positions(segs, w, k, t, seed)
        got = run_binary(s_, w, k, t, seed)
        if exp != got:
            i = next((i for i, (a, b) in enumerate(zip(exp, got)) if a != b), min(len(exp), len(got)))
            return (f"compiled order.cpp vs brute force differ at emission {i} (with_n={with_n}, "
                    f"exp_len={len(exp)} got_len={len(got)})")
    # a sequence and its reverse complement select mirrored loci
    n = len(plain[0])
    rc_s = plain[0][::-1].translate(str.maketrans("ACGT", "TGCA"))
    fwd = set(run_binary(plain[0], w, k, t, seed))
    mirrored = {(n - 1) - (e - k + 1) for e in run_binary(rc_s, w, k, t, seed)}
    # The only allowed exception: a window that is its own reverse complement
    # (possible only for even L) has no preferred strand, so both strands use
    # the forward rule there.
    L = w + k - 1
    codes = plain[1][0][1]
    palin_kmers = set()
    for st in range(len(codes) - L + 1):
        win = codes[st:st + L]
        if win == [3 - b for b in reversed(win)]:
            palin_kmers.update(range(st, st + w))
    unexplained = {e for e in fwd ^ mirrored if (e - k + 1) not in palin_kmers}
    if unexplained:
        return f"{len(unexplained)} loci differ between the strands outside self-reverse-complement windows"
    return None


def main() -> int:
    configs = [
        # (w, k, t, seq_len)
        (10, 15, 5, 4000),   # the benchmark's own k15_w10, t=5
        (10, 19, 9, 20000),  # the benchmark's k19_w10, t=9
        (50, 15, 5, 6000),   # large window, L > 32 (ring path)
        (11, 15, 5, 3000),
        (10, 21, 8, 4000),
        (5, 12, 4, 3000),    # even k: palindromes possible
        (16, 20, 7, 5000),   # L = 35 > 32, odd
        (3, 8, 2, 2000),     # even L: preference ties
        (2, 6, 2, 1500),
        (1, 6, 3, 1200),     # w=1
        (100, 15, 5, 8000),
        (10, 15, 14, 3000),  # t close to k
    ]
    seeds = [0x2545F4914F6CDD1D, 0x9E3779B97F4A7C15, 0x1234567890ABCDEF, 0, 1]
    if not SKETCH_DUMP.is_file():
        print(f"error: {SKETCH_DUMP.name} not built (make sketch_dump)", file=sys.stderr)
        return 1

    n_cases = 0
    failures = []
    rng = random.Random(1234567)
    for (w, k, t, seq_len) in configs:
        for seed in seeds:
            n_cases += 1
            msg = run_case(rng, w, k, t, seed, seq_len)
            if msg:
                failures.append(f"(w={w},k={k},t={t},seed={hex(seed)}): {msg}")

    print(f"ran {n_cases} (w,k,t,seed) cases against the compiled sketch")
    if failures:
        print(f"{len(failures)} FAILURE(S):", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print("all cases match the brute-force reference; selections are strand-symmetric")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
