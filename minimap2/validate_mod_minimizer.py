#!/usr/bin/env python3
"""Brute-force cross-validation of sketch_mod_minimizer() (order.cpp), the
mod-minimizer-xor sketch.

brute_force() scans every t-mer locus of each window directly (no deques, no
incremental state) -- the textbook definition (Groot Koerkamp & Pibiri, WABI
2024): for the k-mer window starting at k-mer ordinal `win`, the t-mer loci
are the tmer_span = w+k-t positions [win, win+tmer_span-1]; pick the smallest
rank (ties -> smallest position x); the sampled k-mer is at win + (x mod w).

Two things are checked against it, on many random (w, k, t, seed) configs:
- incremental(): a Python port of order.cpp's algorithm (one monotone deque
  over all live t-mers, front = leftmost minimum, residues tracked by
  wrap-around counters), kept structurally parallel to the C++ so a mismatch
  points at the same line there;
- the compiled code itself, via the sketch_dump test driver (`make
  sketch_dump`), including sequences with N runs (every N resets the scan) and
  even k (palindromic k-mers are skipped, as in native mm_sketch()). Skipped
  with a note if the driver isn't built.
"""
from __future__ import annotations

import random
import subprocess
import sys
from collections import deque
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
    """Emitted k-mer ordinals (deduped across consecutive windows like the
    C++'s last_emitted check)."""
    tmer_span = w + k - t
    t_raw = rolling_values(seq, t)
    k_raw = rolling_values(seq, k)
    out = []
    last_emitted = -1
    for win in range(len(k_raw) - w + 1):
        best_rank, best_pos = None, None
        for pos in range(win, win + tmer_span):
            r = xor(canon(t_raw[pos], t), seed)
            if best_rank is None or r < best_rank:  # strict: ties keep the smallest pos
                best_rank, best_pos = r, pos
        abs_ord = win + (best_pos - win) % w
        if abs_ord != last_emitted:
            out.append(abs_ord)
            last_emitted = abs_ord
    return out


def incremental(seq: list[int], w: int, k: int, t: int, seed: int) -> list[int]:
    """Port of order.cpp's sketch_mod_minimizer selection."""
    t_raw = rolling_values(seq, t)
    dq: deque[tuple[int, int, int]] = deque()  # (rank, pos, cls), monotone non-decreasing rank
    t_cls = l_cls = 0
    out = []
    last_emitted = -1
    for run_pos in range(len(seq)):
        ot = run_pos - t + 1
        if ot >= 0:
            rank = xor(canon(t_raw[ot], t), seed)
            while dq and dq[-1][0] > rank:
                dq.pop()
            dq.append((rank, ot, t_cls))
            t_cls = 0 if t_cls + 1 == w else t_cls + 1
        ol = run_pos - k + 1
        if ol < 0:
            continue
        slot = l_cls
        l_cls = 0 if l_cls + 1 == w else l_cls + 1
        win = ol - w + 1
        if win < 0:
            continue
        wm = 0 if slot + 1 == w else slot + 1
        while dq[0][1] < win:
            dq.popleft()
        _, _, cls = dq[0]
        best_p = cls - wm
        if best_p < 0:
            best_p += w
        abs_ord = win + best_p
        if abs_ord != last_emitted:
            out.append(abs_ord)
            last_emitted = abs_ord
    return out


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


def run_case(rng: random.Random, w: int, k: int, t: int, seed: int, seq_len: int, use_binary: bool) -> str | None:
    s, segments = make_sequence(rng, seq_len, with_n=False)
    seq = segments[0][1]
    bf = brute_force(seq, w, k, t, seed)
    inc = incremental(seq, w, k, t, seed)
    if bf != inc:
        i = next((i for i, (a, b) in enumerate(zip(bf, inc)) if a != b), min(len(bf), len(inc)))
        return f"port vs brute force differ at emission {i} (bf_len={len(bf)} port_len={len(inc)})"
    if use_binary:
        for with_n in (False, True):
            if with_n:
                s, segments = make_sequence(rng, seq_len, with_n=True)
            exp = expected_end_positions(segments, w, k, t, seed)
            got = run_binary(s, w, k, t, seed)
            if exp != got:
                i = next((i for i, (a, b) in enumerate(zip(exp, got)) if a != b), min(len(exp), len(got)))
                return f"compiled order.cpp vs brute force differ at emission {i} (with_n={with_n}, exp_len={len(exp)} got_len={len(got)})"
    return None


def main() -> int:
    configs = [
        # (w, k, t, seq_len)
        (10, 15, 5, 4000),   # the benchmark's own w10_k15, t=5
        (10, 15, 5, 20000),
        (50, 15, 5, 6000),   # large window
        (11, 15, 5, 3000),
        (10, 21, 8, 4000),
        (5, 12, 4, 3000),    # even k: palindromes possible
        (16, 20, 7, 5000),
        (3, 8, 2, 2000),
        (2, 6, 2, 1500),
        (1, 6, 3, 1200),     # w=1
        (100, 15, 5, 8000),
        (10, 15, 14, 3000),  # t close to k
    ]
    seeds = [0x2545F4914F6CDD1D, 0x9E3779B97F4A7C15, 0x1234567890ABCDEF, 0, 1]
    use_binary = SKETCH_DUMP.is_file()
    if not use_binary:
        print(f"note: {SKETCH_DUMP.name} not built (make sketch_dump) -- checking the Python port only")

    n_cases = 0
    failures = []
    rng = random.Random(1234567)
    for (w, k, t, seq_len) in configs:
        for seed in seeds:
            n_cases += 1
            msg = run_case(rng, w, k, t, seed, seq_len, use_binary)
            if msg:
                failures.append(f"(w={w},k={k},t={t},seed={hex(seed)}): {msg}")

    print(f"ran {n_cases} (w,k,t,seed) cases" + (" (port + compiled binary)" if use_binary else " (port only)"))
    if failures:
        print(f"{len(failures)} FAILURE(S):", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print("all cases match the brute-force reference")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
