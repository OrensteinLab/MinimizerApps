#!/usr/bin/env python3
"""Brute-force cross-validation for sketch_mod_minimizer()'s O(1)-amortized
selection (order.cpp), the kraken2-style expire_a/front_rank[] optimization
that replaced a full O(w) rescan of all residue-class deques on every base
position. Same methodology used to find/fix this function's two earlier bugs
(the tie-break and deque-capacity bugs -- see order.cpp's own comments and
git history): an independent brute-force reference, cross-checked against a
faithful port of the actual algorithm, on many random (w,k,t,seed) configs.

brute_force() scans every t-mer locus in each window directly (no deques, no
caching, no incremental state at all) -- the textbook definition (Groot
Koerkamp & Pibiri, WABI 2024, Algorithm 6): for the k-mer window starting at
l-mer ordinal `win`, the relevant t-mer loci are the tmer_span = w+k-t
absolute positions [win, win+tmer_span-1]; pick the one with the smallest
rank (ties broken by smallest position); the sampled k-mer's window-relative
offset is (winning_pos - win) mod w.

incremental() is a direct Python port of order.cpp's new algorithm (monotone
per-residue deques + a persistent, O(1)-amortized-updated GLOBAL best cache),
kept structurally parallel to the C++ so a mismatch here points at the same
line in order.cpp.
"""
from __future__ import annotations

import random
import sys
from collections import deque

MASK64 = (1 << 64) - 1


def xor_hash(value: int, seed: int) -> int:
    return (value ^ seed) & MASK64


def pcg_rxs_m_xs64(x: int) -> int:
    x &= MASK64
    count = x >> 59
    x ^= x >> (5 + count)
    x = (x * 12605985483714917081) & MASK64
    return x ^ (x >> 43)


def xor_pcg_hash(value: int, seed: int) -> int:
    return pcg_rxs_m_xs64((value ^ seed) & MASK64)


def revcomp(x: int, n: int) -> int:
    r = 0
    for i in range(n):
        r = (r << 2) | (3 - (x & 3))
        x >>= 2
    return r


def canon(x: int, n: int) -> int:
    r = revcomp(x, n)
    return x if x <= r else r


def rolling_values(seq: list[int], length: int) -> list[int]:
    n = len(seq)
    if n < length:
        return []
    mask = (1 << (2 * length)) - 1
    state = 0
    for i in range(length):
        state = (state << 2) | seq[i]
    out = [state]
    for i in range(length, n):
        state = ((state << 2) | seq[i]) & mask
        out.append(state)
    return out


def brute_force(seq: list[int], w: int, k: int, t: int, seed: int, hash_fn) -> list[tuple[int, int]]:
    """Returns [(abs_kmer_ordinal, kmer_value), ...] emitted, deduped like
    the C++'s own last_emitted check."""
    n = len(seq)
    tmer_span = w + k - t
    t_raw = rolling_values(seq, t)
    k_raw = rolling_values(seq, k)
    n_wins = len(k_raw) - w + 1
    if n_wins <= 0 or not t_raw:
        return []
    out = []
    last_emitted = -1
    for win in range(n_wins):
        best_rank, best_pos = None, None
        for pos in range(win, win + tmer_span):
            if pos >= len(t_raw):
                break
            r = hash_fn(canon(t_raw[pos], t), seed)
            if best_rank is None or r < best_rank or (r == best_rank and pos < best_pos):
                best_rank, best_pos = r, pos
        if best_rank is None:
            continue
        best_p = (best_pos - win) % w
        abs_ord = win + best_p
        if abs_ord != last_emitted:
            out.append((abs_ord, k_raw[abs_ord]))
            last_emitted = abs_ord
    return out


def incremental(seq: list[int], w: int, k: int, t: int, seed: int, hash_fn) -> list[tuple[int, int]]:
    """Direct port of order.cpp's new sketch_mod_minimizer selection logic:
    a persistent GLOBAL best (best_a/best_rank/best_pos), updated in true
    O(1) amortized time -- not a per-residue front cache (see order.cpp's
    own comment for why that first attempt barely helped in practice)."""
    tmer_span = w + k - t
    t_raw = rolling_values(seq, t)
    k_raw = rolling_values(seq, k)
    if not t_raw or len(k_raw) < w:
        return []

    deques = [deque() for _ in range(w)]  # each holds (rank, pos), monotone increasing rank

    def push_monotone(a: int, rank: int, pos: int) -> None:
        dq = deques[a]
        while dq and dq[-1][0] > rank:
            dq.pop()
        dq.append((rank, pos))

    def pop_stale_and_peek(a: int, win: int):
        dq = deques[a]
        while dq and dq[0][1] < win:
            dq.popleft()
        return dq[0] if dq else None

    out = []
    last_emitted = -1
    win_start = -1
    lmers_stored = 0
    best_a, best_rank, best_pos = -1, None, None

    for ot in range(len(t_raw)):
        r = hash_fn(canon(t_raw[ot], t), seed)
        push_monotone(ot % w, r, ot)

        ol = ot - (k - t)  # k-mer ordinal whose t-mer window just became complete at this ot
        if ol < 0 or ol >= len(k_raw):
            continue
        if lmers_stored < w:
            lmers_stored += 1
        win = ol - w + 1
        if win < 0 or lmers_stored < w or ot < win + tmer_span - 1:
            continue

        if win_start >= 0 and win == win_start + 1 and best_a >= 0 and best_pos >= win:
            front = pop_stale_and_peek(ot % w, win)
            if front is not None:
                r2, p2 = front
                if r2 < best_rank or (r2 == best_rank and p2 < best_pos):
                    best_a, best_rank, best_pos = ot % w, r2, p2
        else:
            best_a, best_rank, best_pos = -1, None, None
            for a in range(w):
                front = pop_stale_and_peek(a, win)
                if front is None:
                    continue
                r2, p2 = front
                if best_a < 0 or r2 < best_rank or (r2 == best_rank and p2 < best_pos):
                    best_a, best_rank, best_pos = a, r2, p2

        if best_a >= 0:
            wm = win % w
            best_p = ((best_a - wm) % w + w) % w
            abs_ord = win + best_p
            if abs_ord != last_emitted and abs_ord < len(k_raw):
                out.append((abs_ord, k_raw[abs_ord]))
                last_emitted = abs_ord
            win_start = win
    return out


def naive_incremental(seq: list[int], w: int, k: int, t: int, seed: int, hash_fn) -> list[tuple[int, int]]:
    """Deque-based, but rescans+pops ALL w deques every window (no expire_a
    shortcut) -- matches order.cpp's ORIGINAL, already-proven-correct logic
    (pre-optimization). Used only to isolate whether a mismatch is in the
    expire_a caching specifically vs. the shared deque/window machinery."""
    tmer_span = w + k - t
    t_raw = rolling_values(seq, t)
    k_raw = rolling_values(seq, k)
    if not t_raw or len(k_raw) < w:
        return []
    deques = [deque() for _ in range(w)]

    def push_monotone(a: int, rank: int, pos: int) -> None:
        dq = deques[a]
        while dq and dq[-1][0] > rank:
            dq.pop()
        dq.append((rank, pos))

    out = []
    last_emitted = -1
    lmers_stored = 0
    for ot in range(len(t_raw)):
        r = hash_fn(canon(t_raw[ot], t), seed)
        push_monotone(ot % w, r, ot)
        ol = ot - (k - t)
        if ol < 0 or ol >= len(k_raw):
            continue
        if lmers_stored < w:
            lmers_stored += 1
        win = ol - w + 1
        if win < 0 or lmers_stored < w or ot < win + tmer_span - 1:
            continue
        best_a, best_rank, best_pos = -1, None, None
        for a in range(w):
            dq = deques[a]
            while dq and dq[0][1] < win:
                dq.popleft()
            if not dq:
                continue
            r2, p2 = dq[0]
            if best_a < 0 or r2 < best_rank or (r2 == best_rank and p2 < best_pos):
                best_a, best_rank, best_pos = a, r2, p2
        if best_a >= 0:
            wm = win % w
            best_p = ((best_a - wm) % w + w) % w
            abs_ord = win + best_p
            if abs_ord != last_emitted and abs_ord < len(k_raw):
                out.append((abs_ord, k_raw[abs_ord]))
                last_emitted = abs_ord
    return out


def run_case(rng: random.Random, w: int, k: int, t: int, seed: int, hash_fn, seq_len: int) -> str | None:
    seq = [rng.randrange(4) for _ in range(seq_len)]
    bf = brute_force(seq, w, k, t, seed, hash_fn)
    naive = naive_incremental(seq, w, k, t, seed, hash_fn)
    inc = incremental(seq, w, k, t, seed, hash_fn)
    if bf != naive:
        for i, (b, nv) in enumerate(zip(bf, naive)):
            if b != nv:
                return f"BRUTE-FORCE vs NAIVE mismatch at index {i}: brute_force={b} naive={nv} (bf_len={len(bf)} naive_len={len(naive)})"
        return f"BRUTE-FORCE vs NAIVE length mismatch: bf_len={len(bf)} naive_len={len(naive)}"
    if naive != inc:
        for i, (nv, ic) in enumerate(zip(naive, inc)):
            if nv != ic:
                return f"NAIVE vs OPTIMIZED mismatch at index {i}: naive={nv} incremental={ic} (naive_len={len(naive)} inc_len={len(inc)})"
        return f"NAIVE vs OPTIMIZED length mismatch: naive_len={len(naive)} inc_len={len(inc)}"
    return None


def main() -> int:
    configs = [
        # (w, k, t, seq_len)
        (10, 15, 5, 4000),
        (50, 15, 5, 6000),   # the large-window regression case
        (11, 15, 5, 3000),
        (10, 21, 8, 4000),
        (5, 12, 4, 3000),
        (16, 20, 7, 5000),
        (3, 8, 2, 2000),     # tiny w, exercises frequent expire_a wraparound
        (2, 6, 2, 1500),     # w=2: pushed_a and expire_a alternate every step
        (1, 6, 3, 1200),     # w=1: a single residue class is both expire_a and pushed_a every step
        (100, 15, 5, 8000),  # w > tmer_span/2-ish, stress the cache further
        (10, 15, 14, 3000),  # t close to k (num t-mer loci per k-mer small)
        (10, 15, 5, 20000),  # longer sequence, same as our real w10_k15 config
    ]
    hash_fns = [("xor-hash", xor_hash), ("xor-pcg-hash", xor_pcg_hash)]
    seeds = [0x2545F4914F6CDD1D, 0x9E3779B97F4A7C15, 0x1234567890ABCDEF, 0, 1]

    n_cases = 0
    failures = []
    rng = random.Random(1234567)
    for (w, k, t, seq_len) in configs:
        for name, hash_fn in hash_fns:
            for seed in seeds:
                n_cases += 1
                msg = run_case(rng, w, k, t, seed, hash_fn, seq_len)
                if msg:
                    failures.append(f"(w={w},k={k},t={t},hash={name},seed={hex(seed)}): {msg}")

    print(f"ran {n_cases} (w,k,t,hash,seed) cases")
    if failures:
        print(f"{len(failures)} FAILURE(S):", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    print("all cases: incremental matches brute_force exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
