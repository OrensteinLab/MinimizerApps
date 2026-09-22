"""Minimizer density helpers: generic sliding-window "particular" density
counting, plus pure-Python re-implementations of sshash's four hash/order
formulas used by ``expected_density.py`` to Monte Carlo simulate the
*theoretical* expected density_factor of each order on i.i.d. random DNA
(mirrors ``benchmark/kmc/scripts/density.py``'s approach and terminology,
ported for sshash's own orders since kmc's file re-derives kmc-specific
comparator/denylist logic that doesn't apply here).

density_factor = density * (w + 1), where density = n_selected / n_kmers and
w = k - m + 1 is the window size (see orders/GreedyMini/w21_k11_canonical.*:
the "w21" in that filename *is* this w, for k=31/m=11). A uniformly random
total order has expected density_factor -> 2 as w grows (Roberts et al.
2004); GreedyMini's whole point is to push it below that.
"""

from __future__ import annotations

import struct
from collections import deque
from pathlib import Path

MASK64 = (1 << 64) - 1
EVENS64 = 0x5555555555555555
ODDS64 = 0xAAAAAAAAAAAAAAAA
INF_RANK = 0xFFFF


def density_factor(n_selected: int, n_kmers: int, w: int) -> float | None:
    if n_kmers <= 0 or n_selected <= 0:
        return None
    return round(n_selected / n_kmers * (w + 1), 6)


def _sliding_runs(keys: list[int], w: int) -> tuple[int, int, int]:
    """Slide a length-w window over `keys` (one hash value per m-mer position);
    return (n_windows, n_super, n_selected_positions). n_super counts how many
    times the window's argmin *value* changes (super-k-mer boundaries); this
    is the generic "particular density" building block re-derivable from any
    (dataset or simulated) key stream, ported verbatim from kmc's
    scripts/density.py (that helper has no kmc-specific content)."""
    n = len(keys)
    n_win = max(0, n - w + 1)
    dq: deque[int] = deque()
    n_super = 0
    n_sel = 0
    prev_key = None
    prev_pos = None
    for i in range(n):
        while dq and keys[dq[-1]] >= keys[i]:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - w:
            dq.popleft()
        if i + 1 < w:
            continue
        pos = dq[0]
        min_key = keys[pos]
        if prev_key is None or min_key != prev_key:
            n_super += 1
            prev_key = min_key
        if prev_pos != pos:
            n_sel += 1
            prev_pos = pos
    return n_win, n_super, n_sel


# ---------------------------------------------------------------------------
# .gm file I/O (matches common/minimizer-orders/greedymini_rank.h exactly)
# ---------------------------------------------------------------------------


def read_gm_ranks(path: Path) -> list[int]:
    """Read a classic .gm file ([u64 len][len x u64] raw ranks; entries equal
    to `len` are the "not in UHS" sentinel) and return the dense uint16-style
    rank table (INF_RANK for non-UHS), matching GreedyMiniTable::load_from_raw."""
    data = path.read_bytes()
    (n,) = struct.unpack_from("<Q", data, 0)
    raw = struct.unpack_from(f"<{n}Q", data, 8)
    return [INF_RANK if r >= n else r for r in raw]


# ---------------------------------------------------------------------------
# Bit helpers (Python port of common/minimizer-orders/greedymini_rank.h,
# assuming the traditional ACGT=0,1,2,3 packing this fork always builds with)
# ---------------------------------------------------------------------------


def dna_reverse_complement(mmer: int, m: int) -> int:
    """Reverse complement of an m-base 2-bit-per-base DNA m-mer, ACGT=0,1,2,3
    (complement(x) = 3-x under this mapping)."""
    rev = 0
    shift = m * 2 - 2
    x = mmer
    for _ in range(m):
        rev |= (3 - (x & 3)) << shift
        x >>= 2
        shift -= 2
    return rev


def gm_acgt_to_atcg(x: int) -> int:
    x &= MASK64
    msb = (x >> 1) & EVENS64
    lsb = x & EVENS64
    return ((msb ^ lsb) << 1 | msb) & MASK64


def _pext(x: int, mask: int) -> int:
    res = 0
    bitpos = 0
    m = mask
    while m:
        lsb = m & (-m)
        if x & lsb:
            res |= 1 << bitpos
        bitpos += 1
        m &= m - 1
    return res


def gm_key_and_gamma(atcg: int, m: int) -> tuple[int, int]:
    """gm_k == m (no Extend-k, this fork's only configuration): key = all
    Sigma (odd) bits, gamma = all Gamma (even) bits, of the low 2*m bits."""
    low_mask = (1 << (2 * m)) - 1
    key = _pext(atcg, ODDS64 & low_mask)
    gamma = _pext(atcg, EVENS64 & low_mask)
    return key, gamma


# ---------------------------------------------------------------------------
# Order hash functions (pure-Python ports; see
# sshash/include/custom_order_hasher.hpp for the authoritative C++)
# ---------------------------------------------------------------------------


def pcg_rxs_m_xs64(x: int) -> int:
    x &= MASK64
    count = x >> 59
    x ^= x >> (5 + count)
    x = (x * 12605985483714917081) & MASK64
    return x ^ (x >> 43)


def xor_hash(value: int, seed: int) -> int:
    return (value ^ seed) & MASK64


def xor_pcg_hash(value: int, seed: int) -> int:
    return pcg_rxs_m_xs64((value ^ seed) & MASK64)


def native_hash(value: int, seed: int) -> int:
    """Approximates sshash's NATIVE order ((value*mult)^magic, magic =
    pthash::xxhash_64::hash(seed,0)). We use splitmix64(seed) as a stand-in
    for the exact xxhash64-derived magic: for i.i.d. random DNA, expected
    density depends only on the order being an unstructured (roughly
    uniform-random) total order, not on which specific avalanching mix
    produced it, so this is statistically equivalent for this Monte Carlo's
    purpose without needing a full xxhash64 port."""
    z = (seed + 0x9E3779B97F4A7C15) & MASK64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    magic = z ^ (z >> 31)
    return ((value * 0x517CC1B727220A95) ^ magic) & MASK64


def greedymini_hash(value: int, m: int, gm_ranks: list[int], seed: int) -> int:
    """Full port of custom_order_hasher::rank_binary_dna for gm_k == m (no
    Extend-k): dual-orientation lookup, raw Sigma key folded below rank as a
    tie-break (see that file's comment for why)."""
    key_bits = m
    gamma_bits = m
    toggle = seed & ((1 << gamma_bits) - 1)

    def probe(x: int) -> tuple[int, int, int]:
        atcg = gm_acgt_to_atcg(x)
        key, gamma = gm_key_and_gamma(atcg, m)
        return gm_ranks[key], key, gamma

    def finish(rank: int, key: int, gamma: int) -> int:
        mid = (rank << key_bits) | key
        g = gamma ^ toggle
        return (mid << gamma_bits) | g

    ra, ka, ga = probe(value)
    rc = dna_reverse_complement(value, m)
    if rc == value:
        return finish(ra, ka, ga)
    rb, kb, gb = probe(rc)
    if ra < rb:
        return finish(ra, ka, ga)
    if rb < ra:
        return finish(rb, kb, gb)
    va, vb = finish(ra, ka, ga), finish(rb, kb, gb)
    return va if va < vb else vb


def canonical_mmer(x: int, m: int) -> int:
    """sshash's own pre-canonicalization: numeric min(forward, revcomp)."""
    rc = dna_reverse_complement(x, m)
    return x if x <= rc else rc
