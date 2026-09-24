"""Minimizer density helpers: generic sliding-window "particular" density
counting, plus pure-Python re-implementations of sshash's four hash/order
formulas used by ``expected_density.py`` to Monte Carlo simulate the
*theoretical* expected density_factor of each order on i.i.d. random DNA
(mirrors ``benchmark/kmc/scripts/density.py``'s approach and terminology,
ported for sshash's own orders since kmc's file re-derives kmc-specific
comparator/denylist logic that doesn't apply here).

density_factor = density * (w + 1), where density = n_selected / n_kmers and
w = k - m + 1 is the window size (see orders/GreedyMini/w11_k15_canonical.*:
the "w11" in that filename *is* this w, for k=31/m=21). A uniformly random
total order has expected density_factor -> 2 as w grows (Roberts et al.
2004); GreedyMini's whole point is to push it below that, and mod-minimizer
(Groot Koerkamp & Pibiri, WABI 2024) pushes it lower still as m grows relative
to w -- see mod_sampling_runs below.
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


def _gm_k_from_table(gm_ranks: list[int]) -> int:
    """Binary-alphabet GreedyMini tables have exactly 2**gm_k entries."""
    n = len(gm_ranks)
    gm_k = n.bit_length() - 1
    assert (1 << gm_k) == n, f"GreedyMini table size {n} is not a power of 2"
    return gm_k


def gm_extract_dna_key_atcg(atcg: int, l: int, gm_k: int) -> int:
    """Sigma of the TOP gm_k bases of the l-length atcg value (Extend-k: the
    key is the prefix, the remaining l-gm_k bases are the suffix below)."""
    top = (atcg >> (2 * (l - gm_k))) & MASK64
    mask = (1 << (2 * gm_k)) - 1
    top &= mask
    return _pext(top, ODDS64 & mask)


def gm_extract_suffix_sigma_atcg(atcg: int, l: int, gm_k: int) -> int:
    """Sigma of the low (l-gm_k) bases (the Extend-k suffix); 0 when gm_k==l."""
    suf = l - gm_k
    if suf <= 0:
        return 0
    mask = (1 << (2 * suf)) - 1
    return _pext(atcg & mask, ODDS64 & mask)


def gm_extract_dna_tb_atcg(atcg: int, l: int) -> int:
    """Gamma (tie-break) of all l bases."""
    mask = (1 << (2 * l)) - 1
    return _pext(atcg & mask, EVENS64 & mask)


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


def greedymini_hash(value: int, l: int, gm_ranks: list[int], seed: int) -> int:
    """Full port of custom_order_hasher::rank_binary_dna, including Extend-k
    (gm_k, derived from the table's own size, may be < l -- e.g. plain
    greedymini uses gm_k=15 at l=m=21; mod-minimizer-greedymini uses gm_k=11
    at l=t=11, ext_bits=0): dual-orientation lookup, raw Sigma key (+ suffix,
    if extended) folded below rank as a tie-break (see that file's comment
    for why)."""
    gm_k = _gm_k_from_table(gm_ranks)
    ext_bits = l - gm_k
    key_bits = gm_k
    gamma_bits = l  # 1 bit/base (DNA, bpc-1=1), all l bases
    toggle = seed & ((1 << gamma_bits) - 1)

    def probe(x: int) -> tuple[int, int, int]:
        atcg = gm_acgt_to_atcg(x)
        key = gm_extract_dna_key_atcg(atcg, l, gm_k)
        return gm_ranks[key], key, atcg

    def finish(rank: int, key: int, atcg: int) -> int:
        mid = ((rank << key_bits) | key) << ext_bits
        if ext_bits > 0:
            mid |= gm_extract_suffix_sigma_atcg(atcg, l, gm_k)
        gamma = gm_extract_dna_tb_atcg(atcg, l) ^ toggle
        return (mid << gamma_bits) | gamma

    ra, ka, aa = probe(value)
    rc = dna_reverse_complement(value, l)
    if rc == value:
        return finish(ra, ka, aa)
    rb, kb, ab = probe(rc)
    if ra < rb:
        return finish(ra, ka, aa)
    if rb < ra:
        return finish(rb, kb, ab)
    va, vb = finish(ra, ka, aa), finish(rb, kb, ab)
    return va if va < vb else vb


def canonical_mmer(x: int, m: int) -> int:
    """sshash's own pre-canonicalization: numeric min(forward, revcomp)."""
    rc = dna_reverse_complement(x, m)
    return x if x <= rc else rc


# ---------------------------------------------------------------------------
# Mod-sampling (Groot Koerkamp & Pibiri, WABI 2024, Algorithm 6 / Definition 7)
# ---------------------------------------------------------------------------


def _rolling_values(dna_stream: list[int], length: int) -> list[int]:
    """Raw (uncanonicalized) packed integer value of each `length`-window of
    `dna_stream` (2-bit bases), one entry per valid start position."""
    n = len(dna_stream)
    if n < length:
        return []
    mask = (1 << (2 * length)) - 1
    state = 0
    for i in range(length):
        state = (state << 2) | dna_stream[i]
    out = [state]
    for i in range(length, n):
        state = ((state << 2) | dna_stream[i]) & mask
        out.append(state)
    return out


def mod_sampling_runs(
    dna_stream: list[int], k: int, m: int, t: int, hash_fn, seed: int
) -> tuple[int, int]:
    """Mod-sampling particular-density counting on a single (uncanonicalized)
    strand -- matches sshash's own "regular parsing" construction for
    mod-minimizer (see sshash/include/util.hpp's
    compute_mod_minimizer_single_strand doc comment): for each k-mer window,
    scan its k-t+1 t-mer loci with `hash_fn` (= O_t), apply the paper's exact
    tie rule (equal hash -> prefer the smaller i mod w), and count super-kmer
    boundaries (how many times the selected absolute position's raw m-mer
    value changes) -- the same (n_windows, n_super) shape `_sliding_runs`
    returns for the classic scheme, so `density_factor(n_super, n_windows, w)`
    applies unchanged.

    `hash_fn` is one of xor_hash/xor_pcg_hash/native_hash (called as
    hash_fn(t_mer_value, seed)) or a `lambda v, s: greedymini_hash(v, t,
    gm_ranks, s)` closure for mod-minimizer-greedymini (whose GreedyMini
    table has gm_k == t, so no Extend-k is needed here either).
    """
    w = k - m + 1
    num_t_loci = k - t + 1
    n_win = max(0, len(dna_stream) - k + 1)
    if n_win <= 0:
        return 0, 0
    t_raw = _rolling_values(dna_stream, t)
    m_raw = _rolling_values(dna_stream, m)

    n_super = 0
    prev_val = None
    for start in range(n_win):
        o_min = None
        p = 0
        for i in range(num_t_loci):
            o = hash_fn(t_raw[start + i], seed)
            pos_mod_w = i % w
            if o_min is None or o < o_min:
                o_min = o
                p = pos_mod_w
            elif o == o_min and pos_mod_w < p:
                p = pos_mod_w
        val = m_raw[start + p]
        if prev_val is None or val != prev_val:
            n_super += 1
            prev_val = val
    return n_win, n_super
