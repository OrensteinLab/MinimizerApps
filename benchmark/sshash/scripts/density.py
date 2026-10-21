"""Minimizer density helpers: generic sliding-window "particular" density
counting, plus pure-Python re-implementations of sshash's hash/order
formulas used by ``expected_density.py`` to Monte Carlo simulate the
*theoretical* expected density_factor of each order on i.i.d. random DNA
(mirrors ``benchmark/kmc/scripts/density.py``'s approach and terminology,
ported for sshash's own orders since kmc's file re-derives kmc-specific
comparator/denylist logic that doesn't apply here).

density_factor = density * (w + 1), where density = n_selected / n_kmers and
w = k - m + 1 is the window size (see orders/GreedyMini/w11_k15_canonical.*:
the "w11" in that filename *is* this w, for k=31/m=21). A uniformly random
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
    return round(n_selected / n_kmers * (w + 1), 10)


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


def mod_sampling_runs(dna_stream: list[int], k: int, m: int, t: int, hash_fn,
                      seed: int) -> tuple[int, int]:
    """(n_kmers, n_super) of the strand-symmetric mod-minimizer
    (util::compute_mod_minimizer in sshash/include/util.hpp): per k-mer, the
    T = k - t + 1 canonical t-mers are hashed by `hash_fn`, and the centre rule of
    the shared reference (benchmark/common/density.py mod_centre_selections) picks
    the sampled m-mer, two equally close offsets by the numeric order of the packed
    k-mer and its reverse complement (sshash's !(kmer_rc < kmer)). n_super counts
    changes of the sampled position along the stream -- sshash's super-k-mers.

    t-mers and k-mers are packed as sshash packs them (kmer.hpp: drop_char shifts
    the first base out of the low bits), first base least significant -- not this
    module's MSB-first packing of the plain orders -- so the canonical t-mer (and
    hence its xor hash) and the strand tie match the compiled order exactly."""
    from density import mod_centre_selections  # benchmark/common

    n = len(dna_stream)
    w = k - m + 1
    T = k - t + 1
    if n < k:
        return 0, 0
    top = 2 * (t - 1)
    state = 0
    ranks: list[int] = []
    for i, b in enumerate(dna_stream):
        state = (state >> 2) | (b << top)  # first base in the low bits (sshash's packing)
        if i + 1 >= t:
            ranks.append(hash_fn(canonical_mmer(state, t), seed))

    def le_rc(s: int) -> bool:
        kv = 0
        for j, b in enumerate(dna_stream[s:s + k]):
            kv |= b << (2 * j)
        return kv <= dna_reverse_complement(kv, k)  # !(kmer_rc < kmer)

    n_super = 0
    prev = None
    for pos in mod_centre_selections(ranks, w, T, le_rc):
        if pos != prev:
            n_super += 1
            prev = pos
    return n - k + 1, n_super


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
# Bit helpers (Python port of sshash/include/custom_order_hasher.hpp, on sshash's
# default packing A=0,C=1,G=3,T=2: bit 0 of a base is Sigma (AT vs CG), bit 1 is
# Gamma, and complementing a base flips bit 1)
# ---------------------------------------------------------------------------


def dna_reverse_complement(mmer: int, m: int) -> int:
    """Reverse complement of an m-base 2-bit-per-base DNA m-mer, A=0,C=1,G=3,T=2
    (complement(x) = x ^ 2 under this mapping)."""
    rev = 0
    shift = m * 2 - 2
    x = mmer
    for _ in range(m):
        rev |= ((x & 3) ^ 2) << shift
        x >>= 2
        shift -= 2
    return rev


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


def gm_centre_valid(l: int, gm_k: int) -> bool:
    """The extension of the table to the m-mer needs the same number of extra bases on both sides."""
    return l >= gm_k and (l - gm_k) % 2 == 0


def gm_table_mirror_symmetric(gm_ranks: list[int]) -> bool:
    gm_k = _gm_k_from_table(gm_ranks)
    return all(
        gm_ranks[key] == gm_ranks[int(format(key, f"0{gm_k}b")[::-1], 2)] for key in range(len(gm_ranks))
    )


def greedymini_hash(value: int, l: int, gm_ranks: list[int], seed: int, sigma_flip: bool = False) -> int:
    """Port of custom_order_hasher::hash_from_planes (the GreedyMini hash).

    View the m-mer as x y z with |x| = |z| = i = (l - gm_k) / 2, |y| = gm_k.
    Primary key: rank of Sigma(y) (one probe; strand-invariant because the
    table is mirror-symmetric and rc maps y to its mirror). Ties: the smaller
    of (Sigma, Gamma^toggle) of the two strands, compared as (Sigma << l) |
    Gamma -- i.e. min(Sigma(xyz), reverse Sigma(xyz)) first, then Gamma -- so
    the packed value is injective on canonical m-mers."""
    gm_k = _gm_k_from_table(gm_ranks)
    assert gm_centre_valid(l, gm_k), f"centre extension needs l - gm_k even (l={l}, gm_k={gm_k})"
    i = (l - gm_k) // 2
    msk = (1 << (2 * l)) - 1
    flip = ((1 << l) - 1) if sigma_flip else 0
    toggle = seed & ((1 << l) - 1)

    def planes(v: int) -> tuple[int, int]:
        return _pext(v, EVENS64 & msk) ^ flip, _pext(v, ODDS64 & msk) ^ toggle

    s1, g1 = planes(value)
    s2, g2 = planes(dna_reverse_complement(value, l))
    key = (s1 >> (l - i - gm_k)) & ((1 << gm_k) - 1)
    return (gm_ranks[key] << (2 * l)) | min((s1 << l) | g1, (s2 << l) | g2)


def canonical_mmer(x: int, m: int) -> int:
    """sshash's own pre-canonicalization: numeric min(forward, revcomp)."""
    rc = dna_reverse_complement(x, m)
    return x if x <= rc else rc
