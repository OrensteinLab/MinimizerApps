"""Tool-agnostic minimizer/signature density helpers, shared by every app's
own ``scripts/density.py`` (sshash, kmc, minimap2). Extracted from what was
previously duplicated near-verbatim between ``benchmark/kmc/scripts/density.py``
and ``benchmark/sshash/scripts/density.py`` -- see this repo's git history
for the pre-extraction copies.

density_factor = density * (w + 1), where density = n_selected / n_kmers and
w is the window size (see orders/GreedyMini/w11_k15_canonical.*: the "w11" in
that filename *is* this w). A uniformly random total order has expected
density_factor -> 2 as w grows (Roberts et al. 2004); GreedyMini's whole point
is to push it below that, and mod-minimizer (Groot Koerkamp & Pibiri, WABI
2024) pushes it lower still as the sampled substring grows relative to w.

Kept out of this module (still per-app, in each app's own scripts/density.py)
because they encode an app-specific packing/parsing choice, not a shared
algorithm:
  - the app's own "native" order approximation (each app's default order is
    its own, unshared, hash construction).
  - the GreedyMini rank's *packing* into a return value: kraken2/kmc/sshash
    fold Extend-k + a gamma tiebreak below the UHS primary (see
    common/minimizer-orders/greedymini_rank.h's own doc comment and
    orders/README.md's "Extend-k" section); minimap2's simplified ranker
    (minimap2/order.cpp) instead packs a PCG-scrambled tiebreak with no
    Extend-k, since it never needs the result to double as an index key.
    Build each app's own greedymini_hash on top of the bit-math primitives
    below, matching that app's actual ranker.
"""

from __future__ import annotations

import struct
from collections import deque
from pathlib import Path

MASK64 = (1 << 64) - 1
EVENS64 = 0x5555555555555555
ODDS64 = 0xAAAAAAAAAAAAAAAA
INF_RANK = 0xFFFF


def density_factor(n_selected: int, n_kmers: int, w: int, decimals: int = 6) -> float | None:
    if n_kmers <= 0 or n_selected <= 0:
        return None
    return round(n_selected / n_kmers * (w + 1), decimals)


def sliding_runs(keys: list[int], w: int) -> tuple[int, int, int]:
    """Slide a length-w window over `keys` (one hash value per m-mer
    position); return (n_windows, n_super, n_selected_positions). n_super
    counts how many times the window's argmin *value* changes (super-k-mer
    boundaries) -- the generic "particular density" building block,
    re-derivable from any (dataset or simulated) key stream, independent of
    which order produced `keys`."""
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


def gm_k_from_table(gm_ranks: list[int]) -> int:
    """Binary-alphabet GreedyMini tables have exactly 2**gm_k entries."""
    n = len(gm_ranks)
    gm_k = n.bit_length() - 1
    assert (1 << gm_k) == n, f"GreedyMini table size {n} is not a power of 2"
    return gm_k


# ---------------------------------------------------------------------------
# Bit helpers (Python port of common/minimizer-orders/greedymini_rank.h,
# assuming the traditional ACGT=0,1,2,3 packing every app in this repo builds
# with)
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


def canonical_mmer(x: int, m: int) -> int:
    """Numeric min(forward, revcomp) -- the "canonical parsing" convention
    (as opposed to a single-strand/"regular parsing" convention some apps use
    for mod-minimizer specifically; see each app's own density.py)."""
    rc = dna_reverse_complement(x, m)
    return x if x <= rc else rc


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
# Order hash functions (pure-Python ports of common/minimizer-orders/
# pcg_rxs_m_xs.h; xor-hash and xor-pcg-hash are identical across every app)
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


# ---------------------------------------------------------------------------
# Mod-sampling (Groot Koerkamp & Pibiri, WABI 2024, Algorithm 6 / Definition 7)
# ---------------------------------------------------------------------------


def rolling_values(dna_stream: list[int], length: int) -> list[int]:
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
    dna_stream: list[int], k: int, m: int, t: int, hash_fn, seed: int,
    canonical: bool = False,
) -> tuple[int, int]:
    """Mod-sampling particular-density counting: for each k-mer window, scan
    its k-t+1 t-mer loci with `hash_fn` (= O_t), apply the paper's exact tie
    rule (equal hash -> prefer the smaller i mod w), and count super-kmer
    boundaries (how many times the selected absolute position's raw m-mer
    value changes) -- the same (n_windows, n_super) shape `sliding_runs`
    returns for the classic scheme, so `density_factor(n_super, n_windows, w)`
    applies unchanged.

    `canonical=False` (default) matches a single-strand/"regular parsing"
    construction (sshash's own mod-minimizer, per sshash/include/util.hpp's
    compute_mod_minimizer_single_strand doc comment: canonicalizing per-locus
    isn't mirror-equivariant for mod-sampling's tie rule, so sshash keeps
    whichever strand construction walked). `canonical=True` matches a
    canonical (min(fwd,rc)) per-t-mer-locus construction (kraken2::ModSampler
    and minimap2/order.cpp's sketch_mod_minimizer, which only ever need a
    single forward-genome pass and canonicalize each t-mer locus directly).

    `hash_fn` is called as hash_fn(t_mer_value, seed); pass a closure over
    each app's own greedymini_hash for the *-greedymini variant.
    """
    w = k - m + 1
    num_t_loci = k - t + 1
    n_win = max(0, len(dna_stream) - k + 1)
    if n_win <= 0:
        return 0, 0
    t_raw = rolling_values(dna_stream, t)
    m_raw = rolling_values(dna_stream, m)

    n_super = 0
    prev_val = None
    for start in range(n_win):
        o_min = None
        p = 0
        for i in range(num_t_loci):
            v = t_raw[start + i]
            if canonical:
                v = canonical_mmer(v, t)
            o = hash_fn(v, seed)
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
