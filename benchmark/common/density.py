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
# GreedyMini centre extension, injective packing (port of
# common/minimizer-orders/greedymini_centre.h gm_centre::Order with
# injective=true and TwoBitPlanes<true>, i.e. the ACGT = 0..3 encoding that
# kraken2 and minimap2 use; benchmark/sshash/scripts/density.py has the same
# order on sshash's own encoding). l = gm_k + 2i; first base = most
# significant bit of each plane.
# ---------------------------------------------------------------------------


def _reverse_low_bits(v: int, n: int) -> int:
    return int(format(v, f"0{n}b")[::-1], 2) if n else 0


def greedymini_centre_hash(lmer: int, l: int, gm_ranks: list[int], gm_k: int, flip: bool = False) -> int:
    """Value of an ACGT=0123 l-mer: rank of the centre key's Sigma plane in the
    table, then min over the two strands of (Sigma << l | Gamma). Sigma = AT vs
    CG (bit0 ^ bit1, complemented when `flip`: --gm-projection at1-cg0),
    Gamma = bit 1. Strand invariant; a total order on canonical l-mers."""
    i = (l - gm_k) // 2
    mask = (1 << l) - 1
    sigma = _pext(lmer ^ (lmer >> 1), EVENS64 & ((1 << (2 * l)) - 1))
    if flip:
        sigma ^= mask
    gamma = _pext(lmer, ODDS64 & ((1 << (2 * l)) - 1))
    rank = gm_ranks[(sigma >> i) & ((1 << gm_k) - 1)]
    sigma_rc = _reverse_low_bits(sigma, l)
    gamma_rc = _reverse_low_bits(gamma, l) ^ mask
    return (rank << (2 * l)) | min((sigma << l) | gamma, (sigma_rc << l) | gamma_rc)


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


def mod_sampling_density(
    dna_stream: list[int], k: int, w: int, t: int, hash_fn, seed: int,
) -> tuple[int, int]:
    """(n_kmers, n_selected) of the strand-symmetric mod-minimizer as
    minimap2's mod-minimizer-xor and kraken2's mod-minimizer define it
    (minimap2/order.cpp, kraken2/src/mod_minimizer.h):
    per window of L = w+k-1 bases the preferred strand has more G/T (tie: the
    window <= its reverse complement); forward-preferred -> k-mer offset
    x_L mod w (x_L = leftmost minimal canonical t-mer), reverse-preferred ->
    (w-1) - (((T-1) - x_R) mod w) (x_R = rightmost), T = w+k-t. Not a forward
    scheme, so this counts distinct selected positions instead of
    super-k-mer boundaries; density_factor(n_selected, n_kmers, w)."""
    n = len(dna_stream)
    L, T = w + k - 1, w + k - t
    ranks = [hash_fn(canonical_mmer(v, t), seed) for v in rolling_values(dna_stream, t)]
    gt = sum(1 for b in dna_stream[:L] if b >= 2)
    selected: set[int] = set()
    for s in range(n - L + 1):
        if s:
            gt += (dna_stream[s + L - 1] >= 2) - (dna_stream[s - 1] >= 2)
        if 2 * gt != L:
            fwd = 2 * gt > L
        else:
            win = dna_stream[s:s + L]
            fwd = win <= [3 - b for b in reversed(win)]
        r = ranks[s:s + T]
        m = min(r)
        if fwd:
            off = r.index(m) % w
        else:
            off = (w - 1) - (((T - 1) - (T - 1 - r[::-1].index(m))) % w)
        selected.add(s + off)
    return max(0, n - k + 1), len(selected)
