"""sshash-specific density helpers: its own NATIVE-order approximation and
its GreedyMini rank packing (Extend-k + gamma tiebreak, matching
sshash/include/custom_order_hasher.hpp), built on the tool-agnostic
primitives in ``benchmark/common/density.py`` (density_factor, sliding_runs,
.gm file I/O, the GreedyMini bit-math helpers, pcg_rxs_m_xs64/xor_hash/
xor_pcg_hash, and mod_sampling_runs) -- extracted there since kmc's own
density.py duplicated all of these near-verbatim (that file keeps its own
comparator/denylist logic that doesn't apply here; see that module's
docstring for exactly what stayed per-app and why).
"""

from __future__ import annotations

from density import (  # noqa: F401 -- re-exported for existing importers
    MASK64,
    EVENS64,
    ODDS64,
    INF_RANK,
    density_factor,
    read_gm_ranks,
    dna_reverse_complement,
    canonical_mmer,
    gm_acgt_to_atcg,
    gm_extract_dna_key_atcg,
    gm_extract_suffix_sigma_atcg,
    gm_extract_dna_tb_atcg,
    pcg_rxs_m_xs64,
    xor_hash,
    xor_pcg_hash,
    rolling_values,
)
from density import gm_k_from_table as _gm_k_from_table
from density import sliding_runs as _sliding_runs
from density import mod_sampling_runs as _shared_mod_sampling_runs


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


# ---------------------------------------------------------------------------
# Mod-sampling: thin wrapper over the shared implementation, pinned to
# canonical=False (sshash's own single-strand/"regular parsing" construction
# -- see the shared function's docstring in benchmark/common/density.py).
# ---------------------------------------------------------------------------


def mod_sampling_runs(
    dna_stream: list[int], k: int, m: int, t: int, hash_fn, seed: int
) -> tuple[int, int]:
    return _shared_mod_sampling_runs(dna_stream, k, m, t, hash_fn, seed, canonical=False)
