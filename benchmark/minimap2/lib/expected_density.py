#!/usr/bin/env python3
"""Theoretical expected density_factor for every minimap2 order (native,
xor-hash, xor-pcg-hash, greedymini, mod-minimizer-xor-hash,
mod-minimizer-pcg-xor-hash), via Monte Carlo simulation on an i.i.d.
uniform-random DNA stream -- not measured on any real dataset. Mirrors
benchmark/sshash/scripts/expected_density.py's own approach and naming,
ported to minimap2's own order formulas (minimap2/order.cpp), reusing
benchmark/common/density.py's shared, tool-agnostic primitives (the same
ones order.cpp's own C++ is built on -- gm_acgt_to_atcg/
gm_extract_dna_key_atcg/dna_reverse_complement/pcg_rxs_m_xs64/xor_hash/
xor_pcg_hash) rather than re-deriving them a third time.

xor-hash/xor-pcg-hash/greedymini (and their mod-minimizer-* counterparts)
are seed-dependent, so this reports mean and std across the benchmark's
predefined seeds (see config.yaml), not a single value; native has no seed
at all in minimap2 (mm_hash64 takes none), so it's computed once and
repeated n_seeds times as a formality for a uniform table shape.

Feeds the `orders` xlsx sheet's expected_density_factor_mean/std columns
(see lib/export.py); cached to minimap2_density_expected_orders.csv since
it's independent of any actual benchmark run (a pure function of k/w/t +
the seed list) and each simulation costs real Monte Carlo time.
mod-minimizer costs O(k-t+1) more per window than the 4 base orders (it
scans every t-mer locus in each k-mer window, not just one canonicalized
k-mer), so its Monte Carlo uses a smaller default n_bases (see
DEFAULT_MOD_N_BASES) to keep the one-time cache-build tractable.
"""
from __future__ import annotations

import csv
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from density import (  # noqa: E402
    canonical_mmer,
    density_factor,
    dna_reverse_complement,
    gm_acgt_to_atcg,
    gm_extract_dna_key_atcg,
    mod_sampling_runs,
    pcg_rxs_m_xs64,
    read_gm_ranks,
    sliding_runs,
    xor_hash,
    xor_pcg_hash,
)

DEFAULT_N_BASES = 300_000
DEFAULT_MOD_N_BASES = 40_000
DEFAULT_SEED = 0xC0FFEE01

ORDERS = [
    "native", "xor-hash", "xor-pcg-hash", "greedymini",
    "mod-minimizer-xor-hash", "mod-minimizer-pcg-xor-hash",
]
MOD_ORDERS = frozenset({"mod-minimizer-xor-hash", "mod-minimizer-pcg-xor-hash"})

# Which base order each order uses as its hash/rank function -- mod-minimizer-*
# reuses xor-hash/xor-pcg-hash as O_t, exactly as order.cpp's own
# order_rank_tmer() dispatches (see that function's own switch statement).
_BASE_ORDER = {
    "native": "native",
    "xor-hash": "xor-hash",
    "xor-pcg-hash": "xor-pcg-hash",
    "greedymini": "greedymini",
    "mod-minimizer-xor-hash": "xor-hash",
    "mod-minimizer-pcg-xor-hash": "xor-pcg-hash",
}


def mm_hash64(key: int, mask: int) -> int:
    """Direct port of sketch.c's mm_hash64 (Thomas Wang's 64-bit integer
    mix) -- minimap2's own 'native' order, called on the already-
    canonicalized (min of fwd/rc) k-mer, same convention as every other
    order here. `~key` relies on Python's infinite-precision two's
    complement: `(~key) & mask` gives the same bit pattern C's fixed-width
    `~key` would, for any mask width -- a standard, safe idiom, not
    something specific to 64 bits."""
    key = (~key + (key << 21)) & mask
    key = key ^ (key >> 24)
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ^ (key >> 14)
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ^ (key >> 28)
    key = (key + (key << 31)) & mask
    return key


def greedymini_hash(lmer: int, l: int, gm_ranks: list[int], seed: int) -> int:
    """Direct port of order.cpp's greedymini_rank(): dual-orientation UHS
    rank probe (gm_k == l always -- minimap2 disallows Extend-k, checked
    at load time) packed above a PCG-scrambled tiebreak, exactly as
    minimap2 computes and stores it. `seed` is unused (greedymini's ranks
    come from its table, not a seed) -- kept only for hash_fn signature
    parity with xor_hash/xor_pcg_hash."""
    gm_k = l
    atcg = gm_acgt_to_atcg(lmer)
    best = gm_ranks[gm_extract_dna_key_atcg(atcg, l, gm_k)]
    rc = dna_reverse_complement(lmer, l)
    if rc != lmer:
        atcg_rc = gm_acgt_to_atcg(rc)
        p_rc = gm_ranks[gm_extract_dna_key_atcg(atcg_rc, l, gm_k)]
        if p_rc < best:
            best = p_rc
    tiebreak = pcg_rxs_m_xs64(lmer) & 0xFFFFFFFFFFFF
    return (best << 48) | tiebreak


def _hash_fn(base_order: str, length: int, gm_ranks: list[int] | None):
    """Return f(value, seed) -> rank for `base_order`, given an already-
    canonicalized `length`-mer -- used directly as k-mer keys for the 4
    base orders, or as t-mer keys (O_t) for mod-minimizer's mod_sampling_runs."""
    if base_order == "native":
        mask = (1 << (2 * length)) - 1
        return lambda value, seed: mm_hash64(value, mask)
    if base_order == "xor-hash":
        return xor_hash
    if base_order == "xor-pcg-hash":
        return xor_pcg_hash
    if base_order == "greedymini":
        assert gm_ranks is not None
        return lambda value, seed: greedymini_hash(value, length, gm_ranks, seed)
    raise ValueError(f"unknown base order: {base_order}")


def hpc_compress(dna_stream: list[int]) -> list[int]:
    """Collapse consecutive equal bases to one -- the input transformation
    `-H`/is_hpc applies before k-merizing (sketch.c's own tiny_queue_t run-
    length tracking additionally records each k-mer's *span* in bases, but
    density_factor only depends on window count and selected-position
    count, both fully determined by the compressed *values* alone)."""
    return [b for i, b in enumerate(dna_stream) if i == 0 or b != dna_stream[i - 1]]


def expected_density_factor(
    order: str,
    k: int,
    w: int,
    *,
    seed: int,
    t: int | None = None,
    n_bases: int = DEFAULT_N_BASES,
    gm_ranks: list[int] | None = None,
    rng_seed: int = DEFAULT_SEED,
    is_hpc: bool = True,
) -> float | None:
    rng = random.Random(rng_seed)
    dna_stream = [rng.randrange(4) for _ in range(n_bases)]
    if is_hpc:
        dna_stream = hpc_compress(dna_stream)
    base_order = _BASE_ORDER[order]

    if order in MOD_ORDERS:
        assert t is not None
        hash_fn = _hash_fn(base_order, t, gm_ranks)
        # mod_sampling_runs's own (k, m) naming is (outer window span in
        # bases, final sampled length); minimap2's own (w, k) mod-minimizer
        # samples k-mers (m_param = k) over a window of w consecutive
        # k-mers (k_param = w + k - 1, so that its internally-derived
        # `w = k_param - m_param + 1` comes back out to minimap2's own w).
        k_param = w + k - 1
        n_win, n_super = mod_sampling_runs(dna_stream, k_param, k, t, hash_fn, seed, canonical=True)
        return density_factor(n_super, n_win, w)

    # classic scheme: one hash per canonicalized k-mer position, plain
    # sliding-window argmin (see sliding_runs).
    hash_fn = _hash_fn(base_order, k, gm_ranks)
    n = len(dna_stream)
    mask = (1 << (2 * k)) - 1
    state = 0
    for i in range(k):
        state = (state << 2) | dna_stream[i]
    keys = [hash_fn(canonical_mmer(state, k), seed)]
    for i in range(k, n):
        state = ((state << 2) | dna_stream[i]) & mask
        keys.append(hash_fn(canonical_mmer(state, k), seed))
    n_win, n_super, _n_sel = sliding_runs(keys, w)
    return density_factor(n_super, n_win, w)


def compute_all_orders(
    k: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    n_bases: int = DEFAULT_N_BASES,
    mod_n_bases: int = DEFAULT_MOD_N_BASES,
    mod_t: int = 5,
    is_hpc: bool = True,
) -> list[dict]:
    gm_ranks = read_gm_ranks(gm_file) if gm_file and gm_file.is_file() else None

    out: list[dict] = []
    for order in ORDERS:
        is_mod = order in MOD_ORDERS
        needs_gm = _BASE_ORDER[order] == "greedymini"
        t = mod_t if is_mod else None
        if needs_gm and gm_ranks is None:
            out.append({
                "order": order, "t": t, "k": k, "w": w, "n_seeds": 0,
                "density_factor_avg": None, "density_factor_std": None,
            })
            continue
        n_bases_for_order = mod_n_bases if is_mod else n_bases
        vals: list[float] = []
        for seed in seeds:
            v = expected_density_factor(
                order, k, w, seed=seed, t=t, n_bases=n_bases_for_order, gm_ranks=gm_ranks, is_hpc=is_hpc
            )
            if v is not None:
                vals.append(v)
        if not vals:
            out.append({
                "order": order, "t": t, "k": k, "w": w, "n_seeds": 0,
                "density_factor_avg": None, "density_factor_std": None,
            })
            continue
        avg = round(statistics.mean(vals), 4)
        std = round(statistics.pstdev(vals), 4) if len(vals) > 1 else 0.0
        out.append({
            "order": order, "t": t, "k": k, "w": w, "n_seeds": len(vals),
            "density_factor_avg": avg, "density_factor_std": std,
        })
    return out


EXPECTED_ORDERS_FIELDS = ("order", "t", "k", "w", "n_seeds", "density_factor_avg", "density_factor_std")


def expected_orders_csv_path(results: Path) -> Path:
    return results / "minimap2_density_expected_orders.csv"


def load_expected_orders(results: Path) -> dict[str, dict]:
    path = expected_orders_csv_path(results)
    if not path.is_file():
        return {}
    return {row["order"]: row for row in csv.DictReader(path.open()) if row.get("order")}


def ensure_expected_orders(
    results: Path,
    k: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
    mod_n_bases: int = DEFAULT_MOD_N_BASES,
    mod_t: int = 5,
    is_hpc: bool = True,
) -> dict[str, dict]:
    """is_hpc must match whatever the actual benchmark runs use (lib/run.sh
    always passes -x map-pb, which enables -H/HPC k-mers -- see that
    file's own comment) for 'particular' (measured) and 'expected'
    (this function) density_factor to be a fair comparison; not tracked
    in the cache key below since this harness has no non-HPC mode today."""
    path = expected_orders_csv_path(results)
    if not force and path.is_file():
        cached = load_expected_orders(results)
        if cached and all(int(r.get("k") or 0) == k and int(r.get("w") or 0) == w for r in cached.values()):
            return cached
    rows = compute_all_orders(k, w, seeds, gm_file=gm_file, n_bases=n_bases, mod_n_bases=mod_n_bases, mod_t=mod_t, is_hpc=is_hpc)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(EXPECTED_ORDERS_FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    return {r["order"]: {k2: ("" if v is None else str(v)) for k2, v in r.items()} for r in rows}


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(
        description="Theoretical expected density_factor per minimap2 order (Monte Carlo "
        "on i.i.d. random DNA; mean/std over the benchmark's predefined seeds). Feeds the "
        "`orders` xlsx sheet's expected_density_factor_mean/std columns."
    )
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--w", type=int, default=10)
    ap.add_argument("--mod-t", type=int, default=5)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--gm-file", default=None)
    ap.add_argument("--n-bases", type=int, default=DEFAULT_N_BASES)
    ap.add_argument("--mod-n-bases", type=int, default=DEFAULT_MOD_N_BASES)
    ap.add_argument("--seeds", default="0x2545F4914F6CDD1D,0x9E3779B97F4A7C15,0xDEADBEEFCAFEBABE,0x1234567890ABCDEF,0xA5A5A5A5A5A5A5A5")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)

    seeds = [int(s, 16) if s.lower().startswith("0x") else int(s) for s in args.seeds.split(",")]
    rows = ensure_expected_orders(
        Path(args.results_dir), args.k, args.w, seeds,
        gm_file=Path(args.gm_file) if args.gm_file else None,
        force=args.force, n_bases=args.n_bases, mod_n_bases=args.mod_n_bases, mod_t=args.mod_t,
    )
    for order, row in rows.items():
        print(order, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
