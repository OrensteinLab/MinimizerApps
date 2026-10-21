#!/usr/bin/env python3
"""Theoretical expected density_factor for every minimap2 order (native,
xor, xor-pcg-hash, greedymini in each projection, mod-minimizer-xor), via
Monte Carlo simulation on an i.i.d.
uniform-random DNA stream -- not measured on any real dataset. Thin adapter over
benchmark/common/expected_density_driver.py (loop, mean/std, cache, CLI):
this file only has minimap2's variants, order hash functions (minimap2/order.cpp:
mm_hash64 'native', the GreedyMini centre hash, mod-minimizer-xor) and HPC
compression, on benchmark/common/density.py's shared primitives.

xor/xor-pcg-hash/mod-minimizer-xor are seed-dependent, so this reports mean and std across the benchmark's
predefined seeds (see config.yaml), not a single value; native has no seed
at all in minimap2 (mm_hash64 takes none), so it's computed once and
repeated n_seeds times as a formality for a uniform table shape.

Feeds the `orders` xlsx sheet's expected_density_factor_mean/std columns
(see lib/export.py); cached to minimap2_density_expected_k<k>_w<w>_t<mod_t>_<hpc|nohpc>.csv since
it's independent of any actual benchmark run (a pure function of k/w/t/HPC,
the GreedyMini table, the seed list and the Monte Carlo sizes) and each
simulation costs real Monte Carlo time. Every cached row carries a
cache_key over everything but k/w/t/HPC (which are in the file name): the
GreedyMini table's path and content hash, the seeds, n_bases/mod_n_bases
and the Monte Carlo RNG seed, so changing any of them recomputes the cache.
(Both GreedyMini projections get their own rows: each row carries its
projection, gm_k and gm_table columns.) A cache
written before cache_key existed simply misses once.
mod-minimizer costs O(k-t+1) more per window than the 4 base orders (it
scans every t-mer locus in each k-mer window, not just one canonicalized
k-mer), so its Monte Carlo uses a smaller default n_bases (see
DEFAULT_MOD_N_BASES) to keep the one-time cache-build tractable.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import expected_density_driver as shared  # noqa: E402
from order_names import gm_fields, projection_from_flag  # noqa: E402
from density import (  # noqa: E402
    canonical_mmer,
    gm_k_from_table,
    greedymini_centre_hash,
    mod_sampling_density,
    read_gm_ranks,
    xor_hash,
    xor_pcg_hash,
)

DEFAULT_N_BASES = 300_000
DEFAULT_MOD_N_BASES = 40_000
DEFAULT_SEED = shared.DEFAULT_SEED
# density_factor of one draw is rounded to 6 decimals (shared default), mean/std to 4
ROUNDING = shared.Rounding(sample=6, avg=4, std=4, single_std=0.0)

# (order, projection) variants: greedymini once per projection, the others
# have no projection (None).
VARIANTS = [
    ("native", None), ("xor", None), ("xor-pcg-hash", None),
    ("greedymini", "at0-cg1"), ("greedymini", "at1-cg0"),
    ("mod-minimizer-xor", None),
]
MOD_ORDERS = frozenset({"mod-minimizer-xor"})

# Which base order each order uses as its hash/rank function --
# mod-minimizer-xor ranks t-mers by xor, as order.cpp's sketch_mod_minimizer.
_BASE_ORDER = {
    "native": "native",
    "xor": "xor",
    "xor-pcg-hash": "xor-pcg-hash",
    "greedymini": "greedymini",
    "mod-minimizer-xor": "xor",
}


# App-specific (not in common/density.py): minimap2's own 'native' hash and the HPC
# input transformation exist only in minimap2.


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


def greedymini_hash(lmer: int, l: int, gm_ranks: list[int], seed: int, flip: bool = False) -> int:
    """order.cpp's RankState<MM_ORDER_GREEDYMINI>: the shared GreedyMini centre extension
    of the gm_k table to length l (l = gm_k + 2i), injective packing --
    benchmark/common/density.py greedymini_centre_hash, checked against the
    compiled order by minimap2/validate_greedymini.py. `seed` = --order-seed,
    the Gamma tie-break toggle. `flip` = --gm-projection at1-cg0."""
    return greedymini_centre_hash(lmer, l, gm_ranks, gm_k_from_table(gm_ranks), flip, seed)


def _hash_fn(base_order: str, length: int, gm_ranks: list[int] | None, flip: bool = False):
    """Return f(value, seed) -> rank for `base_order`, given an already-
    canonicalized `length`-mer -- used directly as k-mer keys for the 4
    base orders, or as t-mer keys (O_t) for mod-minimizer's mod_sampling_density."""
    if base_order == "native":
        mask = (1 << (2 * length)) - 1
        return lambda value, seed: mm_hash64(value, mask)
    if base_order == "xor":
        return xor_hash
    if base_order == "xor-pcg-hash":
        return xor_pcg_hash
    if base_order == "greedymini":
        assert gm_ranks is not None
        return lambda value, seed: greedymini_hash(value, length, gm_ranks, seed, flip)
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
    projection: str | None = None,
) -> float | None:
    dna_stream = shared.random_dna(n_bases, rng_seed)
    if is_hpc:
        dna_stream = hpc_compress(dna_stream)
    base_order = _BASE_ORDER[order]
    flip = projection_from_flag(projection) == "at1-cg0" if base_order == "greedymini" else False

    if order in MOD_ORDERS:
        # strand-symmetric mod-sampling: distinct selected positions (not a
        # forward scheme, so no super-k-mer boundary count)
        assert t is not None
        n_kmers, n_sel = mod_sampling_density(dna_stream, k, w, t, _hash_fn(base_order, t, gm_ranks, flip), seed)
        return ROUNDING.factor(n_sel, n_kmers, w)

    # classic scheme: one hash per canonicalized k-mer position, plain sliding-window argmin
    keys = shared.canonical_keys(dna_stream, k, _hash_fn(base_order, k, gm_ranks, flip), canonical_mmer, seed)
    return shared.window_density(keys, w, ROUNDING)


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
    gm_info = gm_fields(table=gm_file.name, projection=None) if gm_file else {"gm_k": None, "gm_table": None}

    def row_id(v):
        order, projection = v
        needs_gm = _BASE_ORDER[order] == "greedymini"
        return {
            "order": order, "gm_k": gm_info["gm_k"] if needs_gm else None, "projection": projection,
            "gm_table": gm_info["gm_table"] if needs_gm else None, "mod_t": mod_t if order in MOD_ORDERS else None,
            "k": k, "w": w, "L": k + w - 1,
        }

    def prepare(v):
        if _BASE_ORDER[v[0]] == "greedymini" and gm_ranks is None:
            raise shared.Unavailable  # no table: n_seeds 0
        return gm_ranks

    return shared.compute_rows(
        VARIANTS,
        seeds_of=lambda v: seeds,
        prepare=prepare,
        sample=lambda v, ranks, seed: expected_density_factor(
            v[0], k, w, seed=seed, t=mod_t if v[0] in MOD_ORDERS else None,
            n_bases=mod_n_bases if v[0] in MOD_ORDERS else n_bases, gm_ranks=ranks, is_hpc=is_hpc,
            projection=v[1]),
        row_id=row_id,
        rounding=ROUNDING,
    )


EXPECTED_ORDERS_FIELDS = ("order", "gm_k", "projection", "gm_table", "mod_t", "k", "w", "L", *shared.STAT_FIELDS, "cache_key")


def expected_cache_key(
    seeds: list[int],
    *,
    gm_file: Path | None,
    n_bases: int,
    mod_n_bases: int,
    rng_seed: int = DEFAULT_SEED,
) -> str:
    """Identity of every cache input not already in the file name (k, w, t,
    HPC): the GreedyMini table (resolved path + content sha256; 'missing' when
    absent, so building the table later invalidates a cache computed without
    it), the seed list, the Monte Carlo sizes and RNG seed, and a version of
    the order hashes."""
    gm = "none"
    if gm_file is not None:
        gm_path = Path(gm_file).resolve()
        digest = hashlib.sha256(gm_path.read_bytes()).hexdigest() if gm_path.is_file() else "missing"
        gm = f"{gm_path}:{digest}"
    ident = {
        # bump when an order's hash or the row layout changes (2: greedymini --order-seed
        # toggles Gamma; 3: projection/gm_k/gm_table columns; 4: mod_t and L columns;
        # 5: the mod-minimizer's centre tie rule replaces the G/T-strand rule)
        "hash_version": 5,
        "gm": gm,
        "seeds": [hex(s) for s in seeds],
        "n_bases": n_bases,
        "mod_n_bases": mod_n_bases,
        "rng_seed": hex(rng_seed),
    }
    return hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]


def expected_orders_csv_path(results: Path, k: int, w: int, t: int, is_hpc: bool) -> Path:
    return results / "csvs" / f"minimap2_density_expected_k{k}_w{w}_t{t}_{'hpc' if is_hpc else 'nohpc'}.csv"


def variant_key(row: dict) -> tuple[str, str]:
    """(order, projection or '') -- the key of one expected-density row."""
    return (row["order"], row.get("projection") or "")


def _row_key(row: dict) -> tuple[str, str] | None:
    return variant_key(row) if row.get("order") else None


def load_expected_orders(results: Path, k: int, w: int, t: int, is_hpc: bool) -> dict[tuple[str, str], dict]:
    return shared.read_csv_rows(expected_orders_csv_path(results, k, w, t, is_hpc), _row_key)


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
) -> dict[tuple[str, str], dict]:
    """Keyed by (order, projection or ''). is_hpc must match the dataset's preset (lib/constants.sh:
    dataset_preset; map-pb enables -H/HPC k-mers, map-ont does not) for
    'particular' (measured) and 'expected' (this function) density_factor
    to be a fair comparison: map-pb datasets use the HPC value, map-ont ones
    the non-HPC one; each is cached in its own file."""
    key = expected_cache_key(seeds, gm_file=gm_file, n_bases=n_bases, mod_n_bases=mod_n_bases)
    return shared.ensure_cache(
        expected_orders_csv_path(results, k, w, mod_t, is_hpc),
        EXPECTED_ORDERS_FIELDS,
        VARIANTS,
        key_of_variant=lambda v: (v[0], v[1] or ""),
        key_of_row=_row_key,
        load=lambda: load_expected_orders(results, k, w, mod_t, is_hpc),
        compute=lambda _todo: compute_all_orders(k, w, seeds, gm_file=gm_file, n_bases=n_bases,
                                                 mod_n_bases=mod_n_bases, mod_t=mod_t, is_hpc=is_hpc),
        accept=lambda cached: all(r.get("cache_key") == key for r in cached.values()),
        stamp={"cache_key": key},
        force=force,
    )


def main(argv: list[str] | None = None) -> int:
    ap = shared.base_parser(
        None,
        "Theoretical expected density_factor per minimap2 order (Monte Carlo "
        "on i.i.d. random DNA; mean/std over the benchmark's predefined seeds). Feeds the "
        "`orders` xlsx sheet's expected_density_factor_mean/std columns.",
        default_results=None,
        default_n_bases=DEFAULT_N_BASES,
    )
    ap.add_argument("--k", type=int, default=19, help="19 = map-pb's k; 15 for map-ont")
    ap.add_argument("--w", type=int, default=10)
    ap.add_argument("--mod-t", type=int, default=None, help="default 4 + (k-4) %% w")
    ap.add_argument("--no-hpc", action="store_true", help="raw k-mers (map-ont) instead of HPC (map-pb)")
    ap.add_argument("--gm-file", default=None)
    ap.add_argument("--mod-n-bases", type=int, default=DEFAULT_MOD_N_BASES)
    ap.add_argument("--seeds", default="0x2545F4914F6CDD1D,0x9E3779B97F4A7C15,0xDEADBEEFCAFEBABE,0x1234567890ABCDEF,0xA5A5A5A5A5A5A5A5")
    args = ap.parse_args(argv)

    seeds = [int(s, 16) if s.lower().startswith("0x") else int(s) for s in args.seeds.split(",")]
    rows = ensure_expected_orders(
        Path(args.results_dir), args.k, args.w, seeds,
        gm_file=Path(args.gm_file) if args.gm_file else None,
        force=args.force, n_bases=args.n_bases, mod_n_bases=args.mod_n_bases,
        mod_t=args.mod_t if args.mod_t is not None else 4 + (args.k - 4) % args.w,
        is_hpc=not args.no_hpc,
    )
    for key, row in rows.items():
        print(key, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
