"""Theoretical expected DNA-alphabet density_factor for every signature-order
scheme (default, lex = true canonical lex, xor, xor-pcg-hash, greedymini), via Monte
Carlo simulation on an i.i.d. uniform-random DNA stream -- not measured on any
real dataset. xor/xor-pcg-hash/greedymini are seed-dependent (see benchmark/README.md:
"Run order and the lex baseline" for why seed matters), so this reports mean
and std across the benchmark's predefined seeds, not a single value.

Thin adapter over benchmark/common/expected_density_driver.py (loop, mean/std, cache,
CLI); this file only has kmc's own variants, DNA comparator keys and random
stream. Feeds the `orders` xlsx sheet (see scripts/export.py); cached to
`kmc3_density_expected_orders.csv` since it's independent of any actual
benchmark run (a pure function of k/p/w and the seed list) and each
simulation costs real Monte Carlo time.
"""

from __future__ import annotations

import random
import struct
from pathlib import Path

import expected_density_driver as shared
from density import dna_reverse_complement, flip_gm_ranks, pcg_rxs_m_xs64, sliding_runs  # noqa: F401
from scripts.config import seeds as bench_seeds
from scripts.density import unify_param_row, _cmp_keys_dna, resolve_gm_ranks

DEFAULT_N_BASES = 4_000_000
DEFAULT_SEED = 0xC0FFEE01
N_ORDER_SEEDS = 5  # matches this benchmark's default -n (see README)

ORDER_LABELS = ("default", "lex", "xor", "xor-pcg-hash", "greedymini")
PROJECTION_AT0 = "at0-cg1"
PROJECTION_AT1 = "at1-cg0"  # kmc -oat1
# One expected-density row per (order, projection); projection is "" for every
# order but greedymini, which has both projections (kmc -oat1 = at1-cg0).
ORDER_VARIANTS = tuple((o, "" if o != "greedymini" else PROJECTION_AT0) for o in ORDER_LABELS) + (
    ("greedymini", PROJECTION_AT1),
)


# Shared, identical implementations (common/density.py); the old local names stay as
# aliases (custom_order.cpp's rxs_m_xs64 / dna_reverse_complement).
rxs_m_xs64 = pcg_rxs_m_xs64
dna_rev = dna_reverse_complement


def read_dna_gm(path: Path) -> tuple[int, int, list[int]]:
    """Read a KMC _dna.gm file (see custom_order.cpp's KMC_DNA_GM_MAGIC
    format); return (k, seed, sort_key). This is loading a table already
    built by the C++ gm_sigma_extend tool, not reimplementing it."""
    with path.open("rb") as f:
        _magic, k, seed, n = struct.unpack("<QQQQ", f.read(32))
        sort_key = list(struct.unpack(f"<{n}I", f.read(4 * n)))
    return k, seed, sort_key


def _is_allowed(mmer: int, length: int) -> bool:
    """Port of CMmer::is_allowed (kmc_api/mmer.h) -- KMC's default-scheme denylist."""
    if (mmer & 0x3F) == 0x3F:
        return False
    if (mmer & 0x3F) == 0x3B:
        return False
    if (mmer & 0x3C) == 0x3C:
        return False
    m = mmer
    for _ in range(length - 3):
        if (m & 0xF) == 0:
            return False
        m >>= 2
    if m == 0:
        return False
    if m == 0x04:
        return False
    if (m & 0xF) == 0:
        return False
    return True


def _dna_sort_key_from_rank(rank_by_dna: list[int], p: int, seed: int) -> list[int]:
    """Canonical pick + total order, exactly matching custom_order.cpp's
    build_gm_dna_luts: every (p-mer, revcomp) pair collapses to one sort_key
    (canonical, both-strand counting), seed only breaks true rank ties."""
    n = 1 << (2 * p)
    order = sorted(range(n), key=lambda a: (rank_by_dna[a], a ^ seed))
    packed = [0] * n
    for pos, dna in enumerate(order):
        packed[dna] = pos

    def canonical(i: int) -> int:
        r = dna_rev(i, p)
        ri, rr = rank_by_dna[i], rank_by_dna[r]
        if ri != rr:
            return i if ri < rr else r
        ti, tr = i ^ seed, r ^ seed
        if ti != tr:
            return i if ti < tr else r
        return min(i, r)

    return [packed[canonical(i)] for i in range(n)]


def greedymini_dna_keys(
    k: int, p: int, seed: int, projection: str = PROJECTION_AT0, gm_dir: Path | None = None
) -> list[int]:
    """DNA comparator of a kmc greedymini run with tie-break `seed`: the binary
    .gm ranks (flipped for at1-cg0) extended to DNA by the same gm_sigma_extend
    tool kmc's table build uses. For at0-cg1 this reproduces the real
    `_dna.gm` sort_key bit for bit (tests/test_expected_density.py)."""
    if projection not in (PROJECTION_AT0, PROJECTION_AT1):
        raise ValueError(f"unknown projection: {projection}")
    try:
        _kb, ranks, _src = resolve_gm_ranks(k, p, gm_dir)
    except SystemExit as exc:  # die() on a missing table
        raise FileNotFoundError(str(exc)) from exc
    if projection == PROJECTION_AT1:
        ranks = flip_gm_ranks(ranks)
    return _cmp_keys_dna(ranks, p, seed & 0xFFFFFFFF)


def build_dna_keys(
    order: str,
    k: int,
    p: int,
    *,
    seed: int | None = None,
    gm_dir: Path | None = None,
    projection: str = PROJECTION_AT0,
) -> list[int]:
    """Build the per-DNA-p-mer comparator array (size 4^p) a real KMC run would
    use for `order`, matching kmc_api/mmer.h + custom_order.cpp exactly.
    `k` (the outer k-mer length) only matters for greedymini -- it's encoded
    in the GreedyMini table's on-disk filename (w{k-p+1}_k{p}_canonical...),
    even though the table itself, like every other order here, operates over
    p-mers (4^p entries)."""
    n = 1 << (2 * p)
    if order == "lex":
        # True canonical lex, no is_allowed denylist: MIN(i, revcomp(i)).
        return [min(i, dna_rev(i, p)) for i in range(n)]
    if order == "default":
        # KMC's shipped lex: same, but is_allowed()-excluded p-mers (and their
        # revcomp) are pushed to a sentinel so they can never be the minimizer.
        special = n
        keys = [0] * n
        for i in range(n):
            r = dna_rev(i, p)
            sv = i if _is_allowed(i, p) else special
            rv = r if _is_allowed(r, p) else special
            keys[i] = min(sv, rv)
        return keys
    if order == "xor":
        if seed is None:
            raise ValueError("xor order requires a seed")
        rank_by_dna = [(i ^ seed) & 0xFFFFFFFFFFFFFFFF for i in range(n)]
        return _dna_sort_key_from_rank(rank_by_dna, p, seed & 0xFFFFFFFF)
    if order == "xor-pcg-hash":
        if seed is None:
            raise ValueError("xor-pcg-hash order requires a seed")
        rank_by_dna = [rxs_m_xs64(i ^ seed) for i in range(n)]
        return _dna_sort_key_from_rank(rank_by_dna, p, seed & 0xFFFFFFFF)
    if order == "greedymini":
        # Both LUT backends are the same order (verified identical measured
        # density elsewhere). seed=None = the table's own embedded seed
        # (seeds[0]); the at0 result equals the real _dna.gm sort_key.
        return greedymini_dna_keys(k, p, bench_seeds()[0] if seed is None else seed, projection, gm_dir)
    raise ValueError(f"unknown order: {order}")


ROUNDING = shared.Rounding(sample=4, avg=4, std=4, single_std=None)


def expected_density(
    dna_keys: list[int], p: int, w: int, *, n_bases: int = DEFAULT_N_BASES, seed: int = DEFAULT_SEED
) -> float | None:
    """Monte Carlo expected density_factor for a precomputed comparator on an
    i.i.d. uniform-random DNA stream (no real sequence content, no read
    boundaries -- the theoretical floor a scheme's construction targets).
    kmc-specific: the stream is drawn as one 2p-bit seed state then 2 bits per
    base via getrandbits (not the shared `random_dna`), and the keys are a
    lookup of the p-mer state in the comparator `dna_keys` (already canonical)."""
    rng = random.Random(seed)
    mask = (1 << (2 * p)) - 1
    state = rng.getrandbits(2 * p) & mask
    keys = [dna_keys[state]]
    for _ in range(n_bases - p):
        state = ((state << 2) | rng.getrandbits(2)) & mask
        keys.append(dna_keys[state])
    return shared.window_density(keys, w, ROUNDING)


def _order_seeds(order: str) -> list[int | None]:
    if order in ("default", "lex"):
        return [None]  # deterministic, no seed dependence
    return list(bench_seeds()[:N_ORDER_SEEDS])


def compute_all_orders(
    k: int, p: int, w: int, *, gm_dir: Path | None = None, n_bases: int = DEFAULT_N_BASES
) -> list[dict]:
    """Expected density_factor (mean/std over seeds where applicable) for every
    (order, projection) variant. density_factor_std is None for deterministic
    orders (default, lex) since there's no seed to vary; the export shows 0."""

    def sample(v, _ctx, seed):
        order, projection = v
        try:
            keys = build_dna_keys(order, k, p, seed=seed, gm_dir=gm_dir, projection=projection or PROJECTION_AT0)
        except FileNotFoundError:
            raise shared.Unavailable from None
        return expected_density(keys, p, w, n_bases=n_bases)

    return shared.compute_rows(
        ORDER_VARIANTS,
        seeds_of=lambda v: _order_seeds(v[0]),
        sample=sample,
        row_id=lambda v: {"order": v[0], "projection": v[1], "k": p, "w": w, "L": k},
        rounding=ROUNDING,
    )


EXPECTED_ORDERS_FIELDS = ("order", "projection", "k", "w", "L", "n_seeds", "density_factor_mean", "density_factor_std")


def expected_orders_csv_path(results: Path) -> Path:
    return results / "kmc3_density_expected_orders.csv"


def _key(row: dict) -> tuple[str, str] | None:
    return (row["order"], row.get("projection") or "") if row.get("order") else None


def load_expected_orders(results: Path) -> dict[tuple[str, str], dict]:
    """Load the cache written by ensure_expected_orders, keyed by
    (order, projection); {} if not built yet."""
    return shared.read_csv_rows(expected_orders_csv_path(results), _key, unify_param_row)


def ensure_expected_orders(
    results: Path,
    k: int,
    p: int,
    w: int,
    *,
    gm_dir: Path | None = None,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
) -> dict[tuple[str, str], dict]:
    """Return the cached per-variant expected-density table, computing (and
    caching) it on first use. Independent of any actual benchmark run, so it's
    safe/cheap to reuse across every export once built. A cache with another
    variant set (e.g. pre-projection) or another k/p is stale."""
    return shared.ensure_cache(
        expected_orders_csv_path(results),
        EXPECTED_ORDERS_FIELDS,
        ORDER_VARIANTS,
        key_of_variant=lambda v: v,
        key_of_row=_key,
        load=lambda: load_expected_orders(results),
        compute=lambda _todo: compute_all_orders(k, p, w, gm_dir=gm_dir, n_bases=n_bases),
        accept=lambda cached: all(int(r.get("L") or 0) == k and int(r.get("k") or 0) == p for r in cached.values()),
        force=force,
    )


def main(argv: list[str] | None = None) -> int:
    import sys

    from scripts.common import load_dotenv
    from paths import bench_paths

    load_dotenv()
    default_results = bench_paths("kmc")["results"]

    ap = shared.base_parser(
        "scripts expected-density",
        "Theoretical expected density_factor per signature-order scheme "
        "(Monte Carlo on i.i.d. random DNA; mean/std over predefined seeds for "
        "xor/xor-pcg-hash/greedymini). Feeds the `orders` xlsx sheet's expected_density_factor_mean/std "
        "columns -- run this once (or after changing k/p) to (re)build that cache.",
        default_results=str(default_results),
        default_n_bases=DEFAULT_N_BASES,
    )
    ap.add_argument("--k", type=int, default=28)
    ap.add_argument("-p", type=int, default=9)
    ap.add_argument("--gm-dir", default=None)
    args = ap.parse_args(argv)

    w = args.k - args.p + 1
    results = Path(args.results_dir).resolve() / "csvs"  # --results-dir is the results root
    gm_dir = Path(args.gm_dir) if args.gm_dir else None

    cached = ensure_expected_orders(
        results, args.k, args.p, w, gm_dir=gm_dir, force=args.force, n_bases=args.n_bases
    )
    print(f"wrote {expected_orders_csv_path(results)} ({len(cached)} orders)", file=sys.stderr)
    for label, projection in ORDER_VARIANTS:
        r = cached.get((label, projection))
        if r:
            print(shared.format_row(f"{label} {projection}".strip(), r, 18, "expected_density"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
