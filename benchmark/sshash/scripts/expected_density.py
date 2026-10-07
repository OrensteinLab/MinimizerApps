"""Theoretical expected density_factor for every sshash order (native,
xor-hash, xor-pcg-hash, greedymini-<label>[-at1]), via Monte Carlo simulation
on an i.i.d. uniform-random DNA stream -- not measured on any real dataset.
Mirrors benchmark/kmc/scripts/expected_density.py's approach and naming,
ported to sshash's own order formulas (see scripts/density.py).

xor-hash/xor-pcg-hash/greedymini are seed-dependent, so this reports mean and
std across the benchmark's predefined seeds (see config.yaml), not a single
value; native's own seed only changes an internal mixing constant, not its
statistical behavior as a random order (see density.py:native_hash's
docstring), but it's still varied across seeds here for a fair
apples-to-apples comparison.

Feeds the `orders` xlsx sheet's expected_density_avg/std columns (see
scripts/export.py); cached to sshash_density_expected_orders.csv since it's
independent of any actual benchmark run (a pure function of m/w + the seed
list) and each simulation costs real Monte Carlo time.
"""

from __future__ import annotations

import csv
import random
import statistics
from collections.abc import Iterable
from pathlib import Path

from scripts.density import (
    canonical_mmer,
    density_factor,
    _gm_k_from_table,
    gm_centre_valid,
    greedymini_hash,
    mod_sampling_runs,
    native_hash,
    read_gm_ranks,
    xor_hash,
    xor_pcg_hash,
    _sliding_runs,
)
from scripts.orders import ORDERS, mod_t, parse_greedymini_order, parse_mod_order

DEFAULT_N_BASES = 300_000
# the mod-minimizer's Monte Carlo is O(k) per k-mer in pure Python
DEFAULT_MOD_N_BASES = 60_000
DEFAULT_SEED = 0xC0FFEE01


def _hash_fn(order: str, m: int, gm_ranks: list[int] | None, sigma_flip: bool = False):
    """Return `f(value, seed) -> hash` for `order`, hashing a raw
    (uncanonicalized) m-mer."""
    if order == "native":
        return native_hash
    if order == "xor-hash":
        return xor_hash
    if order == "xor-pcg-hash":
        return xor_pcg_hash
    assert gm_ranks is not None
    return lambda value, seed: greedymini_hash(value, m, gm_ranks, seed, sigma_flip)


def expected_density_factor(
    order: str,
    m: int,
    w: int,
    *,
    seed: int,
    n_bases: int = DEFAULT_N_BASES,
    gm_ranks: list[int] | None = None,
    rng_seed: int = DEFAULT_SEED,
) -> float | None:
    rng = random.Random(rng_seed)
    mod = parse_mod_order(order)
    if mod is not None:
        n_bases = min(n_bases, DEFAULT_MOD_N_BASES)
    dna_stream = [rng.randrange(4) for _ in range(n_bases)]
    if mod is not None:
        k = w + m - 1
        t = mod_t(order, k, m)
        if mod[0] == "xor":
            fn = xor_hash
        else:
            assert gm_ranks is not None
            fn = lambda value, seed: greedymini_hash(value, t, gm_ranks, seed)  # noqa: E731
        n_kmers, n_super = mod_sampling_runs(dna_stream, k, m, t, fn, seed)
        return density_factor(n_super, n_kmers, w)
    plain = parse_greedymini_order(order)
    sigma_flip = plain is not None and plain[1] == "at1-cg0"
    if plain is not None and gm_ranks is not None and not gm_centre_valid(m, _gm_k_from_table(gm_ranks)):
        return None  # needs m - gm_k even (sshash rejects it too)

    # one hash per canonicalized m-mer position, plain sliding-window argmin
    # (see _sliding_runs).
    hash_fn = _hash_fn(order, m, gm_ranks, sigma_flip)
    n = len(dna_stream)
    mask = (1 << (2 * m)) - 1
    state = 0
    for i in range(m):
        state = (state << 2) | dna_stream[i]
    keys = [hash_fn(canonical_mmer(state, m), seed)]
    for i in range(m, n):
        state = ((state << 2) | dna_stream[i]) & mask
        keys.append(hash_fn(canonical_mmer(state, m), seed))
    n_win, n_super, _n_sel = _sliding_runs(keys, w)
    return density_factor(n_super, n_win, w)


def compute_all_orders(
    m: int,
    w: int,
    seeds: list[int],
    *,
    n_bases: int = DEFAULT_N_BASES,
    orders: Iterable[str] | None = None,
) -> list[dict]:
    """Monte Carlo rows for `orders` (default: every order in ORDERS)."""
    from scripts.config import greedymini_label_file, mod_greedymini_label_file

    wanted = set(orders) if orders is not None else None
    out: list[dict] = []
    for order in ORDERS:
        if wanted is not None and order not in wanted:
            continue
        plain = parse_greedymini_order(order)
        mod = parse_mod_order(order)
        gm_ranks = None
        if plain is not None or (mod is not None and mod[0] == "greedymini"):
            path = greedymini_label_file(plain[0]) if plain else mod_greedymini_label_file(mod[1])
            gm_ranks = read_gm_ranks(path) if path.is_file() else None
        needs_table = plain is not None or (mod is not None and mod[0] == "greedymini")
        vals: list[float] = []
        if not needs_table or gm_ranks is not None:
            for seed in seeds:
                v = expected_density_factor(order, m, w, seed=seed, n_bases=n_bases, gm_ranks=gm_ranks)
                if v is not None:
                    vals.append(v)
        if not vals:
            out.append({"order": order, "m": m, "w": w, "n_seeds": 0,
                        "density_factor_avg": None, "density_factor_std": None})
            continue
        out.append({"order": order, "m": m, "w": w, "n_seeds": len(vals),
                    "density_factor_avg": round(statistics.mean(vals), 8),
                    "density_factor_std": round(statistics.pstdev(vals), 10) if len(vals) > 1 else 0.0})
    return out


EXPECTED_ORDERS_FIELDS = ("order", "m", "w", "n_seeds", "density_factor_avg", "density_factor_std")


def expected_orders_csv_path(results: Path) -> Path:
    return results / "sshash_density_expected_orders.csv"


def load_expected_orders(results: Path) -> dict[str, dict]:
    path = expected_orders_csv_path(results)
    if not path.is_file():
        return {}
    return {row["order"]: row for row in csv.DictReader(path.open()) if row.get("order")}


def ensure_expected_orders(
    results: Path,
    m: int,
    w: int,
    seeds: list[int],
    *,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
    orders: Iterable[str] | None = None,
) -> dict[str, dict]:
    """Cached expected density per order. With `orders` given, only those are
    computed (the pure-Python Monte Carlo is slow), merged into whatever the
    cache already holds for this m/w."""
    path = expected_orders_csv_path(results)
    wanted_orders = set(orders) if orders is not None else set(ORDERS)
    kept: list[dict] = []
    if not force and path.is_file():
        cached = load_expected_orders(results)
        if cached and all(int(r.get("m") or 0) == m and int(r.get("w") or 0) == w for r in cached.values()):
            if wanted_orders <= set(cached):
                return cached
            kept = [dict(r) for r in cached.values()]
    have = {r["order"] for r in kept}
    rows = kept + compute_all_orders(m, w, seeds, n_bases=n_bases, orders=wanted_orders - have)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(EXPECTED_ORDERS_FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    return {r["order"]: {k: ("" if v is None else str(v)) for k, v in r.items()} for r in rows}


def main(argv: list[str] | None = None) -> int:
    import argparse
    import sys

    from scripts.common import bench_paths, load_dotenv
    from scripts.config import seeds as bench_seeds

    load_dotenv()
    _, default_results = bench_paths()

    ap = argparse.ArgumentParser(
        prog="scripts expected-density",
        description="Theoretical expected density_factor per sshash order (Monte Carlo "
        "on i.i.d. random DNA; mean/std over the benchmark's predefined seeds). Feeds the "
        "`orders` xlsx sheet's expected_density_avg/std columns.",
    )
    ap.add_argument("--m", type=int, default=21)
    ap.add_argument("--k", type=int, default=31, help="only used to derive w=k-m+1")
    ap.add_argument("--results-dir", default=str(default_results))
    ap.add_argument("--n-bases", type=int, default=DEFAULT_N_BASES)
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)

    w = args.k - args.m + 1
    results = Path(args.results_dir).resolve()
    seeds = list(bench_seeds())[: args.n_seeds]

    cached = ensure_expected_orders(results, args.m, w, seeds, force=args.force, n_bases=args.n_bases)
    print(f"wrote {expected_orders_csv_path(results)} ({len(cached)} rows)", file=sys.stderr)
    for order in ORDERS:
        r = cached.get(order)
        if not r:
            continue
        std = r.get("density_factor_std") or None
        std_str = f" +/- {std}" if std else ""
        print(f"  {order:36s} n_seeds={r.get('n_seeds')}  expected_density_factor={r.get('density_factor_avg')}{std_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
