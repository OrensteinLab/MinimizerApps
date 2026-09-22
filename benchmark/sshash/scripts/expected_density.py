"""Theoretical expected density_factor for every sshash order (native,
xor-hash, xor-pcg-hash, greedymini, and the 4 mod-minimizer-* orders), via
Monte Carlo simulation on an i.i.d. uniform-random DNA stream -- not measured
on any real dataset. Mirrors benchmark/kmc/scripts/expected_density.py's
approach and naming, ported to sshash's own order formulas (see
scripts/density.py).

xor-hash/xor-pcg-hash/greedymini (and their mod-minimizer-* counterparts) are
seed-dependent, so this reports mean and std across the benchmark's
predefined seeds (see config.yaml), not a single value; native's own seed
only changes an internal mixing constant, not its statistical behavior as a
random order (see density.py:native_hash's docstring), but it's still varied
across seeds here for a fair apples-to-apples comparison.

Feeds the `orders` xlsx sheet's expected_density_avg/std columns (see
scripts/export.py); cached to sshash_density_expected_orders.csv since it's
independent of any actual benchmark run (a pure function of k/m/w/t + the
seed list) and each simulation costs real Monte Carlo time. Mod-minimizer
rows cost O(k-t+1) more than the 4 base orders (Algorithm 6 scans every t-mer
locus in each k-mer window, not just m), so their Monte Carlo uses a smaller
default n_bases (see DEFAULT_MOD_N_BASES) to keep the one-time cache-build
tractable.
"""

from __future__ import annotations

import csv
import random
import statistics
from pathlib import Path

from scripts.density import (
    canonical_mmer,
    density_factor,
    greedymini_hash,
    mod_sampling_runs,
    native_hash,
    read_gm_ranks,
    xor_hash,
    xor_pcg_hash,
    _sliding_runs,
)
from scripts.orders import MOD_ORDERS, ORDERS, order_t_values

DEFAULT_N_BASES = 300_000
DEFAULT_MOD_N_BASES = 40_000
DEFAULT_SEED = 0xC0FFEE01

# Base order each order uses as its hash/order function (mod-minimizer-*
# reuses the 4 base orders as O_t -- see orders/README.md).
_BASE_ORDER = {
    "native": "native",
    "xor-hash": "xor-hash",
    "xor-pcg-hash": "xor-pcg-hash",
    "greedymini": "greedymini",
    "mod-minimizer-native-hash": "native",
    "mod-minimizer-xor-hash": "xor-hash",
    "mod-minimizer-pcg-xor-hash": "xor-pcg-hash",
    "mod-minimizer-greedymini": "greedymini",
}


def _hash_fn(base_order: str, length: int, gm_ranks: list[int] | None):
    """Return `f(value, seed) -> hash` for `base_order`, hashing a raw
    (uncanonicalized) `length`-mer -- used directly as m-mer keys for the 4
    base orders, or as t-mer keys (O_t) for mod-minimizer's `mod_sampling_runs`."""
    if base_order == "native":
        return native_hash
    if base_order == "xor-hash":
        return xor_hash
    if base_order == "xor-pcg-hash":
        return xor_pcg_hash
    if base_order == "greedymini":
        assert gm_ranks is not None
        return lambda value, seed: greedymini_hash(value, length, gm_ranks, seed)
    raise ValueError(f"unknown base order: {base_order}")


def expected_density_factor(
    order: str,
    m: int,
    w: int,
    *,
    seed: int,
    t: int | None = None,
    n_bases: int = DEFAULT_N_BASES,
    gm_ranks: list[int] | None = None,
    rng_seed: int = DEFAULT_SEED,
) -> float | None:
    rng = random.Random(rng_seed)
    dna_stream = [rng.randrange(4) for _ in range(n_bases)]
    base_order = _BASE_ORDER[order]
    k = m + w - 1

    if order in MOD_ORDERS:
        assert t is not None
        hash_fn = _hash_fn(base_order, t, gm_ranks)
        n_win, n_super = mod_sampling_runs(dna_stream, k, m, t, hash_fn, seed)
        return density_factor(n_super, n_win, w)

    # classic scheme: one hash per canonicalized m-mer position, plain
    # sliding-window argmin (see _sliding_runs).
    hash_fn = _hash_fn(base_order, m, gm_ranks)
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
    gm_file: Path | None = None,
    mod_t_gm_file: Path | None = None,
    n_bases: int = DEFAULT_N_BASES,
    mod_n_bases: int = DEFAULT_MOD_N_BASES,
) -> list[dict]:
    gm_ranks = read_gm_ranks(gm_file) if gm_file and gm_file.is_file() else None
    mod_gm_ranks = read_gm_ranks(mod_t_gm_file) if mod_t_gm_file and mod_t_gm_file.is_file() else None

    out: list[dict] = []
    for order in ORDERS:
        is_mod = order in MOD_ORDERS
        ranks_for_order = mod_gm_ranks if is_mod else gm_ranks
        needs_gm = _BASE_ORDER[order] == "greedymini"
        for t in order_t_values(order):
            if needs_gm and ranks_for_order is None:
                out.append(
                    {"order": order, "t": t, "m": m, "w": w, "n_seeds": 0,
                     "density_factor_avg": None, "density_factor_std": None}
                )
                continue
            n_bases_for_order = mod_n_bases if is_mod else n_bases
            vals: list[float] = []
            for seed in seeds:
                v = expected_density_factor(
                    order, m, w, seed=seed, t=t, n_bases=n_bases_for_order, gm_ranks=ranks_for_order
                )
                if v is not None:
                    vals.append(v)
            if not vals:
                out.append(
                    {"order": order, "t": t, "m": m, "w": w, "n_seeds": 0,
                     "density_factor_avg": None, "density_factor_std": None}
                )
                continue
            avg = round(statistics.mean(vals), 4)
            std = round(statistics.pstdev(vals), 4) if len(vals) > 1 else 0.0
            out.append(
                {"order": order, "t": t, "m": m, "w": w, "n_seeds": len(vals),
                 "density_factor_avg": avg, "density_factor_std": std}
            )
    return out


EXPECTED_ORDERS_FIELDS = ("order", "t", "m", "w", "n_seeds", "density_factor_avg", "density_factor_std")


def expected_orders_csv_path(results: Path) -> Path:
    return results / "sshash_density_expected_orders.csv"


def _expected_key(order: str, t) -> tuple[str, int | None]:
    return (order, int(t) if t not in (None, "") else None)


def load_expected_orders(results: Path) -> dict[tuple[str, int | None], dict]:
    path = expected_orders_csv_path(results)
    if not path.is_file():
        return {}
    out = {}
    for row in csv.DictReader(path.open()):
        label = row.get("order")
        if label:
            out[_expected_key(label, row.get("t"))] = row
    return out


def ensure_expected_orders(
    results: Path,
    m: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    mod_t_gm_file: Path | None = None,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
    mod_n_bases: int = DEFAULT_MOD_N_BASES,
) -> dict[tuple[str, int | None], dict]:
    path = expected_orders_csv_path(results)
    if not force and path.is_file():
        cached = load_expected_orders(results)
        if cached and all(int(r.get("m") or 0) == m and int(r.get("w") or 0) == w for r in cached.values()):
            return cached
    rows = compute_all_orders(
        m, w, seeds, gm_file=gm_file, mod_t_gm_file=mod_t_gm_file, n_bases=n_bases, mod_n_bases=mod_n_bases
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(EXPECTED_ORDERS_FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    return {_expected_key(r["order"], r["t"]): {k: ("" if v is None else str(v)) for k, v in r.items()} for r in rows}


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
    ap.add_argument("--gm-file", default=None)
    ap.add_argument("--mod-t-gm-file", default=None)
    ap.add_argument("--n-bases", type=int, default=DEFAULT_N_BASES)
    ap.add_argument("--mod-n-bases", type=int, default=DEFAULT_MOD_N_BASES)
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)

    w = args.k - args.m + 1
    results = Path(args.results_dir).resolve()
    gm_file = Path(args.gm_file) if args.gm_file else None
    mod_t_gm_file = Path(args.mod_t_gm_file) if args.mod_t_gm_file else None
    seeds = list(bench_seeds())[: args.n_seeds]

    cached = ensure_expected_orders(
        results, args.m, w, seeds, gm_file=gm_file, mod_t_gm_file=mod_t_gm_file,
        force=args.force, n_bases=args.n_bases, mod_n_bases=args.mod_n_bases,
    )
    print(f"wrote {expected_orders_csv_path(results)} ({len(cached)} rows)", file=sys.stderr)
    for order in ORDERS:
        for t in order_t_values(order):
            r = cached.get(_expected_key(order, t))
            if not r:
                continue
            std = r.get("density_factor_std") or None
            std_str = f" +/- {std}" if std else ""
            label = order if t is None else f"{order} (t={t})"
            print(f"  {label:36s} n_seeds={r.get('n_seeds')}  expected_density_factor={r.get('density_factor_avg')}{std_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
