"""Theoretical expected density_factor for every sshash order (native,
xor-hash, xor-pcg-hash, greedymini), via Monte Carlo simulation on an i.i.d.
uniform-random DNA stream -- not measured on any real dataset. Mirrors
benchmark/kmc/scripts/expected_density.py's approach and naming, ported to
sshash's own order formulas (see scripts/density.py).

xor-hash/xor-pcg-hash/greedymini are seed-dependent, so this reports mean
and std across the benchmark's predefined seeds (see config.yaml), not a
single value; native's own seed only changes an internal mixing constant, not
its statistical behavior as a random order (see density.py:native_hash's
docstring), but it's still varied across seeds here for a fair apples-to-
apples comparison.

Feeds the `orders` xlsx sheet's expected_density_avg/std columns (see
scripts/export.py); cached to sshash_density_expected_orders.csv since it's
independent of any actual benchmark run (a pure function of k/m/w + the seed
list) and each simulation costs real Monte Carlo time.
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
    native_hash,
    read_gm_ranks,
    xor_hash,
    xor_pcg_hash,
    _sliding_runs,
)

DEFAULT_N_BASES = 1_000_000
DEFAULT_SEED = 0xC0FFEE01

ORDER_LABELS = ("native", "xor-hash", "xor-pcg-hash", "greedymini")


def _keys_for_order(order: str, dna_stream: list[int], m: int, seed: int, gm_ranks: list[int] | None) -> list[int]:
    """One hash value per m-mer position of a precomputed 2-bit/base DNA
    stream (`dna_stream[i]` = the i-th base, 0..3): canonicalize each
    m-window then apply the order's hash, exactly matching what
    util::canonical_mmer_at + hasher_type::hash do for a real sshash build."""
    n = len(dna_stream)
    mask = (1 << (2 * m)) - 1
    state = 0
    for i in range(m):
        state = (state << 2) | dna_stream[i]
    keys: list[int] = []

    def hash_value(value: int) -> int:
        if order == "native":
            return native_hash(value, seed)
        if order == "xor-hash":
            return xor_hash(value, seed)
        if order == "xor-pcg-hash":
            return xor_pcg_hash(value, seed)
        if order == "greedymini":
            assert gm_ranks is not None
            return greedymini_hash(value, m, gm_ranks, seed)
        raise ValueError(f"unknown order: {order}")

    keys.append(hash_value(canonical_mmer(state, m)))
    for i in range(m, n):
        state = ((state << 2) | dna_stream[i]) & mask
        keys.append(hash_value(canonical_mmer(state, m)))
    return keys


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
    dna_stream = [rng.randrange(4) for _ in range(n_bases)]
    keys = _keys_for_order(order, dna_stream, m, seed, gm_ranks)
    n_win, n_super, _n_sel = _sliding_runs(keys, w)
    return density_factor(n_super, n_win, w)


def compute_all_orders(
    m: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    n_bases: int = DEFAULT_N_BASES,
) -> list[dict]:
    gm_ranks = read_gm_ranks(gm_file) if gm_file and gm_file.is_file() else None
    out = []
    for order in ORDER_LABELS:
        if order == "greedymini" and gm_ranks is None:
            out.append(
                {"order": order, "m": m, "w": w, "n_seeds": 0,
                 "density_factor_avg": None, "density_factor_std": None}
            )
            continue
        vals: list[float] = []
        for seed in seeds:
            v = expected_density_factor(order, m, w, seed=seed, n_bases=n_bases, gm_ranks=gm_ranks)
            if v is not None:
                vals.append(v)
        if not vals:
            out.append(
                {"order": order, "m": m, "w": w, "n_seeds": 0,
                 "density_factor_avg": None, "density_factor_std": None}
            )
            continue
        avg = round(statistics.mean(vals), 4)
        std = round(statistics.pstdev(vals), 4) if len(vals) > 1 else 0.0
        out.append(
            {"order": order, "m": m, "w": w, "n_seeds": len(vals),
             "density_factor_avg": avg, "density_factor_std": std}
        )
    return out


EXPECTED_ORDERS_FIELDS = ("order", "m", "w", "n_seeds", "density_factor_avg", "density_factor_std")


def expected_orders_csv_path(results: Path) -> Path:
    return results / "sshash_density_expected_orders.csv"


def load_expected_orders(results: Path) -> dict[str, dict]:
    path = expected_orders_csv_path(results)
    if not path.is_file():
        return {}
    out = {}
    for row in csv.DictReader(path.open()):
        label = row.get("order")
        if label:
            out[label] = row
    return out


def ensure_expected_orders(
    results: Path,
    m: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
) -> dict[str, dict]:
    path = expected_orders_csv_path(results)
    if not force and path.is_file():
        cached = load_expected_orders(results)
        if cached and all(int(r.get("m") or 0) == m and int(r.get("w") or 0) == w for r in cached.values()):
            return cached
    rows = compute_all_orders(m, w, seeds, gm_file=gm_file, n_bases=n_bases)
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
    ap.add_argument("--m", type=int, default=11)
    ap.add_argument("--k", type=int, default=31, help="only used to derive w=k-m+1")
    ap.add_argument("--results-dir", default=str(default_results))
    ap.add_argument("--gm-file", default=None)
    ap.add_argument("--n-bases", type=int, default=DEFAULT_N_BASES)
    ap.add_argument("--n-seeds", type=int, default=5)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)

    w = args.k - args.m + 1
    results = Path(args.results_dir).resolve()
    gm_file = Path(args.gm_file) if args.gm_file else None
    seeds = list(bench_seeds())[: args.n_seeds]

    cached = ensure_expected_orders(results, args.m, w, seeds, gm_file=gm_file, force=args.force, n_bases=args.n_bases)
    print(f"wrote {expected_orders_csv_path(results)} ({len(cached)} orders)", file=sys.stderr)
    for label in ORDER_LABELS:
        r = cached.get(label)
        if not r:
            continue
        std = r.get("density_factor_std") or None
        std_str = f" +/- {std}" if std else ""
        print(f"  {label:14s} n_seeds={r.get('n_seeds')}  expected_density_factor={r.get('density_factor_avg')}{std_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
