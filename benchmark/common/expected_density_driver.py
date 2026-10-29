"""Shared Monte Carlo driver for the theoretical expected density_factor of a
scheme on an i.i.d. uniform-random DNA stream (kmc, sshash, minimap2; kraken2
keeps its own C++ tool). Each app's ``expected_density.py`` is a thin adapter
that supplies only:

  (a) the variants to compute for its config,
  (b) ``sample(variant, ctx, seed) -> float | None``: one Monte Carlo draw
      (builds the app's rank/hash function and applies its windowing rule;
      the helpers `random_dna`, `canonical_keys`, `Rounding.factor` cover the
      common stream -> canonical hash -> sliding-window skeleton),
  (c) ``row_id(variant) -> dict`` (the id/parameter columns of the cache row),
  (d) the cache file path and column tuple.

This module owns the loop over variants x seeds, mean/std and its rounding,
the cache read/validate/merge/write (`ensure_cache`) and the shared CLI
arguments. The cache layout (a CSV with the id columns, k, w, L, n_seeds,
density_factor_mean, density_factor_std [, cache_key]) is whatever the adapter
declares in its `fields`; exporters read it through the adapters'
`load_expected_schemes` / `ensure_expected_schemes`.

(Not named `expected_density.py`: minimap2's adapter is imported as the
top-level module `expected_density` and would be shadowed.)
"""

from __future__ import annotations

import csv
import random
import statistics
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from density import density_factor, sliding_runs

DEFAULT_SEED = 0xC0FFEE01  # the Monte Carlo RNG seed (not an order seed)
STAT_FIELDS = ("n_seeds", "density_factor_mean", "density_factor_std")


class Unavailable(Exception):
    """Raised by `prepare`/`sample` when a variant cannot be computed (e.g. its
    GreedyMini table is missing): it contributes no values (n_seeds 0)."""


@dataclass(frozen=True)
class Rounding:
    """Per-app precision, kept as each app's established one so caches stay
    bit-identical: `sample` = decimals of one draw's density_factor, `avg`/`std`
    = decimals of the mean / population std over seeds, `single_std` = std
    reported when only one seed produced a value."""

    sample: int
    avg: int
    std: int
    single_std: float | None

    def factor(self, n_selected: int, n_kmers: int, w: int) -> float | None:
        return density_factor(n_selected, n_kmers, w, self.sample)

    def aggregate(self, vals: Sequence[float]) -> dict:
        out = {"n_seeds": len(vals), "density_factor_mean": None, "density_factor_std": None}
        if vals:
            out["density_factor_mean"] = round(statistics.mean(vals), self.avg)
            out["density_factor_std"] = (
                round(statistics.pstdev(vals), self.std) if len(vals) > 1 else self.single_std
            )
        return out


# ---------------------------------------------------------------------------
# Stream -> keys -> density skeleton
# ---------------------------------------------------------------------------


def random_dna(n_bases: int, rng_seed: int = DEFAULT_SEED) -> list[int]:
    """i.i.d. uniform bases 0..3 (`random.Random(rng_seed).randrange(4)`)."""
    rng = random.Random(rng_seed)
    return [rng.randrange(4) for _ in range(n_bases)]


def canonical_keys(stream: list[int], length: int, hash_fn: Callable[[int, int], int],
                   canonical: Callable[[int, int], int], seed: int) -> list[int]:
    """One `hash_fn(canonical(l-mer), seed)` per l-mer of `stream`, l-mers packed
    first base most significant (the app's `canonical` knows its base encoding)."""
    mask = (1 << (2 * length)) - 1
    state = 0
    for i in range(length):
        state = (state << 2) | stream[i]
    keys = [hash_fn(canonical(state, length), seed)]
    for i in range(length, len(stream)):
        state = ((state << 2) | stream[i]) & mask
        keys.append(hash_fn(canonical(state, length), seed))
    return keys


def window_density(keys: list[int], w: int, rounding: Rounding) -> float | None:
    """density_factor of a forward scheme: super-k-mer boundaries (changes of the
    window's minimal key) per window over `keys`."""
    n_win, n_super, _n_sel = sliding_runs(keys, w)
    return rounding.factor(n_super, n_win, w)


# ---------------------------------------------------------------------------
# Loop over variants x seeds
# ---------------------------------------------------------------------------


def compute_rows(
    variants: Iterable,
    *,
    seeds_of: Callable[[object], Sequence],
    sample: Callable[[object, object, object], float | None],
    row_id: Callable[[object], dict],
    rounding: Rounding,
    prepare: Callable[[object], object] | None = None,
) -> list[dict]:
    """One cache row per variant: `row_id(v)` + n_seeds/mean/std of
    `sample(v, ctx, seed)` over `seeds_of(v)`; `ctx = prepare(v)` once per
    variant (e.g. the loaded GreedyMini table)."""
    out = []
    for v in variants:
        vals: list[float] = []
        try:
            ctx = prepare(v) if prepare else None
            for seed in seeds_of(v):
                try:
                    x = sample(v, ctx, seed)
                except Unavailable:
                    continue
                if x is not None:
                    vals.append(x)
        except Unavailable:
            pass
        out.append({**row_id(v), **rounding.aggregate(vals)})
    return out


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------


def read_csv_rows(path: Path, key_of_row: Callable[[dict], object | None],
                  normalize: Callable[[dict], dict] | None = None) -> dict:
    """{key: row} of the cache (rows whose key is None are ignored); {} if absent."""
    if not path.is_file():
        return {}
    out = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            if normalize:
                row = normalize(row)
            key = key_of_row(row)
            if key is not None:
                out[key] = row
    return out


def write_csv_rows(path: Path, fields: Sequence[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fields))
        writer.writeheader()
        for r in rows:
            writer.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in fields})


def stringify(rows: Iterable[dict], key_of_row: Callable[[dict], object]) -> dict:
    return {key_of_row(r): {k: ("" if v is None else str(v)) for k, v in r.items()} for r in rows}


def ensure_cache(
    path: Path,
    fields: Sequence[str],
    wanted: Sequence,
    *,
    key_of_variant: Callable[[object], object],
    key_of_row: Callable[[dict], object | None],
    load: Callable[[], dict],
    compute: Callable[[list], list[dict]],
    accept: Callable[[dict], bool] = lambda cached: True,
    merge: bool = False,
    stamp: dict | None = None,
    force: bool = False,
) -> dict:
    """The cached {key: row}, computing (and rewriting the CSV) on a miss.
    A cache hits when `load()` is non-empty, `accept(cached)` holds and it covers
    `wanted`: exactly (default; any other variant set is a stale cache) or, with
    `merge`, as a superset -- then a partial miss keeps the cached rows and
    computes only the missing variants. `stamp` (e.g. {"cache_key": ...}) is
    added to every computed row."""
    wanted_keys = {key_of_variant(v) for v in wanted}
    kept: list[dict] = []
    if not force and path.is_file():
        cached = load()
        if cached and accept(cached):
            if merge:
                if wanted_keys <= set(cached):
                    return cached
                kept = [dict(r) for r in cached.values()]
            elif wanted_keys == set(cached):
                return cached
    have = {key_of_row(r) for r in kept}
    todo = [v for v in wanted if key_of_variant(v) not in have]
    new = compute(todo)
    if stamp:
        for r in new:
            r.update(stamp)
    rows = kept + new
    write_csv_rows(path, fields, rows)
    return stringify(rows, key_of_row)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def base_parser(prog: str | None, description: str, *, default_results: str | None, default_n_bases: int):
    """argparse parser with the arguments every adapter's CLI shares."""
    import argparse

    ap = argparse.ArgumentParser(prog=prog, description=description)
    ap.add_argument("--results-dir", default=default_results, required=default_results is None)
    ap.add_argument("--n-bases", type=int, default=default_n_bases)
    ap.add_argument("--force", action="store_true", help="recompute even if a matching cache exists")
    return ap


def format_row(label: str, row: dict, width: int = 36, name: str = "expected_density_factor") -> str:
    std = row.get("density_factor_std") or None
    std_str = f" +/- {std}" if std else ""
    return f"  {label:{width}s} n_seeds={row.get('n_seeds')}  {name}={row.get('density_factor_mean')}{std_str}"
