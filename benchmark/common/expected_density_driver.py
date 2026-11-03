"""Shared driver for the expected density_factor of a scheme (kmc, sshash,
minimap2, kraken2). The density is counted by one C++ tool per app
(benchmark/common/expected_density/, see its driver.h): on DEFAULT_N_BASES
uniformly random bases, a string of its own per seed, with the app's own
selection code, so that the tie-breaking rule and what is counted (the app's
particular density) are the app's. Each app's ``expected_density.py`` is a thin
adapter that supplies only:

  (a) the variants to compute for its config,
  (b) ``sample(variant, ctx, seed) -> float | None``: the tool's density_factor
      of the variant for one seed (`tool_density_factor`),
  (c) ``row_id(variant) -> dict`` (the id/parameter columns of the cache row),
  (d) the cache file path and column tuple.

This module owns the tools, the loop over variants x seeds, mean/std and its
rounding, the cache read/validate/merge/write (`ensure_cache`) and the shared
CLI arguments. The cache layout (a CSV with the id columns, k, w, L, n_seeds,
density_factor_mean, density_factor_std, n_bases [, cache_key]) is whatever
the adapter declares in its `fields`; exporters read it through the adapters'
`load_expected_schemes` / `ensure_expected_schemes`.

(Not named `expected_density.py`: minimap2's adapter is imported as the
top-level module `expected_density` and would be shadowed.)
"""

from __future__ import annotations

import csv
import functools
import os
import statistics
import subprocess
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from density import density_factor

DEFAULT_N_BASES = 1_000_000_000  # random bases per seed, the same for every app and scheme
STAT_FIELDS = ("n_seeds", "density_factor_mean", "density_factor_std")
TOOL_DIR = Path(__file__).resolve().parent / "expected_density"


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
# The C++ tools
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=None)
def tool_path(app: str) -> Path:
    """The app's expected-density tool, brought up to date by `make` in TOOL_DIR once per
    process, so a change to the app's code rebuilds it (SSH_EXT in the environment
    points sshash's build at another pthash checkout)."""
    target = f"build/{app}_expected_density"
    subprocess.run(["make", "-C", str(TOOL_DIR), target], check=True, stdout=subprocess.DEVNULL)
    return TOOL_DIR / target


def tool_density_factor(app: str, options: dict, seed: int, n_bases: int, w: int, rounding: Rounding) -> float | None:
    """density_factor of one variant and seed: `options` are the tool's own
    (--name value; a None value passes the flag alone). Threads: the
    EXPECTED_DENSITY_THREADS environment variable (default 16)."""
    cmd = [str(tool_path(app)), "--w", str(w), "--seeds", str(seed), "--n-bases", str(n_bases),
           "--threads", os.environ.get("EXPECTED_DENSITY_THREADS", "16")]
    for name, value in options.items():
        cmd += [f"--{name}"] + ([] if value is None else [str(value)])
    run = subprocess.run(cmd, capture_output=True, text=True)
    if run.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}: {run.stderr.strip()}")
    out = run.stdout.splitlines()
    row = dict(zip(out[0].split(","), out[-1].split(",")))
    return rounding.factor(int(row["n_selected"]), int(row["n_kmers"]), w)


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
