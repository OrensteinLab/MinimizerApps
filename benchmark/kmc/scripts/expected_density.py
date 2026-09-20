"""Theoretical expected DNA-alphabet density_factor for every signature-order
scheme (default/lex, lexicographic/true-lex, hash, pcg, greedymini), via Monte
Carlo simulation on an i.i.d. uniform-random DNA stream -- not measured on any
real dataset. hash/pcg/greedymini are seed-dependent (see benchmark/README.md:
"Run order and the lex baseline" for why seed matters), so this reports mean
and std across the benchmark's predefined seeds, not a single value.

Feeds the `orders` xlsx sheet (see scripts/export.py); cached to
`kmc3_density_expected_orders.csv` since it's independent of any actual
benchmark run (a pure function of k/p/w and the seed list) and each
simulation costs real Monte Carlo time.
"""

from __future__ import annotations

import csv
import random
import statistics
import struct
from pathlib import Path

from scripts.config import seeds as bench_seeds
from scripts.density import _sliding_runs, density_factor
from scripts.gm import dna_gm_path

DEFAULT_N_BASES = 4_000_000
DEFAULT_SEED = 0xC0FFEE01
N_ORDER_SEEDS = 5  # matches this benchmark's default -n (see README)

ORDER_LABELS = ("default", "lexicographic", "hash", "pcg", "greedymini_dna", "greedymini_binary")


def rxs_m_xs64(x: int) -> int:
    """PCG's RXS-M-XS64 output mixing (must match custom_order.cpp exactly)."""
    x &= (1 << 64) - 1
    count = x >> 59
    x ^= x >> (5 + count)
    x = (x * 12605985483714917081) & ((1 << 64) - 1)
    return x ^ (x >> 43)


def dna_rev(mmer: int, p: int) -> int:
    """Reverse complement of a p-base DNA p-mer (must match
    common/minimizer-orders/dna_revcomp.h's dna_reverse_complement exactly --
    this file's Monte Carlo simulation is an independent from-scratch check,
    not a build of the real order, so it doesn't shell out to the C++ side)."""
    rev = 0
    shift = p * 2 - 2
    m = mmer
    for _ in range(p):
        rev += (3 - (m & 3)) << shift
        m >>= 2
        shift -= 2
    return rev


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


def build_dna_keys(
    order_label: str, k: int, p: int, *, seed: int | None = None, gm_dir: Path | None = None
) -> list[int]:
    """Build the per-DNA-p-mer comparator array (size 4^p) a real KMC run would
    use for `order_label`, matching kmc_api/mmer.h + custom_order.cpp exactly.
    `k` (the outer k-mer length) only matters for greedymini_* -- it's encoded
    in the GreedyMini table's on-disk filename (w{k-p+1}_k{p}_canonical...),
    even though the table itself, like every other order here, operates over
    p-mers (4^p entries)."""
    n = 1 << (2 * p)
    if order_label == "lexicographic":
        # True canonical lex, no is_allowed denylist: MIN(i, revcomp(i)).
        return [min(i, dna_rev(i, p)) for i in range(n)]
    if order_label == "default":
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
    if order_label == "hash":
        if seed is None:
            raise ValueError("hash order requires a seed")
        rank_by_dna = [(i ^ seed) & 0xFFFFFFFFFFFFFFFF for i in range(n)]
        return _dna_sort_key_from_rank(rank_by_dna, p, seed & 0xFFFFFFFF)
    if order_label == "pcg":
        if seed is None:
            raise ValueError("pcg order requires a seed")
        rank_by_dna = [rxs_m_xs64(i ^ seed) for i in range(n)]
        return _dna_sort_key_from_rank(rank_by_dna, p, seed & 0xFFFFFFFF)
    if order_label in ("greedymini_dna", "greedymini_binary"):
        # Both LUT backends are the same order (verified identical measured
        # density elsewhere) -- read the real DNA sort_key straight from the
        # actual _dna.gm file KMC uses, no need to rebuild it.
        path = dna_gm_path(k, p, seed, gm_dir)
        if not path.is_file():
            raise FileNotFoundError(path)
        _k, _seed, sort_key = read_dna_gm(path)
        return sort_key
    raise ValueError(f"unknown order_label: {order_label}")


def expected_density(
    dna_keys: list[int], p: int, w: int, *, n_bases: int = DEFAULT_N_BASES, seed: int = DEFAULT_SEED
) -> float | None:
    """Monte Carlo expected density_factor for a precomputed comparator on an
    i.i.d. uniform-random DNA stream (no real sequence content, no read
    boundaries -- the theoretical floor a scheme's construction targets)."""
    rng = random.Random(seed)
    mask = (1 << (2 * p)) - 1
    state = rng.getrandbits(2 * p) & mask
    keys = [dna_keys[state]]
    for _ in range(n_bases - p):
        state = ((state << 2) | rng.getrandbits(2)) & mask
        keys.append(dna_keys[state])
    n_win, n_super, _n_sel = _sliding_runs(keys, w)
    return density_factor(n_super, n_win, w)


def _order_seeds(order_label: str) -> list[int | None]:
    if order_label in ("default", "lexicographic"):
        return [None]  # deterministic, no seed dependence
    seeds = list(bench_seeds()[:N_ORDER_SEEDS])
    if order_label.startswith("greedymini"):
        # dna_gm_path's first (seeds[0]) table is the unsuffixed file
        # (w{w}_k{p}_canonical_dna.gm) -- seed=None, not seed=seeds[0] -- see
        # README: "_dna.gm embeds the tie-break seed in its header". Only
        # seeds[1:] get an explicit _s{seed} suffix on disk.
        return [None, *seeds[1:]]
    return seeds


def compute_all_orders(
    k: int, p: int, w: int, *, gm_dir: Path | None = None, n_bases: int = DEFAULT_N_BASES
) -> list[dict]:
    """Expected density_factor (mean/std over seeds where applicable) for every
    scheme. One row per order_label; density_factor_std is None for
    deterministic orders (default, lexicographic) since there's no seed to vary."""
    out = []
    for order_label in ORDER_LABELS:
        vals: list[float] = []
        for seed in _order_seeds(order_label):
            try:
                keys = build_dna_keys(order_label, k, p, seed=seed, gm_dir=gm_dir)
            except FileNotFoundError:
                continue
            v = expected_density(keys, p, w, n_bases=n_bases)
            if v is not None:
                vals.append(v)
        if not vals:
            out.append(
                {"order_label": order_label, "k": k, "p": p, "w": w, "n_seeds": 0,
                 "density_factor_avg": None, "density_factor_std": None}
            )
            continue
        avg = round(statistics.mean(vals), 4)
        std = round(statistics.pstdev(vals), 4) if len(vals) > 1 else None
        out.append(
            {"order_label": order_label, "k": k, "p": p, "w": w, "n_seeds": len(vals),
             "density_factor_avg": avg, "density_factor_std": std}
        )
    return out


EXPECTED_ORDERS_FIELDS = ("order_label", "k", "p", "w", "n_seeds", "density_factor_avg", "density_factor_std")


def expected_orders_csv_path(results: Path) -> Path:
    return results / "kmc3_density_expected_orders.csv"


def load_expected_orders(results: Path) -> dict[str, dict]:
    """Load the cache written by ensure_expected_orders; {} if not built yet."""
    path = expected_orders_csv_path(results)
    if not path.is_file():
        return {}
    out = {}
    for row in csv.DictReader(path.open()):
        label = row.get("order_label")
        if not label:
            continue
        out[label] = row
    return out


def ensure_expected_orders(
    results: Path,
    k: int,
    p: int,
    w: int,
    *,
    gm_dir: Path | None = None,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
) -> dict[str, dict]:
    """Return the cached per-order expected-density table, computing (and
    caching) it on first use. Independent of any actual benchmark run, so it's
    safe/cheap to reuse across every export once built."""
    path = expected_orders_csv_path(results)
    if not force and path.is_file():
        cached = load_expected_orders(results)
        if cached and all(int(r.get("k") or 0) == k and int(r.get("p") or 0) == p for r in cached.values()):
            return cached
    rows = compute_all_orders(k, p, w, gm_dir=gm_dir, n_bases=n_bases)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(EXPECTED_ORDERS_FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    return {r["order_label"]: {k: str(v) if v is not None else "" for k, v in r.items()} for r in rows}


def main(argv: list[str] | None = None) -> int:
    import argparse
    import sys

    from scripts.common import load_dotenv
    from scripts.export import bench_paths

    load_dotenv()
    _, default_results = bench_paths()

    ap = argparse.ArgumentParser(
        prog="scripts expected-density",
        description="Theoretical expected density_factor per signature-order scheme "
        "(Monte Carlo on i.i.d. random DNA; mean/std over predefined seeds for "
        "hash/pcg/greedymini). Feeds the `orders` xlsx sheet's expected_density_avg/std "
        "columns -- run this once (or after changing k/p) to (re)build that cache.",
    )
    ap.add_argument("--k", type=int, default=28)
    ap.add_argument("-p", type=int, default=9)
    ap.add_argument("--results-dir", default=str(default_results))
    ap.add_argument("--gm-dir", default=None)
    ap.add_argument("--n-bases", type=int, default=DEFAULT_N_BASES)
    ap.add_argument("--force", action="store_true", help="recompute even if a matching cache exists")
    args = ap.parse_args(argv)

    w = args.k - args.p + 1
    results = Path(args.results_dir).resolve()
    gm_dir = Path(args.gm_dir) if args.gm_dir else None

    cached = ensure_expected_orders(results, args.k, args.p, w, gm_dir=gm_dir, force=args.force)
    print(f"wrote {expected_orders_csv_path(results)} ({len(cached)} orders)", file=sys.stderr)
    for label in ORDER_LABELS:
        r = cached.get(label)
        if not r:
            continue
        std = r.get("density_factor_std") or None
        std_str = f" +/- {std}" if std else ""
        print(f"  {label:18s} n_seeds={r.get('n_seeds')}  expected_density={r.get('density_factor_avg')}{std_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
