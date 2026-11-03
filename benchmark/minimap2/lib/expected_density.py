#!/usr/bin/env python3
"""Theoretical expected density_factor for every minimap2 scheme (native,
xor, xor-pcg-hash, greedymini and mod-minimizer-greedymini in each projection,
mod-minimizer-xor), counted on uniformly random DNA by minimap2's own
mm_sketch / mm_sketch_scheme (the minimap2 tool of benchmark/common/expected_density/),
homopolymer-compressed as -H does -- not measured on any real dataset. Thin
adapter over benchmark/common/expected_density_driver.py (tool, loop, mean/std,
cache, CLI): this file only has minimap2's variants and their tool options.

xor/xor-pcg-hash/mod-minimizer-xor are seed-dependent, so this reports mean and std across the benchmark's
predefined seeds (see config.yaml), not a single value; native has no seed
at all in minimap2 (mm_hash64 takes none), so its std over the seeds' random
strings is that of the string alone.

Feeds the `schemes` xlsx sheet's expected_density_factor_mean/std columns
(see lib/export.py); cached to minimap2_density_expected_k<k>_w<w>_t<mod_t>_<hpc|nohpc>.csv since
it's independent of any actual benchmark run (a pure function of k/w/t/HPC,
the GreedyMini table, the seed list and the string length) and each
count costs real time. Every cached row carries a
cache_key over everything but k/w/t/HPC (which are in the file name): the
GreedyMini tables' paths and content hashes (greedymini's table and
mod-minimizer-greedymini's, which applies only where it fits t = gm_k + 2i), the seeds
and n_bases, so changing any of them recomputes the cache.
(Both GreedyMini projections get their own rows: each row carries its
projection, gm_k and gm_table columns.) A cache
written before cache_key existed simply misses once.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import expected_density_driver as shared  # noqa: E402
from scheme_names import MOD_SCHEMES, gm_fields, gm_table_fits, is_greedymini  # noqa: E402

DEFAULT_N_BASES = shared.DEFAULT_N_BASES
# density_factor of one seed is rounded to 8 decimals, the mean to 6 and the std to 8
ROUNDING = shared.Rounding(sample=8, avg=6, std=8, single_std=0.0)

# (scheme, projection) variants: greedymini once per projection, the others
# have no projection (None).
VARIANTS = [
    ("native", None), ("xor", None), ("xor-pcg-hash", None),
    ("greedymini", "at0-cg1"), ("greedymini", "at1-cg0"),
    ("mod-minimizer-xor", None),
    ("mod-minimizer-greedymini", "at0-cg1"), ("mod-minimizer-greedymini", "at1-cg0"),
]

def variants_for(k: int, mod_t: int, gm_file: Path | None, mod_gm_file: Path | None) -> list[tuple[str, str | None]]:
    """VARIANTS that apply at this (k, t), as lib/run.sh runs them: a scheme with a GreedyMini
    table only where its table fits what it ranks (scheme_names.gm_table_fits: greedymini's
    table the k-mer, mod-minimizer-greedymini's the t-mer; w20_k9 at k=19, w=10, t=9). A scheme
    whose table is not given stays (its rows get n_seeds 0)."""
    files = {"greedymini": gm_file, "mod-minimizer-greedymini": mod_gm_file}

    def applies(scheme: str) -> bool:
        f = files.get(scheme)
        return f is None or gm_table_fits(scheme, gm_fields(table=f.name)["gm_k"], k, mod_t)
    return [v for v in VARIANTS if applies(v[0])]


def compute_all_schemes(
    k: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    mod_gm_file: Path | None = None,
    n_bases: int = DEFAULT_N_BASES,
    mod_t: int = 5,
    is_hpc: bool = True,
) -> list[dict]:
    files = {"greedymini": gm_file, "mod-minimizer-greedymini": mod_gm_file}
    info = {s: gm_fields(table=f.name) if f else {"gm_k": None, "gm_table": None} for s, f in files.items()}

    def row_id(v):
        scheme, projection = v
        gm = info[scheme] if is_greedymini(scheme) else {"gm_k": None, "gm_table": None}
        return {
            "scheme": scheme, "gm_k": gm["gm_k"], "projection": projection,
            "gm_table": gm["gm_table"], "mod_t": mod_t if scheme in MOD_SCHEMES else None,
            "k": k, "w": w, "L": k + w - 1,
        }

    def options(v):
        scheme, projection = v
        opts = {"k": k, "scheme": scheme}
        if is_hpc:
            opts["hpc"] = None
        if scheme in MOD_SCHEMES:
            opts["t"] = mod_t
        if is_greedymini(scheme):
            f = files[scheme]
            if not (f and f.is_file()):
                raise shared.Unavailable  # no table: n_seeds 0
            opts.update({"gm": f, "projection": projection})
        return opts

    return shared.compute_rows(
        variants_for(k, mod_t, gm_file, mod_gm_file),
        seeds_of=lambda v: seeds,
        prepare=options,
        sample=lambda v, opts, seed: shared.tool_density_factor("minimap2", opts, seed, n_bases, w, ROUNDING),
        row_id=row_id,
        rounding=ROUNDING,
    )


EXPECTED_SCHEMES_FIELDS = ("scheme", "gm_k", "projection", "gm_table", "mod_t", "k", "w", "L", *shared.STAT_FIELDS, "cache_key")


def expected_cache_key(
    seeds: list[int],
    *,
    gm_file: Path | None,
    mod_gm_file: Path | None = None,
    n_bases: int,
) -> str:
    """Identity of every cache input not already in the file name (k, w, t,
    HPC): the GreedyMini tables (resolved path + content sha256; 'missing' when
    absent, so building a table later invalidates a cache computed without
    it), the seed list, the string length, and a version of the order hashes
    and of the counting."""
    def table_id(f: Path | None) -> str:
        if f is None:
            return "none"
        path = Path(f).resolve()
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "missing"
        return f"{path}:{digest}"

    ident = {
        # bump when an order's hash or the row layout changes (2: greedymini --order-seed
        # toggles Gamma; 3: projection/gm_k/gm_table columns; 4: mod_t and L columns;
        # 5: the mod-minimizer's centre tie rule; 6: mod-minimizer-greedymini;
        # 7: counted by minimap2's own sketching, a random string per seed;
        # 8: mean/std to 6/8 decimals)
        "hash_version": 8,
        "gm": table_id(gm_file),
        "mod_gm": table_id(mod_gm_file),
        "seeds": [hex(s) for s in seeds],
        "n_bases": n_bases,
    }
    return hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]


def expected_schemes_csv_path(results: Path, k: int, w: int, t: int, is_hpc: bool) -> Path:
    return results / "csvs" / f"minimap2_density_expected_k{k}_w{w}_t{t}_{'hpc' if is_hpc else 'nohpc'}.csv"


def variant_key(row: dict) -> tuple[str, str]:
    """(scheme, projection or '') -- the key of one expected-density row."""
    return (row["scheme"], row.get("projection") or "")


def _row_key(row: dict) -> tuple[str, str] | None:
    return variant_key(row) if row.get("scheme") else None


def load_expected_schemes(results: Path, k: int, w: int, t: int, is_hpc: bool) -> dict[tuple[str, str], dict]:
    return shared.read_csv_rows(expected_schemes_csv_path(results, k, w, t, is_hpc), _row_key)


def ensure_expected_schemes(
    results: Path,
    k: int,
    w: int,
    seeds: list[int],
    *,
    gm_file: Path | None = None,
    mod_gm_file: Path | None = None,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
    mod_t: int = 5,
    is_hpc: bool = True,
) -> dict[tuple[str, str], dict]:
    """Keyed by (scheme, projection or ''). is_hpc must match the dataset's preset (lib/constants.sh:
    dataset_preset; map-pb enables -H/HPC k-mers, map-ont does not) for
    'particular' (measured) and 'expected' (this function) density_factor
    to be a fair comparison: map-pb datasets use the HPC value, map-ont ones
    the non-HPC one; each is cached in its own file."""
    variants = variants_for(k, mod_t, gm_file, mod_gm_file)
    # a table is part of the key only where a variant uses it (the k=15 cache ignores w20_k9)
    used = {v[0] for v in variants}
    key = expected_cache_key(seeds, gm_file=gm_file if "greedymini" in used else None,
                             mod_gm_file=mod_gm_file if "mod-minimizer-greedymini" in used else None,
                             n_bases=n_bases)
    return shared.ensure_cache(
        expected_schemes_csv_path(results, k, w, mod_t, is_hpc),
        EXPECTED_SCHEMES_FIELDS,
        variants,
        key_of_variant=lambda v: (v[0], v[1] or ""),
        key_of_row=_row_key,
        load=lambda: load_expected_schemes(results, k, w, mod_t, is_hpc),
        compute=lambda _todo: compute_all_schemes(k, w, seeds, gm_file=gm_file, mod_gm_file=mod_gm_file, n_bases=n_bases,
                                                 mod_t=mod_t, is_hpc=is_hpc),
        accept=lambda cached: all(r.get("cache_key") == key for r in cached.values()),
        stamp={"cache_key": key},
        force=force,
    )


def main(argv: list[str] | None = None) -> int:
    ap = shared.base_parser(
        None,
        "Theoretical expected density_factor per minimap2 scheme (minimap2's own "
        "sketching on uniformly random DNA; mean/std over the benchmark's predefined seeds). Feeds the "
        "`schemes` xlsx sheet's expected_density_factor_mean/std columns.",
        default_results=None,
        default_n_bases=DEFAULT_N_BASES,
    )
    ap.add_argument("--k", type=int, default=19, help="19 = map-pb's k; 15 for map-ont")
    ap.add_argument("--w", type=int, default=10)
    ap.add_argument("--mod-t", type=int, default=None, help="default 4 + (k-4) %% w")
    ap.add_argument("--no-hpc", action="store_true", help="raw k-mers (map-ont) instead of HPC (map-pb)")
    ap.add_argument("--gm-file", default=None, help="greedymini's table of k-mers")
    ap.add_argument("--mod-gm-file", default=None, help="mod-minimizer-greedymini's table (gm_k-mers, t = gm_k + 2i)")
    ap.add_argument("--seeds", default="0x2545F4914F6CDD1D,0x9E3779B97F4A7C15,0xDEADBEEFCAFEBABE,0x1234567890ABCDEF,0xA5A5A5A5A5A5A5A5")
    args = ap.parse_args(argv)

    seeds = [int(s, 16) if s.lower().startswith("0x") else int(s) for s in args.seeds.split(",")]
    rows = ensure_expected_schemes(
        Path(args.results_dir), args.k, args.w, seeds,
        gm_file=Path(args.gm_file) if args.gm_file else None,
        mod_gm_file=Path(args.mod_gm_file) if args.mod_gm_file else None,
        force=args.force, n_bases=args.n_bases,
        mod_t=args.mod_t if args.mod_t is not None else 4 + (args.k - 4) % args.w,
        is_hpc=not args.no_hpc,
    )
    for key, row in rows.items():
        print(key, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
