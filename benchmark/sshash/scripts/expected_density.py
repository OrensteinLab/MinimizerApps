"""Theoretical expected density_factor for every sshash scheme variant (native,
xor, xor-pcg-hash, greedymini and mod-minimizer-* per table/projection),
counted on uniformly random DNA by sshash's own minimizer iterators (the sshash
tool of benchmark/common/expected_density/) -- not measured on any real dataset.
Thin adapter over benchmark/common/expected_density_driver.py (tool, loop,
mean/std, cache, CLI): this file only has sshash's variants and their tool options.

xor/xor-pcg-hash/greedymini are seed-dependent, so this reports mean and
std across the benchmark's predefined seeds (see config.yaml), not a single
value; native's own seed (sshash -s) only changes an internal mixing constant,
but it's still varied across seeds here, as in the runs.

Feeds the `schemes` xlsx sheet's expected_density_factor_mean/std columns (see
scripts/export.py); cached to sshash_density_expected_schemes.csv since it's
independent of any actual benchmark run (a pure function of m/w + the seed
list) and each count costs real time.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import expected_density_driver as shared
from scripts.density import gm_centre_valid
import scheme_names as sn
from scripts.config import gm_table_file
from scripts.schemes import SCHEMES, expand_variants, mod_t, scheme_dir_name, variant_key

DEFAULT_N_BASES = shared.DEFAULT_N_BASES
ROUNDING = shared.Rounding(sample=10, avg=8, std=10, single_std=0.0)


def compute_all_schemes(
    m: int,
    w: int,
    seeds: list[int],
    *,
    n_bases: int = DEFAULT_N_BASES,
    variants: Iterable[dict] | None = None,
) -> list[dict]:
    """Expected-density rows for `variants` (id dicts as scripts.schemes.expand_variants
    returns; default: every variant config.yaml registers at k = w + m - 1)."""
    if variants is None:
        variants = expand_variants(SCHEMES, w + m - 1, m)

    def options(v):
        scheme = v["scheme"]
        opts = {"k": w + m - 1, "m": m, "scheme": scheme}
        if sn.is_mod_minimizer(scheme):
            opts["t"] = v.get("mod_t") or mod_t(scheme, w + m - 1, m)
        if sn.is_greedymini(scheme):
            path = v.get("gm_file") or gm_table_file(v["gm_table"])
            if not Path(path).is_file():
                raise shared.Unavailable  # no table: n_seeds 0
            if scheme == "greedymini" and not gm_centre_valid(m, int(v["gm_k"])):
                raise shared.Unavailable  # needs m - gm_k even (sshash rejects it too)
            opts.update({"gm": path, "projection": sn.projection_from_flag(v.get("projection"))})
        return opts

    return shared.compute_rows(
        variants,
        seeds_of=lambda v: seeds,
        prepare=options,
        sample=lambda v, opts, seed: shared.tool_density_factor("sshash", opts, seed, n_bases, w, ROUNDING),
        row_id=lambda v: {**{f: v.get(f) for f in sn.ID_FIELDS}, "k": m, "w": w, "L": m + w - 1},
        rounding=ROUNDING,
    )


EXPECTED_SCHEMES_FIELDS = (*sn.ID_FIELDS, "k", "w", "L", *shared.STAT_FIELDS, "mod_tie_rule", "n_bases")
# The mod-minimizer tie rule the cached rows were computed with; a cache without it
# (the retired G/T-strand rule) is recomputed.
MOD_TIE_RULE = "centre"


def expected_schemes_csv_path(results: Path) -> Path:
    return results / "csvs" / "sshash_density_expected_schemes.csv"


def _row_key(row: dict) -> tuple | None:
    """variant_key of a CSV row; None for a stale row (an old scheme spelling)."""
    try:
        return variant_key(sn.parse_id(row, baselines=("native",)))
    except (ValueError, KeyError):
        return None


def load_expected_schemes(results: Path) -> dict[tuple, dict]:
    """The cached rows by variant_key (rows with an unknown scheme spelling, e.g. from
    before the canonical names, are ignored)."""
    return shared.read_csv_rows(expected_schemes_csv_path(results), _row_key)


def ensure_expected_schemes(
    results: Path,
    m: int,
    w: int,
    seeds: list[int],
    *,
    force: bool = False,
    n_bases: int = DEFAULT_N_BASES,
    variants: Iterable[dict] | None = None,
) -> dict[tuple, dict]:
    """Cached expected density per variant (keyed by variant_key). With `variants`
    given, only those are computed, merged into
    whatever the cache already holds for this m/w."""
    wanted = list(variants) if variants is not None else expand_variants(SCHEMES, w + m - 1, m)
    return shared.ensure_cache(
        expected_schemes_csv_path(results),
        EXPECTED_SCHEMES_FIELDS,
        wanted,
        key_of_variant=variant_key,
        key_of_row=_row_key,
        load=lambda: load_expected_schemes(results),
        compute=lambda todo: compute_all_schemes(m, w, seeds, n_bases=n_bases, variants=todo),
        accept=lambda cached: all(
            int(r.get("k") or r.get("m") or 0) == m and int(r.get("w") or 0) == w
            and r.get("mod_tie_rule") == MOD_TIE_RULE and str(r.get("n_bases")) == str(n_bases)
            for r in cached.values()),
        merge=True,
        stamp={"mod_tie_rule": MOD_TIE_RULE, "n_bases": n_bases},
        force=force,
    )


def main(argv: list[str] | None = None) -> int:
    import sys

    from scripts.common import bench_paths, load_dotenv
    from scripts.config import seeds as bench_seeds

    load_dotenv()
    _, default_results = bench_paths()

    ap = shared.base_parser(
        "scripts expected-density",
        "Theoretical expected density_factor per sshash scheme (sshash's own "
        "minimizer iterators on uniformly random DNA; mean/std over the benchmark's predefined seeds). Feeds the "
        "`schemes` xlsx sheet's expected_density_factor_mean/std columns.",
        default_results=str(default_results),
        default_n_bases=DEFAULT_N_BASES,
    )
    ap.add_argument("--m", type=int, default=21)
    ap.add_argument("--k", type=int, default=31, help="only used to derive w=k-m+1")
    ap.add_argument("--n-seeds", type=int, default=5)
    args = ap.parse_args(argv)

    w = args.k - args.m + 1
    results = Path(args.results_dir).resolve()
    seeds = list(bench_seeds())[: args.n_seeds]

    cached = ensure_expected_schemes(results, args.m, w, seeds, force=args.force, n_bases=args.n_bases)
    print(f"wrote {expected_schemes_csv_path(results)} ({len(cached)} rows)", file=sys.stderr)
    for v in expand_variants(SCHEMES, w + args.m - 1, args.m):
        r = cached.get(variant_key(v))
        if r:
            print(shared.format_row(scheme_dir_name(v["scheme"], v["gm_k"], v["projection"], v["mod_t"]), r))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
