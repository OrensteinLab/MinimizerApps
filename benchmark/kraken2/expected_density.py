#!/usr/bin/env python3
"""Expected density_factor of the kraken2 runs, counted on uniformly random DNA
by Kraken2's own MinimizerScanner (the kraken2 tool of
benchmark/common/expected_density/, which counts as benchmark_minimizers counts
the particular density). The density microbenchmark (minimizer_density.csv)
also estimates it, on its own shorter string; `with_expected` replaces that
estimate by the tool's, so every app's expected density is counted the same way.

One value per density row, i.e. per (scheme, gm table, projection, mod_t) and
seed (the toggle mask, toggle_hex), cached in
<results>/csvs/density_expected_kraken2.csv.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import expected_density_driver as shared  # noqa: E402

GM_DIR = Path(__file__).resolve().parents[2] / "orders" / "GreedyMini"
ROUNDING = shared.Rounding(sample=6, avg=6, std=6, single_std=0.0)
FIELDS = ("scheme", "gm_k", "projection", "gm_table", "mod_t", "toggle_hex", "k", "w", "L", "spaces",
          "n_bases", "density_factor")
KEY = ("scheme", "projection", "gm_table", "mod_t", "toggle_hex", "k", "L", "spaces")


def csv_path(results: Path) -> Path:
    return Path(results) / "csvs" / "density_expected_kraken2.csv"


def _key(row: dict) -> tuple:
    return tuple(str(row.get(f) or "") for f in KEY)


def density_factor(row: dict, n_bases: int) -> float | None:
    """The tool's density_factor for one density row (k = minimizer length, L = k-mer length)."""
    k, L = int(row["k"]), int(row["L"])
    options = {"l": k, "k": L, "spaces": row.get("spaces") or 0, "scheme": row["scheme"]}
    if row.get("gm_table"):
        options.update({"gm": GM_DIR / f"{row['gm_table']}.gm", "projection": row.get("projection") or "at0-cg1"})
    return shared.tool_density_factor("kraken2", options, int(row["toggle_hex"], 16), n_bases, L - k + 1, ROUNDING)


def with_expected(dens_rows: list[dict], results: Path, n_bases: int = shared.DEFAULT_N_BASES) -> list[dict]:
    """dens_rows with fac_expected (and dens_expected) replaced by the tool's
    values, computing and caching the ones the cache lacks for n_bases."""
    path = csv_path(results)
    cached = {k: r for k, r in shared.read_csv_rows(path, _key).items() if str(r.get("n_bases")) == str(n_bases)}
    out, new = [], False
    for d in dens_rows:
        d = dict(d)
        if d.get("scheme") and d.get("toggle_hex") and d.get("L"):
            key = _key(d)
            if key not in cached:
                cached[key] = {**{f: d.get(f) for f in FIELDS}, "n_bases": n_bases,
                               "density_factor": density_factor(d, n_bases)}
                new = True
            fac = float(cached[key]["density_factor"])
            w = int(d["L"]) - int(d["k"]) + 1
            d["fac_expected"] = fac
            d["dens_expected"] = round(fac / (w + 1), 8)
        out.append(d)
    if new:
        shared.write_csv_rows(path, FIELDS, sorted(cached.values(), key=_key))
    return out


def main(argv: list[str] | None = None) -> int:
    import csv

    ap = shared.base_parser(None, "Expected density_factor of the density rows of a kraken2 run "
                                  "(results/runs/<library>/minimizer_density.csv), by Kraken2's own scanner.",
                            default_results=os.environ.get("KRAKEN2_RESULTS_DIR"),
                            default_n_bases=shared.DEFAULT_N_BASES)
    ap.add_argument("density_csv")
    args = ap.parse_args(argv)
    with open(args.density_csv) as f:
        rows = list(csv.DictReader(f))
    if args.force and csv_path(Path(args.results_dir)).is_file():
        csv_path(Path(args.results_dir)).unlink()
    for r in with_expected(rows, Path(args.results_dir), args.n_bases):
        print(r["scheme"], r.get("projection") or "", r.get("toggle_hex"), r["fac_expected"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
