#!/usr/bin/env python3
"""Greedymini runs with the flipped projection (kmc -oat1) are their own variant.

lib/run.sh writes them to a separate ``greedymini/k9_at1-cg0/lut_binary`` run dir
and meta.json records ``scheme=greedymini, projection=at1-cg0``; export must keep
them out of the at0 row (summary aggregation keys on the full scheme identity --
scheme, gm_k, projection, gm_table, gm_lut -- so a shared key would pool both
projections' seeds into one mean/std).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.export import (
    iter_run_dirs,
    load_run_row,
    sort_export_rows,
    summary_xlsx_rows,
)


def _write_run(runs_root: Path, scheme_dir: str, seed_idx: int, *, at1: bool, n_super: int) -> None:
    rid = f"f_vesca/k28_p9/gzip_n512_t12_m12/{scheme_dir}/{seed_idx}"
    run_dir = runs_root / rid
    run_dir.mkdir(parents=True)
    seed = "0xdeadbeef"
    label = f"gm_w20_k9_canonical_s{seed}"
    command = ["/usr/bin/time", "-v", "kmc", "-k28", f"-ol{label}", f"-os{int(seed, 16)}"]
    if at1:
        command.append("-oat1")
    meta = {
        "run_id": rid,
        "dataset": "f_vesca",
        "label": "F. vesca",
        "k": 28,
        "input_format": "gzip",
        "threads": 12,
        "mem_limit_gb": 12,
        "n_bins": 512,
        "strict_memory": False,
        "config_tag": "gzip_n512_t12_m12",
        "signature_len": 9,
        "scheme": "greedymini",
        "gm_k": 9,
        "projection": "at1-cg0" if at1 else "at0-cg1",
        "gm_table": "w20_k9_canonical",
        "gm_lut": "binary",
        "signature_order": "custom",
        "signature_order_mode": "greedymini",
        "signature_order_label": label,
        "signature_order_seed": int(seed, 16),
        "seed_idx": seed_idx,
        "exit_code": 0,
        "python_wall_s": 1.0,
        "command": command,
    }
    (run_dir / "meta.json").write_text(json.dumps(meta))
    kmc = {
        "Total": "10.0s",
        "Stats": {
            "#Unique_k-mers": 500,
            "#Unique_counted_k-mers": 400,
            "#Total no. of k-mers": 1000,
            "#Total_reads": 10,
            "#Total_super-k-mers": n_super,
        },
    }
    (run_dir / "kmc.json").write_text(json.dumps(kmc))


def _rows(tmp_path: Path) -> list[dict]:
    runs_root = tmp_path / "runs"
    for i in range(2):
        _write_run(runs_root, "greedymini/k9_at0-cg1/lut_binary", i, at1=False, n_super=100)
        _write_run(runs_root, "greedymini/k9_at1-cg0/lut_binary", i, at1=True, n_super=110)
    rows = [load_run_row(d, runs_root=runs_root) for d in iter_run_dirs(runs_root)]
    return sort_export_rows([r for r in rows if r])


def test_at1_runs_get_their_own_projection(tmp_path):
    rows = _rows(tmp_path)
    assert {r["scheme"] for r in rows} == {"greedymini"}
    proj = {r["run_id"].split("/")[-3]: r["projection"] for r in rows}
    assert proj == {"k9_at0-cg1": "at0-cg1", "k9_at1-cg0": "at1-cg0"}
    assert {r["gm_k"] for r in rows} == {9}
    assert {r["gm_table"] for r in rows} == {"w20_k9_canonical"}
    assert {r["gm_lut"] for r in rows} == {"binary"}


def test_at1_not_pooled_into_at0_summary(tmp_path):
    summary = {r["projection"]: r for r in summary_xlsx_rows(_rows(tmp_path))}
    assert set(summary) == {"at0-cg1", "at1-cg0"}
    assert summary["at0-cg1"]["n_repeats"] == 2
    assert summary["at1-cg0"]["n_repeats"] == 2
    assert summary["at0-cg1"]["n_super_kmers_mean"] == 100
    assert summary["at1-cg0"]["n_super_kmers_mean"] == 110
    assert summary["at0-cg1"]["n_super_kmers_std"] == 0


def test_export_layout(tmp_path):
    """results/runs -> results/csvs/*.csv + results/*.xlsx."""
    from scripts import export

    results = tmp_path / "results"
    runs_root = results / "runs"
    for i in range(2):
        _write_run(runs_root, "greedymini/k9_at0-cg1/lut_binary", i, at1=False, n_super=100)
    csvs = results / "csvs"
    csvs.mkdir()
    (csvs / "kmc3_density_expected_schemes.csv").write_text(
        "scheme,projection,k,w,L,n_seeds,density_factor_mean,density_factor_std\n"
        "greedymini,at0-cg1,9,20,28,5,1.7,0.01\n"
    )
    export.main(["--results-dir", str(results), "--paper-csv", str(tmp_path / "none.csv")])
    for name in ("kmc3_runs.csv", "kmc3_summary.csv", "kmc3_schemes.csv"):
        assert (csvs / name).is_file()
    assert (results / "kmc3_results.xlsx").is_file()
    assert not list(results.glob("*.csv"))


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        test_at1_runs_get_their_own_projection(Path(td) / "a")
        test_at1_not_pooled_into_at0_summary(Path(td) / "b")
    print("ok")
