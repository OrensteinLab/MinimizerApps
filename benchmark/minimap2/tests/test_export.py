"""Offline tests for lib/export.py (synthetic runs, no minimap2 needed) and
lib/schemes.sh (canonical scheme names; greedymini projection as a parameter)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
REPO_ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE / "lib"))
sys.path.insert(0, str(REPO_ROOT / "benchmark" / "common"))

import export  # noqa: E402
import expected_density  # noqa: E402
import step_record  # noqa: E402
from scheme_names import ID_FIELDS  # noqa: E402
from xlsx_check import check_xlsx  # noqa: E402
from xlsx_util import load_xlsx  # noqa: E402

VARIANTS = [
    ("native", {}),
    ("xor", {}),
    ("greedymini", {"gm_k": 15, "projection": "at0-cg1", "gm_table": "w10_k15_canonical"}),
    ("greedymini", {"gm_k": 15, "projection": "at1-cg0", "gm_table": "w10_k15_canonical"}),
    ("mod-minimizer-xor", {"mod_t": 9}),
]


def _write_run(root: Path, scheme: str, params: dict, seed_idx: int, dens: float) -> None:
    from scheme_names import scheme_dir_name
    sdir = scheme_dir_name(scheme, params.get("gm_k"), params.get("projection"), params.get("mod_t"))
    d = root / "runs" / "dsA" / "k19_w10" / sdir / str(seed_idx)
    d.mkdir(parents=True)
    meta = {"dataset": "dsA", "scheme": scheme, "seed": hex(seed_idx), "seed_idx": seed_idx, "k": 19, "w": 10,
            "gm_k": params.get("gm_k"), "projection": params.get("projection"),
            "gm_table": params.get("gm_table"), "mod_t": params.get("mod_t"),
            "preset": "map-pb", "threads": 4, "distinct_minimizers": 1000, "avg_occurrences": 1.2,
            "density_factor": dens,
            "build": {"wall_time_sec": 2.0 + dens, "peak_rss_kb": 2048},
            "map": {"wall_time_sec": 1.0 + dens, "peak_rss_kb": 4096, "index_load_s": 0.5, "pinned_cpus": "2,6,10"},
            "accuracy": {"n_reads": 10, "frac_mapped": 0.9, "frac_correct": 0.8 - dens / 100,
                         "frac_incorrectly_mapped": 0.1, "mean_identity": 0.88}}
    (d / "meta.json").write_text(json.dumps(meta))


def test_export_roundtrip(tmp_path, monkeypatch):
    for i, (scheme, params) in enumerate(VARIANTS):
        for s in range(2):
            _write_run(tmp_path, scheme, params, s, 1.5 + 0.1 * i)
    fake = {(o, p or ""): {"density_factor_mean": "1.9", "density_factor_std": "0.01"}
            for o, p in expected_density.VARIANTS}
    monkeypatch.setattr(expected_density, "ensure_expected_schemes", lambda *a, **k: fake)
    (tmp_path / "data" / "dsA").mkdir(parents=True)
    assert export.main(["--results-dir", str(tmp_path), "--data-dir", str(tmp_path / "data"),
                        "--samtools", ""]) == 0
    xlsx = tmp_path / "minimap2_results.xlsx"
    assert check_xlsx(str(xlsx)) == []
    book = load_xlsx(str(xlsx))
    # no fixed-threshold block (-f), so no run is paired with native's and there are no density_effect sheets
    assert list(book) == ["overview", "schemes", "datasets", "summary", "runs"]
    import csv
    runs = list(csv.DictReader((tmp_path / "csvs" / "minimap2_runs.csv").open()))
    assert len(runs) == 10
    hdr = list(runs[0])
    assert hdr[:3] == ["dataset", "n_reads", "n_reads_in_repeat"] and hdr[3:6] == ["k", "w", "L"]
    assert hdr[6:11] == list(ID_FIELDS)
    assert hdr[11:15] == ["mid_occ_threshold", "threads", "seed", "seed_idx"]
    assert hdr.index("avg_occurrences") + 1 == hdr.index("index_size_gb") and hdr.index("index_size_gb") + 1 == hdr.index("total_wall_s")
    n_rec = 2 * len(step_record.FIELDS)   # step-record columns of build and map at the end
    assert hdr[-2 - n_rec:-n_rec] == ["build_pinned_cpus", "map_pinned_cpus"]
    assert hdr[-n_rec:] == step_record.names("build") + step_record.names("map")
    assert {"map_load_s", "map_run_s"} <= set(hdr) and not {"index_load_s", "map_only_s", "n_seeds"} & set(hdr)
    assert runs[0]["map_pinned_cpus"] == "2,6,10" and runs[0]["map_load_s"] == "0.5" and runs[0]["map_run_s"] != ""
    greedy = [r for r in runs if r["scheme"] == "greedymini"]
    assert {r["projection"] for r in greedy} == {"at0-cg1", "at1-cg0"}
    assert {r["gm_k"] for r in greedy} == {"15"}
    assert all(r["mod_t"] == "" for r in greedy)
    summary = list(csv.DictReader((tmp_path / "csvs" / "minimap2_summary.csv").open()))
    assert len(summary) == 5 and summary[0]["scheme"] == "native"
    assert float(summary[0]["density_factor_ratio"]) == 1.0
    assert "n_repeats" in summary[0] and "n_seeds" not in summary[0] and "map_run_s_mean" in summary[0]
    assert "build_pinned_cpus" not in summary[0] and "index_size_gb_mean" in summary[0]
    schemes = list(csv.DictReader((tmp_path / "csvs" / "minimap2_schemes.csv").open()))
    assert [(o["scheme"], o["projection"]) for o in schemes] == [
        ("native", ""), ("xor", ""), ("greedymini", "at0-cg1"), ("greedymini", "at1-cg0"),
        ("mod-minimizer-xor", "")]
    assert "--gm-projection at1-cg0" in schemes[3]["cli"]
    # unified layout: one row per (scheme, k, w, L), expected for the row, particular for the dataset
    assert list(schemes[0])[:5] == ["scheme", "k", "w", "L", "gm_k"]
    assert list(schemes[0])[-2:] == ["dsA_particular_density_factor_mean", "dsA_particular_density_factor_std"]
    assert [(o["k"], o["w"], o["L"], o["mod_t"]) for o in schemes] == [("19", "10", "28", "")] * 4 + [("19", "10", "28", "9")]
    assert all(o["expected_density_factor_mean"] == "1.9" and o["dsA_particular_density_factor_std"] != "" for o in schemes)
    assert schemes[0]["dsA_particular_density_factor_std"] in ("0", "0.0")  # native: std 0, not blank


def test_schemes_one_row_per_klw_and_loud_on_holes():
    datasets = [{"dataset": "A", "k": 19, "w": 10, "mod_t": 9, "preset": "map-pb"},
                {"dataset": "B", "k": 15, "w": 10, "mod_t": None, "preset": "map-ont"}]
    def srow(ds, scheme, k, t=None):
        return {"dataset": ds, "scheme": scheme, "mid_occ_cap": export.DEFAULT_MID_OCC_CAP, "k": k, "w": 10, "mod_t": t,
                "density_factor_mean": 1.5, "density_factor_std": None}
    summary = [srow("A", "native", 19), srow("A", "mod-minimizer-xor", 19, 9),
               srow("B", "native", 15), srow("B", "mod-minimizer-xor", 15, 5)]
    exp = {kwt: {(o, ""): {"density_factor_mean": "1.6", "density_factor_std": "0"} for o in ("native", "mod-minimizer-xor")}
           for kwt in ((19, 10, 9), (15, 10, 5))}
    rows = export.schemes_rows(summary, datasets, exp)
    assert [(r["scheme"], r["k"], r["w"], r["L"], r["mod_t"]) for r in rows] == [
        ("native", 19, 10, 28, None), ("native", 15, 10, 24, None),
        ("mod-minimizer-xor", 19, 10, 28, 9), ("mod-minimizer-xor", 15, 10, 24, 5)]
    assert rows[0]["A_particular_density_factor_mean"] == 1.5 and rows[0]["A_particular_density_factor_std"] == 0.0
    assert rows[0]["B_particular_density_factor_mean"] is None   # other (k, w, L): accepted empty
    export.check_schemes_complete(rows, datasets)
    rows[3]["B_particular_density_factor_mean"] = None
    with pytest.raises(SystemExit):
        export.check_schemes_complete(rows, datasets)
    rows[3]["B_particular_density_factor_mean"] = 1.0
    rows[1]["expected_density_factor_mean"] = None
    with pytest.raises(SystemExit):
        export.check_schemes_complete(rows, datasets)


def test_legacy_scheme_is_rejected(tmp_path):
    _write_run(tmp_path, "native", {}, 0, 1.0)
    meta = next((tmp_path / "runs").glob("**/meta.json"))
    m = json.loads(meta.read_text())
    m["scheme"] = "greedymini-at1"
    meta.write_text(json.dumps(m))
    with pytest.raises(SystemExit):
        export.load_runs(tmp_path)


def _sh(snippet: str) -> str:
    script = f'set -u; HERE_DIR={HERE}; source {HERE}/lib/constants.sh; source {HERE}/lib/schemes.sh; {snippet}'
    return subprocess.run(["bash", "-c", script], check=True, capture_output=True, text=True).stdout.strip()


def test_schemes_sh_dirs_and_flags():
    assert _sh("scheme_dir_name xor") == "xor"
    assert _sh("scheme_dir_name greedymini") == "greedymini/k15_at0-cg1"
    assert _sh("PROJECTION=at1-cg0; scheme_dir_name greedymini") == "greedymini/k15_at1-cg0"
    assert _sh("scheme_dir_name mod-minimizer-xor 9") == "mod-minimizer-xor/t9"
    flags = _sh("PROJECTION=at1-cg0; build_scheme_flags greedymini 7 0; echo \"${SCHEME_FLAGS[*]}\"")
    assert "--scheme greedymini" in flags and "--gm-projection at1-cg0" in flags
    assert _sh("scheme_is_known greedymini-at1 && echo yes || echo no") == "no"
    assert _sh("PROJECTION=at1-cg0; scheme_dir_name mod-minimizer-greedymini 9") == "mod-minimizer-greedymini/k9_at1-cg0_t9"
    flags = _sh("PROJECTION=at1-cg0; build_scheme_flags mod-minimizer-greedymini 7 9; echo \"${SCHEME_FLAGS[*]}\"").split()
    assert flags[:4] == ["--scheme", "mod-minimizer-greedymini", "--order-seed", "7"]
    assert flags[5].endswith("/w20_k9_canonical.gm") and flags[4] == "--greedymini-file"
    assert flags[6:] == ["--gm-projection", "at1-cg0", "--mod-t", "9"]


def test_schemes_mod_greedymini_only_where_t_is_gm_k():
    """mod-minimizer-greedymini (a table of 9-mers) gets a schemes row at t = 9 (k=19), none at t = 5 (k=15)."""
    datasets = [{"dataset": "A", "k": 19, "w": 10, "mod_t": 9, "preset": "map-pb"},
                {"dataset": "B", "k": 15, "w": 10, "mod_t": None, "preset": "map-ont"}]
    ident = {"gm_k": 9, "projection": "at0-cg1", "gm_table": "w20_k9_canonical"}
    summary = [{"dataset": "A", "scheme": "mod-minimizer-greedymini", **ident, "mid_occ_cap": export.DEFAULT_MID_OCC_CAP,
                "k": 19, "w": 10, "mod_t": 9, "density_factor_mean": 1.4, "density_factor_std": 0.01}]
    exp = {(19, 10, 9): {("mod-minimizer-greedymini", "at0-cg1"): {"density_factor_mean": "1.45", "density_factor_std": "0.01"}}}
    rows = export.schemes_rows(summary, datasets, exp)
    assert [(r["scheme"], r["k"], r["mod_t"]) for r in rows] == [("mod-minimizer-greedymini", 19, 9)]
    assert "w20_k9_canonical.gm" in rows[0]["cli"] and "--mod-t 9" in rows[0]["cli"]
    export.check_schemes_complete(rows, datasets)
