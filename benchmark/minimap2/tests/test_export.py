"""Offline tests for lib/export.py (synthetic runs, no minimap2 needed) and
lib/orders.sh (canonical order names; greedymini projection as a parameter)."""
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
from order_names import ID_FIELDS  # noqa: E402
from xlsx_check import check_xlsx  # noqa: E402
from xlsx_util import load_xlsx  # noqa: E402

VARIANTS = [
    ("native", {}),
    ("xor", {}),
    ("greedymini", {"gm_k": 15, "projection": "at0-cg1", "gm_table": "w10_k15_canonical"}),
    ("greedymini", {"gm_k": 15, "projection": "at1-cg0", "gm_table": "w10_k15_canonical"}),
    ("mod-minimizer-xor", {"mod_t": 9}),
]


def _write_run(root: Path, order: str, params: dict, seed_idx: int, dens: float) -> None:
    from order_names import order_dir_name
    odir = order_dir_name(order, params.get("gm_k"), params.get("projection"), params.get("mod_t"))
    d = root / "runs" / "dsA" / "k19_w10" / odir / str(seed_idx)
    d.mkdir(parents=True)
    meta = {"dataset": "dsA", "order": order, "seed": hex(seed_idx), "seed_idx": seed_idx, "k": 19, "w": 10,
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
    for i, (order, params) in enumerate(VARIANTS):
        for s in range(2):
            _write_run(tmp_path, order, params, s, 1.5 + 0.1 * i)
    fake = {(o, p or ""): {"density_factor_mean": "1.9", "density_factor_std": "0.01"}
            for o, p in expected_density.VARIANTS}
    monkeypatch.setattr(expected_density, "ensure_expected_orders", lambda *a, **k: fake)
    (tmp_path / "data" / "dsA").mkdir(parents=True)
    assert export.main(["--results-dir", str(tmp_path), "--data-dir", str(tmp_path / "data"),
                        "--samtools", ""]) == 0
    xlsx = tmp_path / "minimap2_results.xlsx"
    assert check_xlsx(str(xlsx)) == []
    book = load_xlsx(str(xlsx))
    assert list(book) == ["overview", "orders", "datasets", "summary", "runs"]
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
    greedy = [r for r in runs if r["order"] == "greedymini"]
    assert {r["projection"] for r in greedy} == {"at0-cg1", "at1-cg0"}
    assert {r["gm_k"] for r in greedy} == {"15"}
    assert all(r["mod_t"] == "" for r in greedy)
    summary = list(csv.DictReader((tmp_path / "csvs" / "minimap2_summary.csv").open()))
    assert len(summary) == 5 and summary[0]["order"] == "native"
    assert float(summary[0]["density_factor_ratio"]) == 1.0
    assert "n_repeats" in summary[0] and "n_seeds" not in summary[0] and "map_run_s_mean" in summary[0]
    assert "build_pinned_cpus" not in summary[0] and "index_size_gb_mean" in summary[0]
    orders = list(csv.DictReader((tmp_path / "csvs" / "minimap2_orders.csv").open()))
    assert [(o["order"], o["projection"]) for o in orders] == [
        ("native", ""), ("xor", ""), ("greedymini", "at0-cg1"), ("greedymini", "at1-cg0"),
        ("mod-minimizer-xor", "")]
    assert "--gm-projection at1-cg0" in orders[3]["cli"]
    # unified layout: one row per (order, k, w, L), expected for the row, particular for the dataset
    assert list(orders[0])[:5] == ["order", "k", "w", "L", "gm_k"]
    assert list(orders[0])[-2:] == ["dsA_particular_density_factor_mean", "dsA_particular_density_factor_std"]
    assert [(o["k"], o["w"], o["L"], o["mod_t"]) for o in orders] == [("19", "10", "28", "")] * 4 + [("19", "10", "28", "9")]
    assert all(o["expected_density_factor_mean"] == "1.9" and o["dsA_particular_density_factor_std"] != "" for o in orders)
    assert orders[0]["dsA_particular_density_factor_std"] in ("0", "0.0")  # native: std 0, not blank


def test_orders_one_row_per_klw_and_loud_on_holes():
    datasets = [{"dataset": "A", "k": 19, "w": 10, "mod_t": 9, "preset": "map-pb"},
                {"dataset": "B", "k": 15, "w": 10, "mod_t": None, "preset": "map-ont"}]
    def srow(ds, order, k, t=None):
        return {"dataset": ds, "order": order, "mid_occ_cap": export.DEFAULT_MID_OCC_CAP, "k": k, "w": 10, "mod_t": t,
                "density_factor_mean": 1.5, "density_factor_std": None}
    summary = [srow("A", "native", 19), srow("A", "mod-minimizer-xor", 19, 9),
               srow("B", "native", 15), srow("B", "mod-minimizer-xor", 15, 5)]
    exp = {kwt: {(o, ""): {"density_factor_mean": "1.6", "density_factor_std": "0"} for o in ("native", "mod-minimizer-xor")}
           for kwt in ((19, 10, 9), (15, 10, 5))}
    rows = export.orders_rows(summary, datasets, exp)
    assert [(r["order"], r["k"], r["w"], r["L"], r["mod_t"]) for r in rows] == [
        ("native", 19, 10, 28, None), ("native", 15, 10, 24, None),
        ("mod-minimizer-xor", 19, 10, 28, 9), ("mod-minimizer-xor", 15, 10, 24, 5)]
    assert rows[0]["A_particular_density_factor_mean"] == 1.5 and rows[0]["A_particular_density_factor_std"] == 0.0
    assert rows[0]["B_particular_density_factor_mean"] is None   # other (k, w, L): accepted empty
    export.check_orders_complete(rows, datasets)
    rows[3]["B_particular_density_factor_mean"] = None
    with pytest.raises(SystemExit):
        export.check_orders_complete(rows, datasets)
    rows[3]["B_particular_density_factor_mean"] = 1.0
    rows[1]["expected_density_factor_mean"] = None
    with pytest.raises(SystemExit):
        export.check_orders_complete(rows, datasets)


def test_legacy_order_is_rejected(tmp_path):
    _write_run(tmp_path, "native", {}, 0, 1.0)
    meta = next((tmp_path / "runs").glob("**/meta.json"))
    m = json.loads(meta.read_text())
    m["order"] = "greedymini-at1"
    meta.write_text(json.dumps(m))
    with pytest.raises(SystemExit):
        export.load_runs(tmp_path)


def _sh(snippet: str) -> str:
    script = f'set -u; HERE_DIR={HERE}; source {HERE}/lib/constants.sh; source {HERE}/lib/orders.sh; {snippet}'
    return subprocess.run(["bash", "-c", script], check=True, capture_output=True, text=True).stdout.strip()


def test_orders_sh_dirs_and_flags():
    assert _sh("order_dir_name xor") == "xor"
    assert _sh("order_dir_name greedymini") == "greedymini/k15_at0-cg1"
    assert _sh("PROJECTION=at1-cg0; order_dir_name greedymini") == "greedymini/k15_at1-cg0"
    assert _sh("order_dir_name mod-minimizer-xor 9") == "mod-minimizer-xor/t9"
    flags = _sh("PROJECTION=at1-cg0; build_order_flags greedymini 7 0; echo \"${ORDER_FLAGS[*]}\"")
    assert "--minimizer-order greedymini" in flags and "--gm-projection at1-cg0" in flags
    assert _sh("order_is_known greedymini-at1 && echo yes || echo no") == "no"
