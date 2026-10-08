#!/usr/bin/env python3
"""GM_PROJECTION=at1-cg0 runs (kmc -oat1) must export as their own scheme.

lib/run.sh writes them to a separate ``greedymini_binary_at1`` scheme dir with
an ``_at1`` order label and ``-oat1`` on the recorded command line; export must
keep them out of the at0 ``greedymini_binary`` row (summary aggregation keys on
order_label, so a shared label would pool both projections' seeds into one
mean/std).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.export import (
    NO_EXPECTED_DENSITY,
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
    label = f"gm_w20_k9_canonical{'_at1' if at1 else ''}_s{seed}"
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
        "signature_order": "custom",
        "signature_order_mode": "table_binary",
        "gm_lut_mode": "binary",
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
        _write_run(runs_root, "greedymini_binary", i, at1=False, n_super=100)
        _write_run(runs_root, "greedymini_binary_at1", i, at1=True, n_super=110)
    rows = [load_run_row(d, runs_root=runs_root) for d in iter_run_dirs(runs_root)]
    return sort_export_rows([r for r in rows if r])


def test_at1_runs_get_their_own_label(tmp_path):
    rows = _rows(tmp_path)
    labels = {r["run_id"].split("/")[-2]: r["order_label"] for r in rows}
    assert labels == {
        "greedymini_binary": "greedymini_binary",
        "greedymini_binary_at1": "greedymini_binary_at1",
    }
    proj = {r["order_label"]: r["gm_projection"] for r in rows}
    assert proj == {"greedymini_binary": "at0-cg1", "greedymini_binary_at1": "at1-cg0"}


def test_at1_not_pooled_into_at0_summary(tmp_path):
    summary = {r["order_label"]: r for r in summary_xlsx_rows(_rows(tmp_path))}
    assert set(summary) == {"greedymini_binary", "greedymini_binary_at1"}
    assert summary["greedymini_binary"]["n_seeds"] == 2
    assert summary["greedymini_binary_at1"]["n_seeds"] == 2
    assert summary["greedymini_binary"]["n_super_kmers_mean"] == 100
    assert summary["greedymini_binary_at1"]["n_super_kmers_mean"] == 110
    assert summary["greedymini_binary"]["n_super_kmers_std"] == 0


def test_at1_expected_density_blank():
    assert "greedymini_binary_at1" in NO_EXPECTED_DENSITY


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        test_at1_runs_get_their_own_label(Path(td) / "a")
        test_at1_not_pooled_into_at0_summary(Path(td) / "b")
    test_at1_expected_density_blank()
    print("ok")
