#!/usr/bin/env python3
"""Tests for the repeated/interleaved lexicographic baseline (see README:
'Run order and the lex baseline').

Once lib/run.sh repeats `default` (lex) the same way it repeats custom
schemes, these generic exporters must treat lex like any other scheme:
aggregate its repeats into a mean/std (not a single sample), and expose it
through `{metric}_mean/_std/_ratio`.
"""

from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.export import (
    XLSX_DATA_START_ROW,
    _xlsx_col,
    make_sheet,
    prepare_agg_xlsx_sheet,
    summary_xlsx_rows,
    vs_paper_rows,
    write_book,
)


def _lex_row(wall: float, seed_idx: int = 0, *, n_bins: int = 512, threads: int = 12) -> dict:
    return {
        "dataset": "f_vesca",
        "k": 28,
        "input_format": "gzip",
        "config_tag": f"gzip_n{n_bins}_t{threads}_m12",
        "order": "default",
        "threads": threads,
        "n_bins": n_bins,
        "mem_limit_gb": 12,
        "seed_idx": seed_idx,
        "count_wall_s": wall,
        "count_stage1_s": wall * 0.6,
        "count_stage2_s": wall * 0.4,
        "count_stage3_s": 0.0,
        "tmp_disk_gb": 4.0,
        "n_super_kmers": 120,
        "avg_super_kmer_len": 11.0,
        "density_factor": 2.0,
        "n_unique_counted": 200,
        "n_unique": 200,
        "n_total_kmers": 1000,
        "n_reads": 100,
        "p": 9,
        "w": 20,
        "count_max_rss_gb": 11.2,
        "strict_memory": False,
        "signature_order": "default",
    }


def _gm_row(wall: float, seed_idx: int, **kw) -> dict:
    return {
        **_lex_row(wall, seed_idx, **kw),
        "order": "greedymini",
        "gm_k": 9,
        "projection": "at0-cg1",
        "gm_table": "w20_k9_canonical",
        "gm_lut": "dna",
        "n_super_kmers": 100,
        "density_factor": 1.6,
        "signature_order": "custom",
    }


def _xor_row(wall: float, seed_idx: int, **kw) -> dict:
    return {**_lex_row(wall, seed_idx, **kw), "order": "xor", "signature_order": "custom"}


def _pcg_row(wall: float, seed_idx: int, **kw) -> dict:
    return {**_lex_row(wall, seed_idx, **kw), "order": "xor-pcg-hash", "signature_order": "custom"}


def test_repeated_lex_gets_real_mean_and_std():
    lex_reps = [_lex_row(30.0 + d, i) for i, d in enumerate([0.0, 0.5, -0.5, 0.2, -0.2])]
    summary = summary_xlsx_rows(lex_reps)
    assert len(summary) == 1
    row = summary[0]
    assert row["n_repeats"] == 5
    assert row["count_wall_s_mean"] == 30.0
    assert row["count_wall_s_std"] is not None and row["count_wall_s_std"] > 0


def test_vs_paper_aggregates_lex_repeats():
    lex_reps = [_lex_row(30.0 + d, i) for i, d in enumerate([0.0, 0.5, -0.5, 0.2, -0.2])]
    out = vs_paper_rows(lex_reps, [], "xeon_hdd")
    assert len(out) == 1, "repeats of the same config must collapse to one vs-paper row"
    row = out[0]
    assert row["n_repeats"] == 5
    assert row["count_wall_s"] == 30.0
    assert row["count_wall_s_std"] is not None and row["count_wall_s_std"] > 0


def test_vs_paper_keeps_configs_separate():
    reps_512 = [_lex_row(30.0, i, n_bins=512) for i in range(3)]
    reps_1024 = [_lex_row(15.0, i, n_bins=1024) for i in range(3)]
    out = vs_paper_rows(reps_512 + reps_1024, [], "xeon_hdd")
    assert len(out) == 2
    bins = sorted(r["n_bins"] for r in out)
    assert bins == [512, 1024]


if __name__ == "__main__":
    test_repeated_lex_gets_real_mean_and_std()
    test_vs_paper_aggregates_lex_repeats()
    test_vs_paper_keeps_configs_separate()
    print("ok")
