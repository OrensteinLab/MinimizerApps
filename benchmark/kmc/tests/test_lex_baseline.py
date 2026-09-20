#!/usr/bin/env python3
"""Tests for the repeated/interleaved lexicographic baseline (see README:
'Run order and the lex baseline').

Once lib/run.sh repeats `default` (lex) the same way it repeats custom
schemes, these generic exporters must treat lex like any other scheme:
aggregate its repeats into a mean/std (not a single sample), and expose a
Welch's-t `{metric}_z` alongside `{metric}_ratio` so a ratio far from 1 can
be told apart from ordinary seed noise.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.export import (
    XLSX_DATA_START_ROW,
    XLSX_MERGE_FIELDS,
    Z_SCORE_METRICS,
    _vertical_merge_plan,
    _xlsx_col,
    prepare_agg_xlsx_sheet,
    summary_xlsx_field_order,
    summary_xlsx_rows,
    vs_paper_rows,
)


def _lex_row(wall: float, seed_idx: int = 0, *, n_bins: int = 512, threads: int = 12) -> dict:
    return {
        "dataset": "f_vesca",
        "k": 28,
        "input_format": "gzip",
        "config_tag": f"gzip_n{n_bins}_t{threads}_m12",
        "order_label": "default",
        "threads": threads,
        "n_bins": n_bins,
        "mem_limit": 12,
        "seed_idx": seed_idx,
        "wall_time": wall,
        "stage1": wall * 0.6,
        "stage2": wall * 0.4,
        "stage3": 0.0,
        "disk": 4.0,
        "n_super_kmers": 120,
        "avg_super_kmer_len": 11.0,
        "density_factor": 2.0,
        "n_unique_counted": 200,
        "n_unique": 200,
        "n_total_kmers": 1000,
        "n_reads": 100,
        "p": 9,
        "w": 20,
        "ram": 11.2,
        "strict_memory": False,
        "signature_order": "default",
    }


def _gm_row(wall: float, seed_idx: int, **kw) -> dict:
    return {
        **_lex_row(wall, seed_idx, **kw),
        "order_label": "greedymini_dna",
        "n_super_kmers": 100,
        "density_factor": 1.6,
        "signature_order": "custom",
    }


def _hash_row(wall: float, seed_idx: int, **kw) -> dict:
    return {**_lex_row(wall, seed_idx, **kw), "order_label": "hash", "signature_order": "custom"}


def _pcg_row(wall: float, seed_idx: int, **kw) -> dict:
    return {**_lex_row(wall, seed_idx, **kw), "order_label": "pcg", "signature_order": "custom"}


def test_repeated_lex_gets_real_mean_and_std():
    lex_reps = [_lex_row(30.0 + d, i) for i, d in enumerate([0.0, 0.5, -0.5, 0.2, -0.2])]
    summary = summary_xlsx_rows(lex_reps)
    assert len(summary) == 1
    row = summary[0]
    assert row["n_seeds"] == 5
    assert row["wall_time_mean"] == 30.0
    assert row["wall_time_std"] is not None and row["wall_time_std"] > 0


def test_z_score_column_present_and_blank_for_lex():
    lex_reps = [_lex_row(30.0 + d, i) for i, d in enumerate([0.0, 0.5, -0.5, 0.2, -0.2])]
    gm_reps = [_gm_row(24.0 + d, i) for i, d in enumerate([0.0, 0.3, -0.3, 0.1, -0.1])]
    summary = summary_xlsx_rows(lex_reps + gm_reps)
    fields = prepare_agg_xlsx_sheet(summary)
    for metric in Z_SCORE_METRICS:
        assert f"{metric}_z" in fields

    lex_row = next(r for r in summary if r["order_label"] == "default")
    gm_row = next(r for r in summary if r["order_label"] == "greedymini_dna")
    assert lex_row["wall_time_z"] == {"f": "=0"}
    assert gm_row["wall_time_z"]["f"].startswith("=IFERROR(")


def test_z_score_blank_when_only_one_side_has_repeats():
    # lex has 5 repeats, gm has only 1 (n=1, std=0) -> IFERROR blanks the divide-by-zero.
    lex_reps = [_lex_row(30.0 + d, i) for i, d in enumerate([0.0, 0.5, -0.5, 0.2, -0.2])]
    gm_reps = [_gm_row(24.0, 0)]
    summary = summary_xlsx_rows(lex_reps + gm_reps)
    fields = prepare_agg_xlsx_sheet(summary)
    gm_row = next(r for r in summary if r["order_label"] == "greedymini_dna")
    assert "wall_time_z" in fields
    assert gm_row["wall_time_z"]["f"].startswith("=IFERROR(")


def test_z_score_resolves_through_vertically_merged_n_seeds():
    """Regression: n_seeds is vertically merged (XLSX_MERGE_FIELDS), so pcg/gm -- sharing
    hash's n_seeds=5 -- only get a *visually* merged cell; hash is the only row with an
    actual <c> for it. A z formula that references n_col{own_row} directly divides by a
    blank cell there and IFERROR silently blanks the result, which is exactly what made
    z come out empty for every non-anchor row in a real export."""
    seeds5 = [0.0, 0.5, -0.5, 0.2, -0.2]
    rows = (
        [_lex_row(30.0 + d, i) for i, d in enumerate(seeds5)]
        + [_hash_row(29.0 + d, i) for i, d in enumerate(seeds5)]
        + [_pcg_row(31.0 + d, i) for i, d in enumerate(seeds5)]
        + [_gm_row(24.0 + d, i) for i, d in enumerate(seeds5)]
    )
    summary = summary_xlsx_rows(rows)
    # sanity: lex/hash/pcg/gm all land on n_seeds=5 (lex is repeated too, post-harness-fix),
    # so the merge groups all four into one run anchored on the first (lex) row -- meaning
    # even hash's own n_seeds cell is merged away, not just pcg/gm's.
    assert [r["n_seeds"] for r in summary] == [5, 5, 5, 5]

    _assert_no_z_formula_hits_merged_away_n_seeds(summary, ("hash", "pcg", "greedymini_dna"))


def test_z_score_resolves_when_lex_has_fewer_repeats_than_custom():
    """Same bug, pre-harness-fix shape: lex n_seeds=1 (its own anchor), hash/pcg/gm share
    n_seeds=5 as one merge run anchored on hash -- this is the exact layout of the real
    export that first exposed the bug (pcg/gm referencing their own blanked n_seeds cell)."""
    seeds5 = [0.0, 0.5, -0.5, 0.2, -0.2]
    rows = (
        [_lex_row(30.0, 0)]
        + [_hash_row(29.0 + d, i) for i, d in enumerate(seeds5)]
        + [_pcg_row(31.0 + d, i) for i, d in enumerate(seeds5)]
        + [_gm_row(24.0 + d, i) for i, d in enumerate(seeds5)]
    )
    summary = summary_xlsx_rows(rows)
    assert [r["n_seeds"] for r in summary] == [1, 5, 5, 5]
    _assert_no_z_formula_hits_merged_away_n_seeds(summary, ("hash", "pcg", "greedymini_dna"))


def _assert_no_z_formula_hits_merged_away_n_seeds(summary: list[dict], orders) -> None:
    fields = prepare_agg_xlsx_sheet(summary)
    merge_fields = XLSX_MERGE_FIELDS["summary"]
    _refs, skip, _anchors = _vertical_merge_plan(summary, fields, merge_fields)
    n_col = _xlsx_col(fields.index("n_seeds"))

    for order in orders:
        row_idx = next(i for i, r in enumerate(summary) if r["order_label"] == order)
        formula = summary[row_idx]["wall_time_z"]["f"]
        assert formula.startswith("=IFERROR("), f"{order}: {formula}"
        # every n_col reference in the formula must be a row that actually has a cell,
        # i.e. must NOT be one of the vertically-merged-away (skipped) rows
        for excel_row in re.findall(rf"(?<![A-Za-z]){n_col}(\d+)", formula):
            ref_idx = int(excel_row) - XLSX_DATA_START_ROW
            assert (ref_idx, "n_seeds") not in skip, (
                f"{order}: {formula} references merged-away n_seeds row {excel_row}"
            )


def test_vs_paper_aggregates_lex_repeats():
    lex_reps = [_lex_row(30.0 + d, i) for i, d in enumerate([0.0, 0.5, -0.5, 0.2, -0.2])]
    out = vs_paper_rows(lex_reps, [], "xeon_hdd")
    assert len(out) == 1, "repeats of the same config must collapse to one vs-paper row"
    row = out[0]
    assert row["n_seeds"] == 5
    assert row["wall_time"] == 30.0
    assert row["wall_time_std"] is not None and row["wall_time_std"] > 0


def test_vs_paper_keeps_configs_separate():
    reps_512 = [_lex_row(30.0, i, n_bins=512) for i in range(3)]
    reps_1024 = [_lex_row(15.0, i, n_bins=1024) for i in range(3)]
    out = vs_paper_rows(reps_512 + reps_1024, [], "xeon_hdd")
    assert len(out) == 2
    bins = sorted(r["n_bins"] for r in out)
    assert bins == [512, 1024]


if __name__ == "__main__":
    test_repeated_lex_gets_real_mean_and_std()
    test_z_score_column_present_and_blank_for_lex()
    test_z_score_blank_when_only_one_side_has_repeats()
    test_z_score_resolves_through_vertically_merged_n_seeds()
    test_z_score_resolves_when_lex_has_fewer_repeats_than_custom()
    test_vs_paper_aggregates_lex_repeats()
    test_vs_paper_keeps_configs_separate()
    print("ok")
