#!/usr/bin/env python3
"""Offline tests for scripts/export.py's aggregation/validation/CSV/XLSX
logic -- no sshash binary, no network, no run data required (synthetic rows)."""

from __future__ import annotations

import sys
import tempfile
import zipfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))

from scripts.export import (
    attach_native_ratios,
    summary_field_order,
    summary_rows,
    validate_orders_equivalence,
    vs_reference_rows,
    write_csv,
    write_xlsx,
    _dataset_last_row_indices,
    _sheet_xml,
    _vertical_merge_plan,
)


def _row(dataset, order, seed_idx, num_kmers=1000, check_ok=True, pos_match=True, neg_match=True, **extra):
    row = {
        "dataset": dataset,
        "k": 31,
        "m": 11,
        "order": order,
        "order_seed": seed_idx,
        "seed_idx": seed_idx,
        "gm_file": "",
        "threads": 8,
        "num_kmers": num_kmers,
        "index_bytes_on_disk": 500,
        "bits_per_kmer": 4.0,
        "build_wall_s": 10.0 + seed_idx,
        "build_max_rss_kb": 1000.0,
        "bench_positive_ns": 300.0,
        "bench_negative_ns": 250.0,
        "bench_access_ns": 150.0,
        "query_wall_s": None,
        "query_elapsed_ms": None,
        "check_ok": check_ok,
        "query_pos_match": pos_match,
        "query_neg_match": neg_match,
        "run_id": f"{dataset}/{order}/{seed_idx}",
    }
    row.update(extra)
    return row


def test_validate_orders_equivalence_ok() -> None:
    rows = [
        _row("cod", "native", 0),
        _row("cod", "native", 1),
        _row("cod", "xor-hash", 0),
        _row("cod", "greedymini", 0),
    ]
    n_cmp, errors = validate_orders_equivalence(rows)
    assert n_cmp == 4
    assert errors == []
    print("ok  validate_orders_equivalence: no errors on matching runs")


def test_validate_orders_equivalence_catches_kmer_mismatch() -> None:
    rows = [
        _row("cod", "native", 0, num_kmers=1000),
        _row("cod", "xor-hash", 0, num_kmers=999),
    ]
    _, errors = validate_orders_equivalence(rows)
    assert any("num_kmers differs" in e for e in errors), errors
    print("ok  validate_orders_equivalence: catches num_kmers mismatch")


def test_validate_orders_equivalence_catches_query_mismatch() -> None:
    rows = [
        _row("cod", "native", 0),
        _row("cod", "xor-hash", 0, pos_match=False),
        _row("cod", "greedymini", 0, neg_match=False),
        _row("cod", "xor-pcg-hash", 0, check_ok=False),
    ]
    _, errors = validate_orders_equivalence(rows)
    assert any("query_pos_match is false" in e for e in errors), errors
    assert any("query_neg_match is false" in e for e in errors), errors
    assert any("EVERYTHING OK" in e for e in errors), errors
    print("ok  validate_orders_equivalence: catches query_pos_match/query_neg_match/check_ok failures")


def test_summary_rows_and_ratios() -> None:
    rows = [
        _row("cod", "native", 0, build_wall_s=10.0, bits_per_kmer=5.0),
        _row("cod", "native", 1, build_wall_s=12.0, bits_per_kmer=5.0),
        _row("cod", "greedymini", 0, build_wall_s=8.0, bits_per_kmer=4.0),
        _row("cod", "greedymini", 1, build_wall_s=9.0, bits_per_kmer=4.0),
    ]
    summary = summary_rows(rows)
    assert len(summary) == 2
    attach_native_ratios(summary)
    gm = next(r for r in summary if r["order"] == "greedymini")
    native = next(r for r in summary if r["order"] == "native")
    assert native["build_wall_s_ratio"] == 1.0
    assert gm["build_wall_s_ratio"] < 1.0  # greedymini built faster in this synthetic example
    assert gm["bits_per_kmer_ratio"] < 1.0
    print("ok  summary_rows + attach_native_ratios: ratios computed vs native baseline")


def test_vs_reference_rows() -> None:
    rows = [_row("cod", "native", 0, num_kmers=502_500_000)]
    out = vs_reference_rows(rows)
    cod = next(r for r in out if r["dataset"] == "cod")
    assert cod["measured_num_kmers"] == 502_500_000
    assert abs(cod["measured_vs_published_ratio"] - 1.0) < 0.01
    print("ok  vs_reference_rows: informational ratio computed")


def test_write_csv_and_xlsx_roundtrip() -> None:
    rows = [_row("cod", "native", 0), _row("cod", "greedymini", 0)]
    summary = summary_rows(rows)
    attach_native_ratios(summary)
    with tempfile.TemporaryDirectory(prefix="sshash_export_test_") as td:
        out_dir = Path(td)
        write_csv(out_dir / "sshash_runs.csv", rows)
        write_csv(out_dir / "sshash_summary.csv", summary, summary_field_order())
        assert (out_dir / "sshash_runs.csv").is_file()
        assert (out_dir / "sshash_summary.csv").is_file()

        xlsx_path = out_dir / "sshash_results.xlsx"
        write_xlsx(
            xlsx_path,
            [
                ("summary", summary_field_order(), summary),
                ("runs", list(rows[0].keys()), rows),
            ],
        )
        assert xlsx_path.is_file()
        with zipfile.ZipFile(xlsx_path) as z:
            names = z.namelist()
            assert "xl/workbook.xml" in names
            assert "xl/worksheets/sheet1.xml" in names
            assert "xl/worksheets/sheet2.xml" in names
    print("ok  write_csv / write_xlsx: produce a valid zip with expected parts")


def test_dataset_last_row_indices() -> None:
    rows = [
        _row("cod", "native", 0), _row("cod", "native", 1),
        _row("kestrel", "native", 0), _row("kestrel", "native", 1), _row("kestrel", "native", 2),
    ]
    assert _dataset_last_row_indices(rows) == {1, 4}
    print("ok  _dataset_last_row_indices: border at each dataset's last row")


def test_vertical_merge_plan_stops_at_dataset_boundary() -> None:
    # k/m/w are constant across the whole benchmark, so without a dataset-
    # boundary check a naive merge would span cod+kestrel as one run.
    rows = [_row("cod", "native", i) for i in range(3)] + [_row("kestrel", "native", i) for i in range(2)]
    fields = ["dataset", "k", "m", "w", "order"]
    refs, skip, anchors = _vertical_merge_plan(rows, fields, frozenset({"dataset", "k", "m", "w"}))
    # Two separate merge runs per merged field (cod: rows 0-2, kestrel: rows 3-4),
    # for each of dataset/k/m/w -> 8 merge refs; "order" is not a merge field.
    assert len(refs) == 8, refs
    assert (0, "dataset") in anchors and (3, "dataset") in anchors
    assert (1, "dataset") in skip and (2, "dataset") in skip and (4, "dataset") in skip
    assert not any(field == "order" for (_, field) in anchors)
    print("ok  _vertical_merge_plan: merges within a dataset, never across one")


def test_sheet_xml_has_category_band_header_and_styles() -> None:
    rows = [_row("cod", "native", 0), _row("cod", "native", 1), _row("cod", "greedymini", 0)]
    fields = ["dataset", "k", "m", "w", "order", "num_kmers", "bits_per_kmer", "check_ok"]
    xml = _sheet_xml("runs", fields, rows)
    # row 1 = category band (merged dataset/k/m/w/order under "id / params",
    # num_kmers/bits_per_kmer under "space", check_ok under "correctness")
    assert '<row r="1">' in xml
    assert "id / params" in xml
    assert "space" in xml
    assert "correctness" in xml
    # row 2 = column headers
    assert '<row r="2">' in xml
    # data starts row 3
    assert '<row r="3">' in xml
    # k/m/w/dataset get vertically merged since all 3 rows share one dataset
    assert "<mergeCells" in xml
    # single dataset here -> exactly one divider border, on the last data row
    assert ' s="1"' in xml
    print("ok  _sheet_xml: category band + column headers + merges + border style present")


def test_write_xlsx_includes_styles_part() -> None:
    rows = [_row("cod", "native", 0)]
    with tempfile.TemporaryDirectory(prefix="sshash_export_test_") as td:
        xlsx_path = Path(td) / "out.xlsx"
        write_xlsx(xlsx_path, [("runs", list(rows[0].keys()), rows)])
        with zipfile.ZipFile(xlsx_path) as z:
            names = z.namelist()
            assert "xl/styles.xml" in names
            styles = z.read("xl/styles.xml").decode()
            assert "cellXfs" in styles
    print("ok  write_xlsx: styles.xml part present (dataset-border/bold styles)")


def main() -> int:
    test_validate_orders_equivalence_ok()
    test_validate_orders_equivalence_catches_kmer_mismatch()
    test_validate_orders_equivalence_catches_query_mismatch()
    test_summary_rows_and_ratios()
    test_vs_reference_rows()
    test_write_csv_and_xlsx_roundtrip()
    test_dataset_last_row_indices()
    test_vertical_merge_plan_stops_at_dataset_boundary()
    test_sheet_xml_has_category_band_header_and_styles()
    test_write_xlsx_includes_styles_part()
    print("all export tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
