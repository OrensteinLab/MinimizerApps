#!/usr/bin/env python3
"""Offline tests for scripts/export.py's aggregation/validation/CSV/XLSX
logic -- no sshash binary, no network, no run data required (synthetic rows)."""

from __future__ import annotations

import sys
import tempfile
import zipfile
from pathlib import Path


def test_log2_bucket_bins_follow_the_skew_partitions():
    # sizes: 1 | 2 | 3,4 | 5..8 | 33..64 (last light bin) | 65..128 (skew partition 0) | 4096
    hist = [(1, 100, 700), (2, 10, 40), (3, 4, 20), (4, 2, 12), (8, 1, 9),
            (33, 1, 33), (64, 2, 130), (65, 1, 65), (128, 1, 128), (4096, 1, 4096)]
    bins = log2_bucket_bins(hist, overflow=(3, 99999), heavy_lower=64)
    by_range = {b["bucket_size_range"]: b for b in bins}
    assert [b["bucket_size_range"] for b in bins] == ["1", "2", "3-4", "5-8", "33-64", "65-128", "2049-4096", ">4096"]
    assert by_range["3-4"]["num_buckets"] == 6 and by_range["3-4"]["num_kmers_in_buckets"] == 32
    assert by_range["33-64"]["num_buckets"] == 3
    assert [b["bucket_type"] for b in bins] == ["singleton", "light", "light", "light", "light", "heavy", "heavy", "heavy"]
    assert by_range["65-128"]["log2_bucket_size"] == 7 and by_range[">4096"]["log2_bucket_size"] == 13
    # nothing lost: totals are preserved
    assert sum(b["num_kmers_in_buckets"] for b in bins) == sum(n for _, _, n in hist) + 99999


def test_correlation_rows_formulas_and_cached_values():
    summary = [
        {"density_factor_ratio": 1.0, "num_buckets_ratio": 1.0, "bits_per_kmer_ratio": 1.0, "const_ratio": 1.0},
        {"density_factor_ratio": 0.9, "num_buckets_ratio": 0.9, "bits_per_kmer_ratio": 0.92, "const_ratio": 1.0},
        {"density_factor_ratio": 0.8, "num_buckets_ratio": 0.8, "bits_per_kmer_ratio": 0.85, "const_ratio": 1.0},
    ]
    fields = ["density_factor_ratio", "num_buckets_ratio", "bits_per_kmer_ratio", "const_ratio"]
    rows = correlation_rows(fields, summary)
    assert [r["label"] for r in rows] == ["correlation with density", "correlation with distinct minimizers"]
    formula, value = rows[0]["cells"][2]
    assert formula == 'IFERROR(CORREL($A$3:$A$5,C$3:C$5),"")'  # 3 data rows: 3..5
    assert value == _pearson([1.0, 0.9, 0.8], [1.0, 0.92, 0.85])
    assert rows[1]["cells"][2][0].startswith("IFERROR(CORREL($B$3:$B$5,")
    assert rows[0]["cells"][3][1] is None  # constant column: no correlation
    xml = _sheet_xml("summary", fields, summary, None, rows)
    assert "<f>IFERROR(CORREL(" in xml
    # no colour formatting is attached to the correlation rows (only the _ratio scales of the table)
    assert xml.count("<conditionalFormatting") == 0 or "IFERROR" not in xml.split("<conditionalFormatting")[1]


def test_space_rows_saving_split():
    def run(order, seed, **c):
        base = {"strings_bits_per_kmer": 2.0, "strings_offsets_bits_per_kmer": 0.1, "mphf_bits_per_kmer": 0.4,
                "control_codewords_bits_per_kmer": 4.0, "mid_load_bits_per_kmer": 0.5,
                "skew_mphf_v_bits_per_kmer": 0.1, "skew_h_bits_per_kmer": 0.05, "begin_buckets_bits_per_kmer": 0.0}
        base.update(c)
        return _row("cod", order, seed, **base)

    rows = [run("native", 0), run("native", 1),
            run("greedymini-k13", 0, control_codewords_bits_per_kmer=3.0, mphf_bits_per_kmer=0.3,
                mid_load_bits_per_kmer=0.4, skew_h_bits_per_kmer=0.07),
            run("greedymini-k13", 1, control_codewords_bits_per_kmer=3.0, mphf_bits_per_kmer=0.3,
                mid_load_bits_per_kmer=0.4, skew_h_bits_per_kmer=0.07)]
    summary = summary_rows(rows)
    attach_native_ratios(summary)
    out = {r["order"]: r for r in space_rows(rows, summary)}
    gm = out["greedymini-k13"]
    assert gm["saved_tags"] == 1.0 and gm["saved_mphf"] == 0.1 and gm["saved_light_L"] == 0.1
    assert gm["saved_skew_H"] == -0.02 and gm["saved_strings"] == 0.0
    assert gm["saved_total"] == 1.18
    assert gm["saved_by_minimizers"] == 1.1 and gm["saved_by_occurrences"] == 0.08
    assert abs(gm["share_by_minimizers_pct"] + gm["share_by_occurrences_pct"] - 100) < 0.05
    assert out["native"]["saved_total"] == 0 and out["native"]["share_tags_pct"] is None
    assert set(space_field_order()) >= set(gm)


if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))

from scripts.export import (
    attach_native_ratios,
    bucket_distribution_metrics,
    histogram_overflow,
    log2_bucket_bins,
    correlation_rows,
    space_rows,
    space_field_order,
    _pearson,
    summary_field_order,
    summary_rows,
    validate_orders_equivalence,
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
        "modality": "canonical",
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


def test_histogram_overflow_makes_the_histogram_exact() -> None:
    # sshash does not tally buckets above 4096 per size: 90 singletons + 5 buckets of size 2
    # list 100 k-mers in 95 buckets, but the run has 110 k-mers in 96 buckets (one big bucket).
    listed = [(1, 90, 90), (2, 5, 10)]
    assert histogram_overflow(listed, 110, 96) == (1, 10)
    assert histogram_overflow(listed, 100, 95) is None  # nothing missing
    assert histogram_overflow(listed, 105, 95) is None  # inconsistent: k-mers but no bucket
    assert histogram_overflow([], 110, 96) is None  # no histogram
    # the shares use the true totals once the overflow is included
    plain = bucket_distribution_metrics(listed)
    exact = bucket_distribution_metrics(listed, (1, 10))
    assert plain["kmers_in_singleton_bucket_pct"] == 90.0
    assert exact["kmers_in_singleton_bucket_pct"] == round(100 * 90 / 110, 4)
    assert exact["expected_bucket_cost"] > plain["expected_bucket_cost"]


def test_validate_orders_equivalence_ok() -> None:
    rows = [
        _row("cod", "native", 0),
        _row("cod", "native", 1),
        _row("cod", "xor-hash", 0),
        _row("cod", "greedymini-k15", 0),
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
        _row("cod", "greedymini-k15", 0, neg_match=False),
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
        _row("cod", "greedymini-k15", 0, build_wall_s=8.0, bits_per_kmer=4.0),
        _row("cod", "greedymini-k15", 1, build_wall_s=9.0, bits_per_kmer=4.0),
    ]
    summary = summary_rows(rows)
    assert len(summary) == 2
    attach_native_ratios(summary)
    gm = next(r for r in summary if r["order"] == "greedymini-k15")
    native = next(r for r in summary if r["order"] == "native")
    assert native["build_wall_s_ratio"] == 1.0
    assert gm["build_wall_s_ratio"] < 1.0  # greedymini built faster in this synthetic example
    assert gm["bits_per_kmer_ratio"] < 1.0
    print("ok  summary_rows + attach_native_ratios: ratios computed vs native baseline")


def test_datasets_csv_rows_are_complete() -> None:
    """Every real dataset row (not the in-repo smoke sample) has what `download` needs: an eulertigs
    file with size, md5 and URL, and a read file with a URL or an ENA accession."""
    import csv

    path = Path(__file__).resolve().parents[1] / "data" / "datasets.csv"
    rows = list(csv.DictReader(path.open()))
    real = [r for r in rows if not r["local_path"]]
    assert {"cod", "kestrel", "human", "ncbi-virus", "se", "hprc"} <= {r["dataset"] for r in real}
    for r in real:
        assert r["eulertigs_filename"] and r["eulertigs_url"].startswith("https://"), r["dataset"]
        assert int(r["eulertigs_size"]) > 0 and len(r["eulertigs_md5"]) == 32, r["dataset"]
        assert r["query_filename"] and len(r["query_md5"]) == 32 and int(r["query_size"]) > 0, r["dataset"]
    print("ok  datasets.csv rows are complete")


def test_write_csv_and_xlsx_roundtrip() -> None:
    rows = [_row("cod", "native", 0), _row("cod", "greedymini-k15", 0)]
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
    rows = [_row("cod", "native", 0), _row("cod", "native", 1), _row("cod", "greedymini-k15", 0)]
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
    test_datasets_csv_rows_are_complete()
    test_write_csv_and_xlsx_roundtrip()
    test_dataset_last_row_indices()
    test_vertical_merge_plan_stops_at_dataset_boundary()
    test_sheet_xml_has_category_band_header_and_styles()
    test_write_xlsx_includes_styles_part()
    test_log2_bucket_bins_follow_the_skew_partitions()
    test_correlation_rows_formulas_and_cached_values()
    test_space_rows_saving_split()
    print("all export tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
