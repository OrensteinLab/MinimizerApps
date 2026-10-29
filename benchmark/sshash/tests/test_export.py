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
    assert by_range["3-4"]["n_buckets"] == 6 and by_range["3-4"]["n_kmers_in_buckets"] == 32
    assert by_range["33-64"]["n_buckets"] == 3
    assert [b["bucket_type"] for b in bins] == ["singleton", "light", "light", "light", "light", "heavy", "heavy", "heavy"]
    assert by_range["65-128"]["log2_bucket_size"] == 7 and by_range[">4096"]["log2_bucket_size"] == 13
    # nothing lost: totals are preserved
    assert sum(b["n_kmers_in_buckets"] for b in bins) == sum(n for _, _, n in hist) + 99999


def test_correlation_rows_formulas_and_cached_values():
    summary = [
        {"density_factor_ratio": 1.0, "distinct_minimizers_ratio": 1.0, "bits_per_kmer_ratio": 1.0, "const_ratio": 1.0},
        {"density_factor_ratio": 0.9, "distinct_minimizers_ratio": 0.9, "bits_per_kmer_ratio": 0.92, "const_ratio": 1.0},
        {"density_factor_ratio": 0.8, "distinct_minimizers_ratio": 0.8, "bits_per_kmer_ratio": 0.85, "const_ratio": 1.0},
    ]
    fields = ["density_factor_ratio", "distinct_minimizers_ratio", "bits_per_kmer_ratio", "const_ratio"]
    rows = correlation_rows(fields, summary)
    assert [r["label"] for r in rows] == ["correlation with density", "correlation with distinct minimizers"]
    cells = {f: (formula, value) for f, formula, value in rows[0]["cells"]}
    formula, value = cells["bits_per_kmer_ratio"]
    assert formula == 'IFERROR(CORREL($A$3:$A$5,C$3:C$5),"")'  # 3 data rows: 3..5
    assert value == _pearson([1.0, 0.9, 0.8], [1.0, 0.92, 0.85])
    assert rows[1]["cells"][2][1].startswith("IFERROR(CORREL($B$3:$B$5,")
    assert cells["const_ratio"][1] is None  # constant column: no correlation (cached value blank)


def test_space_rows_saving_split():
    def run(scheme, seed, **c):
        base = {"strings_bits_per_kmer": 2.0, "strings_offsets_bits_per_kmer": 0.1, "mphf_bits_per_kmer": 0.4,
                "control_codewords_bits_per_kmer": 4.0, "mid_load_bits_per_kmer": 0.5,
                "skew_mphf_v_bits_per_kmer": 0.1, "skew_h_bits_per_kmer": 0.05, "begin_buckets_bits_per_kmer": 0.0}
        base.update(c)
        return _row("cod", scheme, seed, **(GM13 if scheme == "greedymini" else {}), **base)

    rows = [run("native", 0), run("native", 1),
            run("greedymini", 0, control_codewords_bits_per_kmer=3.0, mphf_bits_per_kmer=0.3,
                mid_load_bits_per_kmer=0.4, skew_h_bits_per_kmer=0.07),
            run("greedymini", 1, control_codewords_bits_per_kmer=3.0, mphf_bits_per_kmer=0.3,
                mid_load_bits_per_kmer=0.4, skew_h_bits_per_kmer=0.07)]
    summary = summary_rows(rows)
    attach_native_ratios(summary)
    out = {r["scheme"]: r for r in space_rows(rows, summary)}
    gm = out["greedymini"]
    assert gm["saved_tags"] == 1.0 and gm["saved_mphf"] == 0.1 and gm["saved_light_L"] == 0.1
    assert gm["saved_skew_H"] == -0.02 and gm["saved_strings"] == 0.0
    assert gm["saved_total"] == 1.18
    assert gm["saved_by_minimizers"] == 1.1 and gm["saved_by_occurrences"] == 0.08
    assert abs(gm["share_by_minimizers_pct"] + gm["share_by_occurrences_pct"] - 100) < 0.05
    assert out["native"]["saved_total"] == 0 and out["native"]["share_tags_pct"] is None
    assert set(space_field_order()) >= set(gm)


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
    validate_schemes_equivalence,
    write_csv,
    HISTOGRAM_ID_FIELDS,
    bucket_rows,
    build_sheets,
    schemes_rows,
    schemes_field_order,
)
from xlsx_book import write_book
import xlsx_check


GM13 = {"gm_k": 13, "projection": "at0-cg1", "gm_table": "w11_k13_canonical"}
GM15 = {"gm_k": 15, "projection": "at0-cg1", "gm_table": "w11_k15_canonical"}


def test_column_order():
    from scripts.export import RUN_FIELDS, summary_field_order, bucket_field_order, space_field_order
    for f in (list(RUN_FIELDS), summary_field_order(), bucket_field_order(), space_field_order()):
        i = f.index("modality")
        assert f[i - 1] == "mod_t" or f[:i][-1] == "mod_t"
        assert f[i + 1] in ("threads", "n_repeats")
    r = list(RUN_FIELDS)
    assert r[r.index("modality"):r.index("modality") + 4] == ["modality", "threads", "seed", "seed_idx"]


def _row(dataset, scheme, seed_idx, n_kmers=1000, check_ok=True, pos_match=True, neg_match=True, **extra):
    row = {
        "dataset": dataset,
        "k": 11,
        "w": 21,
        "L": 31,
        "scheme": scheme,
        "gm_k": None,
        "projection": None,
        "gm_table": None,
        "mod_t": None,
        "modality": "canonical",
        "seed": seed_idx,
        "seed_idx": seed_idx,
        "threads": 8,
        "n_kmers": n_kmers,
        "index_disk_gb": 5e-7,
        "bits_per_kmer": 4.0,
        "build_wall_s": 10.0 + seed_idx,
        "build_max_rss_gb": 0.001,
        "bench_positive_ns": 300.0,
        "bench_negative_ns": 250.0,
        "bench_access_ns": 150.0,
        "query_wall_s": None,
        "query_run_s": None,
        "check_ok": check_ok,
        "query_pos_match": pos_match,
        "query_neg_match": neg_match,
        "run_id": f"{dataset}/{scheme}/{seed_idx}",
    }
    if scheme in ("mod-minimizer-xor", "mod-minimizer-greedymini"):
        row["mod_t"] = 10
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


def test_validate_schemes_equivalence_ok() -> None:
    rows = [
        _row("cod", "native", 0),
        _row("cod", "native", 1),
        _row("cod", "xor", 0),
        _row("cod", "greedymini", 0, **GM15),
    ]
    n_cmp, errors = validate_schemes_equivalence(rows)
    assert n_cmp == 4
    assert errors == []
    print("ok  validate_schemes_equivalence: no errors on matching runs")


def test_validate_schemes_equivalence_catches_kmer_mismatch() -> None:
    rows = [
        _row("cod", "native", 0, n_kmers=1000),
        _row("cod", "xor", 0, n_kmers=999),
    ]
    _, errors = validate_schemes_equivalence(rows)
    assert any("n_kmers differs" in e for e in errors), errors
    print("ok  validate_schemes_equivalence: catches num_kmers mismatch")


def test_validate_schemes_equivalence_catches_query_mismatch() -> None:
    rows = [
        _row("cod", "native", 0),
        _row("cod", "xor", 0, pos_match=False),
        _row("cod", "greedymini", 0, neg_match=False, **GM15),
        _row("cod", "xor-pcg-hash", 0, check_ok=False),
    ]
    _, errors = validate_schemes_equivalence(rows)
    assert any("query_pos_match is false" in e for e in errors), errors
    assert any("query_neg_match is false" in e for e in errors), errors
    assert any("EVERYTHING OK" in e for e in errors), errors
    print("ok  validate_schemes_equivalence: catches query_pos_match/query_neg_match/check_ok failures")


def test_summary_rows_and_ratios() -> None:
    rows = [
        _row("cod", "native", 0, build_wall_s=10.0, bits_per_kmer=5.0),
        _row("cod", "native", 1, build_wall_s=12.0, bits_per_kmer=5.0),
        _row("cod", "greedymini", 0, build_wall_s=8.0, bits_per_kmer=4.0, **GM15),
        _row("cod", "greedymini", 1, build_wall_s=9.0, bits_per_kmer=4.0, **GM15),
    ]
    summary = summary_rows(rows)
    assert len(summary) == 2
    attach_native_ratios(summary)
    gm = next(r for r in summary if r["scheme"] == "greedymini")
    native = next(r for r in summary if r["scheme"] == "native")
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


def _sample_book_inputs():
    rows = [_row("cod", "native", 0), _row("cod", "native", 1), _row("cod", "greedymini", 0, **GM15),
            _row("cod", "greedymini", 1, **GM15), _row("kestrel", "native", 0),
            _row("kestrel", "greedymini", 0, **GM15)]
    summary = summary_rows(rows)
    attach_native_ratios(summary, "native")
    hist = [{**{f: r.get(f) for f in HISTOGRAM_ID_FIELDS}, "bucket_type": "singleton", "log2_bucket_size": 0,
             "bucket_size_range": "1", "n_buckets": 5, "n_kmers_in_buckets": 5}
            for r in rows[:4]]
    return rows, summary, hist


def test_write_csv_and_xlsx_roundtrip() -> None:
    rows, summary, hist = _sample_book_inputs()
    with tempfile.TemporaryDirectory(prefix="sshash_export_test_") as td:
        out_dir = Path(td)
        write_csv(out_dir / "sshash_runs.csv", rows)
        write_csv(out_dir / "sshash_summary.csv", summary, summary_field_order())
        assert (out_dir / "sshash_runs.csv").is_file()
        assert (out_dir / "sshash_summary.csv").is_file()
        header = (out_dir / "sshash_summary.csv").read_text().splitlines()[0].split(",")
        assert header[4:10] == ["scheme", "gm_k", "projection", "gm_table", "mod_t", "modality"]

        sheets = build_sheets(rows, summary, schemes_rows(summary), bucket_rows(summary),
                              space_rows(rows, summary), hist, [])
        assert [sh["name"] for sh in sheets] == [
            "overview", "schemes", "datasets", "summary", "runs", "buckets", "space_breakdown",
            "bucket_histogram", "super_kmer_histogram"]
        xlsx_path = out_dir / "sshash_results.xlsx"
        names = write_book(xlsx_path, sheets)
        assert names[0] == "overview" and "runs" in names
        assert xlsx_check.check_xlsx(str(xlsx_path)) == []
        with zipfile.ZipFile(xlsx_path) as z:
            assert "xl/charts/chart1.xml" in z.namelist()  # the bucket-histogram charts
            summary_xml = z.read("xl/worksheets/sheet4.xml").decode()
            assert "<f>IFERROR(CORREL(" in summary_xml  # the correlation footer rows
            assert "<mergeCells" in summary_xml
    print("ok  write_csv / write_book: valid workbook with the expected sheets, charts and footers")


def test_layout_paths(monkeypatch=None) -> None:
    import os
    from scripts.common import bench_paths
    from scripts.expected_density import expected_schemes_csv_path
    old = {k: os.environ.pop(k, None) for k in ("SSHASH_DATA_DIR", "SSHASH_RESULTS_DIR")}
    os.environ["SSHASH_BENCH_ROOT"] = "/r"
    try:
        data, results = bench_paths()
        assert (str(data), str(results)) == ("/r/data", "/r/results")
        assert str(expected_schemes_csv_path(results)) == "/r/results/csvs/sshash_density_expected_schemes.csv"
    finally:
        for k, v in old.items():
            if v is not None:
                os.environ[k] = v


def test_schemes_sheet_unified_klw() -> None:
    import schemes_sheet as ssh
    _, summary, _ = _sample_book_inputs()
    schemes = schemes_rows(summary, {})
    fields = [f for f in schemes_field_order(summary)]
    assert fields[:9] == ["scheme", "k", "w", "L", "gm_k", "projection", "gm_table", "mod_t", "description"]
    for r in schemes:
        assert (r["k"], r["w"], r["L"]) == (11, 21, 31), r
    native = next(r for r in schemes if r["scheme"] == "native")
    assert native["cod_particular_density_factor_std"] == 0.0 or native["cod_particular_density_factor_std"] is not None
    # no expected density given -> the helper reports the holes (export dies on them)
    assert ssh.missing_cells(schemes, ["expected_density_factor_mean"])
    print("ok  schemes sheet: k=m, w, L=k+w-1 and column order")


def test_variants_are_separate_groups() -> None:
    """Two tables of one greedymini scheme stay separate summary / space rows; each ratio is
    against native."""
    rows = [_row("cod", "native", 0), _row("cod", "greedymini", 0, **GM13, bits_per_kmer=2.0),
            _row("cod", "greedymini", 0, **GM15, bits_per_kmer=3.0)]
    summary = summary_rows(rows)
    attach_native_ratios(summary, "native")
    assert [(r["scheme"], r["gm_k"]) for r in summary] == [("native", None), ("greedymini", 13), ("greedymini", 15)]
    assert [r["bits_per_kmer_ratio"] for r in summary] == [1.0, 0.5, 0.75]
    assert [(r["scheme"], r["gm_k"]) for r in space_rows(rows, summary)] == [
        ("native", None), ("greedymini", 13), ("greedymini", 15)]
    print("ok  scheme variants (gm_k) are separate groups")


def main() -> int:
    test_validate_schemes_equivalence_ok()
    test_schemes_sheet_unified_klw()
    test_validate_schemes_equivalence_catches_kmer_mismatch()
    test_validate_schemes_equivalence_catches_query_mismatch()
    test_summary_rows_and_ratios()
    test_datasets_csv_rows_are_complete()
    test_write_csv_and_xlsx_roundtrip()
    test_variants_are_separate_groups()
    test_log2_bucket_bins_follow_the_skew_partitions()
    test_correlation_rows_formulas_and_cached_values()
    test_space_rows_saving_split()
    print("all export tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())


def test_runs_and_summary_bands_are_contiguous_and_renamed():
    from scripts.export import RUN_FIELDS, summary_field_order, sheet_categories
    for name, fields in (('runs', list(RUN_FIELDS)), ('summary', summary_field_order())):
        seen, prev = set(), None
        for f, c in zip(fields, sheet_categories(name, fields)):
            if c != prev:
                assert c not in seen, f"band {c!r} fragmented at {f}"
                seen.add(c)
                prev = c
    for gone in ("n_seeds", "query_elapsed_ms", "python_wall_s"):
        assert gone not in RUN_FIELDS and gone not in summary_field_order()
    assert "n_repeats" in summary_field_order() and "query_run_s_mean" in summary_field_order()
    assert "distinct_minimizers" in RUN_FIELDS and "distinct_minimizers_mean" in summary_field_order()
    assert {"build_pinned_cpus", "bench_pinned_cpus", "query_pinned_cpus"} <= set(RUN_FIELDS)
