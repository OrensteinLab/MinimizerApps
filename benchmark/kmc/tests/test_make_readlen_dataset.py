#!/usr/bin/env python3
"""Tests for the read-length-confound derived-dataset builder (see README:
'Read-length confound check')."""

from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.make_readlen_dataset import (
    build,
    max_output_reads_for_target,
    output_paths_for,
    probe_avg_read_len,
)


def _write_fastq_gz(path: Path, n_records: int, read_len: int = 20) -> None:
    """`read_len` must be a multiple of 4 (repeats of "ACGT")."""
    assert read_len % 4 == 0
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = "ACGT" * (read_len // 4)
    qual = "I" * len(seq)
    with gzip.open(path, "wt") as f:
        for i in range(n_records):
            f.write(f"@r{i}\n{seq}\n+\n{qual}\n")


def test_join_concatenates_and_drops_remainder():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=10, read_len=20)
        out = Path(td) / "out.fastq.gz"
        n_in, n_out, total_bases, out_paths = build([p], out, join_n=3, max_output_reads=None)
        assert n_in == 10
        assert n_out == 3  # 10 // 3, remainder of 1 dropped
        assert total_bases == 3 * 3 * 20
        assert out_paths == [out]
        with gzip.open(out, "rt") as f:
            lines = f.read().splitlines()
        assert len(lines) == 3 * 4
        assert len(lines[1]) == 60


def test_join_spans_multiple_source_files():
    with tempfile.TemporaryDirectory() as td:
        a = Path(td) / "a.fastq.gz"
        b = Path(td) / "b.fastq.gz"
        _write_fastq_gz(a, n_records=3, read_len=8)
        _write_fastq_gz(b, n_records=3, read_len=8)
        out = Path(td) / "out.fastq.gz"
        n_in, n_out, total_bases, _out_paths = build([a, b], out, join_n=4, max_output_reads=None)
        assert n_in == 6
        assert n_out == 1  # one pseudo-read spans the a/b boundary
        assert total_bases == 32


def test_max_output_reads_caps_and_stops_early():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=100, read_len=20)
        out = Path(td) / "out.fastq.gz"
        n_in, n_out, total_bases, _out_paths = build([p], out, join_n=1, max_output_reads=5)
        assert n_out == 5
        assert n_in == 5  # stops consuming input as soon as the cap is hit


def test_probe_avg_read_len():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=1000, read_len=100)
        assert probe_avg_read_len([p], n_probe=50) == 100.0


def test_max_output_reads_for_target_matches_actual_build():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=5000, read_len=100)
        target_kmers = 50_000
        join_n, k = 4, 28
        n_reads = max_output_reads_for_target(
            [p], join_n=join_n, k=k, target_kmers=target_kmers, n_probe=1000
        )
        out = Path(td) / "out.fastq.gz"
        _n_in, n_out, total_bases, _out_paths = build(
            [p], out, join_n=join_n, max_output_reads=n_reads
        )
        assert n_out == n_reads
        approx_kmers = total_bases - n_out * (k - 1)
        # probe-based estimate should land within ~5% of the real target
        assert abs(approx_kmers - target_kmers) / target_kmers < 0.05


def test_join_1_is_passthrough_subsample():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=4, read_len=20)
        out = Path(td) / "out.fastq.gz"
        n_in, n_out, total_bases, _out_paths = build([p], out, join_n=1, max_output_reads=None)
        assert n_in == n_out == 4
        assert total_bases == 4 * 20


def test_output_paths_for_single_file():
    out = Path("/tmp/x/h_sapiens_2_shortsub.fastq.gz")
    assert output_paths_for(out, 1) == [out]


def test_output_paths_for_multi_file_naming():
    out = Path("/tmp/x/h_sapiens_2_shortsub.fastq.gz")
    paths = output_paths_for(out, 11)
    assert len(paths) == 11
    assert [p.name for p in paths] == [
        f"h_sapiens_2_shortsub_part{i:02d}.fastq.gz" for i in range(11)
    ]


def test_num_files_round_robins_and_preserves_total():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=23, read_len=20)
        out = Path(td) / "out.fastq.gz"
        n_in, n_out, total_bases, out_paths = build(
            [p], out, join_n=1, max_output_reads=None, num_files=4
        )
        assert n_in == n_out == 23
        assert len(out_paths) == 4
        assert all(p.is_file() for p in out_paths)

        counts = []
        total_reads_seen = 0
        for part in out_paths:
            with gzip.open(part, "rt") as f:
                lines = f.read().splitlines()
            n_records = len(lines) // 4
            counts.append(n_records)
            total_reads_seen += n_records
        assert total_reads_seen == 23  # every read landed in exactly one file
        # round-robin over 23 reads / 4 files -> sizes differ by at most 1
        assert max(counts) - min(counts) <= 1


def test_num_files_matches_single_file_when_1():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "in.fastq.gz"
        _write_fastq_gz(p, n_records=10, read_len=20)
        out = Path(td) / "out.fastq.gz"
        n_in, n_out, total_bases, out_paths = build(
            [p], out, join_n=1, max_output_reads=None, num_files=1
        )
        assert out_paths == [out]
        assert out.is_file()


if __name__ == "__main__":
    test_join_concatenates_and_drops_remainder()
    test_join_spans_multiple_source_files()
    test_max_output_reads_caps_and_stops_early()
    test_probe_avg_read_len()
    test_max_output_reads_for_target_matches_actual_build()
    test_join_1_is_passthrough_subsample()
    test_output_paths_for_single_file()
    test_output_paths_for_multi_file_naming()
    test_num_files_round_robins_and_preserves_total()
    test_num_files_matches_single_file_when_1()
    print("ok")
