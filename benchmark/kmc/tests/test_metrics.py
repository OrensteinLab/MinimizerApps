#!/usr/bin/env python3
"""Tests for super-k-mer / minimizer density metrics."""

from __future__ import annotations

import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.metrics import super_kmer_metrics


def test_density_factor_lex_f_vesca():
    # f_vesca k=28 p=9 lex run (n_tot, n_sup from kmc.json)
    m = super_kmer_metrics(
        k=28,
        signature_len=9,
        n_total_kmers=4_134_079_185,
        n_super_kmers=393_866_689,
    )
    assert m["minimizer_window_w"] == 20
    assert abs(m["density_factor"] - 2.0) < 0.01
    assert abs(m["avg_super_kmer_len"] - 10.496) < 0.01
    assert m["avg_read_len"] is None


def test_avg_read_len_f_vesca():
    m = super_kmer_metrics(
        k=28,
        signature_len=9,
        n_total_kmers=4_134_079_185,
        n_super_kmers=393_866_689,
        n_reads=12_803_137,
    )
    assert abs(m["avg_read_len"] - 349.9) < 0.1


def test_gm_lower_density_factor():
    m = super_kmer_metrics(
        k=28,
        signature_len=9,
        n_total_kmers=4_134_079_185,
        n_super_kmers=352_943_222,
    )
    assert 1.0 < m["density_factor"] < 2.0


if __name__ == "__main__":
    test_density_factor_lex_f_vesca()
    test_avg_read_len_f_vesca()
    test_gm_lower_density_factor()
    print("ok")
