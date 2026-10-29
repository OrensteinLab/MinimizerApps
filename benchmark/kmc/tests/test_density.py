#!/usr/bin/env python3
"""Tests for GreedyMini expected / particular density."""

from __future__ import annotations

import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.density import (
    _sliding_runs,
    attach_density_refs,
    context_is_charged,
    density_factor,
    expected_binary_density,
    expected_binary_density_stream,
    expected_dna_density,
    particular_on_bases,
    resolve_gm_ranks,
)
from scripts.metrics import super_kmer_metrics


def _ranks():
    p, ranks, src = resolve_gm_ranks(28, 9)
    assert p == 9
    assert len(ranks) == 512
    assert src.exists()
    return p, ranks


def test_density_factor_formula():
    assert density_factor(100, 1050, 20) == round(100 / 1050 * 21, 4)


def test_sliding_runs_constant():
    keys = [5] * 30
    n_win, n_super, n_sel = _sliding_runs(keys, w=20)
    assert n_win == 11
    assert n_super == 1
    # left tie-break slides the argmin index; signature value stays the same
    assert n_sel == n_win


def test_sliding_runs_all_distinct_increasing():
    keys = list(range(25))
    n_win, n_super, n_sel = _sliding_runs(keys, w=5)
    assert n_win == 21
    # left-most is always the min in an increasing window, so the min slides every step
    assert n_super == n_win
    assert n_sel == n_win


def test_charged_changes_when_min_slides():
    p, w = 3, 2
    ranks = list(range(8))
    # increasing k-mers 000,001,010: first window min at 0, second at 1 → charged
    # bits: 000 001 010 → 000001010 = 0b000001010
    assert context_is_charged(0b000001010, ranks, p, w)
    # constant 00000: both windows min at leftmost surviving index 0 vs 1 → charged
    # (left tie-break slides). Density of homopolymers is 1.
    assert context_is_charged(0, ranks, p, w)


def test_expected_binary_near_1_72():
    p, ranks = _ranks()
    w = 20
    got = expected_binary_density(ranks, p, w, n_samples=400_000, seed=1)
    df = got["density_factor"]
    assert 1.62 < df < 1.82, df
    stream = expected_binary_density_stream(ranks, p, w, n_bits=800_000, seed=1)
    assert abs(stream["density_factor_positions"] - df) < 0.05


def test_expected_dna_kmc_not_below_binary():
    p, ranks = _ranks()
    w = 20
    binary = expected_binary_density(ranks, p, w, n_samples=200_000, seed=2)
    dna = expected_dna_density(ranks, p, w, n_bases=400_000, seed=2)
    # Lossy AT/GC projection cannot beat the binary order on random input by much.
    assert dna["density_factor"] >= binary["density_factor"] - 0.08
    assert dna["density_factor"] < 2.15


def test_particular_dna_vs_binary_on_repeat():
    p, ranks = _ranks()
    from scripts.density import _cmp_keys_binary, _cmp_keys_dna

    dna_keys = _cmp_keys_dna(ranks, p, 0)
    bin_keys = _cmp_keys_binary(ranks, p)
    # poly-A: binary is all zeros; DNA is all zeros. Super-k-mer run is 1.
    bases = [0] * 200
    acc = particular_on_bases(bases, p=p, w=20, dna_keys=dna_keys, bin_keys=bin_keys)
    assert acc["n_super_dna"] == 1
    assert acc["n_super_bin"] == 1
    assert acc["n_kmers"] == 200 - 28 + 1 or acc["n_kmers"] == 200 - 9 - 20 + 2
    # windows of w signatures: n = (L-p+1) - w + 1 = L - (p+w-1) + 1 = L - 28 + 1
    assert acc["n_kmers"] == 200 - 28 + 1


def test_attach_density_refs():
    part = {
        ("f_vesca", 28, 9): {
            "density_factor_binary_positions": "1.84",
            "density_factor_dna_positions": "1.82",
        }
    }
    gm = {
        "dataset": "f_vesca",
        "k": 28,
        "p": 9,
        "density_factor": 1.79,
        "scheme": "greedymini",
        "gm_k": 9,
        "projection": "at0-cg1",
        "gm_lut": "dna",
        "signature_order": "custom",
        "signature_order_mode": "greedymini",
        "signature_order_label": "gm_w20_k9",
    }
    lex = {
        "dataset": "f_vesca",
        "k": 28,
        "p": 9,
        "density_factor": 2.0,
        "scheme": "default",
        "signature_order": "default",
    }
    attach_density_refs([gm, lex], particular=part)
    assert gm["binary_density_factor"] == 1.84
    assert lex["binary_density_factor"] is None

    gm_bin = {
        **gm,
        "gm_lut": "binary",
        "signature_order_label": "gm_bin",
        "density_factor": 1.84,
    }
    attach_density_refs([gm_bin], particular=part)
    assert gm_bin["binary_density_factor"] is None


def test_kmc_formula_still_2_for_lex_vesca():
    m = super_kmer_metrics(
        k=28,
        signature_len=9,
        n_total_kmers=4_134_079_185,
        n_super_kmers=393_866_689,
    )
    assert abs(m["density_factor"] - 2.0) < 0.01


if __name__ == "__main__":
    test_density_factor_formula()
    test_sliding_runs_constant()
    test_sliding_runs_all_distinct_increasing()
    test_charged_changes_when_min_slides()
    test_expected_binary_near_1_72()
    test_expected_dna_kmc_not_below_binary()
    test_particular_dna_vs_binary_on_repeat()
    test_attach_density_refs()
    test_kmc_formula_still_2_for_lex_vesca()
    print("ok")
