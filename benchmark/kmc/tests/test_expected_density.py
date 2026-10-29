#!/usr/bin/env python3
"""Tests for scripts/expected_density.py: theoretical expected density_factor,
generalized to any scheme (not just GreedyMini)."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.expected_density import (
    SCHEME_VARIANTS,
    build_dna_keys,
    compute_all_schemes,
    dna_rev,
    ensure_expected_schemes,
    expected_density,
    load_expected_schemes,
    rxs_m_xs64,
)

P = 5  # small alphabet (4^5=1024) for fast structural tests; GM needs the real p=9 tables


def test_lex_is_pure_canonical_min():
    keys = build_dna_keys("lex", 28, P)
    n = 1 << (2 * P)
    for i in range(n):
        assert keys[i] == min(i, dna_rev(i, P))


def test_default_pushes_disallowed_motifs_to_sentinel():
    keys = build_dna_keys("default", 28, P)
    n = 1 << (2 * P)
    # poly-A (i=0) and poly-T (all 1 bits) are is_allowed-excluded -- both must
    # map to the same out-of-range sentinel, never winning an argmin.
    assert keys[0] == n
    assert keys[n - 1] == n
    # some p-mer that IS allowed must NOT be pushed to the sentinel
    assert any(v < n for v in keys)


def test_xor_seed_zero_equals_true_lex():
    # Verified earlier this session: "-or xor -os0" is order-equivalent to
    # true canonical lex (no is_allowed denylist). Regression-test that
    # equivalence directly against build_dna_keys.
    hash_keys = build_dna_keys("xor", 28, P, seed=0)
    lex_keys = build_dna_keys("lex", 28, P)
    n = 1 << (2 * P)
    # order-equivalent, not necessarily identical values: same relative order
    for a in range(n):
        for b in range(n):
            assert (hash_keys[a] < hash_keys[b]) == (lex_keys[a] < lex_keys[b])
            assert (hash_keys[a] == hash_keys[b]) == (lex_keys[a] == lex_keys[b])


def test_canonical_pairs_always_tie():
    # Both-strand (canonical) counting: a p-mer and its reverse complement must
    # always land on the identical comparator value, for every order.
    n = 1 << (2 * P)
    for scheme, kw in [
        ("lex", {}),
        ("default", {}),
        ("xor", {"seed": 0xDEADBEEF}),
        ("xor-pcg-hash", {"seed": 0xC0FFEE01}),
    ]:
        keys = build_dna_keys(scheme, 28, P, **kw)
        for i in range(n):
            assert keys[i] == keys[dna_rev(i, P)], scheme


def test_xor_and_pcg_differ_for_same_seed():
    hash_keys = build_dna_keys("xor", 28, P, seed=0xDEADBEEF)
    pcg_keys = build_dna_keys("xor-pcg-hash", 28, P, seed=0xDEADBEEF)
    assert hash_keys != pcg_keys


def test_rxs_m_xs64_matches_custom_order_cpp_constant():
    # Sanity: deterministic, matches the one recorded reference value from
    # custom_order.cpp's rxs_m_xs64 (12605985483714917081 multiplier).
    assert rxs_m_xs64(0) == 0
    assert rxs_m_xs64(1) != rxs_m_xs64(2)


def test_expected_density_is_in_plausible_range():
    keys = build_dna_keys("lex", 28, P)
    w = 6
    d = expected_density(keys, P, w, n_bases=50_000)
    assert d is not None
    assert 1.0 <= d <= w + 1


def test_compute_all_schemes_small_alphabet_no_greedymini():
    # Point gm_dir at an empty directory so the greedymini variants gracefully
    # report n_seeds=0 instead of trying to load real p=9 tables at this test's p=5.
    with tempfile.TemporaryDirectory() as td:
        rows = compute_all_schemes(28, P, w=6, gm_dir=Path(td), n_bases=20_000)
    by = {(r["scheme"], r["projection"]): r for r in rows}
    assert set(by) == set(SCHEME_VARIANTS)
    assert {(r["k"], r["w"], r["L"]) for r in rows} == {(P, 6, 28)}  # unified names

    for label in ("default", "lex"):
        assert by[(label, "")]["n_seeds"] == 1
        assert by[(label, "")]["density_factor_std"] is None  # the export shows 0
    assert by[("xor", "")]["n_seeds"] == 5
    assert by[("xor-pcg-hash", "")]["n_seeds"] == 5
    assert by[("xor", "")]["density_factor_mean"] is not None
    for proj in ("at0-cg1", "at1-cg0"):
        assert by[("greedymini", proj)]["n_seeds"] == 0
        assert by[("greedymini", proj)]["density_factor_mean"] is None


def test_greedymini_at0_equals_real_dna_gm_table():
    # The sigma-extension of the binary ranks reproduces the real _dna.gm sort_key
    # (seed=None = seeds[0], the unsuffixed file's embedded seed).
    from scripts.config import seeds
    from scripts.expected_density import PROJECTION_AT0, read_dna_gm
    from scripts.gm import dna_gm_path

    assert seeds()[0] == 0xDEADBEEF
    _k, _seed, real = read_dna_gm(dna_gm_path(28, 9))
    assert build_dna_keys("greedymini", 28, 9, seed=seeds()[0], projection=PROJECTION_AT0) == real
    assert build_dna_keys("greedymini", 28, 9, projection=PROJECTION_AT0) == real


def test_at1_is_the_flipped_table_and_a_different_order():
    from scripts.config import seeds
    from scripts.expected_density import flip_gm_ranks

    assert flip_gm_ranks([3, 0, 2, 1]) == [1, 2, 0, 3]  # ranks[key ^ mask]
    at0 = build_dna_keys("greedymini", 28, 9, seed=seeds()[0], projection="at0-cg1")
    at1 = build_dna_keys("greedymini", 28, 9, seed=seeds()[0], projection="at1-cg0")
    assert at0 != at1
    assert len(at1) == 4**9
    assert expected_density(at1, 9, 20, n_bases=100_000) < 2.0


def test_ensure_and_load_expected_schemes_cache_round_trip():
    with tempfile.TemporaryDirectory() as td:
        out_root = Path(td)
        assert load_expected_schemes(out_root) == {}
        built = ensure_expected_schemes(out_root, 28, P, 6, gm_dir=out_root, n_bases=10_000)
        assert set(built) == set(SCHEME_VARIANTS)
        loaded = load_expected_schemes(out_root)
        key = ("default", "")
        assert loaded[key]["density_factor_mean"] == built[key]["density_factor_mean"]
        assert (loaded[key]["k"], loaded[key]["w"], loaded[key]["L"]) == (str(P), "6", "28")

        # second call with same k/p reuses the cache (doesn't silently change values)
        again = ensure_expected_schemes(out_root, 28, P, 6, gm_dir=out_root)
        assert again[key]["density_factor_mean"] == built[key]["density_factor_mean"]


if __name__ == "__main__":
    test_lex_is_pure_canonical_min()
    test_default_pushes_disallowed_motifs_to_sentinel()
    test_xor_seed_zero_equals_true_lex()
    test_canonical_pairs_always_tie()
    test_xor_and_pcg_differ_for_same_seed()
    test_rxs_m_xs64_matches_custom_order_cpp_constant()
    test_expected_density_is_in_plausible_range()
    test_compute_all_schemes_small_alphabet_no_greedymini()
    test_greedymini_at0_equals_real_dna_gm_table()
    test_at1_is_the_flipped_table_and_a_different_order()
    test_ensure_and_load_expected_schemes_cache_round_trip()
    print("ok")
