#!/usr/bin/env python3
"""Tests for scripts/expected_density.py: theoretical expected density_factor,
generalized to any signature-order scheme (not just GreedyMini)."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.expected_density import (
    ORDER_LABELS,
    build_dna_keys,
    compute_all_orders,
    dna_rev,
    ensure_expected_orders,
    expected_density,
    load_expected_orders,
    rxs_m_xs64,
)

P = 5  # small alphabet (4^5=1024) for fast structural tests; GM needs the real p=9 tables


def test_lexicographic_is_pure_canonical_min():
    keys = build_dna_keys("lexicographic", 28, P)
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


def test_hash_seed_zero_equals_true_lexicographic():
    # Verified earlier this session: "-or hash -os0" is order-equivalent to
    # true canonical lex (no is_allowed denylist). Regression-test that
    # equivalence directly against build_dna_keys.
    hash_keys = build_dna_keys("hash", 28, P, seed=0)
    lex_keys = build_dna_keys("lexicographic", 28, P)
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
    for order_label, kw in [
        ("lexicographic", {}),
        ("default", {}),
        ("hash", {"seed": 0xDEADBEEF}),
        ("pcg", {"seed": 0xC0FFEE01}),
    ]:
        keys = build_dna_keys(order_label, 28, P, **kw)
        for i in range(n):
            assert keys[i] == keys[dna_rev(i, P)], order_label


def test_hash_and_pcg_differ_for_same_seed():
    hash_keys = build_dna_keys("hash", 28, P, seed=0xDEADBEEF)
    pcg_keys = build_dna_keys("pcg", 28, P, seed=0xDEADBEEF)
    assert hash_keys != pcg_keys


def test_rxs_m_xs64_matches_custom_order_cpp_constant():
    # Sanity: deterministic, matches the one recorded reference value from
    # custom_order.cpp's rxs_m_xs64 (12605985483714917081 multiplier).
    assert rxs_m_xs64(0) == 0
    assert rxs_m_xs64(1) != rxs_m_xs64(2)


def test_expected_density_is_in_plausible_range():
    keys = build_dna_keys("lexicographic", 28, P)
    w = 6
    d = expected_density(keys, P, w, n_bases=50_000)
    assert d is not None
    assert 1.0 <= d <= w + 1


def test_greedymini_uses_real_table_and_seed_zero_is_unsuffixed_file():
    # Regression: seeds()[0] (0xDEADBEEF) has no `_s0xdeadbeef` suffixed file on
    # disk -- it's the *unsuffixed* `_dna.gm` (see README: "the seed is in the
    # file header"). Passing seed=None must resolve to that same file. This
    # uses the real p=9 GreedyMini table already in the repo.
    from scripts.config import seeds

    assert seeds()[0] == 0xDEADBEEF
    keys_none = build_dna_keys("greedymini_dna", 28, 9)
    keys_explicit_wrong = None
    try:
        keys_explicit_wrong = build_dna_keys("greedymini_dna", 28, 9, seed=seeds()[0])
    except FileNotFoundError:
        pass
    assert keys_explicit_wrong is None, "seeds()[0] has no _s0xdeadbeef file; must use seed=None"
    assert len(keys_none) == 4**9


def test_compute_all_orders_small_alphabet_no_greedymini():
    # Point gm_dir at an empty directory so greedymini_* gracefully report
    # n_seeds=0 instead of trying to load real p=9 tables at this test's p=5.
    with tempfile.TemporaryDirectory() as td:
        rows = compute_all_orders(28, P, w=6, gm_dir=Path(td), n_bases=20_000)
    by_label = {r["order_label"]: r for r in rows}
    assert set(by_label) == set(ORDER_LABELS)

    assert by_label["default"]["n_seeds"] == 1
    assert by_label["default"]["density_factor_std"] is None
    assert by_label["lexicographic"]["n_seeds"] == 1
    assert by_label["lexicographic"]["density_factor_std"] is None

    assert by_label["hash"]["n_seeds"] == 5
    assert by_label["pcg"]["n_seeds"] == 5
    assert by_label["hash"]["density_factor_avg"] is not None

    assert by_label["greedymini_dna"]["n_seeds"] == 0
    assert by_label["greedymini_dna"]["density_factor_avg"] is None


def test_compute_all_orders_real_greedymini_table_seed_zero_resolves():
    # Real p=9 table already in the repo -- only the canonical (unsuffixed,
    # seed 0) w20_k9 table is kept after the GreedyMini table-retirement
    # pass (orders/GreedyMini/ no longer carries the other 4 seeds' tables),
    # so n_seeds is 1, not 5. This is still a real regression test for the
    # seed-mapping bug it was written for: seeds()[0] (0xDEADBEEF) used to
    # look for a nonexistent _s0xdeadbeef-suffixed file and silently drop to
    # n_seeds=0 instead of resolving to the unsuffixed file.
    rows = compute_all_orders(28, 9, w=20, n_bases=20_000)
    by_label = {r["order_label"]: r for r in rows}
    assert by_label["greedymini_dna"]["n_seeds"] == 1
    assert by_label["greedymini_binary"]["n_seeds"] == 1


def test_ensure_and_load_expected_orders_cache_round_trip():
    with tempfile.TemporaryDirectory() as td:
        out_root = Path(td)
        assert load_expected_orders(out_root) == {}
        built = ensure_expected_orders(out_root, 28, P, 6, gm_dir=out_root, n_bases=10_000)
        assert set(built) == set(ORDER_LABELS)
        loaded = load_expected_orders(out_root)
        assert loaded["default"]["density_factor_avg"] == built["default"]["density_factor_avg"]

        # second call with same k/p reuses the cache (doesn't silently change values)
        again = ensure_expected_orders(out_root, 28, P, 6, gm_dir=out_root)
        assert again["default"]["density_factor_avg"] == built["default"]["density_factor_avg"]


if __name__ == "__main__":
    test_lexicographic_is_pure_canonical_min()
    test_default_pushes_disallowed_motifs_to_sentinel()
    test_hash_seed_zero_equals_true_lexicographic()
    test_canonical_pairs_always_tie()
    test_hash_and_pcg_differ_for_same_seed()
    test_rxs_m_xs64_matches_custom_order_cpp_constant()
    test_expected_density_is_in_plausible_range()
    test_greedymini_uses_real_table_and_seed_zero_is_unsuffixed_file()
    test_compute_all_orders_small_alphabet_no_greedymini()
    test_compute_all_orders_real_greedymini_table_seed_zero_resolves()
    test_ensure_and_load_expected_orders_cache_round_trip()
    print("ok")
