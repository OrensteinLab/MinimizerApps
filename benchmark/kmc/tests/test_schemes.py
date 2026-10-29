#!/usr/bin/env python3
"""Tests for scheme identity (canonical names + separate gm columns)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.schemes import (
    SCHEME_ORDER,
    chart_label,
    is_at1_projection,
    is_greedymini,
    scheme_from_json,
    resolve_scheme_label,
    variant_key,
    variant_sort_key,
)


def _gm(**extra):
    row = {
        "scheme": "greedymini",
        "gm_k": 9,
        "projection": "at0-cg1",
        "gm_table": "w20_k9_canonical",
        "gm_lut": "binary",
    }
    row.update(extra)
    return row


def test_scheme_order_is_canonical_names_only():
    assert SCHEME_ORDER == ("default", "lex", "xor", "xor-pcg-hash", "greedymini")


def test_resolve_is_the_scheme_field():
    for scheme in SCHEME_ORDER:
        assert resolve_scheme_label({"scheme": scheme}) == scheme
    assert resolve_scheme_label({}) == "default"


@pytest.mark.parametrize("old", ["hash", "pcg", "lexicographic", "greedymini_dna", "greedymini_binary", "greedymini_binary_at1"])
def test_old_spellings_are_rejected(old):
    with pytest.raises(ValueError):
        resolve_scheme_label({"scheme": old})


def test_is_greedymini():
    assert is_greedymini(_gm())
    assert is_greedymini(_gm(projection="at1-cg0", gm_lut="dna"))
    assert not is_greedymini({"scheme": "default"})
    assert not is_greedymini({"scheme": "xor"})
    assert not is_greedymini({"scheme": "lex"})


def test_scheme_from_json():
    assert scheme_from_json("default", "", "") == "default"
    # true lex is "-orxor -os0 -ollex": kmc.json reports mode xor, label lex
    assert scheme_from_json("custom", "xor", "lex") == "lex"
    assert scheme_from_json("custom", "xor", "xor_s0x1") == "xor"
    assert scheme_from_json("custom", "xor-pcg-hash", "xor-pcg-hash_s0x1") == "xor-pcg-hash"
    assert scheme_from_json("custom", "greedymini", "gm_w20_k9_canonical_s0x1") == "greedymini"
    for old in ("hash", "pcg", "table", "table_dna", "table_binary"):
        with pytest.raises(ValueError):
            scheme_from_json("custom", old, "x")


def test_variants_differ_by_every_gm_column():
    base = _gm()
    keys = {
        variant_key(base),
        variant_key(_gm(projection="at1-cg0")),
        variant_key(_gm(gm_lut="dna")),
        variant_key(_gm(gm_table="w20_k9_other")),
        variant_key(_gm(gm_k=10)),
    }
    assert len(keys) == 5
    # non-greedymini schemes ignore stray gm fields
    assert variant_key({"scheme": "xor", "gm_k": 9, "projection": "at1-cg0"}) == ("xor", None, None, None, None)


def test_at1_projection_is_its_own_variant():
    a0, a1 = _gm(), _gm(projection="at1-cg0")
    assert not is_at1_projection(a0) and is_at1_projection(a1)
    assert variant_sort_key(a0) != variant_sort_key(a1)
    assert resolve_scheme_label(a0) == resolve_scheme_label(a1) == "greedymini"


def test_chart_labels():
    assert chart_label({"scheme": "xor-pcg-hash"}) == "xor-pcg-hash"
    assert chart_label(_gm(gm_lut="dna")) == "greedymini (dna)"
    assert chart_label(_gm(projection="at1-cg0")) == "greedymini (binary, at1-cg0)"


def test_variant_sort_follows_scheme_order():
    rows = [_gm(), {"scheme": "xor"}, {"scheme": "default"}, {"scheme": "lex"}, {"scheme": "xor-pcg-hash"}]
    assert [r["scheme"] for r in sorted(rows, key=variant_sort_key)] == list(SCHEME_ORDER)
