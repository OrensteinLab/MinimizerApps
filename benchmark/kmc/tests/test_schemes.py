#!/usr/bin/env python3
"""Tests for scheme label resolution."""

from __future__ import annotations

import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.schemes import SCHEME_ORDER, is_greedymini, resolve_scheme_label


def test_legacy_greedymini_dir():
    row = {
        "signature_order": "custom",
        "signature_order_label": "gm_w20_k9_canonical_s0xdeadbeef",
        "signature_order_mode": "table",
        "run_id": "f_vesca/k28_p9/gzip_n512_t12_m12/greedymini/0",
    }
    assert resolve_scheme_label(row) == "greedymini_dna"


def test_binary_lut_mode():
    row = {
        "signature_order": "custom",
        "signature_order_label": "gm_w20_k9_canonical_s0x16180339",
        "signature_order_mode": "table_binary",
        "gm_lut_mode": "binary",
        "run_id": "f_vesca/k28_p9/gm_lut_binary/greedymini_binary/0",
    }
    assert resolve_scheme_label(row) == "greedymini_binary"


def test_lex_not_swallowed_by_hash():
    # Real kmc.json for a lex run self-reports signature_order_mode ==
    # "hash" (it's plumbed through -or hash under the hood -- see
    # lib/run.sh) -- must resolve to its own label ("lexicographic": true
    # canonical lex, no is_allowed denylist -- distinct from "default", which
    # is what KMC actually ships), not "hash".
    row = {
        "signature_order": "custom",
        "signature_order_label": "lexicographic",
        "signature_order_mode": "hash",
        "run_id": "f_vesca/k28_p9/gzip_n512_t12_m12/lexicographic/0",
    }
    assert resolve_scheme_label(row) == "lexicographic"
    assert not is_greedymini(row)


def test_is_greedymini():
    assert is_greedymini(
        {
            "signature_order": "custom",
            "signature_order_mode": "table_dna",
            "signature_order_label": "gm_w20_k9",
        }
    )
    assert is_greedymini(
        {
            "signature_order": "custom",
            "signature_order_mode": "table_binary",
            "gm_lut_mode": "binary",
        }
    )
    assert not is_greedymini({"signature_order": "default"})
    assert not is_greedymini(
        {
            "signature_order": "custom",
            "signature_order_mode": "hash",
            "signature_order_label": "rand_hash",
        }
    )


def _at1_row(**extra):
    row = {
        "signature_order": "custom",
        "signature_order_label": "gm_w20_k9_canonical_at1_s0xdeadbeef",
        "signature_order_mode": "table_binary",
        "gm_lut_mode": "binary",
        "run_id": "f_vesca/k28_p9/gzip_n512_t12_m12/greedymini_binary_at1/0",
    }
    row.update(extra)
    return row


def test_at1_projection_is_its_own_scheme():
    # GM_PROJECTION=at1-cg0 (kmc -oat1) is a different order: it must never
    # resolve to (and so pool/aggregate with) the at0 greedymini_binary row.
    assert resolve_scheme_label(_at1_row()) == "greedymini_binary_at1"
    # any single signal is enough: scheme dir, label suffix, explicit field
    assert resolve_scheme_label(_at1_row(signature_order_label="gm_w20_k9_canonical_s0x1")) == (
        "greedymini_binary_at1"
    )
    assert resolve_scheme_label(_at1_row(run_id="f_vesca/x/greedymini_binary/0")) == "greedymini_binary_at1"
    assert resolve_scheme_label(
        _at1_row(run_id="a/b", signature_order_label="gm_x_s0x1", gm_projection="at1-cg0")
    ) == "greedymini_binary_at1"
    # explicit at0 wins over a stale heuristic
    assert resolve_scheme_label(_at1_row(gm_projection="at0-cg1")) == "greedymini_binary"
    assert is_greedymini(_at1_row())
    assert "greedymini_binary_at1" in SCHEME_ORDER


def test_at0_binary_unaffected_by_at1_detection():
    row = _at1_row(
        signature_order_label="gm_w20_k9_canonical_s0xdeadbeef",
        run_id="f_vesca/k28_p9/gzip_n512_t12_m12/greedymini_binary/0",
    )
    assert resolve_scheme_label(row) == "greedymini_binary"
    # non-GM labels never get the suffix
    assert resolve_scheme_label(
        {"signature_order": "custom", "signature_order_label": "rand_hash_s0x1", "gm_projection": "at1-cg0"}
    ) == "hash"


if __name__ == "__main__":
    test_legacy_greedymini_dir()
    test_at1_projection_is_its_own_scheme()
    test_at0_binary_unaffected_by_at1_detection()
    test_binary_lut_mode()
    test_lex_not_swallowed_by_hash()
    test_is_greedymini()
    print("ok")
