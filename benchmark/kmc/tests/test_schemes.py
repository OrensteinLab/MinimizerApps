#!/usr/bin/env python3
"""Tests for scheme label resolution."""

from __future__ import annotations

import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.schemes import is_greedymini, resolve_scheme_label


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


if __name__ == "__main__":
    test_legacy_greedymini_dir()
    test_binary_lut_mode()
    test_lex_not_swallowed_by_hash()
    test_is_greedymini()
    print("ok")
