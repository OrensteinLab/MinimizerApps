"""Signature-order identity for benchmark export and charts.

kmc uses the canonical order names of benchmark/common/order_names.py:
``xor``, ``xor-pcg-hash``, ``greedymini`` plus its two baselines ``default``
(KMC's shipped is_allowed-protected lex) and ``lex`` (true canonical lex, run as
``-orxor -os0``). kmc has no mod-minimizer.

A run's identity is a set of separate columns, never packed into one string:

  order       canonical name or baseline (``default`` / ``lex``)
  gm_k        GreedyMini table k = the signature length p      (greedymini only)
  projection  ``at0-cg1`` | ``at1-cg0`` (kmc ``-oat1``)         (greedymini only)
  gm_table    table stem, e.g. ``w20_k9_canonical``             (greedymini only)
  gm_lut      table LUT backend, ``dna`` | ``binary``           (greedymini only;
              kmc-specific: same order, different storage -- kept as its own column)
"""

from __future__ import annotations

import order_names

BASELINES = ("default", "lex")
SCHEME_ORDER = ("default", "lex", "xor", "xor-pcg-hash", "greedymini")
GM_LUTS = ("dna", "binary")
GM_ID_FIELDS = ("gm_k", "projection", "gm_table", "gm_lut")
ID_FIELDS = ("order", *GM_ID_FIELDS)

PROJECTION_AT0 = "at0-cg1"  # kmc's default DNA->binary projection
PROJECTION_AT1 = "at1-cg0"  # kmc -oat1


def order_from_json(signature_order: str, mode: str, label: str) -> str:
    """Order of a run from kmc.json's self-report (when meta.json has no
    ``order``): ``default``, ``lex`` (true lex = ``-orxor -os0 -ollex``), or
    the canonical mode name (xor / xor-pcg-hash / greedymini)."""
    if (signature_order or "default") == "default":
        return "default"
    if (label or "").strip() == "lex":
        return "lex"
    return order_names.validate_order((mode or "").strip(), BASELINES)


def resolve_scheme_label(row: dict) -> str:
    """The run's order (canonical name or baseline)."""
    return order_names.validate_order((row.get("order") or "default").strip(), BASELINES)


def is_greedymini(row: dict) -> bool:
    return order_names.is_greedymini(resolve_scheme_label(row))


def is_at1_projection(row: dict) -> bool:
    return (row.get("projection") or "").strip().lower() == PROJECTION_AT1


def variant_key(row: dict) -> tuple:
    """(order, gm_k, projection, gm_table, gm_lut): one row per variant in the
    orders / summary sheets; non-greedymini orders carry None for the gm fields."""
    order = resolve_scheme_label(row)
    if not order_names.is_greedymini(order):
        return (order, None, None, None, None)
    gk = row.get("gm_k")
    return (
        order,
        int(gk) if gk not in (None, "") else None,
        row.get("projection") or None,
        row.get("gm_table") or None,
        row.get("gm_lut") or None,
    )


def variant_sort_key(row: dict) -> tuple:
    key = variant_key(row)
    try:
        rank = SCHEME_ORDER.index(key[0])
    except ValueError:
        rank = len(SCHEME_ORDER)
    return (rank, key[1] or 0, key[2] or "", key[3] or "", key[4] or "")


def chart_label(row: dict) -> str:
    """Chart/legend label: the order name; greedymini adds its LUT backend and,
    when flipped, the projection (e.g. ``greedymini (binary, at1-cg0)``)."""
    order, _gk, proj, _tab, lut = variant_key(row)
    if order != "greedymini":
        return order
    extra = [x for x in (lut, proj if proj == PROJECTION_AT1 else None) if x]
    return "greedymini (%s)" % ", ".join(extra) if extra else "greedymini"
