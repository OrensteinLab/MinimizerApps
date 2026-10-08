"""Signature scheme labels for benchmark export and charts."""

from __future__ import annotations

import re

SCHEME_ORDER = (
    "default",
    "lexicographic",
    "hash",
    "pcg",
    "greedymini_dna",
    "greedymini_binary",
    # GM_PROJECTION=at1-cg0 (kmc -oat1): the flipped DNA->binary projection of
    # the binary table -- a different order, so its own scheme row.
    "greedymini_binary_at1",
)

SCHEME_CHART_LABEL = {
    "default": "default",
    "lexicographic": "lex",
    "hash": "hash",
    "pcg": "pcg",
    "greedymini_dna": "GM-DNA",
    "greedymini_binary": "GM-bin",
    "greedymini_binary_at1": "GM-bin-at1",
    # legacy aliases from older runs / exports
    "greedymini": "GM-DNA",
    "greedymini_bin": "GM-bin",
}

_LEGACY_GM = frozenset({"greedymini", "greedymini_bin"})

PROJECTION_AT0 = "at0-cg1"  # kmc's default DNA->binary projection
PROJECTION_AT1 = "at1-cg0"  # kmc -oat1 (lib/run.sh: GM_PROJECTION=at1-cg0)
AT1_SUFFIX = "_at1"

# lib/gm.sh gm_label_for_k appends "_at1" before run_schemes' "_s{seed}".
_AT1_LABEL_RE = re.compile(r"_at1(?:_s[0-9a-fx]+)?$", re.IGNORECASE)


def is_at1_projection(row: dict) -> bool:
    """True for a GM run made with GM_PROJECTION=at1-cg0 (kmc -oat1).

    Detected from (in order) an explicit gm_projection field (export.py sets it
    from meta.json's command line), the scheme dir lib/run.sh writes
    (``greedymini_binary_at1``), or the ``-ol`` label lib/gm.sh builds
    (``gm_w20_k9_canonical_at1_s0x...``)."""
    proj = (row.get("gm_projection") or "").strip().lower()
    if proj:
        return proj == PROJECTION_AT1
    run_id = (row.get("run_id") or "").replace("\\", "/").lower()
    if any(part.startswith("greedymini") and part.endswith(AT1_SUFFIX) for part in run_id.split("/")):
        return True
    raw = (row.get("signature_order_label") or "").strip()
    return bool(raw.lower().startswith("gm_") and _AT1_LABEL_RE.search(raw))


def resolve_scheme_label(row: dict) -> str:
    """Map a run row to a display scheme label (``_at1`` suffix on GM schemes
    run with the flipped projection, so they never pool with the at0 rows)."""
    label = _resolve_base_scheme_label(row)
    if label.startswith("greedymini") and not label.endswith(AT1_SUFFIX) and is_at1_projection(row):
        return label + AT1_SUFFIX
    return label


def _resolve_base_scheme_label(row: dict) -> str:
    sig_order = row.get("signature_order") or "default"
    raw = (row.get("signature_order_label") or "").strip()
    mode = (row.get("signature_order_mode") or "").strip().lower()
    gm_lut = (row.get("gm_lut_mode") or "").strip().lower()
    run_id = (row.get("run_id") or "").strip().lower()

    if sig_order == "default" or (not raw and sig_order != "custom"):
        if "greedymini_binary" in run_id:
            return "greedymini_binary"
        if "greedymini_dna" in run_id or "/greedymini/" in run_id.replace("\\", "/"):
            return "greedymini_dna"
        return "default"

    if gm_lut == "binary":
        return "greedymini_binary"
    if gm_lut == "dna":
        return "greedymini_dna"

    # True canonical lex (no is_allowed denylist -- see README) is implemented
    # as "-or hash -os 0" (order-equivalent to MIN(i, revcomp(i))), so its real
    # kmc.json self-reports signature_order_mode == "hash". Must be checked
    # before the mode == "hash" dispatch below, or it would be misclassified
    # as the random-hash scheme. Reported as "lexicographic" (KMC's shipped,
    # is_allowed-protected order is "default" -- see SCHEME_ORDER above).
    if raw == "lexicographic":
        return "lexicographic"

    low = raw.lower()
    if low.startswith("rand_hash") or mode == "hash":
        return "hash"
    if low.startswith("rand_pcg") or mode == "pcg":
        return "pcg"

    if mode == "table_binary":
        return "greedymini_binary"
    if mode in ("table_dna", "table") or low.startswith("gm_"):
        return "greedymini_dna"

    if "greedymini_binary" in run_id:
        return "greedymini_binary"
    if "greedymini_dna" in run_id:
        return "greedymini_dna"

    return raw


def is_greedymini(row: dict) -> bool:
    """True for GreedyMini DNA or binary-LUT runs (either projection)."""
    label = resolve_scheme_label(row)
    return label.startswith("greedymini") or label in _LEGACY_GM
