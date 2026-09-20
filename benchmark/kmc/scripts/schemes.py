"""Signature scheme labels for benchmark export and charts."""

from __future__ import annotations

SCHEME_ORDER = (
    "default",
    "lexicographic",
    "hash",
    "pcg",
    "greedymini_dna",
    "greedymini_binary",
)

SCHEME_CHART_LABEL = {
    "default": "default",
    "lexicographic": "lex",
    "hash": "hash",
    "pcg": "pcg",
    "greedymini_dna": "GM-DNA",
    "greedymini_binary": "GM-bin",
    # legacy aliases from older runs / exports
    "greedymini": "GM-DNA",
    "greedymini_bin": "GM-bin",
}

_LEGACY_GM = frozenset({"greedymini", "greedymini_bin"})


def resolve_scheme_label(row: dict) -> str:
    """Map a run row to a display scheme label."""
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
    """True for GreedyMini DNA or binary-LUT runs."""
    label = resolve_scheme_label(row)
    return label.startswith("greedymini") or label in _LEGACY_GM
