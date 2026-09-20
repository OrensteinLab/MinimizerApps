"""Low-level .xlsx/.csv writer primitives shared by kraken2's xlsxutil.py and
kmc's scripts/export.py + scripts/density.py. Each app's own writer builds
very different worksheets (kraken2 adds charts, kmc adds styles/merges/
conditional formatting) -- only the genuinely app-independent pieces (cell
reference math, XML escaping, the fixed root .rels part, and CSV writing)
live here; verified byte-identical to both apps' prior independent
implementations before sharing.
"""

from __future__ import annotations

import csv
import os
from xml.sax.saxutils import escape as _sax_escape


def col_letter(n: int) -> str:
    """1-based column index -> Excel column letter ("A", "B", ..., "AA", ...)."""
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def xml_escape(text) -> str:
    """Escape &, <, >, and " for use in XML text/attribute content."""
    return _sax_escape(str(text), {'"': "&quot;"})


def root_rels_xml() -> str:
    """The package-level `_rels/.rels` part: identical for any single-workbook .xlsx."""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )


def write_csv(path, fieldnames, rows) -> None:
    """Write `rows` (dicts) as CSV with `fieldnames` as the header/column
    order. Always creates the parent directory and always writes a header,
    even for zero rows; missing/None values become "" ."""
    os.makedirs(os.path.dirname(str(path)) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fieldnames), extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})
