"""Low-level .xlsx/.csv primitives shared by the kraken2, kmc and sshash exports.

Cell reference math, XML escaping, the fixed root .rels part and CSV writing
are used by every app's own sheet writer (sshash/kmc build their worksheets
in their own export.py). The minimal .xlsx *reader* (load_xlsx) and the
CSV/TSV table loader (load_table) are used by kraken2; the styled
multi-sheet writer used by kraken2 lives in xlsx_book.py.

Python 3.6 compatible (kraken2's run_benchmark.sh calls `python3`).
"""

from __future__ import print_function

import csv
import os
import zipfile
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape as _sax_escape


def col_letter(n):
    """1-based column index -> Excel column letter ("A", "B", ..., "AA", ...)."""
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def xml_escape(text):
    """Escape &, <, >, and " for use in XML text/attribute content."""
    return _sax_escape(str(text), {'"': "&quot;"})


def root_rels_xml():
    """The package-level `_rels/.rels` part: identical for any single-workbook .xlsx."""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )


def write_csv(path, fieldnames, rows):
    """Write `rows` (dicts) as CSV with `fieldnames` as the header/column
    order. Always creates the parent directory and always writes a header,
    even for zero rows; missing/None values become "" ."""
    os.makedirs(os.path.dirname(str(path)) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fieldnames), extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})


# ---------------------------------------------------------------------------
# Reader: minimal .xlsx -> {sheet: {row: {col: text}}}; CSV/TSV table loader
# ---------------------------------------------------------------------------

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def _col_row(ref):
    col = "".join(c for c in ref if c.isalpha())
    row = int("".join(c for c in ref if c.isdigit()))
    n = 0
    for c in col:
        n = n * 26 + (ord(c) - 64)
    return n, row


def _shared_strings(zf):
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    out = []
    for si in root.findall("m:si", NS):
        out.append("".join(t.text or "" for t in si.findall(".//m:t", NS)))
    return out


def _cell_value(cell, shared):
    t = cell.attrib.get("t")
    v = cell.find("m:v", NS)
    if v is None or v.text is None:
        is_el = cell.find("m:is", NS)
        if is_el is not None:
            return "".join(t.text or "" for t in is_el.findall(".//m:t", NS))
        return ""
    if t == "s":
        return shared[int(v.text)]
    return v.text


def workbook_sheets(zf):
    root = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rid_to_target = {}
    for rel in rels:
        rid = rel.attrib.get("Id")
        target = rel.attrib.get("Target")
        if rid and target:
            rid_to_target[rid] = target
    sheets = []
    for sh in root.findall("m:sheets/m:sheet", NS):
        name = sh.attrib.get("name")
        rid = sh.attrib.get(REL)
        target = rid_to_target.get(rid, "")
        if target.startswith("/"):
            path = target.lstrip("/")
        else:
            path = "xl/" + target.lstrip("/")
            if not path.startswith("xl/"):
                path = "xl/" + target
        sheets.append((name, path))
    return sheets


def sheet_rows(zf, sheet_path, shared):
    root = ET.fromstring(zf.read(sheet_path))
    rows = {}
    for c in root.findall(".//m:c", NS):
        ref = c.attrib.get("r")
        if not ref:
            continue
        col, row = _col_row(ref)
        rows.setdefault(row, {})[col] = _cell_value(c, shared)
    return rows


def row_values(cols, min_col=1, max_col=None):
    if not cols:
        return []
    if max_col is None:
        max_col = max(cols)
    return [cols.get(i, "") for i in range(min_col, max_col + 1)]


def load_xlsx(path):
    with zipfile.ZipFile(path) as zf:
        shared = _shared_strings(zf)
        sheets = {}
        for name, spath in workbook_sheets(zf):
            sheets[name] = sheet_rows(zf, spath, shared)
        return sheets


# --- writer ---


def _coerce_cell(text):
    if text is None:
        return ""
    s = str(text).strip()
    if s == "":
        return ""
    try:
        if s.startswith("0") and len(s) > 1 and "." not in s.lower() and "e" not in s.lower():
            return s
        if any(c in s for c in ".eE"):
            return float(s)
        return int(s)
    except ValueError:
        return text


def load_table(path):
    """Read a .csv/.tsv (delimiter sniffed) -> (fieldnames, [row dicts]) with
    numeric-looking cells coerced to int/float."""
    with open(path, newline="", encoding="utf-8") as f:
        first = f.readline()
        f.seek(0)
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        reader = csv.DictReader(f, delimiter=delim)
        fields = list(reader.fieldnames or [])
        rows = []
        for row in reader:
            rows.append({k: _coerce_cell(row.get(k, "")) for k in fields})
        return fields, rows
