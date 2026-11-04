"""Styled multi-sheet .xlsx writer (stdlib only), the one writer all four apps
(kraken2, kmc, sshash, minimap2) use: row 1 = category band (merged across
contiguous columns that share a category; bold, centered, no fill), row 2 =
column headers, data from row 3, frozen below the header; vertically merged identity columns (dataset,
library, ...) inside blocks, a medium bottom border after each block, and
green/white/red colour scales on `*_ratio` columns.

A sheet is a dict:
  name        sheet name (sanitised to Excel's 31-char / no-[]:*?/\\ rules)
  fields      column ids (header text)
  rows        list of dicts (or lists) keyed by field
  categories  optional list, one category label per field (band text)
  merge_fields  fields merged vertically inside a block (identical values)
  merge_scope   optional {field: [fields]}: merge that field across runs of rows
                that agree on those fields instead of inside a block (so a
                library / dataset cell can span several bordered blocks)
  block_fields  fields whose change starts a new block (border after each
                block; merges never cross a block). Default: no blocks.
  ratio_fields  {field: "lower" | "higher"} colour scale; "lower" = below the
                baseline is better (green), "higher" = above is better. Any
                field may be listed (not only `*_ratio`), so a caller restricts
                the scale to a subset or sets per-metric polarity just by what
                it puts in the dict.
  ratio_baseline  white point of the scales (default 1.0); a number for all
                ratio fields or {field: number}
  banner      optional (label, value) written in row 1 (A1, B1) above the
                category band; value may be a number/text or {"f": formula}.
                Everything below shifts down one row (see row_offset()).
  header_labels  optional {field: displayed header text}; `fields` stay the
                  row keys (merge/block/ratio keys always use the field id)
  cell_fill   optional callable (field, row_dict, value) -> "FFRRGGBB" | None,
                background fill of a data cell (merged anchors included)
  merge_key   optional callable (row_dict, field) -> comparable, replaces the
                raw value when deciding whether adjacent rows merge (needed
                when the cell holds a per-row formula)
  charts      optional list of chart dicts, or a callable (sheet_name) -> list
                (sheet_name = the sanitised name, for cell references). A chart:
                {"xml": complete chartSpace xml, "col": 0-based anchor column,
                 "row": 0-based anchor row, "cx": width EMU, "cy": height EMU}.
                The writer adds the drawing/chart parts and relationships.

A cell value may be a formula: {"f": "=A1/B1" (leading "=" optional), "v":
cached value or None}. Formulas, footer formulas and the banner are authored
for the layout WITHOUT a banner (data from DATA_START_ROW); the writer shifts
their row numbers by the banner offset. To get a numeric colour-scale range
for a formula column, put the plain value in row["__<field>_value"] (not
written). Use data_start_row(sheet) when you need the real first data row.
  widths      {field: column width in characters} (default: from content)
  wrap_fields fields with wrapped text (long prose columns)
  freeze_cols number of leading columns to freeze (default 0)
  footer_rows rows written right under the data: [{"label": text (column A),
                "cells": [(field, formula, cached value or None), ...]}]; the
                formula (no leading "=") is a live Excel formula and the cached
                value is stored next to it for viewers that do not recalculate

Numbers are written as numbers, None / "" / NaN / inf as blank cells. Python
3.6 compatible. No header/band fills: bold, centered, wrapped only.
"""

from __future__ import print_function

import math
import os
import re
import zipfile

from xlsx_util import col_letter, root_rels_xml, xml_escape

HEADER_ROW = 2
DATA_START_ROW = 3

_NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_NS_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
_NS_CT = "http://schemas.openxmlformats.org/package/2006/content-types"
_INVALID_SHEET = re.compile(r"[:\\/?*\[\]]")

def sanitize_sheet_name(name, used):
    base = _INVALID_SHEET.sub(" ", str(name or "Sheet")).strip() or "Sheet"
    base = base[:31]
    cand = base
    i = 2
    while cand in used:
        suffix = "_%d" % i
        cand = base[: 31 - len(suffix)] + suffix
        i += 1
    used.add(cand)
    return cand


# ---------------------------------------------------------------------------
# styles: a registry that hands out cellXfs indices on demand
# ---------------------------------------------------------------------------

_NUMFMTS = {"General": 0, "0.0000": 164, "0.000": 165, "#,##0": 3}


class Styles(object):
    def __init__(self):
        self.fills = [None, None]  # 0 none, 1 gray125 (reserved)
        self.xfs = []  # (bold, fill_idx, border, numfmt_id, wrap, halign, valign)
        self._xf_index = {}
        self.get()  # xf 0 = default

    def fill(self, rgb):
        if rgb is None:
            return 0
        if rgb not in self.fills:
            self.fills.append(rgb)
        return self.fills.index(rgb)

    def get(self, bold=False, fill=None, border=False, numfmt="General",
            wrap=False, halign=None, valign=None):
        key = (bool(bold), self.fill(fill), bool(border), _NUMFMTS[numfmt],
               bool(wrap), halign, valign)
        if key not in self._xf_index:
            self._xf_index[key] = len(self.xfs)
            self.xfs.append(key)
        return self._xf_index[key]

    def xml(self):
        fills = ['<fill><patternFill patternType="none"/></fill>',
                 '<fill><patternFill patternType="gray125"/></fill>']
        for rgb in self.fills[2:]:
            fills.append(
                '<fill><patternFill patternType="solid"><fgColor rgb="%s"/>'
                '<bgColor indexed="64"/></patternFill></fill>' % rgb)
        xfs = []
        for bold, fill, border, numfmt, wrap, halign, valign in self.xfs:
            attrs = ['numFmtId="%d"' % numfmt, 'fontId="%d"' % (1 if bold else 0),
                     'fillId="%d"' % fill, 'borderId="%d"' % (1 if border else 0),
                     'xfId="0"']
            if numfmt:
                attrs.append('applyNumberFormat="1"')
            if bold:
                attrs.append('applyFont="1"')
            if fill:
                attrs.append('applyFill="1"')
            if border:
                attrs.append('applyBorder="1"')
            if wrap or halign or valign:
                al = []
                if halign:
                    al.append('horizontal="%s"' % halign)
                if valign:
                    al.append('vertical="%s"' % valign)
                if wrap:
                    al.append('wrapText="1"')
                attrs.append('applyAlignment="1"')
                xfs.append('<xf %s><alignment %s/></xf>' % (" ".join(attrs), " ".join(al)))
            else:
                xfs.append('<xf %s/>' % " ".join(attrs))
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<styleSheet xmlns="%s">'
            '<numFmts count="2"><numFmt numFmtId="164" formatCode="0.0000"/>'
            '<numFmt numFmtId="165" formatCode="0.000"/></numFmts>'
            '<fonts count="2">'
            '<font><sz val="11"/><color theme="1"/><name val="Calibri"/></font>'
            '<font><b/><sz val="11"/><color theme="1"/><name val="Calibri"/></font>'
            '</fonts>'
            '<fills count="%d">%s</fills>'
            '<borders count="2">'
            '<border><left/><right/><top/><bottom/><diagonal/></border>'
            '<border><left/><right/><top/><bottom style="medium"><color auto="1"/></bottom>'
            '<diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="%d">%s</cellXfs>'
            '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
            '</styleSheet>'
        ) % (_NS_MAIN, len(fills), "".join(fills), len(xfs), "".join(xfs))


# ---------------------------------------------------------------------------
# cells and sheets
# ---------------------------------------------------------------------------

def _blank(value):
    if value is None or value == "":
        return True
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return True
    return False


def _cell(value, col, row, style=0):
    """col is 0-based."""
    ref = "%s%d" % (col_letter(col + 1), row)
    s_attr = ' s="%d"' % style if style else ""
    if isinstance(value, dict) and "f" in value:
        return _formula_cell(str(value["f"]).lstrip("="), value.get("v"), col, row, style)
    if _blank(value):
        return '<c r="%s"%s/>' % (ref, s_attr)
    if isinstance(value, bool):
        return '<c r="%s"%s t="b"><v>%d</v></c>' % (ref, s_attr, 1 if value else 0)
    if isinstance(value, (int, float)):
        return '<c r="%s"%s><v>%s</v></c>' % (ref, s_attr, repr(value))
    return '<c r="%s"%s t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>' % (
        ref, s_attr, xml_escape(value))


def _formula_cell(formula, cached, col, row, style=0):
    ref = "%s%d" % (col_letter(col + 1), row)
    s_attr = ' s="%d"' % style if style else ""
    if _blank(cached):
        return '<c r="%s"%s><f>%s</f></c>' % (ref, s_attr, xml_escape(formula))
    if isinstance(cached, (int, float)) and not isinstance(cached, bool):
        return '<c r="%s"%s><f>%s</f><v>%s</v></c>' % (ref, s_attr, xml_escape(formula), repr(cached))
    return '<c r="%s"%s t="str"><f>%s</f><v>%s</v></c>' % (
        ref, s_attr, xml_escape(formula), xml_escape(cached))


def _row_get(row, field, idx):
    if isinstance(row, dict):
        return row.get(field)
    return row[idx] if idx < len(row) else None


def _mergekey(v):
    if _blank(v):
        return None
    return v


def _col_width(field, rows, sheet):
    w = sheet.get("widths", {}).get(field)
    if w:
        return w
    longest = len(str(field))
    for r in rows[:400]:
        v = r.get(field) if isinstance(r, dict) else None
        if _blank(v) or isinstance(v, dict):
            continue
        if isinstance(v, float):
            n = 9
        else:
            n = len(str(v))
        longest = max(longest, n)
    return max(8, min(longest + 2, 48))


_CELL_REF_RE = re.compile(r"(\$?)([A-Z]{1,3})(\$?)(\d+)")


def shift_formula_rows(formula, delta):
    """Shift every cell-reference row number in `formula` by `delta` (what the
    writer applies to formulas authored for the banner-less layout)."""
    if not delta:
        return formula
    return _CELL_REF_RE.sub(
        lambda m: "%s%s%s%d" % (m.group(1), m.group(2), m.group(3), int(m.group(4)) + delta),
        formula)


def _shift_value(v, delta):
    if delta and isinstance(v, dict) and "f" in v:
        out = dict(v)
        out["f"] = shift_formula_rows(str(v["f"]), delta)
        return out
    return v


def row_offset(sheet):
    """1 if the sheet has a banner row above the band, else 0."""
    return 1 if sheet.get("banner") else 0


def header_row(sheet):
    """1-based row of the column headers for this sheet spec."""
    return HEADER_ROW + row_offset(sheet)


def data_start_row(sheet):
    """1-based first data row for this sheet spec (DATA_START_ROW, +1 with a banner)."""
    return DATA_START_ROW + row_offset(sheet)


def _num(v):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or _blank(v):
        return None
    return v


def _scale_bounds(rows, field, baseline):
    vals = []
    for r in rows:
        v = _num(r.get("__%s_value" % field))
        if v is None:
            v = _num(r.get(field))
        if v is not None:
            vals.append(v)
    lo = min([baseline] + vals)
    hi = max([baseline] + vals)
    # keep the baseline a genuine interior point (a data-driven min/max
    # degenerates when every value sits on one side of it)
    if lo == baseline:
        lo = baseline - 0.01
    if hi == baseline:
        hi = baseline + 0.01
    return lo, hi


def _conditional_formats(sheet, fields, rows, first_row, last_row):
    parts = []
    prio = 0
    ratio_fields = sheet.get("ratio_fields") or {}
    base_spec = sheet.get("ratio_baseline", 1.0)
    for field in fields:
        direction = ratio_fields.get(field)
        if not direction:
            continue
        base = base_spec.get(field, 1.0) if isinstance(base_spec, dict) else base_spec
        prio += 1
        letter = col_letter(fields.index(field) + 1)
        lo, hi = _scale_bounds(rows, field, base)
        good, bad = "FF63BE7B", "FFF8696B"
        c_lo, c_hi = (good, bad) if direction == "lower" else (bad, good)
        parts.append(
            '<conditionalFormatting sqref="%s%d:%s%d"><cfRule type="colorScale" '
            'priority="%d"><colorScale><cfvo type="num" val="%s"/><cfvo type="num" val="%s"/>'
            '<cfvo type="num" val="%s"/><color rgb="%s"/><color rgb="FFFFFFFF"/>'
            '<color rgb="%s"/></colorScale></cfRule></conditionalFormatting>'
            % (letter, first_row, letter, last_row, prio, repr(lo), repr(base), repr(hi),
               c_lo, c_hi))
    return "".join(parts)


def _has_formula(sheet):
    if sheet.get("banner") and isinstance(sheet["banner"][1], dict):
        return True
    if sheet.get("footer_rows"):
        return True
    for r in sheet["rows"]:
        vals = r.values() if isinstance(r, dict) else r
        for v in vals:
            if isinstance(v, dict) and "f" in v:
                return True
    return False


def _sheet_xml(sheet, styles, sheet_name="", drawing_rid=None):
    fields = list(sheet["fields"])
    off = row_offset(sheet)
    hdr = HEADER_ROW + off
    d0 = DATA_START_ROW + off
    rows = [r if isinstance(r, dict) else dict(zip(fields, r)) for r in sheet["rows"]]
    if off:
        rows = [dict((k, _shift_value(v, off)) for k, v in r.items()) for r in rows]
    cats = sheet.get("categories") or [""] * len(fields)
    labels = sheet.get("header_labels") or {}
    merge_fields = set(sheet.get("merge_fields") or ())
    block_fields = list(sheet.get("block_fields") or ())
    wrap_fields = set(sheet.get("wrap_fields") or ())
    freeze_cols = int(sheet.get("freeze_cols") or 0)
    cell_fill = sheet.get("cell_fill")
    merge_key = sheet.get("merge_key")
    ncols = len(fields)

    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="%s" xmlns:r="%s">' % (_NS_MAIN, _NS_REL)]
    pane = '<pane %sySplit="%d" topLeftCell="%s%d" activePane="%s" state="frozen"/>' % (
        ('xSplit="%d" ' % freeze_cols) if freeze_cols else "",
        hdr, col_letter(freeze_cols + 1), d0,
        "bottomRight" if freeze_cols else "bottomLeft")
    out.append('<sheetViews><sheetView workbookViewId="0">%s</sheetView></sheetViews>' % pane)
    out.append('<sheetFormatPr defaultRowHeight="15"/>')
    cols = []
    for i, f in enumerate(fields, 1):
        cols.append('<col min="%d" max="%d" width="%s" customWidth="1"/>'
                    % (i, i, _col_width(f, rows, sheet)))
    if cols:
        out.append("<cols>%s</cols>" % "".join(cols))
    out.append("<sheetData>")
    merges = []

    if off:
        label, value = sheet["banner"]
        out.append('<row r="1">%s%s</row>' % (
            _cell(label, 0, 1, styles.get(bold=True)),
            _cell(_shift_value(value, off), 1, 1)))
    # category band
    band = off + 1
    cells = []
    i = 0
    while i < ncols:
        cat = cats[i]
        start = i
        while i < ncols and cats[i] == cat:
            i += 1
        st = styles.get(bold=True, halign="center", valign="center")
        for c in range(start, i):
            cells.append(_cell(cat if c == start else None, c, band, st))
        if i - 1 > start:
            merges.append("%s%d:%s%d" % (col_letter(start + 1), band, col_letter(i), band))
    out.append('<row r="%d">%s</row>' % (band, "".join(cells)))
    # headers
    cells = []
    for c, f in enumerate(fields):
        st = styles.get(bold=True, halign="center", valign="center", wrap=True)
        cells.append(_cell(labels.get(f, f), c, hdr, st))
    out.append('<row r="%d">%s</row>' % (hdr, "".join(cells)))

    # blocks
    def block_of(r):
        return tuple(_mergekey(r.get(f)) for f in block_fields)

    n = len(rows)
    block_ids = []
    cur = -1
    prev = None
    for k, r in enumerate(rows):
        b = block_of(r) if block_fields else ()
        if k == 0 or b != prev:
            cur += 1
            prev = b
        block_ids.append(cur)
    last_in_block = set(k for k in range(n)
                        if k == n - 1 or block_ids[k] != block_ids[k + 1])
    scope_cache = {}

    def ids_for(f):
        scope = (sheet.get("merge_scope") or {}).get(f)
        if not scope:
            return block_ids
        key = tuple(scope)
        if key not in scope_cache:
            ids, c0, p0 = [], -1, None
            for k2, r2 in enumerate(rows):
                b2 = tuple(_mergekey(r2.get(x)) for x in key)
                if k2 == 0 or b2 != p0:
                    c0 += 1
                    p0 = b2
                ids.append(c0)
            scope_cache[key] = ids
        return scope_cache[key]

    def mval(k, f):
        if merge_key is not None:
            v = merge_key(rows[k], f)
        else:
            v = rows[k].get(f)
            if isinstance(v, dict) and "f" in v:
                v = str(v["f"])
        return _mergekey(v)

    # vertical merge plan
    skip = set()
    anchor = set()
    if block_fields:
        for f in fields:
            if f not in merge_fields:
                continue
            c = fields.index(f)
            ids = ids_for(f)
            start = 0
            for k in range(1, n + 1):
                same = (k < n and ids[k] == ids[start]
                        and mval(k, f) == mval(start, f)
                        and mval(start, f) is not None)
                if same:
                    continue
                if k - 1 > start:
                    letter = col_letter(c + 1)
                    merges.append("%s%d:%s%d" % (letter, d0 + start, letter, d0 + k - 1))
                    anchor.add((start, f))
                    for q in range(start + 1, k):
                        skip.add((q, f))
                start = k
    for k, r in enumerate(rows):
        rn = d0 + k
        border = bool(block_fields) and k in last_in_block
        cells = []
        for c, f in enumerate(fields):
            numfmt = "0.0000" if f.endswith("_ratio") else "General"
            if (k, f) in skip:
                cells.append(_cell(None, c, rn, styles.get(border=border, valign="center")))
                continue
            v = r.get(f)
            fill = cell_fill(f, r, v) if cell_fill else None
            if (k, f) in anchor:
                st = styles.get(border=border, valign="center", halign="center", fill=fill)
            else:
                st = styles.get(border=border, numfmt=numfmt, wrap=f in wrap_fields,
                                valign="top" if f in wrap_fields else None, fill=fill)
            cells.append(_cell(v, c, rn, st))
        out.append('<row r="%d">%s</row>' % (rn, "".join(cells)))
    for i, fr in enumerate(sheet.get("footer_rows") or ()):
        rn = d0 + n + i
        cells = [_cell(fr.get("label"), 0, rn, styles.get(bold=True))]
        by_col = sorted((fields.index(f), fm, cv) for f, fm, cv in fr.get("cells", ()) if f in fields)
        for c, fm, cv in by_col:
            if c == 0:
                continue
            cells.append(_formula_cell(shift_formula_rows(fm, off), cv, c, rn,
                                       styles.get(numfmt="0.000")))
        out.append('<row r="%d">%s</row>' % (rn, "".join(cells)))
    out.append("</sheetData>")
    if merges:
        out.append('<mergeCells count="%d">%s</mergeCells>' % (
            len(merges), "".join('<mergeCell ref="%s"/>' % m for m in merges)))
    if rows and sheet.get("ratio_fields"):
        out.append(_conditional_formats(sheet, fields, rows, d0, d0 + n - 1))
    if drawing_rid:
        out.append('<drawing r:id="%s"/>' % drawing_rid)
    out.append("</worksheet>")
    return "\n".join(out)


_DRAWING_NS = "http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing"


def _rels_xml(items):
    body = "".join('<Relationship Id="%s" Type="%s" Target="%s"/>' % i for i in items)
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="%s">%s</Relationships>' % (_NS_PKG, body))


def _drawing_xml(charts):
    anchors = []
    for n, ch in enumerate(charts, 1):
        anchors.append(
            '<oneCellAnchor><from><col>%d</col><colOff>0</colOff><row>%d</row><rowOff>0</rowOff></from>'
            '<ext cx="%d" cy="%d"/><graphicFrame macro=""><nvGraphicFramePr>'
            '<cNvPr id="%d" name="Chart %d"/><cNvGraphicFramePr/></nvGraphicFramePr>'
            '<xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/></xfrm>'
            '<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/chart">'
            '<c:chart xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart" r:id="rId%d"/>'
            '</a:graphicData></a:graphic></graphicFrame><clientData/></oneCellAnchor>'
            % (ch["col"], ch["row"], ch.get("cx", 5400000), ch.get("cy", 2880000), n + 1, n, n))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<wsDr xmlns="%s" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
            'xmlns:r="%s">%s</wsDr>' % (_DRAWING_NS, _NS_REL, "".join(anchors)))


def write_book(path, sheets, recalc_on_load=None):
    """Write the workbook; returns the list of (sanitised) sheet names.

    sheets          list of sheet dicts (see module docstring).
    recalc_on_load  True/False to force `fullCalcOnLoad` (+ A1 refMode); None =
                    on when any sheet has formulas (footer rows, {"f": ...}
                    cells, a formula banner).
    """
    dirname = os.path.dirname(str(path)) or "."
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    styles = Styles()
    used = set()
    names = []
    parts = []
    extra = []
    overrides = [
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]
    wb_sheets = []
    rels = ['<Relationship Id="rIdStyles" Type="%s/styles" Target="styles.xml"/>' % _NS_REL]
    chart_no = 0
    drawing_no = 0
    any_formula = False
    for i, sh in enumerate(sheets, 1):
        name = sanitize_sheet_name(sh["name"], used)
        names.append(name)
        wb_sheets.append('<sheet name="%s" sheetId="%d" r:id="rIdSheet%d"/>' % (xml_escape(name), i, i))
        rels.append('<Relationship Id="rIdSheet%d" Type="%s/worksheet" Target="worksheets/sheet%d.xml"/>'
                    % (i, _NS_REL, i))
        overrides.append('<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % i)
        charts = sh.get("charts")
        if callable(charts):
            charts = charts(name)
        charts = list(charts or ())
        rid = None
        if charts:
            drawing_no += 1
            rid = "rId1"
            drels = []
            for n, ch in enumerate(charts, 1):
                chart_no += 1
                extra.append(("xl/charts/chart%d.xml" % chart_no, ch["xml"]))
                drels.append(("rId%d" % n, _NS_REL + "/chart", "../charts/chart%d.xml" % chart_no))
                overrides.append('<Override PartName="/xl/charts/chart%d.xml" ContentType="application/vnd.openxmlformats-officedocument.drawingml.chart+xml"/>' % chart_no)
            extra.append(("xl/drawings/drawing%d.xml" % drawing_no, _drawing_xml(charts)))
            extra.append(("xl/drawings/_rels/drawing%d.xml.rels" % drawing_no, _rels_xml(drels)))
            extra.append(("xl/worksheets/_rels/sheet%d.xml.rels" % i, _rels_xml(
                [("rId1", _NS_REL + "/drawing", "../drawings/drawing%d.xml" % drawing_no)])))
            overrides.append('<Override PartName="/xl/drawings/drawing%d.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/>' % drawing_no)
        any_formula = any_formula or _has_formula(sh)
        parts.append(("xl/worksheets/sheet%d.xml" % i, _sheet_xml(sh, styles, name, rid)))

    recalc = any_formula if recalc_on_load is None else recalc_on_load
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="%s">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>%s</Types>'
    ) % (_NS_CT, "".join(overrides))
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="%s" xmlns:r="%s"><sheets>%s</sheets>%s</workbook>'
    ) % (_NS_MAIN, _NS_REL, "".join(wb_sheets),
         '<calcPr refMode="A1" fullCalcOnLoad="1"/>' if recalc else "")
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="%s">%s</Relationships>' % (_NS_PKG, "".join(rels))
    )
    with zipfile.ZipFile(str(path), "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", root_rels_xml())
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/styles.xml", styles.xml())
        for part, body in parts + extra:
            z.writestr(part, body)
    return names
