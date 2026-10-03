"""Styled multi-sheet .xlsx writer (stdlib only) in the layout the sshash and
kmc exports use: row 1 = category band (merged across contiguous columns that
share a category, category-coloured), row 2 = column headers, data from row
3, frozen below the header; vertically merged identity columns (dataset,
library, ...) inside blocks, a medium bottom border after each block, and
green/white/red colour scales on `*_ratio` columns.

A sheet is a dict:
  name        sheet name (sanitised to Excel's 31-char / no-[]:*?/\\ rules)
  fields      column ids (header text)
  rows        list of dicts (or lists) keyed by field
  categories  optional list, one category label per field (band + header fill)
  merge_fields  fields merged vertically inside a block (identical values)
  block_fields  fields whose change starts a new block (border after each
                block; merges never cross a block). Default: no blocks.
  ratio_fields  {field: "lower" | "higher"} colour scale; "lower" = below 1
                is better (green), "higher" = above 1 is better
  widths      {field: column width in characters} (default: from content)
  wrap_fields fields with wrapped text (long prose columns)
  freeze_cols number of leading columns to freeze (default 0)
  charts      list of chart specs (see xlsx_charts.py), placed below the table

Numbers are written as numbers, None / "" / NaN / inf as blank cells. Python
3.6 compatible.
"""

from __future__ import print_function

import math
import os
import re
import zipfile

from xlsx_charts import chart_xml, drawing_anchor_xml, CT_CHART, CT_DRAW
from xlsx_util import col_letter, root_rels_xml, xml_escape

HEADER_ROW = 2
DATA_START_ROW = 3

_NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_NS_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
_NS_CT = "http://schemas.openxmlformats.org/package/2006/content-types"
_INVALID_SHEET = re.compile(r"[:\\/?*\[\]]")

# Light fills for category bands / headers; assigned in first-seen order unless
# the category is listed in CATEGORY_FILLS.
CATEGORY_FILLS = {
    "id / params": "FFD9D9D9",
    "library": "FFD9D9D9",
    "dataset": "FFD9D9D9",
    "order": "FFD9D9D9",
    "build": "FFDDEBF7",
    "classify": "FFE2EFDA",
    "density": "FFFFF2CC",
    "nextminimizer": "FFFCE4D6",
    "vs reference": "FFEDEDED",
    "accuracy": "FFE4DFEC",
    "expected (theoretical)": "FFFFF2CC",
    "reference": "FFE7E6E6",
}
_FALLBACK_FILLS = (
    "FFDDEBF7", "FFE2EFDA", "FFFFF2CC", "FFFCE4D6", "FFE4DFEC",
    "FFD0E0E3", "FFF4CCCC", "FFEAD1DC", "FFD9EAD3", "FFFFE599",
)


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
    if _blank(value):
        return '<c r="%s"%s/>' % (ref, s_attr)
    if isinstance(value, bool):
        return '<c r="%s"%s t="b"><v>%d</v></c>' % (ref, s_attr, 1 if value else 0)
    if isinstance(value, (int, float)):
        return '<c r="%s"%s><v>%s</v></c>' % (ref, s_attr, repr(value))
    return '<c r="%s"%s t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>' % (
        ref, s_attr, xml_escape(value))


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
        if _blank(v):
            continue
        if isinstance(v, float):
            n = 9
        else:
            n = len(str(v))
        longest = max(longest, n)
    return max(8, min(longest + 2, 48))


def _ratio_bounds(rows, field):
    vals = [r.get(field) for r in rows
            if isinstance(r.get(field), (int, float)) and not _blank(r.get(field))]
    lo = min([1.0] + vals)
    hi = max([1.0] + vals)
    if lo == 1.0:
        lo = 0.99
    if hi == 1.0:
        hi = 1.01
    return lo, hi


def _color_scales(fields, rows, ratio_fields, last_row):
    parts = []
    prio = 0
    for field in fields:
        direction = ratio_fields.get(field)
        if not direction:
            continue
        prio += 1
        letter = col_letter(fields.index(field) + 1)
        lo, hi = _ratio_bounds(rows, field)
        good, bad = "FF63BE7B", "FFF8696B"
        c_lo, c_hi = (good, bad) if direction == "lower" else (bad, good)
        parts.append(
            '<conditionalFormatting sqref="%s%d:%s%d"><cfRule type="colorScale" '
            'priority="%d"><colorScale><cfvo type="num" val="%s"/><cfvo type="num" val="1"/>'
            '<cfvo type="num" val="%s"/><color rgb="%s"/><color rgb="FFFFFFFF"/>'
            '<color rgb="%s"/></colorScale></cfRule></conditionalFormatting>'
            % (letter, DATA_START_ROW, letter, last_row, prio, repr(lo), repr(hi), c_lo, c_hi))
    return "".join(parts)


def _sheet_xml(sheet, styles, cat_fill, drawing_rid):
    fields = list(sheet["fields"])
    rows = [r if isinstance(r, dict) else dict(zip(fields, r)) for r in sheet["rows"]]
    cats = sheet.get("categories") or [""] * len(fields)
    merge_fields = set(sheet.get("merge_fields") or ())
    block_fields = list(sheet.get("block_fields") or ())
    ratio_fields = sheet.get("ratio_fields") or {}
    wrap_fields = set(sheet.get("wrap_fields") or ())
    freeze_cols = int(sheet.get("freeze_cols") or 0)
    ncols = len(fields)

    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           '<worksheet xmlns="%s" xmlns:r="%s">' % (_NS_MAIN, _NS_REL)]
    pane = '<pane %sySplit="%d" topLeftCell="%s%d" activePane="%s" state="frozen"/>' % (
        ('xSplit="%d" ' % freeze_cols) if freeze_cols else "",
        HEADER_ROW, col_letter(freeze_cols + 1), DATA_START_ROW,
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

    # row 1: category band
    cells = []
    i = 0
    while i < ncols:
        cat = cats[i]
        start = i
        while i < ncols and cats[i] == cat:
            i += 1
        st = styles.get(bold=True, fill=cat_fill(cat), halign="center", valign="center")
        for c in range(start, i):
            cells.append(_cell(cat if c == start else None, c, 1, st))
        if i - 1 > start:
            merges.append("%s1:%s1" % (col_letter(start + 1), col_letter(i)))
    out.append('<row r="1">%s</row>' % "".join(cells))
    # row 2: headers
    cells = []
    for c, f in enumerate(fields):
        st = styles.get(bold=True, fill=cat_fill(cats[c]), halign="center",
                        valign="center", wrap=True)
        cells.append(_cell(f, c, HEADER_ROW, st))
    out.append('<row r="%d">%s</row>' % (HEADER_ROW, "".join(cells)))

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
    # vertical merge plan
    skip = set()
    anchor = set()
    if block_fields:
        for f in fields:
            if f not in merge_fields:
                continue
            c = fields.index(f)
            start = 0
            for k in range(1, n + 1):
                same = (k < n and block_ids[k] == block_ids[start]
                        and _mergekey(rows[k].get(f)) == _mergekey(rows[start].get(f))
                        and _mergekey(rows[start].get(f)) is not None)
                if same:
                    continue
                if k - 1 > start:
                    letter = col_letter(c + 1)
                    merges.append("%s%d:%s%d" % (letter, DATA_START_ROW + start,
                                                 letter, DATA_START_ROW + k - 1))
                    anchor.add((start, f))
                    for q in range(start + 1, k):
                        skip.add((q, f))
                start = k
    for k, r in enumerate(rows):
        rn = DATA_START_ROW + k
        border = block_fields and k in last_in_block
        cells = []
        for c, f in enumerate(fields):
            numfmt = "0.0000" if f.endswith("_ratio") else "General"
            if (k, f) in skip:
                cells.append(_cell(None, c, rn, styles.get(border=border, valign="center")))
                continue
            v = r.get(f)
            if (k, f) in anchor:
                st = styles.get(border=border, valign="center", halign="center")
            else:
                st = styles.get(border=border, numfmt=numfmt, wrap=f in wrap_fields,
                                valign="top" if f in wrap_fields else None)
            cells.append(_cell(v, c, rn, st))
        out.append('<row r="%d">%s</row>' % (rn, "".join(cells)))
    out.append("</sheetData>")
    if merges:
        out.append('<mergeCells count="%d">%s</mergeCells>' % (
            len(merges), "".join('<mergeCell ref="%s"/>' % m for m in merges)))
    if ratio_fields and rows:
        out.append(_color_scales(fields, rows, ratio_fields, DATA_START_ROW + n - 1))
    if drawing_rid:
        out.append('<drawing r:id="%s"/>' % drawing_rid)
    out.append("</worksheet>")
    return "\n".join(out)


def write_book(path, sheets):
    """Write the workbook; returns the list of (sanitised) sheet names."""
    dirname = os.path.dirname(str(path)) or "."
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    styles = Styles()
    assigned = {}

    def cat_fill(cat):
        if not cat:
            return None
        if cat in CATEGORY_FILLS:
            return CATEGORY_FILLS[cat]
        if cat not in assigned:
            assigned[cat] = _FALLBACK_FILLS[len(assigned) % len(_FALLBACK_FILLS)]
        return assigned[cat]

    used = set()
    names = []
    parts = []
    overrides = [
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]
    wb_sheets = []
    rels = ['<Relationship Id="rIdStyles" Type="%s/styles" Target="styles.xml"/>' % _NS_REL]
    chart_seq = 0
    drawing_seq = 0
    for i, sh in enumerate(sheets, 1):
        name = sanitize_sheet_name(sh["name"], used)
        names.append(name)
        wb_sheets.append('<sheet name="%s" sheetId="%d" r:id="rIdSheet%d"/>' % (xml_escape(name), i, i))
        rels.append('<Relationship Id="rIdSheet%d" Type="%s/worksheet" Target="worksheets/sheet%d.xml"/>'
                    % (i, _NS_REL, i))
        overrides.append('<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % i)
        charts = sh.get("charts") or []
        drawing_rid = None
        if charts:
            drawing_seq += 1
            drawing_rid = "rId1"
            drels = []
            rids = []
            for j, spec in enumerate(charts, 1):
                chart_seq += 1
                parts.append(("xl/charts/chart%d.xml" % chart_seq, chart_xml(spec)))
                overrides.append('<Override PartName="/xl/charts/chart%d.xml" ContentType="%s"/>' % (chart_seq, CT_CHART))
                drels.append('<Relationship Id="rId%d" Type="%s/chart" Target="../charts/chart%d.xml"/>' % (j, _NS_REL, chart_seq))
                rids.append("rId%d" % j)
            start_row = DATA_START_ROW + len(sh["rows"]) + 2
            parts.append(("xl/drawings/drawing%d.xml" % drawing_seq, drawing_anchor_xml(start_row, rids)))
            parts.append(("xl/drawings/_rels/drawing%d.xml.rels" % drawing_seq,
                          '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                          '<Relationships xmlns="%s">%s</Relationships>' % (_NS_PKG, "".join(drels))))
            parts.append(("xl/worksheets/_rels/sheet%d.xml.rels" % i,
                          '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                          '<Relationships xmlns="%s"><Relationship Id="rId1" Type="%s/drawing" '
                          'Target="../drawings/drawing%d.xml"/></Relationships>' % (_NS_PKG, _NS_REL, drawing_seq)))
            overrides.append('<Override PartName="/xl/drawings/drawing%d.xml" ContentType="%s"/>' % (drawing_seq, CT_DRAW))
        parts.append(("xl/worksheets/sheet%d.xml" % i, _sheet_xml(sh, styles, cat_fill, drawing_rid)))

    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="%s">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>%s</Types>'
    ) % (_NS_CT, "".join(overrides))
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="%s" xmlns:r="%s"><sheets>%s</sheets></workbook>'
    ) % (_NS_MAIN, _NS_REL, "".join(wb_sheets))
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
        for part, body in parts:
            z.writestr(part, body)
    return names
