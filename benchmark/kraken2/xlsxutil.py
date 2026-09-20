#!/usr/bin/env python3
"""Minimal XLSX reader/writer (stdlib only).

Read Additional file 1, or write CSV + .xlsx workbooks without openpyxl.
"""

from __future__ import print_function

import csv
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

_MONOREPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_MONOREPO_ROOT / "benchmark" / "common"))

from xlsx_util import col_letter, root_rels_xml, write_csv, xml_escape  # noqa: E402

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

_INVALID_SHEET = re.compile(r'[:\\/?*\[\]]')
_NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_NS_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
_NS_CT = "http://schemas.openxmlformats.org/package/2006/content-types"


_xml_escape = xml_escape
_col_letter = col_letter


def sanitize_sheet_name(name, used):
    base = _INVALID_SHEET.sub(" ", str(name or "Sheet")).strip() or "Sheet"
    base = base[:31]
    cand = base
    i = 2
    while cand in used:
        suffix = "_%d" % i
        cand = (base[: 31 - len(suffix)] + suffix)
        i += 1
    used.add(cand)
    return cand


def _is_number(value):
    if isinstance(value, bool) or value is None:
        return False
    if isinstance(value, (int, float)):
        return True
    return False


def _cell_xml(col, row, value):
    ref = "%s%d" % (_col_letter(col), row)
    if isinstance(value, dict) and value.get("formula"):
        f = str(value["formula"])
        if f.startswith("="):
            f = f[1:]
        cached = value.get("value")
        if cached is not None and _is_number(cached):
            try:
                if cached == cached:
                    return '<c r="%s"><f>%s</f><v>%s</v></c>' % (
                        ref,
                        _xml_escape(f),
                        repr(float(cached)),
                    )
            except TypeError:
                pass
        return '<c r="%s"><f>%s</f></c>' % (ref, _xml_escape(f))
    if value is None or value == "":
        return '<c r="%s"/>' % ref
    if _is_number(value):
        try:
            if value != value:  # NaN
                return '<c r="%s" t="inlineStr"><is><t>nan</t></is></c>' % ref
        except TypeError:
            pass
        return '<c r="%s" t="n"><v>%s</v></c>' % (ref, repr(float(value) if isinstance(value, float) else value).rstrip("L"))
    return '<c r="%s" t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>' % (
        ref,
        _xml_escape(value),
    )


def _sheet_xml(fieldnames, rows, drawing_rid=None, extra_cells=None):
    ns_r = ' xmlns:r="%s"' % _NS_REL if drawing_rid else ""
    cells = {}
    for i, name in enumerate(fieldnames, 1):
        cells[(i, 1)] = name
    for ridx, row in enumerate(rows, 2):
        for i, name in enumerate(fieldnames, 1):
            val = row.get(name, "") if isinstance(row, dict) else row[i - 1]
            cells[(i, ridx)] = val
    if extra_cells:
        cells.update(extra_cells)
    rows_map = {}
    for (c, r), val in cells.items():
        rows_map.setdefault(r, []).append((c, val))
    lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="%s"%s><sheetData>' % (_NS_MAIN, ns_r),
    ]
    for ridx in sorted(rows_map):
        cols = sorted(rows_map[ridx])
        xml = "".join(_cell_xml(c, ridx, val) for c, val in cols)
        lines.append('<row r="%d">%s</row>' % (ridx, xml))
    lines.append("</sheetData>")
    if drawing_rid:
        lines.append('<drawing r:id="%s"/>' % drawing_rid)
    lines.append("</worksheet>")
    return "\n".join(lines)


_NS_C = "http://schemas.openxmlformats.org/drawingml/2006/chart"
_NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
_NS_XDR = "http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing"
_CT_CHART = "application/vnd.openxmlformats-officedocument.drawingml.chart+xml"
_CT_DRAW = "application/vnd.openxmlformats-officedocument.drawing+xml"
_CHART_COLORS = (
    "4C78A8",
    "F58518",
    "54A24B",
    "E45756",
    "72B7B2",
    "B279A2",
    "FF9DA6",
    "9D755D",
)


def _chart_num_src(ref, values):
    """Excel numRef (sheet formula) when ref is set; otherwise a numeric literal."""
    if ref:
        cache = ""
        if values:
            pts = "".join(
                '<c:pt idx="%d"><c:v>%s</c:v></c:pt>' % (i, repr(float(v)))
                for i, v in enumerate(values)
                if v is not None
            )
            n = sum(1 for v in values if v is not None)
            cache = (
                "<c:numCache><c:formatCode>General</c:formatCode>"
                '<c:ptCount val="%d"/>%s</c:numCache>' % (n, pts)
            )
        return "<c:numRef><c:f>%s</c:f>%s</c:numRef>" % (_xml_escape(ref), cache)
    return _lit_pts(values or [], True)


def _lit_pts(values, numeric):
    pts = []
    for i, v in enumerate(values):
        if numeric:
            pts.append('<c:pt idx="%d"><c:v>%s</c:v></c:pt>' % (i, repr(float(v))))
        else:
            pts.append(
                '<c:pt idx="%d"><c:v>%s</c:v></c:pt>' % (i, _xml_escape(v))
            )
    tag = "numLit" if numeric else "strLit"
    fmt = '<c:formatCode>General</c:formatCode>' if numeric else ""
    return '<c:%s>%s<c:ptCount val="%d"/>%s</c:%s>' % (
        tag,
        fmt,
        len(values),
        "".join(pts),
        tag,
    )


def _ax_scaling_xml(ymin=None, ymax=None):
    extra = ""
    if ymin is not None:
        extra += '<c:min val="%s"/>' % repr(float(ymin))
    if ymax is not None:
        extra += '<c:max val="%s"/>' % repr(float(ymax))
    return "<c:scaling><c:orientation val='minMax'/>%s</c:scaling>" % extra


def _ax_title_xml(text):
    if not text:
        return ""
    return (
        "<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/>"
        "<a:p><a:r><a:t>%s</a:t></a:r></a:p></c:rich></c:tx>"
        '<c:overlay val="0"/></c:title>' % _xml_escape(text)
    )


def _err_bars_xml(stds, direction="y"):
    if not any(float(x or 0) != 0 for x in (stds or [])):
        return ""
    return (
        "<c:errBars>"
        '<c:errDir val="%s"/>'
        '<c:errBarType val="both"/>'
        '<c:errValType val="cust"/>'
        '<c:noEndCap val="0"/>'
        "<c:plus>%s</c:plus><c:minus>%s</c:minus>"
        "</c:errBars>"
    ) % (direction, _lit_pts(stds, True), _lit_pts(stds, True))


def _chartspace_wrap(title, plot_inner):
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<c:chartSpace xmlns:c="%s" xmlns:a="%s" xmlns:r="%s">'
        "<c:chart>"
        "<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/>"
        "<a:p><a:r><a:rPr sz='1200' b='1'/><a:t>%s</a:t></a:r></a:p>"
        "</c:rich></c:tx><c:overlay val='0'/></c:title>"
        "<c:plotArea><c:layout/>%s</c:plotArea>"
        '<c:legend><c:legendPos val="b"/><c:overlay val="0"/></c:legend>'
        '<c:plotVisOnly val="1"/>'
        "</c:chart></c:chartSpace>"
    ) % (_NS_C, _NS_A, _NS_REL, _xml_escape(title), plot_inner)


def _scatter_chart_xml(spec):
    series = spec.get("series") or []
    title = spec.get("title") or ""
    x_title = spec.get("x_title") or ""
    y_title = spec.get("y_title") or ""
    x_log = bool(spec.get("x_log"))
    y_log = bool(spec.get("y_log"))
    style = spec.get("scatter_style") or "marker"
    msize = int(spec.get("marker_size") or 8)
    line_w = int(spec.get("line_width") or 12700)
    grid_xml = "<c:majorGridlines/>" if spec.get("grid", True) else ""
    ser_xml = []
    for i, ser in enumerate(series):
        color = _CHART_COLORS[i % len(_CHART_COLORS)]
        xs = ser.get("x") or []
        ys = ser.get("y") or []
        stds = ser.get("stds") or [0.0] * len(ys)
        x_src = _chart_num_src(ser.get("x_ref"), xs)
        y_src = _chart_num_src(ser.get("y_ref"), ys)
        has_line = style in ("line", "lineMarker", "smooth", "smoothMarker")
        ln = (
            "<c:spPr><a:ln w='%d'><a:solidFill><a:srgbClr val='%s'/>"
            "</a:solidFill></a:ln></c:spPr>" % (line_w, color)
            if has_line
            else "<c:spPr><a:ln><a:noFill/></a:ln></c:spPr>"
        )
        if style in ("line", "smooth") or int(spec.get("marker_size") or 8) <= 0:
            marker = "<c:marker><c:symbol val='none'/></c:marker>"
        else:
            marker = (
                "<c:marker><c:symbol val='circle'/><c:size val='%d'/>"
                "<c:spPr><a:solidFill><a:srgbClr val='%s'/></a:solidFill>"
                "<a:ln w='12700'><a:solidFill><a:srgbClr val='%s'/></a:solidFill></a:ln>"
                "</c:spPr></c:marker>" % (msize, color, color)
            )
        ser_xml.append(
            "<c:ser>"
            '<c:idx val="%d"/><c:order val="%d"/>'
            "<c:tx><c:v>%s</c:v></c:tx>"
            "%s%s"
            "<c:xVal>%s</c:xVal><c:yVal>%s</c:yVal>%s"
            "</c:ser>"
            % (
                i,
                i,
                _xml_escape(ser.get("name") or "series"),
                ln,
                marker,
                x_src,
                y_src,
                _err_bars_xml(stds, "y"),
            )
        )
    x_scale = (
        "<c:scaling><c:logBase val='10'/><c:orientation val='minMax'/></c:scaling>"
        if x_log
        else _ax_scaling_xml(spec.get("x_min"), spec.get("x_max"))
    )
    y_scale = (
        "<c:scaling><c:logBase val='10'/><c:orientation val='minMax'/></c:scaling>"
        if y_log
        else _ax_scaling_xml(spec.get("y_min"), spec.get("y_max"))
    )
    plot = (
        "<c:scatterChart>"
        '<c:scatterStyle val="%s"/>'
        '<c:varyColors val="0"/>'
        "%s"
        '<c:axId val="1"/><c:axId val="2"/>'
        "</c:scatterChart>"
        '<c:valAx><c:axId val="1"/>'
        "%s"
        '<c:delete val="0"/><c:axPos val="b"/>%s'
        "%s"
        '<c:crossAx val="2"/><c:crosses val="min"/>'
        "</c:valAx>"
        '<c:valAx><c:axId val="2"/>'
        "%s"
        '<c:delete val="0"/><c:axPos val="l"/>%s'
        "%s"
        '<c:crossAx val="1"/><c:crosses val="min"/>'
        "</c:valAx>"
    ) % (
        _xml_escape(style),
        "".join(ser_xml),
        x_scale,
        grid_xml,
        _ax_title_xml(x_title),
        y_scale,
        grid_xml,
        _ax_title_xml(y_title),
    )
    return _chartspace_wrap(title, plot)


def _bar_chart_xml(spec):
    cats = spec.get("cats") or []
    series = spec.get("series") or []
    title = spec.get("title") or ""
    y_title = spec.get("y_title") or ""
    x_title = spec.get("x_title") or ""
    ser_xml = []
    for i, ser in enumerate(series):
        color = _CHART_COLORS[i % len(_CHART_COLORS)]
        means = ser.get("means") or []
        stds = ser.get("stds") or [0.0] * len(means)
        ser_xml.append(
            "<c:ser>"
            '<c:idx val="%d"/><c:order val="%d"/>'
            "<c:tx><c:v>%s</c:v></c:tx>"
            '<c:spPr><a:solidFill><a:srgbClr val="%s"/></a:solidFill></c:spPr>'
            "<c:cat>%s</c:cat><c:val>%s</c:val>%s"
            "</c:ser>"
            % (
                i,
                i,
                _xml_escape(ser.get("name") or "series"),
                color,
                _lit_pts(cats, False),
                _lit_pts(means, True),
                _err_bars_xml(stds, "y"),
            )
        )
    plot = (
        "<c:barChart>"
        '<c:barDir val="col"/><c:grouping val="clustered"/>'
        '<c:varyColors val="0"/>'
        "%s"
        '<c:gapWidth val="80"/>'
        '<c:axId val="1"/><c:axId val="2"/>'
        "</c:barChart>"
        '<c:catAx><c:axId val="1"/>'
        "%s"
        '<c:delete val="0"/><c:axPos val="b"/>'
        "%s"
        '<c:tickLblPos val="nextTo"/><c:crossAx val="2"/>'
        "</c:catAx>"
        '<c:valAx><c:axId val="2"/>'
        "%s"
        '<c:delete val="0"/><c:axPos val="l"/><c:majorGridlines/>'
        "%s"
        '<c:crossAx val="1"/>'
        "</c:valAx>"
    ) % (
        "".join(ser_xml),
        _ax_scaling_xml(),
        _ax_title_xml(x_title),
        _ax_scaling_xml(spec.get("y_min"), spec.get("y_max")),
        _ax_title_xml(y_title),
    )
    return _chartspace_wrap(title, plot)


def _chart_xml(spec):
    if spec.get("type") == "scatter":
        return _scatter_chart_xml(spec)
    return _bar_chart_xml(spec)


def _drawing_xml(n_charts, start_row, chart_rids):
    """Place charts in a 2-column grid below the data table."""
    anchors = []
    width_cols, height_rows = 9, 16
    for i, rid in enumerate(chart_rids):
        col = (i % 2) * width_cols
        row = start_row + (i // 2) * height_rows
        anchors.append(
            "<xdr:twoCellAnchor>"
            "<xdr:from><xdr:col>%d</xdr:col><xdr:colOff>0</xdr:colOff>"
            "<xdr:row>%d</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:from>"
            "<xdr:to><xdr:col>%d</xdr:col><xdr:colOff>0</xdr:colOff>"
            "<xdr:row>%d</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:to>"
            "<xdr:graphicFrame>"
            "<xdr:nvGraphicFramePr>"
            '<xdr:cNvPr id="%d" name="Chart %d"/>'
            "<xdr:cNvGraphicFramePr><a:graphicFrameLocks noGrp='1'/>"
            "</xdr:cNvGraphicFramePr></xdr:nvGraphicFramePr>"
            "<xdr:xfrm><a:off x='0' y='0'/><a:ext cx='0' cy='0'/></xdr:xfrm>"
            "<a:graphic><a:graphicData uri='%s'>"
            '<c:chart xmlns:c="%s" xmlns:r="%s" r:id="%s"/>'
            "</a:graphicData></a:graphic>"
            "</xdr:graphicFrame><xdr:clientData/>"
            "</xdr:twoCellAnchor>"
            % (
                col,
                row,
                col + width_cols,
                row + height_rows - 1,
                i + 2,
                i + 1,
                _NS_C,
                _NS_C,
                _NS_REL,
                rid,
            )
        )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<xdr:wsDr xmlns:xdr="%s" xmlns:a="%s">%s</xdr:wsDr>'
        % (_NS_XDR, _NS_A, "".join(anchors))
    )


def write_xlsx(path, sheets, charts_by_name=None, extra_cells_by_name=None):
    """Write workbook. sheets: list of (name, fieldnames, rows).

    charts_by_name: {sheet_name: [chart_spec, ...]}
    Column charts: {title, y_title, cats, series:[{name,means,stds}]}
    Scatter: {type:'scatter', title, x_title, y_title, series:[{name,x,y,stds}]}
    extra_cells_by_name: {sheet_name: {(col, row): value}} 1-based cells.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    charts_by_name = charts_by_name or {}
    extra_cells_by_name = extra_cells_by_name or {}
    used = set()
    packed = []
    for name, fields, rows in sheets:
        packed.append((sanitize_sheet_name(name, used), fields, rows))
    ct_overrides = []
    wb_sheets = []
    rels = []
    extra_parts = []  # (zip_path, xml)
    chart_seq = 0
    drawing_seq = 0
    for i, (name, fields, rows) in enumerate(packed, 1):
        rid = "rId%d" % i
        ct_overrides.append(
            '<Override PartName="/xl/worksheets/sheet%d.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.'
            'spreadsheetml.worksheet+xml"/>' % i
        )
        wb_sheets.append(
            '<sheet name="%s" sheetId="%d" r:id="%s"/>' % (_xml_escape(name), i, rid)
        )
        rels.append(
            '<Relationship Id="%s" Type="%s/worksheet" Target="worksheets/sheet%d.xml"/>'
            % (rid, "http://schemas.openxmlformats.org/officeDocument/2006/relationships", i)
        )
        specs = charts_by_name.get(name) or []
        extra_cells = extra_cells_by_name.get(name) or {}
        drawing_rid = None
        extra_max_row = max((r for (_, r) in extra_cells), default=0)
        chart_start = max(len(rows) + 2, extra_max_row + 2)
        if specs:
            drawing_seq += 1
            d_idx = drawing_seq
            chart_rids = []
            drawing_rels = []
            for j, spec in enumerate(specs, 1):
                chart_seq += 1
                crid = "rId%d" % j
                chart_rids.append(crid)
                extra_parts.append(
                    ("xl/charts/chart%d.xml" % chart_seq, _chart_xml(spec))
                )
                ct_overrides.append(
                    '<Override PartName="/xl/charts/chart%d.xml" ContentType="%s"/>'
                    % (chart_seq, _CT_CHART)
                )
                drawing_rels.append(
                    '<Relationship Id="%s" Type="%s/chart" Target="../charts/chart%d.xml"/>'
                    % (
                        crid,
                        "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
                        chart_seq,
                    )
                )
            extra_parts.append(
                (
                    "xl/drawings/drawing%d.xml" % d_idx,
                    _drawing_xml(len(specs), chart_start, chart_rids),
                )
            )
            extra_parts.append(
                (
                    "xl/drawings/_rels/drawing%d.xml.rels" % d_idx,
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<Relationships xmlns="%s">%s</Relationships>'
                    % (_NS_PKG, "".join(drawing_rels)),
                )
            )
            extra_parts.append(
                (
                    "xl/worksheets/_rels/sheet%d.xml.rels" % i,
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<Relationships xmlns="%s">'
                    '<Relationship Id="rId1" Type="%s/drawing" '
                    'Target="../drawings/drawing%d.xml"/>'
                    "</Relationships>"
                    % (
                        _NS_PKG,
                        "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
                        d_idx,
                    ),
                )
            )
            ct_overrides.append(
                '<Override PartName="/xl/drawings/drawing%d.xml" ContentType="%s"/>'
                % (d_idx, _CT_DRAW)
            )
            drawing_rid = "rId1"
            extra_parts.append(
                (
                    "xl/worksheets/sheet%d.xml" % i,
                    _sheet_xml(
                        fields, rows, drawing_rid=drawing_rid, extra_cells=extra_cells
                    ),
                )
            )
        else:
            extra_parts.append(
                (
                    "xl/worksheets/sheet%d.xml" % i,
                    _sheet_xml(fields, rows, extra_cells=extra_cells),
                )
            )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="%s">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        "%s</Types>"
    ) % (_NS_CT, "".join(ct_overrides))
    root_rels = root_rels_xml()
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="%s" xmlns:r="%s"><sheets>%s</sheets>'
        '<calcPr calcId="0" fullCalcOnLoad="1"/></workbook>'
    ) % (_NS_MAIN, _NS_REL, "".join(wb_sheets))
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="%s">%s</Relationships>'
    ) % (_NS_PKG, "".join(rels))
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("xl/workbook.xml", workbook)
        zf.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        for zip_path, xml in extra_parts:
            zf.writestr(zip_path, xml)


def write_tables(
    xlsx_path,
    csv_dir,
    sheets,
    xlsx_sheets=None,
    charts_by_name=None,
    extra_cells_by_name=None,
):
    """Write each sheet as CSV under csv_dir and one combined xlsx_path.

    xlsx_sheets: optional list of (tab_name, fields, rows) for the workbook
    (short Excel names). CSVs always use `sheets` names.
    charts_by_name: optional {tab_name: [chart_spec, ...]} column or scatter charts.
    extra_cells_by_name: optional {tab_name: {(col, row): value}} (xlsx only).
    """
    if csv_dir:
        os.makedirs(csv_dir, exist_ok=True)
        for name, fields, rows in sheets:
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "sheet"
            write_csv(os.path.join(csv_dir, safe + ".csv"), fields, rows)
    if xlsx_path:
        write_xlsx(
            xlsx_path,
            xlsx_sheets if xlsx_sheets is not None else sheets,
            charts_by_name=charts_by_name,
            extra_cells_by_name=extra_cells_by_name,
        )


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


def pack_files(xlsx_path, named_paths, csv_dir=None):
    sheets = []
    for name, path in named_paths:
        if not path or not os.path.isfile(path):
            continue
        fields, rows = load_table(path)
        if not fields:
            continue
        sheets.append((name, fields, rows))
    if not sheets:
        sys.exit("no tables to pack")
    write_tables(xlsx_path, csv_dir, sheets)


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print(
            "Usage:\n"
            "  xlsxutil.py write --xlsx OUT.xlsx [--csv-dir DIR] NAME:PATH [NAME:PATH ...]\n"
            "  PATH may be .csv or .tsv"
        )
        return 0
    if argv[0] != "write":
        sys.exit("unknown command (expected write)")
    xlsx = csv_dir = None
    pairs = []
    i = 1
    while i < len(argv):
        if argv[i] == "--xlsx":
            xlsx = argv[i + 1]
            i += 2
        elif argv[i] == "--csv-dir":
            csv_dir = argv[i + 1]
            i += 2
        elif ":" in argv[i]:
            name, path = argv[i].split(":", 1)
            pairs.append((name, path))
            i += 1
        else:
            sys.exit("unexpected argument: " + argv[i])
    if not xlsx or not pairs:
        sys.exit("write requires --xlsx and NAME:PATH arguments")
    pack_files(xlsx, pairs, csv_dir=csv_dir)
    print("wrote", xlsx)
    return 0


if __name__ == "__main__":
    sys.exit(main())

