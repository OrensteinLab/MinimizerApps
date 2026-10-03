"""DrawingML charts (clustered columns with mean+-std error bars, scatter
with error bars) for the hand-rolled .xlsx writer in xlsx_book.py.

Chart spec dicts:
  column: {title, y_title, cats, series:[{name, means, stds}], y_min?, y_max?}
  scatter: {type:'scatter', title, x_title, y_title, series:[{name, x, y, stds}]}
"""

from __future__ import print_function

from xlsx_util import xml_escape

_NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

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
        return "<c:numRef><c:f>%s</c:f>%s</c:numRef>" % (xml_escape(ref), cache)
    return _lit_pts(values or [], True)


def _lit_pts(values, numeric):
    pts = []
    for i, v in enumerate(values):
        if numeric:
            pts.append('<c:pt idx="%d"><c:v>%s</c:v></c:pt>' % (i, repr(float(v))))
        else:
            pts.append(
                '<c:pt idx="%d"><c:v>%s</c:v></c:pt>' % (i, xml_escape(v))
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
        '<c:overlay val="0"/></c:title>' % xml_escape(text)
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
    ) % (_NS_C, _NS_A, _NS_REL, xml_escape(title), plot_inner)


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
                xml_escape(ser.get("name") or "series"),
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
        xml_escape(style),
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
                xml_escape(ser.get("name") or "series"),
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



chart_xml = _chart_xml


CT_CHART = _CT_CHART
CT_DRAW = _CT_DRAW


def drawing_anchor_xml(start_row, chart_rids):
    """Charts in a 2-column grid starting at 1-based sheet row `start_row`."""
    return _drawing_xml(len(chart_rids), start_row - 1, chart_rids)
