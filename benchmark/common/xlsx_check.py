"""Structural checks for the hand-rolled .xlsx files (no Excel/LibreOffice/
openpyxl needed): well-formed XML parts, package relationships resolve, style
indices exist, rows/cells are ascending, merges do not overlap, sheet names
are valid. Returns a list of problem strings (empty = OK)."""

from __future__ import print_function

import posixpath
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REF = re.compile(r"^([A-Z]+)(\d+)$")


def _colnum(letters):
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n


def _rng(ref):
    a, _, b = ref.partition(":")
    ma, mb = _REF.match(a), _REF.match(b or a)
    return (_colnum(ma.group(1)), int(ma.group(2)), _colnum(mb.group(1)), int(mb.group(2)))


def check_xlsx(path):
    problems = []
    z = zipfile.ZipFile(path)
    names = set(z.namelist())
    trees = {}
    for n in sorted(names):
        if n.endswith(".xml") or n.endswith(".rels"):
            try:
                trees[n] = ET.fromstring(z.read(n))
            except ET.ParseError as e:
                problems.append("%s: not well-formed (%s)" % (n, e))
    for need in ("[Content_Types].xml", "_rels/.rels", "xl/workbook.xml", "xl/styles.xml"):
        if need not in names:
            problems.append("missing part " + need)
    if problems:
        return problems
    # relationships resolve
    for n, root in trees.items():
        if not n.endswith(".rels"):
            continue
        base = posixpath.dirname(posixpath.dirname(n))
        for rel in root:
            t = rel.attrib.get("Target", "")
            full = t.lstrip("/") if t.startswith("/") else posixpath.normpath(posixpath.join(base, t))
            if full not in names:
                problems.append("%s: target %s missing" % (n, full))
    # content-type overrides exist
    for o in trees["[Content_Types].xml"]:
        pn = o.attrib.get("PartName")
        if pn and pn.lstrip("/") not in names:
            problems.append("content-type override for missing part " + pn)
    nxf = len(trees["xl/styles.xml"].find(_M + "cellXfs"))
    sheets = trees["xl/workbook.xml"].find(_M + "sheets")
    seen = set()
    for sh in sheets:
        nm = sh.attrib["name"]
        if len(nm) > 31 or re.search(r"[:\\/?*\[\]]", nm) or nm.lower() in seen:
            problems.append("bad/duplicate sheet name %r" % nm)
        seen.add(nm.lower())
    for n in sorted(names):
        if not re.match(r"xl/worksheets/sheet\d+\.xml$", n):
            continue
        root = trees[n]
        prev_row = 0
        for row in root.iter(_M + "row"):
            r = int(row.attrib["r"])
            if r <= prev_row:
                problems.append("%s: row %d out of order" % (n, r))
            prev_row = r
            prev_col = 0
            for c in row:
                m = _REF.match(c.attrib["r"])
                if int(m.group(2)) != r:
                    problems.append("%s: cell %s in row %d" % (n, c.attrib["r"], r))
                col = _colnum(m.group(1))
                if col <= prev_col:
                    problems.append("%s: cell %s out of order" % (n, c.attrib["r"]))
                prev_col = col
                s = int(c.attrib.get("s", 0))
                if s >= nxf:
                    problems.append("%s: style %d >= %d" % (n, s, nxf))
        mc = root.find(_M + "mergeCells")
        if mc is not None:
            rects = [_rng(m.attrib["ref"]) for m in mc]
            if int(mc.attrib["count"]) != len(rects):
                problems.append("%s: mergeCells count mismatch" % n)
            for i in range(len(rects)):
                for j in range(i + 1, len(rects)):
                    a, b = rects[i], rects[j]
                    if not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1]):
                        problems.append("%s: overlapping merges %s %s" % (n, rects[i], rects[j]))
        # child element order of worksheet
        order = ["sheetPr", "dimension", "sheetViews", "sheetFormatPr", "cols", "sheetData",
                 "mergeCells", "conditionalFormatting", "drawing"]
        idx = [order.index(ch.tag.replace(_M, "")) for ch in root if ch.tag.replace(_M, "") in order]
        if idx != sorted(idx):
            problems.append("%s: worksheet children out of schema order" % n)
    return problems


if __name__ == "__main__":
    bad = 0
    for p in sys.argv[1:]:
        pr = check_xlsx(p)
        print(p, "OK" if not pr else "%d problems" % len(pr))
        for x in pr[:20]:
            print("  ", x)
        bad += bool(pr)
    sys.exit(1 if bad else 0)
