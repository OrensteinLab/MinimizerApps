import os
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import xlsx_book as xb  # noqa: E402
from xlsx_check import check_xlsx  # noqa: E402

CHART = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
         '<chartSpace xmlns="http://schemas.openxmlformats.org/drawingml/2006/chart">'
         '<chart><plotArea><layout/></plotArea></chart></chartSpace>')


def sheet(**kw):
    rows = [{"dataset": "d", "scheme": "xor", "t": 1.0, "t_ratio": 1.0, "t_z": 0.5},
            {"dataset": "d", "scheme": "greedymini", "t": 2.0, "t_ratio": 2.0},
            {"dataset": "e", "scheme": "xor", "t": 1.0, "t_ratio": 1.0, "t_z": -2.5}]
    s = {"name": "summary", "fields": ["dataset", "scheme", "t", "t_ratio", "t_z"], "rows": rows,
         "block_fields": ["dataset"], "merge_fields": ["dataset"]}
    s.update(kw)
    return s


class XlsxBook(unittest.TestCase):
    def write(self, sheets, **kw):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "b.xlsx")
        xb.write_book(p, sheets, **kw)
        self.assertEqual(check_xlsx(p), [])
        return zipfile.ZipFile(p)

    def test_no_fills_and_basic(self):
        z = self.write([sheet(categories=["a", "a", "b", "b", "b"], ratio_fields={"t_ratio": "lower"})])
        styles = z.read("xl/styles.xml").decode()
        self.assertNotIn("solid", styles)
        sh = z.read("xl/worksheets/sheet1.xml").decode()
        self.assertIn("colorScale", sh)
        self.assertNotIn("calcPr", z.read("xl/workbook.xml").decode())

    def test_polarity_and_baseline(self):
        z = self.write([sheet(ratio_fields={"t_ratio": "higher"}, ratio_baseline=1.5)])
        sh = z.read("xl/worksheets/sheet1.xml").decode()
        self.assertIn('<cfvo type="num" val="1.5"/>', sh)
        self.assertLess(sh.index("FFF8696B"), sh.index("FF63BE7B"))

    def test_banner_labels_fill_formula(self):
        rows = [{"dataset": "d", "scheme": "xor", "t": 1.0, "t_ratio": {"f": "=C3/C3"}, "__t_ratio_value": 1.0},
                {"dataset": "d", "scheme": "g", "t": 2.0, "t_ratio": {"f": "=C4/C3"}, "__t_ratio_value": 2.0}]
        s = sheet(rows=rows, banner=("x", {"f": "=SUM(C3:C4)"}), header_labels={"t": "time"},
                  ratio_fields={"t_ratio": "lower"},
                  cell_fill=lambda f, r, v: "FFCCFFFF" if f == "t" and v == 2.0 else None,
                  merge_key=lambda r, f: r.get(f))
        z = self.write([s])
        sh = z.read("xl/worksheets/sheet1.xml").decode()
        self.assertIn("SUM(C4:C5)", sh)
        self.assertIn("C4/C4", sh)
        self.assertIn('<row r="3">', sh)
        self.assertIn("time", sh)
        self.assertIn("calcPr", z.read("xl/workbook.xml").decode())
        self.assertIn("FFCCFFFF", z.read("xl/styles.xml").decode())
        self.assertEqual((xb.header_row(s), xb.data_start_row(s)), (3, 4))

    def test_charts(self):
        got = []

        def charts(name):
            got.append(name)
            return [{"xml": CHART, "col": 8, "row": 1}]

        z = self.write([sheet(charts=charts)])
        self.assertEqual(got, ["summary"])
        self.assertIn("xl/charts/chart1.xml", z.namelist())
        self.assertIn("xl/drawings/drawing1.xml", z.namelist())
        self.assertIn("<drawing", z.read("xl/worksheets/sheet1.xml").decode())

    def test_footer(self):
        s = sheet(footer_rows=[{"label": "corr", "cells": [("t_ratio", "CORREL(C3:C5,D3:D5)", 0.5)]}])
        z = self.write([s])
        self.assertIn("calcPr", z.read("xl/workbook.xml").decode())


if __name__ == "__main__":
    unittest.main()
