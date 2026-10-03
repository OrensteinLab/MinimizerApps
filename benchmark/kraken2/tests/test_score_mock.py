#!/usr/bin/env python3
"""Unit tests for score_mock.py on a synthetic Kraken 2 report (no DB needed).

Run: python3 -m unittest discover -s benchmark/kraken2/tests
"""

import csv
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import score_mock  # noqa: E402

# Expected: E. coli (genus Escherichia) 60 %, Salmonella enterica 30 %,
# S. cerevisiae (yeast) 10 % by genomic DNA.
COMPOSITION = """taxid,name,domain,genus_taxid,genus_name,gdna_pct,genome_copy_pct,cell_pct,strain
562,Escherichia coli,bacteria,561,Escherichia,60,0,0,x
28901,Salmonella enterica,bacteria,590,Salmonella,30,0,0,x
4932,Saccharomyces cerevisiae,eukaryote,4930,Saccharomyces,10,0,0,x
"""

# 1000 reads: 100 unclassified, 900 classified.
#   E. coli 500 (incl. 100 on a strain), Escherichia-genus-only 20,
#   Salmonella enterica 200, Shigella flexneri 150 (off target), bacteria-only 30
#   (above genus), yeast 0, tiny off-target Listeria monocytogenes 2.
REPORT = [
    (10.0, 100, 100, "U", 0, "unclassified"),
    (90.0, 900, 0, "R", 1, "root"),
    (90.0, 900, 30, "D", 2, "  Bacteria"),
    (52.0, 520, 20, "G", 561, "    Escherichia"),
    (50.0, 500, 400, "S", 562, "      Escherichia coli"),
    (10.0, 100, 100, "S1", 83333, "        Escherichia coli K-12"),
    (15.0, 150, 0, "G", 620, "    Shigella"),
    (15.0, 150, 150, "S", 623, "      Shigella flexneri"),
    (20.0, 200, 0, "G", 590, "    Salmonella"),
    (20.0, 200, 200, "S", 28901, "      Salmonella enterica"),
    (0.2, 2, 0, "G", 1637, "    Listeria"),
    (0.2, 2, 2, "S", 1639, "      Listeria monocytogenes"),
]
# Root / Bacteria clades are patched to 902 in setUpClass (sum of the parts).


def write_report(path, rows):
    with open(path, "w") as f:
        for pct, clade, taxon, code, taxid, name in rows:
            f.write("{:.2f}\t{}\t{}\t{}\t{}\t{}\n".format(pct, clade, taxon, code, taxid, name))


class ScoreMockTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rows = [list(r) for r in REPORT]
        rows[1][1] = 902  # root clade
        rows[2][1] = 902  # bacteria clade
        rows[0][1] = 98   # unclassified so that totals are 1000
        rows[0][2] = 98
        cls.report = os.path.join(cls.tmp, "report.txt")
        write_report(cls.report, rows)
        cls.comp = os.path.join(cls.tmp, "composition.csv")
        with open(cls.comp, "w") as f:
            f.write(COMPOSITION)

    def metrics(self, **kw):
        rows = score_mock.score_report(self.report, self.comp, scheme="T", **kw)
        return {(r["scope"], r["rank"], r["min_frac"]): r for r in rows}

    def test_parse_report(self):
        unc, nodes = score_mock.parse_report(self.report)
        self.assertEqual(unc, 98)
        self.assertEqual(nodes[1].clade, 902)
        self.assertEqual(nodes[83333].parent.taxid, 562)
        self.assertEqual(nodes[562].parent.taxid, 561)

    def test_species_all(self):
        m = self.metrics()[("all", "species", 1e-3)]
        # detected at >=1e-3 of 1000 reads (>=1 read): E. coli, S. enterica,
        # Shigella flexneri, Listeria (2 reads); yeast missing.
        self.assertEqual((m["tp_taxa"], m["fp_taxa"], m["fn_taxa"]), (2, 2, 1))
        self.assertAlmostEqual(m["precision"], 0.5)
        self.assertAlmostEqual(m["recall"], 2 / 3.0)
        self.assertEqual(m["n_reads"], 1000)

    def test_threshold_removes_tiny_listeria(self):
        # 1e-2 * 1000 = 10 reads: Listeria (2 reads) drops out.
        m = self.metrics(min_fracs=(0.01,))[("all", "species", 0.01)]
        self.assertEqual((m["tp_taxa"], m["fp_taxa"], m["fn_taxa"]), (2, 1, 1))
        self.assertAlmostEqual(m["precision"], 2 / 3.0)

    def test_l1_and_bray_curtis(self):
        m = self.metrics()[("all", "species", 1e-4)]
        # observed over species reads: 500+150+200+2 = 852
        obs = {"ecoli": 500 / 852.0, "senterica": 200 / 852.0, "shig": 150 / 852.0,
               "list": 2 / 852.0}
        exp = {"ecoli": 0.6, "senterica": 0.3, "yeast": 0.1}
        l1 = (abs(exp["ecoli"] - obs["ecoli"]) + abs(exp["senterica"] - obs["senterica"])
              + exp["yeast"] + obs["shig"] + obs["list"])
        self.assertAlmostEqual(m["l1"], l1)
        self.assertAlmostEqual(m["bray_curtis"], l1 / 2)

    def test_read_fractions_sum_to_one(self):
        for (scope, rank, _t), m in self.metrics().items():
            total = (m["on_target_frac"] + m["off_target_frac"] + m["above_rank_frac"]
                     + m["unclassified_frac"])
            self.assertAlmostEqual(total, 1.0, msg=(scope, rank))
        m = self.metrics()[("all", "species", 1e-4)]
        self.assertAlmostEqual(m["classified_frac"], 0.902)
        self.assertAlmostEqual(m["on_target_frac"], 0.7)
        self.assertAlmostEqual(m["off_target_frac"], 0.152)
        self.assertIn("Shigella flexneri", m["top_off_target"])

    def test_genus_rank_uses_report_lineage(self):
        # no nodes.dmp: expected genus of E. coli comes from the report lineage
        # (561); yeast genus absent from the report -> composition-table genus.
        m = self.metrics()[("all", "genus", 1e-4)]
        self.assertEqual((m["tp_taxa"], m["fp_taxa"], m["fn_taxa"]), (2, 2, 1))
        # Escherichia clade (520) counts as on-target at genus rank
        self.assertAlmostEqual(m["on_target_frac"], (520 + 200) / 1000.0)

    def test_bacteria_scope_drops_yeast(self):
        m = self.metrics()[("bacteria", "species", 1e-4)]
        self.assertEqual(m["n_expected"], 2)
        self.assertEqual(m["fn_taxa"], 0)
        self.assertAlmostEqual(m["recall"], 1.0)

    def test_perfect_profile(self):
        tmp = tempfile.mkdtemp()
        rep = os.path.join(tmp, "r.txt")
        write_report(rep, [
            (0.0, 0, 0, "U", 0, "unclassified"),
            (100.0, 1000, 0, "R", 1, "root"),
            (100.0, 1000, 0, "D", 2, "  Bacteria"),
            (60.0, 600, 0, "G", 561, "    Escherichia"),
            (60.0, 600, 600, "S", 562, "      Escherichia coli"),
            (30.0, 300, 0, "G", 590, "    Salmonella"),
            (30.0, 300, 300, "S", 28901, "      Salmonella enterica"),
            (10.0, 100, 0, "G", 4930, "    Saccharomyces"),
            (10.0, 100, 100, "S", 4932, "      Saccharomyces cerevisiae"),
        ])
        rows = score_mock.score_report(rep, self.comp)
        for r in rows:
            self.assertAlmostEqual(r["l1"], 0.0, msg=r["rank"])
            self.assertAlmostEqual(r["f1"], 1.0)
            self.assertAlmostEqual(r["classified_frac"], 1.0)

    def test_bacteria_scope_removes_genus_level_yeast_reads(self):
        tmp = tempfile.mkdtemp()
        rep = os.path.join(tmp, "r.txt")
        write_report(rep, [
            (0.0, 0, 0, "U", 0, "unclassified"),
            (100.0, 1000, 0, "R", 1, "root"),
            (90.0, 900, 0, "D", 2, "  Bacteria"),
            (60.0, 600, 0, "G", 561, "    Escherichia"),
            (60.0, 600, 600, "S", 562, "      Escherichia coli"),
            (30.0, 300, 0, "G", 590, "    Salmonella"),
            (30.0, 300, 300, "S", 28901, "      Salmonella enterica"),
            (10.0, 100, 100, "G", 4930, "  Saccharomyces"),  # genus-level only
        ])
        rows = score_mock.score_report(rep, self.comp)
        m = {(r["scope"], r["rank"], r["min_frac"]): r for r in rows}
        b = m[("bacteria", "species", 1e-4)]
        self.assertEqual(b["n_reads"], 900)  # yeast genus reads removed
        self.assertAlmostEqual(b["classified_frac"], 1.0)
        self.assertAlmostEqual(m[("all", "species", 1e-4)]["above_rank_frac"], 0.1)

    def test_cli_writes_table(self):
        out = os.path.join(self.tmp, "m.tsv")
        self.assertEqual(score_mock.main(["--report", self.report, "--composition",
                                          self.comp, "--out-tsv", out, "--scheme", "S"]), 0)
        with open(out, newline="") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 2 * 2 * 3)
        self.assertEqual(set(rows[0]), set(score_mock.FIELDS))


if __name__ == "__main__":
    unittest.main()
