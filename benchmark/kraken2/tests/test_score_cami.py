#!/usr/bin/env python3
"""Unit tests for score_cami.py / prepare_cami_reads.py / the comparison
summaries on tiny synthetic data (no DB needed).

Run: python3 -m unittest discover -s benchmark/kraken2/tests
"""

import gzip
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import compare_benchmark  # noqa: E402
import prepare_cami_reads  # noqa: E402
import score_cami  # noqa: E402
from evaluate_calls import load_taxonomy  # noqa: E402

# taxid, parent, rank
NODES = [
    (1, 1, "no rank"),
    (2, 1, "superkingdom"),        # Bacteria
    (2157, 1, "superkingdom"),     # Archaea
    (543, 2, "family"),
    (561, 543, "genus"),           # Escherichia
    (562, 561, "species"),         # E. coli
    (83333, 562, "no rank"),       # strain
    (570, 543, "genus"),           # Klebsiella
    (573, 570, "species"),         # K. pneumoniae
    (580, 543, "genus"),           # Citrobacter: NOT in the library
    (581, 580, "species"),
    (2200, 2157, "genus"),         # an archaeal genus in the library
    (2201, 2200, "species"),
    (777, 543, "no rank"),         # CAMI new_genus label (family-level child)
]


def write_nodes(path):
    with open(path, "w") as f:
        for t, p, r in NODES:
            f.write("{}\t|\t{}\t|\t{}\t|\t\t|\n".format(t, p, r))


class ScoreCamiTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.nodes = os.path.join(cls.tmp, "nodes.dmp")
        write_nodes(cls.nodes)
        # library with E. coli strain, K. pneumoniae and the archaeon
        cls.db = os.path.join(cls.tmp, "db")
        os.makedirs(os.path.join(cls.db, "library", "bacteria"))
        with open(os.path.join(cls.db, "library", "bacteria", "prelim_map.txt"), "w") as f:
            for t in (83333, 573, 2201):
                f.write("TAXID\tkraken:taxid|{0}|NC_{0}.1\t{0}\tdesc\turl\n".format(t))
        # truth: 10 reads
        reads = [
            ("r1", 83333), ("r2", 562), ("r3", 573), ("r4", 573), ("r5", 581),
            ("r6", 581), ("r7", 2201), ("r8", 777), ("r9", 99999), ("r10", 1),
        ]
        cls.truth = os.path.join(cls.tmp, "truth.tsv")
        with open(cls.truth, "w") as f:
            f.write("#seqid\ttaxid\tdomain\n")
            for sid, t in reads:
                f.write("{}/1\t{}\tall\n".format(sid, t))  # mate suffix is stripped
        with open(os.path.join(cls.tmp, "truth_extra.tsv"), "w") as f:
            f.write("#seqid\ttaxid\tgenome_id\tnovelty_category\n")
            for sid, t in reads:
                f.write("{}\t{}\tg\t{}\n".format(sid, t, "new_species" if sid == "r8" else "known_strain"))
        calls = {
            "r1": 562,    # strain read called at species: TP (genus+species)
            "r2": 561,    # genus call for a species read: TP genus, VP species
            "r3": 570,    # genus call: species VP, genus TP
            "r4": 561,    # wrong genus: FP
            "r5": 580,    # true genus Citrobacter not in library: excluded, consistent
            "r6": 0,      # excluded, unclassified
            "r7": 0,      # FN
            "r8": 543,    # family-level label, no genus node: 'no'
            "r9": 561,    # taxid not in nodes.dmp: 'no'
            "r10": 561,   # root taxid: 'no'
        }
        cls.kraken = os.path.join(cls.tmp, "out.kraken")
        with open(cls.kraken, "w") as f:
            for sid, t in calls.items():
                f.write("{}\t{}\t{}\t150|150\t0:1\n".format("C" if t else "U", sid, t))
        cls.parent, cls.ranks = load_taxonomy(cls.nodes)

    def run_score(self, db=None):
        out = os.path.join(self.tmp, "m.csv")
        args = ["--nodes", self.nodes, "--truth", self.truth, "--kraken", self.kraken,
                "--out-tsv", out, "--scheme", "T"]
        if db:
            args += ["--library-dir", db]
        score_cami.main(args)
        import csv
        with open(out, newline="") as f:
            return {(r["group"], r["rank"]): r for r in csv.DictReader(f)}

    def test_library_taxids(self):
        self.assertEqual(score_cami.library_taxids(self.db), {83333, 573, 2201})

    def test_genus_buckets(self):
        g = self.run_score(self.db)[("all", "genus")]
        # r1 TP, r2 TP, r3 TP, r4 FP, r7 FN scored; r5 r6 excluded; r8 r9 r10 no
        self.assertEqual((g["tp"], g["fp"], g["fn"], g["vp"]), ("3", "1", "1", "0"))
        self.assertEqual((g["excluded"], g["no"]), ("2", "3"))
        self.assertEqual((g["ex_consistent"], g["ex_unclassified"], g["ex_off"]), ("1", "1", "0"))
        self.assertAlmostEqual(float(g["sensitivity"]), 3 / 5.0)
        self.assertAlmostEqual(float(g["ppv"]), 3 / 4.0)
        self.assertAlmostEqual(float(g["excluded_frac"]), 0.2)

    def test_species_buckets(self):
        s = self.run_score(self.db)[("all", "species")]
        # in library at species: r1 (562 via strain), r2 (562), r3/r4 (573), r7 (2201)
        # r1 TP; r2 VP (genus call); r3 VP; r4 FP; r7 FN
        self.assertEqual((s["tp"], s["fp"], s["fn"], s["vp"]), ("1", "1", "1", "2"))
        # r5 r6: species 581 not in library -> excluded; r8 has no species node
        self.assertEqual(s["excluded"], "2")

    def test_without_library_nothing_is_excluded(self):
        g = self.run_score()[("all", "genus")]
        self.assertEqual(g["excluded"], "0")
        self.assertEqual(int(g["tp"]) + int(g["fp"]) + int(g["fn"]) + int(g["vp"]), 7)

    def test_groups(self):
        res = self.run_score(self.db)
        self.assertIn(("bacteria", "genus"), res)
        self.assertIn(("archaea", "genus"), res)
        self.assertIn(("novelty:new_species", "genus"), res)
        self.assertEqual(res[("archaea", "genus")]["fn"], "1")

    def test_coverage_only(self):
        self.assertEqual(score_cami.main(["--nodes", self.nodes, "--truth", self.truth,
                                          "--library-dir", self.db]), 0)

    def test_prepare_subsample_and_truth(self):
        d = tempfile.mkdtemp()
        reads = os.path.join(d, "anonymous_reads.fq.gz")
        mapping = os.path.join(d, "reads_mapping.tsv.gz")
        meta = os.path.join(d, "metadata.tsv")
        with gzip.open(reads, "wt") as f, gzip.open(mapping, "wt") as m:
            m.write("#anonymous_read_id\tgenome_id\ttax_id\tread_id\n")
            for i in range(20):
                for mate in (1, 2):
                    f.write("@S0R{}/{}\nACGT\n+\nIIII\n".format(i, mate))
                    m.write("S0R{}/{}\tOtu{}\t{}\tx\n".format(i, mate, i % 2, 562 if i % 2 else 573))
        with open(meta, "w") as f:
            f.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\nOtu0\t1\t573\tnew_strain\nOtu1\t2\t562\tknown_strain\n")
        out = os.path.join(d, "out")
        prepare_cami_reads.main(["--reads", reads, "--mapping", mapping, "--metadata", meta,
                                 "--out-dir", out, "--pairs", "5", "--nodes", self.nodes])
        with open(os.path.join(out, "truth.tsv")) as f:
            rows = [l.split("\t") for l in f.read().splitlines() if not l.startswith("#")]
        # stride = 20 // 5 = 4: pairs 0, 4, 8, 12, 16 (original IDs kept)
        self.assertEqual([r[0] for r in rows], ["S0R0", "S0R4", "S0R8", "S0R12", "S0R16"])
        self.assertEqual({r[1] for r in rows}, {"573"})
        self.assertEqual({r[2] for r in rows}, {"bacteria"})
        with open(os.path.join(out, "R2.fastq")) as f:
            self.assertTrue(f.readline().startswith("@S0R0/2"))
        with open(os.path.join(out, "truth_extra.tsv")) as f:
            self.assertIn("new_strain", f.read())

    def test_compare_summaries(self):
        rows = []
        for seed, f1 in ((1, 0.8), (2, 0.6)):
            rows.append({"scheme": "Random", "group": "all", "rank": "genus", "seed_idx": seed,
                         "sensitivity": 0.5, "ppv": 0.9, "f1": f1, "excluded_frac": 0.1, "n_scored": 100})
        fields, out = compare_benchmark.summarize_cami(rows)
        self.assertEqual(out[0]["n_seeds"], 2)
        self.assertAlmostEqual(out[0]["f1_mean"], 0.7)
        self.assertIn("excluded_frac_mean", fields)
        mrows = [{"scheme": "Random", "scope": "all", "rank": "species", "min_frac": 1e-4,
                  "seed_idx": s, **{m: 0.5 for m in compare_benchmark.MOCK_METRICS}} for s in (1, 2)]
        fields, out = compare_benchmark.summarize_mock(mrows)
        self.assertEqual(len(out), 1)
        self.assertAlmostEqual(out[0]["bray_curtis_mean"], 0.5)


if __name__ == "__main__":
    unittest.main()
