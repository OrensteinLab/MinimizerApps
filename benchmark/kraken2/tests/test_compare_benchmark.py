"""compare_benchmark: single-run classify metrics (DB load vs processing), the
<dataset>_mhg<m> file layout and the published paper values.

    python3 -m unittest discover -s benchmark/kraken2/tests
"""

import csv
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = os.path.dirname(HERE)
sys.path.insert(0, K2)
sys.path.insert(0, os.path.join(os.path.dirname(K2), "common"))

import compare_benchmark as C  # noqa: E402
import export_layout as E  # noqa: E402

PERF = """Loading database information... done.
20000 sequences (4.00 Mbp) processed in 0.236s (5086.9 Kseq/m, 1017.38 Mbp/m).
  3806 sequences classified (19.03%)
  16194 sequences unclassified (80.97%)
\tCommand being timed: "k2 classify"
\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:00.49
\tMaximum resident set size (kbytes): 299424
"""


class Compare(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def write(self, name, text):
        with open(os.path.join(self.d, name), "w") as f:
            f.write(text)

    def test_single_run_metrics(self):
        self.write("classification_performance_cami_mhg0.txt", PERF)
        self.write("classification_performance_strain-exclusion_mhg2.txt", PERF)
        self.write("classification_performance_cami_mhg2.txt", PERF)
        self.assertEqual(C.discover_classify_cells(self.d),
                         [("strain-exclusion", 2), ("cami", 0), ("cami", 2)])
        run = {"outdir": self.d, "seed_idx": 1, "context": "lib/xor/1"}
        m = C.collect_run_metrics(run, 1e-6, {1: run}, "cami", 0)
        self.assertAlmostEqual(m["class_s"], 0.49)
        self.assertAlmostEqual(m["class_proc_s"], 0.236)
        self.assertAlmostEqual(m["db_load_s"], 0.254)
        self.assertAlmostEqual(m["class_gb"], 299424 * 1024 / 1e9, places=4)
        self.assertAlmostEqual(m["cls_pct"], 19.03)
        self.assertAlmostEqual(m["mrd_min"], 0.02 / (0.236 / 60.0))
        self.assertEqual(m["minimum_hit_groups"], 0)
        self.assertIsNone(m["build_pinned_cpus"])
        self.write("cpus_build.txt", "0-5,64-69\n")
        self.write("classification_performance_cami_mhg0_cpus.txt", "0,64\n")
        m = C.collect_run_metrics(run, 1e-6, {1: run}, "cami", 0)
        self.assertEqual((m["build_pinned_cpus"], m["classify_pinned_cpus"]), ("0-5,64-69", "0,64"))
        self.assertFalse([k for k in m if "std" in k or k in ("n_classify", "class_walls")])
        # a PROF-less run has no profile
        self.assertIsNone(m["scan_s"])
        # missing file: nothing, no crash
        m = C.collect_run_metrics(run, 1e-6, {1: run}, "nodata", 0)
        self.assertIsNone(m["class_s"])
        self.assertIsNone(m["db_load_s"])

    def test_manifest_five_columns_and_cell_keys(self):
        man = os.path.join(self.d, "m")
        with open(man, "w") as f:
            f.write("xor-pcg-hash\tlib/xor-pcg-hash/1\t%s\t0101\t1\n" % self.d)
            f.write("xor\tlib/xor/1\t%s\t0101\t1\n" % self.d)
        runs = C.load_manifest(man)
        self.assertEqual([r["scheme"] for r in runs], ["xor", "xor-pcg-hash"])
        self.assertEqual(runs[0]["seed_idx"], 1)
        key = C._cell_key("xor-pcg-hash", "fda-argos", 2)
        self.assertEqual(C._scheme_cell(key), ("fda-argos", 2))
        self.assertEqual(C._scheme_base(key), "xor-pcg-hash")

    def test_manifest_identity_columns_and_variant_labels(self):
        man = os.path.join(self.d, "m")
        with open(man, "w") as f:
            f.write("greedymini\tlib/greedymini/k15_at0-cg1/1\t%s\t0101\t1\t15\tat0-cg1\tw5_k15_canonical\t\n" % self.d)
            f.write("greedymini\tlib/greedymini/k15_at1-cg0/1\t%s\t0101\t1\t15\tat1-cg0\tw5_k15_canonical\t\n" % self.d)
            f.write("mod-minimizer-greedymini\tlib/mod-minimizer-greedymini/k6_at0-cg1_t6/1\t%s\t0101\t1\t6\tat0-cg1\tw5_k6_canonical\t6\n" % self.d)
            f.write("xor\tlib/xor/1\t%s\t0101\t1\t\t\t\t\n" % self.d)
        runs = C.load_manifest(man)
        # one parameterisation per order: plain order label; several: "<order> (<params>)"
        self.assertEqual([r["scheme"] for r in runs],
                         ["xor", "greedymini (k15_at0-cg1)", "greedymini (k15_at1-cg0)",
                          "mod-minimizer-greedymini"])
        g = runs[2]
        self.assertEqual((g["order"], g["gm_k"], g["projection"], g["gm_table"], g["mod_t"]),
                         ("greedymini", 15, "at1-cg0", "w5_k15_canonical", None))
        self.assertEqual(runs[3]["mod_t"], 6)
        self.assertEqual(C.scheme_rank(runs[2]["scheme"]), 2)

    def test_old_scheme_names_rejected(self):
        man = os.path.join(self.d, "m")
        with open(man, "w") as f:
            f.write("pcg\tlib/pcg/1\t%s\t0101\t1\n" % self.d)
        with self.assertRaises(ValueError):
            C.load_manifest(man)

    def test_published_values_match_tables(self):
        pub = dict(((p["table"], p["domain"], p["metric"]), p["value"]) for p in E.read_published())
        with open(os.path.join(K2, "tables", "table_s1.csv")) as f:
            k2 = [r for r in csv.DictReader(f) if r["classifier"] == "Kraken 2"][0]
        self.assertAlmostEqual(pub[("S1", "prokaryote", "f1")], float(k2["prok_genus_f1"]))
        self.assertAlmostEqual(pub[("S1", "virus", "sensitivity")], float(k2["viral_genus_sensitivity"]))
        with open(os.path.join(K2, "tables", "table_s5_kraken2_expected.csv")) as f:
            rows = list(csv.DictReader(f))
        for m in ("concordant", "discordant", "ancestor", "unclassified"):
            mean = sum(float(r["genus_" + m]) for r in rows) / len(rows)
            self.assertAlmostEqual(pub[("S5", None, m)], mean, places=5)

    def test_density_csv_and_data_dir_in_new_layout(self):
        root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, root, True)
        run = os.path.join(root, "results", "runs", "strain_exclusion")
        os.makedirs(run)
        os.makedirs(os.path.join(root, "data"))
        with open(os.path.join(run, "minimizer_density.csv"), "w") as f:
            f.write("scheme,seed_idx,k,w,L,fac_expected,fac_particular\nxor,1,15,5,19,1.5,1.6\n")
        rows = C.load_density_rows(run)
        self.assertEqual([(r["scheme"], r["seed_idx"], r["fac_expected"]) for r in rows], [("xor", 1, 1.5)])
        old = dict(os.environ)
        try:
            os.environ.pop("KRAKEN2_DATA_DIR", None)
            os.environ.pop("KRAKEN2_BENCH_ROOT", None)
            self.assertEqual(os.path.realpath(C.resolve_data_dir(run)), os.path.realpath(os.path.join(root, "data")))
        finally:
            os.environ.clear()
            os.environ.update(old)


if __name__ == "__main__":
    unittest.main()
