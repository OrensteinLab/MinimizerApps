"""Export layer: sheet layout, summary / ratio maths, per-library export and the
multi-library merge, on synthetic per-run numbers (no database needed).

    python3 -m unittest discover -s benchmark/kraken2/tests
"""

import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = os.path.dirname(HERE)
sys.path.insert(0, K2)
sys.path.insert(0, os.path.join(os.path.dirname(K2), "common"))

import export_layout as E  # noqa: E402
from xlsx_check import check_xlsx  # noqa: E402
from xlsx_util import load_xlsx  # noqa: E402

ORDERS = ["Random", "RandomPcg", "GreedyMini_w5_k15_canonical", "ModRandom"]


def per_seed(library, dataset, seeds=(1, 2, 3)):
    rows = []
    for oi, order in enumerate(ORDERS):
        for s in seeds:
            walls = [2.0 + oi + 0.1 * s + 0.01 * i for i in range(3)]
            rows.append({
                "scheme": E.family(order), "dataset": dataset, "seed_idx": s,
                "context": "%s/%s/%d [%s]" % (library, order, s, dataset),
                "seed": "1" * 20 + "0" * 44,
                "build_s": 100.0 * (oi + 1) + s, "build_gb": 10.0 + oi,
                "class_s": sum(walls) / 3, "class_s_std": 0.01, "class_gb": 9.0 + oi,
                "mrd_min": 1.0, "cls_pct": 77.0, "n_classify": 3, "n_reads": 1000,
                "class_walls": walls,
                "nextmin_s": 0.3 + 0.01 * oi, "lookups_pm": 2e5,
                "concordance": None if oi == 0 else 95.0 + oi, "jaccard": None if oi == 0 else 98.0,
            })
    return rows


def acc_rows(library, dataset, kind):
    rows = []
    for oi, order in enumerate(ORDERS):
        for s in (1, 2, 3):
            ctx = "%s/%s/%d" % (library, order, s)
            base = {"scheme": E.family(order), "dataset": dataset, "seed_idx": s, "seed": "x", "context": ctx}
            if kind == "labeled":
                for dom in ("all", "virus"):
                    for rank in ("genus", "species"):
                        r = dict(base, domain=dom, rank=rank, sensitivity=0.7 - 0.01 * oi + 0.001 * s,
                                 ppv=0.9, f1=0.8)
                        rows.append(r)
            elif kind == "argos":
                rows.append(dict(base, rank="genus", scope="unweighted_mean", sra_run="",
                                 concordant_frac=0.88 - 0.01 * oi, discordant_frac=0.04,
                                 ancestor_frac=0.03, unclassified_frac=0.05))
                rows.append(dict(base, rank="genus", scope="experiment", sra_run="SRR1",
                                 concordant_frac=0.9, discordant_frac=0.04, ancestor_frac=0.03,
                                 unclassified_frac=0.03))
            elif kind == "cami":
                for rank in ("genus", "species"):
                    rows.append(dict(base, group="all", rank=rank, sensitivity=0.6, ppv=0.9, f1=0.7,
                                     excluded_frac=0.1, n_scored=900))
            elif kind == "mock":
                for rank in ("genus", "species"):
                    rows.append(dict(base, scope="all", rank=rank, min_frac=1e-4, precision=0.5,
                                     recall=0.9, f1=0.64, l1=0.3, bray_curtis=0.15,
                                     classified_frac=0.9, on_target_frac=0.8, off_target_frac=0.1))
    return rows


DENS = [{"scheme": E.family(o), "seed_idx": s, "l": 31.0, "k": 35.0, "w": 5.0,
         "fac_expected": 1.5 + 0.1 * i, "fac_particular": 1.6 + 0.1 * i, "ns_per_minimizer": 90.0}
        for i, o in enumerate(ORDERS) for s in (1, 2)]
QUERY = [{"dataset": "strain-exclusion", "sources": "held-out genomes (Mason)", "n_sources": 50,
          "n_fragments": 1000, "read_length_bp": 100, "paired": "yes", "path": "/x"}]


def make_library(base, lib, datasets):
    run_root = os.path.join(base, "runs", lib)
    os.makedirs(run_root)
    ps, acc = [], {}
    for ds in datasets:
        ps += per_seed(lib, ds)
        kind = {"strain-exclusion": "labeled", "fda-argos": "argos", "cami": "cami", "mock": "mock"}.get(ds)
        if kind:
            acc.setdefault(kind, []).extend(acc_rows(lib, ds, kind))
    out = E.export_library(run_root, os.path.join(run_root, "benchmark_comparison.xlsx"), lib, ps, acc, DENS, QUERY)
    return run_root, out


class ExportLayout(unittest.TestCase):
    def setUp(self):
        self.base = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.base, ignore_errors=True)

    def test_library_export(self):
        run_root, (xlsx, tdir, tables) = make_library(
            self.base, "strain_exclusion", ["strain-exclusion", "fda-argos", "cami", "mock"])
        self.assertEqual(check_xlsx(xlsx), [])
        names = list(load_xlsx(xlsx))
        self.assertEqual(names[:6], list(E.CORE_SHEETS))
        for name in names:
            self.assertTrue(os.path.isfile(os.path.join(tdir, "kraken2_%s.csv" % name)), name)
        runs = tables["runs"][1]
        self.assertEqual(len(runs), 4 * 4 * 3)
        summ = tables["summary"][1]
        self.assertEqual(len(summ), 4 * 4)
        r0 = [r for r in summ if r["dataset"] == "strain-exclusion" and r["order"] == "Random"][0]
        self.assertEqual(r0["build_s_ratio"], 1.0)
        self.assertAlmostEqual(r0["build_s_mean"], 102.0)  # seeds 1..3 -> 101,102,103
        self.assertAlmostEqual(r0["build_s_std"], 1.0)
        pcg = [r for r in summ if r["dataset"] == "strain-exclusion" and r["order"] == "RandomPcg"][0]
        self.assertAlmostEqual(pcg["build_s_ratio"], round(202.0 / 102.0, 4))
        self.assertAlmostEqual(r0["se_sens_genus_mean"], 0.702, places=6)
        self.assertIsNone(r0.get("cami_f1_genus_mean"))
        # accuracy columns only filled for their own dataset
        cami = [r for r in summ if r["dataset"] == "cami" and r["order"] == "Random"][0]
        self.assertEqual(cami["cami_f1_genus_mean"], 0.7)
        self.assertEqual(cami.get("se_f1_genus_mean"), None)
        mock = [r for r in summ if r["dataset"] == "mock" and r["order"] == "Random"][0]
        self.assertEqual(mock["mock_f1_species_mean"], 0.64)
        for sheet in ("cami_summary", "mock_summary", "fda_argos_summary", "strain_exclusion_summary"):
            self.assertIn(sheet, tables)
        self.assertEqual(len(tables["classify_repeats"][1]), 4 * 4 * 3 * 3)
        # density joined per seed (seeds 1,2 only), mean over those
        self.assertAlmostEqual(r0["fac_expected_mean"], 1.5)
        # std rounding keeps >= 4 significant digits
        self.assertEqual(E.round_std(1.23456e-5), 1.235e-05)
        self.assertEqual(E.round_std(0.0), 0.0)
        self.assertEqual(E.round_std(123.456789), 123.4568)
        # orders: one row per order with params
        orders = tables["orders"][1]
        self.assertEqual([o["order"] for o in orders], ORDERS)
        self.assertEqual((orders[0]["l"], orders[0]["k"], orders[0]["w"]), (31, 35, 5))
        self.assertEqual(orders[2]["table"], "w5_k15_canonical")
        self.assertIn("strain_exclusion_particular_density_factor_mean", tables["orders"][0])
        # legacy tables removed, results merged workbook written by export_results
        res = E.export_results(self.base)
        self.assertEqual(check_xlsx(res[0]), [])

    def test_merge_libraries(self):
        make_library(self.base, "standard_p10", ["speed"])
        make_library(self.base, "strain_exclusion", ["strain-exclusion", "cami"])
        xlsx, rdir, n = E.export_results(self.base)
        self.assertEqual(n, 2)
        self.assertEqual(check_xlsx(xlsx), [])
        _f, summ = E.read_csv(os.path.join(rdir, "kraken2_summary.csv"))
        libs = [r["library"] for r in summ]
        self.assertEqual(libs, sorted(libs, key=E.library_key))
        self.assertEqual(set(libs), set(["standard_p10", "strain_exclusion"]))
        # ratio is per library x dataset
        for r in summ:
            if r["order"] == "Random":
                self.assertEqual(r["build_s_ratio"], 1.0)
        _f, runs = E.read_csv(os.path.join(rdir, "kraken2_runs.csv"))
        self.assertEqual(len(runs), 3 * 4 * 3)
        self.assertEqual(runs[0]["seed"][:2], "0x")
        _f, ds = E.read_csv(os.path.join(rdir, "kraken2_datasets.csv"))
        self.assertIn("speed", [d["dataset"] for d in ds])
        sheets = load_xlsx(xlsx)
        self.assertEqual(list(sheets)[:6], list(E.CORE_SHEETS))

    def test_overview_documents_every_column(self):
        run_root, (xlsx, tdir, tables) = make_library(
            self.base, "strain_exclusion", ["strain-exclusion", "fda-argos", "cami", "mock"])
        documented = " ".join(r["column"] for r in E.overview_rows())
        import re
        for sheet in ("runs", "summary", "orders", "libraries", "datasets", "vs_reference",
                      "timing_correlation", "classify_repeats"):
            for f in tables[sheet][0]:
                base = re.sub(r"_(mean|std|ratio)$", "", f)
                base = re.sub(r"^[a-z_0-9.]+_particular_density_factor", "<library>_particular_density_factor", base)
                self.assertTrue(base in documented or f in documented, "%s.%s undocumented" % (sheet, f))


if __name__ == "__main__":
    unittest.main()
