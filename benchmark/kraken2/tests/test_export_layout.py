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

ORDERS = ["Random", "RandomPcg", "GreedyMini_w5_k15_canonical", "ModRandom"]   # run directory names
SHORT = ["Random", "RandomPcg", "GreedyMini", "ModRandom"]               # displayed names


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
                "mrd_min": 1.0, "cls_pct": 77.0, "class_gb_std": 0.0, "n_classify": 3, "n_reads": 1000,
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
            base = {"scheme": E.family(order), "dataset": dataset, "seed_idx": s, "seed": "1" * 20 + "0" * 44, "context": ctx}
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
    if not os.path.isdir(run_root):
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

    def lib(self, name="strain_exclusion", datasets=("strain-exclusion", "fda-argos", "cami", "mock")):
        return make_library(self.base, name, list(datasets))

    def test_sheets_and_files(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        self.assertEqual(check_xlsx(xlsx), [])
        self.assertEqual(list(tables), list(E.SHEETS))
        self.assertEqual(list(load_xlsx(xlsx)), list(E.SHEETS))
        for name in tables:
            self.assertTrue(os.path.isfile(os.path.join(tdir, "kraken2_%s.csv" % name)), name)
        # no charts in the workbook
        import zipfile
        self.assertFalse([n for n in zipfile.ZipFile(xlsx).namelist() if "chart" in n or "drawing" in n])

    def test_build_once_per_library_order(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        build_runs, build = tables["build_runs"][1], tables["build"][1]
        self.assertEqual(len(build_runs), 4 * 3)          # order x seed, not x dataset
        self.assertEqual(len(build), 4)
        r0 = [r for r in build if r["order"] == "Random"][0]
        self.assertAlmostEqual(r0["build_wall_s_mean"], 102.0)   # seeds 1..3 -> 101,102,103
        self.assertAlmostEqual(r0["build_wall_s_std"], 1.0)
        self.assertEqual(r0["build_wall_s_ratio"], 1.0)
        pcg = [r for r in build if r["order"] == "RandomPcg"][0]
        self.assertAlmostEqual(pcg["build_wall_s_ratio"], round(202.0 / 102.0, 4))
        self.assertAlmostEqual(r0["expected_density_factor_mean"], 1.5)   # density seeds 1,2
        self.assertFalse([f for f in tables["build"][0] if "ns_per_minimizer" in f])
        self.assertEqual(sorted(r["order"] for r in build), sorted(SHORT))
        # classify / classify_runs: classify phase only
        for sheet in ("classify", "classify_runs"):
            for f in tables[sheet][0]:
                self.assertFalse(f.startswith(("build_", "expected_", "particular_", "rank_")), (sheet, f))

    def test_classify_and_classify_runs(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        runs, summ = tables["classify_runs"][1], tables["classify"][1]
        self.assertEqual(len(runs), 4 * 4 * 3)
        self.assertEqual(len(summ), 4 * 4)
        r0 = [r for r in summ if r["dataset"] == "strain-exclusion" and r["order"] == "Random"][0]
        self.assertAlmostEqual(r0["classify_wall_s_mean"], 2.21, places=6)
        self.assertEqual(r0["classify_wall_s_ratio"], 1.0)
        self.assertEqual(runs[0]["seed"][:2], "0x")
        self.assertIn("classify_wall_s_repeat_std", tables["classify_runs"][0])
        self.assertIn("n_classify_repeats", tables["classify_runs"][0])
        # no dataset prefix; strain-exclusion and cami share the read-level names
        self.assertFalse([f for f in tables["classify_runs"][0] if f.startswith(("se_", "argos_", "cami_", "mock_"))])
        self.assertAlmostEqual(r0["genus_sensitivity_mean"], 0.702, places=6)
        cami = [r for r in summ if r["dataset"] == "cami" and r["order"] == "Random"][0]
        self.assertEqual(cami["genus_f1_mean"], 0.7)
        self.assertEqual(cami["genus_excluded_fraction_mean"], 0.1)
        self.assertIsNone(cami.get("genus_concordant_fraction_mean"))
        mock = [r for r in summ if r["dataset"] == "mock" and r["order"] == "Random"][0]
        self.assertEqual(mock["species_taxa_f1_mean"], 0.64)
        self.assertIsNone(mock.get("genus_f1_mean"))   # taxa_f1 is not the read-level f1
        argos = [r for r in summ if r["dataset"] == "fda-argos" and r["order"] == "RandomPcg"][0]
        self.assertAlmostEqual(argos["genus_concordant_fraction_mean"], 0.87)
        # old cryptic names are gone everywhere
        gone = ("class_s", "class_gb", "mrd_min", "cls_pct", "mape_g", "nextmin_", "fac_", "build_s", "build_gb")
        for name, (fields, _rows) in tables.items():
            if name == "overview":
                continue
            for f in fields:
                self.assertFalse(f.startswith(gone), (name, f))

    def test_ratio_rules_and_bands(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields = tables["classify"][0]
        for neutral in ("classified_pct", "genus_ancestor_fraction", "genus_unclassified_fraction",
                        "genus_excluded_fraction", "species_excluded_fraction"):
            self.assertIn(neutral + "_mean", fields)
            self.assertNotIn(neutral + "_ratio", fields)
        for directional in ("genus_f1", "species_ppv", "genus_discordant_fraction", "species_bray_curtis",
                            "genus_taxa_f1", "genus_off_target_fraction", "classify_mreads_per_min"):
            self.assertIn(directional + "_ratio", fields)
        spec = E.sheet_spec("classify", fields, tables["classify"][1])
        self.assertEqual(spec["ratio_fields"]["genus_f1_ratio"], "higher")
        self.assertEqual(spec["ratio_fields"]["genus_discordant_fraction_ratio"], "lower")
        self.assertEqual(spec["ratio_fields"]["classify_wall_s_ratio"], "lower")
        self.assertNotIn("classified_pct_ratio", spec["ratio_fields"])
        cats = E._categories("classify", fields)
        self.assertEqual(cats[fields.index("genus_f1_mean")], E.BAND_READ)
        self.assertEqual(cats[fields.index("genus_excluded_fraction_mean")], E.BAND_CAMI)
        self.assertEqual(cats[fields.index("genus_concordant_fraction_mean")], E.BAND_ARGOS)
        self.assertEqual(cats[fields.index("species_bray_curtis_mean")], E.BAND_MOCK)
        # a dataset-less library drops the unused accuracy columns
        _rr, (_x, _t, t2) = make_library(self.base, "standard_p10", ["speed"])
        self.assertFalse([f for f in t2["classify"][0] if f.startswith(("genus_", "species_"))])

    def test_ratio_blank_when_reference_missing_or_zero(self):
        rows = [
            {"library": "L", "dataset": "d", "order": "Random", "classify_wall_s_mean": 0.0},
            {"library": "L", "dataset": "d", "order": "RandomPcg", "classify_wall_s_mean": 5.0},
            {"library": "L", "dataset": "e", "order": "RandomPcg", "classify_wall_s_mean": 5.0},  # no reference
            {"library": "L", "dataset": "f", "order": "Random", "classify_wall_s_mean": 4.0},
            {"library": "L", "dataset": "f", "order": "ModRandom", "classify_wall_s_mean": None},
        ]
        E.attach_ratios(rows, ["classify_wall_s", "classified_pct"])
        self.assertIsNone(rows[0]["classify_wall_s_ratio"])   # reference 0: blank, not an error
        self.assertIsNone(rows[1]["classify_wall_s_ratio"])
        self.assertIsNone(rows[2]["classify_wall_s_ratio"])
        self.assertEqual(rows[3]["classify_wall_s_ratio"], 1.0)
        self.assertIsNone(rows[4]["classify_wall_s_ratio"])
        self.assertNotIn("classified_pct_ratio", rows[0])

    def test_profile_vs_reference_detail(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        prof = tables["profile"][1]
        self.assertEqual(len(prof), 4 * 4)
        pr = [r for r in prof if r["dataset"] == "strain-exclusion" and r["order"] == "ModRandom"][0]
        self.assertAlmostEqual(pr["nextminimizer_s_ratio"], round(0.33 / 0.3, 4))
        vr = tables["vs_reference"][1]
        self.assertEqual(len(vr), 3 * 4)                       # no reference order
        self.assertIn("jaccard_pct_mean", tables["vs_reference"][0])
        # one accuracy sheet per dataset, aggregated over seeds, only the columns that apply
        self.assertEqual(list(tables)[-4:], ["accuracy_strain_exclusion", "accuracy_fda_argos",
                                             "accuracy_cami", "accuracy_mock"])
        self.assertNotIn("accuracy_detail", tables)
        for sheet, dataset, has, hasnt in (
                ("accuracy_strain_exclusion", "strain-exclusion", "sensitivity_mean", "concordant_fraction_mean"),
                ("accuracy_fda_argos", "fda-argos", "concordant_fraction_mean", "sensitivity_mean"),
                ("accuracy_cami", "cami", "excluded_fraction_mean", "concordant_fraction_mean"),
                ("accuracy_mock", "mock", "taxa_f1_mean", "sensitivity_mean")):
            fields, rows = tables[sheet]
            self.assertEqual(set(r["dataset"] for r in rows), set([dataset]))
            self.assertIn(has, fields)
            self.assertNotIn(hasnt, fields)
            self.assertIn("n_seeds", fields)
            self.assertNotIn("seed_idx", fields)
            self.assertNotIn("seed", fields)
            for f in fields:   # no column that is blank for the whole dataset
                self.assertTrue(f in ("library", "dataset", "order", "n_seeds")
                                or any(r.get(f) not in (None, "") for r in rows), (sheet, f))
            self.assertEqual(sorted(set(r["order"] for r in rows)), sorted(SHORT))
        self.assertIn("min_frac", tables["accuracy_mock"][0])
        se = tables["accuracy_strain_exclusion"][1]
        self.assertEqual(len(se), 4 * 4)                    # order x (domain x rank): seeds collapsed
        self.assertEqual(se[0]["n_seeds"], 3)

    def test_accuracy_aggregation(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields, rows = tables["accuracy_strain_exclusion"]
        # slice first, order last: the 4 orders of a slice are consecutive, Random first
        for i in range(0, len(rows), 4):
            chunk = rows[i:i + 4]
            self.assertEqual(len(set((r["domain"], r["rank"]) for r in chunk)), 1)
            self.assertEqual([r["order"] for r in chunk], SHORT)
        # mean of the per-seed values (0.7 - 0.01 * oi + 0.001 * s, s = 1..3) and ratio to Random
        r = [x for x in rows if x["order"] == "ModRandom" and x["domain"] == "all" and x["rank"] == "genus"][0]
        ref = [x for x in rows if x["order"] == "Random" and x["domain"] == "all" and x["rank"] == "genus"][0]
        self.assertAlmostEqual(ref["sensitivity_mean"], 0.702, places=6)
        self.assertAlmostEqual(r["sensitivity_mean"], 0.672, places=6)
        self.assertAlmostEqual(r["sensitivity_std"], 0.001, places=6)
        self.assertAlmostEqual(r["sensitivity_ratio"], round(0.672 / 0.702, 4))
        self.assertEqual(ref["sensitivity_ratio"], 1.0)
        # column order: mean, std, ratio adjacent per metric; counts have only a mean
        i = fields.index("sensitivity_mean")
        self.assertEqual(fields[i:i + 3], ["sensitivity_mean", "sensitivity_std", "sensitivity_ratio"])
        # directions / neutral
        cf, crows = tables["accuracy_cami"]
        self.assertIn("f1_ratio", cf)
        self.assertNotIn("excluded_fraction_ratio", cf)
        af, _ar = tables["accuracy_fda_argos"]
        self.assertIn("discordant_fraction_ratio", af)
        self.assertNotIn("ancestor_fraction_ratio", af)
        self.assertNotIn("unclassified_fraction_ratio", af)
        mf, _mr = tables["accuracy_mock"]
        self.assertIn("bray_curtis_ratio", mf)
        self.assertIn("off_target_fraction_ratio", mf)
        self.assertNotIn("top_off_target", " ".join(mf))
        self.assertNotIn("tp_ratio", " ".join(af + mf))
        # spec: merged slice columns, borders per slice, direction of the colour scale
        spec = E.sheet_spec("accuracy_strain_exclusion", fields, rows)
        self.assertEqual(spec["ratio_fields"]["sensitivity_ratio"], "higher")
        self.assertEqual(E.sheet_spec("accuracy_mock", mf, _mr)["ratio_fields"]["bray_curtis_ratio"], "lower")
        self.assertIn("domain", spec["block_fields"])
        self.assertIn("rank", spec["merge_fields"])
        self.assertNotIn("order", spec["merge_fields"])

    def test_build_measured_compact(self):
        a = "previous session (db built 2026-08-18 with the earlier binaries on a differently loaded host; x)"
        b = a.replace("08-18", "08-19")
        self.assertEqual(E.short_measured(a), "previous session (2026-08-18)")
        self.assertEqual(E.join_measured([E.short_measured(x) for x in (a, a, b, a)]),
                         "previous session (2026-08-18/19)")
        self.assertEqual(E.join_measured(["this run"] * 5), "this run")
        self.assertEqual(E.join_measured(["this run", "previous session (2026-08-31)", "previous session (2026-09-01)"]),
                         "this run; previous session (2026-08-31/2026-09-01)")
        self.assertIsNone(E.join_measured([None, ""]))

    def test_orders_density(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields, orders = tables["orders"]
        self.assertIn("expected_density_factor_mean", fields)
        self.assertIn("strain_exclusion_particular_density_factor_mean", fields)
        o = dict((r["order"], r) for r in orders)
        self.assertAlmostEqual(o["Random"]["expected_density_factor_mean"], 1.5)
        self.assertAlmostEqual(o["Random"]["expected_density_factor_std"], 0.0)
        self.assertAlmostEqual(o["RandomPcg"]["strain_exclusion_particular_density_factor_mean"], 1.7)

    def test_merge_libraries_and_stale_files(self):
        make_library(self.base, "standard_p10", ["speed"])
        self.lib("strain_exclusion", ("strain-exclusion", "cami"))
        # previous layout (accuracy_detail sheet): skipped with a message, no crash
        prev = os.path.join(self.base, "runs", "standard_prev", "tables")
        os.makedirs(prev)
        for n in ("classify_runs", "accuracy_detail"):
            with open(os.path.join(prev, "kraken2_%s.csv" % n), "w") as f:
                f.write("library,dataset,order,classify_wall_s\nstandard_p0.5,speed,GreedyMini_w5_k15_canonical,1\n")
        # a library exported by an earlier layout is skipped, not mixed in
        old = os.path.join(self.base, "runs", "standard", "tables")
        os.makedirs(old)
        with open(os.path.join(old, "kraken2_runs.csv"), "w") as f:
            f.write("library,dataset,order,class_s\nstandard,speed,Random,1\n")
        # ... and so is one exported with the previous sheet names (`summary` / `runs`)
        old2 = os.path.join(self.base, "runs", "standard_p0.5", "tables")
        os.makedirs(old2)
        with open(os.path.join(old2, "kraken2_runs.csv"), "w") as f:
            f.write("library,dataset,order,classify_wall_s\nstandard_p0.5,speed,Random,1\n")
        with open(os.path.join(old2, "kraken2_summary.csv"), "w") as f:
            f.write("library,dataset,order,classify_wall_s_mean\nstandard_p0.5,speed,Random,1\n")
        # db dirs of libraries without runs never list
        os.makedirs(os.path.join(self.base, "db", "standard", "library"))
        os.makedirs(os.path.join(self.base, "db", "standard_p0.5", "library"))
        stale = os.path.join(self.base, "results")
        if not os.path.isdir(stale):
            os.makedirs(stale)
        with open(os.path.join(stale, "kraken2_cami_summary.csv"), "w") as f:
            f.write("x\n")
        xlsx, rdir, n = E.export_results(self.base)
        self.assertEqual(n, 2)
        self.assertFalse(os.path.exists(os.path.join(rdir, "kraken2_cami_summary.csv")))
        self.assertEqual(check_xlsx(xlsx), [])
        _f, summ = E.read_csv(os.path.join(rdir, "kraken2_classify.csv"))
        libs = [r["library"] for r in summ]
        self.assertEqual(libs, sorted(libs, key=E.library_key))
        self.assertEqual(set(libs), set(["standard_p10", "strain_exclusion"]))
        _f, lrows = E.read_csv(os.path.join(rdir, "kraken2_libraries.csv"))
        self.assertEqual(sorted(r["library"] for r in lrows), ["standard_p10", "strain_exclusion"])
        _f, orows = E.read_csv(os.path.join(rdir, "kraken2_overview.csv"))
        for name in ("build", "build_runs", "classify", "classify_runs"):
            self.assertIn(name, [r["column"] for r in orows if r["sheet"] == "sheets"])
        self.assertFalse([f for f in os.listdir(rdir) if f in ("kraken2_summary.csv", "kraken2_runs.csv")])
        for r in summ:
            if r["order"] == "Random":
                self.assertEqual(r["classify_wall_s_ratio"], 1.0)
        _f, build = E.read_csv(os.path.join(rdir, "kraken2_build.csv"))
        self.assertEqual(len(build), 2 * 4)
        _f, runs = E.read_csv(os.path.join(rdir, "kraken2_classify_runs.csv"))
        self.assertEqual(len(runs), 3 * 4 * 3)
        _f, ds = E.read_csv(os.path.join(rdir, "kraken2_datasets.csv"))
        self.assertIn("speed", [d["dataset"] for d in ds])
        self.assertEqual(list(load_xlsx(xlsx)), [s for s in E.SHEETS if s in load_xlsx(xlsx)])

    def test_build_measured_provenance(self):
        run_root = os.path.join(self.base, "runs", "strain_exclusion")
        d = os.path.join(run_root, "RandomPcg", "2", "results")
        os.makedirs(d)
        with open(os.path.join(d, "build_measured.txt"), "w") as f:
            f.write("previous session (2026-08-18)\n")
        _rr, (_x, _t, tables) = self.lib("strain_exclusion", ("strain-exclusion",))
        br = tables["build_runs"][1]
        pcg2 = [r for r in br if r["order"] == "RandomPcg" and r["seed_idx"] == 2][0]
        self.assertEqual(pcg2["build_measured"], "previous session (2026-08-18)")
        rnd = [r for r in br if r["order"] == "Random"][0]
        self.assertEqual(rnd["build_measured"], "this run")
        b = [r for r in tables["build"][1] if r["order"] == "RandomPcg"][0]
        self.assertEqual(b["build_measured"], "this run; previous session (2026-08-18)")

    def test_overview_documents_every_column(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        import re
        documented = " ".join(r["column"] for r in E.overview_rows())
        for sheet, (fields, _rows) in tables.items():
            if sheet == "overview":
                continue
            for f in fields:
                base = re.sub(r"_(mean|std|ratio)$", "", f)
                base = re.sub(r"^.*_particular_density_factor", "<library>_particular_density_factor", base)
                self.assertTrue(base in documented or f in documented, "%s.%s undocumented" % (sheet, f))

    def test_short_order_names(self):
        self.assertEqual(E.short_order("GreedyMini_w5_k15_canonical"), "GreedyMini")
        self.assertEqual(E.short_order("ModGreedyMini_w5_k6_canonical"), "ModGreedyMini")
        self.assertEqual(E.short_order("ModRandom"), "ModRandom")
        self.assertEqual(E.short_order("GreedyMini"), "GreedyMini")      # already short
        self.assertEqual(E.gm_table("ModGreedyMini_w5_k6_canonical"), "w5_k6_canonical")
        self.assertIsNone(E.gm_table("Random"))
        self.assertEqual(E.family("GreedyMini"), "GreedyMini")
        # two tables of one scheme in a library keep their long names (no silent mixing)
        m = E.display_orders(["GreedyMini_w5_k15_canonical", "GreedyMini_w4_k15_canonical", "Random"])
        self.assertEqual(m["GreedyMini_w5_k15_canonical"], "GreedyMini_w5_k15_canonical")
        self.assertEqual(m["Random"], "Random")

    def test_orders_sheet_and_joins(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        orders = tables["orders"][1]
        self.assertEqual([o["order"] for o in orders], SHORT)
        self.assertEqual(orders[2]["gm_table"], "w5_k15_canonical")
        self.assertIn("w5_k15_canonical", orders[2]["cli"])
        self.assertNotIn("table", tables["orders"][0])
        # density rows keyed by family ("GreedyMini"), by long name and by short name all join
        for key in ("GreedyMini", "GreedyMini_w5_k15_canonical"):
            dens = [dict(d, scheme=key) if d["scheme"] == "GreedyMini" else d for d in DENS]
            ps = per_seed("L", "speed")
            raw = E.build_tables("L", ps, {}, dens)
            b = [r for r in raw["build_runs"] if r["order"] == "GreedyMini" and r["seed_idx"] == 1][0]
            self.assertAlmostEqual(b["expected_density_factor"], 1.7)

    def test_round_std(self):
        self.assertEqual(E.round_std(1.23456e-5), 1.235e-05)
        self.assertEqual(E.round_std(0.0), 0.0)
        self.assertEqual(E.round_std(123.456789), 123.4568)


if __name__ == "__main__":
    unittest.main()
