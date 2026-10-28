"""Export layer: sheet layout, summary / ratio maths, per-library export and the
multi-library merge, on synthetic per-run numbers (no database needed).

    python3 -m unittest discover -s benchmark/kraken2/tests
"""

import os
import re
import shutil
from collections import OrderedDict
import struct
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = os.path.dirname(HERE)
sys.path.insert(0, K2)
sys.path.insert(0, os.path.join(os.path.dirname(K2), "common"))
import step_record  # noqa: E402

import export_layout as E  # noqa: E402
import order_names  # noqa: E402
from xlsx_check import check_xlsx  # noqa: E402
from xlsx_util import load_xlsx  # noqa: E402

ORDERS = ["xor", "xor-pcg-hash", "greedymini/k15_at0-cg1", "mod-minimizer-xor/t6"]   # run directory names
SHORT = ["xor", "xor-pcg-hash", "greedymini", "mod-minimizer-xor"]                  # canonical `order` values
TABLES = {"greedymini": "w5_k15_canonical"}


def ident(order_dir):
    """The five identity fields of a run directory + its table (as the manifest columns)."""
    idd = order_names.parse_order_dir(order_dir)
    idd["gm_table"] = TABLES.get(idd["order"])
    return idd
MHGS = (0, 2)


def per_seed(library, dataset, seeds=(1, 2, 3)):
    rows = []
    for mhg in MHGS:
        for oi, order in enumerate(ORDERS):
            for sd in seeds:
                wall = 2.0 + oi + 0.1 * sd + 10 * mhg
                rows.append(dict(ident(order), **{
                    "scheme": order.split("/")[0], "dataset": dataset, "minimum_hit_groups": mhg,
                    "seed_idx": sd,
                    "context": "%s/%s/%d [%s, mhg%d]" % (library, order, sd, dataset, mhg),
                    "seed": "1" * 20 + "0" * 44,
                    "build_s": 100.0 * (oi + 1) + sd, "build_gb": 10.0 + oi,
                    "class_s": wall, "db_load_s": 0.5 * wall, "class_proc_s": 0.5 * wall,
                    "class_gb": 9.0 + oi, "mrd_min": 1.0, "cls_pct": 77.0 - mhg, "n_reads": 1000,
                    "scan_s": 0.3 + 0.01 * oi, "lookups_pm": 2e5,
                    "concordance": None if oi == 0 else 95.0 + oi, "jaccard": None if oi == 0 else 98.0,
                }))
    return rows


def acc_rows(library, dataset, kind):
    rows = []
    for mhg, oi, order, s in [(m, o, ORDERS[o], sd) for m in MHGS for o in range(len(ORDERS)) for sd in (1, 2, 3)]:
        if True:
            ctx = "%s/%s/%d" % (library, order, s)
            base = {"scheme": order.split("/")[0], "dataset": dataset, "minimum_hit_groups": mhg,
                    "seed_idx": s, "seed": "1" * 20 + "0" * 44, "context": ctx}
            if kind == "labeled":
                for dom in ("all", "prokaryote", "virus"):
                    for rank in ("genus", "species"):
                        r = dict(base, domain=dom, rank=rank, sensitivity=0.7 - 0.01 * oi + 0.001 * s - 0.1 * mhg,
                                 ppv=0.9, f1=0.8, n_truth=1000, n_calls=900, tp=700, fp=50, fn=100, vp=150,
                                 no=0, n=1000)
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
                    for grp, tp in (("all", 90), ("prokaryote", 80), ("archaea", 10), ("bacteria", 70), ("virus", 4)):
                        t = tp + s
                        rows.append(dict(base, group=grp, rank=rank, n_truth=200, n_calls=150, no=2, excluded=20,
                                         excluded_frac=0.1, ex_unclassified=15, ex_consistent=2, ex_off=3,
                                         n_scored=t + 20 + 10 + 5, tp=t, fp=10, fn=20, vp=5,
                                         sensitivity=0.6, ppv=0.9, f1=0.7))
    return rows


DENS = [dict(dict((k, v) for k, v in ident(o).items() if k != "order"), scheme=o.split("/")[0],
             seed_idx=s, k=31.0, w=5.0, L=35.0, fac_expected=1.5 + 0.1 * i,
             fac_particular=1.6 + 0.1 * i, ns_per_minimizer=90.0)
        for i, o in enumerate(ORDERS) for s in (1, 2)]
QUERY = [{"dataset": "strain-exclusion", "sources": "held-out genomes (Mason)", "n_sources": 50,
          "n_fragments": 1000, "read_length_bp": 100, "paired": "yes", "path": "/x"}]


def make_library(base, lib, datasets):
    run_root = os.path.join(base, "results", "runs", lib)
    if not os.path.isdir(run_root):
        os.makedirs(run_root)
    ps, acc = [], {}
    # fake hash.k2d headers (capacity, size, key_bits, value_bits): size = 1000 * (order index + 1)
    for oi, order in enumerate(ORDERS):
        for sd in (1, 2, 3):
            d = os.path.join(run_root, order, str(sd), "db")
            if not os.path.isdir(d):
                os.makedirs(d)
            with open(os.path.join(d, "hash.k2d"), "wb") as f:
                f.write(struct.pack("<4Q", 5000, 1000 * (oi + 1) + sd, 18, 14))
                f.write(b"\0" * 64)
            with open(os.path.join(d, "opts.k2d"), "wb") as f:   # k, l, spaced-seed mask
                f.write(struct.pack("<QQQ", 35, 31, 0x3ffffffff3333333))
    for ds in datasets:
        ps += per_seed(lib, ds)
        kind = {"strain-exclusion": "labeled", "fda-argos": "argos", "cami": "cami"}.get(ds)
        if kind:
            acc.setdefault(kind, []).extend(acc_rows(lib, ds, kind))
    out = E.export_library(run_root, lib, ps, acc, DENS, QUERY)
    return run_root, out


class ExportLayout(unittest.TestCase):
    def setUp(self):
        self.base = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.base, ignore_errors=True)

    def fresh_base(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        return d

    def lib(self, name="strain_exclusion", datasets=("strain-exclusion", "fda-argos", "cami")):
        return make_library(self.base, name, list(datasets))

    def test_sheets_and_files(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        self.assertEqual(check_xlsx(xlsx), [])
        # profile is present here (fixture has timing data); paper_check only for strain_exclusion
        self.assertEqual(list(tables), list(E.SHEETS))
        # the export appends the density-effect sheets (benchmark/common/density_effect.py)
        self.assertEqual(list(load_xlsx(xlsx)), list(E.SHEETS) + ["density_effect", "density_effect_points", "significance"])
        for name in tables:
            self.assertTrue(os.path.isfile(os.path.join(tdir, "kraken2_%s.csv" % name)), name)
        import zipfile
        self.assertFalse([n for n in zipfile.ZipFile(xlsx).namelist() if "chart" in n or "drawing" in n])

    def test_no_old_names_or_repeat_columns(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        banned = ("_repeat_", "build_measured", "n_classify")
        for name, (fields, rows) in tables.items():
            for f in fields:
                self.assertFalse(any(b in f for b in banned), (name, f))
            if name in ("overview", "paper_check"):
                continue
            for r in rows:
                for k, v in r.items():
                    if k in ("order", "scheme"):
                        self.assertIn(v, SHORT)
        for r in tables["overview"][1]:
            text = " ".join(str(x) for x in r.values())
            for old in ("build_measured", "_repeat_std", "n_classify", "RandomXor", "ModRandom", "GreedyMini_w"):
                self.assertNotIn(old, text)

    def test_legend_matches_sheets(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        leg = tables["overview"][1]
        for sheet in ("build", "build_runs", "summary", "runs", "profile", "vs_reference"):
            if sheet not in tables:
                continue
            for f in tables[sheet][0]:
                base = re.sub(r"_(mean|std|ratio)$", "", f)
                if f in E.ID or base in ("library", "dataset", "minimum_hit_groups", "n_repeats", "run_id",
                                         "seed_idx", "seed", "gm_k", "projection", "gm_table", "mod_t"):
                    continue
                self.assertTrue(any(sheet in r["sheet"].split("/") and base in re.split(r" / |, ", r["column"])
                                    for r in leg), (sheet, f))
        # no legend row for a sheet that is not in the workbook
        for r in leg:
            for sh in r["sheet"].split("/"):
                if sh in E.SHEETS:
                    self.assertIn(sh, tables, r["column"])
        # no pre-spec names left in the legend
        text = " ".join(" ".join(map(str, r.values())) for r in leg)
        for old in ("n_seeds", "db_load_s", "classify_processing_s", "least-busy"):
            self.assertNotIn(old, text)

    def test_build_once_per_library_order(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        build_runs, build = tables["build_runs"][1], tables["build"][1]
        self.assertEqual(len(build_runs), 4 * 3)          # order x seed, not x dataset / mhg
        self.assertEqual(len(build), 4)
        self.assertNotIn("minimum_hit_groups", tables["build"][0])
        self.assertNotIn("minimum_hit_groups", tables["build_runs"][0])
        r0 = [r for r in build if r["order"] == "xor"][0]
        self.assertAlmostEqual(r0["build_wall_s_mean"], 102.0)
        self.assertAlmostEqual(r0["build_wall_s_std"], 1.0)
        self.assertEqual(r0["build_wall_s_ratio"], 1.0)
        pcg = [r for r in build if r["order"] == "xor-pcg-hash"][0]
        self.assertAlmostEqual(pcg["build_wall_s_ratio"], round(202.0 / 102.0, 4))
        self.assertAlmostEqual(r0["expected_density_factor_mean"], 1.5)
        self.assertEqual(sorted(r["order"] for r in build), sorted(SHORT))
        for sheet in ("summary", "runs"):
            for f in tables[sheet][0]:
                self.assertFalse(f.startswith(("build_", "index_size", "rank_")) and f != "build_pinned_cpus", (sheet, f))
        n_rec = len(step_record.FIELDS)   # step-record columns follow the CPU lists
        self.assertEqual(tables["runs"][0][-2 - n_rec:], ["build_pinned_cpus", "classify_pinned_cpus"]
                         + step_record.names("classify"))
        self.assertEqual(tables["build_runs"][0][-1 - n_rec:], ["build_pinned_cpus"] + step_record.names("build"))
        for sheet in ("build", "summary"):
            self.assertFalse([f for f in tables[sheet][0] if "pagefaults" in f or "read_gb" in f], sheet)

    def test_column_order_density_first(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        ids = {"build": 7, "build_runs": 9, "summary": 9, "runs": 11}
        for sheet in ("build", "build_runs", "summary", "runs"):
            fields = tables[sheet][0]
            i = ids[sheet]
            names = [f.replace("_mean", "") for f in fields[i:] if not f.endswith(("_std", "_ratio"))]
            self.assertEqual(names[:3], ["expected_density_factor", "density_factor", "distinct_minimizers"], sheet)
            if sheet.startswith("build"):
                self.assertEqual(names[3:6], ["index_size_gb", "build_wall_s", "build_max_rss_gb"])
            else:
                self.assertEqual(names[3], "classify_wall_s")
            cats = E._categories(sheet, fields)
            self.assertEqual(cats[fields.index(fields[i])], "density")
            self.assertIn("minimizers", cats)
        # runs carries the per-seed order properties joined by seed
        runs = tables["runs"][1]
        r = [x for x in runs if x["order"] == "greedymini" and x["seed_idx"] == 2][0]
        self.assertEqual(r["distinct_minimizers"], 3002)
        self.assertEqual(tables["orders"][0][-4:-2], ["expected_density_factor_mean", "expected_density_factor_std"])
        self.assertEqual(tables["orders"][0][:10], ["order", "k", "w", "L", "gm_k", "projection", "gm_table",
                                                     "mod_t", "description", "cli"])

    def test_distinct_minimizers_header(self):
        d = tempfile.mkdtemp()
        try:
            os.makedirs(os.path.join(d, "xor", "1", "db"))
            p = os.path.join(d, "xor", "1", "db", "hash.k2d")
            with open(p, "wb") as f:
                f.write(struct.pack("<4Q", 10, 7, 18, 14) + b"x" * 100)
            self.assertEqual(E._distinct_minimizers(d, "xor", 1), 7)
            self.assertIsNone(E._distinct_minimizers(d, "xor", 2))          # no db
            with open(p, "wb") as f:
                f.write(b"short")
            self.assertIsNone(E._distinct_minimizers(d, "xor", 1))           # truncated header
            self.assertIsNone(E._distinct_minimizers(None, "xor", 1))
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_distinct_minimizers_columns(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        build = dict((r["order"], r) for r in tables["build"][1])
        self.assertAlmostEqual(build["xor"]["distinct_minimizers_mean"], 1002.0)
        self.assertEqual(build["xor"]["distinct_minimizers_ratio"], 1.0)
        self.assertAlmostEqual(build["greedymini"]["distinct_minimizers_ratio"], round(3002 / 1002, 4))
        self.assertIn("distinct_minimizers_ratio", tables["summary"][0])
        spec = E.sheet_spec("build", tables["build"][0], tables["build"][1])
        self.assertEqual(spec["ratio_fields"]["distinct_minimizers_ratio"], "lower")

    def test_correlation_rows(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        import zipfile
        for sheet in ("build", "summary"):
            fields, rows = tables[sheet]
            footer, flat = E.correlation_rows(sheet, fields, rows)
            self.assertEqual([f["label"].split(" (")[0] for f in footer],
                             ["correlation with density", "correlation with distinct minimizers"])
            n = len(rows)
            dcol = E.col_letter(fields.index("density_factor_ratio") + 1)
            cells = dict((f, (fm, cv)) for f, fm, cv in footer[0]["cells"])
            self.assertEqual(set(cells), set(f for f in fields if f.endswith("_ratio")))
            fm, cv = cells["density_factor_ratio"]
            self.assertEqual(fm, 'IFERROR(CORREL($%s$3:$%s$%d,%s$3:%s$%d),"")' % (dcol, dcol, n + 2, dcol, dcol, n + 2))
            self.assertEqual(cv, 1.0)
            # cached value: Pearson r over the pairs where both cells are numbers
            col = "distinct_minimizers_ratio"
            pairs = [(r["density_factor_ratio"], r[col]) for r in rows
                     if r.get("density_factor_ratio") is not None and r.get(col) is not None]
            self.assertAlmostEqual(cells[col][1], round(E._pearson([a for a, _ in pairs], [b for _, b in pairs]), 4), places=4)
            for fr in flat:
                self.assertEqual(fr["sheet"], sheet)
        self.assertIn("(n = 4)", E.correlation_rows("build", *tables["build"])[0][0]["label"])
        # blanks are skipped pairwise; constant column -> no value
        fields = ["library", "density_factor_ratio", "x_ratio", "y_ratio"]
        rows = [{"density_factor_ratio": 1.0, "x_ratio": 1.0, "y_ratio": 2.0},
                {"density_factor_ratio": 2.0, "x_ratio": None, "y_ratio": 2.0},
                {"density_factor_ratio": 3.0, "x_ratio": 3.0, "y_ratio": 2.0},
                {"density_factor_ratio": 4.0, "x_ratio": 4.0, "y_ratio": 2.0}]
        _f, flat = E.correlation_rows("summary", fields, rows)
        v = dict((r["column"], r) for r in flat if r["row_label"] == "correlation with density")
        self.assertEqual(v["x_ratio"]["n"], 3)
        self.assertAlmostEqual(v["x_ratio"]["value"], 1.0)
        self.assertIsNone(v["y_ratio"]["value"])
        # workbook XML: formula + cached value; csv mirror has data only + correlations csv
        z = zipfile.ZipFile(xlsx)
        names = list(load_xlsx(xlsx))
        sheetxml = z.read("xl/worksheets/sheet%d.xml" % (names.index("build") + 1)).decode()
        self.assertIn("CORREL($", sheetxml)
        self.assertIn("<v>1.0</v>", sheetxml.split("CORREL(")[1])
        self.assertEqual(check_xlsx(xlsx), [])
        with open(os.path.join(tdir, "kraken2_build.csv")) as f:
            self.assertNotIn("correlation", f.read())
        crows = E.read_csv(os.path.join(tdir, "kraken2_correlations.csv"))[1]
        self.assertTrue(crows)
        self.assertEqual(set(r["sheet"] for r in crows), set(["build", "summary"]))
        self.assertEqual(set(E.read_csv(os.path.join(tdir, "kraken2_correlations.csv"))[0]),
                         set(E.CORREL_FIELDS))

    def test_classify_and_runs(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        runs, summ = tables["runs"][1], tables["summary"][1]
        self.assertEqual(len(runs), 3 * 4 * 2 * 3)        # dataset x order x mhg x seed
        self.assertEqual(len(summ), 3 * 4 * 2)
        f = tables["runs"][0]
        self.assertEqual(f[:9], ["run_id", "library", "dataset", "order", "gm_k",
                                 "projection", "gm_table", "mod_t", "minimum_hit_groups"])
        self.assertEqual(tables["summary"][0][:9], ["library", "dataset", "order", "gm_k",
                                                    "projection", "gm_table", "mod_t", "minimum_hit_groups", "n_repeats"])
        g = [r for r in runs if r["order"] == "greedymini"][0]
        self.assertEqual((g["gm_k"], g["projection"], g["gm_table"], g["mod_t"]),
                         (15, "at0-cg1", "w5_k15_canonical", None))
        self.assertEqual(g["run_id"].split("/")[1:3], ["greedymini", "k15_at0-cg1"])
        m = [r for r in runs if r["order"] == "mod-minimizer-xor"][0]
        self.assertEqual((m["gm_k"], m["projection"], m["mod_t"]), (None, None, 6))
        for c in ("classify_load_s", "classify_run_s", "classify_wall_s", "classify_max_rss_gb"):
            self.assertIn(c, f)
        self.assertLess(f.index("classify_wall_s"), f.index("classify_load_s"))
        r0 = [r for r in summ if r["dataset"] == "strain-exclusion" and r["order"] == "xor"
              and r["minimum_hit_groups"] == 0][0]
        self.assertAlmostEqual(r0["classify_wall_s_mean"], 2.2, places=6)    # seeds 1..3 -> 2.1 2.2 2.3
        self.assertAlmostEqual(r0["classify_wall_s_std"], 0.1, places=6)
        self.assertAlmostEqual(r0["classify_load_s_mean"], 1.1, places=6)
        self.assertAlmostEqual(r0["classify_run_s_mean"], 1.1, places=6)
        self.assertEqual(r0["classify_wall_s_ratio"], 1.0)
        self.assertEqual(r0["n_repeats"], 3)
        r2 = [r for r in summ if r["dataset"] == "strain-exclusion" and r["order"] == "xor"
              and r["minimum_hit_groups"] == 2][0]
        self.assertAlmostEqual(r2["classify_wall_s_mean"], 22.2, places=6)
        self.assertAlmostEqual(r2["classified_pct_mean"], 75.0)
        self.assertEqual(runs[0]["seed"][:2], "0x")
        self.assertFalse([c for c in f if c.startswith(("se_", "argos_", "cami_", "mock_"))])
        self.assertAlmostEqual(r0["genus_sensitivity_mean"], 0.702, places=6)
        cami = [r for r in summ if r["dataset"] == "cami" and r["order"] == "xor"][0]
        self.assertEqual(cami["genus_f1_mean"], 0.7)
        self.assertIsNone(cami.get("genus_concordant_fraction_mean"))
        argos = [r for r in summ if r["dataset"] == "fda-argos" and r["order"] == "xor-pcg-hash"][0]
        self.assertAlmostEqual(argos["genus_concordant_fraction_mean"], 0.87)

    def test_sorting_and_ratio_reference_per_dataset_mhg(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        summ = tables["summary"][1]
        keys = [(E.dataset_key(r["dataset"]), r["minimum_hit_groups"], E.order_key(r)) for r in summ]
        self.assertEqual(keys, sorted(keys))
        # the orders of one (dataset, mhg) are consecutive, xor first
        for i in range(0, len(summ), 4):
            chunk = summ[i:i + 4]
            self.assertEqual([r["order"] for r in chunk], SHORT)
            self.assertEqual(len(set((r["dataset"], r["minimum_hit_groups"]) for r in chunk)), 1)
        # ratio = order mean / xor mean at the same dataset + mhg (not across mhg)
        for r in summ:
            ref = [x for x in summ if x["order"] == "xor" and x["dataset"] == r["dataset"]
                   and x["minimum_hit_groups"] == r["minimum_hit_groups"]][0]
            self.assertAlmostEqual(r["classify_wall_s_ratio"],
                                   round(r["classify_wall_s_mean"] / ref["classify_wall_s_mean"], 4))
        pcg2 = [r for r in summ if r["dataset"] == "cami" and r["order"] == "xor-pcg-hash"
                and r["minimum_hit_groups"] == 2][0]
        self.assertAlmostEqual(pcg2["classify_wall_s_ratio"], round(23.2 / 22.2, 4))
        runs = tables["runs"][1]
        rk = [(E.dataset_key(r["dataset"]), r["minimum_hit_groups"], E.order_key(r), r["seed_idx"]) for r in runs]
        self.assertEqual(rk, sorted(rk))
        # sheet spec: borders per (dataset, mhg) block, dataset cell spans both mhg blocks
        spec = E.sheet_spec("summary", tables["summary"][0], summ)
        self.assertEqual(spec["block_fields"], ("library", "dataset", "minimum_hit_groups"))
        self.assertIn("minimum_hit_groups", spec["merge_fields"])
        self.assertEqual(spec["merge_scope"]["dataset"], ("library", "dataset"))

    def test_ratio_rules_and_bands(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields = tables["summary"][0]
        for neutral in ("classified_pct", "genus_ancestor_fraction", "genus_unclassified_fraction",
                        "genus_excluded_fraction", "species_excluded_fraction"):
            self.assertIn(neutral + "_mean", fields)
            self.assertNotIn(neutral + "_ratio", fields)
        for directional in ("genus_f1", "species_ppv", "genus_discordant_fraction",
                            "classify_mreads_per_min",
                            "classify_load_s", "classify_run_s"):
            self.assertIn(directional + "_ratio", fields)
        spec = E.sheet_spec("summary", fields, tables["summary"][1])
        self.assertEqual(spec["ratio_fields"]["genus_f1_ratio"], "higher")
        self.assertEqual(spec["ratio_fields"]["genus_discordant_fraction_ratio"], "lower")
        self.assertEqual(spec["ratio_fields"]["classify_load_s_ratio"], "lower")
        self.assertNotIn("classified_pct_ratio", spec["ratio_fields"])
        cats = E._categories("summary", fields)
        self.assertEqual(cats[fields.index("genus_f1_mean")], E.BAND_READ)
        self.assertEqual(cats[fields.index("genus_excluded_fraction_mean")], E.BAND_CAMI)
        self.assertEqual(cats[fields.index("genus_concordant_fraction_mean")], E.BAND_ARGOS)
        _rr, (_x, _t, t2) = make_library(self.fresh_base(), "standard_p10", ["speed"])
        self.assertFalse([f for f in t2["summary"][0] if f.startswith(("genus_", "species_"))])

    def test_ratio_blank_when_reference_missing_or_zero(self):
        def row(ds, mhg, order, v):
            return {"library": "L", "dataset": ds, "minimum_hit_groups": mhg, "order": order,
                    "classify_wall_s_mean": v}
        rows = [
            row("d", 0, "xor", 0.0), row("d", 0, "xor-pcg-hash", 5.0),
            row("e", 0, "xor-pcg-hash", 5.0),                       # no reference
            row("f", 0, "xor", 4.0), row("f", 0, "mod-minimizer-xor", None),
            row("g", 2, "xor-pcg-hash", 6.0),                       # reference exists only at mhg 0
            row("g", 0, "xor", 3.0),
        ]
        E.attach_ratios(rows, ["classify_wall_s", "classified_pct"])
        self.assertIsNone(rows[0]["classify_wall_s_ratio"])
        self.assertIsNone(rows[1]["classify_wall_s_ratio"])
        self.assertIsNone(rows[2]["classify_wall_s_ratio"])
        self.assertEqual(rows[3]["classify_wall_s_ratio"], 1.0)
        self.assertIsNone(rows[4]["classify_wall_s_ratio"])
        self.assertIsNone(rows[5]["classify_wall_s_ratio"])    # never against another mhg
        self.assertNotIn("classified_pct_ratio", rows[0])

    def test_profile_vs_reference_detail(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        prof = tables["profile"][1]
        self.assertEqual(len(prof), 3 * 4 * 2)
        for sheet in ("profile", "vs_reference"):
            self.assertIn("minimum_hit_groups", tables[sheet][0])
        pr = [r for r in prof if r["dataset"] == "strain-exclusion" and r["order"] == "mod-minimizer-xor"
              and r["minimum_hit_groups"] == 2][0]
        self.assertAlmostEqual(pr["scan_s_ratio"], round(0.33 / 0.3, 4))
        vr = tables["vs_reference"][1]
        self.assertEqual(len(vr), 3 * 3 * 2)
        self.assertIn("jaccard_pct_mean", tables["vs_reference"][0])
        self.assertEqual(list(tables)[-4:], ["accuracy_strain_exclusion", "accuracy_fda_argos",
                                             "accuracy_cami", "paper_check"])
        for sheet, dataset, has, hasnt in (
                ("accuracy_strain_exclusion", "strain-exclusion", "sensitivity_mean", "concordant_fraction_mean"),
                ("accuracy_fda_argos", "fda-argos", "concordant_fraction_mean", "sensitivity_mean"),
                ("accuracy_cami", "cami", "excluded_fraction_mean", "concordant_fraction_mean")):
            fields, rows = tables[sheet]
            self.assertEqual(set(r["dataset"] for r in rows), set([dataset]))
            self.assertEqual(set(r["minimum_hit_groups"] for r in rows), set(MHGS))
            self.assertEqual(fields[:2], ["library", "dataset"])
            self.assertLess(fields.index("order"), fields.index("minimum_hit_groups"))
            self.assertEqual(fields.index("n_repeats"), fields.index("minimum_hit_groups") + 1)
            self.assertIn(has, fields)
            self.assertNotIn(hasnt, fields)
            self.assertNotIn("seed_idx", fields)
            self.assertEqual(sorted(set(r["order"] for r in rows)), sorted(SHORT))
        se = tables["accuracy_strain_exclusion"][1]
        self.assertEqual(len(se), 4 * 2 * 3 * 2)           # order x mhg x group x rank
        self.assertEqual(se[0]["n_repeats"], 3)

    def test_accuracy_sheets_share_layout(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        sf, sr = tables["accuracy_strain_exclusion"]
        cf, cr = tables["accuracy_cami"]
        ident = ["library", "dataset", "group", "rank"] + list(order_names.ID_FIELDS) + ["minimum_hit_groups", "n_repeats"]
        self.assertEqual(sf[:len(ident)], ident)
        self.assertEqual(cf[:len(ident)], ident)
        for f in ("domain", "n_mean"):
            self.assertNotIn(f, sf)
        self.assertNotIn("domain", cf)
        cami_only = [f for f in cf if f.startswith(("ex_unclassified", "ex_consistent", "ex_off"))]
        self.assertTrue(cami_only)
        self.assertEqual([f for f in cf if f not in cami_only], sf)     # identical up to the cami-only columns
        base = [f.replace("_mean", "") for f in sf[len(ident):] if not f.endswith(("_std", "_ratio"))]
        self.assertEqual(base, ["n_truth", "n_calls", "n_scored", "tp", "vp", "fp", "fn", "no", "excluded",
                                "excluded_fraction", "sensitivity", "ppv", "f1"])
        # strain-exclusion: nothing excluded, n_scored = n_truth
        for r in sr:
            self.assertEqual(r["excluded_mean"], 0)
            self.assertEqual(r["excluded_fraction_mean"], 0)
            self.assertEqual(r["n_scored_mean"], r["n_truth_mean"])
        # same column groups (bands) for the same columns in both sheets
        sc, cc = E._categories("accuracy_strain_exclusion", sf), E._categories("accuracy_cami", cf)
        self.assertEqual(dict(zip(sf, sc)), dict((f, c) for f, c in zip(cf, cc) if f in sf))
        self.assertEqual(sc[sf.index("tp_mean")], E.BAND_COUNTS)
        self.assertEqual(cc[cf.index("ex_off_mean")], E.BAND_EXCLUDED_CALLS)
        # group order in both sheets: all, prokaryote, virus first, then bacteria, archaea
        self.assertEqual(list(OrderedDict.fromkeys(r["group"] for r in sr)), ["all", "prokaryote", "virus"])
        self.assertEqual(list(OrderedDict.fromkeys(r["group"] for r in cr)),
                         ["all", "prokaryote", "virus", "bacteria", "archaea"])
        self.assertEqual(len(cr), 5 * 4 * 2 * 2)          # group x order x mhg x rank

    def test_accuracy_aggregation(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields, rows = tables["accuracy_strain_exclusion"]
        for i in range(0, len(rows), 4):
            chunk = rows[i:i + 4]
            self.assertEqual(len(set((r["minimum_hit_groups"], r["group"], r["rank"]) for r in chunk)), 1)
            self.assertEqual([r["order"] for r in chunk], SHORT)
        # mhg is the outermost sort key after dataset
        mh = [r["minimum_hit_groups"] for r in rows]
        self.assertEqual(mh, sorted(mh))

        def pick(order, mhg):
            return [x for x in rows if x["order"] == order and x["group"] == "all" and x["rank"] == "genus"
                    and x["minimum_hit_groups"] == mhg][0]
        ref0, ref2 = pick("xor", 0), pick("xor", 2)
        self.assertAlmostEqual(ref0["sensitivity_mean"], 0.702, places=6)
        self.assertAlmostEqual(ref2["sensitivity_mean"], 0.502, places=6)
        r = pick("mod-minimizer-xor", 2)
        self.assertAlmostEqual(r["sensitivity_mean"], 0.472, places=6)
        self.assertAlmostEqual(r["sensitivity_std"], 0.001, places=6)
        self.assertAlmostEqual(r["sensitivity_ratio"], round(0.472 / 0.502, 4))   # against xor at mhg 2
        self.assertEqual(ref0["sensitivity_ratio"], 1.0)
        i = fields.index("sensitivity_mean")
        self.assertEqual(fields[i:i + 3], ["sensitivity_mean", "sensitivity_std", "sensitivity_ratio"])
        cf, _c = tables["accuracy_cami"]
        self.assertIn("f1_ratio", cf)
        self.assertNotIn("excluded_fraction_ratio", cf)
        af, _ar = tables["accuracy_fda_argos"]
        self.assertIn("discordant_fraction_ratio", af)
        self.assertNotIn("ancestor_fraction_ratio", af)
        spec = E.sheet_spec("accuracy_strain_exclusion", fields, rows)
        self.assertEqual(spec["ratio_fields"]["sensitivity_ratio"], "higher")
        self.assertIn("minimum_hit_groups", spec["block_fields"])
        self.assertIn("group", spec["block_fields"])
        self.assertNotIn("order", spec["merge_fields"])
        self.assertEqual(spec["merge_scope"]["minimum_hit_groups"], ("library", "dataset", "minimum_hit_groups"))

    def test_paper_check(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields, rows = tables["paper_check"]
        self.assertEqual(fields, E.PAPER_FIELDS)
        self.assertEqual(len(rows), 10)                    # 6 x S1 + 4 x S5
        by = dict(((r["table"], r["domain"], r["metric"]), r) for r in rows)
        r = by[("S1", "prokaryote", "sensitivity")]
        self.assertAlmostEqual(r["ours"], 0.702, places=6)      # xor, mhg 0 only (not mhg 2: 0.502)
        self.assertEqual(r["published"], 0.9396)
        self.assertAlmostEqual(r["difference"], round(0.702 - 0.9396, 6))
        self.assertEqual(r["within_tolerance"], "no")
        self.assertEqual(r["tolerance"], E.PAPER_TOLERANCE)
        self.assertEqual(r["n_repeats"], 3)
        r = by[("S5", None, "concordant")]
        self.assertAlmostEqual(r["ours"], 0.88)
        self.assertEqual(r["n_experiments"], 1)
        self.assertIn("experiments", r["note"])
        # published values come from the checked-in tables
        pub = dict(((p["table"], p["domain"], p["metric"]), p["value"]) for p in E.read_published())
        self.assertEqual(pub[("S1", "prokaryote", "f1")], 0.9657)
        self.assertEqual(pub[("S1", "virus", "f1")], 0.7536)
        self.assertAlmostEqual(sum(pub[("S5", None, m)] for m in
                                   ("concordant", "discordant", "ancestor", "unclassified")), 1.0, places=5)
        # within tolerance / not in repo / missing
        base = {"accuracy_strain_exclusion": [dict(library="strain_exclusion", order="xor", minimum_hit_groups=0,
                                                   rank="genus", group="prokaryote", n_repeats=5,
                                                   f1_mean=0.9650, f1_std=0.001)]}
        pubs = [dict(table="S1", row="Kraken 2", domain="prokaryote", rank="genus", metric="f1", value=0.9657),
                dict(table="S1", row="Kraken 2", domain="virus", rank="genus", metric="f1", value=None),
                dict(table="S1", row="Kraken 2", domain="prokaryote", rank="genus", metric="ppv", value=0.9)]
        out = E.paper_check_rows(base, pubs)
        self.assertEqual(out[0]["within_tolerance"], "yes")
        self.assertEqual(out[1]["within_tolerance"], "not in repo")
        self.assertIsNone(out[2]["within_tolerance"])
        self.assertIsNone(E.paper_check_rows({"accuracy_strain_exclusion": [dict(
            base["accuracy_strain_exclusion"][0], minimum_hit_groups=2)]}, pubs) or None)
        # no sheet without the xor / mhg 0 / strain_exclusion cell
        _rr, (_x, _t, t2) = make_library(self.fresh_base(), "standard_p10", ["strain-exclusion"])
        self.assertNotIn("paper_check", t2)

    def test_orders_density(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        fields, orders = tables["orders"]
        self.assertIn("expected_density_factor_mean", fields)
        self.assertIn("strain_exclusion_particular_density_factor_mean", fields)
        o = dict((r["order"], r) for r in orders)
        self.assertAlmostEqual(o["xor"]["expected_density_factor_mean"], 1.5)
        self.assertAlmostEqual(o["xor-pcg-hash"]["strain_exclusion_particular_density_factor_mean"], 1.7)
        self.assertEqual(o["mod-minimizer-xor"]["sketching"], "mod-minimizer")
        self.assertEqual(o["mod-minimizer-xor"]["rank_order"], "xor")
        self.assertEqual(o["xor-pcg-hash"]["sketching"], "minimizer")
        self.assertEqual(o["xor-pcg-hash"]["rank_order"], "xor-pcg-hash")
        self.assertEqual(o["greedymini"]["cli"].split()[:4], ["--sketching", "minimizer", "--order", "greedymini"])
        # unified names: k = minimizer length (old l), w = window, L = k + w - 1
        self.assertEqual(fields[10:14], ["sketching", "rank_order", "minimizer_spaces", "spaced_seed_mask"])
        self.assertNotIn("l", fields)
        for r in orders:
            self.assertEqual((r["k"], r["w"], r["L"]), (31, 5, 35))

    def test_orders_empty_cells_fail(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        rows = [dict(r) for r in tables["orders"][1]]
        E.check_orders(rows)
        rows[1]["expected_density_factor_mean"] = None
        with self.assertRaises(SystemExit) as cm:
            E.check_orders(rows)
        self.assertIn("xor-pcg-hash/expected_density_factor_mean", str(cm.exception))
        # an order without density rows makes the whole export fail
        nodens = [d for d in DENS if d["scheme"] != "xor-pcg-hash"]
        with self.assertRaises(SystemExit):
            E.export_library(run_root, "L", per_seed("L", "speed"), {}, nodens, QUERY)

    def test_merge_libraries_single_copy(self):
        make_library(self.base, "standard_p10", ["speed"])
        self.lib("strain_exclusion", ("strain-exclusion", "cami"))
        os.makedirs(os.path.join(self.base, "db", "standard", "library"))
        results = os.path.join(self.base, "results")
        # a stale csv of a sheet that no longer exists is removed
        with open(os.path.join(results, "csvs", "kraken2_cami_summary.csv"), "w") as f:
            f.write("x\n")
        # re-exporting a library replaces its rows (no duplicates)
        self.lib("strain_exclusion", ("strain-exclusion", "cami"))
        xlsx, rdir, n, _t = E.export_results(self.base)
        self.assertEqual(n, 2)
        self.assertEqual(xlsx, os.path.join(results, "kraken2_results.xlsx"))
        self.assertEqual(rdir, os.path.join(results, "csvs"))
        self.assertFalse(os.path.exists(os.path.join(rdir, "kraken2_cami_summary.csv")))
        self.assertEqual(check_xlsx(xlsx), [])
        # csvs / workbook exist only under results/csvs and results/
        for root, _d, files in os.walk(results):
            for f in files:
                if f.endswith((".csv", ".xlsx")):
                    self.assertIn(os.path.basename(root), ("csvs", os.path.basename(results)), os.path.join(root, f))
        self.assertTrue(os.path.isfile(xlsx))
        _f, summ = E.read_csv(os.path.join(rdir, "kraken2_summary.csv"))
        libs = [r["library"] for r in summ]
        self.assertEqual(libs, sorted(libs, key=E.library_key))
        self.assertEqual(set(libs), set(["standard_p10", "strain_exclusion"]))
        self.assertTrue(all(r["minimum_hit_groups"] in (0, 2) for r in summ))
        _f, lrows = E.read_csv(os.path.join(rdir, "kraken2_libraries.csv"))
        self.assertEqual(sorted(r["library"] for r in lrows), ["standard_p10", "strain_exclusion"])
        for r in summ:
            if r["order"] == "xor":
                self.assertEqual(r["classify_wall_s_ratio"], 1.0)
        _f, build = E.read_csv(os.path.join(rdir, "kraken2_build.csv"))
        self.assertEqual(len(build), 2 * 4)
        _f, runs = E.read_csv(os.path.join(rdir, "kraken2_runs.csv"))
        self.assertEqual(len(runs), 3 * 4 * 2 * 3)
        _f, pc = E.read_csv(os.path.join(rdir, "kraken2_paper_check.csv"))
        self.assertEqual(len(pc), 10)
        # re-rendering without a new library gives the same tables
        _x, _r, n2, t2 = E.export_results(self.base)
        self.assertEqual(n2, 2)
        self.assertEqual(len(t2["runs"][1]), len(runs))
        # a base dir with no csvs (or csvs without minimum_hit_groups) returns None
        empty = tempfile.mkdtemp()
        try:
            self.assertIsNone(E.export_results(empty))
            os.makedirs(os.path.join(empty, "results", "csvs"))
            with open(os.path.join(empty, "results", "csvs", "kraken2_runs.csv"), "w") as f:
                f.write("library,dataset,order,classify_wall_s\nstandard,speed,xor,1\n")
            self.assertIsNone(E.export_results(empty))
        finally:
            shutil.rmtree(empty, ignore_errors=True)

    def test_no_per_library_copy(self):
        run_root, _out = self.lib()
        self.assertEqual(sorted(n for n in os.listdir(run_root) if not os.path.isdir(os.path.join(run_root, n))), [])
        self.assertFalse(os.path.exists(os.path.join(run_root, "tables")))

    def test_resolve_layout(self):
        lay = E.resolve_layout("/x/kraken2/results/runs/strain_exclusion")
        self.assertEqual((lay["root"], lay["data"], lay["csvs"], lay["xlsx"]),
                         ("/x/kraken2", "/x/kraken2/data", "/x/kraken2/results/csvs", "/x/kraken2/results"))
        self.assertIsNone(E.resolve_layout("/x/somewhere/else"))

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

    def test_order_identity_helpers(self):
        # the order is the canonical name; the params are separate fields (never in the name)
        r = E.id_from_context("lib/greedymini/k15_at1-cg0/3 [cami, mhg2]", {("greedymini", 15, "at1-cg0", None): "w5_k15_canonical"})
        self.assertEqual(r, {"order": "greedymini", "gm_k": 15, "projection": "at1-cg0",
                             "gm_table": "w5_k15_canonical", "mod_t": None})
        self.assertEqual(E.order_dir(r), "greedymini/k15_at1-cg0")
        r = E.id_from_context("lib/mod-minimizer-greedymini/k6_at0-cg1_t6/1")
        self.assertEqual((r["order"], r["gm_k"], r["projection"], r["mod_t"]), ("mod-minimizer-greedymini", 6, "at0-cg1", 6))
        self.assertEqual(E.id_from_context("lib/xor-pcg-hash/2")["order"], "xor-pcg-hash")
        self.assertEqual(E.split_context("lib/xor-pcg-hash/3 [cami, mhg2]"), ("lib", "xor-pcg-hash", "3"))
        self.assertEqual(E.split_context("lib/greedymini/k15_at0-cg1/3"), ("lib", "greedymini/k15_at0-cg1", "3"))
        with self.assertRaises(ValueError):
            E.id_from_context("lib/greedymini_w5_k15_canonical/1")
        # sort: canonical order, then gm_k / projection
        a = {"order": "greedymini", "gm_k": 15, "projection": "at0-cg1"}
        b = {"order": "greedymini", "gm_k": 15, "projection": "at1-cg0"}
        self.assertLess(E.order_key(a), E.order_key(b))
        self.assertLess(E.order_key("xor"), E.order_key("xor-pcg-hash"))
        self.assertLess(E.order_key("xor-pcg-hash"), E.order_key(a))

    def test_two_projections_are_separate_rows(self):
        rows = []
        for proj in ("at0-cg1", "at1-cg0"):
            for sd in (1, 2):
                rows.append({"scheme": "greedymini", "dataset": "speed", "minimum_hit_groups": 0, "seed_idx": sd,
                             "context": "L/greedymini/k15_%s/%d" % (proj, sd), "seed": "1" * 64, "gm_table": "w5_k15_canonical",
                             "build_s": 10.0 + sd + (5 if proj == "at1-cg0" else 0), "class_s": 2.0})
        rows.append({"scheme": "xor", "dataset": "speed", "minimum_hit_groups": 0, "seed_idx": 1,
                     "context": "L/xor/1", "seed": "1" * 64, "build_s": 10.0, "class_s": 2.0})
        raw = E.build_tables("L", rows, {}, [])
        self.assertEqual(len(raw["build_runs"]), 5)
        build = E.aggregate(raw["build_runs"], E.BUILD_NAMES, E.AGG_KEYS["build"])
        self.assertEqual([(r["order"], r["projection"]) for r in build],
                         [("xor", None), ("greedymini", "at0-cg1"), ("greedymini", "at1-cg0")])
        self.assertAlmostEqual(build[1]["build_wall_s_mean"], 11.5)
        self.assertAlmostEqual(build[2]["build_wall_s_mean"], 16.5)
        spec = E.sheet_spec("build_runs", E.build_runs_field_order(), raw["build_runs"])
        for f in order_names.ID_FIELDS:
            self.assertIn(f, spec["merge_fields"])
        self.assertEqual(E.order_dir(raw["build_runs"][1]), "greedymini/k15_at0-cg1")

    def test_orders_sheet_and_joins(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        orders = tables["orders"][1]
        self.assertEqual([o["order"] for o in orders], SHORT)
        self.assertEqual(orders[2]["gm_table"], "w5_k15_canonical")
        self.assertEqual((orders[2]["gm_k"], orders[2]["projection"], orders[3]["mod_t"]), (15, "at0-cg1", 6))
        self.assertIn("w5_k15_canonical", orders[2]["cli"])
        self.assertIn("--gm-projection at0-cg1", orders[2]["cli"])
        raw = E.build_tables("L", per_seed("L", "speed"), {}, DENS)
        b = [r for r in raw["build_runs"] if r["order"] == "greedymini" and r["seed_idx"] == 1][0]
        self.assertAlmostEqual(b["expected_density_factor"], 1.7)
        # a density row of another projection does not match
        other = [dict(d, projection="at1-cg0") if d["scheme"] == "greedymini" else d for d in DENS]
        raw = E.build_tables("L", per_seed("L", "speed"), {}, other)
        b = [r for r in raw["build_runs"] if r["order"] == "greedymini" and r["seed_idx"] == 1][0]
        self.assertIsNone(b["expected_density_factor"])
        # xor-pcg-hash density rows never leak into xor
        raw = E.build_tables("L", per_seed("L", "speed"), {}, DENS)
        bx = [r for r in raw["build_runs"] if r["order"] == "xor" and r["seed_idx"] == 1][0]
        self.assertAlmostEqual(bx["expected_density_factor"], 1.5)

    def test_no_mock_dataset(self):
        run_root, (xlsx, tdir, tables) = self.lib()
        self.assertNotIn("accuracy_mock", tables)
        self.assertNotIn("accuracy_mock", E.SHEETS)
        self.assertNotIn("mock", E.DATASET_ORDER)
        import zipfile
        for name, (fields, rows) in tables.items():
            self.assertFalse([f for f in fields if "mock" in f or "bray" in f or "taxa_" in f], (name, fields))
            for r in rows:
                for v in r.values():
                    self.assertFalse(isinstance(v, str) and "mock" in v.lower() and name != "overview", (name, v))
        ov = " ".join(str(v) for r in tables["overview"][1] for v in r.values()).lower()
        self.assertNotIn("mock", ov)
        self.assertNotIn("zymo", ov)
        self.assertFalse([n for n in os.listdir(tdir) if "mock" in n])
        self.assertNotIn(b"mock", zipfile.ZipFile(xlsx).read("xl/workbook.xml"))

    def test_round_std(self):
        self.assertEqual(E.round_std(1.23456e-5), 1.235e-05)
        self.assertEqual(E.round_std(0.0), 0.0)
        self.assertEqual(E.round_std(123.456789), 123.4568)


if __name__ == "__main__":
    unittest.main()
