"""Scheme identity helpers shared by the shell harness (bench_const.py) and the
density table merge (density_csv.py): canonical orders only, params as fields.

    python3 -m unittest discover -s benchmark/kraken2/tests
"""

import contextlib
import csv
import io
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = os.path.dirname(HERE)
sys.path.insert(0, K2)
sys.path.insert(0, os.path.join(os.path.dirname(K2), "common"))

import bench_const as B  # noqa: E402
import density_csv as D  # noqa: E402
import order_names  # noqa: E402


class BenchConst(unittest.TestCase):
    def test_scheme_order_is_canonical(self):
        self.assertEqual(sorted(B.SCHEME_ORDER, key=B.SCHEME_ORDER.get), list(order_names.CANONICAL))
        self.assertEqual(B.scheme_rank("greedymini (k15_at1-cg0)"), B.SCHEME_ORDER["greedymini"])
        self.assertEqual(B.scheme_rank("greedymini [cami, mhg0]"), B.SCHEME_ORDER["greedymini"])
        self.assertEqual(B.scheme_rank("nonsense"), 99)

    def test_mod_t(self):
        self.assertEqual(B.mod_t(31, 35), 6)        # r + ((l - r) % w) = 4 + (27 % 5)

    def test_scheme_id_and_run_dir(self):
        i = B.scheme_id("greedymini", "/x/w5_k15_canonical.gm", "at1-cg0")
        self.assertEqual((i["order"], i["gm_k"], i["projection"], i["gm_table"], i["mod_t"]),
                         ("greedymini", 15, "at1-cg0", "w5_k15_canonical", None))
        self.assertEqual(B.scheme_run_dir(i), "greedymini/k15_at1-cg0")
        i = B.scheme_id("mod-minimizer-greedymini", "w5_k6_canonical.gm", None)
        self.assertEqual(B.scheme_run_dir(i), "mod-minimizer-greedymini/k6_at0-cg1_t6")
        self.assertEqual(B.scheme_run_dir(B.scheme_id("xor")), "xor")
        self.assertEqual(B.scheme_run_dir(B.scheme_id("mod-minimizer-xor")), "mod-minimizer-xor/t6")
        with self.assertRaises(ValueError):
            B.scheme_id("pcg")

    def test_cli_scheme_id(self):
        out = subprocess.check_output(
            [sys.executable, os.path.join(K2, "bench_const.py"), "scheme-id", "greedymini",
             "/x/w5_k15_canonical.gm", "at0-cg1", "31", "35"]).decode().strip()
        self.assertEqual(out, "greedymini/k15_at0-cg1|15|at0-cg1|w5_k15_canonical|")


class DensityCsv(unittest.TestCase):
    def test_merge_adds_identity_columns_and_keeps_projections_apart(self):
        d = tempfile.mkdtemp()
        out = os.path.join(d, "density.csv")
        head = "scheme,seed_idx,toggle_hex,l,k,w,fac_expected,fac_particular\n"
        for proj, val in (("at0-cg1", "1.4"), ("at1-cg0", "1.5")):
            new = os.path.join(d, "new.csv")
            with open(new, "w") as f:
                f.write(head + "greedymini,1,0x1,31,35,5,%s,%s\nxor,1,0x1,31,35,5,2.0,2.0\n" % (val, val))
            D.cmd_merge([out, new, "7", "31", "35", proj, "/x/w5_k15_canonical.gm", "/x/w5_k6_canonical.gm"])
        with open(out) as f:
            rows = list(csv.DictReader(f))
        gm = [r for r in rows if r["scheme"] == "greedymini"]
        self.assertEqual(sorted((r["projection"], r["fac_expected"]) for r in gm),
                         [("at0-cg1", "1.4"), ("at1-cg0", "1.5")])
        self.assertEqual(gm[0]["gm_k"], "15")
        self.assertEqual(gm[0]["gm_table"], "w5_k15_canonical")
        x = [r for r in rows if r["scheme"] == "xor"]
        self.assertEqual(len(x), 1)
        self.assertEqual(x[0]["gm_k"], "")
        # the tool's l / k header is rewritten with the unified names k, w, L
        self.assertEqual(list(rows[0].keys())[:6], ["scheme", "gm_k", "projection", "gm_table", "mod_t", "seed_idx"])
        self.assertNotIn("l", rows[0])
        self.assertEqual((x[0]["k"], x[0]["w"], x[0]["L"]), ("31", "5", "35"))
        fields = list(rows[0].keys())
        self.assertEqual(fields[fields.index("k"):fields.index("k") + 3], ["k", "w", "L"])


    def test_missing_recomputes_mod_rows_without_the_tie_rule(self):
        d = tempfile.mkdtemp()
        out = os.path.join(d, "density.csv")
        args = ["1", "0x1", "31", "35", "7", "at0-cg1", "/x/w5_k15_canonical.gm", "/x/w5_k6_canonical.gm"]
        new = os.path.join(d, "new.csv")
        with open(new, "w") as f:
            f.write("scheme,seed_idx,toggle_hex,l,k,w,fac_expected,fac_particular\n"
                    "xor,1,0x1,31,35,5,2.0,2.0\nmod-minimizer-xor,1,0x1,31,35,5,1.37,1.37\n")
        D.cmd_merge([out, new] + ["7", "31", "35", "at0-cg1", "/x/w5_k15_canonical.gm", "/x/w5_k6_canonical.gm"])

        def missing():
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                D.cmd_missing([out] + args + ["xor", "mod-minimizer-xor"])
            return buf.getvalue().split()

        self.assertEqual(missing(), [])
        with open(out) as f:
            rows = list(csv.DictReader(f))
        self.assertEqual({r["scheme"]: r["mod_tie_rule"] for r in rows},
                         {"xor": "", "mod-minimizer-xor": D.MOD_TIE_RULE})
        # a row of the old rule (no stamp) is missing; a plain order is not
        for r in rows:
            r["mod_tie_rule"] = ""
        with open(out, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        self.assertEqual(missing(), ["mod-minimizer-xor"])


if __name__ == "__main__":
    unittest.main()
