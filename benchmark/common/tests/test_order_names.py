import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import order_names as on  # noqa: E402


class OrderNames(unittest.TestCase):
    def test_validate(self):
        self.assertEqual(on.validate_order("xor"), "xor")
        with self.assertRaises(ValueError):
            on.validate_order("native")
        self.assertEqual(on.validate_order("native", on.BASELINES), "native")
        with self.assertRaises(ValueError):
            on.validate_order("xor-hash")

    def test_projection(self):
        self.assertEqual(on.projection_from_flag(True), "at1-cg0")
        self.assertEqual(on.projection_from_flag("at1"), "at1-cg0")
        self.assertEqual(on.projection_from_flag(None), "at0-cg1")
        with self.assertRaises(ValueError):
            on.projection_from_flag("zz")

    def test_table(self):
        p = on.parse_table_name("/x/orders/GreedyMini/w11_k15_canonical.gm")
        self.assertEqual((p["gm_table"], p["gm_k"], p["gm_w"]), ("w11_k15_canonical", 15, 11))
        self.assertEqual(on.parse_table_name("w20_k9_canonical_dna.gm")["gm_k"], 9)
        self.assertEqual(on.gm_fields("w5_k6_canonical.gm", "at1")["projection"], "at1-cg0")
        with self.assertRaises(ValueError):
            on.parse_table_name("junk.gm")

    def test_make_id_roundtrip(self):
        i = on.make_id("greedymini", gm_table="w11_k15_canonical.gm", projection="at1")
        self.assertEqual(i, {"order": "greedymini", "gm_k": 15, "projection": "at1-cg0",
                             "gm_table": "w11_k15_canonical", "mod_t": None})
        self.assertEqual(on.parse_id(on.format_id(i)), i)
        x = on.make_id("xor", gm_k=3, projection="at1", mod_t=4)
        self.assertEqual(tuple(x[f] for f in on.ID_FIELDS[1:]), (None,) * 4)
        m = on.make_id("mod-minimizer-greedymini", gm_k=6, mod_t="6")
        self.assertEqual((m["projection"], m["mod_t"]), ("at0-cg1", 6))

    def test_dir_name(self):
        self.assertEqual(on.order_dir_name("xor"), "xor")
        self.assertEqual(on.order_dir_name("greedymini", 15, "at1-cg0"), "greedymini/k15_at1-cg0")
        self.assertEqual(on.order_dir_name("mod-minimizer-xor", mod_t=10), "mod-minimizer-xor/t10")
        for args in [("greedymini", 13, "at0-cg1", None),
                     ("mod-minimizer-greedymini", 6, "at1-cg0", 6),
                     ("mod-minimizer-xor", None, None, 9), ("xor-pcg-hash", None, None, None)]:
            d = on.parse_order_dir(on.order_dir_name(*args))
            self.assertEqual(d["order"], args[0])
            self.assertEqual(d["gm_k"], args[1])
            self.assertEqual(d["mod_t"], args[3])


if __name__ == "__main__":
    unittest.main()
