import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import orders_sheet as os_  # noqa: E402


class T(unittest.TestCase):
    def test_fields(self):
        f = os_.orders_fields(["a", "b"])
        self.assertEqual(f[:4], ("order", "k", "w", "L"))
        self.assertNotIn("gm_lut", f)
        self.assertIn("gm_lut", os_.orders_fields(["a"], with_gm_lut=True))
        self.assertEqual(f[-2:], ("b_particular_density_factor_mean", "b_particular_density_factor_std"))

    def test_klw(self):
        r = os_.with_klw({}, 21, 11)
        self.assertEqual((r["k"], r["w"], r["L"]), (21, 11, 31))

    def test_missing(self):
        rows = [{"a": 1, "b": ""}, {"a": None, "b": 2}]
        self.assertEqual(os_.missing_cells(rows, ["a", "b"]), [(0, "b"), (1, "a")])
        self.assertEqual(os_.missing_cells(rows, ["a", "b"], optional=["b"]), [(1, "a")])


if __name__ == "__main__":
    unittest.main()
