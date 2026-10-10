import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from correlation import correlation_rows, pearson  # noqa: E402


class CorrelationTest(unittest.TestCase):
    def test_pearson(self):
        self.assertEqual(pearson([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertEqual(pearson([1, 2, 3], [3, 2, 1]), -1.0)
        self.assertIsNone(pearson([1, 2], [1, 2]))          # too few points
        self.assertIsNone(pearson([1, 1, 1], [1, 2, 3]))    # constant series

    def test_rows(self):
        fields = ["order", "d_ratio", "x_ratio", "y"]
        rows = [{"d_ratio": 1.0, "x_ratio": 1.0}, {"d_ratio": 2.0, "x_ratio": 2.0},
                {"d_ratio": 3.0, "x_ratio": 3.5}, {"d_ratio": None, "x_ratio": 9.0}]
        out = correlation_rows(fields, rows, (("correlation with density", "d_ratio"),
                                              ("absent", "nope_ratio")))
        self.assertEqual([r["label"] for r in out], ["correlation with density"])
        cells = {f: (fm, v) for f, fm, v in out[0]["cells"]}
        self.assertEqual(sorted(cells), ["d_ratio", "x_ratio"])          # only *_ratio columns
        self.assertIn("CORREL($B$3:$B$6,C$3:C$6)", cells["x_ratio"][0])
        self.assertAlmostEqual(cells["x_ratio"][1], 0.993, places=3)     # rows with a missing driver skipped
        self.assertEqual(correlation_rows(fields, [], (("a", "d_ratio"),)), [])


if __name__ == "__main__":
    unittest.main()
