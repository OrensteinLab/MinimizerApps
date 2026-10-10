import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import paths  # noqa: E402


class T(unittest.TestCase):
    def test_layout(self):
        p = paths.bench_paths("kmc", {"KMC3_BENCH_ROOT": "/r"})
        self.assertEqual(str(p["data"]), "/r/data")
        self.assertEqual(str(p["runs"]), "/r/results/runs")
        self.assertEqual(str(p["csvs"]), "/r/results/csvs")
        self.assertEqual(str(p["xlsx"]), "/r/results")

    def test_override_and_missing(self):
        p = paths.bench_paths("sshash", {"SSHASH_BENCH_ROOT": "/r", "SSHASH_DATA_DIR": "/d"})
        self.assertEqual(str(p["data"]), "/d")
        with self.assertRaises(SystemExit):
            paths.bench_paths("minimap2", {})


if __name__ == "__main__":
    unittest.main()
