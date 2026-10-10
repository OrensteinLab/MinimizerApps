import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import units  # noqa: E402


class T(unittest.TestCase):
    def test_conv(self):
        self.assertEqual(units.gb_from_bytes(2_500_000_000), 2.5)
        self.assertEqual(units.gb_from_kib(1_000_000), 1.024)
        self.assertIsNone(units.gb_from_bytes(None))
        self.assertIsNone(units.gb_from_kib(None))


if __name__ == "__main__":
    unittest.main()
