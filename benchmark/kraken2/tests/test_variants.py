"""bench_expand_variants / bench_configure_scheme (benchmark_interface.sh): every
GreedyMini scheme runs once per projection, the others once."""

import os
import subprocess
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = os.path.dirname(HERE)


def sh(body):
    script = 'source "%s/benchmark_interface.sh"; CODEBASE=/x; bench_iface_init; %s' % (K2, body)
    env = dict(os.environ, KRAKEN2_BENCH_ROOT="/nonexistent")
    env.pop("GM_PROJECTION", None)
    r = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True, env=env)
    return r.returncode, r.stdout.split(), r.stderr


class Variants(unittest.TestCase):
    def test_default_is_both_projections(self):
        rc, out, _ = sh('bench_expand_variants; echo "${VARIANTS[@]}"')
        self.assertEqual(rc, 0)
        self.assertEqual(out, ["xor", "xor-pcg-hash", "greedymini@at0-cg1", "greedymini@at1-cg0",
                               "mod-minimizer-xor", "mod-minimizer-greedymini@at0-cg1",
                               "mod-minimizer-greedymini@at1-cg0"])

    def test_single_projection_and_subset(self):
        rc, out, _ = sh('GM_PROJECTION=at1-cg0; DEFAULT_SCHEMES=(xor greedymini); bench_expand_variants; echo "${VARIANTS[@]}"')
        self.assertEqual(out, ["xor", "greedymini@at1-cg0"])

    def test_invalid_projection(self):
        rc, _out, err = sh('GM_PROJECTION=at0-cg1,bad; bench_expand_variants')
        self.assertNotEqual(rc, 0)
        self.assertIn("bad", err)

    def test_run_dirs(self):
        rc, out, _ = sh('bench_expand_variants; for v in "${VARIANTS[@]}"; do bench_configure_scheme "$v"; '
                        'echo "$SCHEME_BASE|$SCHEME_PROJECTION"; done')
        self.assertEqual(rc, 0)
        self.assertEqual(out, ["xor|", "xor-pcg-hash|", "greedymini/k15_at0-cg1|at0-cg1",
                               "greedymini/k15_at1-cg0|at1-cg0", "mod-minimizer-xor/t6|",
                               "mod-minimizer-greedymini/k6_at0-cg1_t6|at0-cg1",
                               "mod-minimizer-greedymini/k6_at1-cg0_t6|at1-cg0"])


if __name__ == "__main__":
    unittest.main()
