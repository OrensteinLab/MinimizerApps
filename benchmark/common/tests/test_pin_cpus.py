import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pin_cpus import choose, fmt_list, parse_list  # noqa: E402


def topology(nodes=2, cores=4, smt=2):
    """cpu ids as on the dev host: node = id % nodes, siblings = id and id + nodes*cores."""
    ncpu = nodes * cores * smt
    allowed = set(range(ncpu))
    node = {c: c % nodes for c in allowed}
    core = {c: c % (nodes * cores) for c in allowed}
    return allowed, node, core


class PinCpusTest(unittest.TestCase):
    def test_lists(self):
        self.assertEqual(parse_list("0-2,5,7-8"), [0, 1, 2, 5, 7, 8])
        self.assertEqual(fmt_list([8, 0, 1, 2, 5, 7]), "0-2,5,7-8")

    def test_single_node_one_thread_per_core(self):
        allowed, node, core = topology()
        nd, cpus = choose(3, 1, allowed, node, core, {c: 0.0 for c in allowed})
        self.assertEqual(len({node[c] for c in cpus}), 1)           # one NUMA node
        self.assertEqual(len({core[c] for c in cpus}), len(cpus))   # no two siblings
        self.assertEqual(len(cpus), 4)

    def test_avoids_busy_node_and_busy_sibling(self):
        allowed, node, core = topology()
        busy = {c: 0.0 for c in allowed}
        for c in allowed:
            if node[c] == 0:
                busy[c] = 0.9                                       # node 0 loaded
        busy[1] = 0.9                                               # one core of node 1 loaded ...
        nd, cpus = choose(3, 0, allowed, node, core, busy)
        self.assertEqual(nd, 1)
        self.assertNotIn(1, cpus)
        self.assertNotIn(9, cpus)                                   # ... through either hardware thread
        busy2 = {c: 0.0 for c in allowed}
        busy2[9] = 0.9                                              # sibling (9 = 1 + 8) loaded only
        _, cpus2 = choose(3, 0, allowed, node, core, busy2)
        self.assertNotIn(1, cpus2)

    def test_more_threads_than_cores_uses_second_hardware_threads_same_node(self):
        allowed, node, core = topology()
        nd, cpus = choose(6, 0, allowed, node, core, {c: 0.0 for c in allowed})
        self.assertEqual(len(cpus), 6)
        self.assertEqual(len({node[c] for c in cpus}), 1)           # node has 4 cores x 2 threads = 8 cpus

    def test_spills_to_other_node_only_when_needed(self):
        allowed, node, core = topology()
        nd, cpus = choose(10, 0, allowed, node, core, {c: 0.0 for c in allowed})
        self.assertEqual(len(cpus), 10)
        self.assertEqual(len({node[c] for c in cpus}), 2)

    def test_respects_allowed_set(self):
        allowed, node, core = topology()
        allowed = {c for c in allowed if node[c] == 1}
        nd, cpus = choose(2, 0, allowed, node, core, {c: 0.0 for c in allowed})
        self.assertEqual(nd, 1)
        self.assertTrue(set(cpus) <= allowed)


if __name__ == "__main__":
    unittest.main()
