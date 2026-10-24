"""Measurement-protocol helpers: free 2 MB blocks from buddyinfo, the node gate's decision, the step record's
columns, and bench_step.py running a program pinned, memory-bound and recorded."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import node_gate  # noqa: E402
import pin_cpus  # noqa: E402
import step_record  # noqa: E402

BUDDY = """Node 0, zone      DMA      1      1      1      0      2      1      1      0      1      1      3
Node 0, zone   Normal    100     50     20     10      5      2      1      1      1      4      8
Node 1, zone   Normal      0      0      0      0      0      0      0      0      0      0      0
"""


def test_free_2mb_gb_counts_blocks_of_order_9_and_up():
    got = pin_cpus.free_2mb_gb(BUDDY)
    # node 0: DMA 1 x 2 MB + 3 x 4 MB, Normal 4 x 2 MB + 8 x 4 MB = 54 MB
    assert got[0] == pytest.approx(54 / 1024.0)
    assert got[1] == 0


def test_gate_decision():
    st = {"busy": 2.0, "quiet_cores": 20, "avail_gb": 100.0, "free_2mb_gb": 5.0}
    r = node_gate.ready
    assert r(st, 4, 18, 50, 0, False, False, False) == (True, True)
    assert r(st, 4, 26, 50, 0, False, False, False) == (False, False)     # too few quiet cores
    assert r(st, 4, 26, 50, 0, True, False, False) == (True, True)        # ... until the load wait times out
    assert r(dict(st, busy=6.0), 4, 0, 50, 0, False, True, True) == (False, False)
    assert r(st, 4, 0, 50, 10, False, False, False) == (True, False)      # huge pages short
    assert r(st, 4, 0, 50, 10, False, False, True) == (True, True)        # ... until that wait is waived
    assert r(dict(st, avail_gb=10.0), 4, 0, 50, 0, False, False, False) == (False, False)
    assert r(dict(st, avail_gb=10.0), 4, 0, 50, 0, False, True, False) == (True, True)


def test_gate_clocks_count_only_the_time_the_conditions_before_held():
    busy = {"busy": 9.0, "quiet_cores": 0, "avail_gb": 100.0, "free_2mb_gb": 0.0}
    calm = dict(busy, busy=1.0, quiet_cores=30)
    low = dict(calm, avail_gb=10.0)
    c = node_gate.Clocks(0.0, load_wait=3600, mem_wait=1800, huge_wait=1800, huge_gb=10)
    step = lambda t, st: c.step(t, st, 8, 4, 50)
    # a node busy for 2500 s: neither the memory nor the huge-page wait has run
    assert step(2500, busy) == (False, False, (False, False, False))
    # calm but memory short for 1000 s: the memory wait runs, the huge-page one does not
    assert step(3500, low) == (False, False, (False, False, False))
    assert c.mem_t == 1000 and c.huge_t == 0
    # busy again (load wait not yet over): the memory wait stops, it is not used up by the load
    assert step(3590, busy)[2] == (False, False, False) and c.mem_t == 1000
    # memory available from here: the huge-page wait starts
    assert step(3600, calm) == (True, False, (True, False, False)) and c.huge_t == 10
    assert c.compaction_due(True, calm)                                     # first compaction at its start
    c.compactions += 1
    assert step(3700, calm)[1] is False and not c.compaction_due(True, calm)   # the next after a third
    assert step(4210, calm)[1] is False and c.compaction_due(True, calm)
    c.compactions += 2
    assert not c.compaction_due(True, calm)                                 # at most 3
    assert step(5389, calm)[1] is False
    assert step(5390, calm) == (True, True, (True, True, True))             # 1800 s of huge-page wait
    # an intermittent load: only the calm stretches count for the memory wait
    c = node_gate.Clocks(0.0, 10000, 1800, 1800, 0)                       # load wait not over in this part
    t = 0
    for _ in range(4):
        t += 500; c.step(t, low, 8, 4, 50)       # calm, memory short: +500
        t += 500; c.step(t, busy, 8, 4, 50)      # busy: +0
    assert c.mem_t == 2000 and t == 4000                                    # 4 x 500 s calm, 2000 s busy
    assert c.step(4001, low, 8, 4, 50)[:2] == (True, True)                  # waived (1800 s of calm reached)


def test_cpu_clock_is_weighted_by_time():
    import bench_step
    # three early samples 1 s apart at 1000 MHz, one at 3000 MHz covering the last 7 s
    assert bench_step.time_weighted([(0, 1000), (1, 1000), (2, 1000), (3, 3000)], 10) == 2400
    assert bench_step.time_weighted([], 10) is None


def test_quiet_cores_needs_both_threads_idle():
    busy = {0: 0.0, 64: 0.0, 2: 0.0, 66: 0.9, 4: 0.3, 68: 0.0}
    core = {0: 0, 64: 0, 2: 2, 66: 2, 4: 4, 68: 4}
    assert node_gate.quiet_cores(busy, core) == 1


def test_node_available_counts_reclaimable_cache():
    m = pin_cpus.node_meminfo(0)
    if not m:
        pytest.skip("no NUMA sysfs")
    assert pin_cpus.node_available_gb(0) >= m["MemFree"] / 1048576.0


def test_step_record_columns(tmp_path):
    assert step_record.columns(tmp_path / "missing.json", "build") == {c: None for c in step_record.names("build")}
    p = tmp_path / "step_build.json"
    p.write_text(json.dumps({"node": 1, "huge_fraction": 0.0, "remote_frac": None, "extra": 5}))
    cols = step_record.columns(p, "build")
    assert cols["build_node"] == 1 and cols["build_huge_fraction"] == 0.0 and "build_extra" not in cols
    assert set(cols) == set(step_record.names("build"))


@pytest.mark.skipif(not os.path.isdir("/sys/devices/system/node/node0"), reason="no NUMA sysfs")
def test_bench_step_pins_binds_and_records(tmp_path):
    cpu = min(pin_cpus.node_cpus(0) & pin_cpus.allowed_cpus())
    rec = tmp_path / "step.json"
    prog = ("import os; print(sorted(os.sched_getaffinity(0))); "
            "print(open('/proc/self/numa_maps').read().count('bind:0') > 0)")
    out = subprocess.run([sys.executable, str(HERE / "bench_step.py"), "--cpus", str(cpu), "--node", "0",
                          "--record", str(rec), "--", sys.executable, "-c", prog],
                         stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.split("\n")
    assert out[0] == "[%d]" % cpu and out[1] == "True"
    r = json.loads(rec.read_text())
    assert r["node"] == 0 and r["membind"] is True and r["rc"] == 0 and r["cpus"] == str(cpu)
    # load and disk figures are only recorded for steps of 1 s or more
    short_only = {"foreign_busy_cores", "intruder_busy_cores", "disk_busy_frac", "remote_frac"}
    assert set(step_record.FIELDS) - short_only <= set(r)
    # a short step is sampled after it has allocated: 0.5 s samples at the start
    rec2 = tmp_path / "short.json"
    prog = "import time; b = bytearray(100 << 20); time.sleep(0.8)"   # ends before the old first sample at 1 s
    subprocess.run([sys.executable, str(HERE / "bench_step.py"), "--cpus", str(cpu), "--record", str(rec2), "--",
                    sys.executable, "-c", prog], check=True)
    r2 = json.loads(rec2.read_text())
    assert r2["peak_rss_gb"] >= 0.1 and r2["huge_fraction"] is not None
    # the program's exit status is passed through
    assert subprocess.run([sys.executable, str(HERE / "bench_step.py"), "--cpus", str(cpu), "--",
                           "sh", "-c", "exit 3"]).returncode == 3
