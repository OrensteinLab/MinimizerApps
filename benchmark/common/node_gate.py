#!/usr/bin/env python3
"""Wait until one NUMA node is ready for a timed step (used by wait_low_load.sh when BENCH_NODE is set).

    node_gate.py --node K [--cores C] [--max-busy B] [--mem-gb M] [--huge-gb H]
                 [--load-wait S] [--mem-wait S] [--huge-wait S] [--window W]

Ready =
  * at least C quiet physical cores on the node (both hardware threads under 25% busy: the cores the
    step will be pinned to, one thread per core), and fewer than B busy cores over all the node's
    CPUs (other people's load: the call is made while the slot is idle, and the other slot runs on
    the other node);
  * at least M GB available on the node: free memory plus reclaimable page cache and slab (the
    node's own page cache fills up with a long run's input reads and is reclaimed on demand);
  * for memory-latency-bound steps, at least H GB free in blocks of 2 MB or more, the memory
    transparent huge pages come from. On this host THP defrag is "madvise", so a program that does
    not madvise never compacts memory itself; when the node is short, the gate compacts it: a
    helper bound to the node maps H GB with MADV_HUGEPAGE, touches and unmaps it (page cache is
    reclaimed, memory compacted, no privileges needed), then the condition is checked again.
Each condition is waited for at most its own time: the CPU one --load-wait seconds from the start
(default 3600), the memory one --mem-wait seconds of the time the CPU one held, and the huge-page one
--huge-wait seconds of the time both held (default 1800 each), so a long or intermittent CPU wait does
not use up the others; the up to 3 compactions are spread over the huge-page wait (at its start, after
a third and after two thirds of it). After that the step runs anyway and the line says so (the step record keeps the load and memory state of
every run, so such steps can be told apart).
Prints one "[load] ..." line to stderr. Python >= 3.6, standard library only.
"""
import argparse
import ctypes
import os
import subprocess
import sys
import time

import pin_cpus

MADV_HUGEPAGE = 14
PROT_RW = 0x1 | 0x2
MAP_PRIVATE_ANON = 0x02 | 0x20


def busy_cores(cpus, window):
    """Busy cores among `cpus` over `window` seconds (sum of per-CPU busy fractions)."""
    return sum(pin_cpus.busy_fractions(cpus, window).values())


QUIET = 0.25


def quiet_cores(busy, core):
    """Physical cores whose every hardware thread is under QUIET busy."""
    worst = {}
    for c, b in busy.items():
        worst[core[c]] = max(worst.get(core[c], 0.0), b)
    return sum(1 for b in worst.values() if b < QUIET)


def state(node, window):
    cpus = pin_cpus.node_cpus(node) & pin_cpus.allowed_cpus()
    busy = pin_cpus.busy_fractions(cpus, window)
    return {"busy": sum(busy.values()), "quiet_cores": quiet_cores(busy, pin_cpus.core_of(cpus)),
            "avail_gb": pin_cpus.node_available_gb(node),
            "free_2mb_gb": pin_cpus.free_2mb_gb().get(node, 0.0)}


def cpu_quiet(st, max_busy, cores):
    return st["busy"] < max_busy and st["quiet_cores"] >= cores


def ready(st, max_busy, cores, mem_gb, huge_gb, load_waived, mem_waived, huge_waived):
    """(cpu+memory ready, huge pages ready) -- pure function of a state sample."""
    cpu_ok = load_waived or cpu_quiet(st, max_busy, cores)
    base = cpu_ok and (mem_waived or st["avail_gb"] >= mem_gb)
    return base, base and (huge_waived or st["free_2mb_gb"] >= huge_gb)


class Clocks:
    """The waits of the three conditions; the memory and huge-page ones count only the time during which the
    conditions before them held (a node whose load comes and goes does not use them up)."""
    MAX_COMPACTIONS = 3

    def __init__(self, t0, load_wait, mem_wait, huge_wait, huge_gb):
        self.t0, self.load_wait, self.mem_wait, self.huge_wait, self.huge_gb = t0, load_wait, mem_wait, huge_wait, huge_gb
        self.last = t0
        self.mem_t = self.huge_t = 0.0    # time the CPU condition / the CPU and memory ones held
        self.compactions = 0

    def step(self, now, st, max_busy, cores, mem_gb):
        """(base, ok, (load_waived, mem_waived, huge_waived)) for the sample st taken at `now`; the time since
        the previous sample counts for a wait when the conditions before it hold at this one."""
        dt, self.last = now - self.last, now
        load_waived = now - self.t0 >= self.load_wait
        if load_waived or cpu_quiet(st, max_busy, cores):
            self.mem_t += dt
        mem_waived = self.mem_t >= self.mem_wait
        base, _ = ready(st, max_busy, cores, mem_gb, 0, load_waived, mem_waived, True)
        if base:
            self.huge_t += dt
        huge_waived = self.huge_gb > 0 and self.huge_t >= self.huge_wait
        base, ok = ready(st, max_busy, cores, mem_gb, self.huge_gb, load_waived, mem_waived, huge_waived)
        return base, ok, (load_waived, mem_waived, huge_waived)

    def compaction_due(self, base, st):
        """A compaction at the start of the huge-page wait, after a third and after two thirds of it."""
        if not base or st["free_2mb_gb"] >= self.huge_gb or self.compactions >= self.MAX_COMPACTIONS:
            return False
        return self.huge_t >= self.compactions * self.huge_wait / self.MAX_COMPACTIONS


def compact_in_process(gb):
    """Map gb GB of private anonymous memory with MADV_HUGEPAGE, touch every 2 MB, unmap."""
    libc = ctypes.CDLL(None, use_errno=True)
    libc.mmap.restype = ctypes.c_void_p
    libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
    n = int(gb * (1 << 30))
    addr = libc.mmap(None, n, PROT_RW, MAP_PRIVATE_ANON, -1, 0)
    if addr in (None, ctypes.c_void_p(-1).value):
        raise OSError(ctypes.get_errno(), "mmap")
    libc.madvise(ctypes.c_void_p(addr), ctypes.c_size_t(n), MADV_HUGEPAGE)
    step = 2 << 20
    for off in range(0, n, step):
        ctypes.c_char.from_address(addr + off).value = b"\x01"
    libc.munmap(ctypes.c_void_p(addr), ctypes.c_size_t(n))


def compact(node, gb):
    """Run compact_in_process(gb) in a child bound to the node (CPUs and memory), untimed."""
    cpus = sorted(pin_cpus.node_cpus(node) & pin_cpus.allowed_cpus())
    here = os.path.dirname(os.path.abspath(__file__))
    cmd = [sys.executable, os.path.join(here, "bench_step.py"), "--cpus", ",".join(map(str, cpus)),
           "--node", str(node), "--", sys.executable, os.path.abspath(__file__), "--compact-gb", str(gb)]
    subprocess.call(cmd, stdout=subprocess.DEVNULL)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--node", type=int)
    ap.add_argument("--cores", type=int, default=0, help="quiet physical cores the step needs")
    ap.add_argument("--max-busy", type=float, default=8.0)
    ap.add_argument("--load-wait", type=float, default=3600.0)
    ap.add_argument("--mem-gb", type=float, default=0.0)
    ap.add_argument("--huge-gb", type=float, default=0.0)
    ap.add_argument("--mem-wait", type=float, default=1800.0)
    ap.add_argument("--huge-wait", type=float, default=1800.0)
    ap.add_argument("--window", type=float, default=5.0)
    ap.add_argument("--compact-gb", type=float, default=None, help=argparse.SUPPRESS)  # the helper itself
    a = ap.parse_args(argv)
    if a.compact_gb is not None:
        compact_in_process(a.compact_gb)
        return 0
    if a.node is None:
        ap.error("--node is required")
    t0 = time.time()
    clocks = Clocks(t0, a.load_wait, a.mem_wait, a.huge_wait, a.huge_gb)
    first = None
    while True:
        st = state(a.node, a.window)
        first = first or st
        now = time.time()
        base, ok, (load_waived, mem_waived, huge_waived) = clocks.step(now, st, a.max_busy, a.cores, a.mem_gb)
        if ok:
            break
        if clocks.compaction_due(base, st):
            compact(a.node, a.huge_gb)   # then sample again
            clocks.compactions += 1
    compactions = clocks.compactions
    notes = []
    if load_waived and not cpu_quiet(st, a.max_busy, a.cores):
        notes.append("load wait timed out")
    if mem_waived and st["avail_gb"] < a.mem_gb:
        notes.append("memory wait timed out")
    if huge_waived and st["free_2mb_gb"] < a.huge_gb:
        notes.append("huge-page wait timed out")
    if compactions:
        notes.append("%d compaction(s)" % compactions)
    print("[load] node %d: %.1f busy cores (first %.1f, max %g), %d quiet cores (need %d), %.0f GB available "
          "(need %g), %.0f GB in 2 MB blocks (need %g), waited %.0fs%s" % (
              a.node, st["busy"], first["busy"], a.max_busy, st["quiet_cores"], a.cores, st["avail_gb"], a.mem_gb,
                                                st["free_2mb_gb"], a.huge_gb, time.time() - t0,
                                                (" (" + ", ".join(notes) + ")") if notes else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
