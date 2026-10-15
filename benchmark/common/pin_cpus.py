#!/usr/bin/env python3
"""Choose the CPUs a timed benchmark step is pinned to (used by pin_cpus.sh).

    pin_cpus.py N [--spare S] [--window SECONDS]   ->  prints "<numa node> <cpu list>"

One NUMA node, one logical CPU per physical core:
  * topology from sysfs (node cpulists, thread_siblings_list), load from two /proc/stat samples;
    only CPUs in this process's Cpus_allowed_list are considered. Needs no privileges.
  * a core's busyness is the MAX over its siblings, so a core is "quiet" only if both hardware
    threads are, and the sibling of a pinned CPU is not running someone else's job;
  * N + S cores are taken (S spare cores for helper threads: readers, IO, main thread) from the
    single node whose best N + S cores are least busy (ties -> lowest node id); only one thread
    of each core is used. If the best node has fewer than N + S cores, its second hardware
    threads are added (same node); if that is still not enough, the remaining nodes are added.
Linux's default memory policy allocates on the node of the CPU that touches a page, so confining
the threads to one node also keeps the step's own memory local (no numactl needed).
Python >= 3.6, standard library only.
"""
import argparse
import os
import re
import sys
import time


def parse_list(text):
    out = []
    for part in text.strip().split(","):
        if not part:
            continue
        a, _, b = part.partition("-")
        out.extend(range(int(a), int(b or a) + 1))
    return out


def fmt_list(cpus):
    cpus = sorted(cpus)
    runs, start, prev = [], None, None
    for c in cpus:
        if start is None:
            start = prev = c
        elif c == prev + 1:
            prev = c
        else:
            runs.append((start, prev)); start = prev = c
    if start is not None:
        runs.append((start, prev))
    return ",".join(str(a) if a == b else "%d-%d" % (a, b) for a, b in runs)


def allowed_cpus():
    try:
        return set(os.sched_getaffinity(0))
    except AttributeError:
        text = open("/proc/self/status").read()
        return set(parse_list(re.search(r"Cpus_allowed_list:\s*(\S+)", text).group(1)))


def read(path, default=""):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return default


def node_of(allowed):
    """cpu -> numa node (0 for all if sysfs has no node info)."""
    node = {}
    base = "/sys/devices/system/node"
    for d in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        m = re.fullmatch(r"node(\d+)", d)
        if m:
            for c in parse_list(read("%s/%s/cpulist" % (base, d))):
                node[c] = int(m.group(1))
    return {c: node.get(c, 0) for c in allowed}


def core_of(allowed):
    """cpu -> core id = lowest CPU id of its hardware-thread sibling group (itself if unknown)."""
    core = {}
    for c in allowed:
        sib = parse_list(read("/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list" % c, str(c)))
        core[c] = min(sib) if sib else c
    return core


def cpu_times():
    t = {}
    for line in open("/proc/stat"):
        m = re.match(r"cpu(\d+)\s+(.*)", line)
        if m:
            v = [int(x) for x in m.group(2).split()]
            idle = v[3] + (v[4] if len(v) > 4 else 0)       # idle + iowait
            t[int(m.group(1))] = (sum(v[:8]), idle)         # exclude guest columns (already in user/nice)
    return t


def busy_fractions(allowed, window):
    a = cpu_times(); time.sleep(window); b = cpu_times()
    out = {}
    for c in allowed:
        if c in a and c in b:
            dt = b[c][0] - a[c][0]
            out[c] = (1.0 - (b[c][1] - a[c][1]) / dt) if dt > 0 else 0.0
        else:
            out[c] = 1.0                                    # offline / unknown: never pick it
    return out


def choose(n, spare, allowed, node, core, busy):
    """Return (node id, sorted cpu list) -- pure function of the topology and the load."""
    want = n + spare
    cores = {}                                              # (node, core) -> cpus sorted by id
    for c in sorted(allowed):
        cores.setdefault((node[c], core[c]), []).append(c)
    cbusy = {k: max(busy[c] for c in v) for k, v in cores.items()}
    nodes = sorted({k[0] for k in cores})

    def best(nd):
        ks = sorted((k for k in cores if k[0] == nd), key=lambda k: (cbusy[k], k[1]))
        return ks

    scored = []
    for nd in nodes:
        ks = best(nd)
        top = ks[:want]
        scored.append((sum(cbusy[k] for k in top) + (want - len(top)), nd))   # missing cores count as busy
    nd = min(scored)[1]
    ks = best(nd)
    picked = [cores[k][0] for k in ks[:want]]               # one thread per core
    if len(picked) < want:                                  # second hardware threads, same node
        for k in ks:
            for c in cores[k][1:]:
                if len(picked) < want:
                    picked.append(c)
    if len(picked) < want:                                  # still short: remaining nodes
        for other in nodes:
            if other != nd:
                for k in best(other):
                    for c in cores[k]:
                        if len(picked) < want:
                            picked.append(c)
    return nd, sorted(picked)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("n", type=int, help="threads of the step")
    ap.add_argument("--spare", type=int, default=2, help="extra cores for helper threads (default 2)")
    ap.add_argument("--window", type=float, default=1.0, help="load sampling window in seconds (default 1)")
    a = ap.parse_args(argv)
    allowed = allowed_cpus()
    nd, cpus = choose(a.n, a.spare, allowed, node_of(allowed), core_of(allowed), busy_fractions(allowed, a.window))
    print(nd, fmt_list(cpus))
    return 0


if __name__ == "__main__":
    sys.exit(main())
