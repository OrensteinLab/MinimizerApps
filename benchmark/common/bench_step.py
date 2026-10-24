#!/usr/bin/env python3
"""Run one timed benchmark step pinned to CPUs and memory of one NUMA node, and record its memory state.

    bench_step.py --cpus LIST [--node K] [--no-membind] [--record STEP.json] [--poll S] -- PROG ARGS...

Put it where `taskset -c LIST` was, inside /usr/bin/time (pin_cpus.sh's BENCH_PIN does this):
the child gets the CPU affinity LIST and, unless --no-membind, the memory policy MPOL_BIND to node K
(set_mempolicy(2) through libc, no privileges or numactl needed; affinity and policy are inherited
by the step's own children and survive exec). Its exit status is the program's.

While the program runs, the step's process tree is sampled 0.5, 1, 2, 4 and 8 s after the start (so that a
short step, e.g. a classification of a few seconds, is sampled after loading too) and then every S
seconds (default 10; each sample walks the page tables of big processes, so it stays rare; the next one is
scheduled from the end of the previous one) from
/proc/<pid>/smaps_rollup -- the peak of the summed Rss, with the AnonHugePages of that sample -- and once, after 10 s,
from /proc/<pid>/numa_maps for the share of its anonymous pages on another node. The record
(--record, JSON) holds: node, cpus, membind, epoch_start/epoch_end (to line up the two slots), the node's
free and available memory and free 2 MB blocks before the step, peak_rss_gb, peak_anon_huge_gb,
huge_fraction (= AnonHugePages / Rss at the peak-Rss sample: ~1 when the step's memory came in huge pages, ~0
for 4 KB pages, which makes memory-latency-bound steps markedly slower), remote_frac, the mean CPU
clock of the pinned CPUs, and for steps of 1 s or more: foreign_busy_cores (other processes' busy
cores on the node's CPUs outside LIST), intruder_busy_cores (busy time on the pinned CPUs beyond
the step's own CPU time, i.e. other users' threads the pinning cannot keep out), and disk_busy_frac /
disk_read_gb / disk_write_gb of the device holding the results (device-wide: other users and the
other slot included). The wrapper adds ~50 ms of start-up to the step's wall time, and its few
samples' CPU time (milliseconds) to the CPU time /usr/bin/time reports. Python >= 3.6, standard library only.
"""
import argparse
import ctypes
import json
import os
import re
import resource
import signal
import subprocess
import sys
import time

import pin_cpus

MPOL_BIND = 2
SYS_SET_MEMPOLICY = {"x86_64": 238, "aarch64": 237}


FIRST_POLL = 0.5   # first memory sample (s); the interval then doubles up to --poll (about 4 extra samples per step)

def bind_memory(node):
    """set_mempolicy(MPOL_BIND, {node}) for the calling process (raises OSError on failure)."""
    nr = SYS_SET_MEMPOLICY.get(os.uname().machine)
    if nr is None:
        raise OSError("set_mempolicy: unsupported architecture %s" % os.uname().machine)
    libc = ctypes.CDLL(None, use_errno=True)
    mask = (ctypes.c_ulong * 4)()
    mask[node // 64] = 1 << (node % 64)
    if libc.syscall(nr, MPOL_BIND, mask, ctypes.c_ulong(256 + 1)) != 0:
        raise OSError(ctypes.get_errno(), "set_mempolicy")


def tree(pid):
    """pid and all its descendants (from /proc/<pid>/task/*/children)."""
    out, todo = [], [pid]
    while todo:
        p = todo.pop()
        out.append(p)
        try:
            for t in os.listdir("/proc/%d/task" % p):
                with open("/proc/%d/task/%s/children" % (p, t)) as f:
                    todo.extend(int(x) for x in f.read().split())
        except OSError:
            pass
    return out


def smaps(pid):
    """(Rss kB, AnonHugePages kB) of one process, (0, 0) if it is gone."""
    rss = huge = 0
    try:
        with open("/proc/%d/smaps_rollup" % pid) as f:
            for line in f:
                if line.startswith("Rss:"):
                    rss = int(line.split()[1])
                elif line.startswith("AnonHugePages:"):
                    huge = int(line.split()[1])
    except OSError:
        pass
    return rss, huge


def numa_pages(pid, node):
    """(anonymous pages on node, on other nodes) of one process from numa_maps (file pages, e.g.
    shared libraries already in another node's page cache, are not the step's allocations)."""
    local = remote = 0
    try:
        with open("/proc/%d/numa_maps" % pid) as f:
            for line in f:
                if " anon=" not in line:
                    continue
                for n, c in re.findall(r"\bN(\d+)=(\d+)", line):
                    if int(n) == node:
                        local += int(c)
                    else:
                        remote += int(c)
    except OSError:
        pass
    return local, remote


def block_device(path):
    """Kernel name (e.g. dm-1) of the block device holding path, or None."""
    try:
        st = os.stat(path)
        link = os.path.realpath("/sys/dev/block/%d:%d" % (os.major(st.st_dev), os.minor(st.st_dev)))
        return os.path.basename(link) if os.path.isdir(link) else None
    except OSError:
        return None


def diskstats(dev):
    """(sectors read, sectors written, io_ticks ms) of a block device, or None."""
    if not dev:
        return None
    try:
        with open("/proc/diskstats") as f:
            for line in f:
                v = line.split()
                if len(v) > 12 and v[2] == dev:
                    return int(v[5]), int(v[9]), int(v[12])
    except OSError:
        pass
    return None


def time_weighted(samples, t_end):
    """Mean of (time, value) samples, each weighted by the time until the next one (the last: until t_end)."""
    if not samples:
        return None
    ts = [t for t, _ in samples[1:]] + [max(t_end, samples[-1][0])]
    w = [max(b - a, 1e-3) for (a, _), b in zip(samples, ts)]
    return round(sum(wi * v for wi, (_, v) in zip(w, samples)) / sum(w))


def cpu_mhz(cpus):
    """Mean current clock (MHz) of `cpus` from /proc/cpuinfo, or None."""
    vals, cur = [], None
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("processor"):
                    cur = int(line.split(":")[1])
                elif line.startswith("cpu MHz") and cur in cpus:
                    vals.append(float(line.split(":")[1]))
    except (OSError, ValueError):
        pass
    return sum(vals) / len(vals) if vals else None


def busy_seconds(t_before, t_after, cpus):
    """Busy CPU seconds over `cpus` between two cpu_times() samples (USER_HZ ticks / 100)."""
    tot = 0.0
    hz = float(os.sysconf("SC_CLK_TCK"))
    for c in cpus:
        if c in t_before and c in t_after:
            tot += ((t_after[c][0] - t_before[c][0]) - (t_after[c][1] - t_before[c][1])) / hz
    return tot


def write_record(path, rec):
    if not path:
        return
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    with open(path, "w") as f:
        json.dump(rec, f, sort_keys=True)
        f.write("\n")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--" not in argv:
        print("usage: bench_step.py --cpus LIST [--node K] [--record F] -- PROG ARGS...", file=sys.stderr)
        return 2
    i = argv.index("--")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cpus", required=True)
    ap.add_argument("--node", type=int, default=None)
    ap.add_argument("--no-membind", action="store_true")
    ap.add_argument("--record", default=None)
    ap.add_argument("--poll", type=float, default=10.0, help="seconds between memory samples (default 10)")
    ap.add_argument("--io-path", default=None, help="file system whose disk I/O is recorded (default: the record's)")
    a = ap.parse_args(argv[:i])
    cmd = argv[i + 1:]
    cpus = set(pin_cpus.parse_list(a.cpus))
    node = a.node
    if node is None:
        nodes = set(pin_cpus.node_of(cpus).values())
        node = nodes.pop() if len(nodes) == 1 else None
    membind = node is not None and not a.no_membind

    def child_setup():
        os.sched_setaffinity(0, cpus)
        if membind:
            bind_memory(node)

    rec = {"node": node, "cpus": a.cpus, "membind": membind, "epoch_start": round(time.time(), 3)}
    if node is not None:
        rec["mem_free_gb_before"] = round(pin_cpus.node_meminfo(node).get("MemFree", 0) / 1048576.0, 2)
        rec["mem_avail_gb_before"] = round(pin_cpus.node_available_gb(node), 2)
        rec["free_2mb_gb_before"] = round(pin_cpus.free_2mb_gb().get(node, 0.0), 2)
    dev = block_device(a.io_path or (os.path.dirname(os.path.abspath(a.record)) if a.record else os.getcwd()))
    io_before = diskstats(dev)
    others = (pin_cpus.node_cpus(node) - cpus) if node is not None else set()
    t_before = pin_cpus.cpu_times()
    t0 = time.monotonic()
    try:
        proc = subprocess.Popen(cmd, preexec_fn=child_setup)
    except (OSError, subprocess.SubprocessError) as e:   # e.g. the memory binding failed
        rec.update({"rc": 127, "error": str(e)})
        write_record(a.record, rec)
        print("bench_step: cannot start %s: %s" % (cmd[0] if cmd else "?", e), file=sys.stderr)
        return 127
    # SIGINT reaches the whole process group already; forward SIGTERM (sent to this process only)
    signal.signal(signal.SIGTERM, lambda s, f: proc.send_signal(s))
    peak_rss = peak_huge = huge_at_peak = 0
    remote = None
    mhz = []
    next_mem = t0 + FIRST_POLL                # memory samples at 0.5, 1, 2, 4, 8 s, then every --poll s
    numa_at = t0 + max(10.0, a.poll)          # one numa_maps sample, once the program runs
    while True:
        now = time.monotonic()
        if now >= next_mem:
            pids = tree(proc.pid)
            rss = huge = 0
            for p in pids:
                r, h = smaps(p)
                rss += r
                huge += h
            if rss > peak_rss:
                peak_rss, huge_at_peak = rss, huge
            peak_huge = max(peak_huge, huge)
            m = cpu_mhz(cpus)
            if m:
                mhz.append((now, m))
            if node is not None and remote is None and now >= numa_at:
                lr = [numa_pages(p, node) for p in pids]
                l, r = sum(x[0] for x in lr), sum(x[1] for x in lr)
                remote = round(r / float(l + r), 4) if l + r else None
            done = time.monotonic()
            next_mem = done + min(a.poll, max(FIRST_POLL, now - t0))
        try:
            proc.wait(timeout=max(0.05, min(0.5, next_mem - time.monotonic())))
            break
        except subprocess.TimeoutExpired:
            pass
    wall = time.monotonic() - t0
    t_after = pin_cpus.cpu_times()
    io_after = diskstats(dev)
    own = resource.getrusage(resource.RUSAGE_CHILDREN)
    own_cpu = own.ru_utime + own.ru_stime
    rec.update({"rc": proc.returncode, "wall_s": round(wall, 3), "epoch_end": round(time.time(), 3),
                "peak_rss_gb": round(peak_rss / 1048576.0, 3), "peak_anon_huge_gb": round(peak_huge / 1048576.0, 3),
                "huge_fraction": round(huge_at_peak / float(peak_rss), 4) if peak_rss else None,
                "remote_frac": remote,
                "cpu_mhz": time_weighted(mhz, t0 + wall),
                "io_device": dev})
    if wall >= 1.0:   # load figures are meaningless over a fraction of a second
        rec["foreign_busy_cores"] = round(busy_seconds(t_before, t_after, others) / wall, 2) if others else None
        rec["intruder_busy_cores"] = round(max(0.0, busy_seconds(t_before, t_after, cpus) - own_cpu) / wall, 2)
        if io_before and io_after:
            rec["disk_busy_frac"] = round(min(1.0, (io_after[2] - io_before[2]) / 1000.0 / wall), 3)
            rec["disk_read_gb"] = round((io_after[0] - io_before[0]) * 512 / 1e9, 2)
            rec["disk_write_gb"] = round((io_after[1] - io_before[1]) * 512 / 1e9, 2)
    write_record(a.record, rec)
    return proc.returncode if proc.returncode >= 0 else 128 - proc.returncode


if __name__ == "__main__":
    sys.exit(main())
