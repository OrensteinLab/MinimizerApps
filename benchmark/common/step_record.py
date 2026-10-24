"""Columns of a timed step's record (bench_step.py --record), shared by the four exporters.

    columns(path, step) -> {"<step>_node": ..., "<step>_huge_fraction": ..., ...}

Every key is present (None when the record is missing, e.g. runs made before records existed), so
the run tables keep a fixed column set. DOC describes each column for the workbooks' doc sheets.
"""
import json
import os

FIELDS = ("node", "membind", "huge_fraction", "remote_frac", "foreign_busy_cores", "intruder_busy_cores",
          "disk_busy_frac", "cpu_mhz", "free_2mb_gb_before", "mem_avail_gb_before", "epoch_start", "epoch_end")

DOC = {
    "node": "NUMA node the step ran on (CPUs and, with membind, memory)",
    "membind": "whether the step's memory was bound to that node (set_mempolicy MPOL_BIND)",
    "huge_fraction": "AnonHugePages / Rss of the step's processes at their peak Rss: ~1 = its memory came in 2 MB "
                     "transparent huge pages, ~0 = 4 KB pages (memory-latency-bound steps are markedly slower)",
    "remote_frac": "share of the step's anonymous pages on another NUMA node, sampled once after 10 s (0 with "
                   "membind; blank for shorter steps)",
    "foreign_busy_cores": "other processes' busy cores on the node's CPUs outside the pinned set, averaged over "
                          "the step (contention the pinning cannot exclude: shared L3 and memory bandwidth)",
    "intruder_busy_cores": "busy cores on the pinned CPUs beyond the step's own CPU time, averaged over the step: "
                           "other users' threads the pinning cannot keep off those CPUs",
    "disk_busy_frac": "fraction of the step's wall time the results disk was busy (/proc/diskstats io_ticks; "
                      "device-wide, so it includes the other slot and other users)",
    "cpu_mhz": "mean clock of the pinned CPUs over the step (MHz, /proc/cpuinfo; the frequency is not controlled)",
    "free_2mb_gb_before": "GB free on the node in blocks of 2 MB or more before the step (what huge pages can use)",
    "mem_avail_gb_before": "GB available on the node before the step (free + reclaimable page cache and slab)",
    "epoch_start": "start of the step (Unix seconds), to line up the two slots' steps",
    "epoch_end": "end of the step (Unix seconds)",
}


def names(step):
    return ["%s_%s" % (step, f) for f in FIELDS]


def columns(path, step):
    rec = {}
    if path and os.path.isfile(str(path)):
        try:
            with open(str(path)) as f:
                rec = json.load(f)
        except (OSError, ValueError):
            rec = {}
    return {"%s_%s" % (step, f): rec.get(f) for f in FIELDS}
