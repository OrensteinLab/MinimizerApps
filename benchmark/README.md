# Benchmark harnesses: shared protocol and tools

Four apps are benchmarked across minimizer orders (random xor / xor-pcg-hash, GreedyMini, mod-minimizer
variants, and each app's own baseline): `minimap2/` (read mapping), `sshash/` (k-mer dictionary),
`kraken2/` (taxonomic classification) and `kmc/` (k-mer counting). Each folder has its own README (data,
`.env`, usage, outputs); this file describes what they share.

## Measurement protocol

Timings on a shared host are only comparable if every timed step is made under the same conditions. Pinning
the CPUs is not enough on its own: single-threaded, memory-latency-bound steps (sshash query and bench,
kraken2 classify) also depend on whether their memory comes in 2 MB transparent huge pages and stays on their
own NUMA node. Every timed step of every app therefore follows the same protocol
(`common/compare_result_sets.py` compares the run-to-run stability of two result sets, to check a protocol).

1. **Slots: one NUMA node per app** (`full_run_all.sh`, `BENCH_NODE`). The apps run in two slots, one per node,
   in parallel; the apps of a slot run one after the other. Every timed step of an app runs on its slot's node,
   so all orders and seeds of an app are measured on the same node.
   - The two slots share no cores, L3 cache, memory controller or memory.
   - They do share the disk (`/localdata`: kmc's temp files, kraken2's library reads, sshash's external sort)
     and the inter-socket link. Interleaving the orders (point 5) keeps this from favouring one order, and the
     step records keep the disk's busy share and each step's start and end, so overlaps can be checked afterwards.
   - Default assignment: node 0 runs minimap2 then kmc, node 1 runs kraken2 then sshash. Node 1 usually has more
     free 2 MB blocks, which the huge-page-sensitive steps need.
   - Without `BENCH_NODE` (an app run by hand), each step takes the node whose best cores are least busy.
2. **Node gate before each step** (`common/wait_low_load.sh` -> `common/node_gate.py`). The harness waits, while
   its slot is idle, until its node has:
   - **CPUs**: as many quiet physical cores as the step will be pinned to (both hardware threads under 25% busy),
     and fewer than `BENCH_NODE_LOAD_MAX` (default 8) busy cores over the node's CPUs. Measured from `/proc/stat`
     over `BENCH_LOAD_WINDOW` s (default 5), so it is other users' load: the other slot is on the other node.
   - **Memory**: the step's memory available (`wait_low_load MEM_GB`, set above each step's peak RSS). Available
     means free plus reclaimable page cache and slab: with the memory bound to the node, a long run fills the
     node's page cache with its own input reads.
   - **Huge pages**, for the memory-latency-bound steps (sshash bench and query, kraken2 classify): `HUGE_GB` free
     in blocks of 2 MB or more, the memory huge pages come from (`wait_low_load MEM_GB HUGE_GB`). On this host
     THP is `enabled=always` but `defrag=madvise`, so a program that does not madvise never compacts memory.
     When the node is short, the gate compacts it, untimed and without privileges: a helper bound to the node
     maps `HUGE_GB` with `MADV_HUGEPAGE`, touches it and unmaps it.

   Each wait is bounded, so a busy host slows the run down but never stalls it: CPUs `BENCH_LOAD_WAIT` (default
   3600 s), memory `BENCH_MEM_WAIT` and huge pages `BENCH_HUGE_WAIT` (default 1800 s each), each counted from
   when the conditions before it were first met, so a long CPU wait does not use up the memory and huge-page
   waits. The up to 3 compactions are spread over the huge-page wait (at its start, after a third and after
   two thirds of it). After a timeout the
   step runs anyway; the `[load] node N: ...` line in the run log says so, and the step record keeps the load
   and memory state, so such steps can be told apart. Without `BENCH_NODE`, a machine-wide gate (fewer than
   `BENCH_LOAD_MAX` busy cores) is used instead. `BENCH_WAIT_LOW_LOAD=0` disables the gate.
3. **CPU pinning, always** (`common/pin_cpus.sh`, selection in `common/pin_cpus.py`). The CPUs are chosen per
   step on the slot's node, needing no privileges:
   - **one logical CPU per physical core**: a core counts as busy if *either* hardware thread is, so a pinned
     thread never shares its core with its own sibling or with someone else's job;
   - **N threads + spare cores** for helper threads (readers, I/O, the main thread). Multi-threaded steps get
     N + 2 (`BENCH_PIN_SPARE`). The single-threaded sshash bench and query get 1 + 1, so the OS cannot move
     them across many cores; kraken2 classify (`--threads 1`, worker + reader) gets 2 + 2.

   Pinning is not optional: if the selection fails, the benchmark stops instead of timing unpinned.
   `BENCH_CPUS=<list>` forces an explicit list.
4. **Memory bound to the node** (`common/bench_step.py`, *inside* `/usr/bin/time -v`). It sets the CPU affinity
   and `set_mempolicy(MPOL_BIND)` to the slot's node, through libc: no `numactl`, no privileges. The policy and
   affinity are inherited by the step's children and kept across `exec`, so the program and its page-cache
   reads allocate on its own node only. The gate's memory check is what makes the binding safe.
   `BENCH_MEMBIND=0` falls back to the kernel's local-first policy.
   - The wrapper's Python start-up and its memory samples (at 0.5, 1, 2, 4 and 8 s, then one every 10 s) add a
     small amount to the step's wall and CPU time; a sample reads the page tables, so it costs about 6 ms per GB
     of 4 KB pages.
   - The kernel's automatic NUMA balancing skips a process bound with `MPOL_BIND`, so bound steps run without
     balancing scans.
5. **Repeats for every order, baseline included.** Each order runs for each seed in `config.yaml`; minimap2
   `native`, sshash `native` and kmc `default` are repeated like the others (their repeats only measure timing
   noise), so every comparison has replicates on both sides. kraken2 has no native run, so its reference order
   is `xor`. Orders are interleaved within a repeat, so drift hits all orders alike: kmc runs one `-O all -n N`
   call, sshash and kraken2 classify loop over the orders inside each seed, and the minimap2 `sweep` is
   seed-major (all orders for seed 0, then seed 1, ...).
6. **What is recorded for every step**:
   - from `/usr/bin/time -v`: wall seconds, CPU seconds (`*_cpu_s` = user + system, all threads), peak RSS
     (`*_max_rss_gb`), major/minor page faults and file-system reads (`*_major_pagefaults`, `*_read_gb`);
   - the load/run split where available (`map_load_s`/`map_run_s`, `classify_load_s`/`classify_run_s`,
     `query_load_s`/`query_run_s`);
   - the CPU list (`{step}_pinned_cpus`);
   - the **step record** (`step_*.json` in the run directory, from `bench_step.py`; columns
     `{step}_node`, `_membind`, `_huge_fraction`, `_remote_frac`, `_foreign_busy_cores`, `_intruder_busy_cores`,
     `_disk_busy_frac`, `_cpu_mhz`, `_free_2mb_gb_before`, `_mem_avail_gb_before`, `_epoch_start`, `_epoch_end`):
     - `huge_fraction`: AnonHugePages / Rss of the step's processes at their peak Rss (sampled at 0.5, 1, 2,
       4 and 8 s, then every 10 s, so a step of a few seconds is sampled after loading), ~1 = huge pages;
     - `remote_frac`: share of its anonymous pages on another node, 0 with the binding;
     - `foreign_busy_cores`: other processes' load on the node's unpinned CPUs during the step;
     - `intruder_busy_cores`: busy time on the pinned CPUs beyond the step's own CPU time, i.e. other users'
       threads that pinning cannot keep off them;
     - `disk_busy_frac`: share of the step the results disk was busy (device-wide);
     - `cpu_mhz`: mean clock of the pinned CPUs (the frequency is not controllable here);
     - `epoch_start` / `epoch_end`;
     - the node's available memory and free 2 MB blocks before the step.

     `<results>/host_settings.txt` keeps the kernel settings the timings depend on (THP, NUMA balancing, zone
     reclaim) and the commit.
7. **Headline metrics.**
   - Wall time of each step.
   - For kraken2 classify, `classify_run_s`: the read processing without the database load, which depends on
     the page cache much more than the processing does.
   - For sshash, `query_run_s` and the `bench` lookup times, read together with `query_huge_fraction` /
     `bench_huge_fraction`.
8. **A fresh results directory per protocol or code change** (`full_run_all.sh` refuses a non-empty one). The
   harnesses skip runs that already have results, so resuming into an older directory would mix conditions
   (another protocol, or another version of an order such as the mod-minimizer's tie rule).

### Wall time or CPU time

Wall time is the headline (what a user waits for). CPU time is the diagnostic: it mostly excludes waiting
(contention, input reads), so a wall-time difference **without** a CPU-time difference is load or cache, not
the order. A step with `*_read_gb` > 0 read from disk (cold page cache) and is not comparable with warm runs;
the page cache is recorded, not controlled. CPU time is not immune to contention either (hyperthread sharing
and memory bandwidth inflate it, which is why pinning uses one thread per core).

### Significance tests (`common/order_tests.py`)

```bash
python common/order_tests.py --app {minimap2,sshash,kraken2,kmc} --results-dir <results dir> [--exclude-cold-gb 1]
```

Per context (dataset, repeat-filter condition, ...) and timed metric (wall and CPU seconds of each step) it
compares every GreedyMini / mod-minimizer variant, xor and xor-pcg-hash with the baseline order and with xor:
Welch's t-test on log times over the repeats (ratio of geometric means, 90% CI, p, Holm-adjusted over the
candidates of a block), plus TOST equivalence at +-5%. No scipy (`common/stats_util.py`). Output:
`<results dir>/csvs/order_tests_<app>.csv` and a per-metric verdict count. `--exclude-cold-gb G` drops the
steps that read more than G GB from disk, a sensitivity check for the page-cache confound.

## Running everything

```bash
benchmark/full_run_all.sh                                    # two slots, default results directory name
SLOTS="0:minimap2,kmc 1:kraken2,sshash" RESULTS_NAME=<name> benchmark/full_run_all.sh   # the default slots, explicit
SLOTS="1:sshash" RESULTS_NAME=<name> benchmark/full_run_all.sh                         # one app, on node 1
RESUME=1 RESULTS_NAME=<name> benchmark/full_run_all.sh       # continue an interrupted run (finished runs are skipped)
```

- **Refusals and exit status:** it refuses to start in a results directory that already holds runs, refuses an
  app listed in two slots, and exits 1 if any app failed. Within an app's `experiments/full_run.sh` a failed step
  is reported and the rest still runs (a re-run resumes); the script then exits 1, which counts as that app failing.
- **Logs:** `<app root>/<RESULTS_NAME>.log`, appended to on `RESUME=1`. Each app's `experiments/full_run.sh` ends
  with its export. The kmc export reuses its expected and per-dataset density caches from `OLD_RESULTS` (default
  `<kmc root>/results`); point it at an archived summary's `csvs/` parent when that directory is gone.
- **Disk:** the peak comes from kmc's temporary files (largest for the human datasets), kraken2's databases (all
  kept until classification ends) and minimap2's indexes. minimap2 keeps only seed 0's index of every variant,
  for the scheme analysis; the others are deleted after mapping and their size kept in `index_bytes.txt`
  (`KEEP_ALL_INDEXES=1` keeps them all).

## Workbook and column naming (all four exporters)

- Units are suffixes: `_s` seconds, `_gb` decimal GB (1 GB = 10^9 bytes), `_pct` percent, `_ns`; counts start
  with `n_`; aggregated sheets use `<metric>_mean`, `<metric>_std`, `<metric>_ratio` and `n_repeats` (number of
  runs aggregated; for the baselines the repeats only measure timing noise).
- Step columns are `{step}_wall_s`, `{step}_cpu_s`, `{step}_max_rss_gb`, `{step}_major_pagefaults`,
  `{step}_read_gb`, `{step}_pinned_cpus`; the split is `{step}_load_s` / `{step}_run_s`. I/O and pinning
  columns are in the runs sheets only.
- The xlsx headers are the CSV headers; the overview sheet has one legend row per column (checked by each
  app's tests).

## How each metric follows the density (`common/density_effect.py`)

Each app's export adds the sheets `density_effect` and `density_effect_points` to its workbook and writes
`<results dir>/csvs/density_effect_<app>.csv` (by hand: `python common/density_effect.py APP=RESULTS_DIR ...
[--out WORKBOOK.xlsx]`).
- **Pairing.** Every run of a non-default scheme is paired with the default scheme's run of the same setting
  (dataset and its parameters) and seed, i.e. of the same round of the protocol; when the default has several
  runs there, their mean is used.
- **Points.** Each pair gives the density factor and the metric as ratios to the default's; the ratios make
  the datasets comparable, whereas raw values would correlate mostly through the size of the dataset.
- **Correlation.** Per (app, metric), the Pearson r of the two ratios over the runs of all settings (`all`)
  and of each setting, with the number of runs and of schemes: r near 1 means that the metric follows the
  density, near 0 that it does not depend on it.
- **minimap2.** Only the fixed-threshold block (`-f` = native's threshold) is paired, found by position as in
  `order_tests.py`: a derived threshold differs between schemes, and one that happens to equal native's does
  not make the run comparable. A result set without a fixed block has no `density_effect` sheets.

## Whether each scheme changes each metric (`common/significance.py`)

Each app's export adds the sheet `significance` to its workbook and writes
`<results dir>/csvs/significance_<app>.csv` (by hand: `python common/significance.py APP=RESULTS_DIR ...`), for
the metrics of the paper's figures (build and processing times, KMC3 stages, peak memory, index size,
temporary disk, bits per k-mer).
- **Pairing.** As in `density_effect`: every run is paired with the default scheme's run of the same setting
  and seed (the same round), and gives d = ln(metric / default's).
- **Test.** Per (setting, metric, scheme), over the seeds: ratio = exp(mean d), the geometric mean of the
  ratios, with its 95% CI; a paired t-test of mean d = 0 (p); Holm-adjusted over the schemes of the same
  setting and metric (p_holm); TOST equivalence at +-5% (the 90% CI of the ratio in [1/1.05, 1.05]).
- **Verdict.** lower / higher when p_holm < 0.05, else equivalent (TOST) or inconclusive.
- Not tested: kraken2 at minimum hit groups 0 and the sshash 13-mer GreedyMini tables (not in the paper).
  Unlike `order_tests.py` (Welch's test, wall and CPU time, run by hand), the test is paired by round.

## Profiling (`<app>/experiments/profile_run.sh`)

To explain the runtime, each app has a profiling script for a finished full run. It first mirrors the
result set `SOURCE_RESULTS` (default `results`) into `RESULTS_NAME` (default `results-prof`) with
`common/mirror_results.py` -- files of 8 MiB or more hard-linked (databases, indexes: only read), smaller ones
copied (logs, caches, tables an export rewrites) -- so the source result set is never modified. It builds an
instrumented binary in a separate directory (the timed binaries are not touched), runs it on the full run's
databases or indexes for the five seeds under the measurement protocol, writes the profiles into the mirror,
and re-exports the mirror, which adds the profile columns and the `significance` sheet. The instrumentation
compiles to nothing in the default builds.

| app | flag | what it measures | file per run | columns |
|---|---|---|---|---|
| kraken2 | `make ENABLE_MINIMIZER_PROF=1` (classify) | time of the three phases of classifying a read pair (selected k-mers into a token stream, hash-table lookups, resolving the hits), lookups, lookups that miss, GreedyMini ranks on the finite (UHS) and infinite path | `results/classify_profile_<ds>_mhg<m>.txt` | `profile` sheet: `scan_s`, `lookup_s`, `resolve_s`, `scan_pct_of_phases`, `lookups_per_mbp`, `ns_per_lookup`, `lookup_miss_pct`, `uhs_pct` |
| sshash | `cmake -DSSHASH_PROFILE=On` | streaming query: lookups that decode a bucket and the positions they decode, lookups that verify positions and the positions they verify (index rebuilt from the run's recorded command, temporary) | `query_profile.log` | `query_prof_*` |
| minimap2 | `CFLAGS+=-DMM_PHASE_PROF` | one-thread mapping of the first 3000 reads (fixed threshold): time of computing the selected k-mers, collecting their hits, chaining, aligning | `map_profile.log` | `map_prof_*_pct` |
| kmc | `tools/bench_signature_impls` | kmc's signature code with the default order and the GreedyMini order (DNA-alphabet and binary table), with the window minimum, per base, one core | `<results>/profile/signature_<ds>_<seed>.txt` | `signature_cost` sheet |

Two build phases come from the timed runs' own logs, without profiling: kraken2's `build_estimate_s` /
`build_db_s` (its log's time stamps) and minimap2's `build_collect_s` / `build_collect_pct` (its "collected
minimizers" stamp).

## Shared code (`common/`)

| File | Purpose |
|---|---|
| `wait_low_load.sh`, `node_gate.py` | wait for a quiet machine / NUMA node before a timed step (CPUs, memory, free 2 MB blocks, compaction) |
| `bench_step.py`, `step_record.py` | run a timed step pinned and memory-bound, and record its node, huge pages, load and disk; the record's columns |
| `pin_cpus.sh`, `pin_cpus.py` | CPU pinning: one NUMA node, one thread per core, spare cores |
| `order_tests.py`, `stats_util.py` | significance tests of order differences; scipy-free t / Welch / Holm |
| `density_effect.py` | how closely each metric follows the density: paired ratios to the default scheme, Pearson r |
| `compare_result_sets.py` | run-to-run stability of two result sets |
| `units.py` | decimal-GB conversions (bytes, KiB, 512-byte blocks) |
| `xlsx_book.py`, `xlsx_util.py`, `xlsx_check.py` | workbook writing and validation |
| `order_names.py`, `orders_sheet.py`, `density.py`, `expected_density_driver.py`, `correlation.py`, `paths.py`, `yaml_util.py`, `experiments_lib.py` | order identities and the shared analysis helpers |
