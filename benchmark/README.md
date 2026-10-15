# Benchmark harnesses: shared protocol and tools

Four apps are benchmarked across minimizer orders (random xor / xor-pcg-hash, GreedyMini, mod-minimizer
variants, and each app's own baseline): `minimap2/` (read mapping), `sshash/` (k-mer dictionary),
`kraken2/` (taxonomic classification) and `kmc/` (k-mer counting). Each folder has its own README (data,
`.env`, usage, outputs); this file describes what they share.

## Measurement protocol

Timings on a shared host are only comparable if every timed step is made under the same conditions.
Every timed step of every app follows the same protocol:

1. **Quiet machine first** (`common/wait_low_load.sh`). Before each timed step the harness waits until fewer
   than `BENCH_LOAD_MAX` (default 4) cores are busy, measured from `/proc/stat` over `BENCH_LOAD_WINDOW`
   seconds (default 5) while the benchmark itself is idle, so it sees other users' load, not its own.
   The busy-core figure is printed (`[load] ...`) and lands in the run log. `BENCH_WAIT_LOW_LOAD=0` disables it.
2. **CPU pinning, always** (`common/pin_cpus.sh`, selection in `common/pin_cpus.py`). The step runs under
   `taskset -c <list>` placed *inside* `/usr/bin/time -v`, so time and peak RSS are the program's own.
   The list is chosen per step, needs no privileges (only `taskset`, `/proc`, `/sys`), and follows one rule set:
   - **one NUMA node**: the node whose best cores are least busy (a step never spans sockets, and Linux's
     local-first memory policy then keeps its memory on that node);
   - **one logical CPU per physical core**: a core counts as busy if *either* hardware thread is, so a pinned
     thread never shares its core with its own sibling or with someone else's job;
   - **N threads + 2 spare cores** (`BENCH_PIN_SPARE`) for helper threads (readers, I/O, the main thread).
   Pinning is not optional: if `taskset` or the selection is unavailable the benchmark stops instead of timing
   unpinned. `BENCH_CPUS=<list>` forces an explicit list. The list used is stored with the run
   (`cpus*.txt` in the run directory, exported as `{step}_pinned_cpus`).
   What is pinned, per app: minimap2 build and map (and `analysis/run_seed_experiments.sh`, `-t 1`); sshash
   build, bench and query; kraken2 build and classify (classify is `--threads 1`: 2 cores + 2 spare);
   kmc the counting call.
3. **Repeats for every order, baseline included.** Each order runs for each seed in `config.yaml` (5): minimap2
   `native` and kmc `default` are repeated like the others (their repeats only measure timing noise), so
   every comparison has replicates on both sides. sshash `native` was already run per seed; kraken2 has no
   native run, so its reference order is `xor`. The orders of a repeat are interleaved (kmc: one
   `-O all -n 5` call), not run as separate passes, so cache and load drift hits all orders alike.
4. **What is recorded for every step** (`/usr/bin/time -v`): wall seconds, CPU seconds (`*_cpu_s` = user +
   system, all threads), peak RSS (`*_max_rss_gb`), major page faults and file-system reads
   (`*_major_pagefaults`, `*_read_gb`, decimal GB from 512-byte blocks), and where available the load/run
   split (`map_load_s`/`map_run_s`, `classify_load_s`/`classify_run_s`, `query_load_s`/`query_run_s`).

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
benchmark/full_run_all.sh                        # minimap2, sshash, kraken2, kmc -> <app root>/results_r2
RESULTS_NAME=results_r3 APPS="kmc sshash" benchmark/full_run_all.sh
RESUME=1 benchmark/full_run_all.sh               # continue an interrupted run (finished runs are skipped)
```

One app after the other (never two at once). It refuses to start in a results directory that already holds
runs, so runs made under other conditions (e.g. unpinned) cannot mix with new ones. Per-app logs:
`<app root>/<RESULTS_NAME>.log`; each app's `experiments/full_run.sh` ends with its export. Expect days of
timed work (kraken2 builds and the kmc `h_sapiens_3` runs dominate) and roughly 0.9 TB of disk: minimap2 keeps
every `index.mmi` (there is no index trimming) and kraken2 every built database.

## Workbook and column naming (all four exporters)

- Units are suffixes: `_s` seconds, `_gb` decimal GB (1 GB = 10^9 bytes), `_pct` percent, `_ns`; counts start
  with `n_`; aggregated sheets use `<metric>_mean`, `<metric>_std`, `<metric>_ratio` and `n_repeats` (number of
  runs aggregated; for the baselines the repeats only measure timing noise).
- Step columns are `{step}_wall_s`, `{step}_cpu_s`, `{step}_max_rss_gb`, `{step}_major_pagefaults`,
  `{step}_read_gb`, `{step}_pinned_cpus`; the split is `{step}_load_s` / `{step}_run_s`. I/O and pinning
  columns are in the runs sheets only.
- The xlsx headers are the CSV headers; the overview sheet has one legend row per column (checked by each
  app's tests).

## Shared code (`common/`)

| File | Purpose |
|---|---|
| `wait_low_load.sh` | wait for a quiet machine before a timed step |
| `pin_cpus.sh`, `pin_cpus.py` | CPU pinning: one NUMA node, one thread per core, spare cores |
| `order_tests.py`, `stats_util.py` | significance tests of order differences; scipy-free t / Welch / Holm |
| `units.py` | decimal-GB conversions (bytes, KiB, 512-byte blocks) |
| `xlsx_book.py`, `xlsx_util.py`, `xlsx_check.py` | workbook writing and validation |
| `order_names.py`, `orders_sheet.py`, `density.py`, `expected_density_driver.py`, `correlation.py`, `paths.py`, `yaml_util.py`, `experiments_lib.py` | order identities and the shared analysis helpers |
