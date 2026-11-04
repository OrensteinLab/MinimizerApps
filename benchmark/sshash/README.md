# SSHash minimizer-scheme benchmark

Benchmarks [SSHash](https://github.com/jermp/sshash) (a compact k-mer dictionary index) across
its **minimizer scheme**: `native` (sshash's own baseline), `xor`, `xor-pcg-hash`,
`greedymini` (the GreedyMini rank tables, center-extended to the m-mer),
and the strand-symmetric mod-minimizer `mod-minimizer-xor` / `mod-minimizer-greedymini`
(see below and `../../orders/README.md`). Scheme names are the canonical ones shared by all
four benchmarks (`benchmark/common/scheme_names.py`) and never carry parameters: the GreedyMini
table, projection and the mod-minimizer t-mer length are the separate result fields `gm_k`,
`projection`, `gm_table` and `mod_t`. Unlike kmc's `lex`
baseline, sshash has no seed-less scheme — every scheme, including `native`, takes a seed, so this
harness runs every scheme once per seed with no special-casing.

## What is measured

For each dataset, k=31, m=21 (fixed across every scheme -- see "Choice of m" below), and
scheme:

- **construction time + peak RSS**: `sshash build` wrapped in `/usr/bin/time -v`
- **correctness**: `sshash check` (must print `EVERYTHING OK!`)
- **query time, primary**: `sshash bench`'s own internal micro-benchmark (`positive lookup`,
  `negative lookup`, `access`, `iterator`, all in `avg_nanosec_per_kmer`) -- this is the same
  headline query-time metric the SSHash paper itself reports
- **query time, secondary**: `sshash query` against the dataset's real paired FASTQ reads,
  also `/usr/bin/time -v`-wrapped
  (streaming queries extend along the stored string after each hit, for every scheme, so the
  report's `query_n_searches`/`query_n_extensions` split is comparable across schemes)
- **space**: the built `.sshash` file's size on disk (`os.path.getsize`, which should match
  `build`'s own reported `index_size_in_bytes` (exported as `index_size_gb`) up to a small serialization-estimate delta) and
  derived bits/kmer
- **scheme equivalence** (hard-fail gate, not just a warning -- see below): every scheme must
  report the identical `n_kmers` for the same dataset, and must answer a fixed set of
  positive/negative k-mer queries identically

## Dataset

Six collections (k=31 eulertigs) from Zenodo record
[17582116](https://zenodo.org/records/17582116): three whole genomes and three larger collections
(a virus set, a bacterial pangenome and a human pangenome), each paired with the real-read query
FASTQ the paper's benchmark scripts (jermp/kmer_sets_benchmark, `query_sshash_indexes.sh`) use for
it. The first three FASTQs come from ENA (species/accession mapping resolved by hand via the ENA
API's `scientific_name`); the last three are files of the record itself:

| id | collection | eulertigs | query FASTQ |
|---|---|---|---|
| `cod` | Gadus morhua | `cod.k31.eulertigs.fa.gz` (~177 MB) | `SRR12858649.fastq.gz` (~54 MB) |
| `kestrel` | Falco tinnunculus | `kestrel.k31.eulertigs.fa.gz` (~342 MB) | `SRR11449743_1.fastq.gz` (~478 MB) |
| `human` | Homo sapiens | `human.k31.eulertigs.fa.gz` (~883 MB) | `SRR5833294.fastq.gz` (~1.67 GB) |
| `ncbi-virus` | NCBI RefSeq viruses (18,836 assemblies) | `ncbi-virus.k31.eulertigs.fa.gz` (~136 MB) | `ncbi-queries.fastq.gz` (~6 MB) |
| `se` | Salmonella enterica pangenome (AllTheBacteria, 534,751 genomes) | `se.k31.eulertigs.fa.gz` (~510 MB) | `SRR27871075_1.fastq.gz` (~439 MB) |
| `hprc` | Homo sapiens pangenome (HPRC) | `hprc.k31.eulertigs.fa.gz` (~1.88 GB) | `SRR5833294.fastq.gz` (human's file; `download` links it instead of fetching it twice) |

Only the record's k=31 eulertigs files are used here (its k=63 files are out of scope: they need a
second sshash build with `SSHASH_USE_MAX_KMER_LENGTH_63` and GreedyMini tables for another window).
Exact filenames, sizes, and md5s are in `data/datasets.csv` (committed; download URLs are the Zenodo
API's file endpoints and, for the first three query FASTQs, ENA's standard `vol1/fastq/...` path
convention; `download` verifies every md5). `scripts/zenodo.py` can re-fetch the record's file
catalog (`python -m scripts zenodo list-files`) to re-verify or extend those values later -- routine
`download`/`run`/`export` never needs network access to Zenodo itself beyond the one-time dataset
download. m stays 21 for every dataset (the paper tunes m per collection; here one m keeps the
GreedyMini w11 tables applicable and the datasets comparable).

There is also a `smoke` row in `data/datasets.csv` pointing at the in-repo
`sshash/data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz` sample (no download, no query
FASTQ) -- see "Smoke testing" below.

## `.env` setup

```bash
cp benchmark/sshash/.env.example benchmark/sshash/.env
# edit SSHASH_BENCH_ROOT — use a directory outside the repo for large downloads
```

| variable | meaning |
|---|---|
| `SSHASH_BENCH_ROOT` | on-disk root: `data/<id>/` (datasets) and `results/{runs,csvs}/` underneath |
| `SSHASH_DATA_DIR` | optional; default `$SSHASH_BENCH_ROOT/data` (datasets live at `$SSHASH_DATA_DIR/<id>/`) |
| `SSHASH_RESULTS_DIR` | optional; default `$SSHASH_BENCH_ROOT/results` |

CLI `--data-dir`/`--results-dir` override the env file. Existing environment variables are not
overwritten by `.env`.

**Git:** the harness (`scripts/`, `lib/`, `data/*.csv`, `data/queries/*.txt`, `experiments/`,
…) is meant to be committed. `benchmark/.gitignore` excludes `.env`, downloaded eulertigs/FASTQ
under `data/<id>/`, and `results/`.

## Layout

```
benchmark/sshash/
  config.yaml         # seeds, paths, dataset list, GreedyMini table labels
  sshash_bench.sh     # main entry (download / run / export / experiment)
  Makefile            # `make venv` / `make test`
  lib/*.sh            # bash helpers
  scripts/            # python -m scripts …
  experiments/*.yaml  # multi-config sweep definitions
  data/datasets.csv   # dataset catalog (committed)
  data/queries/       # committed small pos/neg k-mer sets for the offline sanity test
  tests/
  README.md
```

| Path | Role |
|---|---|
| `config.yaml` | Shared seeds (same 10 hex literals as `benchmark/kmc/config.yaml`, for cross-tool consistency), paths, dataset list, GreedyMini table labels |
| `sshash_bench.sh` | Main entry (download / run / export / experiment) |
| `lib/*.sh` | Bash helpers: `constants.sh` (shared defaults + `config.yaml` loading), `env.sh` (.env loading), `datasets.sh` (download + md5), `schemes.sh` (scheme -> CLI flags), `run.sh` (build/check/bench/query + equivalence, `/usr/bin/time -v` wrapping), `experiment.sh` (YAML sweep dispatch) |
| `scripts/config.py` | Load `config.yaml` |
| `scripts/schemes.py` | Scheme name -> sshash CLI flags (Python twin of `lib/schemes.sh`) |
| `scripts/equivalence.py` | Fixed positive/negative k-mer query-set generation + per-run equivalence checks |
| `scripts/export.py` | CSV / XLSX export + the equivalence hard-fail gate |
| `scripts/zenodo.py` | Zenodo record 17582116 lookup (maintainers only; not on the normal run path) |
| `scripts/meta.py` | Per-run `meta.json` writer |
| `scripts/common.py` | Shared helpers (`.env` loading, `die()`) |
| `scripts/density.py` | Sliding-window "particular" density + pure-Python re-implementations of sshash's order formulas |
| `scripts/expected_density.py` | Theoretical `density_factor` per scheme, counted on uniformly random DNA by sshash's own minimizer iterators (the `sshash` tool of `benchmark/common/expected_density/`; feeds the `schemes` xlsx sheet); thin adapter over the shared `benchmark/common/expected_density_driver.py` (tool, variant x seed loop, mean/std, cache) |
| `scripts/experiments.py` | Load/expand `experiments/*.yaml` (shared engine in `benchmark/common/experiments_lib.py`) |
| `data/` | Dataset catalog (`datasets.csv`) + committed sanity-test query sets (`queries/`) |
| `experiments/*.yaml` | Multi-config sweep definitions |
| `tests/` | Offline sanity checks: `test_scheme_sanity.py` (scheme/equivalence correctness against the in-repo sample) + `test_export.py` (`export.py`'s aggregation/validation/CSV/XLSX logic, synthetic rows, no sshash binary needed); `make venv` once, then `make test` |

Python CLI (from the repo root, `BENCHMARK_DIR=benchmark/sshash` set by `sshash_bench.sh`):

```bash
python3 -m scripts export --results-dir …
python3 -m scripts experiment list
python3 -m scripts schemes flags xor 0xDEADBEEF
python3 -m scripts schemes variants all --k 31 --m 21   # the (scheme, table, projection) runs
python3 -m scripts zenodo list-files
```

## Usage

Requires **Python 3.9+** (a `.venv` is auto-preferred if present at `benchmark/sshash/.venv` --
kept inside the harness directory, not under `sshash/` itself, since that's upstream-tracked
source this harness doesn't touch; `make -C benchmark/sshash venv` creates it) and **Bash**.

```bash
make -C sshash all                                  # build sshash (cmake + make); do not edit its sources
cp benchmark/sshash/.env.example benchmark/sshash/.env   # edit SSHASH_BENCH_ROOT
benchmark/sshash/sshash_bench.sh download -d cod
benchmark/sshash/sshash_bench.sh run -d cod -n 1     # quick smoke: 1 seed, all 9 schemes
# → $SSHASH_BENCH_ROOT/results/csvs/sshash_{runs,summary,...}.csv
# → $SSHASH_BENCH_ROOT/results/sshash_results.xlsx
```

### `sshash_bench.sh` (main entry)

```bash
benchmark/sshash/sshash_bench.sh -h
benchmark/sshash/sshash_bench.sh download -d cod
benchmark/sshash/sshash_bench.sh run -d cod -n 1                     # 1 seed, all schemes
benchmark/sshash/sshash_bench.sh run -d cod,kestrel,human,ncbi-virus,se,hprc -S all -n 5 -t 16   # full sweep
benchmark/sshash/sshash_bench.sh export                              # re-export CSV/XLSX only
benchmark/sshash/sshash_bench.sh all -d cod -n 1                      # download + run + export
```

Defaults: `cod,kestrel,human,ncbi-virus,se,hprc`, k=31, m=21, `-S all` (native, xor, xor-pcg-hash, greedymini,
mod-minimizer-xor, mod-minimizer-greedymini), `-M canonical`, `-n 5` seeds, 8 threads.

The command name (`download`/`run`/`export`/`all`/`experiment`) must come **first**, before any
options -- e.g. `sshash_bench.sh run --dataset-file PATH -n 1`, not the reverse.

`-M/--modality` selects `canonical` (default, single-probe) or `regular` (two-probe, `sshash
build --regular`), a separate axis from the scheme, e.g. `-S native -M regular`. `--force` re-runs a config even if its `meta.json` already exists (default:
skip unchanged configs).

**Experiments** (multi-config sweeps from YAML):

```bash
benchmark/sshash/sshash_bench.sh experiment list
benchmark/sshash/sshash_bench.sh experiment sweep_schemes   # primary: all schemes x 5 seeds, all 6 datasets
benchmark/sshash/sshash_bench.sh experiment sweep_m        # secondary/optional: native only, m in {7,9,11,13}
```

Experiment files live in `experiments/`. Each has shared `defaults` plus a `configs` list; keys
are `datasets`, `k`, `m`, `schemes`, `n_seeds`, `threads`, `force`, `run_tag`.

### Smoke testing

`--dataset-file PATH` registers an ad hoc local FASTA under dataset id `adhoc` (no download, no
query FASTQ, bypasses `data/datasets.csv` entirely):

```bash
benchmark/sshash/sshash_bench.sh run --dataset-file sshash/data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz -n 1
```

Equivalently, `data/datasets.csv` has a committed `smoke` row pointing at the same in-repo
sample via its `local_path` column:

```bash
benchmark/sshash/sshash_bench.sh run -d smoke -n 1
```

Both are useful when a real Zenodo dataset is too slow/expensive to build in the moment --
the whole `run`/`export` pipeline (build, check, bench, query-equivalence, CSV/XLSX, the
hard-fail gate) runs unmodified against either.

### Choice of m

**m is fixed at 21 for every scheme** in the primary comparison (`sweep_schemes.yaml` /
`run -S all`), by deliberate design: the point of this benchmark is isolating the *scheme's*
effect on space/time at a fixed minimizer length, not conflating it with an m sweep. m=21 gives
w=k-m+1=11, the window the GreedyMini tables `w11_k13`/`w11_k15` were built for, and using the same
(k, m) for every scheme keeps the comparison apples-to-apples.
`greedymini` runs `orders/GreedyMini/w11_k13_canonical.gm` and `w11_k15_canonical.gm` (two variants),
extended to `l=21` around the center `gm_k`-mer (`m - gm_k` even). Varying
m *is* useful as its own, separate question -- `experiments/sweep_m.yaml` covers that, explicitly
at a single scheme (`native`) so it isn't mistaken for a second scheme comparison, and is
optional/secondary to the primary `sweep_schemes` results.

### Equivalence checking

Two independent hard-fail gates run at `export` time (mirroring kmc's
`validate_custom_vs_lex`, generalized here since sshash has no seed-less baseline to compare
everything else *against* -- every scheme, including `native`, is just another group member):

1. **`n_kmers` must be identical** across every scheme/seed run of the same dataset. sshash
   already asserts this internally in a sense (`build --check`/`sshash check` verify the index
   is self-consistent), but this is the cross-scheme check that the *set of k-mers indexed*
   didn't silently change because of the scheme.
2. **Fixed positive/negative k-mer query answers must be identical.** Once per dataset,
   `scripts/equivalence.py` builds a throwaway `native`-scheme reference index, samples ~10,000
   real k-mers from the input as positives and ~10,000 random 31-mers verified absent (via
   `sshash query` against that reference index) as negatives, and caches them under
   `$SSHASH_DATA_DIR/<id>/queries/{pos,neg}.txt`. Every subsequent scheme/seed run re-queries
   both fixed sets against its own index and records `query_pos_match`/`query_neg_match`
   booleans (all-present / all-absent, from `sshash query`'s own aggregate
   `query_n_positive_kmers`/`query_n_negative_kmers` counts -- it doesn't report per-k-mer answers, so
   the aggregate count against a fixed all-positive or all-negative file is the check).
   `sshash check`'s own `EVERYTHING OK!` is also folded into this gate (`check_ok`).

`export` (`scripts/export.py:validate_schemes_equivalence`) reads these fields back from every
run's `meta.json` and `die()`s hard, printing every violation, if any is false or if
`n_kmers` disagrees -- not a warning, a failed export.

The offline sanity test (`tests/test_scheme_sanity.py::test_sshash_integration`) exercises the
same logic end-to-end against the in-repo sample dataset, using the **committed**
`data/queries/salmonella_enterica_k31_{pos,neg}.txt` (3000 k-mers each, generated the same way
against a native-scheme build of that sample) instead of generating a fresh set, so the test has
no network dependency and a fixed expected `n_kmers` (4787534) to assert against.

## Outputs

On-disk layout (no wrapping folder under the root):

```
$SSHASH_BENCH_ROOT/
  data/<dataset>/      # eulertigs, query FASTQ, queries/{pos,neg}.txt
  results/runs/        # per-run artifacts (below)
  results/csvs/        # exported tables (sshash_*.csv, expected-density cache)
  results/sshash_results.xlsx   # workbook
```

Per-run artifacts live under `$SSHASH_BENCH_ROOT/results/runs/` in a nested layout:

```
runs/
  {dataset}/
    {modality}/            # canonical | regular -- see -M/--modality
      {scheme}[/params]/    # native | xor | xor-pcg-hash | mod-minimizer-xor/t10 |
                            # greedymini/k15_at0-cg1 | mod-minimizer-greedymini/k10_at0-cg1_t10
                            # (scheme_names.scheme_dir_name: k<gm_k>_<projection>_t<mod_t>)
        {seed_idx}/          # 0…n-1
          build.log          # sshash build --verbose --bucket-stats stdout+stderr (trailing JSON blob)
          build_time.txt     # GNU /usr/bin/time -v
          check.log          # sshash check --verbose (EVERYTHING OK! or diagnostics)
          bench.log          # sshash bench (trailing JSON blob)
          query.log          # sshash query vs the real dataset FASTQ (if configured)
          query_time.txt     # GNU /usr/bin/time -v for query.log
          equivalence.json   # {"pos": {...}, "neg": {...}} fixed-set check results
          meta.json          # run identity + config + exit codes + command (k/m recorded here, not in the path)
          # index.sshash is deleted after check/bench/query/equivalence unless --keep-index
```

`export` (also run at the end of `run`/`all`) writes CSVs to
`$SSHASH_BENCH_ROOT/results/csvs/` and the workbook to `$SSHASH_BENCH_ROOT/results/`:

| file | rows | purpose | in `.xlsx` |
|---|---|---|---|
| `sshash_runs.csv` | one per run | all metrics + telemetry | `runs` |
| `sshash_schemes.csv` | one per scheme variant (k, w, L) | description, CLI, expected + particular density_factor | `schemes` |
| `sshash_summary.csv` | one per dataset x scheme variant | mean/std across seeds, ratio vs native | `summary` |
| `sshash_space_breakdown.csv` | one per dataset x scheme variant | bits/k-mer by component, saving vs native, share of the saving | `space_breakdown` |
| `significance_sshash.csv` | one per metric x dataset x scheme | paired test of the change from native (`benchmark/common/significance.py`, see `benchmark/README.md`) | `significance` |
| `density_effect_sshash.csv` | one per metric x (all, dataset) | Pearson r of the metric and the density factor, both as ratios to native's run of the same seed (`benchmark/common/density_effect.py`) | `density_effect` (points: `density_effect_points`) |
| `sshash_results.xlsx` | — | `overview`, `schemes`, `datasets`, `summary`, `runs`, `buckets`, `space_breakdown`, `bucket_histogram`, `super_kmer_histogram`, `density_effect`, `density_effect_points`, `significance` (benchmark/common/xlsx_book.py) | — |

Re-export only: `benchmark/sshash/sshash_bench.sh export` (no network needed once local run
data exists).

### CSV columns (main fields)

**Shared keys:** `dataset`, `k` (minimizer length, sshash's m), `w` (window), `L` (= k+w-1, sshash's k-mer length), the scheme variant `scheme` (`native`/`xor`/`xor-pcg-hash`/`greedymini`/
`mod-minimizer-xor`/`mod-minimizer-greedymini`), `gm_k`, `projection` (`at0-cg1`/`at1-cg0`),
`gm_table` (table file stem), `mod_t` (blank where they do not apply), then the app parameters `modality`
(`canonical`/`regular`, see `-M/--modality`) and, in runs, `threads`, then the seed columns (`n_repeats`, the number of runs aggregated into the row's mean/std, in
aggregated sheets; `seed`, `seed_idx` in runs). Summary, bucket and space rows are grouped by all of these.

**Space:** `n_kmers`, `index_size_gb` (from `build`'s own JSON), `index_disk_gb`
(`os.path.getsize`, the ground truth), `bits_per_kmer` (`index_disk_gb * 1e9 * 8 / num_kmers`). All sizes and memory are decimal GB (1 GB = 10^9 bytes, `benchmark/common/units.py`); meta.json and the run dirs keep the raw bytes / KiB, converted at export.

**Density:** `density_factor` and `distinct_minimizers` (M, same value as `n_buckets`), `n_minimizer_positions` (runs).

**Build:** `build_wall_s`/`build_max_rss_gb` (GNU time, KiB -> GB), `build_steps_total_s` and
per-step `build_step{1..7}_s` (from `build --verbose`'s own microsecond timings).

**Query (primary, `sshash bench`):** `bench_positive_ns`, `bench_negative_ns`,
`bench_access_ns`, `bench_iterator_ns` (all `avg_nanosec_per_kmer`).

**Query (secondary, real FASTQ):** `query_wall_s`/`query_max_rss_gb` (GNU time),
`query_run_s` (sshash's own timed query phase, seconds), `query_load_s` (= `query_wall_s - query_run_s`), `query_n_kmers`/`query_n_positive_kmers`/`query_n_negative_kmers`.

**Correctness/equivalence:** `check_ok`, `query_pos_match`, `query_neg_match` (all booleans;
`export` hard-fails if any is false for any run).

### XLSX sheets

| sheet | content |
|---|---|
| `overview` | column reference: what every non-obvious column means |
| `schemes` | one row per (scheme variant, k, w, L) with the unified names shared by all apps (k = minimizer length = 21, w = window = 11, L = k + w - 1 = 31 = sshash's k-mer length): scheme, k, w, L, gm_k, projection, gm_table, mod_t, description, cli, expected + per-dataset particular density_factor mean/std (no empty density cell, `export` fails otherwise; also written as `sshash_schemes.csv`; layout in `benchmark/common/schemes_sheet.py`) |
| `datasets` | one row per dataset: k, w, L, measured n_kmers |
| `summary` | mean/std per dataset x scheme across seeds, plus `{metric}_ratio` = this variant's mean / the same-modality native's mean (native = 1); green (better) → red (worse) color scale on the ratio-scale metrics. Two rows right under the table hold, for every `{metric}_ratio` column, its Excel `CORREL` with the density ratio and with the distinct-minimizer ratio (`num_buckets_ratio`) over the dataset x scheme rows. |
| `space_breakdown` | one row per dataset x scheme: where the bits/k-mer go (strings, offsets, MPHF, tags, light array L, skew index f_i+V_i and H), each component's saving over native and its share of the total saving, and the saving split into what is sized by the distinct minimizers (MPHF + tags) and what by the occurrences (L + skew index) |
| `runs` | one row per run (every metric + telemetry, unfiltered); after `experiments/profile_run.sh` also the streaming-query profile `query_prof_full_lookups`, `query_prof_positions_per_full_lookup`, `query_prof_verify_calls`, `query_prof_positions_per_verify` (the two ratios also on `summary`) |

## Measurement protocol

Build, bench and query are timed under the protocol in [`../README.md`](../README.md): node-quiet wait (bench and
query also wait for free huge-page memory), pinning (build `$THREADS` cores + 2 spare, the single-threaded bench
and query 1 + 1; one NUMA node, one thread per core; lists in `cpus.txt`, `cpus_bench.txt`, `cpus_query.txt`,
exported as `build_/bench_/query_pinned_cpus`; with its memory bound to the slot's node and a step record (`step_*.json`, exported as `build|bench|query_node`, `_huge_fraction`, `_foreign_busy_cores`, ... -- see the protocol)) and per-step wall and CPU
seconds, RSS, major page faults and reads (runs sheets only). Query wall time splits into `query_load_s`
(index load and start-up) and `query_run_s` (sshash's own timed queries, seconds). `native` runs once per seed
like every other scheme. `experiments/full_run.sh` needs the default-encoding binary
(`sshash/build/sshash`, `SSHASH_BIN` overrides; a build with the traditional nucleotide encoding cannot run
GreedyMini).

## Notes

- Every scheme (including `native`) takes `-s SEED`; there is no seed-less baseline the way
  kmc's `lex` is, so this benchmark never special-cases a "run once" scheme the way kmc does.
- `greedymini` (sshash `--scheme greedymini --greedy-mini-file TABLE`, the center extension of the
  table to the m-mer, see `orders/README.md`) runs once per variant, a (table, projection)
  pair from `config.yaml`'s `greedymini` section (`tables`: `w11_k13_canonical.gm`,
  `w11_k15_canonical.gm`; `projections`: `at0-cg1`). Add `at1-cg0` there (or pass
  `--gm-projections at0-cg1,at1-cg0`; `--gm-tables` likewise selects tables) for the flipped
  DNA->binary polarity (`sshash build --gm-projection at1-cg0`: AT->1, CG->0 instead of the
  default AT->0, CG->1). Each variant is its own run directory (`greedymini/k15_at1-cg0/...`),
  summary/schemes row and `meta.json` entry with its own `gm_k`, `projection` and `gm_table`.
  The center extension needs `m - gm_k` even, so e.g. a table with an even `gm_k` is rejected at
  m=21; `sshash_bench.sh run` checks this up front (`python -m scripts schemes check-gm greedymini M
  TABLE`). Example:
  `./sshash_bench.sh run -d cod,kestrel -n 5 -S native,greedymini --gm-projections at0-cg1,at1-cg0`.
- Strand-symmetric mod-minimizer schemes (sshash `--scheme mod-minimizer-xor` / `--scheme
  mod-minimizer-greedymini` with `--mod-t T`, see `orders/README.md`): `mod-minimizer-xor` ranks
  the k-mer's t-mers by xor, `T = 4 + (m - 4) mod w` (10 at k=31, m=21);
  `mod-minimizer-greedymini` by a GreedyMini table of t-mers from `config.yaml`'s
  `mod_minimizer_greedymini` section (`w22_k10`, t=10), `T` = the table's gm_k. Both are in
  `all`; `--regular` does not apply to them. With
  t = m mod w (both here) the scheme is forward, which the build checks; any other t makes the
  two strands' offset rules disagree (a t=11 table gave density ~2.2, worse than native).
  The schemes sheet's expected density is counted on uniformly random DNA (equal for both polarities
  up to sampling noise by symmetry); the particular density per dataset is what can differ.
- The built `.sshash` index is deleted after each run's check/bench/query/equivalence steps
  finish (to keep disk usage bounded across cod/kestrel/human x 9 schemes x N seeds x up to 2
  `t` values); pass `--keep-index` to retain it.
- `sshash build`/`bench`/`query --verbose` all print one trailing single-line JSON blob (values
  are JSON *strings*, parsed with `json.loads` then cast) -- `lib/run.sh` captures each tool's
  combined stdout+stderr to its own `.log` file, and `scripts/export.py` finds that trailing
  `{...}` via regex rather than assuming it's the literal last line of a particular stream (it
  was observed on stderr in hand-testing, not stdout).
- Query FASTQ download uses ENA's standard `vol1/fastq/<prefix6>[/<subfolder>]/<accession>/<file>`
  path convention, generalized (`lib/datasets.sh:ena_url`) beyond `benchmark/kmc`'s version to
  handle accessions needing the extra numeric subfolder (7/8/9-digit run numbers).
