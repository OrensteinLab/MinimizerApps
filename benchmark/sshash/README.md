# SSHash minimizer-scheme benchmark

Benchmarks [SSHash](https://github.com/jermp/sshash) (a compact *k*-mer dictionary index) across its minimizer scheme: `native` (SSHash's own), `xor`, `xor-pcg-hash`, `greedymini` (the GreedyMini rank tables, center-extended to the selected *k*-mer), and the canonical mod-minimizers `mod-minimizer-xor` / `mod-minimizer-greedymini` (see `../../orders/README.md`). Scheme names are the canonical ones shared by all four benchmarks (`benchmark/common/scheme_names.py`) and never carry parameters: the GreedyMini table, the polarity of Φ and the mod-minimizer *t*-mer length are the separate result fields `gm_k`, `projection`, `gm_table` and `mod_t`. Every scheme, including `native`, takes a seed, so every scheme variant runs once per seed.

**The paper's runs** are `experiments/full_run.sh`: all six datasets, all schemes and variants (`-S all`), 5 seeds, 16 threads, canonical modality, *L* = 31, *k* = 21, *w* = 11. The paper reports `greedymini` with the table `w11_k15` and polarity of Φ `at0-cg1`; the `w11_k13` table is an additional run. The defaults of `sshash_bench.sh` match these settings.

## Notation

The paper's notation is used throughout, and the CSVs and xlsx use it as columns `k, w, L`:

| paper | meaning | SSHash name |
|---|---|---|
| `k` | selected *k*-mer length (21) | minimizer length `m` (`-m`) |
| `w` | number of *k*-mers per window, `w = L - k + 1` (11) | — |
| `L` | SSHash's *k*-mer (*L*-mer) length, `L = k + w - 1` (31) | `k` (`-k`) |
| `t` | mod-minimizer *t*-mer length (10) | `--mod-t` |

SSHash's flags, `meta.json` and the `sshash` commands keep SSHash's own names (`-k 31 -m 21`); the exporter renames on output.

## What is measured

For each dataset and scheme variant, at *L* = 31, *k* = 21 (the same for every scheme, see "Choice of *k*"):

- **construction time + peak RSS**: `sshash build` wrapped in `/usr/bin/time -v`
- **correctness**: `sshash check` (must print `EVERYTHING OK!`)
- **query time, primary**: `sshash bench`'s micro-benchmark (`positive lookup`, `negative lookup`, `access`, `iterator`, all in `avg_nanosec_per_kmer`), the query-time metric the SSHash paper reports
- **query time, secondary**: `sshash query` against the dataset's real FASTQ reads, also `/usr/bin/time -v`-wrapped (streaming queries extend along the stored string after each hit, for every scheme, so `query_n_searches` / `query_n_extensions` are comparable across schemes)
- **space**: the built `.sshash` file's size on disk (`os.path.getsize`; matches `build`'s reported `index_size_in_bytes`, exported as `index_size_gb`, up to a small serialization-estimate delta) and bits/*k*-mer
- **scheme equivalence** (hard-fail gate, see below): every scheme must report the same `n_kmers` for a dataset, and answer a fixed set of positive/negative *k*-mer queries identically

## Datasets

Six collections (*L* = 31 eulertigs) from Zenodo record [17582116](https://zenodo.org/records/17582116): three whole genomes and three larger collections (a virus set, a bacterial pangenome and a human pangenome), each paired with the real-read query FASTQ that the SSHash paper's benchmark scripts (jermp/kmer_sets_benchmark, `query_sshash_indexes.sh`) use for it. The first three FASTQs come from ENA (species/accession mapping via the ENA API's `scientific_name`); the last three are files of the record itself:

| id | collection | eulertigs | query FASTQ |
|---|---|---|---|
| `cod` | Gadus morhua | `cod.k31.eulertigs.fa.gz` (~177 MB) | `SRR12858649.fastq.gz` (~54 MB) |
| `kestrel` | Falco tinnunculus | `kestrel.k31.eulertigs.fa.gz` (~342 MB) | `SRR11449743_1.fastq.gz` (~478 MB) |
| `human` | Homo sapiens | `human.k31.eulertigs.fa.gz` (~883 MB) | `SRR5833294.fastq.gz` (~1.67 GB) |
| `ncbi-virus` | NCBI RefSeq viruses (18,836 assemblies) | `ncbi-virus.k31.eulertigs.fa.gz` (~136 MB) | `ncbi-queries.fastq.gz` (~6 MB) |
| `se` | Salmonella enterica pangenome (AllTheBacteria, 534,751 genomes) | `se.k31.eulertigs.fa.gz` (~510 MB) | `SRR27871075_1.fastq.gz` (~439 MB) |
| `hprc` | Homo sapiens pangenome (HPRC) | `hprc.k31.eulertigs.fa.gz` (~1.88 GB) | `SRR5833294.fastq.gz` (human's file; `download` links it instead of fetching it twice) |

Only the record's *L* = 31 eulertigs are used (its *L* = 63 files would need a second SSHash build with `SSHASH_USE_MAX_KMER_LENGTH_63` and GreedyMini tables for another window). Exact filenames, sizes and md5s are in `data/datasets.csv` (download URLs are the Zenodo API's file endpoints and, for the first three query FASTQs, ENA's `vol1/fastq/...` paths; `download` verifies every md5). `python -m scripts zenodo list-files` (`scripts/zenodo.py`) re-fetches the record's file catalog to re-verify or extend those values; `download`/`run`/`export` need no other access to Zenodo. *k* = 21 for every dataset (the SSHash paper tunes m per collection; one *k* keeps the GreedyMini *w* = 11 tables applicable and the datasets comparable).

`data/datasets.csv` also has a `smoke` row pointing at the in-repo `sshash/data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz` sample (no download, no query FASTQ); see "Smoke testing".

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

CLI `--data-dir`/`--results-dir` override the env file. Existing environment variables are not overwritten by `.env`. `benchmark/.gitignore` excludes `.env`, the downloaded eulertigs/FASTQ under `data/<id>/`, and `results/`.

## Layout

```
benchmark/sshash/
  config.yaml         # seeds, paths, dataset list, GreedyMini tables and projections
  sshash_bench.sh     # main entry (download / run / export / all / experiment)
  Makefile            # make venv / make test
  lib/*.sh            # bash helpers
  scripts/            # python -m scripts …
  experiments/        # full_run.sh, profile_run.sh, *.yaml sweep definitions
  data/datasets.csv   # dataset catalog
  data/queries/       # small pos/neg k-mer sets for the offline sanity test
  tests/
```

| Path | Role |
|---|---|
| `config.yaml` | Seeds (the same 10 hex literals as `benchmark/kmc/config.yaml`), paths, dataset list, GreedyMini tables and projections |
| `sshash_bench.sh` | Main entry |
| `lib/*.sh` | Bash helpers: `constants.sh` (run defaults + `config.yaml` loading), `env.sh` (.env loading), `datasets.sh` (download + md5), `schemes.sh` (scheme -> CLI flags), `run.sh` (build/check/bench/query + equivalence, `/usr/bin/time -v` wrapping), `experiment.sh` (YAML sweep dispatch) |
| `scripts/config.py` | Load `config.yaml` |
| `scripts/schemes.py` | Scheme name -> sshash CLI flags (Python twin of `lib/schemes.sh`) |
| `scripts/equivalence.py` | Fixed positive/negative *k*-mer query sets + per-run equivalence checks |
| `scripts/export.py` | CSV / XLSX export + the equivalence hard-fail gate |
| `scripts/zenodo.py` | Zenodo record 17582116 lookup (maintainers only; not on the run path) |
| `scripts/meta.py` | Per-run `meta.json` writer |
| `scripts/common.py` | Shared helpers (`.env` loading, `die()`) |
| `scripts/density.py` | Sliding-window particular density + Python re-implementations of SSHash's scheme formulas |
| `scripts/expected_density.py` | Expected `density_factor` per scheme variant, counted on uniformly random DNA by SSHash's own minimizer iterators (the `sshash` tool of `benchmark/common/expected_density/`; feeds the `schemes` sheet); thin adapter over `benchmark/common/expected_density_driver.py` |
| `scripts/experiments.py` | Load/expand `experiments/*.yaml` (shared engine in `benchmark/common/experiments_lib.py`) |
| `experiments/full_run.sh` | The paper's runs, then export |
| `experiments/profile_run.sh` | Streaming-query profile of a finished full run (the `query_prof_*` columns) |
| `experiments/*.yaml` | Multi-config sweep definitions |
| `tests/` | Offline checks: `test_scheme_sanity.py` (scheme/equivalence correctness on the in-repo sample) + `test_export.py` (export aggregation/validation/CSV/XLSX on synthetic rows, no sshash binary needed); `make venv` once, then `make test`, from `benchmark/sshash` |

Python CLI (from `benchmark/sshash`; `sshash_bench.sh` runs it there):

```bash
python3 -m scripts export --results-dir …
python3 -m scripts experiment list
python3 -m scripts schemes flags xor 0xDEADBEEF
python3 -m scripts schemes variants greedymini,mod-minimizer-greedymini --k 31 --m 21   # the (scheme, table, projection) runs
python3 -m scripts zenodo list-files
```

## Usage

Requires Python >= 3.9, Bash, and `curl` and `md5sum` for the downloads. The scripts use `benchmark/sshash/.venv` when it exists, otherwise `python3` (override with `PYTHON=`); `make -C benchmark/sshash venv` creates it with pytest. They never write into the `sshash/` fork.

The GreedyMini schemes read the `.gm` tables (`w11_k13`, `w11_k15`, `w22_k10` under `orders/GreedyMini/`), which these scripts do not build: build them first as in the root [README](../../README.md) (the `gm_serializer` loop over `orders/GreedyMini/*.txt`).

From the repo root:

```bash
make -C sshash all                                       # sshash/build/sshash, the default binary of `sshash_bench.sh`
make -C benchmark/sshash venv                            # benchmark/sshash/.venv with pytest
cp benchmark/sshash/.env.example benchmark/sshash/.env   # set SSHASH_BENCH_ROOT
benchmark/sshash/sshash_bench.sh download -d cod
benchmark/sshash/sshash_bench.sh run -d cod -n 1         # quick smoke: 1 seed, all 6 schemes (10 variants)
# → $SSHASH_BENCH_ROOT/results/csvs/sshash_*.csv and results/sshash_results.xlsx

benchmark/sshash/experiments/full_run.sh                 # the paper's runs; RESULTS_NAME=<dir under SSHASH_BENCH_ROOT>, default results
```

### `sshash_bench.sh`

From `benchmark/sshash`:

```bash
./sshash_bench.sh -h
./sshash_bench.sh download -d cod
./sshash_bench.sh run -d cod -n 1                       # 1 seed, all schemes
./sshash_bench.sh run -S all -n 5 -t 16                 # all six datasets, as experiments/full_run.sh
./sshash_bench.sh export                                # re-export CSV/XLSX only
./sshash_bench.sh all -d cod -n 1                       # download + run + export
```

Defaults: `cod,kestrel,human,ncbi-virus,se,hprc`, `-k 31 -m 21`, `-S all` (native, xor, xor-pcg-hash, greedymini, mod-minimizer-xor, mod-minimizer-greedymini), `-M canonical`, `-n 5` seeds, 16 threads. `-S all` runs 10 variants: native, xor, xor-pcg-hash, mod-minimizer-xor, greedymini × (w11_k13, w11_k15) × (at0-cg1, at1-cg0), and mod-minimizer-greedymini (w22_k10) × (at0-cg1, at1-cg0).

The command (`download`/`run`/`export`/`all`/`experiment`) must come first, before any options, e.g. `sshash_bench.sh run --dataset-file PATH -n 1`.

`-M/--modality` selects `canonical` (default, single-probe) or `regular` (two-probe, `sshash build --regular`), a separate axis from the scheme, e.g. `-S native -M regular`. `--force` re-runs a config even if its `meta.json` exists (default: skip).

**Experiments** (multi-config sweeps from YAML), from `benchmark/sshash`:

```bash
./sshash_bench.sh experiment list
./sshash_bench.sh experiment sweep_schemes   # all schemes x 5 seeds, all 6 datasets, 16 threads (the runs of full_run.sh)
./sshash_bench.sh experiment sweep_m         # optional: native only, m in {7,9,11,13}, cod
```

The paper's entry point is `experiments/full_run.sh`. Each experiment file has shared `defaults` plus a `configs` list; keys are `datasets`, `k`, `m`, `schemes`, `n_seeds`, `threads`, `force`, `run_tag`.

### Smoke testing

`--dataset-file PATH` registers an ad hoc local FASTA under dataset id `adhoc` (no download, no query FASTQ, bypasses `data/datasets.csv`). From `benchmark/sshash`:

```bash
./sshash_bench.sh run --dataset-file ../../sshash/data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz -n 1
./sshash_bench.sh run -d smoke -n 1     # the same sample via data/datasets.csv's smoke row (local_path)
```

The whole `run`/`export` pipeline (build, check, bench, query equivalence, CSV/XLSX, the hard-fail gate) runs unmodified against either.

### Choice of *k*

*k* (SSHash's m) is 21 for every scheme in the scheme comparison (`full_run.sh`, `run -S all`), so the comparison isolates the scheme's effect at a fixed selected *k*-mer length. *k* = 21 gives *w* = *L* - *k* + 1 = 11, the window of the GreedyMini tables `w11_k13`/`w11_k15`. `greedymini` runs `orders/GreedyMini/w11_k15_canonical.gm` (reported in the paper) and `w11_k13_canonical.gm` (additional), extended to the 21-mer around the center `gm_k`-mer (needs `21 - gm_k` even). The effect of *k* alone is a separate question: `experiments/sweep_m.yaml` sweeps it for `native` only.

### Equivalence checking

Two hard-fail gates run at `export` time:

1. **`n_kmers` must be identical** across every scheme/seed run of the same dataset. `build --check`/`sshash check` verify that each index is self-consistent; this is the cross-scheme check that the set of indexed *k*-mers did not change with the scheme.
2. **Fixed positive/negative *k*-mer query answers must be identical.** Once per dataset, `scripts/equivalence.py` builds a throwaway `native` reference index, samples ~10,000 real *k*-mers from the input as positives and ~10,000 random 31-mers verified absent (via `sshash query` against the reference) as negatives, and caches them under `$SSHASH_DATA_DIR/<id>/queries/{pos,neg}.txt`. Every scheme/seed run re-queries both sets against its own index and records `query_pos_match`/`query_neg_match` (all present / all absent, from `sshash query`'s aggregate `query_n_positive_kmers`/`query_n_negative_kmers`; it does not report per-k-mer answers, so the count against an all-positive or all-negative file is the check). `sshash check`'s `EVERYTHING OK!` is part of this gate (`check_ok`).

`export` (`scripts/export.py:validate_schemes_equivalence`) reads these fields from every run's `meta.json` and fails, printing every violation, if any is false or if `n_kmers` disagrees.

The offline sanity test (`tests/test_scheme_sanity.py::test_sshash_integration`) runs the same logic end-to-end on the in-repo sample with the committed `data/queries/salmonella_enterica_k31_{pos,neg}.txt` (3000 *k*-mers each, generated the same way against a `native` build of that sample), so it needs no network and asserts a fixed `n_kmers` (4787534).

## Outputs

On-disk layout:

```
$SSHASH_BENCH_ROOT/
  data/<dataset>/      # eulertigs, query FASTQ, queries/{pos,neg}.txt
  results/runs/        # per-run artifacts (below)
  results/csvs/        # exported tables (sshash_*.csv, expected-density cache)
  results/sshash_results.xlsx   # workbook
```

Per-run artifacts under `$SSHASH_BENCH_ROOT/results/runs/`:

```
runs/
  {dataset}/
    {modality}/            # canonical | regular (-M/--modality)
      {scheme}[/params]/    # native | xor | xor-pcg-hash | mod-minimizer-xor/t10 |
                            # greedymini/k15_at0-cg1 | mod-minimizer-greedymini/k10_at0-cg1_t10
                            # (scheme_names.scheme_dir_name: k<gm_k>_<projection>_t<mod_t>)
        {seed_idx}/          # 0…n-1
          build.log          # sshash build --verbose --bucket-stats stdout+stderr (trailing JSON line)
          build_time.txt     # GNU /usr/bin/time -v
          check.log          # sshash check --verbose (EVERYTHING OK! or diagnostics)
          bench.log          # sshash bench (trailing JSON line)
          query.log          # sshash query vs the dataset's FASTQ (if configured)
          query_time.txt     # GNU /usr/bin/time -v for query.log
          equivalence.json   # {"pos": {...}, "neg": {...}} fixed-set check results
          meta.json          # run identity + config + exit codes + command (L and k are recorded here, not in the path)
          # index.sshash is deleted after check/bench/query/equivalence unless --keep-index
```

`export` (also run at the end of `run`/`all`) writes CSVs to `$SSHASH_BENCH_ROOT/results/csvs/` and the workbook to `$SSHASH_BENCH_ROOT/results/`. It computes the expected density itself (cached in `sshash_density_expected_schemes.csv`) when the cache is missing or does not match.

| file | rows | purpose | in `.xlsx` |
|---|---|---|---|
| `sshash_runs.csv` | one per run | all metrics + telemetry | `runs` |
| `sshash_schemes.csv` | one per scheme variant (*k*, *w*, *L*) | description, CLI, expected + particular density_factor | `schemes` |
| `sshash_summary.csv` | one per dataset × scheme variant | mean/std over seeds, ratio to native | `summary` |
| `sshash_buckets.csv` | one per dataset × scheme variant | bucket metrics (mean/std over seeds, ratio to native) | `buckets` |
| `sshash_space_breakdown.csv` | one per dataset × scheme variant | bits/*k*-mer by component, saving vs native, share of the saving | `space_breakdown` |
| `sshash_bucket_histogram.csv` | one per run × log2 bucket-size bin | bucket-size histogram | `bucket_histogram` |
| `sshash_super_kmer_histogram.csv` | one per run × superwindow length | superwindow-length histogram | `super_kmer_histogram` |
| `significance_sshash.csv` | one per metric × dataset × scheme | paired test of the change from native (`benchmark/common/significance.py`, see `benchmark/README.md`) | `significance` |
| `density_effect_sshash.csv` | one per metric × (all, dataset) | Pearson r of the metric and the density factor, both as ratios to native's run of the same seed (`benchmark/common/density_effect.py`) | `density_effect` (points: `density_effect_points`) |
| `sshash_density_expected_schemes.csv` | one per scheme variant | expected-density cache | — |
| `sshash_results.xlsx` | — | `overview`, `schemes`, `datasets`, `summary`, `runs`, `buckets`, `space_breakdown`, `bucket_histogram`, `super_kmer_histogram`, `density_effect`, `density_effect_points`, `significance` (`benchmark/common/xlsx_book.py`) | — |

Re-export only, from `benchmark/sshash`: `./sshash_bench.sh export` (no network needed).

### CSV columns (main fields)

**Shared keys:** `dataset`, `k` (selected *k*-mer length, SSHash's m), `w`, `L` (= *k* + *w* - 1, SSHash's *k*), the scheme variant `scheme` (`native`/`xor`/`xor-pcg-hash`/`greedymini`/`mod-minimizer-xor`/`mod-minimizer-greedymini`), `gm_k`, `projection` (`at0-cg1`/`at1-cg0`), `gm_table` (table file stem), `mod_t` (blank where they do not apply), then the method's parameters `modality` (`canonical`/`regular`) and, in runs, `threads`, then the seed columns (`n_repeats`, the number of runs in the row's mean/std, in aggregated sheets; `seed`, `seed_idx` in runs). Summary, bucket and space rows are grouped by all of these.

**Space:** `n_kmers`, `index_size_gb` (from `build`'s JSON), `index_disk_gb` (`os.path.getsize`, the ground truth), `bits_per_kmer` (`index_disk_gb * 1e9 * 8 / n_kmers`). All sizes and memory are decimal GB (1 GB = 10^9 bytes, `benchmark/common/units.py`); `meta.json` and the run dirs keep the raw bytes / KiB, converted at export.

**Density:** `density_factor` and `distinct_minimizers` (the number of distinct selected *k*-mers, same value as `n_buckets`), `n_minimizer_positions` (runs).

**Build:** `build_wall_s`/`build_max_rss_gb` (GNU time, KiB -> GB), `build_steps_total_s` and per-step `build_step{1..7}_s` (from `build --verbose`'s microsecond timings).

**Query (primary, `sshash bench`):** `bench_positive_ns`, `bench_negative_ns`, `bench_access_ns`, `bench_iterator_ns` (all `avg_nanosec_per_kmer`).

**Query (secondary, real FASTQ):** `query_wall_s`/`query_max_rss_gb` (GNU time), `query_run_s` (SSHash's own timed query phase, seconds), `query_load_s` (= `query_wall_s - query_run_s`), `query_n_kmers`/`query_n_positive_kmers`/`query_n_negative_kmers`.

**Correctness/equivalence:** `check_ok`, `query_pos_match`, `query_neg_match` (booleans; `export` fails if any is false for any run).

### XLSX sheets

| sheet | content |
|---|---|
| `overview` | column reference: what every non-obvious column means |
| `schemes` | one row per (scheme variant, *k*, *w*, *L*) (*k* = 21, *w* = 11, *L* = 31): scheme, *k*, *w*, *L*, gm_k, projection, gm_table, mod_t, description, cli, expected + per-dataset particular density_factor mean/std (`export` fails on an empty density cell; also written as `sshash_schemes.csv`; layout in `benchmark/common/schemes_sheet.py`) |
| `datasets` | one row per dataset: *k*, *w*, *L*, measured n_kmers |
| `summary` | mean/std per dataset × scheme variant over seeds, plus `{metric}_ratio` = this variant's mean / the same-modality native's mean (native = 1); green (better) → red (worse) color scale on the ratio-scale metrics. Two rows under the table hold, for every `{metric}_ratio` column, its Excel `CORREL` with the density ratio (`density_factor_ratio`) and with the distinct-minimizer ratio (`distinct_minimizers_ratio`) over the dataset × scheme rows. |
| `buckets` | one row per dataset × scheme variant: bucket metrics (bucket sizes, skew-index share, bits/*k*-mer of the tag array and the light array L, …), mean/std over seeds and ratio to native; a row under the table holds the `CORREL` of `n_buckets_ratio` with the tag-array and light-array ratios |
| `space_breakdown` | one row per dataset × scheme variant: where the bits/*k*-mer go (strings, offsets, MPHF, tags, light array L, skew index f_i+V_i and H), each component's saving over native and its share of the total saving, and the saving split into what is sized by the distinct minimizers (MPHF + tags) and what by the occurrences (L + skew index) |
| `bucket_histogram` | one row per run × log2 bucket-size bin (singleton / light / heavy buckets, number of buckets and of *k*-mers in them; from `build --bucket-stats`) |
| `super_kmer_histogram` | one row per run × superwindow length (*k*-mers per superwindow): number of superwindows of that length |
| `runs` | one row per run (every metric + telemetry, unfiltered); after `experiments/profile_run.sh` also the streaming-query profile `query_prof_full_lookups`, `query_prof_positions_per_full_lookup`, `query_prof_verify_calls`, `query_prof_positions_per_verify` (the two ratios also on `summary`) |

## Measurement protocol

Build, bench and query are timed under the protocol in [`../README.md`](../README.md): node-quiet wait (bench and query also wait for free huge-page memory), pinning (build `$THREADS` cores + 2 spare, the single-threaded bench and query 1 + 1; one NUMA node, one thread per core; lists in `cpus.txt`, `cpus_bench.txt`, `cpus_query.txt`, exported as `build_/bench_/query_pinned_cpus`), memory bound to the slot's node and a step record (`step_*.json`, exported as `build|bench|query_node`, `_huge_fraction`, `_foreign_busy_cores`, …; see the protocol), and per-step wall and CPU seconds, RSS, major page faults and reads (runs sheets only). Query wall time splits into `query_load_s` (index load and start-up) and `query_run_s` (SSHash's own timed queries, seconds). `experiments/full_run.sh` needs the default-encoding binary (`sshash/build/sshash`, `SSHASH_BIN` overrides; a build with the traditional nucleotide encoding cannot run GreedyMini).

## Notes

- `greedymini` (`sshash build --scheme greedymini --greedy-mini-file TABLE`, the center extension of the table to the selected *k*-mer, see `orders/README.md`) runs once per variant, a (table, projection) pair from `config.yaml`'s `greedymini` section (`tables`: `w11_k13_canonical.gm`, `w11_k15_canonical.gm`; `projections`: `at0-cg1`, `at1-cg0`). `--gm-tables` / `--gm-projections` replace those lists for one run. `at1-cg0` flips the polarity of Φ (`sshash build --gm-projection at1-cg0`: A,T -> 1, C,G -> 0 instead of A,T -> 0, C,G -> 1). Each variant is its own run directory (`greedymini/k15_at1-cg0/...`), summary/schemes row and `meta.json` entry with its own `gm_k`, `projection` and `gm_table`. The center extension needs `k - gm_k` even, so a table with an even `gm_k` is rejected at *k* = 21; `sshash_bench.sh run` checks this up front (`python -m scripts schemes check-gm greedymini M TABLE`). Example: `./sshash_bench.sh run -d cod,kestrel -n 5 -S native,greedymini --gm-tables w11_k15_canonical --gm-projections at0-cg1`.
- Strand-symmetric mod-minimizers (`sshash build --scheme mod-minimizer-xor` / `--scheme mod-minimizer-greedymini` with `--mod-t t`, see `orders/README.md`): `mod-minimizer-xor` ranks the *L*-mer's *t*-mers by xor, `t = 4 + (k - 4) mod w` (10 at *L* = 31, *k* = 21); `mod-minimizer-greedymini` by a GreedyMini table of *t*-mers from `config.yaml`'s `mod_minimizer_greedymini` section (`w22_k10`, *t* = 10 = the table's `gm_k`). Both are in `all`; `--regular` does not apply to them. With *t* ≡ *k* (mod *w*) (both here) the scheme is forward, which the build checks. The schemes sheet's expected density is counted on uniformly random DNA (equal for both polarities up to sampling noise, by symmetry); the particular density per dataset is what can differ.
- The built `.sshash` index is deleted after each run's check/bench/query/equivalence steps (6 datasets × 10 variants × 5 seeds would otherwise accumulate); pass `--keep-index` to keep it.
- `sshash build`/`bench`/`query --verbose` each print one trailing single-line JSON object (values are JSON strings, parsed with `json.loads` then cast); `lib/run.sh` captures each tool's combined stdout+stderr to its own `.log`, and `scripts/export.py` finds the trailing `{...}` by regex rather than assuming a particular stream.
- Query FASTQ download uses ENA's `vol1/fastq/<prefix6>[/<subfolder>]/<accession>/<file>` path convention (`lib/datasets.sh:ena_url`), including the extra numeric subfolder of 7/8/9-digit run numbers.
