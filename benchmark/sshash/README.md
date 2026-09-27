# SSHash minimizer-order benchmark

Benchmarks [SSHash](https://github.com/jermp/sshash) (a compact k-mer dictionary index) across
its **minimizer order**: `native` (sshash's own baseline), `xor-hash`, `xor-pcg-hash`,
`greedymini` (the GreedyMini rank table integration), and their 4 `mod-minimizer-*` counterparts
(Groot Koerkamp & Pibiri, WABI 2024 — see `../../orders/README.md`). Unlike kmc's `lex`
baseline, sshash has no seed-less order — every order, including `native`, takes a seed, so this
harness runs every order once per seed with no special-casing; the 3 hash-based
`mod-minimizer-*` orders additionally run once per `t` in `{10, 11}` (see "Choice of m" below).

## What is measured

For each dataset, k=31, m=21 (fixed across every order -- see "Choice of m" below), and
order (x `t`, for `mod-minimizer-*`):

- **construction time + peak RSS**: `sshash build` wrapped in `/usr/bin/time -v`
- **correctness**: `sshash check` (must print `EVERYTHING OK!`)
- **query time, primary**: `sshash bench`'s own internal micro-benchmark (`positive lookup`,
  `negative lookup`, `access`, `iterator`, all in `avg_nanosec_per_kmer`) -- this is the same
  headline query-time metric the SSHash paper itself reports
- **query time, secondary**: `sshash query` against the dataset's real paired FASTQ reads,
  also `/usr/bin/time -v`-wrapped
- **space**: the built `.sshash` file's size on disk (`os.path.getsize`, which should match
  `build`'s own reported `index_size_in_bytes` up to a small serialization-estimate delta) and
  derived bits/kmer
- **order equivalence** (hard-fail gate, not just a warning -- see below): every order must
  report the identical `num_kmers` for the same dataset, and must answer a fixed set of
  positive/negative k-mer queries identically

## Dataset

Three whole-genome collections (k=31 eulertigs) from Zenodo record
[17582116](https://zenodo.org/records/17582116), each paired with a real-read query FASTQ from
ENA (species/accession mapping resolved by hand via the ENA API's `scientific_name`, since the
record's own description text doesn't spell it out):

| id | species | eulertigs | query FASTQ | published distinct 31-mers (approx) |
|---|---|---|---|---|
| `cod` | Gadus morhua | `cod.k31.eulertigs.fa.gz` (~177 MB) | `SRR12858649.fastq.gz` (~54 MB) | ~502.5M |
| `kestrel` | Falco tinnunculus | `kestrel.k31.eulertigs.fa.gz` (~342 MB) | `SRR11449743_1.fastq.gz` (~478 MB) | ~1.15B |
| `human` | Homo sapiens | `human.k31.eulertigs.fa.gz` (~883 MB) | `SRR5833294.fastq.gz` (~1.67 GB) | ~2.51B |

Only the record's k=31 eulertigs files are used here (its k=63 files and its three pangenome
collections -- ncbi-virus/se/hprc -- are out of scope for this benchmark). Exact filenames,
sizes, and md5s are in `data/datasets.csv` (committed; download URLs are the Zenodo API's file
endpoints and, for the query FASTQ, ENA's standard `vol1/fastq/...` path convention). `scripts/zenodo.py`
can re-fetch the record's file catalog (`python -m scripts zenodo list-files`) to re-verify or
extend those values later -- routine `download`/`run`/`export` never needs network access to
Zenodo itself beyond the one-time dataset download.

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
| `SSHASH_BENCH_ROOT` | dataset root (`data/<id>/` underneath) and default results parent |
| `SSHASH_DATA_DIR` | optional override of the dataset root |
| `SSHASH_RESULTS_DIR` | optional; default `$SSHASH_BENCH_ROOT/results` |

CLI `--data-dir`/`--results-dir` override the env file. Existing environment variables are not
overwritten by `.env`.

**Git:** the harness (`scripts/`, `lib/`, `data/*.csv`, `data/queries/*.txt`, `experiments/`,
…) is meant to be committed. `benchmark/.gitignore` excludes `.env`, downloaded eulertigs/FASTQ
under `data/<id>/`, and `results/`.

## Layout

```
benchmark/sshash/
  config.yaml         # seeds, paths, defaults (single source of truth)
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
| `config.yaml` | Shared seeds (same 10 hex literals as `benchmark/kmc/config.yaml`, for cross-tool consistency), paths, defaults |
| `sshash_bench.sh` | Main entry (download / run / export / experiment) |
| `lib/*.sh` | Bash helpers: `constants.sh` (shared defaults + `config.yaml` loading), `env.sh` (.env loading), `datasets.sh` (download + md5), `orders.sh` (order -> CLI flags), `run.sh` (build/check/bench/query + equivalence, `/usr/bin/time -v` wrapping), `experiment.sh` (YAML sweep dispatch) |
| `scripts/config.py` | Load `config.yaml` |
| `scripts/orders.py` | Order name -> sshash CLI flags (Python twin of `lib/orders.sh`) |
| `scripts/equivalence.py` | Fixed positive/negative k-mer query-set generation + per-run equivalence checks |
| `scripts/export.py` | CSV / XLSX export + the equivalence hard-fail gate |
| `scripts/zenodo.py` | Zenodo record 17582116 lookup (maintainers only; not on the normal run path) |
| `scripts/meta.py` | Per-run `meta.json` writer |
| `scripts/common.py` | Shared helpers (`.env` loading, `die()`) |
| `scripts/density.py` | Sliding-window "particular" density + pure-Python re-implementations of sshash's order formulas |
| `scripts/expected_density.py` | Monte Carlo theoretical `density_factor` per order (feeds the `orders` xlsx sheet) |
| `scripts/experiments.py` | Load/expand `experiments/*.yaml` (shared engine in `benchmark/common/experiments_lib.py`) |
| `data/` | Dataset catalog (`datasets.csv`) + committed sanity-test query sets (`queries/`) |
| `experiments/*.yaml` | Multi-config sweep definitions |
| `tests/` | Offline sanity checks: `test_order_sanity.py` (order/equivalence correctness against the in-repo sample) + `test_export.py` (`export.py`'s aggregation/validation/CSV/XLSX logic, synthetic rows, no sshash binary needed); `make venv` once, then `make test` |

Python CLI (from the repo root, `BENCHMARK_DIR=benchmark/sshash` set by `sshash_bench.sh`):

```bash
python3 -m scripts export --results-dir …
python3 -m scripts experiment list
python3 -m scripts orders flags xor-hash 0xDEADBEEF
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
benchmark/sshash/sshash_bench.sh run -d cod -n 1     # quick smoke: 1 seed, all 9 orders
# → $SSHASH_BENCH_ROOT/results/sshash_{runs,summary,vs_reference}.csv
# → $SSHASH_BENCH_ROOT/results/sshash_results.xlsx
```

### `sshash_bench.sh` (main entry)

```bash
benchmark/sshash/sshash_bench.sh -h
benchmark/sshash/sshash_bench.sh download -d cod
benchmark/sshash/sshash_bench.sh run -d cod -n 1                     # 1 seed, all orders
benchmark/sshash/sshash_bench.sh run -d cod,kestrel,human -O all -n 5 -t 16   # full sweep
benchmark/sshash/sshash_bench.sh export                              # re-export CSV/XLSX only
benchmark/sshash/sshash_bench.sh all -d cod -n 1                      # download + run + export
```

Defaults: `cod,kestrel,human`, k=31, m=21, `-O all` (native, xor-hash, xor-pcg-hash, greedymini,
mod-minimizer-native-hash, mod-minimizer-xor-hash, mod-minimizer-pcg-xor-hash,
mod-minimizer-greedymini, mod-minimizer-greedymini-v2), `-M canonical`, `-n 5` seeds, 8 threads.

The command name (`download`/`run`/`export`/`all`/`experiment`) must come **first**, before any
options -- e.g. `sshash_bench.sh run --dataset-file PATH -n 1`, not the reverse.

`-M/--modality` selects `canonical` (default, single-probe) or `regular` (two-probe, `sshash
build --regular`); a `mod-minimizer-*` order is always `regular` regardless of this flag, so e.g.
`-O native -M regular` is the apples-to-apples comparison against mod-minimizer's own (mandatory)
two-probe path. `--force` re-runs a config even if its `meta.json` already exists (default:
skip unchanged configs).

**Experiments** (multi-config sweeps from YAML):

```bash
benchmark/sshash/sshash_bench.sh experiment list
benchmark/sshash/sshash_bench.sh experiment sweep_orders   # primary: all 9 orders (x t for mod-minimizer-*) x 5 seeds, all 3 datasets
benchmark/sshash/sshash_bench.sh experiment sweep_m        # secondary/optional: native only, m in {7,9,11,13}
```

Experiment files live in `experiments/`. Each has shared `defaults` plus a `configs` list; keys
are `datasets`, `k`, `m`, `orders`, `gm_file`, `n_seeds`, `threads`, `force`, `run_tag`.

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

**m is fixed at 21 for every order** in the primary comparison (`sweep_orders.yaml` /
`run -O all`), by deliberate design: the point of this benchmark is isolating the *order's*
effect on space/time at a fixed minimizer length, not conflating it with an m sweep. m=21 (up
from an earlier m=11) matches k=31/m=21/w=11 — the mod-minimizer paper's own tested SSHash
regime (Groot Koerkamp & Pibiri, WABI 2024): mod-minimizer needs m large relative to w to show
its benefit, and using the same (k, m) for every order keeps the comparison apples-to-apples.
`greedymini` uses `orders/GreedyMini/w11_k15_canonical.gm` (`gm_k=15`, Extend-k'd up to `l=21`);
`mod-minimizer-greedymini`'s O_t table is `orders/GreedyMini/w21_k11_canonical.gm` at `t=11`
(`gm_k=11=t`, no extend needed there) and `w22_k10_canonical.gm` at `t=10` (`gm_k=10=t`). Varying
m *is* useful as its own, separate question -- `experiments/sweep_m.yaml` covers that, explicitly
at a single order (`native`) so it isn't mistaken for a second order comparison, and is
optional/secondary to the primary `sweep_orders` results.

The 3 hash-based `mod-minimizer-*` orders, and `mod-minimizer-greedymini` itself, sweep `t` (the
length of the inner order O_t): `t=10` reproduces the mod-minimizer paper's own published
Human/GRCh38 numbers (its own formula, Definition 16, gives `t=10` for this regime); `t=11`
instead matches `w21_k11_canonical.gm`'s own length exactly (the paper never used GreedyMini as
O_t, so there's no "official" t for it either way). `mod-minimizer-greedymini-v2` is a separate
order: an A/B comparison against an alternate/updated `t=11` GreedyMini table
(`w21_k11_canonical_v2.gm`), fixed at `t=11`. Both GreedyMini tables' `t=10` case and any further
labeled table variant are looked up from `config.yaml`'s `mod_t_gm_file_t10`/`mod_t_gm_labels`
keys -- adding another labeled variant later (`mod-minimizer-greedymini-v3`, say) needs one new
config.yaml line, not a code change (see `lib/orders.sh::order_gm_file`). See
`../../orders/README.md#mod-minimizer`.

### Equivalence checking

Two independent hard-fail gates run at `export` time (mirroring kmc's
`validate_custom_vs_lex`, generalized here since sshash has no seed-less baseline to compare
everything else *against* -- every order, including `native`, is just another group member):

1. **`num_kmers` must be identical** across every order/seed run of the same dataset. sshash
   already asserts this internally in a sense (`build --check`/`sshash check` verify the index
   is self-consistent), but this is the cross-order check that the *set of k-mers indexed*
   didn't silently change because of the order.
2. **Fixed positive/negative k-mer query answers must be identical.** Once per dataset,
   `scripts/equivalence.py` builds a throwaway `native`-order reference index, samples ~10,000
   real k-mers from the input as positives and ~10,000 random 31-mers verified absent (via
   `sshash query` against that reference index) as negatives, and caches them under
   `$SSHASH_DATA_DIR/data/<id>/queries/{pos,neg}.txt`. Every subsequent order/seed run re-queries
   both fixed sets against its own index and records `query_pos_match`/`query_neg_match`
   booleans (all-present / all-absent, from `sshash query`'s own aggregate
   `num_positive_kmers`/`num_negative_kmers` counts -- it doesn't report per-k-mer answers, so
   the aggregate count against a fixed all-positive or all-negative file is the check).
   `sshash check`'s own `EVERYTHING OK!` is also folded into this gate (`check_ok`).

`export` (`scripts/export.py:validate_orders_equivalence`) reads these fields back from every
run's `meta.json` and `die()`s hard, printing every violation, if any is false or if
`num_kmers` disagrees -- not a warning, a failed export.

The offline sanity test (`tests/test_order_sanity.py::test_sshash_integration`) exercises the
same logic end-to-end against the in-repo sample dataset, using the **committed**
`data/queries/salmonella_enterica_k31_{pos,neg}.txt` (3000 k-mers each, generated the same way
against a native-order build of that sample) instead of generating a fresh set, so the test has
no network dependency and a fixed expected `num_kmers` (4787534) to assert against.

## Outputs

Per-run artifacts live under `$SSHASH_BENCH_ROOT/results/runs/` in a nested layout:

```
runs/
  {dataset}/
    {modality}/            # canonical | regular -- see -M/--modality
      {order}/              # native | xor_hash | xor_pcg_hash | greedymini |
                             # mod_minimizer_{native_hash,xor_hash,pcg_xor_hash,greedymini} |
                             # mod_minimizer_greedymini_<label> (e.g. _v2)
        [t{t}/]              # mod-minimizer-* only: t10 or t11 (fixed t11 for a labeled variant)
        {seed_idx}/          # 0…n-1
          build.log          # sshash build --verbose stdout+stderr (trailing JSON blob)
          build_time.txt     # GNU /usr/bin/time -v
          check.log          # sshash check --verbose (EVERYTHING OK! or diagnostics)
          bench.log          # sshash bench (trailing JSON blob)
          query.log          # sshash query vs the real dataset FASTQ (if configured)
          query_time.txt     # GNU /usr/bin/time -v for query.log
          equivalence.json   # {"pos": {...}, "neg": {...}} fixed-set check results
          meta.json          # run identity + config + exit codes + command (k/m recorded here, not in the path)
          # index.sshash is deleted after check/bench/query/equivalence unless --keep-index
```

`export` (also run at the end of `run`/`all`) writes CSVs and a workbook to
`$SSHASH_BENCH_ROOT/results/`:

| file | rows | purpose | in `.xlsx` |
|---|---|---|---|
| `sshash_runs.csv` | one per run | all metrics + telemetry | `runs` |
| `sshash_summary.csv` | one per dataset x order | mean/std across seeds, ratio vs native | `summary` |
| `sshash_vs_reference.csv` | one per dataset | informational: measured `num_kmers` vs the Zenodo record's own published distinct-31-mer count (not order-specific) | `vs_reference` |
| `sshash_results.xlsx` | — | `index`, `orders`, `datasets`, `summary`, `runs`, `vs_reference` | — |

Re-export only: `benchmark/sshash/sshash_bench.sh export` (no network needed once local run
data exists).

### CSV columns (main fields)

**Shared keys:** `dataset`, `k`, `m`, `w` (`k-m+1`), `modality` (`canonical`/`regular`, see
`-M/--modality`), `scheme` (`minimizer`/`mod-minimizer`) + `order` (the *normalized* order name,
identical across a scheme's base and `mod-minimizer-*` variants: `native-hash`/`xor-hash`/
`pcg-xor-hash`/`greedymini` -- e.g. `-O native` and `-O mod-minimizer-native-hash` both write
`order=native-hash`, distinguished only by `scheme`), `order_seed`, `seed_idx`, `t` (blank except
for `scheme=mod-minimizer`), `gm_file`, `threads`.

**Space:** `num_kmers`, `index_size_in_bytes` (from `build`'s own JSON), `index_bytes_on_disk`
(`os.path.getsize`, the ground truth), `bits_per_kmer` (`index_bytes_on_disk * 8 / num_kmers`).

**Build:** `build_wall_s`/`build_max_rss_kb` (GNU time), `build_total_build_time_s` and
per-step `build_step{1..7}_s` (from `build --verbose`'s own microsecond timings).

**Query (primary, `sshash bench`):** `bench_positive_ns`, `bench_negative_ns`,
`bench_access_ns`, `bench_iterator_ns` (all `avg_nanosec_per_kmer`).

**Query (secondary, real FASTQ):** `query_wall_s`/`query_max_rss_kb` (GNU time),
`query_elapsed_ms` (sshash's own), `query_num_kmers`/`query_num_positive_kmers`/`query_num_negative_kmers`.

**Correctness/equivalence:** `check_ok`, `query_pos_match`, `query_neg_match` (all booleans;
`export` hard-fails if any is false for any run).

### XLSX sheets

| sheet | content |
|---|---|
| `index` | column reference: what every non-obvious column means |
| `orders` | one row per (order, t): what it is, its sshash CLI flags, particular + expected density |
| `datasets` | one row per dataset: k, m, measured num_kmers |
| `summary` | mean/std per dataset x order across seeds, plus `{metric}_ratio` = this order's mean / the same-modality native's mean (native = 1); green (better) → red (worse) color scale on the ratio-scale metrics. Since every `mod-minimizer-*` order is unconditionally `regular` modality, its ratio needs a `native`-under-`regular` run to divide by -- `run_schemes` (`lib/run.sh`) builds one automatically per dataset/seed whenever a sweep includes any `mod-minimizer-*` order, regardless of `-M`, so these columns are never silently blank |
| `runs` | one row per run (every metric + telemetry, unfiltered) |
| `vs_reference` | informational: measured vs the Zenodo record's own published k-mer counts |

## Notes

- Every order (including `native`) takes `-s SEED`; there is no seed-less baseline the way
  kmc's `lex` is, so this benchmark never special-cases a "run once" order the way kmc does.
- `greedymini`'s `--greedy-mini-file` defaults to `orders/GreedyMini/w11_k15_canonical.gm`
  (`gm_k=15`, Extend-k'd up to `-m 21`). Override with `--greedy-mini-file PATH` or the
  `gm_file` experiment-YAML key. `mod-minimizer-greedymini`'s own O_t tables are fixed to
  `config.yaml`'s `mod_t_gm_file`/`mod_t_gm_file_t10` (not overridable from the CLI, to avoid
  ambiguity when several such orders run in the same sweep); a labeled table variant (e.g.
  `mod-minimizer-greedymini-v2`) instead reads `config.yaml`'s top-level `mod_t_gm_labels` map --
  add a new label there to register another variant, no code change needed.
- The built `.sshash` index is deleted after each run's check/bench/query/equivalence steps
  finish (to keep disk usage bounded across cod/kestrel/human x 9 orders x N seeds x up to 2
  `t` values); pass `--keep-index` to retain it.
- `sshash build`/`bench`/`query --verbose` all print one trailing single-line JSON blob (values
  are JSON *strings*, parsed with `json.loads` then cast) -- `lib/run.sh` captures each tool's
  combined stdout+stderr to its own `.log` file, and `scripts/export.py` finds that trailing
  `{...}` via regex rather than assuming it's the literal last line of a particular stream (it
  was observed on stderr in hand-testing, not stdout).
- Query FASTQ download uses ENA's standard `vol1/fastq/<prefix6>[/<subfolder>]/<accession>/<file>`
  path convention, generalized (`lib/datasets.sh:ena_url`) beyond `benchmark/kmc`'s version to
  handle accessions needing the extra numeric subfolder (7/8/9-digit run numbers).
