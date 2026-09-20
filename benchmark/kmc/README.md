# KMC 3 counting benchmark (KMC only)

Reproduces the **k-mer counting tables** from Kokot, Długosz & Deorowicz, *KMC 3: counting and manipulating k-mer statistics*, Bioinformatics 33(17):2759–2761 (2017). Competitors are not run.

Dataset URLs, flags, and published KMC 3 times/RAM/disk are taken from the paper supplement ([arXiv:1701.08022](https://arxiv.org/abs/1701.08022) `supp.pdf`).

## What is measured

For each dataset and \(k \in \{28,55\}\) (optional 40 and 65 from the supplement):

- wall time (KMC JSON `Total`, plus GNU `time`)
- peak RAM (GNU `/usr/bin/time -v` max RSS)
- temp disk (KMC JSON `Tmp_size`)
- unique / counted k-mer stats

Paper command (KMC 3.0.0):

```text
./kmc3 -ci2 -v -t12 -k<k> <fastq|@list> res ./bins
```

H. sapiens 3 used `-m32`. Other sets used the default 12 GB (`-m12`). Gzip input used `@fv.txt` / `@gg.txt` / … file lists.

## Disk (order of magnitude, gzip)

| id | organism | gzip FASTQ |
|---|---|---|
| `f_vesca` | F. vesca | ~3 GB |
| `g_gallus` | G. gallus | ~20–26 GB |
| `m_balbisiana` | M. balbisiana | ~49 GB |
| `h_sapiens_2` | ERA015743 ERR024163–186 | ~106 GB |
| `h_sapiens_3` | ERR174324–341 | ~600 GB |

Uncompressed copies are several times larger. Prefer `--input gzip` first.

Paths come from `benchmark/.env` (see `.env.example`). Copy the example if `.env` is missing:

```bash
cp benchmark/.env.example benchmark/.env
# edit KMC3_BENCH_ROOT — use a directory outside the repo for large downloads
```

| variable | meaning |
|---|---|
| `KMC3_BENCH_ROOT` | FASTQ root (`data/<id>/` underneath) and default results parent |
| `KMC3_DATA_DIR` | optional override of the FASTQ root |
| `KMC3_RESULTS_DIR` | optional; default `$KMC3_BENCH_ROOT/results` |

CLI `--data-dir` / `--results-dir` override the env file. Existing environment variables are not overwritten by `.env`.

**Git:** the harness (`scripts/`, `lib/`, `data/*.csv`, `experiments/`, …) is meant to be committed. `.gitignore` excludes `.env`, downloaded FASTQ under `data/*/`, generated `.gm` tables, and `results/`.

## Usage

Requires **Python 3** for CSV/XLSX export and experiment YAML. **Bash** handles download, prepare, and kmc runs.

From the repo root:

```bash
python3 -m venv .venv
source .venv/bin/activate
make -j
benchmark/kmc3_bench.sh download -d f_vesca
benchmark/kmc3_bench.sh run -d f_vesca -k 28,55 --kmc bin/kmc -t 12
# → $KMC3_BENCH_ROOT/results/kmc3_{runs,summary,vs_paper}.csv
# → $KMC3_BENCH_ROOT/results/kmc3_results.xlsx  (summary + runs + vs_paper sheets)
```

### Layout

```
benchmark/
  config.yaml         # seeds, paths, defaults (single source of truth)
  kmc3_bench.sh       # main entry (download / run / export / experiment)
  run_config_sweep.sh
  lib/*.sh            # bash helpers
  scripts/            # python -m scripts …
  experiments/*.yaml  # multi-config sweep definitions
  minimizers/         # GreedyMini rank tables (.txt, .gm)
  data/               # datasets.csv, paper reference (catalog)
  tests/
  README.md
```

| Path | Role |
|---|---|
| `config.yaml` | Shared seeds, paths, defaults (loaded by bash + Python) |
| `kmc3_bench.sh` | Main entry (download / run / export / experiment) |
| `lib/*.sh` | Bash helpers (env, python, gm, datasets, runs, experiments) |
| `scripts/` | Python package (`python -m scripts …`) |
| `scripts/config.py` | Load `config.yaml`; `config shell-env` for bash |
| `scripts/export.py` | CSV / XLSX export |
| `scripts/experiments.py` | Load / expand experiment YAML |
| `scripts/gm.py` | GreedyMini path helpers + table build/ensure (shells out to the C++ tools below) |
| `tools/gm_sigma_extend.cpp` | Binary `.gm` → DNA `_dna.gm` (C++, built via `tools/Makefile`) |
| `scripts/meta.py` | Per-run `meta.json` writer |
| `data/` | Dataset catalog (`datasets.csv`, `paper_kmc3.csv`) |
| `experiments/*.yaml` | Multi-config sweep definitions |
| `../../orders/GreedyMini/` | Rank lists (`.txt`, shared with kraken2) and tables (`.gm`, `_dna.gm`) |
| `../../common/minimizer-orders/gm_serializer.cpp` | Rank list (`.txt`) → binary `.gm` (C++, shared with kraken2, built via that directory's `Makefile`) |
| `tests/` | Sanity checks — order/density correctness, not report formatting; `make venv` once, then `make test` (pytest) |
| `run_config_sweep.sh` | Wrapper → `experiment sweep_config` |

Python CLI (from repo root, `BENCHMARK_DIR=benchmark` set by kmc3_bench.sh):

```bash
python3 -m scripts export --results-dir …
python3 -m scripts experiment list
python3 -m scripts gm ensure --k 28 -p 9 --seed 0xDEADBEEF
```

### kmc3_bench.sh (main entry)

```bash
benchmark/kmc3_bench.sh -h
benchmark/kmc3_bench.sh download -d f_vesca
benchmark/kmc3_bench.sh run -d f_vesca -k 28 -n 1      # quick smoke (all schemes, 1 rep)
benchmark/kmc3_bench.sh run -O all -n 10 --force       # full multi-seed benchmark
benchmark/kmc3_bench.sh export                           # re-export CSV/XLSX only
benchmark/kmc3_bench.sh all -d f_vesca -k 28 -n 1       # download + run + export
```

Defaults: all five datasets, `k=28,55`, gzip input, 12 threads, `-O all`, `-n 1` repeat of each scheme
(including lex — see "Run order and the lex baseline" below).

**Experiments** (multi-config sweeps from YAML; KMC allows `-n` bins in **64..2000**):

```bash
benchmark/kmc3_bench.sh experiment list
benchmark/kmc3_bench.sh experiment sweep_config                 # all datasets, -O all, 5 seeds, 8 configs
benchmark/kmc3_bench.sh experiment sweep_gm_lut_fixed_dataset   # f_vesca, GM-DNA vs GM-bin, 6 layouts
benchmark/kmc3_bench.sh experiment sweep_gm_lut_fixed_config    # all datasets, GM-DNA vs GM-bin, t=24 n=1024
benchmark/run_config_sweep.sh                                   # same as sweep_config
```

`sweep_config` runs **800 KMC jobs** (5 datasets × 8 configs × [4 schemes × 5 repeats], lex repeated the same as the custom schemes — see "Run order and the lex baseline") at k=28 gzip only. Smoke-test with `run -d f_vesca -k 28 -n 1`. Full sweep including `h_sapiens_3` (~600 GB gzip) needs multi‑TB disk and days of wall time.

GM LUT experiments compare the native DNA table (`_dna.gm`) vs the binary rank LUT (`gm_lut_mode: dna|binary` + `run_tag`); the binary LUT is kmc's default backend for a binary `.gm` source, use `-ode` to opt into the DNA-expanded table instead:

| experiment | scope | KMC jobs | ~wall |
|---|---|---|---|
| `sweep_gm_lut_fixed_dataset` | f_vesca, bins 512/1024/1536 × t12/t24 | 12 | ~4 min |
| `sweep_gm_lut_fixed_config` | five datasets, t=24 bins=1024 | 10 | ~1.5 h |

Experiment files live in `experiments/`. Each file has shared `defaults` plus a `configs` list.
Per-config keys: `threads`, `bins`, `mem`, `strict_mem`, `datasets`, `k`, `input`,
`signature_order`, `order_label`, `n_seeds`, `force`, `gm_lut_mode`, `run_tag`, …

Example (`experiments/sweep_config.yaml` excerpt):

```yaml
# sweep_config.yaml (excerpt)
defaults:
  datasets: f_vesca,g_gallus,m_balbisiana,h_sapiens_2,h_sapiens_3
  k: 28
  input: gzip
  signature_order: all
  n_seeds: 5
  force: true

configs:
  - threads: 12
    bins: 512
  - threads: 24
    bins: 1536
    mem: 8
    strict_mem: true
  # … remaining 512/1024/1536 × threads 12/24
```

Run directories are nested so dataset / params / config / scheme / seed do not overwrite each other.
Use `run_tag` in YAML or `--run-tag` to replace the auto-generated **config** segment.

Custom signature order (GreedyMini table, random hash, or random PCG):

GreedyMini UHS rank lists live in `../../orders/GreedyMini/*.txt` (shared
with kraken2). Build tools:

| Tool | Role |
|---|---|
| `../../common/minimizer-orders/gm_serializer` (C++, shared with kraken2) | Rank list (`.txt`) → binary `.gm` |
| `tools/gm_sigma_extend` (C++) | Binary `.gm` → DNA `_dna.gm` |
| `scripts/gm.py` | Path helpers + `ensure` (auto-build on benchmark runs; shells out to the tools above, doesn't reimplement them) |

The pipeline **auto-builds** tables on GM runs:

1. `{stem}.txt` → `{stem}.gm` (binary ranks, seed-independent)
2. `{stem}.gm` → `{stem}_dna.gm` or `{stem}_dna_s0x{seed}.gm` (σ-extended LUTs)

The `_dna.gm` file **embeds the tie-break seed** in its header (magic `KMCANDM1`); sort/sig tables are
baked for that seed. KMC reads the seed from the file — **do not pass `-os`** with a pre-built `_dna.gm`
(the benchmark omits it automatically; `meta.json` records the embedded seed). Pass `-os` only when
supplying a raw binary `.gm` via `--custom-order` (unusual; the default pipeline always builds `_dna.gm`).
For `-n > 1` GM seeds, the pipeline builds one `_dna.gm` per seed. Re-run after editing `.txt` is automatic (mtime check).

Manual build (optional):

```bash
python3 -m scripts gm ensure --k 28 -p 9 --seed 0xDEADBEEF

# or step-by-step, using the C++ tools directly:
make -C ../common/minimizer-orders
../common/minimizer-orders/gm_serializer --canonical -k 9 \
    ../orders/GreedyMini/w20_k9_canonical.txt \
    -o ../orders/GreedyMini/w20_k9_canonical.gm

make -C tools
tools/gm_sigma_extend \
    ../orders/GreedyMini/w20_k9_canonical.gm \
    --seed 0xDEADBEEF \
    -o ../orders/GreedyMini/w20_k9_canonical_dna.gm

# All schemes + N seeds (default lex + hash + pcg + gm)
benchmark/kmc3_bench.sh \
  -d f_vesca,g_gallus,m_balbisiana,h_sapiens_2,h_sapiens_3 \
  --kmc bin/kmc -k 28,55 -t 12 \
  -O all -n 10
```

`-O all` runs lex, hash (`rand_hash`), pcg (`rand_pcg`), and gm (`gm_w{k-p+1}_k{p}`), each custom scheme with `-n` seeds from `config.yaml`. Optional `--custom-order` overrides the per-k `.gm` path.

Predefined seeds are in `config.yaml` (10 fixed hex values). Custom runs use `seed_idx` `0…n-1` in the results tree; `signature_order_label` in exports still includes the hex seed (e.g. `rand_hash_s0xDEADBEEF`).

### Run order and the lex baseline

`-O all` **interleaves** lex/hash/pcg/gm within each repeat `i` (lex, then hash, then pcg, then
gm, then repeat for `i+1`) instead of running all of one scheme back-to-back before starting the
next. It also repeats lex `-n` times, the same as the custom schemes, instead of running it once.

Both changes fix the same problem: with lex run once, first, and each custom scheme run as an
uninterrupted 5-seed block afterward, any time-of-day or machine-state drift over the course of a
multi-hour run lands as a step change *between* schemes — which reads exactly like a real
ordering effect but isn't one. Interleaving spreads that drift evenly across all four schemes
instead of loading it onto lex's single sample, and repeating lex gives it a real mean/std instead
of one unrepeatable data point. `-O default -n N` also accepts `N > 1` now, for the same reason,
when running lex on its own.

This changes `n_seeds` and adds `{metric}_std` for `order_label == default` rows in every
export (previously always `n_seeds=1`, `std=0`) — `kmc3_vs_paper.csv` and the `runs`/`summary`
sheets need no code changes for this, the existing per-config aggregation already treats lex like
any other scheme once it has real repeats.

**`{metric}_z` (summary/gm_lut sheets):** alongside `{metric}_ratio` (`= this / lex`), the summary
and gm_lut xlsx sheets now carry `wall_time_z`, `stage1_z`, `stage2_z` — Welch's t-statistic of
`(this mean − lex mean)` against both sides' pooled seed/repeat std. `|z| ≥ 2` is a reasonable
"probably a real difference, not seed noise" cutoff (highlighted red/green in the sheet); a ratio
far from 1 with `|z| < 2` is more likely noise, especially on short/fast configs where a handful of
seeds gives a poor std estimate. The cell is blank (`IFERROR`) when either side has too few repeats
(e.g. `n_seeds=1`) to estimate a std from.

All five datasets, gzip then uncompressed:

```bash
benchmark/kmc3_bench.sh all \
  --kmc bin/kmc \
  -k 28,40,55,65 \
  --input gzip,plain \
  -t 12
```

`prepare` (called from `all` when `plain` is requested) gunzips into `$KMC3_BENCH_ROOT/data/<id>/fastq_plain/`. Delete those files after the plain runs if you need space.

## Outputs

Per-run artifacts live under `$KMC3_BENCH_ROOT/results/runs/` in a nested layout:

```
runs/
  {dataset}/                 # e.g. f_vesca
    k{k}_p{p}/               # e.g. k28_p9  (-k and -p / signature length)
      {config}/              # e.g. gzip_n512_t12_m12  or  gzip_n1536_t24_m8_sm
        {scheme}/            # default | lexicographic | hash | pcg | greedymini
          {seed_idx}/        # 0…n-1 for every scheme, including lex (-n repeats)
            kmc.json         # KMC -j summary
            time.txt         # GNU /usr/bin/time -v
            meta.json        # run config + command
            out.kmc_*        # sidecars (removed unless --keep-db)
```

**Config segment** (unless `--run-tag`): `{input}_n{bins}_t{threads}_m{mem}` with `_sm` when `--strict-mem` (e.g. `gzip_n512_t12_m12`, `gzip_n1536_t24_m8_sm`). This string is also `config_tag` in CSV exports.

**run_id** in CSVs is the relative path under `runs/`, e.g. `f_vesca/k28_p9/gzip_n512_t12_m12/hash/2`.

`export` (also at the end of `run` / `all`) writes CSVs and a workbook to `$KMC3_BENCH_ROOT/results/`:

| file | rows | purpose | in `.xlsx` |
|---|---|---|---|
| `kmc3_runs.csv` | one per run | Metrics + telemetry | `runs` |
| `kmc3_summary.csv` | one per run (per seed) | Raw headline metrics | — |
| `kmc3_vs_paper.csv` | one per lex run | vs `data/paper_kmc3.csv` (`wall_time_ratio` = ours ÷ paper) | `vs_paper` |
| `kmc3_density_expected_orders.csv` | one per scheme | theoretical expected density, any order (see `python -m scripts expected-density` below) | `orders` |
| `kmc3_results.xlsx` | — | `index`, `orders`, `datasets`, `summary`, `readlen_check`, `gm_lut`, `runs`, `vs_paper` | — |

Re-export only: `benchmark/kmc3_bench.sh export` (or `python3 -m scripts export --results-dir … --paper-platform xeon_hdd|xeon_ssd`).

`kmc3_vs_paper.csv` is omitted if `paper_kmc3.csv` is missing. Custom vs lex ratios live on the `summary` sheet (`{metric}_ratio` columns); export still **fails** if custom k-mer counts differ from lex.

### CSV columns (main fields)

**Shared keys** (all result CSVs): `dataset`, `k`, `input_format` (`gzip` \| `plain`), `config_tag` (e.g. `gzip_n512_t12_m12`, `gzip_n1536_t24_m8_sm`), `n_bins`, `threads`, `strict_memory`, `mem_limit_gb`.

**Schemes:** `order_label` is one of `default`, `lexicographic`, `hash`, `pcg`, `greedymini_dna` / `greedymini_binary` (from `-O default|lex|hash|pcg|gm`). Custom runs also have `signature_order_label` (e.g. `rand_hash_s0xDEADBEEF`, `gm_w20_k9_canonical_s0xdeadbeef`).

| file | additional columns |
|---|---|
| `kmc3_runs.csv` | `run_id`, `label`, `seed_idx`, metrics, telemetry (`user_s`, `cpu_pct`, `host`, …) |
| `kmc3_summary.csv` | per-run `wall_time_s`, `ram_gb`, … (same as one row per seed) |

Both CSVs carry every dataset unfiltered, including read-length-check derived ones — only the xlsx `summary` sheet (below) excludes them in favor of its own `readlen_check` sheet.

### XLSX sheets

| sheet | content |
|---|---|
| `index` | metrics reference: what every non-obvious column means, its units, and how to read it (first sheet; formerly named `legend`, then briefly `metrics`) |
| `orders` | one row per signature-order scheme: what it is, whether it has the `is_allowed` denylist, a small-sample offline particular-density diagnostic (GreedyMini rows only), theoretical `expected_density_avg`/`expected_density_std` on i.i.d. random DNA (any order, not just GreedyMini -- mean/std across this benchmark's predefined seeds for hash/pcg/greedymini; blank std for the deterministic default/lexicographic), and real measured `density_factor` per dataset |
| `datasets` | one row per dataset: size / read-length / config stats (from its own `default` run) plus published KMC3-paper reference numbers where available |
| `summary` | mean / std / ratio / z formulas; charts below (paper datasets only, see `readlen_check`) |
| `readlen_check` | read-length-confound-check derived datasets (`*_shortsub`/`*_longsub`, …), aggregated like `summary` but kept out of it since they aren't paper datasets — present only when such runs exist |
| `gm_lut` | GreedyMini DNA-vs-binary LUT sweep, aggregated like `summary` |
| `runs` | one row per seed (metrics + telemetry) — every dataset, unfiltered; rows grouped in the same dataset → bins → threads → mem → strict_mem → scheme order as `summary` (seed_idx orders repeats within a scheme block), and its "id / params" columns share `summary`'s field order as far as per-run vs per-config granularity allows |
| `vs_paper` | default (lex) runs vs published KMC3 numbers |

### XLSX `summary` sheet (aggregated)

Row 1 is a one-off banner (not part of the table): `Correlation (density_factor_ratio vs disk_ratio):` and a live Excel `=CORREL(...)` formula over the whole sheet — a single sheet-wide statistic, shown once instead of repeated down every row. The category-band row, column headers, and data table sit below it (rows 2/3/4+), and rows 1-3 are frozen together.

One row per **dataset / k / input / config / scheme**, sorted by **dataset → bins → threads → memory → scheme**, with metrics aggregated across seeds. `{metric}_ratio` for wall time, RAM, and disk use a green (better) → yellow (≈1) → red (worse) color scale.

| column pattern | meaning |
|---|---|
| `n_seeds` | number of seed runs in the group |
| `{metric}_mean`, `{metric}_std` | mean and sample std across seeds (`wall_time_s`, `ram_gb`, `disk_gb`, stages, `n_super_kmers`, …) |
| `{metric}_ratio` | Excel formula `=mean / lex_mean` (lex row = **1**); vs the `default` (lex) baseline for the same **dataset / k / input / config sweep** |

### `python -m scripts expected-density` (theoretical density, any order)

Monte Carlo simulation of density_factor on i.i.d. uniform-random DNA (no real sequence content, no read boundaries) for **every** signature-order scheme — `default`, `lexicographic`, `hash`, `pcg`, `greedymini_dna`/`_binary` — generalized beyond the original GreedyMini-only `scripts density expected`. hash/pcg/greedymini are seed-dependent, so this reports mean **and** std across this benchmark's predefined seeds (`default`/`lexicographic` are deterministic — no seed, std blank). Feeds the `orders` sheet's `expected_density_avg`/`expected_density_std` columns.

```bash
python3 -m scripts expected-density                    # k=28 p=9 (this benchmark's defaults)
python3 -m scripts expected-density --force             # recompute even if a matching cache exists
```

Writes `kmc3_density_expected_orders.csv` to the results dir. Like `kmc3_density_particular.csv`, this is a **separate, explicit step** — `export`/`orders_rows()` only *read* the cache (warning to stderr if missing), they never trigger the (one-time, ~tens-of-seconds) computation themselves, so routine `run`/`export` stays fast.

**Charts** (below the summary table): one row per **dataset**, **ratio vs lex** per config sweep — x = scheme (hash / pcg / GM), **legend = config**. Bar height = mean paired-seed ratio vs lex for that dataset and sweep; error bars = std across seeds. One chart per metric (wall time, RAM, disk).

| file | additional columns |
|---|---|
| `kmc3_vs_paper.csv` | `paper_platform`, `n_seeds`, `wall_time_s`, `wall_time_std`, `paper_time_s`, `wall_time_ratio`, `stage1_s`, `stage2_s`, `ram_gb`, `paper_ram_gb`, `disk_gb`, `paper_disk_gb`, `n_unique_counted` |

One row per **lex config** (bins/threads/memory), averaged across its `n_seeds` repeats — not one row per repeat.

**Units:** `ram_gb` is max RSS in **GiB** (KiB ÷ 1024²). `disk_gb` is KMC temp size in **decimal GB** (bytes ÷ 1000³, matching KMC JSON `Tmp_size` MB × 10⁶). Paper reference RAM/disk in `paper_kmc3.csv` uses the same GB convention as the supplement tables.

**Validation:** `export` **fails** if a custom-order run’s unique / counted / total k-mer counts differ from the matching lex run for the same dataset, k, input, and config. Super-k-mer counts may differ. Custom runs without a lex baseline for that key also fail.

Times will not match 2017 hardware; use `wall_time_ratio` in `kmc3_vs_paper.csv`. Gzip vs plain **unique counted k-mers** for the same k should match (export warns on mismatch).

## Read-length confound check

GreedyMini's wall-time win only reproduced reliably on `f_vesca`, which is both the smallest
dataset (~4.1B k-mers) *and* the one with by far the longest reads (avg ~350bp vs ~64-101bp for
the other four — see `avg_read_len` in `kmc3_summary.csv`). With only 5 datasets, size and read
length are confounded, so it isn't possible to tell from the existing runs alone which one
actually drives the effect.

`kmc3_bench.sh readlen-dataset` / `readlen-check` build a derived dataset from an
already-downloaded one by concatenating every `--rl-join` consecutive FASTQ records into one
pseudo-read (quality strings are just concatenated too; this is a mechanical transform to control
read length going into KMC's splitter/signature path, not simulated long-read sequencing — the
sequence at each join point has no biological meaning). Both commands target a k-mer count
automatically: they probe a sample of the source's native read length and compute how many
pseudo-reads to emit, so there's no manual arithmetic. Output is registered as a new dataset id
under `$KMC3_DATA_DIR/data/`, so the normal harness works on it unmodified (pass `--mem`
explicitly — derived datasets aren't in `data/datasets.csv`, so there's no `mem_gb` to look up).

**One-shot (recommended):** `readlen-check` builds *both* size-matched variants from one source —
`join=1` (native read length, isolates **size**) and `join=<N>` (isolates **length**, the direct
"does the win reappear" replication) — then benchmarks both and exports, in one command:

```bash
# needs f_vesca already benchmarked once (for its lex n_total_kmers); or pass --rl-target-kmers directly
benchmark/kmc3_bench.sh readlen-check \
  --rl-source h_sapiens_2 --rl-against f_vesca --rl-join 4 \
  -k 28 -n 5 --mem 12 --bins 512 -t 12
```

This creates `h_sapiens_2_shortsub` (join=1) and `h_sapiens_2_longsub` (join=4), each capped to
f_vesca's `n_total_kmers`, runs `-O all -n 5` on both, and re-exports — `--rl-join 4` on
h_sapiens_2's ~100bp native reads lands at ~400bp, close to f_vesca's ~350bp. Both derived
datasets only read as much of `h_sapiens_2`'s FASTQ as needed (a few percent of it), so building
them is cheap even though the source is ~106 GB; the two `-O all -n 5` benchmark runs are the real
cost (comparable to the `f_vesca` sweep).

**Step by step**, or to build just one variant, use `readlen-dataset` directly:

```bash
benchmark/kmc3_bench.sh readlen-dataset \
  --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_shortsub --rl-join 1 --rl-against f_vesca
benchmark/kmc3_bench.sh readlen-dataset \
  --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_longsub --rl-join 4 --rl-against f_vesca

benchmark/kmc3_bench.sh run -d h_sapiens_2_shortsub -k 28 -O all -n 5 --mem 12 --bins 512 -t 12 --force
benchmark/kmc3_bench.sh run -d h_sapiens_2_longsub  -k 28 -O all -n 5 --mem 12 --bins 512 -t 12 --force
```

`--rl-target-kmers <N>` overrides the `--rl-against` lookup directly (e.g. if `f_vesca` hasn't been
benchmarked yet: f_vesca's own `n_total_kmers` is ≈4,134,079,185 at k=28).

Compare GreedyMini's `wall_time_z`/`wall_time_ratio` against lex for both datasets in the resulting
`summary` sheet:

- Win reappears on `h_sapiens_2_longsub` but not `h_sapiens_2_shortsub` → **read length** drives it.
- Win reappears on `h_sapiens_2_shortsub` (regardless of `h_sapiens_2_longsub`) → **size** drives it.
- Win reappears on neither → something specific to `f_vesca` beyond size/length (e.g. its
  particular k-mer/GC composition), not a general size or length effect.

Both commands are thin wrappers around `python -m scripts make-readlen-dataset`
(`scripts/make_readlen_dataset.py`) — reach for that directly only if scripting something the bash
commands don't cover (e.g. `--max-output-reads` for a literal read count instead of
`--target-kmers`).

## Dataset files

Catalog and reference CSVs in `benchmark/data/` (committed; not under per-dataset download dirs):

| file | columns | role |
|---|---|---|
| `datasets.csv` | `dataset`, `label`, `mem_gb`, `filename`, `url` | FASTQ list for download; `mem_gb` → kmc `-m` unless `--mem` overrides; empty `url` → ENA path from `filename` |
| `paper_kmc3.csv` | `dataset`, `k`, `platform`, `input_format`, `paper_ram_gb`, `paper_disk_gb`, `paper_time_s` | Published KMC 3 numbers; `platform` is `xeon_hdd` (default) or `xeon_ssd`; joined in `kmc3_vs_paper.csv` |

## Notes

- Current `kmc` default signature length is `-p9`; KMC 3.0.0 used the then-default (7). Pass `--signature-len 7` for closer algorithmic match.
- Custom signature order: `kmc -co<table.gm> -ol<label> [-os<seed>]` (raw binary `.gm`; seed breaks rank ties at runtime), `kmc -co<table_dna.gm> -ol<label>` (pre-built `_dna.gm`; seed is in the file header), `-or hash|pcg -ol<label> [-os<seed>]` (random xor or xor+PCG). Default without these flags is unchanged KMC lex signatures. Sidecar `{output}.kmc_co` / `.kmc_col` persist mode, seed, and table.
- Ratio color scales (`gm_lut`, `summary`, `readlen_check`, `runs`) clamp their min/max to always straddle `1.0` (`_ratio_color_scale_bounds`, driven by shadow `__{field}_value` floats captured alongside each ratio formula) instead of relying on the data's raw min/max as the color-scale stops. Fixes a case in `gm_lut` where `disk_ratio`'s data max was exactly `1.0` (GreedyMini never loses to lex in that sweep), which made the white/`1.0` midpoint coincide with the max stop and suppressed white/red entirely.
- `runs` rows are sorted in the same per-config blocks as `summary` (dataset -> bins -> threads -> mem_limit -> strict_memory -> scheme, then seed_idx within a scheme), and its "id / params" columns follow `summary`'s field order for every field the two sheets share.
- **`-O lex`**: KMC's default `-O default` lexicographic order is not *pure* canonical lex — `CMmer::is_allowed()` (`kmc_api/mmer.h`) excludes poly-A/poly-T-type motifs from ever being a minimizer (a denylist from the original KMC2 paper, to bound the largest bin's size; see §2.2 there). `-O lex` runs true canonical lex with no denylist, for isolating that specific effect — implemented as `-or hash -os0` under the hood (verified order-equivalent to `MIN(kmer, revcomp(kmer))` with no exclusions; no KMC core change needed), but reported end-to-end as its own `lexicographic` scheme (not merged into the `hash` scheme it's plumbed through, and distinct from `default`, which is what KMC actually ships). Standalone comparison scheme, always seed 0 (deterministic; `-n` just repeats it for wall-time noise stats, like `-O default`) — not part of `-O all`.
- ENA now publishes G. gallus / M. balbisiana as paired `_1`/`_2` FASTQ (the 2015 DDBJ unpaired `.fastq.bz2` names 404). `data/datasets.csv` lists the current ENA files.
