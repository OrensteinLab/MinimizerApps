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

Paths come from `benchmark/kmc/.env` (see `.env.example`). Copy the example if `.env` is missing:

```bash
cp benchmark/kmc/.env.example benchmark/kmc/.env
# edit KMC3_BENCH_ROOT — use a directory outside the repo for large downloads
```

| variable | meaning |
|---|---|
| `KMC3_BENCH_ROOT` | root of all generated data (see layout below) |
| `KMC3_DATA_DIR` | optional; default `$KMC3_BENCH_ROOT/data` (datasets live at `$KMC3_DATA_DIR/<id>/`) |
| `KMC3_RESULTS_DIR` | optional; default `$KMC3_BENCH_ROOT/results` |

CLI `--data-dir` / `--results-dir` override the env file. Existing environment variables are not overwritten by `.env`.

On-disk layout under `$KMC3_BENCH_ROOT` (no extra `benchmark/` folder):

```
$KMC3_BENCH_ROOT/
  data/<dataset>/              # fastq_gz/, fastq_plain/, files_gz.lst, files_plain.lst
  results/
    runs/                      # per-run outputs (see Outputs)
    csvs/                      # kmc3_*.csv, incl. kmc3_density_{expected_schemes,particular}.csv
    kmc3_results.xlsx          # workbook
```

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
# → $KMC3_BENCH_ROOT/results/csvs/kmc3_{runs,summary,vs_paper}.csv
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
| `../../common/minimizer-schemes/gm_serializer.cpp` | Rank list (`.txt`) → binary `.gm` (C++, shared with kraken2, built via that directory's `Makefile`) |
| `tests/` | Sanity checks — scheme/density correctness, not report formatting; `make venv` once, then `make test` (pytest) |
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
benchmark/kmc3_bench.sh run -S all -n 10 --force       # full multi-seed benchmark
benchmark/kmc3_bench.sh export                           # re-export CSV/XLSX only
benchmark/kmc3_bench.sh all -d f_vesca -k 28 -n 1       # download + run + export
```

Defaults: all five datasets, `k=28,55`, gzip input, 12 threads, `-S all`, `-n 1` repeat of each scheme
(including lex — see "Run order and the lex baseline" below).

**Experiments** (multi-config sweeps from YAML; KMC allows `-n` bins in **64..2000**):

```bash
benchmark/kmc3_bench.sh experiment list
benchmark/kmc3_bench.sh experiment sweep_config                 # all datasets, -S all, 5 seeds, 8 configs
benchmark/kmc3_bench.sh experiment sweep_gm_lut_fixed_dataset   # f_vesca, GM-DNA vs GM-bin, 6 layouts
benchmark/kmc3_bench.sh experiment sweep_gm_lut_fixed_config    # all datasets, GM-DNA vs GM-bin, t=24 n=1024
benchmark/run_config_sweep.sh                                   # same as sweep_config
```

`sweep_config` runs **800 KMC jobs** (5 datasets × 8 configs × [4 schemes × 5 repeats], lex repeated the same as the custom schemes — see "Run order and the lex baseline") at k=28 gzip only. Smoke-test with `run -d f_vesca -k 28 -n 1`. Full sweep including `h_sapiens_3` (~600 GB gzip) needs multi‑TB disk and days of wall time.

GM LUT experiments compare the native DNA table (`_dna.gm`) vs the binary rank LUT (`gm_lut: dna|binary` + `run_tag`); the binary LUT is kmc's default backend for a binary `.gm` source, use `-ode` to opt into the DNA-expanded table instead:

| experiment | scope | KMC jobs | ~wall |
|---|---|---|---|
| `sweep_gm_lut_fixed_dataset` | f_vesca, bins 512/1024/1536 × t12/t24 | 12 | ~4 min |
| `sweep_gm_lut_fixed_config` | five datasets, t=24 bins=1024 | 10 | ~1.5 h |

Experiment files live in `experiments/`. Each file has shared `defaults` plus a `configs` list.
Per-config keys: `threads`, `bins`, `mem`, `strict_mem`, `datasets`, `k`, `input`,
`scheme`, `order_label` (the kmc `-ol` label), `n_seeds`, `force`, `gm_lut`, `run_tag`, …

Schemes are the canonical names of `benchmark/common/scheme_names.py` / `orders/README.md`:
`-S default|lex|xor|xor-pcg-hash|greedymini|all` (kmc has no mod-minimizer; `default` and `lex`
are kmc's own baselines; no other spellings are accepted). The GreedyMini identity is separate
parameters, never part of the scheme name: `--gm-lut dna|binary` (or `GM_LUT`, yaml `gm_lut`; the
table LUT backend, same order, different storage) and `--projection at0-cg1|at1-cg0` (or
`GM_PROJECTION`; `at1-cg0` = kmc `-oat1`, needs `gm_lut=binary`: the flipped DNA->binary
projection, a different order). Run directory: `greedymini/k<p>_<projection>/lut_<gm_lut>/<seed_idx>`.

Example (`experiments/sweep_config.yaml` excerpt):

```yaml
# sweep_config.yaml (excerpt)
defaults:
  datasets: f_vesca,g_gallus,m_balbisiana,h_sapiens_2,h_sapiens_3
  k: 28
  input: gzip
  scheme: all
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

Custom signature order (greedymini table, xor, or xor-pcg-hash):

GreedyMini UHS rank lists live in `../../orders/GreedyMini/*.txt` (shared
with kraken2). Build tools:

| Tool | Role |
|---|---|
| `../../common/minimizer-schemes/gm_serializer` (C++, shared with kraken2) | Rank list (`.txt`) → binary `.gm` |
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
make -C ../common/minimizer-schemes
../common/minimizer-schemes/gm_serializer --canonical -k 9 \
    ../orders/GreedyMini/w20_k9_canonical.txt \
    -o ../orders/GreedyMini/w20_k9_canonical.gm

make -C tools
tools/gm_sigma_extend \
    ../orders/GreedyMini/w20_k9_canonical.gm \
    --seed 0xDEADBEEF \
    -o ../orders/GreedyMini/w20_k9_canonical_dna.gm

# All schemes + N seeds (default + xor + xor-pcg-hash + greedymini)
benchmark/kmc3_bench.sh \
  -d f_vesca,g_gallus,m_balbisiana,h_sapiens_2,h_sapiens_3 \
  --kmc bin/kmc -k 28,55 -t 12 \
  -S all -n 10
```

`-S all` runs default, xor (`xor_s<seed>`), xor-pcg-hash (`xor-pcg-hash_s<seed>`), and greedymini (`gm_w{k-p+1}_k{p}_canonical_s<seed>`), each custom scheme with `-n` seeds from `config.yaml`. Optional `--custom-order` overrides the per-k `.gm` path.

Predefined seeds are in `config.yaml` (10 fixed hex values). Custom runs use `seed_idx` `0…n-1` in the results tree; `signature_order_label` in exports still includes the hex seed (e.g. `xor_s0xDEADBEEF`).

### Run order and the lex baseline

`-S all` **interleaves** default/xor/xor-pcg-hash/greedymini within each repeat `i` (default, then xor, then xor-pcg-hash, then
greedymini, then repeat for `i+1`) instead of running all of one scheme back-to-back before starting the
next. It also repeats lex `-n` times, the same as the custom schemes, instead of running it once.

Both changes fix the same problem: with lex run once, first, and each custom scheme run as an
uninterrupted 5-seed block afterward, any time-of-day or machine-state drift over the course of a
multi-hour run lands as a step change *between* schemes — which reads exactly like a real
ordering effect but isn't one. Interleaving spreads that drift evenly across all four schemes
instead of loading it onto lex's single sample, and repeating lex gives it a real mean/std instead
of one unrepeatable data point. `-S default -n N` also accepts `N > 1` now, for the same reason,
when running lex on its own.

This changes `n_repeats` and adds `{metric}_std` for `scheme == default` rows in every
export (previously always `n_repeats=1`, `std=0`) — `kmc3_vs_paper.csv` and the `runs`/`summary`
sheets need no code changes for this, the existing per-config aggregation already treats lex like
any other scheme once it has real repeats.

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
        {scheme}/            # default | lex | xor | xor-pcg-hash | greedymini/k<p>_<projection>/lut_<gm_lut>
          {seed_idx}/        # 0…n-1 for every scheme, including lex (-n repeats)
            kmc.json         # KMC -j summary
            time.txt         # GNU /usr/bin/time -v
            meta.json        # run config + command
            out.kmc_*        # sidecars (removed unless --keep-db)
```

**Config segment** (unless `--run-tag`): `{input}_n{bins}_t{threads}_m{mem}` with `_sm` when `--strict-mem` (e.g. `gzip_n512_t12_m12`, `gzip_n1536_t24_m8_sm`). This string is also `config_tag` in CSV exports.

**run_id** in CSVs is the relative path under `runs/`, e.g. `f_vesca/k28_p9/gzip_n512_t12_m12/xor/2`.

`export` (also at the end of `run` / `all`) writes CSVs to `$KMC3_BENCH_ROOT/results/csvs/` and the workbook to `$KMC3_BENCH_ROOT/results/`:

| file | rows | purpose | in `.xlsx` |
|---|---|---|---|
| `kmc3_runs.csv` | one per run | Metrics + telemetry | `runs` |
| `kmc3_summary.csv` | one per run (per seed) | Raw headline metrics | — |
| `kmc3_vs_paper.csv` | one per lex run | vs `data/paper_kmc3.csv` (`count_wall_s_ratio` = ours ÷ paper) | `vs_paper` |
| `kmc3_density_expected_schemes.csv` | one per scheme | theoretical expected density, any scheme (see `python -m scripts expected-density` below) | `schemes` |
| `kmc3_signature_cost.csv` | one per dataset x signature variant | `experiments/profile_run.sh`: ns per base of the signature and of the window minimum, mean / std over seeds, ratio to the default order | `signature_cost` |
| `significance_kmc.csv` | one per metric x setting x scheme | paired test of the change from the default scheme (`benchmark/common/significance.py`, see `benchmark/README.md`) | `significance` |
| `density_effect_kmc.csv` | one per metric x (all, setting = `dataset/n_bins/threads/mem_limit_gb`) | Pearson r of the metric and the density factor, both as ratios to the default scheme's run of the same seed (`benchmark/common/density_effect.py`) | `density_effect` (points: `density_effect_points`) |
| `kmc3_results.xlsx` | — | `overview`, `schemes`, `datasets`, `summary`, `readlen_check`, `gm_lut`, `runs`, `vs_paper`, `signature_cost` (after `experiments/profile_run.sh`), `density_effect`, `density_effect_points`, `significance` | — |

Re-export only: `benchmark/kmc3_bench.sh export` (or `python3 -m scripts export --results-dir … --paper-platform xeon_hdd|xeon_ssd`).

`kmc3_vs_paper.csv` is omitted if `paper_kmc3.csv` is missing. Custom vs lex ratios live on the `summary` sheet (`{metric}_ratio` columns); export still **fails** if custom k-mer counts differ from lex.

### CSV columns (main fields)

**Shared keys** (all result CSVs): `dataset`, `k`, `w`, `L`, `input_format` (`gzip` \| `plain`), `config_tag` (e.g. `gzip_n512_t12_m12`, `gzip_n1536_t24_m8_sm`), `n_bins`, `threads`, `strict_memory`, `mem_limit_gb`.

**Parameter names (all sheets and CSVs, shared with the other apps):** `k` = minimizer (signature) length (kmc `-p`, old column `p`), `w` = window in k-mers, `L` = `k + w - 1` = kmc's k-mer length (kmc `-k`, old column `k`) -- note the swap. They appear together as `k, w, L`. The kmc CLI flags, run-dir names (`k28_p9`), meta.json and `data/paper_kmc3.csv` (whose `k` is the k-mer length `L`) keep kmc's own names; the exporter renames on output.

**Schemes:** `scheme` is one of `default`, `lex`, `xor`, `xor-pcg-hash`, `greedymini` (from `-S default|lex|xor|xor-pcg-hash|greedymini`). GreedyMini runs carry the separate columns `gm_k` (= minimizer length `k` of the table), `projection` (`at0-cg1` / `at1-cg0`; `-oat1` = `at1-cg0`), `gm_table` (e.g. `w20_k9_canonical`) and `gm_lut` (`dna` | `binary`), blank for the other schemes. A summary/schemes row is one full identity (scheme, gm_k, projection, gm_table, gm_lut), so the two projections and the two LUT backends are never pooled; expected density is modelled for both projections (the at1 one flips the binary table, as `kmc -oat1` does). Custom runs also have `signature_order_label` (e.g. `xor_s0xDEADBEEF`, `gm_w20_k9_canonical_s0xdeadbeef`). Runs written by `scripts/meta.py` store these fields in meta.json; kmc.json reports `signature_order_mode` (`xor` / `xor-pcg-hash` / `greedymini`) and `signature_order_gm_lut`.

| file | additional columns |
|---|---|
| `kmc3_runs.csv` | `run_id`, `label`, `seed_idx`, metrics, telemetry (`user_s`, `cpu_pct`, `host`, …) |
| `kmc3_summary.csv` | per-run `count_wall_s`, `count_max_rss_gb`, … (same as one row per seed) |

Both CSVs carry every dataset unfiltered, including read-length-check derived ones — only the xlsx `summary` sheet (below) excludes them in favor of its own `readlen_check` sheet.

### XLSX sheets

| sheet | content |
|---|---|
| `overview` | metrics reference: what every non-obvious column means, its units, and how to read it (first sheet; formerly named `index`, `legend`, `metrics`) |
| `schemes` (+ `kmc3_schemes.csv`) | one row per (scheme variant, `k`, `w`, `L`): `scheme, k, w, L, gm_k, projection, gm_table, gm_lut`, `description`, `cli`, `is_allowed_denylist` (no `mod_t`: not a kmc parameter), `expected_density_factor_mean/std` (theoretical, i.i.d. random DNA; mean/std over the benchmark's seeds for xor/xor-pcg-hash/greedymini, std 0 for deterministic default/lex), then per dataset `<dataset>_particular_density_factor_mean/std` (mean/std over seeds of the kmc runs of that variant on that dataset, any bins/threads/input config; where a variant has no run on a dataset, a density-only scan -- std 0 -- from `kmc3_density_particular.csv`). Every variant that appears in `kmc3_runs.csv` gets a row, and `export` **fails** if any density cell would be empty. |
| `datasets` | one row per dataset: size / read-length / config stats (from its own `default` run) plus published KMC3-paper reference numbers where available |
| `summary` | mean / std / ratio formulas; charts below (paper datasets only, see `readlen_check`) |
| `readlen_check` | read-length-confound-check derived datasets (`*_shortsub`/`*_longsub`, …), aggregated like `summary` but kept out of it since they aren't paper datasets — present only when such runs exist |
| `gm_lut` | GreedyMini DNA-vs-binary LUT sweep, aggregated like `summary` |
| `runs` | one row per seed (metrics + telemetry) — every dataset, unfiltered; rows grouped in the same dataset → bins → threads → mem → strict_mem → scheme order as `summary` (seed_idx orders repeats within a scheme block), and its "id / params" columns share `summary`'s field order as far as per-run vs per-config granularity allows |
| `vs_paper` | default (lex) runs vs published KMC3 numbers |
| `signature_cost` | only after `experiments/profile_run.sh`: cost per base of KMC3's signature with the default order and the GreedyMini order (DNA and binary table), alone and with the window minimum, and superwindows per kb; mean / std over seeds, ratio to the default |

### XLSX `summary` sheet (aggregated)

**Column order** (every sheet/CSV): dataset columns (`dataset`, `n_total_kmers`, `avg_read_len`, `n_reads`, `input_format`, `label`), then scheme columns (`k`, `w`, `L`, `scheme`, `gm_k`, `projection`, `gm_table`, `gm_lut`, `signature_order_label`), then the app parameters (`n_bins`, `threads`, `strict_memory`, `mem_limit_gb`), then the seed columns (`n_repeats` / `seed_idx` / `seed`), then the metrics (`vs_paper` has no scheme columns).

The category-band row and column headers are rows 1-2, the data table starts at row 3 (rows 1-2 frozen). A medium border separates each setting (dataset + input format + bins + threads + memory). Right under the table, a `correlation with density` row holds, for every `{metric}_ratio` column, a live Excel `=CORREL(density_factor_ratio range, that column's range)` over all summary rows (ratios = scheme / lex of the same setting, so rows of different settings are comparable; the value is also cached for viewers that do not recalculate).

One row per **dataset / k / input / config / scheme**, sorted by **dataset → bins → threads → memory → scheme**, with metrics aggregated across seeds. `{metric}_ratio` for wall time, RAM, and disk use a green (better) → yellow (≈1) → red (worse) color scale.

| column pattern | meaning |
|---|---|
| `n_repeats` | number of seed runs in the group |
| `{metric}_mean`, `{metric}_std` | mean and sample std across seeds (`count_wall_s`, `count_max_rss_gb`, `tmp_disk_gb`, stages, `n_super_kmers`, …) |
| `{metric}_ratio` | Excel formula `=mean / lex_mean` (lex row = **1**); vs the `default` (lex) baseline for the same **dataset / k / input / config sweep** |

### `python -m scripts expected-density` (theoretical density, any scheme)

(`scripts/expected_density.py` is a thin adapter over the shared driver `benchmark/common/expected_density_driver.py`: the C++ tool, variant x seed loop, mean/std, cache.)

density_factor on uniformly random DNA (no real sequence content, no read boundaries) for **every** scheme — `default`, `lex`, `xor`, `xor-pcg-hash`, `greedymini` (both projections) — counted by KMC's own splitter loop and signature orders (`benchmark/common/expected_density/`, the `kmc3` tool: superwindows and k-mers exactly as `kmc` counts them, 10^9 bases, a string of its own per seed). xor/xor-pcg-hash/greedymini are seed-dependent, so this reports mean **and** std across this benchmark's predefined seeds (`default`/`lex` are deterministic — no seed, std 0 on the sheet). Feeds the `schemes` sheet's `expected_density_factor_mean`/`expected_density_factor_std` columns.

```bash
python3 -m scripts expected-density                    # k=28 p=9 (this benchmark's defaults)
python3 -m scripts expected-density --force             # recompute even if a matching cache exists
```

Writes `kmc3_density_expected_schemes.csv` to `results/csvs/` (`--results-dir` is the results root). Like `kmc3_density_particular.csv`, this is a **separate, explicit step** — `export`/`schemes_rows()` only *read* the cache (warning to stderr if missing), they never trigger the (one-time, ~minutes) computation themselves, so routine `run`/`export` stays fast.

### `python -m scripts density particular-schemes` (particular density without a kmc run)

Density-only scan (no kmc counting; same `scan_fastq` as `density particular`, first `--max-kmers` 12M k-mers) of one scheme variant's DNA comparator on FASTQs under `$KMC3_DATA_DIR` (`$KMC3_BENCH_ROOT/data`); upserts rows keyed `(scheme, projection, dataset, k, L)` into `kmc3_density_particular.csv` (legacy per-dataset rows have a blank `scheme`). Used for the `schemes` cells with no kmc run, e.g. lex and the greedymini binary/at1 variants on the sub-sampled human sets:

```bash
python3 -m scripts density particular-schemes --scheme lex --datasets h_sapiens_2_longsub,h_sapiens_2_shortsub,h_sapiens_2_shortsub_mf11
python3 -m scripts density particular-schemes --scheme greedymini --projection at1-cg0 --datasets ...
```

The scan's super-k-mer count equals kmc's on the same reads (checked for default, lex and greedymini at0/at1 on a 100k-read sample).

**Charts** (below the summary table): one row per **dataset**, **ratio vs lex** per config sweep — x = scheme (xor / xor-pcg-hash / greedymini), **legend = config**. Bar height = mean paired-seed ratio vs lex for that dataset and sweep; error bars = std across seeds. One chart per metric (wall time, RAM, disk).

| file | additional columns |
|---|---|
| `kmc3_vs_paper.csv` | `paper_platform`, `n_repeats`, `count_wall_s`, `count_wall_s_std`, `paper_count_wall_s`, `count_wall_s_ratio`, `count_stage1_s`, `count_stage2_s`, `count_max_rss_gb`, `paper_count_max_rss_gb`, `tmp_disk_gb`, `paper_disk_gb`, `n_unique_counted` |

One row per **lex config** (bins/threads/memory), averaged across its `n_repeats` repeats — not one row per repeat.

**Units:** all memory and space metrics are **decimal GB** (1 GB = 10^9 bytes) via `benchmark/common/units.py`; columns end in `_gb`. `count_max_rss_gb` = peak RSS (GNU time reports KiB: KiB x 1024 / 10^9; meta/kmc raw values are unchanged, the conversion happens at export), `tmp_disk_gb` / `tmp_disk_strict_mem_gb` = KMC temp size (bytes / 10^9), `paper_count_max_rss_gb` / `paper_disk_gb` = the supplement's GB (`data/paper_kmc3.csv` keeps its own `paper_ram_gb` column, mapped on load). Old names: `ram`, `ram_gb`, `max_rss_kb`, `disk`, `disk_strict_mem`, `paper_ram`, `paper_disk` and the `{ram,disk}_{mean,std,ratio,z}` columns (now `count_max_rss_gb_*`, `tmp_disk_gb_*`). `mem_limit_gb` stays the user's `kmc -m` setting (GB).

**Headers:** the xlsx uses exactly the CSV column names (no abbreviations); every column has a row in the `overview` sheet. `runs` telemetry: `user`, `sys`, `cpu_pct`, `major_pagefaults`, `minor_pagefaults`, `fs_inputs` (512-byte blocks), `read_gb` (decimal GB), `fs_outputs`, `count_pinned_cpus` (CPU list from the run's `cpus.txt`, kept even when empty), `ctx_voluntary`, `ctx_involuntary`, `harness_wall_s` (was `python_wall_s`), `kmc_version`, `host`, `exit_code`. `disk_gb` is now `tmp_disk_gb`, `n_seeds` is `n_repeats`, `mem_limit` is `mem_limit_gb`.

**Validation:** `export` **fails** if a custom-scheme run’s unique / counted / total k-mer counts differ from the matching lex run for the same dataset, k, input, and config. Super-k-mer counts may differ. Custom runs without a lex baseline for that key also fail.

Times will not match 2017 hardware; use `count_wall_s_ratio` in `kmc3_vs_paper.csv`. Gzip vs plain **unique counted k-mers** for the same k should match (export warns on mismatch).

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
under `$KMC3_DATA_DIR/` (`$KMC3_BENCH_ROOT/data`), so the normal harness works on it unmodified (pass `--mem`
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
f_vesca's `n_total_kmers`, runs `-S all -n 5` on both, and re-exports — `--rl-join 4` on
h_sapiens_2's ~100bp native reads lands at ~400bp, close to f_vesca's ~350bp. Both derived
datasets only read as much of `h_sapiens_2`'s FASTQ as needed (a few percent of it), so building
them is cheap even though the source is ~106 GB; the two `-S all -n 5` benchmark runs are the real
cost (comparable to the `f_vesca` sweep).

**Step by step**, or to build just one variant, use `readlen-dataset` directly:

```bash
benchmark/kmc3_bench.sh readlen-dataset \
  --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_shortsub --rl-join 1 --rl-against f_vesca
benchmark/kmc3_bench.sh readlen-dataset \
  --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_longsub --rl-join 4 --rl-against f_vesca

benchmark/kmc3_bench.sh run -d h_sapiens_2_shortsub -k 28 -S all -n 5 --mem 12 --bins 512 -t 12 --force
benchmark/kmc3_bench.sh run -d h_sapiens_2_longsub  -k 28 -S all -n 5 --mem 12 --bins 512 -t 12 --force
```

`--rl-target-kmers <N>` overrides the `--rl-against` lookup directly (e.g. if `f_vesca` hasn't been
benchmarked yet: f_vesca's own `n_total_kmers` is ≈4,134,079,185 at k=28).

Compare GreedyMini's `count_wall_s_ratio` against lex for both datasets in the resulting
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

## Measurement protocol

Every kmc call is timed under the protocol in [`../README.md`](../README.md): it waits for its node to be quiet,
runs pinned (`$THREADS` cores + 2 spare, one NUMA node, one thread per core; the list is in
`cpus.txt` of the run directory and exported as `count_pinned_cpus`; with its memory bound to the slot's node and a step record (`step_count.json`, exported as `count_node`, `_huge_fraction`, `_foreign_busy_cores`, ... -- see the protocol)), and records wall and CPU seconds
(`count_wall_s`, `count_cpu_s`), the KMC stages, peak RSS, major page faults and file-system reads
(`major_pagefaults`, `fs_inputs` in 512-byte blocks, `read_gb`). Stage 1 reads the input, so `read_gb` > 0
means a cold input and a slower stage 1; stage 2 and `count_cpu_s` are the cache-independent numbers.
`default` is repeated like the other schemes and interleaved with them (`experiments/full_run.sh` runs
`-S all -n 5`). `experiments/full_run.sh` is the main-configuration run (150 runs: 5 datasets x 6 schemes x 5 repeats; bins 1024, threads 24).

## Notes

- Current `kmc` default signature length is `-p9`; KMC 3.0.0 used the then-default (7). Pass `--signature-len 7` for closer algorithmic match.
- Custom signature order: `kmc -co<table.gm> -ol<label> [-os<seed>]` (raw binary `.gm`; seed breaks rank ties at runtime), `kmc -co<table_dna.gm> -ol<label>` (pre-built `_dna.gm`; seed is in the file header), `-or xor|xor-pcg-hash -ol<label> [-os<seed>]` (random xor, or xor then PCG hash); `kmc.json` reports `signature_order_mode` = `xor` / `xor-pcg-hash` / `greedymini` and, for tables, `signature_order_gm_lut` = `dna` (4^p DNA table) | `binary` (2^p rank LUT; `-oat1` flips the projection). Default without these flags is unchanged KMC lex signatures. Sidecar `{output}.kmc_co` / `.kmc_col` persist mode, seed, and table.
- Ratio color scales (`gm_lut`, `summary`, `readlen_check`, `runs`) clamp their min/max to always straddle `1.0` (`_ratio_color_scale_bounds`, driven by shadow `__{field}_value` floats captured alongside each ratio formula) instead of relying on the data's raw min/max as the color-scale stops. Fixes a case in `gm_lut` where `tmp_disk_gb_ratio`'s data max was exactly `1.0` (GreedyMini never loses to lex in that sweep), which made the white/`1.0` midpoint coincide with the max stop and suppressed white/red entirely.
- `runs` rows are sorted in the same per-config blocks as `summary` (dataset -> bins -> threads -> mem_limit_gb -> strict_memory -> scheme, then seed_idx within a scheme), and its "id / params" columns follow `summary`'s field order for every field the two sheets share.
- **`-S lex`**: KMC's default `-S default` lexicographic order is not *pure* canonical lex — `CMmer::is_allowed()` (`kmc_api/mmer.h`) excludes poly-A/poly-T-type motifs from ever being a minimizer (a denylist from the original KMC2 paper, to bound the largest bin's size; see §2.2 there). `-S lex` runs true canonical lex with no denylist, for isolating that specific effect — implemented as `-or xor -os0` under the hood (verified order-equivalent to `MIN(kmer, revcomp(kmer))` with no exclusions; no KMC core change needed), but reported end-to-end as its own `lex` scheme (not merged into the `xor` scheme it's plumbed through, and distinct from `default`, which is what KMC actually ships). Standalone comparison scheme, always seed 0 (deterministic; `-n` just repeats it for wall-time noise stats, like `-S default`) — not part of `-S all`.
- ENA now publishes G. gallus / M. balbisiana as paired `_1`/`_2` FASTQ (the 2015 DDBJ unpaired `.fastq.bz2` names 404). `data/datasets.csv` lists the current ENA files.
