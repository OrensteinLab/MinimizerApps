# KMC3 counting benchmark

Runs KMC3 on the k-mer counting datasets of Kokot, Długosz & Deorowicz, *KMC 3: counting and manipulating k-mer statistics*, Bioinformatics 33(17):2759–2761 (2017), with each scheme as the signature order, and compares time, RAM and temp disk with the default scheme and with the published numbers. Competitors are not run. Dataset URLs, flags, and published KMC3 times/RAM/disk are from the paper supplement ([arXiv:1701.08022](https://arxiv.org/abs/1701.08022) `supp.pdf`).

**The paper's runs** are `experiments/full_run.sh`: L = 28 (k = 9, w = 20), gzip input, 1024 bins, 24 threads, all five datasets, 5 repeats of each scheme. It reports `default`, `xor` and `xor-pcg-hash` (one `-S all -n 5` call) and GreedyMini with the binary table and polarity of Φ `at0-cg1` (run tag `gm_proj_n1024_t24`). The same script also runs GreedyMini with the DNA table (`at0-cg1`) and with the binary table and polarity `at1-cg0`: 4 schemes in 6 variants, 150 runs. The defaults of `kmc3_bench.sh` match these settings.

## Notation

The paper's notation is used throughout, and the CSVs and xlsx use it as columns `k, w, L`:

| paper | meaning | KMC3 name |
|---|---|---|
| `k` | selected k-mer length | signature length `p` (`-p`, default 9) |
| `w` | number of k-mers per window, `w = L - k + 1` | — |
| `L` | KMC3's k-mer (L-mer) length, `L = k + w - 1` | `-k` |

KMC3's flags, the run directories (`k28_p9`), `meta.json` and `data/paper_kmc3.csv` (whose `k` is L) keep KMC3's own names; the exporter renames on output. GreedyMini has a table only for L = 28, k = 9 (`w20_k9_canonical`), so other `-k` values work only with `-S default|lex|xor|xor-pcg-hash`.

## What is measured

For each dataset and L:

- wall time (KMC JSON `Total`, plus GNU `time`)
- peak RAM (GNU `/usr/bin/time -v` max RSS)
- temp disk (KMC JSON `Tmp_size`)
- unique / counted k-mer stats

Command in the KMC3 paper (KMC 3.0.0):

```text
./kmc3 -ci2 -v -t12 -k<k> <fastq|@list> res ./bins
```

H. sapiens 3 used `-m32`, the other sets the default 12 GB (`-m12`). Gzip input used `@fv.txt` / `@gg.txt` / … file lists.

## Datasets and setup

| id | organism | gzip FASTQ (GB) |
|---|---|---|
| `f_vesca` | F. vesca | 3.4 |
| `g_gallus` | G. gallus | 24 |
| `m_balbisiana` | M. balbisiana | 51 |
| `h_sapiens_2` | H. sapiens 2 (ERA015743, ERR024163–186) | 106 |
| `h_sapiens_3` | H. sapiens 3 (ERR174324–341) | 614 |

Uncompressed copies are several times larger; the paper's runs use `--input gzip`.

Paths come from `benchmark/kmc/.env` (see `.env.example`):

| variable | meaning |
|---|---|
| `KMC3_BENCH_ROOT` | root of all generated data (see layout below); use a directory outside the repo |
| `KMC3_DATA_DIR` | optional; default `$KMC3_BENCH_ROOT/data` (datasets live at `$KMC3_DATA_DIR/<id>/`) |
| `KMC3_RESULTS_DIR` | optional; default `$KMC3_BENCH_ROOT/results` |

CLI `--data-dir` / `--results-dir` override the env file. Existing environment variables are not overwritten by `.env`.

On-disk layout under `$KMC3_BENCH_ROOT`:

```
$KMC3_BENCH_ROOT/
  data/<dataset>/              # fastq_gz/, fastq_plain/, files_gz.lst, files_plain.lst
  results/
    runs/                      # per-run outputs (see Outputs)
    csvs/                      # kmc3_*.csv, incl. the density caches
    kmc3_results.xlsx          # workbook
```

`.gitignore` excludes `.env`, downloaded FASTQ, generated `.gm` tables and results.

## Usage

Requires Python >= 3.9 (standard library only) for export and experiment YAML; Bash handles download (`curl`), prepare (`gzip`) and the kmc runs. The scripts use `kmc/.venv/bin/python` when it exists, otherwise `python3` (override with `PYTHON=`).

From the repo root (after building the shared scheme code and the `.gm` tables as in the root [README](../../README.md)):

```bash
make -C kmc all                                   # kmc/bin/kmc, the default binary of `kmc3_bench.sh`
make -C benchmark/kmc venv                        # kmc/.venv with pytest (for make test)
cp benchmark/kmc/.env.example benchmark/kmc/.env  # set KMC3_BENCH_ROOT
benchmark/kmc/kmc3_bench.sh download -d f_vesca
benchmark/kmc/kmc3_bench.sh run -d f_vesca -n 1   # quick smoke: all schemes, 1 repeat
# → $KMC3_BENCH_ROOT/results/csvs/kmc3_*.csv and results/kmc3_results.xlsx
```

`export` (also run at the end of `run`) needs the expected-density cache; build it once before the first export (see [expected density](#python--m-scripts-expected-density-theoretical-density-any-scheme)):

```bash
cd benchmark/kmc && python3 -m scripts expected-density
```

The paper's runs (needs all five datasets downloaded):

```bash
benchmark/kmc/experiments/full_run.sh             # RESULTS_NAME=<dir under KMC3_BENCH_ROOT>, default results
```

### Layout

```
benchmark/kmc/
  config.yaml          # seeds, paths, defaults (loaded by bash + Python)
  kmc3_bench.sh        # main entry (download / prepare / run / export / experiment / readlen-*)
  run_config_sweep.sh  # wrapper → experiment sweep_config
  Makefile             # venv, test
  lib/*.sh             # bash helpers
  scripts/             # python -m scripts …
  tools/               # C++ helpers + Makefile
  experiments/         # full_run.sh, profile_run.sh, *.yaml sweep definitions
  data/                # datasets.csv, paper_kmc3.csv
  tests/
```

| Path | Role |
|---|---|
| `config.yaml` | Shared seeds, paths, defaults |
| `kmc3_bench.sh` | Main entry |
| `lib/*.sh` | Bash helpers (env, python, gm, datasets, runs, experiments, read length) |
| `scripts/` | Python package (`python -m scripts …`) |
| `scripts/config.py` | Load `config.yaml` |
| `scripts/export.py` | CSV / XLSX export |
| `scripts/experiments.py` | Load / expand experiment YAML |
| `scripts/gm.py` | GreedyMini path helpers + table build/ensure (shells out to the C++ tools) |
| `scripts/meta.py` | Per-run `meta.json` writer |
| `tools/gm_sigma_extend.cpp` | Binary `.gm` → DNA `_dna.gm` |
| `tools/bench_signature_impls.cpp` | Signature-cost profile (`experiments/profile_run.sh`) |
| `experiments/full_run.sh` | The paper's runs, then export |
| `experiments/profile_run.sh` | Signature-cost profile (`signature_cost` sheet) |
| `experiments/*.yaml` | Multi-config sweep definitions |
| `data/` | Dataset catalog (`datasets.csv`, `paper_kmc3.csv`) |
| `../../orders/GreedyMini/` | Rank lists (`.txt`, shared with the other methods) and tables (`.gm`, `_dna.gm`) |
| `../../common/minimizer-schemes/gm_serializer.cpp` | Rank list (`.txt`) → binary `.gm` (shared with the other methods) |
| `tests/` | Scheme/density correctness checks: `make venv` once, then `make test` (pytest), from `benchmark/kmc` |

Python CLI (from `benchmark/kmc`; `kmc3_bench.sh` runs it there):

```bash
python3 -m scripts export --results-dir …
python3 -m scripts experiment list
python3 -m scripts gm ensure --k 28 -p 9 --seed 0xDEADBEEF
```

### kmc3_bench.sh

From `benchmark/kmc`:

```bash
./kmc3_bench.sh -h
./kmc3_bench.sh download -d f_vesca
./kmc3_bench.sh run -d f_vesca -n 1          # quick smoke (all schemes, 1 repeat)
./kmc3_bench.sh run -S all -n 5 --force       # all datasets, 5 repeats
./kmc3_bench.sh export                        # re-export CSV/XLSX only
./kmc3_bench.sh all -d f_vesca -n 1           # download + run + export
```

Defaults: all five datasets, `-k 28`, gzip input, 24 threads, 1024 bins, `-S all`, `-n 1`, GreedyMini with the binary table (`--gm-lut binary`) and polarity of Φ `at0-cg1`. `-S all` runs `default`, `xor`, `xor-pcg-hash` and `greedymini`; `-S lex` is a separate scheme (see Notes).

`-S all` interleaves the schemes within each repeat `i` (default, xor, xor-pcg-hash, greedymini, then `i+1`), and repeats `default` `-n` times like the seeded schemes, so time-of-day or machine drift spreads evenly over the schemes and every scheme has a mean and std. `-S default -n N` and `-S lex -n N` repeat that scheme N times.

**Experiments** (multi-config sweeps from YAML; KMC allows `-n` bins in 64..2000), from `benchmark/kmc`:

```bash
./kmc3_bench.sh experiment list
./kmc3_bench.sh experiment sweep_config                 # all datasets, -S all, 5 seeds, 7 configs
./kmc3_bench.sh experiment sweep_gm_lut_fixed_dataset   # f_vesca, GreedyMini DNA vs binary table, 6 layouts
./kmc3_bench.sh experiment sweep_gm_lut_fixed_config    # all datasets, GreedyMini DNA vs binary table, t24 n1024
./run_config_sweep.sh                                   # same as sweep_config
```

`sweep_config` runs 700 KMC jobs (5 datasets × 7 configs × [4 schemes × 5 repeats]) at L = 28, gzip only, with the GreedyMini DNA table. Smoke-test with `run -d f_vesca -n 1`. A full sweep including `h_sapiens_3` (614 GB gzip) needs multi-TB disk and days of wall time.

The GreedyMini table experiments compare the DNA table (`_dna.gm`) with the binary rank table (`gm_lut: dna|binary` + `run_tag`); for a binary `.gm` KMC uses the binary table unless `-ode` asks for the DNA-expanded one:

| experiment | scope | KMC jobs | ~wall |
|---|---|---|---|
| `sweep_gm_lut_fixed_dataset` | f_vesca, bins 512/1024/1536 × t12/t24 | 12 | ~4 min |
| `sweep_gm_lut_fixed_config` | five datasets, t=24 bins=1024 | 10 | ~1.5 h |

Each experiment file has shared `defaults` plus a `configs` list. Per-config keys: `threads`, `bins`, `mem`, `strict_mem`, `datasets`, `k`, `input`, `scheme`, `order_label` (the kmc `-ol` label), `n_seeds`, `force`, `gm_lut`, `run_tag`, …

Schemes are the canonical names of `benchmark/common/scheme_names.py` / `orders/README.md`: `-S default|lex|xor|xor-pcg-hash|greedymini|all` (KMC3 has no mod-minimizer; no other spellings are accepted). The GreedyMini variant is given by separate parameters, never by the scheme name: `--gm-lut binary|dna` (or `GM_LUT`, yaml `gm_lut`; the table backend, same order, different storage) and `--projection at0-cg1|at1-cg0` (or `GM_PROJECTION`; the polarity of Φ; `at1-cg0` = kmc `-oat1`, needs `gm_lut=binary`, and is a different order). Run directory: `greedymini/k<p>_<projection>/lut_<gm_lut>/<seed_idx>`.

`experiments/sweep_config.yaml` excerpt:

```yaml
defaults:
  datasets: f_vesca,g_gallus,m_balbisiana,h_sapiens_2,h_sapiens_3
  k: 28
  input: gzip
  scheme: all
  gm_lut: dna
  n_seeds: 5
  force: true

configs:
  - threads: 12
    bins: 512
  - threads: 24
    bins: 512
  # … 1024/1536 × threads 12/24, plus threads 12, bins 1536, mem 8, strict_mem
```

Run directories are nested so dataset / parameters / config / scheme / seed do not overwrite each other. `run_tag` in YAML or `--run-tag` replaces the auto-generated config segment.

### GreedyMini tables

GreedyMini rank lists live in `../../orders/GreedyMini/*.txt` (shared with the other methods). Build tools:

| Tool | Role |
|---|---|
| `../../common/minimizer-schemes/gm_serializer` | Rank list (`.txt`) → binary `.gm` |
| `tools/gm_sigma_extend` | Binary `.gm` → DNA `_dna.gm` |
| `scripts/gm.py` | Path helpers + `ensure` (shells out to the tools above) |

GreedyMini runs build missing or outdated tables (mtime check against the `.txt`):

1. `{stem}.txt` → `{stem}.gm` (binary ranks, seed-independent; the default table)
2. `{stem}.gm` → `{stem}_dna.gm` or `{stem}_dna_s0x{seed}.gm` (DNA tables, one per seed)

A `_dna.gm` embeds the tie-break seed in its header (magic `KMCANDM1`), and KMC reads it from there, so `kmc3_bench.sh` does not pass `-os` with a `_dna.gm` (`meta.json` records the embedded seed). With the binary `.gm` it passes the seed with `-os`.

Manual build (optional), from `benchmark/kmc`:

```bash
python3 -m scripts gm ensure --k 28 -p 9 --seed 0xDEADBEEF

# or step by step, with the C++ tools directly:
make -C ../../common/minimizer-schemes
../../common/minimizer-schemes/gm_serializer --canonical -k 9 \
    ../../orders/GreedyMini/w20_k9_canonical.txt \
    -o ../../orders/GreedyMini/w20_k9_canonical.gm

make -C tools
tools/gm_sigma_extend \
    ../../orders/GreedyMini/w20_k9_canonical.gm \
    --seed 0xDEADBEEF \
    -o ../../orders/GreedyMini/w20_k9_canonical_dna.gm
```

`-S all` runs `default`, xor (`xor_s<seed>`), xor-pcg-hash (`xor-pcg-hash_s<seed>`) and greedymini (`gm_w{L-k+1}_k{k}_canonical_s<seed>`), each seeded scheme with `-n` seeds from `config.yaml`. `--custom-order` overrides the `.gm` path.

`config.yaml` has 10 fixed hex seeds. Runs use `seed_idx` `0…n-1` in the results tree; `signature_order_label` in exports includes the hex seed (e.g. `xor_s0xDEADBEEF`).

Gzip then uncompressed input, from `benchmark/kmc`:

```bash
./kmc3_bench.sh all --input gzip,plain
```

`prepare` (called from `all` when `plain` is requested) gunzips into `$KMC3_BENCH_ROOT/data/<id>/fastq_plain/`. Delete those files after the plain runs if you need space.

## Outputs

Per-run artifacts live under `$KMC3_BENCH_ROOT/results/runs/`:

```
runs/
  {dataset}/                 # e.g. f_vesca
    k{L}_p{k}/               # KMC3's names: -k and -p, e.g. k28_p9
      {config}/              # e.g. gzip_n1024_t24_m12, gzip_n1536_t12_m8_sm, gm_proj_n1024_t24
        {scheme}/            # default | lex | xor | xor-pcg-hash | greedymini/k<p>_<projection>/lut_<gm_lut>
          {seed_idx}/        # 0…n-1 for every scheme (-n repeats)
            kmc.json         # KMC -j summary
            time.txt         # GNU /usr/bin/time -v
            meta.json        # run config + command
            out.kmc_*        # databases (removed unless --keep-db)
```

**Config segment** (unless `--run-tag`): `{input}_n{bins}_t{threads}_m{mem}`, with `_sm` when `--strict-mem`. This string is also `config_tag` in the CSVs.

**run_id** in the CSVs is the relative path under `runs/`, e.g. `f_vesca/k28_p9/gzip_n1024_t24_m12/xor/2`.

`export` (also at the end of `run` / `all`) writes CSVs to `$KMC3_BENCH_ROOT/results/csvs/` and the workbook to `$KMC3_BENCH_ROOT/results/`:

| file | rows | purpose | in `.xlsx` |
|---|---|---|---|
| `kmc3_runs.csv` | one per run | metrics + telemetry | `runs` |
| `kmc3_summary.csv` | one per dataset / L / input / config / scheme variant | metrics aggregated over repeats (mean / std / ratio) | `summary` |
| `kmc3_schemes.csv` | one per scheme variant | expected and particular density | `schemes` |
| `kmc3_gm_lut.csv` | one per GreedyMini table-sweep config | the `gm_lut_*` run tags, aggregated like `summary` (only when such runs exist) | `gm_lut` |
| `kmc3_vs_paper.csv` | one per `default` config | vs `data/paper_kmc3.csv` (`count_wall_s_ratio` = ours ÷ paper) | `vs_paper` |
| `kmc3_signature_cost.csv` | one per dataset × signature variant | `experiments/profile_run.sh`: ns per base of the signature and of the window minimum, mean / std over seeds, ratio to the default order | `signature_cost` |
| `significance_kmc.csv` | one per metric × setting × scheme | paired test of the change from the default scheme (`benchmark/common/significance.py`, see `benchmark/README.md`) | `significance` |
| `density_effect_kmc.csv` | one per metric × (all, setting = `dataset/n_bins/threads/mem_limit_gb`) | Pearson r of the metric and the density factor, both as ratios to the default scheme's run of the same seed (`benchmark/common/density_effect.py`) | `density_effect` (points: `density_effect_points`) |
| `kmc3_density_expected_schemes.csv` | one per scheme variant | cache of `python -m scripts expected-density` (read by export, never written by it) | — |
| `kmc3_density_particular.csv` | one per scheme variant × dataset | cache of `python -m scripts density particular-schemes` (density-only scans) | — |
| `kmc3_results.xlsx` | — | `overview`, `schemes`, `datasets`, `summary`, `readlen_check`, `gm_lut`, `runs`, `vs_paper`, `signature_cost` (after `experiments/profile_run.sh`), `density_effect`, `density_effect_points`, `significance` | — |

Re-export only, from `benchmark/kmc`: `./kmc3_bench.sh export` (or `python3 -m scripts export --results-dir … --paper-platform xeon_hdd|xeon_ssd`).

`kmc3_vs_paper.csv` is omitted if `paper_kmc3.csv` is missing.

### CSV columns (main fields)

**Shared keys** (all result CSVs): `dataset`, `k`, `w`, `L`, `input_format` (`gzip` \| `plain`), `config_tag` (e.g. `gzip_n1024_t24_m12`), `n_bins`, `threads`, `strict_memory`, `mem_limit_gb`.

**Schemes:** `scheme` is one of `default`, `lex`, `xor`, `xor-pcg-hash`, `greedymini`. GreedyMini runs carry the separate columns `gm_k` (= selected k-mer length `k` of the table), `projection` (polarity of Φ: `at0-cg1` / `at1-cg0`; `-oat1` = `at1-cg0`), `gm_table` (e.g. `w20_k9_canonical`) and `gm_lut` (`dna` | `binary`), blank for the other schemes. A summary/schemes row is one full variant (scheme, gm_k, projection, gm_table, gm_lut), so the two polarities and the two table backends are never pooled. Seeded runs also have `signature_order_label` (e.g. `xor_s0xDEADBEEF`, `gm_w20_k9_canonical_s0xdeadbeef`). `meta.json` stores these fields; `kmc.json` reports `signature_order_mode` (`xor` / `xor-pcg-hash` / `greedymini`) and `signature_order_gm_lut`.

| file | additional columns |
|---|---|
| `kmc3_runs.csv` | `run_id`, `label`, `seed_idx`, metrics, telemetry (`user_s`, `cpu_pct`, `host`, …) |
| `kmc3_summary.csv` | `n_repeats`, `{metric}_mean`, `{metric}_std`, `{metric}_ratio` |
| `kmc3_vs_paper.csv` | `paper_platform`, `n_repeats`, `count_wall_s`, `count_wall_s_std`, `paper_count_wall_s`, `count_wall_s_ratio`, `count_stage1_s`, `count_stage2_s`, `count_max_rss_gb`, `paper_count_max_rss_gb`, `tmp_disk_gb`, `paper_disk_gb`, `n_unique_counted` |

`kmc3_runs.csv` and `kmc3_summary.csv` carry every dataset, including the read-length-check datasets; only the xlsx `summary` sheet excludes those, in favor of its `readlen_check` sheet. `kmc3_vs_paper.csv` has one row per `default` config (bins / threads / memory), averaged over its repeats.

**Units:** all memory and space metrics are decimal GB (1 GB = 10^9 bytes, `benchmark/common/units.py`); columns end in `_gb`. `count_max_rss_gb` = peak RSS (GNU time reports KiB: KiB × 1024 / 10^9; the conversion happens at export), `tmp_disk_gb` / `tmp_disk_strict_mem_gb` = KMC temp size (bytes / 10^9), `paper_count_max_rss_gb` / `paper_disk_gb` = the supplement's GB (`data/paper_kmc3.csv` has its own `paper_ram_gb` column, mapped on load). `mem_limit_gb` is the `kmc -m` setting (GB).

**Headers:** the xlsx uses the CSV column names; every column has a row in the `overview` sheet. `runs` telemetry: `user`, `sys`, `cpu_pct`, `major_pagefaults`, `minor_pagefaults`, `fs_inputs` (512-byte blocks), `read_gb` (decimal GB), `fs_outputs`, `count_pinned_cpus` (CPU list from the run's `cpus.txt`, kept even when empty), `ctx_voluntary`, `ctx_involuntary`, `harness_wall_s`, `kmc_version`, `host`, `exit_code`.

**Validation:** `export` fails if a seeded scheme's unique / counted / total k-mer counts differ from the matching `default` run (same dataset, L, input and config), or if that run is missing. Superwindow counts may differ. Gzip vs plain unique counted k-mers for the same L should match (export warns on mismatch).

Times will not match 2017 hardware; use `count_wall_s_ratio` in `kmc3_vs_paper.csv`.

### XLSX sheets

| sheet | content |
|---|---|
| `overview` | metrics reference: what every non-obvious column means, its units, and how to read it (first sheet) |
| `schemes` (+ `kmc3_schemes.csv`) | one row per (scheme variant, `k`, `w`, `L`): `scheme, k, w, L, gm_k, projection, gm_table, gm_lut`, `description`, `cli`, `is_allowed_denylist`, `expected_density_factor_mean/std` (i.i.d. random DNA; mean/std over the seeds for xor/xor-pcg-hash/greedymini, std 0 for the deterministic `default`/`lex`), then per dataset `<dataset>_particular_density_factor_mean/std` (mean/std over the seeds of the kmc runs of that variant on that dataset, any config; where a variant has no run on a dataset, a density-only scan, std 0, from `kmc3_density_particular.csv`). Every variant in `kmc3_runs.csv` gets a row, and `export` fails if any density cell would be empty. |
| `datasets` | one row per dataset: size / read-length / config stats (from its `default` run) plus the KMC3 paper's numbers where available |
| `summary` | mean / std / ratio formulas; charts below (paper datasets only, see `readlen_check`) |
| `readlen_check` | read-length-check datasets (`*_shortsub` / `*_longsub`, …), aggregated like `summary` — only when such runs exist |
| `gm_lut` | GreedyMini DNA vs binary table sweep, aggregated like `summary` |
| `runs` | one row per run (metrics + telemetry), every dataset; rows grouped like `summary` (dataset → bins → threads → mem → strict_mem → scheme, then `seed_idx`), id / params columns in `summary`'s field order |
| `vs_paper` | `default` runs vs the published KMC3 numbers |
| `signature_cost` | only after `experiments/profile_run.sh`: cost per base of KMC3's signature with the default order and the GreedyMini order (DNA and binary table), alone and with the window minimum, and superwindows per kb; mean / std over seeds, ratio to the default |

### XLSX `summary` sheet

**Column order** (every sheet/CSV): dataset columns (`dataset`, `n_total_kmers`, `avg_read_len`, `n_reads`, `input_format`, `label`), then scheme columns (`k`, `w`, `L`, `scheme`, `gm_k`, `projection`, `gm_table`, `gm_lut`, `signature_order_label`), then the method's parameters (`n_bins`, `threads`, `strict_memory`, `mem_limit_gb`), then the seed columns (`n_repeats` / `seed_idx` / `seed`), then the metrics (`vs_paper` has no scheme columns).

Header rows 1–2 (category band and column names, frozen); data from row 3. A medium border separates each setting (dataset + input format + bins + threads + memory). Under the table, a `correlation with density` row holds, for every `{metric}_ratio` column, a live `=CORREL(density_factor_ratio range, that column's range)` over all summary rows (ratios to the `default` scheme of the same setting, so settings are comparable; the value is also cached for viewers that do not recalculate).

One row per dataset / L / input / config / scheme variant, sorted by dataset → bins → threads → memory → scheme, with metrics aggregated over the repeats. `{metric}_ratio` for wall time, RAM and disk use a green (better) → white (1.0) → red (worse) color scale, with 1.0 always inside the scale.

| column pattern | meaning |
|---|---|
| `n_repeats` | number of runs in the group |
| `{metric}_mean`, `{metric}_std` | mean and sample std over the repeats (`count_wall_s`, `count_max_rss_gb`, `tmp_disk_gb`, stages, `n_super_kmers`, …) |
| `{metric}_ratio` | Excel formula `=mean / default_mean` (`default` row = 1) for the same dataset / L / input / config |

**Charts** (below the summary table): one row per dataset, ratio to `default` per config: x = scheme (xor / xor-pcg-hash / greedymini), legend = config. Bar height = mean paired-seed ratio; error bars = std over seeds. One chart per metric (wall time, RAM, disk).

### `python -m scripts expected-density` (theoretical density, any scheme)

Density factor on uniformly random DNA (no read boundaries) for every scheme — `default`, `lex`, `xor`, `xor-pcg-hash`, `greedymini` (both polarities) — counted by KMC3's own splitter loop and signature orders (`benchmark/common/expected_density/`, the `kmc3` tool: superwindows and k-mers as `kmc` counts them, 10^9 bases, a string of its own per seed). The seeded schemes report mean and std over the predefined seeds (`default`/`lex` are deterministic, std 0). Feeds the `schemes` sheet's `expected_density_factor_mean/std` columns. `scripts/expected_density.py` is a thin adapter over the shared driver `benchmark/common/expected_density_driver.py`.

From `benchmark/kmc`:

```bash
python3 -m scripts expected-density            # L=28 k=9 (this benchmark's defaults)
python3 -m scripts expected-density --force    # recompute even if a matching cache exists
```

Writes `kmc3_density_expected_schemes.csv` to `results/csvs/` (`--results-dir` is the results root). This is a separate step, needed once before the first export: `export` only reads the cache and fails on the empty density cells if it is missing. `experiments/full_run.sh` copies the caches from `OLD_RESULTS` (default `results`) when they exist there.

### `python -m scripts density particular-schemes` (particular density without a kmc run)

Density-only scan (no kmc counting; first `--max-kmers` 12M k-mers) of one scheme variant on the FASTQs under `$KMC3_DATA_DIR`; upserts rows keyed `(scheme, projection, dataset, k, L)` into `kmc3_density_particular.csv`. Fills the `schemes` cells of variants with no kmc run on a dataset, e.g. `lex` and the GreedyMini binary/`at1-cg0` variants on the read-length-check datasets. From `benchmark/kmc`:

```bash
python3 -m scripts density particular-schemes --scheme lex --datasets h_sapiens_2_longsub,h_sapiens_2_shortsub,h_sapiens_2_shortsub_mf11
python3 -m scripts density particular-schemes --scheme greedymini --projection at1-cg0 --datasets ...
```

The scan's superwindow count equals kmc's on the same reads (checked for `default`, `lex` and greedymini `at0-cg1`/`at1-cg0` on a 100k-read sample).

## Read-length confound check

`f_vesca` is both the smallest dataset (~4.1B k-mers) and the one with by far the longest reads (avg ~350 bp vs ~64–101 bp for the other four; `avg_read_len` in `kmc3_summary.csv`), so across the five datasets size and read length are confounded. `kmc3_bench.sh readlen-dataset` / `readlen-check` separate them.

They build a derived dataset from an already-downloaded one by concatenating every `--rl-join` consecutive FASTQ records into one pseudo-read (quality strings are concatenated too; a mechanical transform to control read length, not simulated long reads — the sequence at each join point has no biological meaning). Both target a k-mer count: they sample the source's read length and compute how many pseudo-reads to emit. The output is registered as a new dataset id under `$KMC3_DATA_DIR/`, so the usual commands work on it (pass `--mem`: derived datasets are not in `data/datasets.csv`).

**One-shot:** `readlen-check` builds both size-matched variants from one source — `join=1` (native read length, isolates size) and `join=<N>` (isolates length) — then benchmarks both and exports. From `benchmark/kmc`:

```bash
# needs f_vesca benchmarked once (for its default n_total_kmers), or pass --rl-target-kmers
./kmc3_bench.sh readlen-check \
  --rl-source h_sapiens_2 --rl-against f_vesca --rl-join 4 \
  -k 28 -n 5 --mem 12 --bins 512 -t 12
```

This creates `h_sapiens_2_shortsub` (join=1) and `h_sapiens_2_longsub` (join=4), each capped to f_vesca's `n_total_kmers`, runs `-S all -n 5` on both, and re-exports. `--rl-join 4` on h_sapiens_2's ~100 bp reads gives ~400 bp, close to f_vesca's ~350 bp. Both datasets read only a few percent of `h_sapiens_2`'s FASTQ, so building them is cheap; the two benchmark runs are the real cost (comparable to `f_vesca`).

**Step by step**, or to build one variant, use `readlen-dataset`:

```bash
./kmc3_bench.sh readlen-dataset \
  --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_shortsub --rl-join 1 --rl-against f_vesca
./kmc3_bench.sh readlen-dataset \
  --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_longsub --rl-join 4 --rl-against f_vesca

./kmc3_bench.sh run -d h_sapiens_2_shortsub -k 28 -S all -n 5 --mem 12 --bins 512 -t 12 --force
./kmc3_bench.sh run -d h_sapiens_2_longsub  -k 28 -S all -n 5 --mem 12 --bins 512 -t 12 --force
```

`--rl-target-kmers <N>` replaces the `--rl-against` lookup (f_vesca's `n_total_kmers` is ≈4,134,079,185 at L=28).

Compare each scheme's `count_wall_s_ratio` on the two datasets in the `readlen_check` sheet: an effect on `h_sapiens_2_longsub` only points to read length, one on `h_sapiens_2_shortsub` to size, and one on neither to something specific to `f_vesca`'s sequence.

Both commands wrap `python -m scripts make-readlen-dataset` (`scripts/make_readlen_dataset.py`), which also takes e.g. `--max-output-reads` for a literal read count instead of `--target-kmers`.

## Dataset files

Catalog and reference CSVs in `benchmark/kmc/data/`:

| file | columns | role |
|---|---|---|
| `datasets.csv` | `dataset`, `label`, `mem_gb`, `filename`, `url` | FASTQ list for download; `mem_gb` → kmc `-m` unless `--mem` overrides; empty `url` → ENA path from `filename` |
| `paper_kmc3.csv` | `dataset`, `k`, `platform`, `input_format`, `paper_ram_gb`, `paper_disk_gb`, `paper_time_s` | Published KMC3 numbers (`k` is L); `platform` is `xeon_hdd` (default) or `xeon_ssd`; joined in `kmc3_vs_paper.csv` |

## Measurement protocol

Every kmc call is timed under the protocol in [`../README.md`](../README.md): it waits for its node to be quiet, runs pinned (`$THREADS` cores + 2 spare, one NUMA node, one thread per core; the list is in `cpus.txt` of the run directory and exported as `count_pinned_cpus`), with its memory bound to the slot's node and a step record (`step_count.json`, exported as `count_node`, `_huge_fraction`, `_foreign_busy_cores`, …; see the protocol), and records wall and CPU seconds (`count_wall_s`, `count_cpu_s`), the KMC stages, peak RSS, major page faults and file-system reads (`major_pagefaults`, `fs_inputs` in 512-byte blocks, `read_gb`). Stage 1 reads the input, so `read_gb` > 0 means a cold input and a slower stage 1; stage 2 and `count_cpu_s` are the cache-independent numbers. `default` is repeated like the other schemes and interleaved with them.

## Notes

- KMC3's default signature length is `-p9`; KMC 3.0.0 used 7. Pass `--signature-len 7` for a closer match (there is no GreedyMini table for k = 7).
- Signature order flags of the fork: `kmc -co<table.gm> -ol<label> [-os<seed>]` (binary `.gm`; the seed breaks rank ties at runtime), `kmc -co<table_dna.gm> -ol<label>` (`_dna.gm`; seed in the file header), `-or xor|xor-pcg-hash -ol<label> [-os<seed>]` (xor, or xor then PCG hash). `kmc.json` reports `signature_order_mode` = `xor` / `xor-pcg-hash` / `greedymini` and, for tables, `signature_order_gm_lut` = `dna` (4^k DNA table) | `binary` (2^k rank table; `-oat1` flips the polarity of Φ). Without these flags KMC3 uses its default order. Sidecars `{output}.kmc_co` / `.kmc_col` store mode, seed and table.
- **`-S lex`**: KMC3's `default` order is not pure canonical lexicographic: `CMmer::is_allowed()` (`kmc_api/mmer.h`) excludes poly-A/poly-T-type motifs from ever being selected (a denylist from the KMC2 paper, to bound the largest bin; §2.2 there). `-S lex` runs canonical lexicographic order with no denylist, to isolate that effect. It is run as `-or xor -os0` (order-equivalent to `MIN(kmer, revcomp(kmer))` with no exclusions) but reported as its own `lex` scheme. It is deterministic (seed 0; `-n` repeats it for noise statistics) and not part of `-S all`.
- ENA publishes G. gallus / M. balbisiana as paired `_1`/`_2` FASTQ; `data/datasets.csv` lists those files.
