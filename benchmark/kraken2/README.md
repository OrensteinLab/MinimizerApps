# Benchmark harness (fork-only)

This is the lab's benchmark harness for kraken2, kept in the monorepo's
top-level `benchmark/` tree (outside `kraken2/`, which stays upstream-adjacent)
so future `kraken2` upstream syncs never touch it: host setup, data prep, and
one benchmark that times every minimizer scheme and scores them on labeled
and real reads. Upstream Kraken 2 wrappers stay in `kraken2/scripts/` (`k2`,
`kraken2-build`, …).

`run_benchmark.sh` is the driver. `run_eval.sh` is an alias.

## Libraries

Shared FASTA lives under `db/<id>/`. Scheme indexes and results live under
`runs/<id>/`. The four ids:

| Id | Flag | Shared library | Fraction |
|----|------|----------------|----------|
| `standard` | default or `-N` | `db/standard` | 100% NCBI |
| `standard_p10` | `-m` | `db/standard_p10` | ~10% (`SAMPLE_FRAC=0.1`) |
| `standard_p0.5` | `-z` | `db/standard_p0.5` | ~0.5% (`SAMPLE_SMOKE_FRAC=0.005`) |
| `strain_exclusion` | `-E` | `db/strain_exclusion` | Wood 2019 held-out genomes |

`-E` cannot be combined with `-m`/`-z`/`-N`.

Wood et al. used **one** library for accuracy: the strain-exclusion set
(50 Table S3 taxa held out). Both query sets are classified against that
library (`run_benchmark.sh -E`):

| Dataset | Path | Library | Metrics |
|---|---|---|---|
| **strain-exclusion** | `data/strain-exclusion/` | `-E` only | genus/species sensitivity / PPV / F1 (Table S1) |
| **fda-argos** | `data/fda-argos/` | `-E` only | genus (a) concordant (b) discordant (c) ancestor (d) unclassified (Table S5) — **no** sensitivity/PPV |

NCBI `-m`/`-z`/`-N` are a fork speed experiment on an unlabeled corpus
(`data/speed/` or legacy `data/benchmark_R*.fastq`). They do **not** score
Mason or ARGOS (classifying held-out genomes against full NCBI is the
optimistic case the paper avoided; ARGOS concordance in the paper also used
the exclusion DBs).

**Strain-exclusion** (Wood et al. 2019) is that *library*: 50 Table S3 taxa
held out, then Mason reads from those genomes. Scripts named `prepare_strex_*`
mean **str**ain-**ex**clusion. FASTA+taxonomy go to `db/strain_exclusion/`;
`run_benchmark.sh -E` writes indexes under `runs/strain_exclusion/`.

Accuracy (`evaluate_calls.py`): TP = true rank taxon or descendant; VP =
ancestor; FP = off-lineage; FN = unclassified. Strain-exclusion reports
sensitivity = TP/(TP+VP+FN+FP) and PPV = TP/(TP+FP). Table S1 Kraken 2 genus
targets: prokaryote F1 ≈ 0.966, viral F1 ≈ 0.754.

**Classify flags for paper accuracy (`-E`):** `k2 classify` defaults to
`--minimum-hit-groups 2` (Kraken 2 after 2.0.8-beta). Wood et al. / 2.0.8-beta
`kraken2` did not apply that threshold. `run_benchmark.sh` passes
`--minimum-hit-groups 0` automatically for `-E` so Table S1/S5 scores
match the supplementary tables. For manual runs on the 208 index or Mason reads,
always add `--minimum-hit-groups 0`.

FDA-ARGOS is scored only as Table S5 genus concordance on SRA-assigned taxa
(percent of fragments: concordant, discordant, ancestor, unclassified).
Sensitivity and PPV are not reported for these real reads. Headline numbers
are the unweighted mean of the 25 SRR experiments; per-run rows and a pooled
fragment row are in `tables/accuracy_fda_argos.csv`.

Official extras: [Zenodo 3365797](https://doi.org/10.5281/zenodo.3365797),
[kraken2-experiment-code](https://github.com/DerrickWood/kraken2-experiment-code).
Tables S3/S5: `benchmark/kraken2/tables/`.

## Prerequisites

- Linux; `g++` with OpenMP and `zlib` (`./benchmark/kraken2/setup.sh` installs via apt)
- Conda env with [bioconda `kraken2`](https://anaconda.org/channels/bioconda/packages/kraken2/overview)
  (`setup.sh` can install Miniconda and create `kraken_env` from `environment.yml`)
- Disk for NCBI `db/standard` (tens of GB) and/or the Zenodo strain-exclusion tarball (~9 GB)
- Rank tables under `../../orders/GreedyMini/` (in-repo; from
  [GreedyMini](https://github.com/OrensteinLab/GreedyMini))

`make all` builds this fork’s binaries. Bioconda is the toolchain; repo
`scripts/` and `src/` take precedence on `PATH` via `benchmark_setup.sh`.

## Quick start

```bash
./benchmark/kraken2/setup.sh
conda activate kraken_env
make -j"$(nproc)" all

cp benchmark/kraken2/.env.example benchmark/kraken2/.env
# set KRAKEN2_BENCHMARK_DIR, CONDA_ENV
source benchmark/kraken2/.env

make -C src test
make -C src all ENABLE_MINIMIZER_PROF=1

# Strain-exclusion Mason + Table S5 ARGOS reads
./benchmark/kraken2/download_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"

# NCBI libraries + standard_p10 / standard_p0.5
./benchmark/kraken2/download_library.sh -o "${KRAKEN2_BENCHMARK_DIR}"

# Optional: strain-exclusion FASTA → db/strain_exclusion
./benchmark/kraken2/prepare_strex_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"

# 10% NCBI library, unlabeled speed FASTQs only (not Tables S1/S5)
./benchmark/kraken2/run_benchmark.sh -o "${KRAKEN2_BENCHMARK_DIR}" -m -S Random GreedyMini -n1 -r1

# Paper: held-out library × Mason + FDA-ARGOS
./benchmark/kraken2/run_benchmark.sh -o "${KRAKEN2_BENCHMARK_DIR}" -E -S Random GreedyMini -n1 -r1
```

`-a` skips FDA-ARGOS. `-c` re-summarizes an existing run (auto-picks the only
tree with a manifest, or pass `-m`/`-z`/`-N`/`-E`).

`./benchmark/kraken2/run_benchmark.sh -m -n8` → `runs/standard_p10`.  
`./benchmark/kraken2/run_benchmark.sh -z -S Random -n1 -r1` → `runs/standard_p0.5`.

`db/standard_p10` and `db/standard_p0.5` are created by `download_library.sh` unless `-S`.
Rebuild only smoke from an existing `db/standard`:  
`./benchmark/kraken2/download_library.sh -Z -o "${KRAKEN2_BENCHMARK_DIR}"`.

## Directory layout (`KRAKEN2_BENCHMARK_DIR`)

```
data/
  strain-exclusion/   # Mason PE + truth.tsv  (download_data.sh)
    R1.fastq  R2.fastq  truth.tsv
  fda-argos/          # Table S5 SRRs         (download_data.sh)
    R1.fastq  R2.fastq  [truth.tsv]
db/
  standard/           # full NCBI library
  strain_exclusion/   # Wood 2019 held-out FASTA (prepare_strex_data.sh)
  standard_p10/       # ~10%
  standard_p0.5/      # ~0.5%
runs/
  {standard,strain_exclusion,standard_p10,standard_p0.5}/
    {scheme}/{seed}/{db,results}
    benchmark_comparison.txt|.xlsx
    tables/
    figures/
```

- `tables/*.csv` — `contents`, `scheme_summary`, `vs_ref`, `per_seed`,
  plus `accuracy_strain_exclusion_summary` /
  `accuracy_strain_exclusion` and `accuracy_fda_argos_summary` /
  `accuracy_fda_argos` when those queries ran. XLSX tabs: `0 contents`,
  `1 timing summary`, `2 vs Random`, `3 timing per seed`,
  `4 strain-exclusion summary`, `5 strain-exclusion per seed`,
  `6 FDA-ARGOS summary`, `7 FDA-ARGOS per experiment`.
- `figures/` — histograms prefixed `strain-exclusion_` / `fda-argos_` when those
  sets ran; `fda-argos_genus_concordance.png` is the four Table S5 percentages

Existing `data/benchmark_R1.fastq` is still used as the NCBI unlabeled speed
set. Prefer `data/speed/R1.fastq`. Older trees (`sample_db/`, top-level `p10/`,
`db/` with `library/` directly inside, `db_strain_exclusion/`) are still
used if the new path is missing.

## Schemes (`-S`)

| Label | CLI |
|-------|-----|
| `Random` | `--sketching-scheme Minimizer --order RandomXor` |
| `RandomPcg` | `--sketching-scheme Minimizer --order RandomPcg` |
| `GreedyMini` | `--sketching-scheme Minimizer --order GreedyMini` + `w5_k15_*.gm` |
| `ModRandom` | `--sketching-scheme Mod --order RandomXor` |
| `ModGreedyMini` | `--sketching-scheme Mod --order GreedyMini` + `w5_k6_*.gm` |

GreedyMini tables (`../../orders/GreedyMini/`): `-g` / `GM_FILE` (`w5_k15_canonical.gm`),
`-G` / `GM_DNA_FILE`, `-M` / `MOD_GM_FILE` (`w5_k6_canonical.gm`).

## Script map

| File | Role |
|---|---|
| `setup.sh` | apt + conda env |
| `.env` | `KRAKEN2_BENCHMARK_DIR`, threads |
| `benchmark_setup.sh` | conda + `scripts/k2` on PATH |
| `benchmark_interface.sh` | `-S` schemes, `-m/-z/-N` NCBI size, `-E` strain-exclusion library |
| `xlsxutil.py` | read/write CSV + XLSX (no openpyxl) |
| `create_sample_db.py` | NCBI `library.fna` subsample → `db/standard_p10` / `db/standard_p0.5` |
| `download_library.sh` | NCBI taxonomy + standard libraries + subsamples |
| `download_data.sh` | `data/strain-exclusion` + Table S5 `data/fda-argos` |
| `prepare_strex_data.sh` | Zenodo FASTA → `db/strain_exclusion/` |
| `prepare_fda_argos_reads.py` | ENA download of Table S5 `sra_run` accessions |
| `run_benchmark.sh` | build + classify + score + CSV/XLSX |
| `bench_const.py` | shared `DEFAULT_SEEDS`, scheme order, density metrics |
| `compare_benchmark.py` | stats workbook |
| `score_kraken.py` / `score_argos.py` | strain-exclusion sens/PPV/F1; ARGOS Table S5 genus (a–d) |
| `benchmark_minimizers.sh` | density microbench (CSV) |

## Minimizer density (`benchmark_minimizers.sh`)

Density microbench uses NCBI FASTAs, not query reads. Default `db/standard` →
`runs/standard`; `-m` / `-z` match `standard_p10` / `standard_p0.5`. Writes
`tables/minimizer_density.csv`.

Per scheme:

- `dens_theory` — Minimizer `2/(w+1)`, Mod `1/w`
- `dens_expected` — i.i.d. random DNA
- `dens_particular` — library FASTAs
- `fac_* = dens_* × (w+1)` — density factor
- `ns/minimizer`, super-mer `sm_mean` / `sm_max` (identical-key run lengths)

```bash
make -C kraken2/src benchmark_minimizers
./benchmark/kraken2/benchmark_minimizers.sh -m -S Random ModRandom
```

Selected-ℓ-mer **binning** (Nyström-Persson Table 2, Discount reads, occupancy
histograms) lives in the sibling repo `../minimizer-binning/` — not in this
fork.

## Notes

- `data/simulator.pl` is not Mason 2. Prefer Zenodo FASTQs.
- FDA-ARGOS runs are the 25 `sra_run` values in `benchmark/kraken2/tables/table_s5.csv`.
- Table S5 paper averages are unweighted per-experiment; `score_argos.py`
  writes that mean plus a pooled-fragment row. Per-SRR percentages need
  `truth.tsv` (`download_data.sh` after taxonomy is present).
