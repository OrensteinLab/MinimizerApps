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
| **cami** | `data/cami/` | `-E` (or `-Q` with `-m`/`-N`) | genus/species sens / PPV / F1 over reads whose true taxon is in the library, plus the excluded fraction (see [CAMI II marine](#cami-ii-marine-datacami)) |
| **mock** | `data/mock/` | `-E` (or `-Q`) | profile level: taxon precision / recall, L1, Bray-Curtis vs the known composition (see [mock community](#mock-community-datamock)) |

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
fragment row are in `kraken2_fda_argos_runs.csv` (sheet `fda_argos_runs`).

## Extra query sets (CAMI and mock community)

The two Wood et al. sets resemble each other (simulated reads from NCBI
genomes; isolate SRA runs). Two more query sets test the minimizer schemes on
other kinds of reads. Both are optional: `run_benchmark.sh` classifies them
only when prepared (`-C` / `-K` skip them; `-Q` adds them to a `-m`/`-z`/`-N`
run; `-E` includes them automatically). The existing sets and flags are
unchanged.

| Dataset | Reads | Size on disk | Truth |
|---|---|---|---|
| `data/cami/` | CAMI II marine short-read sample 0, 2x150 bp, 250 000 of 16 647 395 pairs | ~170 MB | per read (CAMI gold standard) |
| `data/mock/` | ZymoBIOMICS D6300 mock, SRR12324251, 2x150 bp NovaSeq, 250 000 of 29 386 307 pairs | ~175 MB | composition only |

Subsampling is systematic (every k-th pair, k = pairs in the run // 250 000,
no RNG), so every machine gets the same reads.

### CAMI II marine (`data/cami/`)

```bash
./benchmark/kraken2/prepare_cami_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"   # ~6 GB download, deleted after sampling
```

- Source: CAMI II toy marine dataset (Meyer et al. 2022, *Nat Methods*
  19:429, doi:10.1038/s41592-022-01431-4), PUBLISSO repository
  <https://frl.publisso.de/data/frl:6425521/marine/short_read/>:
  `marmgCAMI2_sample_0_reads.tar.gz` (md5 `3af16716a9613490d3a57688e2fd51e1`)
  and `marmgCAMI2_setup.tar.gz` (md5 `52e0f2834ae826568b089db40347a59e`).
  Real CAMI data was used (no CAMISIM fallback needed): a complex community of
  977 genomes with differential abundances, new strains / species / genera, plasmids
  and viruses; simulated with CAMISIM/ART (profile `mbarc`, 270 bp fragments).
- Files: `R1.fastq R2.fastq` (original IDs `S0R<i>/1`, `/2`), `truth.tsv`
  (`#seqid taxid domain`, one row per pair; kraken2 reports the pair as
  `S0R<i>`, `score_kraken.py` strips `/1`), `truth_extra.tsv` (+ genome id and
  CAMI `novelty_category`), `taxonomic_profile.txt` (CAMI gold-standard profile of
  the whole sample), `SOURCE.txt`.
- Truth: `reads_mapping.tsv.gz` `tax_id` of the source genome (NCBI taxdump
  2017-02). New-species / new-genus genomes carry a genus / family taxid, so
  they have no species / genus node: they are counted as `no` at that rank.
  Plasmid (45202) and unclassified (32644) reads likewise.
- Scoring (`score_cami.py`, library-aware): at each rank every read is `no`
  (no node at the rank), `excluded` (true taxon at the rank has no sequence in
  the run's library, taken from the library's `prelim_map*.txt`) or scored
  (TP/VP/FP/FN as in `evaluate_calls.py`). Sensitivity / PPV / F1 use scored
  reads only; excluded reads are reported with what the classifier did with them
  (`ex_unclassified`, `ex_consistent` = call on the true lineage, `ex_off` =
  spurious call). Groups: `all`, domain, `novelty:<CAMI category>`.
  On `data/strain-exclusion` this scorer reproduces `score_kraken.py` exactly.
- **Library.** Share of reads (with a genus / species node) whose true taxon is
  absent from the library, measured on this subsample with
  `score_cami.py --nodes <db>/taxonomy/nodes.dmp --truth data/cami/truth.tsv --library-dir <db>`:

  | Library | Genus excluded | Species excluded |
  |---|---|---|
  | `strain_exclusion` (`-E`) | 9.6 % | 22.6 % |
  | `standard` (`-N`) | 0.2 % | 8.0 % |
  | `standard_p10` (`-m`) | 41.9 % | 74.9 % |
  | `standard_p0.5` (`-z`) | 66.3 % | 96.7 % |

  Recommended: `-E` (same library and `--minimum-hit-groups 0` as the other
  query sets, 90 % of reads keep a reference at genus level). Full `standard`
  is the best-covered but needs the 741 GB DB (`-N -Q`). The subsampled
  `standard_p10` / `standard_p0.5` drop most genera: use them only for timing.
- Caveats: simulated reads (ART error model), 2017 taxonomy vs the library's
  taxonomy (taxids missing from `nodes.dmp` count as `no`), strain-level
  novelty is not scored. Only the marine sample is validated; `-s strain|plant`
  and `-n SAMPLE` select other CAMI II short-read samples
  (`strmgCAMI2` / `rhimgCAMI2`; plant has host reads).

### Mock community (`data/mock/`)

```bash
./benchmark/kraken2/prepare_mock_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"   # streams ~3.5 GB, keeps the sample
```

- Source: ZymoBIOMICS Microbial Community Standard (D6300: 8 bacteria + 2
  yeasts), BioProject PRJNA648136, run **SRR12324251** (ENA/SRA), Illumina
  NovaSeq 6000, 2x150 bp, 29 386 307 pairs, shotgun (`WGS`, `RANDOM`).
- Expected composition: Zymo D6300 product sheet Table 1, stored in
  `tables/mock_zymo_d6300_expected.csv` (copied to `data/mock/composition.csv`).
  Reads scale with genomic-DNA mass, so `gdna_pct` (12 % per bacterium, 2 % per
  yeast) is the default expected abundance; `genome_copy_pct` and `cell_pct` are
  available with `score_mock.py --abundance`.
- Scoring (`score_mock.py`, needs no DB, works from the `kraken2 --report`):
  genus / species precision, recall and F1 of the detected taxa (a taxon is
  detected at clade reads / reads >= `min_frac`; 0, 1e-4 (headline), 1e-3),
  L1 and Bray-Curtis between expected and observed relative abundances,
  classified fraction, on-target / off-target / above-rank / unclassified
  fractions (sum to 1) and the top off-target taxa. Scopes: `all` and
  `bacteria` (yeast removed from expected and from the read denominator; the
  Kraken standard libraries have no fungi, so `all` recall is capped at 0.8).
  Unit tests: `python3 -m unittest discover -s benchmark/kraken2/tests`.
- Caveats: the lot-specific composition deviates from theory by up to ~15 %
  relative; extraction/library-prep bias is not removed (the sampled reads
  have mean GC 53.6 %, vs ~46 % expected for the nominal DNA mix, so this run
  is GC-rich-biased; the bias is the same for every scheme, but absolute L1 /
  Bray-Curtis include it, compare schemes rather than reading them as error); E. coli and Shigella,
  Salmonella and other close relatives make species-level off-target calls
  expected; D6300 strains and reference genomes were updated in lot ZRC190633.
- Library: all eight bacterial species have sequence in every library
  (`strain_exclusion`, `standard`, `standard_p10`, `standard_p0.5`); the two
  yeasts are in none (no fungi), so their 4 % of the DNA can only be
  unclassified or off-target. `-E` recommended (same library as the other sets).

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

# Optional extra query sets (CAMI II marine, ZymoBIOMICS D6300 mock)
./benchmark/kraken2/prepare_cami_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"
./benchmark/kraken2/prepare_mock_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"

# NCBI libraries + standard_p10 / standard_p0.5
./benchmark/kraken2/download_library.sh -o "${KRAKEN2_BENCHMARK_DIR}"

# Optional: strain-exclusion FASTA → db/strain_exclusion
./benchmark/kraken2/prepare_strex_data.sh -o "${KRAKEN2_BENCHMARK_DIR}"

# 10% NCBI library, unlabeled speed FASTQs only (not Tables S1/S5)
./benchmark/kraken2/run_benchmark.sh -o "${KRAKEN2_BENCHMARK_DIR}" -m -S Random GreedyMini -n1 -r1

# Paper: held-out library × Mason + FDA-ARGOS
./benchmark/kraken2/run_benchmark.sh -o "${KRAKEN2_BENCHMARK_DIR}" -E -S Random GreedyMini -n1 -r1
```

`-a` skips FDA-ARGOS, `-C` skips CAMI, `-K` skips the mock community; `-Q` adds
CAMI and the mock to a `-m`/`-z`/`-N` run. `-c` re-summarizes an existing run (auto-picks the only
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
  cami/               # CAMI II marine        (prepare_cami_data.sh)
    R1.fastq  R2.fastq  truth.tsv  truth_extra.tsv  taxonomic_profile.txt  SOURCE.txt
  mock/               # ZymoBIOMICS D6300     (prepare_mock_data.sh)
    R1.fastq  R2.fastq  composition.csv  SOURCE.txt
db/
  standard/           # full NCBI library
  strain_exclusion/   # Wood 2019 held-out FASTA (prepare_strex_data.sh)
  standard_p10/       # ~10%
  standard_p0.5/      # ~0.5%
runs/
  {standard,strain_exclusion,standard_p10,standard_p0.5}/
    {scheme}/{seed}/{db,results}
    benchmark_comparison.txt|.xlsx   # this library only (same sheets)
    tables/                          # kraken2_<sheet>.csv of this library
    figures/
results/                             # all libraries merged
  kraken2_results.xlsx
  kraken2_<sheet>.csv
```

- **Results workbook** (same layout as the sshash / kmc apps: `export_layout.py`,
  helpers in `../common/xlsx_book.py`). Each `compare_benchmark.py` call (the
  one at the end of `run_benchmark.sh`, or `-c`) writes this library's
  `kraken2_<sheet>.csv` under `runs/<id>/tables/` plus
  `runs/<id>/benchmark_comparison.xlsx`, then merges every library that has
  been exported into `<base>/results/kraken2_results.xlsx` and
  `<base>/results/kraken2_<sheet>.csv` (also runnable alone:
  `python3 benchmark/kraken2/export_layout.py --base-dir <base>`). One CSV per
  sheet, same rows and columns as the sheet. Sheets, in order:

  | Sheet | Rows |
  |---|---|
  | `overview` | column reference: what every sheet and column means, by category |
  | `libraries` | one per Kraken library: sequences, taxa, FASTA GB, role, query sets and orders run |
  | `datasets` | one per query set: reads, read length, truth type, provenance, libraries |
  | `orders` | one per scheme/order: CLI, l / k / w, spaced seed, expected and per-library particular density factor |
  | `summary` | library × dataset × order: mean / std over seeds (std rounded to ≥ 4 significant digits) and `_ratio` to Random of the same library and dataset (green = better, red = worse); build, classify, density, NextMinimizer, vs-Random agreement and the headline accuracy of the dataset (`se_*` strain-exclusion, `argos_*`, `cami_*`, `mock_*`) |
  | `runs` | one per seed × library × dataset × order: the raw per-run numbers (library / dataset / order cells merged, a border after each library × dataset block) |
  | `vs_reference` | paired agreement with Random (MAPE, concordance, Jaccard, new classified / unclassified) |
  | `timing_correlation` | Pearson r over orders of density factor vs build / classify time and RSS |
  | `classify_repeats` | every timed classify repeat |
  | `strain_exclusion_summary` / `_runs`, `fda_argos_summary` / `_runs`, `cami_summary` / `_runs`, `mock_summary` / `_runs` | accuracy detail (domain × rank, per-SRR, per group, per detection threshold); the `_summary` sheets carry the charts (metric vs scheme, metric vs classify time), `summary` carries the timing / density charts, per library |

  Detail sheets appear only for the query sets that ran. These replace the
  former `tables/{contents,query_sets,scheme_summary,per_seed,vs_ref,
  classify_repeats,accuracy_*}.csv` and the numbered `0 contents` … `11 mock
  per seed` tabs; `compare_benchmark.py` deletes those legacy CSV names from
  `tables/` when it writes the new ones. `tables/density_per_seed.csv` /
  `minimizer_density.csv` (density microbenchmark) are inputs and are kept.
- `figures/` — histograms prefixed `strain-exclusion_` / `fda-argos_` / `cami_` / `mock_` when those
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
| `export_layout.py` | results workbook + `kraken2_<sheet>.csv` (sheets above); `../common/xlsx_book.py` / `xlsx_charts.py` / `xlsx_util.py` write and read .xlsx without openpyxl, `../common/xlsx_check.py` validates the package |
| `create_sample_db.py` | NCBI `library.fna` subsample → `db/standard_p10` / `db/standard_p0.5` |
| `download_library.sh` | NCBI taxonomy + standard libraries + subsamples |
| `download_data.sh` | `data/strain-exclusion` + Table S5 `data/fda-argos` |
| `prepare_strex_data.sh` | Zenodo FASTA → `db/strain_exclusion/` |
| `prepare_cami_data.sh` / `prepare_cami_reads.py` | CAMI II marine sample → `data/cami/` (subsample + truth) |
| `prepare_mock_data.sh` / `prepare_mock_reads.py` | ZymoBIOMICS D6300 SRR12324251 → `data/mock/` |
| `score_cami.py` / `score_mock.py` | CAMI library-aware sens/PPV/F1; mock profile precision/recall/L1/Bray-Curtis |
| `tests/` | `python3 -m unittest discover -s benchmark/kraken2/tests` (synthetic data, no DB) |
| `prepare_fda_argos_reads.py` | ENA download of Table S5 `sra_run` accessions |
| `run_benchmark.sh` | build + classify + score + CSV/XLSX |
| `bench_const.py` | shared `DEFAULT_SEEDS`, scheme order, density metrics |
| `compare_benchmark.py` | statistics text report + per-run numbers fed to `export_layout.py` |
| `score_kraken.py` / `score_argos.py` | strain-exclusion sens/PPV/F1; ARGOS Table S5 genus (a–d) |
| `benchmark_minimizers.sh` | density microbench (CSV) |
| `density_sweep.sh` / `density_summary.py` | density factor vs ℓ over several random streams (calls `benchmark_minimizers`), mean ± std table |

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

Density **sweeps** over ℓ and streams (mean ± std per ℓ, old vs new revision,
`Random` baselines) use the same tool through `density_sweep.sh` +
`density_summary.py`; the exact commands and the measured GreedyMini
prefix-vs-centre tables are in
[`../../kraken2/minimizers/README.md`](../../kraken2/minimizers/README.md#density).
Raw sweep rows: `tables/gm_centre_density_*.csv`.

## Notes

- `data/simulator.pl` is not Mason 2. Prefer Zenodo FASTQs.
- FDA-ARGOS runs are the 25 `sra_run` values in `benchmark/kraken2/tables/table_s5.csv`.
- Table S5 paper averages are unweighted per-experiment; `score_argos.py`
  writes that mean plus a pooled-fragment row. Per-SRR percentages need
  `truth.tsv` (`download_data.sh` after taxonomy is present).
