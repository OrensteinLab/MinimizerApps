# Kraken2 benchmark

Host setup, data preparation and one benchmark that builds a Kraken2 database
per scheme and seed, times the build and the classification, and scores the
calls on simulated and real reads. The benchmark scripts live outside `kraken2/`, so
upstream syncs of the fork never touch it. `run_benchmark.sh` is the driver;
`experiments/full_run.sh` runs the paper's experiments.

**Notation.** *k* = length of the selected *k*-mer, *w* = number of *k*-mers in a
window, *L* = *k* + *w* − 1, *t* = *t*-mer length of the mod-minimizer. Kraken2's own
names: its minimizer length ℓ (`-l`, 31) is *k*, its *k*-mer length (`-k`, 35) is
*L*, so *w* = 5. Kraken2's flags and the columns written by the code keep their
names (`distinct_minimizers`, `--minimum-hit-groups`, ...).

## Libraries

Shared library FASTA lives under `db/<id>/`; scheme databases and results under
`results/runs/<id>/`.

| Id | Flag | Library | Use |
|----|------|---------|-----|
| `strain_exclusion` | `-E` | Wood et al. 2019 strain-exclusion library (50 Table S3 taxa held out) | the paper's runs (accuracy and timing) |
| `standard` | default or `-N` | 100% NCBI standard | timing only |
| `standard_p10` | `-m` | ~10% of `standard` (`SAMPLE_FRAC=0.1`) | timing only |
| `standard_p0.5` | `-z` | ~0.5% of `standard` (`SAMPLE_SMOKE_FRAC=0.005`) | smoke tests |

`-E` cannot be combined with `-m`/`-z`/`-N`. The NCBI libraries classify only an
unlabeled speed set (paired FASTQs you put in `data/speed/R1.fastq`, `R2.fastq`), never the
labeled sets: classifying held-out genomes against full NCBI is the optimistic
case Wood et al. avoided.

## Query sets

All three are classified against the strain-exclusion library (`-E`):

| Dataset | Path | Metrics |
|---|---|---|
| strain-exclusion | `data/strain-exclusion/` | genus / species sensitivity, PPV, F1 (Table S1) |
| fda-argos | `data/fda-argos/` | genus concordant / discordant / ancestor / unclassified (Table S5); no sensitivity / PPV |
| cami (not in Wood et al.) | `data/cami/` | genus / species sensitivity, PPV, F1 over reads whose true taxon is in the library, plus the excluded fraction |

`-a` skips FDA-ARGOS and `-C` skips CAMI; `-Q` adds CAMI to a `-m`/`-z`/`-N`
run. CAMI is classified only when `data/cami/` has been prepared.

**Strain-exclusion** (`prepare_strex_*` = **str**ain-**ex**clusion): Mason
reads from the 50 held-out genomes, classified against the library without
them. Accuracy (`evaluate_calls.py`, `score_kraken.py`): TP = call is the true
taxon at the rank or a descendant; VP = ancestor; FP = off-lineage; FN =
unclassified; sensitivity = TP/(TP+VP+FN+FP), PPV = TP/(TP+FP). Table S1
Kraken 2 genus targets: prokaryote F1 ≈ 0.966, viral F1 ≈ 0.754.

**FDA-ARGOS**: the 25 `sra_run` values of `tables/table_s5.csv`, scored as
Table S5 genus concordance (percent of fragments). Headline numbers are the
unweighted mean over the 25 runs, as in the paper; per-run rows and a pooled
row are on the `accuracy_fda_argos` sheet (`scope` = `unweighted_mean` /
`pooled` / `experiment`).

**Hit-group settings.** `k2 classify` defaults to `--minimum-hit-groups 2`;
Wood et al. (Kraken 2.0.8-beta) applied no threshold, so
`--minimum-hit-groups 0` is the paper setting. `run_benchmark.sh` classifies
every (seed, dataset, scheme) cell once with each; output files carry the
suffix `_mhg0` / `_mhg2`. For manual `k2 classify` runs that should match the
paper, add `--minimum-hit-groups 0`.

### CAMI II marine (`data/cami/`)

CAMI II marine short-read sample 0, 2x150 bp, 250 000 of 16 647 395 pairs
(~170 MB), with per-read truth from the CAMI gold standard. It tests the
schemes on a different kind of reads than the two Wood et al. sets.

```bash
benchmark/kraken2/prepare_cami_data.sh -o "$KRAKEN2_BENCH_ROOT"   # ~6 GB download, deleted after sampling
```

- **Source:** CAMI II toy marine dataset (Meyer et al. 2022, *Nat Methods*
  19:429, doi:10.1038/s41592-022-01431-4), PUBLISSO
  <https://frl.publisso.de/data/frl:6425521/marine/short_read/>:
  `marmgCAMI2_sample_0_reads.tar.gz` (md5 `3af16716a9613490d3a57688e2fd51e1`)
  and `marmgCAMI2_setup.tar.gz` (md5 `52e0f2834ae826568b089db40347a59e`). 977
  genomes with differential abundances, new strains / species / genera,
  plasmids and viruses; simulated with CAMISIM/ART (profile `mbarc`, 270 bp
  fragments).
- **Subsampling:** systematic (every *k*-th pair, *k* = pairs in the run // 250 000,
  no RNG), so every machine gets the same reads.
- **Files:** `R1.fastq R2.fastq` (IDs `S0R<i>/1`, `/2`), `truth.tsv`
  (`#seqid taxid domain`, one row per pair), `truth_extra.tsv` (+ genome id and
  CAMI `novelty_category`), `taxonomic_profile.txt`, `SOURCE.txt`.
- **Truth:** the `tax_id` of the source genome (`reads_mapping.tsv.gz`, NCBI
  taxdump 2017-02). Reads with no node at a rank (new-species / new-genus
  genomes, plasmid and unclassified reads) count as `no` at that rank.
- **Scoring** (`score_cami.py`, library-aware): at each rank a read is `no`,
  `excluded` (its true taxon has no sequence in the library, from the
  library's `prelim_map*.txt`) or scored (TP/VP/FP/FN as above). Sensitivity /
  PPV / F1 use scored reads only; for excluded reads the scorer reports
  `ex_unclassified`, `ex_consistent` (call on the true lineage) and `ex_off`
  (off-lineage call). Groups: `all`, domain, `novelty:<CAMI category>`. On
  `data/strain-exclusion` it reproduces `score_kraken.py` exactly.
- **Library coverage** (share of reads with a node at the rank whose true taxon
  is absent from the library; `score_cami.py --nodes <db>/taxonomy/nodes.dmp
  --truth data/cami/truth.tsv --library-dir <db>`):

  | Library | Genus excluded | Species excluded |
  |---|---|---|
  | `strain_exclusion` (`-E`) | 9.6 % | 22.6 % |
  | `standard` (`-N`) | 0.2 % | 8.0 % |
  | `standard_p10` (`-m`) | 41.9 % | 74.9 % |
  | `standard_p0.5` (`-z`) | 66.3 % | 96.7 % |

  CAMI sensitivities are therefore over a library-covered subset and are not
  directly comparable to the paper's definition over all reads.
- **Caveats:** simulated reads (ART error model); 2017 taxonomy vs the library's
  taxonomy (taxids missing from `nodes.dmp` count as `no`); strain-level
  novelty is not scored. Only the marine sample is validated (`-s strain|plant`
  and `-n SAMPLE` select other CAMI II short-read samples).

Wood et al. extras: [Zenodo 3365797](https://doi.org/10.5281/zenodo.3365797),
[kraken2-experiment-code](https://github.com/DerrickWood/kraken2-experiment-code).
Tables S1, S3, S5: `benchmark/kraken2/tables/`. The Zenodo and SRA downloads
have no published checksums.

## Prerequisites

- Linux, `g++` with OpenMP, zlib, Python >= 3.9; `curl` or `wget`, `tar` and
  `md5sum` for the downloads.
- The `kraken_env` conda env with
  [bioconda `kraken2`](https://anaconda.org/channels/bioconda/packages/kraken2/overview)
  (`setup.sh` installs the apt packages, Miniconda if needed, and the env from
  `kraken2/environment.yml`). Bioconda supplies the library and taxonomy
  tools; `benchmark_setup.sh` puts the fork's `kraken2/scripts/` and
  `kraken2/src/` first on `PATH`.
- The GreedyMini tables `orders/GreedyMini/w5_k15_canonical.gm` and
  `w5_k6_canonical.gm` (built from the `.txt` rank lists, see the root README).
- Disk for the paper's runs: ~43 GB for the Zenodo download (34 GB), the
  strain-exclusion library and taxonomy, and the reads, plus ~290 GB for
  the 35 databases (up to 10 GB each; all are kept, the export reads them). The optional full NCBI
  `db/standard` takes 795 GB (727 GB library, 68 GB taxonomy).

## Quick start

From the repo root:

```bash
benchmark/kraken2/setup.sh
conda activate kraken_env
make -C kraken2 all && make -C kraken2 test

cp benchmark/kraken2/.env.example benchmark/kraken2/.env
# set KRAKEN2_BENCH_ROOT and CONDA_ENV, then:
source benchmark/kraken2/.env

# Strain-exclusion library + taxonomy -> db/strain_exclusion (needed by -E and
# by download_data.sh, whose ARGOS truth uses the taxonomy)
benchmark/kraken2/prepare_strex_data.sh -o "$KRAKEN2_BENCH_ROOT"
# Strain-exclusion Mason reads + the 25 FDA-ARGOS runs -> data/
benchmark/kraken2/download_data.sh -o "$KRAKEN2_BENCH_ROOT"
# CAMI II marine -> data/cami (part of the paper's runs)
benchmark/kraken2/prepare_cami_data.sh -o "$KRAKEN2_BENCH_ROOT"

# The paper's runs: run_benchmark.sh -E -n 5 -t 16 -j 1 -P at0-cg1,at1-cg0
# -S <all five schemes>, into $KRAKEN2_BENCH_ROOT/$RESULTS_NAME (default results)
RESULTS_NAME=results benchmark/kraken2/experiments/full_run.sh

# A quick check: two schemes, one seed
benchmark/kraken2/run_benchmark.sh -o "$KRAKEN2_BENCH_ROOT" -E -S xor greedymini -n 1
```

`run_benchmark.sh`'s defaults are the paper's: all five schemes, both
polarities of Φ (`-P at0-cg1,at1-cg0`), 5 seeds (`-n`), 16 build threads
(`-t`), one build at a time (`-j`); only the library flag `-E` must be given.
`full_run.sh` resumes when re-run (databases with a `hash.k2d` are not
rebuilt). `-c` re-exports an existing run (pass the library flag, or it picks
the only tree with a manifest). The classify profile (the `profile` sheet) is
a separate run on a separate instrumented build:
`benchmark/kraken2/experiments/profile_run.sh`.

Optional, not used by the paper: the NCBI standard library and its subsamples,

```bash
benchmark/kraken2/download_library.sh -o "$KRAKEN2_BENCH_ROOT"      # db/standard, standard_p10, standard_p0.5
benchmark/kraken2/download_library.sh -Z -o "$KRAKEN2_BENCH_ROOT"   # only re-make standard_p0.5 from db/standard
benchmark/kraken2/run_benchmark.sh -o "$KRAKEN2_BENCH_ROOT" -m -S xor greedymini -n 1   # speed set on standard_p10
```

## Directory layout (`KRAKEN2_BENCH_ROOT`)

```
data/
  strain-exclusion/   R1.fastq R2.fastq truth.tsv                  (download_data.sh)
  fda-argos/          R1.fastq R2.fastq truth.tsv                  (download_data.sh)
  cami/               R1.fastq R2.fastq truth.tsv truth_extra.tsv ... (prepare_cami_data.sh)
  speed/              R1.fastq R2.fastq, unlabeled, for -m/-z/-N   (yours)
db/
  strain_exclusion/   library + taxonomy                           (prepare_strex_data.sh)
  standard/  standard_p10/  standard_p0.5/                         (download_library.sh)
zenodo/               Wood et al. 2019 reference tarball
results/
  runs/<library>/
    <scheme dir>/<seed>/{db,results}   # scheme dir: xor | greedymini/k15_at0-cg1 | ...
    .benchmark_run_manifest
    benchmark_comparison.txt           # text summary
    minimizer_density.csv              # density microbench rows (an input of the export)
  csvs/kraken2_<sheet>.csv             # one CSV per sheet, all libraries
  kraken2_results.xlsx                 # the workbook, all libraries
```

Several `run_benchmark.sh -E -S <schemes>` invocations into the same results
dir add up: each replaces only the manifest lines of its own run dirs, and the
export covers every scheme. A scheme / seed whose `db/hash.k2d` exists is not
rebuilt, only classified (`db/` may hold symlinks to a database built
elsewhere, with its `results/build_performance.txt` copied next to it).

## Results workbook

`compare_benchmark.py` (run at the end of `run_benchmark.sh`, or with `-c`)
replaces this library's rows in `results/csvs/` and rewrites
`results/kraken2_results.xlsx` (layout in `export_layout.py`, writer
`benchmark/common/xlsx_book.py`; same layout and column naming as the other
benchmarks). `python3 benchmark/kraken2/export_layout.py --base-dir <ROOT>`
re-renders the workbook from the CSVs. The `overview` sheet is the column
reference; in short:

- Raw columns carry their unit (`build_wall_s`, `classify_max_rss_gb`,
  `classified_pct`, `lookups_per_mbp`). Aggregated sheets add `<metric>_mean`,
  `<metric>_std` (across seeds; replicates are the seeds, one classify each)
  and, for metrics with a direction, `<metric>_ratio` = mean / the `xor` mean
  of the same library, dataset and `minimum_hit_groups` (green = better, red =
  worse). The reference scheme is `xor` (Kraken2's default order).
- Identity columns: `library, dataset, scheme, gm_k, projection, gm_table,
  mod_t, minimum_hit_groups, n_repeats` / `seed_idx`. Classify-side sheets carry
  `minimum_hit_groups` and are sorted dataset → `minimum_hit_groups` → scheme;
  build sheets do not.
- Memory and size columns are decimal GB (`_gb`, 10^9 bytes,
  `benchmark/common/units.py`).
- The classify wall time splits into `classify_load_s` (database load) and
  `classify_run_s` (the reads); compare schemes on `classify_run_s`.

| Sheet | Rows |
|---|---|
| `overview` | what every sheet and column means |
| `libraries` | one per library: sequences, taxa, FASTA GB, query sets and schemes run |
| `datasets` | one per query set: reads, read length, truth, scoring, provenance |
| `schemes` | one per scheme variant: `k`, `w`, `L`, GreedyMini table, `mod_t`, CLI, expected and particular density factors |
| `build` / `build_runs` | library × scheme, mean / std / ratio over seeds / per seed: density factors, `distinct_minimizers` (the `hash.k2d` header), build time and RSS, `index_size_gb` |
| `summary` / `runs` | library × dataset × `minimum_hit_groups` × scheme, mean / std / ratio / per seed: classify time, RSS, throughput, `classified_pct`, accuracy; `runs` also has I/O and pinning telemetry |
| `profile` | classify phase profile (`experiments/profile_run.sh`); only when its files exist |
| `vs_reference` | paired agreement with `xor` of the same seed, dataset and `minimum_hit_groups` |
| `accuracy_strain_exclusion`, `accuracy_fda_argos`, `accuracy_cami` | accuracy per slice (`minimum_hit_groups` × group × rank, or SRA run): mean / std / ratio over seeds |
| `paper_check` | `xor` at `minimum_hit_groups` 0 on `strain_exclusion` next to the published Table S1 / S5 Kraken 2 values (`tables/paper_kraken2_published.csv`), absolute tolerance 0.01 |
| `density_effect`, `density_effect_points` | Pearson r between each metric's ratio to `xor` and the density factor's (`benchmark/common/density_effect.py`) |
| `significance` | paired-by-seed tests against `xor`: *t*-test on log ratios, Bonferroni, TOST at ±5% (`benchmark/common/significance.py`) |

`build` and `summary` end with two rows, `correlation with density` and
`correlation with distinct minimizers` (live `CORREL` formulas over the
`*_ratio` columns; values also in `kraken2_correlations.csv`). Accuracy
columns are per dataset and never averaged across datasets; a row fills only
its own dataset's group.

## Schemes (`-S`)

| Scheme | Run directory (`results/runs/<lib>/<dir>/<seed>/`) | Kraken2 CLI |
|-------|------|-----|
| `xor` | `xor` | `--sketching minimizer --order xor` |
| `xor-pcg-hash` | `xor-pcg-hash` | `--sketching minimizer --order xor-pcg-hash` |
| `greedymini` | `greedymini/k15_<projection>` | `--sketching minimizer --order greedymini --greedymini-file w5_k15_canonical.gm --gm-projection <projection>` |
| `mod-minimizer-xor` | `mod-minimizer-xor/t6` | `--sketching mod-minimizer --order xor` |
| `mod-minimizer-greedymini` | `mod-minimizer-greedymini/k6_<projection>_t6` | `--sketching mod-minimizer --order greedymini --greedymini-file w5_k6_canonical.gm --gm-projection <projection>` |

The names are the schemes shared by all four methods
(`benchmark/common/scheme_names.py`) and are exact (any other spelling is
rejected). Parameters are separate fields (manifest, `minimizer_density.csv`
and sheet columns `gm_k, projection, gm_table, mod_t`):

- `projection` is the polarity of Φ, the GreedyMini DNA-to-binary projection:
  `at0-cg1` (A,T → 0, C,G → 1) or `at1-cg0`. `-P` / env `GM_PROJECTION` takes a
  comma list (default both). Each GreedyMini scheme runs once per polarity, so
  the default run has 7 variants. `--dry-run` lists the variants and which are
  done.
- `mod_t` = *t* = r + ((*k* − r) mod *w*), r = 4: 6 for *k* = 31, *w* = 5.
- GreedyMini tables (`orders/GreedyMini/`): `-g` / env `GM_FILE`
  (`w5_k15_canonical.gm`, `greedymini`), `-M` / env `MOD_GM_FILE`
  (`w5_k6_canonical.gm`, `mod-minimizer-greedymini`).

`mod-minimizer` is the canonical mod-minimizer (`kraken2/src/mod_minimizer.h`):
mod-sampling (Groot Koerkamp & Pibiri, 2024) with the center rule shared with
Minimap2 and SSHash (`common/minimizer-schemes/mod_minimizer_rule.h`): a tie at
the minimum takes the offset closest to the *L*-mer's center, two equally close
ones by the *L*-mer ≤ its reverse complement, so a read and its reverse
complement give the same selected *k*-mers.

## Run output layout (`results/runs/<lib>/<scheme dir>/<seed>/results/`)

One classify per (seed, dataset, scheme, `--minimum-hit-groups` m ∈ {0, 2});
`<ds>` = `strain-exclusion` | `fda-argos` | `cami` | `speed`.

| File | Content |
|------|---------|
| `build_performance.txt` | `k2 build` stderr, then the GNU `time -v` block |
| `<ds>_classification_mhg<m>.txt` | `k2 classify --output` |
| `<ds>_report_mhg<m>.txt` | `k2 classify --report` |
| `classification_performance_<ds>_mhg<m>.txt` | classify stderr (`processed in`, classified %) then GNU `time -v`; DB load = wall − processing |
| `classify_profile_<ds>_mhg<m>.txt` | profiling runs only (`experiments/profile_run.sh`) |
| `{strain_exclusion,fda_argos,cami}_metrics_mhg<m>.tsv` | scorer output |

`results/runs/<lib>/.benchmark_run_manifest` has one tab-separated line per
scheme and seed, written when that seed finished classifying: `scheme`,
`context` (`<lib>/<scheme dir>/<seed>`), results dir, toggle (bit string),
seed index, `gm_k`, `projection`, `gm_table`, `mod_t`.

## Script map

| File | Role |
|---|---|
| `setup.sh` | apt packages + conda env |
| `.env` | `KRAKEN2_BENCH_ROOT` (optional `KRAKEN2_DATA_DIR`, `KRAKEN2_RESULTS_DIR`), threads |
| `benchmark_setup.sh` | conda + the fork's `scripts/` and `src/` on `PATH` |
| `benchmark_interface.sh` | shared options: `-S`, `-g/-M/-P`, `-m/-z/-N/-E` |
| `prepare_strex_data.sh` | Zenodo FASTA + taxonomy → `db/strain_exclusion/` |
| `download_data.sh` | `data/strain-exclusion` + `data/fda-argos` |
| `prepare_cami_data.sh` / `prepare_cami_reads.py` | CAMI II marine → `data/cami/` |
| `prepare_fda_argos_reads.py` | ENA download of the Table S5 runs |
| `download_library.sh` / `create_sample_db.py` | NCBI taxonomy + standard library + subsamples |
| `run_benchmark.sh` | build + classify + score + export |
| `experiments/full_run.sh`, `experiments/profile_run.sh` | the paper's run; the classify profile |
| `bench_const.py` | `DEFAULT_SEEDS`, scheme order, run-dir ids, density metrics |
| `score_kraken.py` / `score_argos.py` / `score_cami.py` | strain-exclusion / FDA-ARGOS / CAMI scoring |
| `compare_benchmark.py` | text report + per-run numbers for `export_layout.py` |
| `export_layout.py` | `results/kraken2_results.xlsx` + `results/csvs/kraken2_<sheet>.csv` |
| `benchmark_minimizers.sh` | density microbench (CSV) |
| `density_sweep.sh` / `density_summary.py` | density factor vs *k* over several random streams |
| `tests/` | `python3 -m unittest discover -s benchmark/kraken2/tests` (synthetic data, no DB) |

## Density microbench (`benchmark_minimizers.sh`)

Measures the density of every scheme at the benchmark's configuration (*k* = 31,
*w* = 5, *L* = 35, minimizer spaces *s* = 7, the same GreedyMini tables as the
databases) and the `NextMinimizer` cost:

- **expected** (`dens_expected`, `fac_expected`): random i.i.d. ACGT (`-N`
  bases, default 100 Mbp; Monte Carlo, varies slightly with `--seed`). The
  export replaces it with the count of `expected_density.py`
  (`benchmark/common/expected_density/`: the same scanner, 10^9 bases), as for
  the other methods.
- **particular** (`dens_particular`, `fac_particular`): the library FASTAs of
  the shared DB, read in full (8.7 GB, ~9.3 Gbp for `strain_exclusion`).
- `fac_* = dens_* × (w+1)` (random minimizer ≈ 2); `dens_theory` = 2/(*w*+1)
  (minimizer) or 1/w (mod-minimizer); `ns_per_minimizer`; super-mer lengths
  `sm_mean_*` / `sm_max_*`.

Seed *i* uses the toggle mask `DEFAULT_SEEDS[i-1]` (`bench_const.py`, passed
as `-T` with `--seed-idx i`). One row per scheme and seed goes into
`results/runs/<lib>/minimizer_density.csv` (columns `scheme, gm_k, projection,
gm_table, mod_t, seed_idx, toggle_hex, k, w, L, spaces, dens_*, fac_*,
ns_per_minimizer, sm_*`; `k` is the tool's `-l`, `L` its `-k`). Existing rows
are skipped unless `-f`. Cost: ~8 min per seed (all schemes) with `-j 2` on the
9 Gbp library. `ns_per_minimizer` is a single-thread timing: run it on an idle
machine.

`run_benchmark.sh` runs it after the builds (and with `-c`) for the run's
schemes, seeds, tables and library, logging to
`results/runs/<lib>/density_microbench.log`; `-D` skips it, `DENSITY_JOBS`
sets its `-j` (default 0 = all cores). It needs
`make -C kraken2/src benchmark_minimizers`.

```bash
make -C kraken2/src benchmark_minimizers
benchmark/kraken2/benchmark_minimizers.sh -o "$KRAKEN2_BENCH_ROOT" -E -S xor greedymini -n 5 -j 2
benchmark/kraken2/benchmark_minimizers.sh -E -S xor -n 5 -f                          # redo
benchmark/kraken2/benchmark_minimizers.sh -m -S xor mod-minimizer-xor -T 0xe37e28c4271b5a2d   # one toggle, seed_idx 0
```

Here `-N` is the random-DNA length, not the full-NCBI flag; the default library
is `db/standard`, `-m` / `-z` / `-E` pick the others (or set `SHARED_DB`).

Density sweeps over *k* (`density_sweep.sh` + `density_summary.py`) and the
measured GreedyMini numbers are in
[`kraken2/minimizers/README.md`](../../kraken2/minimizers/README.md#density);
raw rows in `tables/gm_center_density_*.csv`.

## Measurement protocol

Every build and classify follows the protocol in
[`benchmark/README.md`](../README.md): wait for a quiet node (classify also
waits for free huge-page memory), pin to one NUMA node with one thread per
physical core (build: `-t` cores + 2 spare; classify, one thread: 2 + 2 spare)
with memory bound to that node, and write a step record. The pinned CPU lists
(`cpus_build.txt`, `classification_performance_<ds>_mhg<m>_cpus.txt`) and the
step records are exported as `build_*` / `classify_*` telemetry columns of the
runs sheets.

## Notes

- `kraken2/data/simulator.pl` is not Mason 2; use the Zenodo reads.
- The FDA-ARGOS truth (`truth.tsv`) needs the taxonomy, so run
  `download_data.sh` after `prepare_strex_data.sh`.
