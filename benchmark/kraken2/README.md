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
`results/runs/<id>/`. The four ids:

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
| **cami** (supporting set, not in the original paper) | `data/cami/` | `-E` (or `-Q` with `-m`/`-N`) | genus/species sens / PPV / F1 over reads whose true taxon is in the library, plus the excluded fraction (see [CAMI II marine](#cami-ii-marine-datacami)) |

NCBI `-m`/`-z`/`-N` are a fork speed experiment on an unlabeled corpus
(`data/speed/` or legacy `data/benchmark_R*.fastq`). They do **not** score
Mason or ARGOS (classifying held-out genomes against full NCBI is the
optimistic case the paper avoided; ARGOS concordance in the paper also used
the exclusion DBs).

**Strain-exclusion** (Wood et al. 2019) is that *library*: 50 Table S3 taxa
held out, then Mason reads from those genomes. Scripts named `prepare_strex_*`
mean **str**ain-**ex**clusion. FASTA+taxonomy go to `db/strain_exclusion/`;
`run_benchmark.sh -E` writes indexes under `results/runs/strain_exclusion/`.

Accuracy (`evaluate_calls.py`): TP = true rank taxon or descendant; VP =
ancestor; FP = off-lineage; FN = unclassified. Strain-exclusion reports
sensitivity = TP/(TP+VP+FN+FP) and PPV = TP/(TP+FP). Table S1 Kraken 2 genus
targets: prokaryote F1 ≈ 0.966, viral F1 ≈ 0.754.

**Hit-group settings.** `k2 classify` defaults to `--minimum-hit-groups 2`
(Kraken 2 after 2.0.8-beta). Wood et al. / 2.0.8-beta `kraken2` did not apply
that threshold, so `--minimum-hit-groups 0` is the paper-equivalent setting
(Table S1/S5 scores). `run_benchmark.sh` always classifies every (seed,
dataset, scheme) cell exactly once with `--minimum-hit-groups 0` and once with
`2`; output files carry the suffix `_mhg0` / `_mhg2` (see "Run output layout").
For manual runs on the 208 index or Mason reads, add `--minimum-hit-groups 0`.

FDA-ARGOS is scored only as Table S5 genus concordance on SRA-assigned taxa
(percent of fragments: concordant, discordant, ancestor, unclassified).
Sensitivity and PPV are not reported for these real reads. Headline numbers
are the unweighted mean of the 25 SRR experiments; per-run rows and a pooled
fragment row are on the `accuracy_fda_argos` sheet (`scope` = `unweighted_mean` / `pooled` / `experiment`).

## Extra query set (CAMI)

The two Wood et al. sets resemble each other (simulated reads from NCBI
genomes; isolate SRA runs). CAMI II marine reads are a supporting set (not in
the original paper) that tests the minimizer schemes on a different kind of
reads. It is optional: `run_benchmark.sh` classifies it only when prepared
(`-C` skips it; `-Q` adds it to a `-m`/`-z`/`-N` run; `-E` includes it
automatically). The existing sets and flags are unchanged. A ZymoBIOMICS D6300
mock community was tried and dropped on 2026-10-07 (little information on
scheme differences); its harness, scorer and results are gone.

| Dataset | Reads | Size on disk | Truth |
|---|---|---|---|
| `data/cami/` | CAMI II marine short-read sample 0, 2x150 bp, 250 000 of 16 647 395 pairs | ~170 MB | per read (CAMI gold standard) |

Subsampling is systematic (every k-th pair, k = pairs in the run // 250 000,
no RNG), so every machine gets the same reads.

### CAMI II marine (`data/cami/`)

```bash
./benchmark/kraken2/prepare_cami_data.sh -o "${KRAKEN2_BENCH_ROOT}"   # ~6 GB download, deleted after sampling
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

  Recommended: `-E` (same library as the other
  query sets, 90 % of reads keep a reference at genus level). Full `standard`
  is the best-covered but needs the 741 GB DB (`-N -Q`). The subsampled
  `standard_p10` / `standard_p0.5` drop most genera: use them only for timing.
- Caveats: simulated reads (ART error model), 2017 taxonomy vs the library's
  taxonomy (taxids missing from `nodes.dmp` count as `no`), strain-level
  novelty is not scored. Only the marine sample is validated; `-s strain|plant`
  and `-n SAMPLE` select other CAMI II short-read samples
  (`strmgCAMI2` / `rhimgCAMI2`; plant has host reads).

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
# set KRAKEN2_BENCH_ROOT, CONDA_ENV
source benchmark/kraken2/.env

make -C src test
make -C src all ENABLE_MINIMIZER_PROF=1

# Strain-exclusion Mason + Table S5 ARGOS reads
./benchmark/kraken2/download_data.sh -o "${KRAKEN2_BENCH_ROOT}"

# Optional supporting query set (CAMI II marine)
./benchmark/kraken2/prepare_cami_data.sh -o "${KRAKEN2_BENCH_ROOT}"

# NCBI libraries + standard_p10 / standard_p0.5
./benchmark/kraken2/download_library.sh -o "${KRAKEN2_BENCH_ROOT}"

# Optional: strain-exclusion FASTA → db/strain_exclusion
./benchmark/kraken2/prepare_strex_data.sh -o "${KRAKEN2_BENCH_ROOT}"

# 10% NCBI library, unlabeled speed FASTQs only (not Tables S1/S5)
./benchmark/kraken2/run_benchmark.sh -o "${KRAKEN2_BENCH_ROOT}" -m -S xor greedymini -n1

# Paper: held-out library × Mason + FDA-ARGOS
./benchmark/kraken2/run_benchmark.sh -o "${KRAKEN2_BENCH_ROOT}" -E -S xor greedymini -n1
```

`-a` skips FDA-ARGOS, `-C` skips CAMI; `-Q` adds
CAMI to a `-m`/`-z`/`-N` run. `-c` re-summarizes an existing run (auto-picks the only
tree with a manifest, or pass `-m`/`-z`/`-N`/`-E`).

`./benchmark/kraken2/run_benchmark.sh -m -n8` → `results/runs/standard_p10`.  
`./benchmark/kraken2/run_benchmark.sh -z -S xor -n1` → `results/runs/standard_p0.5`.

`db/standard_p10` and `db/standard_p0.5` are created by `download_library.sh` unless `-S`.
Rebuild only smoke from an existing `db/standard`:  
`./benchmark/kraken2/download_library.sh -Z -o "${KRAKEN2_BENCH_ROOT}"`.

## Directory layout (`KRAKEN2_BENCH_ROOT`)

```
data/
  strain-exclusion/   # Mason PE + truth.tsv  (download_data.sh)
    R1.fastq  R2.fastq  truth.tsv
  fda-argos/          # Table S5 SRRs         (download_data.sh)
    R1.fastq  R2.fastq  [truth.tsv]
  cami/               # CAMI II marine        (prepare_cami_data.sh)
    R1.fastq  R2.fastq  truth.tsv  truth_extra.tsv  taxonomic_profile.txt  SOURCE.txt
db/
  standard/           # full NCBI library
  strain_exclusion/   # Wood 2019 held-out FASTA (prepare_strex_data.sh)
  standard_p10/       # ~10%
  standard_p0.5/      # ~0.5%
zenodo/               # Wood 2019 reference tarball (download_data.sh, prepare_strex_data.sh)
results/
  runs/
    {standard,strain_exclusion,standard_p10,standard_p0.5}/   # one dir per library
      {order dir}/{seed}/{db,results}   # order dir: xor | greedymini/k15_at0-cg1 | ...
      benchmark_comparison.txt          # run summary text
      minimizer_density.csv             # INPUT: density microbenchmark rows
  csvs/kraken2_<sheet>.csv              # the only copy of the exported tables (all libraries)
  kraken2_results.xlsx                  # the only workbook (all libraries)
```

- **Results workbook** (same layout and column naming as the sshash / kmc apps:
  `export_layout.py`, writer in `../common/xlsx_book.py`; no charts). Each
  `compare_benchmark.py` call (the one at the end of `run_benchmark.sh`, or `-c`)
  builds this library's tables in memory, replaces its rows in the tables already
  in `results/csvs/` and rewrites `results/kraken2_results.xlsx`; there is no
  per-library copy (`results/runs/<id>/` keeps only `benchmark_comparison.txt`).
  `python3 benchmark/kraken2/export_layout.py --base-dir <ROOT>` re-renders the
  workbook from `results/csvs/`. One CSV per sheet, same rows and columns as the sheet. Column names:
  the raw value carries its unit (`build_wall_s`, `classify_max_rss_gb`,
  `classified_pct`, `lookups_per_mbp`); the aggregated sheets add
  `<metric>_mean`, `<metric>_std` and, for metrics with a direction,
  `<metric>_ratio` = mean / the `xor` mean of the same library, dataset and
  `minimum_hit_groups`, coloured green = better, red = worse. Every cell is
  classified once per `--minimum-hit-groups` value (0 = paper setting, 2 = k2
  default), so the classify-side sheets (`summary`, `runs`, `profile`,
  `vs_reference`, the four `accuracy_*`) carry `minimum_hit_groups` as an
  app-parameter column (after `mod_t`, before `n_repeats` / `seed_idx`; column order: dataset, order, app parameter, seed, metrics) / slice key and are sorted dataset → minimum_hit_groups →
  order; `build` / `build_runs` do not (building is independent of it). The
  replicates are the seeds (one classify run each; no repeats, no cold round),
  so `_std` is across seeds; `n_repeats` = number of runs (repeats) aggregated into the row's mean / std. `runs` and `build_runs` also carry I/O telemetry (`build_major_pagefaults`, `build_read_gb`, `classify_major_pagefaults`, `classify_read_gb`; runs-only, not in `build` / `summary`) and, in a `telemetry` band at the end, the text columns `build_pinned_cpus` (`cpus_build.txt`) and `classify_pinned_cpus` (`classification_performance_<dataset>_mhg<n>_cpus.txt`): every timed step is pinned to one NUMA node, one thread per physical core plus 2 spare cores, with its memory bound to that node (`benchmark/common/pin_cpus.sh`, `bench_step.py`); blank when the file is absent. The `overview` legend only lists sheets and columns that exist in the workbook. Sheets, in order:

  | Sheet | Rows |
  |---|---|
  | `overview` | column reference: what every sheet and column means, by category |
  | `libraries` | one per Kraken library: sequences, taxa, FASTA GB, role, query sets and orders run |
  | `datasets` | one per query set: reads, read length, `truth` (simulated / real), `scoring`, provenance, libraries |
  | `orders` | one per (order variant, `k`, `w`, `L`), columns `order, k, w, L, gm_k, projection, gm_table, mod_t, description, cli`, then `sketching, rank_order, minimizer_spaces, spaced_seed_mask`, then the density factors; unified names: `k` = minimizer length (the kraken2 option `l`, 31), `w` = window (5), `L` = k + w - 1 = kraken2's own k-mer length (35); no empty cells (the export fails if an order lacks density rows): CLI, density factors (`expected_density_factor_mean/_std`, and per library `<library>_particular_density_factor_mean/_std`; also on `build`), `gm_table` (`w5_k15_canonical` / `w5_k6_canonical`), spaced seed. Orders are shown by short name (`greedymini`, `mod-minimizer-greedymini`; the table is in `gm_table`) |
  | `build` | library × order, mean / std / ratio over seeds, columns in the SSHash order (id / params, density, minimizers, space, build): `expected_density_factor`, `density_factor`, `distinct_minimizers` (the `size` field of the `hash.k2d` header), then `build_wall_s`, `build_max_rss_gb`, `index_size_gb` and the density microbenchmark (`expected_density_factor`, `density_factor`). Build does not depend on the query set, so it appears once here, not on every dataset row |
  | `build_runs` | the per-seed rows behind `build` (same split as `summary` / `runs` in the sshash / kmc apps) |
  | `summary` | library × dataset × `minimum_hit_groups` × order, classify phase only; right after the identity columns (`library, dataset, order, gm_k, projection, gm_table, mod_t, minimum_hit_groups, n_repeats`) the order's density factors and `distinct_minimizers` (repeated for every dataset / `minimum_hit_groups` block), then the classify columns: `classify_wall_s` (= `classify_load_s` + `classify_run_s`; processing = the classify `processed in` seconds, DB load = wall − processing), `classify_max_rss_gb`, `classify_mreads_per_min` (per processing time), `classified_pct` and the accuracy of the dataset (std rounded to ≥ 4 significant digits) |
  | `runs` | one per seed × library × dataset × order, the rows behind `summary`: the raw numbers of the single classify run of that seed |
  | `profile` | classify phase profile (`scan_s`, `lookup_s`, `resolve_s`, `scan_pct_of_phases`, `lookups_per_mbp`, `ns_per_lookup`, `lookup_miss_pct`, `uhs_pct`, ...; `experiments/profile_run.sh`, a classify built with `ENABLE_MINIMIZER_PROF=1`), library × dataset × `minimum_hit_groups` × order, mean / std / ratio; only present when `classify_profile_<dataset>_mhg<m>.txt` files exist |
  | `vs_reference` | paired agreement with `xor` of the same seed, dataset and `minimum_hit_groups` (`genus_mape_pct`, `concordance_pct`, `jaccard_pct`, ...), mean / std |
  | `accuracy_strain_exclusion`, `accuracy_fda_argos`, `accuracy_cami` | accuracy of that dataset per slice (`minimum_hit_groups` × `group` × `rank`, or SRA run for FDA-ARGOS): one row per order, `<metric>_mean`, `<metric>_std` over seeds and `<metric>_ratio` to `xor` of the same slice; sorted dataset → `minimum_hit_groups` → slice → order, slice cells merged, border between slices. `accuracy_strain_exclusion` and `accuracy_cami` have the same columns in the same order (identity, counts `n_truth n_calls n_scored tp vp fp fn no`, `excluded excluded_fraction`, then CAMI-only `ex_unclassified ex_consistent ex_off`, then `sensitivity ppv f1`; see the shared definitions below). Groups: `all`, `prokaryote`, `virus` (both sheets), then CAMI `bacteria`, `archaea`, `other`, `novelty:*`; CAMI `prokaryote` = the archaea + bacteria reads scored as one group by `score_cami.py` (as `prokaryote` is in the strain-exclusion truth). Higher-is-better ratios: sensitivity, ppv, f1, concordant; lower: discordant; neutral (no ratio): ancestor, unclassified, excluded, counts. Per-seed rows are not on these sheets (the raw numbers stay in each run's `results/*_metrics_mhg<m>.tsv`). Omitted when the dataset did not run |
  | `paper_check` | paper equivalence for `xor` at `minimum_hit_groups` 0 on the `strain_exclusion` library: our strain-exclusion genus sensitivity / PPV / F1 (prokaryote, virus; Table S1) and FDA-ARGOS genus concordant / discordant / ancestor / unclassified (unweighted mean over experiments; Table S5) next to the published Kraken 2 values (`tables/paper_kraken2_published.csv`, derived from `table_s1.csv` and `table_s5_kraken2_expected.csv`), `difference` and `within_tolerance` (absolute tolerance 0.01, `PAPER_TOLERANCE` in `export_layout.py`). Also `kraken2_paper_check.csv`. Differences are reported as they are |
  | `density_effect`, `density_effect_points` | how closely each metric of `build_runs` / `runs` follows the density: Pearson r of the metric and the density factor, both as ratios to `xor`'s run of the same seed and setting (library, or dataset and `minimum_hit_groups`), over all settings and per setting, and the paired points (`benchmark/common/density_effect.py`; also `csvs/density_effect_kraken2.csv`) |
  | `significance` | whether each scheme changes each metric of the paper's figures: paired by seed with `xor` of the same setting (library for the build; dataset, `minimum_hit_groups` and library for classify), t-test on the log ratios, Holm, TOST at +-5% (`benchmark/common/significance.py`; also `csvs/significance_kraken2.csv`) |

  Under the data of `build` and `summary` two rows, as in the SSHash summary
  sheet: `correlation with density` and `correlation with distinct minimizers`.
  Every `*_ratio` column holds a live `IFERROR(CORREL(<reference ratio column>,
  <this column>),"")` over all data rows (reference: `density_factor_ratio`
  resp. `distinct_minimizers_ratio`), with the computed value cached next to the
  formula. `build` has 5 rows (n = 5, in the label); in `summary` the density
  and minimizer values repeat per dataset / `minimum_hit_groups` block, so the
  number of independent points is the number of orders. The CSVs hold data only;
  the values are in `kraken2_correlations.csv` (`sheet, row_label, column,
  value, n`).

  **Accuracy definitions shared by `accuracy_strain_exclusion` and `accuracy_cami`**
  (Wood et al. 2019, both scorers identical): TP = the call is the true taxon at
  the rank or a descendant, VP = the call is an ancestor of it, FP = the call is
  off-lineage, FN = unclassified; `sensitivity = TP/(TP+VP+FN+FP)` (= `tp / n_scored`),
  `ppv = TP/(TP+FP)`, `f1 = 2·sens·ppv/(sens+ppv)` (0 when undefined);
  `n_scored = tp+vp+fp+fn`. **Denominator:** sensitivity is computed over
  `n_scored` reads. Strain-exclusion: `n_scored = n_truth` (every true genus and
  species is in the library, `excluded` = 0, measured with `score_cami.py`
  coverage mode); CAMI: reads whose true genus / species has no sequence in the
  library are excluded (reported, not scored: 9.6 % genus / 22.6 % species of the
  reads that have a node at the rank), so CAMI sensitivities are over a
  different, library-covered subset and are not directly comparable to the
  paper's definition over all reads. `ex_unclassified` / `ex_consistent` /
  `ex_off` (CAMI only) say what the classifier did with the excluded reads (no
  call / call on the true lineage / off-lineage call).

  Accuracy columns have no dataset prefix because each `summary` / `runs` row is
  one dataset. The kinds differ and are not comparable across datasets (never
  averaged): strain-exclusion and CAMI share the read-level names
  `genus_` / `species_` `sensitivity` / `ppv` / `f1` (CAMI over reads whose taxon
  is in the library, plus `*_excluded_fraction`); FDA-ARGOS has
  `genus_{concordant,discordant,ancestor,unclassified}_fraction`. A row fills only its
  own dataset's column group; columns no row uses are dropped. Neutral columns
  without ratio: `classified_pct`, `genus_ancestor_fraction`,
  `genus_unclassified_fraction`, `*_excluded_fraction`, the profile counts and
  path split. If `results/csvs/kraken2_runs.csv` lacks the `minimum_hit_groups` column
  (not produced by the current `compare_benchmark.py`) it is ignored with a
  notice; `results/runs/<id>/minimizer_density.csv` (density microbenchmark) is
  an input and is kept. `results/runs/.libraries_cache.json` caches library
  FASTA statistics. `kraken2_*` CSVs of sheets
  that no longer exist (e.g. `profile` without timing files) are deleted when
  the new ones are written.
  Only libraries that have runs are listed (`libraries`, `datasets`, every sheet):
  `db/<id>` directories without runs (e.g. the standard libraries) never appear.
- **One unified benchmark dir, several invocations.** Run
  `run_benchmark.sh -E -S <schemes> ...` several times into the same `-o` dir
  (e.g. `-S xor greedymini`, then `-S xor-pcg-hash mod-minimizer-xor
  mod-minimizer-greedymini`): an invocation replaces only the
  `.benchmark_run_manifest` lines of ITS run dirs (order + params; other lines are
  kept), so the final
  `compare_benchmark.py` / merged workbook cover every order. A scheme/seed whose
  `results/runs/<lib>/<order dir>/<seed>/db/hash.k2d` exists is not rebuilt, only classified
  (`db/` may be a directory of symlinks to a database built elsewhere, with
  `results/build_performance.txt` copied next to it).

Existing `data/benchmark_R1.fastq` is still used as the NCBI unlabeled speed
set. Prefer `data/speed/R1.fastq`. Older trees (`sample_db/`, top-level `p10/`,
`db/` with `library/` directly inside, `db_strain_exclusion/`) are still
used if the new path is missing.

## Schemes (`-S`)

| Order (`-S`) | Run directory (`results/runs/<lib>/<dir>/<seed>/`) | CLI |
|-------|------|-----|
| `xor` | `xor` | `--sketching minimizer --order xor` |
| `xor-pcg-hash` | `xor-pcg-hash` | `--sketching minimizer --order xor-pcg-hash` |
| `greedymini` | `greedymini/k<gm_k>_<projection>` (e.g. `greedymini/k15_at0-cg1`) | `--sketching minimizer --order greedymini --greedymini-file w5_k15_*.gm --gm-projection <projection>` |
| `mod-minimizer-xor` | `mod-minimizer-xor/t<mod_t>` (`t6`) | `--sketching mod-minimizer --order xor` |
| `mod-minimizer-greedymini` | `mod-minimizer-greedymini/k<gm_k>_<projection>_t<mod_t>` (`k6_at0-cg1_t6`) | `--sketching mod-minimizer --order greedymini --greedymini-file w5_k6_*.gm --gm-projection <projection>` |

The five names are the canonical orders shared by all four apps
(`benchmark/common/order_names.py`). An order never carries parameters: the
GreedyMini `gm_k`, the `projection` (`at0-cg1` = A,T -> 0, C,G -> 1, the
`at1-cg0` = flipped; `-P at0-cg1,at1-cg0` / env `GM_PROJECTION` (comma list,
default both) / `--gm-projection` of the tools),
the table stem `gm_table` and the mod-minimizer `mod_t` (= r + ((l - r) mod w),
r = 4; 6 for l = 31, k = 35) are separate fields: manifest columns, columns of
`minimizer_density.csv`, and the `gm_k, projection, gm_table, mod_t` columns of
every sheet (part of every block / merge / group key). Run directories are
`order_names.order_dir_name`. One invocation runs every GreedyMini scheme
(`greedymini`, `mod-minimizer-greedymini`) for every projection in the list
(default both, so the default run is 7 variants: xor, xor-pcg-hash,
mod-minimizer-xor once; each GreedyMini scheme x 2 projections); each
(scheme, projection) is its own run directory, manifest row and density row,
and a variant whose `db/hash.k2d` exists is not rebuilt (density rows already in
`minimizer_density.csv` are skipped likewise). `run_benchmark.sh --dry-run`
lists the variants and which are already done without building anything;
`-P at0-cg1` restricts to one projection. The report labels read
`greedymini (k15_at0-cg1)` / `greedymini (k15_at1-cg0)`.

`mod-minimizer` is strand-symmetric mod-sampling (this fork, not a published
scheme; `kraken2/src/mod_minimizer.h`): mod-sampling (Groot Koerkamp & Pibiri,
2024) with t == l (mod w) and the centre tie rule shared with minimap2 and
sshash (`common/minimizer-orders/mod_minimizer_rule.h`): a tie at the minimum
takes the offset closest to the k-mer's centre, two equally close ones by the
k-mer <= its reverse complement, so a read and its reverse complement give the
same minimizers. Databases built with the earlier forward-only mod-minimizer,
or with the G/T-strand tie rule, are rejected (rebuild them).

Names are exact (no aliases: every other spelling, old or new, is rejected). The
name is `sketching[-]order` (`scheme_label` in `kraken2/src/sketching_schemes.h`);
the same names are the `scheme` column of the manifest, all metrics TSVs and
`minimizer_density.csv`, and the `order` column of the sheets. The C++ tools
`scan_dump` / `benchmark_minimizers` take the same names with `--gm-projection`.
(sshash uses `--greedy-mini-file`, kmc its own options.)

greedymini tables (`../../orders/GreedyMini/`): `-g` / env `GM_FILE`
(`w5_k15_canonical.gm`, used by `greedymini`), `-M` / env `MOD_GM_FILE`
(`w5_k6_canonical.gm`, used by `mod-minimizer-greedymini`).

## Run output layout (`results/runs/<lib>/<order dir>/<seed>/results/`)

One classify per (seed, dataset, scheme, `--minimum-hit-groups` m in {0, 2});
`<ds>` = `strain-exclusion` | `fda-argos` | `cami` | `speed`.

| File | Content |
|------|---------|
| `build_performance.txt` | `time -v` of `k2 build` (stderr of the build, then GNU time block) |
| `<ds>_classification_mhg<m>.txt` | `k2 classify --output` |
| `<ds>_report_mhg<m>.txt` | `k2 classify --report` |
| `classification_performance_<ds>_mhg<m>.txt` | classify stderr (`Loading database information... done.`, `N sequences (X Mbp) processed in Ys (...)`, `N sequences classified (P%)`, `N sequences unclassified (P%)`) followed by the GNU `time -v` block (wall clock, max RSS); DB-load time = wall - "processed in" |
| `classify_profile_<ds>_mhg<m>.txt` | only for `ENABLE_MINIMIZER_PROF` classify builds (via `KRAKEN2_MINIMIZER_PROF`; `experiments/profile_run.sh`); absent otherwise |
| `strain_exclusion_metrics_mhg<m>.tsv`, `fda_argos_metrics_mhg<m>.tsv`, `cami_metrics_mhg<m>.tsv` | scorer output; second column `minimum_hit_groups` (= m) |

`results/runs/<lib>/.benchmark_run_manifest` has one tab-separated line per scheme and
seed, written after that seed finished classifying: `scheme` (the order), `context`
(`<lib>/<order dir>/<seed>`), `results dir`, `toggle` (bit string), `seed index`,
`gm_k`, `projection`, `gm_table`, `mod_t` (blank when not applicable).
There is no repeat/cold-round mechanism: replicates are the seeds only.

## Script map

| File | Role |
|---|---|
| `setup.sh` | apt + conda env |
| `.env` | `KRAKEN2_BENCH_ROOT` (<ROOT>; optional `KRAKEN2_DATA_DIR`, `KRAKEN2_RESULTS_DIR`), threads |
| `benchmark_setup.sh` | conda + `scripts/k2` on PATH |
| `benchmark_interface.sh` | `-S` orders, `-g/-M/-P` GreedyMini tables + projection, `-m/-z/-N` NCBI size, `-E` strain-exclusion library |
| `export_layout.py` | `results/kraken2_results.xlsx` + `results/csvs/kraken2_<sheet>.csv` (sheets above); `../common/xlsx_book.py` / `xlsx_charts.py` / `xlsx_util.py` write and read .xlsx without openpyxl, `../common/xlsx_check.py` validates the package |
| `create_sample_db.py` | NCBI `library.fna` subsample → `db/standard_p10` / `db/standard_p0.5` |
| `download_library.sh` | NCBI taxonomy + standard libraries + subsamples |
| `download_data.sh` | `data/strain-exclusion` + Table S5 `data/fda-argos` |
| `prepare_strex_data.sh` | Zenodo FASTA → `db/strain_exclusion/` |
| `prepare_cami_data.sh` / `prepare_cami_reads.py` | CAMI II marine sample → `data/cami/` (subsample + truth) |
| `score_cami.py` | CAMI library-aware sens/PPV/F1 |
| `tests/` | `python3 -m unittest discover -s benchmark/kraken2/tests` (synthetic data, no DB) |
| `prepare_fda_argos_reads.py` | ENA download of Table S5 `sra_run` accessions |
| `run_benchmark.sh` | build + classify + score + CSV/XLSX |
| `bench_const.py` | shared `DEFAULT_SEEDS`, canonical scheme order, `scheme-id` (run dir + gm_k/projection/gm_table/mod_t), density metrics |
| `compare_benchmark.py` | statistics text report + per-run numbers fed to `export_layout.py` |
| `score_kraken.py` / `score_argos.py` | strain-exclusion sens/PPV/F1; ARGOS Table S5 genus (a–d) |
| `benchmark_minimizers.sh` | density microbench (CSV) |
| `density_sweep.sh` / `density_summary.py` | density factor vs ℓ over several random streams (calls `benchmark_minimizers`), mean ± std table |

**Units.** All memory / size columns are `_gb` in decimal GB (1 GB = 10^9 bytes, `../common/units.py`): peak RSS from GNU time (KiB) is converted via bytes (`gb_from_kib`), file sizes with `gb_from_bytes` (`index_size_gb`, `build_max_rss_gb`, `classify_max_rss_gb`, `fasta_gb`, Table S1 `max_rss_gb`). Earlier exports divided KiB by 1024^2 (GiB); re-export to refresh.

## Density microbench (`benchmark_minimizers.sh`)

Measures the minimizer density of every scheme at the benchmark's own
configuration (ℓ = 31, k = 35 so w = 5, minimizer spaces s = 7, greedymini table
`-g`, same as the built databases) and the `NextMinimizer` cost. Two densities:

- **expected** (`dens_expected`, `fac_expected`) — random i.i.d. ACGT of
  `-N` bases (default 100 Mbp in the script; the tool is Monte Carlo, so it varies
  slightly with the RNG seed `--seed`);
- **particular** (`dens_particular`, `fac_particular`) — the library FASTAs of
  the shared DB (`db/<lib>/library/**/*.f{na,a,asta}`), read in full (8.7 GB, ~9.3
  Gbp for `strain_exclusion`), no sampling;
- `fac_* = dens_* × (w+1)` (random minimizer ⇒ 2), `dens_theory` is
  `2/(w+1)` (Minimizer) / `1/w` (Mod); `ns_per_minimizer`, super-mer lengths
  `sm_mean_*` / `sm_max_*`.

**Per seed.** Seed *i* uses the toggle mask `DEFAULT_SEEDS[i-1]`
(`bench_const.py`, bit string converted to `0x…`, the value `k2 build
--toggle-mask` gets), passed as `-T` with `--seed-idx i`. One row per scheme and
seed goes into `results/runs/<lib>/minimizer_density.csv` (columns `scheme,
gm_k, projection, gm_table, mod_t, seed_idx, toggle_hex, k, w, L, spaces, dens_*, fac_*, ns_per_minimizer, sm_*`).
Names are the unified ones: `k` = minimizer length (the tool's `-l`), `w`, `L` = k + w - 1 (the tool's `-k`); a file with the tool's old `l, k, w` header is mapped on read and rewritten. Re-running replaces the rows of the same scheme (+ gm fields) + seed and keeps the rest; a
scheme/seed whose row (same toggle, ℓ, k, s) exists is skipped unless `-f`.
Cost: ~8 min per (seed, all schemes) with `-j 2` on the 9 Gbp library, dominated
by the particular pass. `ns_per_minimizer` is a single-thread micro-timing: run
on an idle machine for comparable numbers.

**In the pipeline.** `run_benchmark.sh` runs it after the builds of all seeds
(and also with `-c`), for the same `-S` schemes, `-n` seeds, GM tables and library
as the run, writing the log to `results/runs/<lib>/density_microbench.log`; `-D` skips it
(`DENSITY_JOBS` sets `-j`, default 0 = all cores). `compare_benchmark.py` then
reads the CSV (`load_density_rows`) into the `build` sheet. The binary must exist:
`make -C kraken2/src benchmark_minimizers` (otherwise a warning, no density).

```bash
make -C kraken2/src benchmark_minimizers
# all seeds of a run, e.g. the strain-exclusion library (-m / -z / default: other libraries)
./benchmark/kraken2/benchmark_minimizers.sh -o "$KRAKEN2_BENCH_ROOT" -E \
    -S xor greedymini -n 5 -j 2
./benchmark/kraken2/benchmark_minimizers.sh -E -S xor -n 5 -f   # redo
./benchmark/kraken2/benchmark_minimizers.sh -m -S xor mod-minimizer-xor -T 0xe37e28c4271b5a2d  # one toggle, seed_idx 0
```

(`-N` here is the random-DNA length, not the full-NCBI library flag; the default
library is `db/standard`, `-m` / `-z` / `-E` pick the others, or set `SHARED_DB`.)

Density **sweeps** over ℓ and streams (mean ± std per ℓ, `xor` baselines) use
the same tool through `density_sweep.sh` + `density_summary.py`; the exact
commands and the measured greedymini tables are in
[`../../kraken2/minimizers/README.md`](../../kraken2/minimizers/README.md#density).
Raw sweep rows: `tables/gm_centre_density_*.csv`.

## Measurement protocol

Every build and classify is timed under the protocol in [`../README.md`](../README.md): node-quiet wait (classify
also waits for free huge-page memory), pinning (build: `--threads` cores + 2 spare; classify, which is
`--threads 1`: 2 + 2 spare = 4 cores; one NUMA node, one thread per core; lists in `cpus_build.txt` and
`classification_performance_<dataset>_mhg<n>_cpus.txt`, exported as `build_pinned_cpus` /
`classify_pinned_cpus`; with its memory bound to the slot's node and a step record (`step_*.json`, exported as `build|classify_node`, `_huge_fraction`, `_foreign_busy_cores`, ... -- see the protocol)). The headline classify time is
`classify_run_s` (without the database load) and per-step wall, CPU (`build_cpu_s`,
`classify_cpu_s`), RSS, major page faults and reads (`*_major_pagefaults`, `*_read_gb`, runs sheets only).
Classify wall time splits into `classify_load_s` (database load) and `classify_run_s` (the reads); compare
orders on `classify_run_s`. kraken2 has no native run: the reference order is `xor`. `experiments/full_run.sh`
is the full run (5 seeds, 35 databases of about 9.4 GB each).

## Notes

- `data/simulator.pl` is not Mason 2. Prefer Zenodo FASTQs.
- FDA-ARGOS runs are the 25 `sra_run` values in `benchmark/kraken2/tables/table_s5.csv`.
- Table S5 paper averages are unweighted per-experiment; `score_argos.py`
  writes that mean plus a pooled-fragment row. Per-SRR percentages need
  `truth.tsv` (`download_data.sh` after taxonomy is present).
