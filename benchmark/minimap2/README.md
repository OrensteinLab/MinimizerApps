# Minimap2 benchmark

Compares Minimap2's minimizer schemes (native, xor, xor-pcg-hash, greedymini,
mod-minimizer-xor and mod-minimizer-greedymini) the way Jain et al. 2020
("Weighted minimizer sampling improves long read mapping" -- Winnowmap,
Bioinformatics 36(Supp_1):i111-i118) compared Winnowmap's weighted minimizers
with Minimap2's own: density, index build and mapping time and memory, and
mapping accuracy, overall and on reads from inside a long near-identical repeat.
The GreedyMini table (`gm_table`, whose k is `gm_k`) and the polarity of Φ
(`--gm-projection at0-cg1|at1-cg0`) are parameters of the greedymini schemes,
not part of the scheme name; each polarity is its own set of runs and rows (see
[`orders/README.md`](../../orders/README.md)).

Notation: k = selected k-mer length (`-k`), w = number of k-mers per window
(`-w`), L = k + w - 1, t = the mod-minimizer's t-mer length.

## What is measured

Per (dataset, scheme) run, from `minimap2 -d` (index build) and one mapping pass:

- **Density**: `density_factor = total_occurrences * (w+1) / n_sketch_kmers`.
  `total_occurrences` = distinct minimizers x average occurrences (both from
  `-d`'s stderr); `n_sketch_kmers` counts the k-mers Minimap2 sketches
  (`lib/sketch_kmers.py`: length - k + 1 per N-free segment, homopolymer-compressed
  under `map-pb`, raw under `map-ont`). Total occurrences, not distinct
  minimizers, are the numerator: distinct minimizers collapse a repeat's copies
  and undercount selected positions (`distinct_minimizers`/`avg_occurrences` are
  reported, informational only). `-d`'s "average spacing" divides by the raw
  length and is kept only as `density_factor_raw_bp`. Same definition as
  `benchmark/common/density.py`.
- **Index build and mapping time/peak memory** (`/usr/bin/time -v`, as in the
  other benchmarks).
- **Mapping accuracy**: fraction of simulated reads whose primary alignment
  overlaps the true origin by at least 10% of the union of the two intervals
  (`lib/score_accuracy.py`'s `frac_correct`, the criterion of the Winnowmap and
  Minimap2 papers), overall and for reads inside the repeat
  (`frac_correct_in_repeat`), plus mean identity (PAF `dv:f:`/`de:f:`) and a
  stricter 50%-of-true-span variant (`frac_correct_strict50`). Real reads (D3)
  have no truth: mapped fraction and identity only, as in the paper.

## Datasets

**`smoke`** (default, `$DATASET` unset): a synthetic ~300kb genome with a
20kb region duplicated verbatim at a distant location (`lib/make_dataset.py`,
no download dependency, fast to regenerate) plus simulated long reads with
known ground truth (`truth.tsv`), 30% deliberately sampled from inside the
repeat. Good for fast iteration and CI-scale sanity checks; too small and
too artificial for real accuracy conclusions.

The three benchmark datasets are the Winnowmap paper's (Jain et al. 2020,
Table 1), built and evaluated the way the paper describes:

| Dataset | Genome | Reads | Preset (k) | Evaluation |
|---|---|---|---|---|
| `chm13_chrX` (D1) | T2T CHM13 chrX v0.7, 154.3 Mbp | PBSIM 1.0.3 CLR, 10% error, 15 kb mean (100-25,000 bp), 5x: 51,459 reads, N50 15,170 (paper: 51,424, N50 15,184) | `map-pb` (19) | unmapped / incorrectly mapped (truth); centromere = DXZ1 array chrX:57,828,561-60,934,693 (the paper's Table 3), 1,090 reads |
| `grch38_wg` (D2) | GRCh38 with ALT contigs (UCSC hg38.fa: 455 sequences, 3.21 Gbp, of which 261 ALT contigs / 109.5 Mb) | same simulation, 0.5x, reads containing "NN" dropped (paftools `pbsim2fq`'s rule): ~101.4k reads (paper: 101,571) | `map-pb` (19) | as D1; plus a chr9 segmental-duplication repeat set (our addition, not in the paper) |
| `chm13_rel3_ont` (D3) | same GRCh38 with ALT contigs (symlink to D2's) | 1,000 real CHM13 rel3 ultra-long ONT reads, uniform random sample of reads >= 1 kbp with mean quality >= 7 (the ONT pass threshold; our addition: the paper states no quality filter, and 31% of rel3 reads >= 1 kbp are below Q7, almost all unmappable) | `map-ont` (15) | no truth: unmapped %, time, memory, density only (as in the paper) |

k is each preset's own default, as the paper ran Minimap2 with `-cx map-pb` /
`-cx map-ont` (`lib/constants.sh:preset_default_k`; `K=` overrides), and w = 10:
k = 19, L = 28 on D1/D2 and k = 15, L = 24 on D3. Every scheme uses that k
(greedymini's `w10_k15_canonical` table is center-extended to k = 19), and the
mod-minimizer's t = 4 + (k-4) mod w (9 at k = 19, 5 at k = 15).
mod-minimizer-greedymini ranks the t-mers by the shared
`orders/GreedyMini/w20_k9_canonical` table (`MOD_GM_TABLE`; its k = t = 9, its
w = w + k - t = 20), so it runs only on D1 and D2; at k = 15 (t = 5) `run`
skips it (`scheme_table_fits`). Each dataset's preset is a one-word
`data/<dataset>/preset` file (absent = `map-pb`; `$PRESET` overrides it),
recorded in every run's `meta.json`.

To regenerate, from `benchmark/minimap2/` (`$E` = the `minimap2_bench` conda
env; `setup.sh` builds PBSIM 1.0.3 into `$E/opt/pbsim1`; the downloads need `curl`, `md5sum` and, for D3, `pigz`):

```bash
P1=$E/opt/pbsim1
# D1
mkdir -p "$MINIMAP2_DATA_DIR/chm13_chrX" && cd "$MINIMAP2_DATA_DIR/chm13_chrX"
curl -sL -o genome.fa.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.chrX_v0.7.fasta.gz
gunzip genome.fa.gz  # verify md5 89b3dd61db66177dd830527b920956fa first
printf 'chrX_v0.7\t57828561\t60934693\tcentromere_dxz1_array\n' > repeats.bed  # the paper's Table 3 coordinates
cd - && $E/bin/python3 lib/pbsim_reads.py --simulator pbsim1 --dataset-dir "$MINIMAP2_DATA_DIR/chm13_chrX" \
    --pbsim-bin $P1/bin/pbsim --model $P1/data/model_qc_clr \
    --depth 5 --length-mean 15000 --accuracy-mean 0.90 --seed 1

# D2: UCSC hg38.fa.gz (GRCh38 with ALT contigs; md5 from the same directory's md5sum.txt) as
#   genome.fa (+ the chr9 segdup pair as repeats.bed), then the same command with
#   --dataset-dir .../grch38_wg --depth 0.5
#   curl -sL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# D3: random 1,000 reads >= 1 kbp, mean Q >= 7, of the 355 GB rel3 FASTQ (md5 92026d97a898c2f5b65074048a1caabf),
#   streamed once (a parallel ranged download into parts is much faster than one connection)
mkdir -p "$MINIMAP2_DATA_DIR/chm13_rel3_ont" && cd "$MINIMAP2_DATA_DIR/chm13_rel3_ont"
ln -s ../grch38_wg/genome.fa ../grch38_wg/genome.fa.fai . && echo map-ont > preset
curl -sL https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/nanopore/rel3/rel3.fastq.gz \
  | tee >(md5sum > rel3.fastq.gz.md5) | pigz -dc \
  | $E/bin/python3 "$OLDPWD/lib/sample_reads.py" --n 1000 --min-len 1000 --min-mean-q 7 --seed 1 --out-dir .
```

Any other dataset: provide `data/<dataset>/{genome.fa,repeats.bed}`, then
`DATASET=<dataset> ./minimap2_bench.sh pbsim-reads [depth]` fills in
`reads.fa`/`truth.tsv` (PBSIM 1.0.3 CLR reads as for D1; depth defaults to 5).
`environment.yml`/`setup.sh` set up the `minimap2_bench` conda env with `pbsim2`
and build PBSIM 1.0.3 (the paper's simulator); `lib/pbsim_reads.py` drives
either (`--simulator`) and parses the truth from the simulator's `.maf` output.
`make_dataset.py`'s read simulator is a lightweight stand-in, used for `smoke`.

## `.env` setup

```bash
cp .env.example .env   # from benchmark/minimap2/; then set MINIMAP2_BENCH_ROOT
```

| Variable | Default | Purpose |
|---|---|---|
| `MINIMAP2_BENCH_ROOT` | *(required)* | Root for `data/` and `results/` (SSD, off-repo -- e.g. `/localdata/<user>/minimap2`) |
| `MINIMAP2_DATA_DIR` | `$MINIMAP2_BENCH_ROOT/data` | Datasets |
| `MINIMAP2_RESULTS_DIR` | `$MINIMAP2_BENCH_ROOT/results` | `runs/` (per-run outputs), `csvs/` (exported CSVs, analysis tables, expected-density caches) and the workbook `minimap2_results.xlsx`; `experiments/full_run.sh` sets it to `$MINIMAP2_BENCH_ROOT/${RESULTS_NAME:-results}` |

CLI/env overrides take precedence; a variable already set in your shell is
never clobbered by `.env`.

## Layout

| Path | Role |
|---|---|
| `minimap2_bench.sh` | Entry point: `make-dataset`, `pbsim-reads [depth]`, `run <scheme> [seed] [seed_idx]`, `sweep [schemes...]`, `export [--datasets a,b]`, `list` |
| `lib/env.sh` | Loads `.env`, requires `MINIMAP2_BENCH_ROOT` |
| `lib/constants.sh` | Defaults: k = the preset's (`preset_default_k`), `W=10`, `MOD_T` = 4 + (k-4) mod w, `THREADS=16`, `DATASET=smoke`, optional `MID_OCC` (fixed map-time `-f` cap; results under `k<k>_w<w>_f<cap>/`), `$PYTHON`; the seeds (`$ALL_SEEDS`) from `config.yaml` |
| `lib/schemes.sh` | Scheme -> CLI flags (`--scheme`, `--order-seed`, ...); the default sweep is every scheme. Table parameters: `GM_TABLE` (greedymini, default `w10_k15_canonical`), `MOD_GM_TABLE` (mod-minimizer-greedymini, default `w20_k9_canonical`); polarity of Φ: `PROJECTION` for a single `run` (default `at0-cg1`), `PROJECTIONS` for `sweep` (default both). greedymini's seed toggles its Gamma tie-break |
| `lib/make_dataset.py` | Synthetic genome + repeat + simulated reads generator |
| `lib/run.sh` | One (scheme, seed, seed_idx) run: build, map, score, `meta.json`, under `runs/<dataset>/k<k>_w<w>[_f<cap>]/<scheme>[/<params>]/<seed_idx>/` (`scheme_names.scheme_dir_name`: e.g. `greedymini/k15_at1-cg0`, `mod-minimizer-xor/t9`, `mod-minimizer-greedymini/k9_at0-cg1_t9`) |
| `lib/score_accuracy.py` | PAF vs. `truth.tsv` accuracy/identity scoring (mapped fraction only for real reads) |
| `lib/pbsim_reads.py` | Simulated reads + truth for a dataset (PBSIM 1.0.3 or PBSIM2) |
| `lib/sample_reads.py` | Uniform reservoir sample of N reads of a minimum length (optionally mean quality) from a FASTQ stream (D3) |
| `lib/write_meta.py` | Assembles one run's `meta.json` from `/usr/bin/time -v` logs + accuracy JSON |
| `lib/export.py` | Aggregates `runs/**/meta.json` into `csvs/minimap2_{overview,datasets,summary,schemes,runs}.csv` and `minimap2_results.xlsx` (those sheets plus `density_effect`, `density_effect_points` and `significance`, from `benchmark/common/`). `summary` has `{metric}_ratio` vs. native; the `overview` sheet documents every column. `smoke` is excluded |
| `lib/sketch_kmers.py` | Counts the k-mers Minimap2 sketches per genome, HPC or raw (density_factor's denominator), cached as `genome.fa.{hpc,raw}_k<k>.json` |
| `lib/expected_density.py` | Expected `density_factor` per scheme variant for a dataset's k, w, t and HPC setting, counted on random DNA by Minimap2's own `mm_sketch` / `mm_sketch_scheme` ([`benchmark/common/expected_density/`](../common/expected_density/README.md)); feeds the `schemes` sheet, cached as `csvs/minimap2_density_expected_k<k>_w<w>_t<mod_t>_{hpc,nohpc}.csv` |
| `experiments/full_run.sh` | The paper's runs (see Usage) |
| `experiments/profile_run.sh` | Mapping-phase profile of a finished full run (adds `map_prof_*_pct` to the runs sheet) |
| `experiments/sweep_schemes.yaml` | Documents the default sweep's shape (not read by `minimap2_bench.sh`) |
| `tests/` | Offline pytest: schemes build and map with sane densities, export, and the greedymini / mod-minimizer sketches against Python references |
| `environment.yml`, `setup.sh` | Conda env for `pbsim2`/`samtools`, plus PBSIM 1.0.3 built from source |
| `analysis/` | Scheme analysis (`run_seed_experiments.sh`, `seed_stats.py`, `accuracy_paired.py`, `stats.py`): per-read seeds, anchors and time, timing replicates, paired accuracy tests and statistics; writes `csvs/analysis_*.csv` and `csvs/analysis/<dataset>/` |

## Usage

Run from `benchmark/minimap2/`:

```bash
./minimap2_bench.sh make-dataset            # regenerate data/smoke/ (built-in simulator)
./minimap2_bench.sh list                    # scheme variants the default sweep covers
./minimap2_bench.sh sweep                   # smoke: every scheme x every seed (skips runs with a meta.json)
DATASET=chm13_chrX ./minimap2_bench.sh sweep             # one paper dataset, Minimap2's own mid_occ
DATASET=chm13_chrX MID_OCC=1244 ./minimap2_bench.sh sweep  # same, map-time repeat cap fixed (-f 1244)
DATASET=chm13_chrX ./minimap2_bench.sh run xor-pcg-hash 0x123   # one scheme/seed
FORCE=1 ./minimap2_bench.sh sweep           # rerun complete runs
./minimap2_bench.sh export                  # -> $MINIMAP2_RESULTS_DIR/csvs/*.csv + minimap2_results.xlsx

python3 -m pytest tests/ -q                 # offline tests (Python >= 3.9)
```

`MOD_T=...` overrides t; t is part of the run directory name (`.../t<t>/`), so
runs with different t never collide.

**The paper's runs** are `experiments/full_run.sh`:

```bash
RESULTS_NAME=results THREADS=16 experiments/full_run.sh
```

Defaults: `DATASETS="chm13_chrX:3000 grch38_wg:500 chm13_rel3_ont:1000"`
(dataset:size of the read sample of the scheme analysis), `THREADS=16`,
`REPS=5`. Per dataset it runs `sweep` with Minimap2's own mid_occ, then `sweep`
with `MID_OCC` fixed to native's effective threshold (every scheme filters
repeats by the same rule), then `analysis/run_seed_experiments.sh`; finally
`analysis/accuracy_paired.py`, `analysis/stats.py` and `export`. It sets
`MINIMAP2_RESULTS_DIR=$MINIMAP2_BENCH_ROOT/${RESULTS_NAME:-results}`, overriding
`.env`. Re-running skips runs that already have a `meta.json`.

## Measurement protocol

Index build and mapping are timed under the protocol in
[`../README.md`](../README.md): node-quiet wait, pinning to `$THREADS` cores + 2
spare on one NUMA node, one thread per core, memory bound to that node (lists in
`cpus_build.txt`/`cpus_map.txt`, exported as `build_pinned_cpus` /
`map_pinned_cpus`; step records `step_*.json`, exported as `build|map_node`,
`_huge_fraction`, `_foreign_busy_cores`, ...). `sweep` is seed-major (all
schemes for seed 0, then seed 1, ...) and runs every scheme for every seed,
`native` included (its repeats build identical indexes and give the baseline
replicates). Per step: wall and CPU seconds (`build_cpu_s`, `map_cpu_s`), RSS,
major page faults; mapping time splits into `map_load_s` (Minimap2's "loaded/built
the index" stamp) and `map_run_s`. The per-read seed experiments use `-t 1`.

Indexes are 5-8 GB each on the human genome. Only seed 0's `index.mmi` of each
run is kept (the scheme analysis reads it; the paper's runs keep 44, about
165 GB); the other seeds' are deleted after mapping, their size kept in
`index_bytes.txt`. `KEEP_ALL_INDEXES=1` keeps them all (about 5x the space).

## Known limitations

**The mod-minimizer helps only when w < k** (here k = 19 or 15, w = 10).

**mod-minimizer-xor is canonical**: plain mod-sampling (the k-mer at offset x
mod w, x the minimal t-mer's position) is not reverse-complement symmetric, so
Minimap2, Kraken2 and SSHash share the canonical rule of
`common/minimizer-schemes/mod_minimizer_rule.h`: with t ≡ k (mod w), a unique
minimal t-mer selects offset x mod w; a tie at the minimum takes the offset
closest to the window's center, and of two equally close ones the left if the
window is <= its reverse complement, else the right (the center rule of
Cologni & Pibiri, as in SSHash). A read and its reverse complement then select
the same loci, and the selection never moves backwards (the scheme is forward).
`analysis/canonical_mod.py` measures density and strand agreement of plain
mod-sampling, a G/T-strand canonical variant and simd-minimizers' canonical
minimizer.

`minimap2/scheme.cpp`'s `sketch_mod_minimizer` keeps the window's t-mer ranks
in a ring buffer with the leftmost/rightmost minimum, rescanning only when the
leftmost minimum leaves the window; `tests/test_mod_minimizer_sketch.py`
(run by `make -C minimap2 test` from the repository root) checks the compiled
sketch against a brute force, including N gaps and palindromic k-mers.
