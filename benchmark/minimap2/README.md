# minimap2 benchmark harness

Compares minimap2's minimizer orders (native, xor, xor-pcg-hash,
greedymini and mod-minimizer-xor; greedymini's AT/CG projection
(`--gm-projection at0-cg1|at1-cg0`) and GreedyMini table (`gm_table`, whose k
is `gm_k`) are separate parameters, not part of the order name, and each
projection is its own set of runs/rows; see
[`orders/README.md`](../../orders/README.md)) the way Jain et al. 2020
("Weighted minimizer sampling improves long read mapping" -- Winnowmap,
Bioinformatics 36(Supp_1):i111-i118) evaluated Winnowmap's weighted
minimizers against minimap2's stock (uniform-hash) ones: minimizer density,
index build/mapping time and memory, and mapping accuracy specifically on
reads originating inside a long, near-identical repeat.

## What is measured

Per (dataset, order) run, from `minimap2 -d` (index build) and a plain
mapping pass:

- **Minimizer density**: `density_factor = total_occurrences * (w+1) /
  n_sketch_kmers`, on minimap2's own sketching scale: `n_sketch_kmers`
  counts the k-mers minimap2 sketches (`lib/sketch_kmers.py`: length minus
  k-1 per N-free segment), homopolymer-compressed under an HPC preset
  (`-x map-pb`) and raw otherwise (`-x map-ont`). minimap2 `-d`'s own "average spacing"
  divides by the *raw* length instead, understating density by the genome's
  compression ratio (~0.70 on human) plus N gaps -- kept only as the runs
  sheet's `density_factor_raw_bp`. total_occurrences (distinct minimizers x
  average occurrences, both from minimap2's `-d` stderr) is total *occurrences*, not
  `distinct minimizers` (also parsed, reported separately as
  `distinct_minimizers`/`avg_occurrences`, informational only) -- distinct
  minimizers collapses every repeated occurrence of the same minimizer key
  to one count, so on a genome with any repeat content it always undercounts
  total selected positions, and undercounts *more* for an order that's
  better at consistently re-selecting the same minimizer across a repeat's
  copies (exactly the property these custom orders are meant to have). That
  undercount can push a distinct-minimizer-based ratio below the proven
  theoretical minimum of `density_factor = 1` (every scheme selects at least
  one minimizer per `w+1` window) -- using total occurrences avoids this.
  Matches `benchmark/common/density.py`'s own `density_factor` definition
  used for kmc/sshash; a uniformly random order has expected
  density_factor -> 2 as w grows, GreedyMini/mod-minimizer push it below
  that -- see that module's docstring.
- **Index build time/peak memory** and **mapping time/peak memory**
  (`/usr/bin/time -v`, same wrapper the other three harnesses use).
- **Mapping accuracy**: fraction of simulated reads whose primary PAF
  alignment overlaps their known true origin by at least 10% of the *union*
  of the true and reported intervals (`lib/score_accuracy.py`'s
  `frac_correct`) -- the exact criterion both Winnowmap and minimap2's own
  papers use, so these numbers are comparable to their published figures.
  Reported both overall and restricted to reads originating inside the
  repeat (`frac_correct_in_repeat` -- the papers' own headline comparison),
  plus mean per-base identity from PAF's `dv:f:`/`de:f:` tag. A stricter
  50%-of-true-span threshold is also reported (`frac_correct_strict50`) for
  continuity with this harness's own earlier (pre-literature-alignment)
  numbers. Real reads (D3) have no truth: only the mapped fraction (the
  paper's "unmapped reads") and identity, as in the paper.

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

k is each preset's own default, as the paper ran minimap2 with `-cx map-pb`
/ `-cx map-ont` (`lib/constants.sh:preset_default_k`; `K=` overrides); every
order uses that k (GreedyMini's k15 table is centre-extended to k=19), and
mod-minimizer's t = 4 + (k-4) % w (9 at k=19, 5 at k=15). Each dataset's
preset is a one-word `data/<dataset>/preset` file (absent = `map-pb`;
`$PRESET` overrides it), recorded in every run's `meta.json`.

To regenerate (`$E` = the `minimap2_bench` conda env; `setup.sh` builds PBSIM
1.0.3 into `$E/opt/pbsim1`):

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
  | python3 "$OLDPWD/lib/sample_reads.py" --n 1000 --min-len 1000 --min-mean-q 7 --seed 1 --out-dir .
```

Any other dataset works the same way: point `minimap2_bench.sh run`/`sweep`
at `$DATASET` and provide `data/<dataset>/{genome.fa,repeats.bed}` yourself,
then run `pbsim-reads` (or `make_dataset.py`'s own simulator, genome/reads
only, no real download) to fill in `reads.fa`/`truth.tsv`.

`make_dataset.py`'s own read simulation is a lightweight, dependency-free
stand-in for a real long-read simulator. `environment.yml`/`setup.sh` set up
a conda env (`minimap2_bench`) with `pbsim2` and build PBSIM 1.0.3 (the
paper's simulator); `lib/pbsim_reads.py` drives either (`--simulator`) and
parses ground truth from the simulator's `.maf` output -- its docstring has
the details.

## `.env` setup

```bash
cp benchmark/minimap2/.env.example benchmark/minimap2/.env   # then edit if needed
```

| Variable | Default | Purpose |
|---|---|---|
| `MINIMAP2_BENCH_ROOT` | *(required)* | Root for `data/` and `results/` (SSD, off-repo -- e.g. `/localdata/<user>/minimap2`) |
| `MINIMAP2_DATA_DIR` | `$MINIMAP2_BENCH_ROOT/data` | Datasets |
| `MINIMAP2_RESULTS_DIR` | `$MINIMAP2_BENCH_ROOT/results` | `runs/` (per-run outputs), `csvs/` (exported CSVs, analysis tables, expected-density caches), and the workbook `minimap2_results.xlsx` directly in it |

CLI/env overrides take precedence; a variable already set in your shell is
never clobbered by `.env`.

## Layout

| Path | Role |
|---|---|
| `minimap2_bench.sh` | Entry point: `make-dataset`, `run <order> [seed] [seed_idx]`, `sweep [orders...]`, `export [--datasets a,b]`, `list` |
| `lib/env.sh` | Loads `.env`, requires `MINIMAP2_BENCH_ROOT` |
| `lib/constants.sh` | Defaults (k = the dataset preset's own via `preset_default_k`, `w=10`, mod_t = 4 + (k-4) % w; `$PYTHON` resolution; optional `MID_OCC` = fixed map-time `-f` cap, results under `k<k>_w<w>_f<cap>/`) + `config.yaml` seed loading (`$ALL_SEEDS`) |
| `lib/orders.sh` | Order name -> CLI flags (`ALL_ORDERS` == `DEFAULT_ORDERS`); `sweep` runs every order across the full seed list, `native` included (its repeats build identical indexes and give the baseline replicates; greedymini's seed toggles its Gamma tie-break). greedymini's parameters are `GM_TABLE` (default `w10_k15_canonical`) and `PROJECTION` (`at0-cg1`/`at1-cg0`, for a single `run`); `sweep` runs greedymini once per projection in `PROJECTIONS` (default both). Run-dir names come from `order_names.order_dir_name` |
| `lib/make_dataset.py` | Synthetic genome + repeat + simulated reads generator |
| `lib/run.sh` | Single (order, seed, seed_idx) run: build, map, score, `meta.json` (results land under `runs/<dataset>/k<k>_w<w>/<order>[/<params>]/<seed_idx>/`, `<params>` = `k<gm_k>_<projection>` for greedymini, `t<mod_t>` for mod-minimizer, e.g. `greedymini/k15_at1-cg0`; `meta.json` carries the separate `order, gm_k, projection, gm_table, mod_t` fields) |
| `lib/score_accuracy.py` | PAF vs. `truth.tsv` accuracy/identity scoring (mapped fraction only for real reads) |
| `lib/pbsim_reads.py` | Simulated reads + truth for a dataset (PBSIM 1.0.3 or PBSIM2) |
| `lib/sample_reads.py` | Uniform reservoir sample of N reads >= L bp (optionally mean quality >= Q) from a FASTQ stream (the paper's D3) |
| `lib/write_meta.py` | Assembles one run's `meta.json` from `/usr/bin/time -v` logs + accuracy JSON |
| `lib/export.py` | Aggregates every `runs/**/meta.json` into one CSV per sheet (`csvs/minimap2_{overview,datasets,summary,orders,runs}.csv`) plus `minimap2_results.xlsx` written by `benchmark/common/xlsx_book.write_book` (same 5 sheets, then `density_effect` / `density_effect_points`, `significance` (`csvs/significance_minimap2.csv`, `benchmark/common/significance.py`) and `csvs/density_effect_minimap2.csv` from `benchmark/common/density_effect.py`: Pearson r of each metric and the density factor as ratios to native's run of the same seed, over the fixed-threshold block only; `overview` = per-column legend, with a category-band header (no fills), dataset dividers/merges (`order, gm_k, projection, gm_table, mod_t` are separate id columns), and `summary`'s `{metric}_ratio` vs. native with green/white/red conditional coloring, plus, as in sshash's export, correlation rows right under `summary`'s table (Excel `CORREL` of every `{metric}_ratio` column with the density and distinct-minimizer ratios); the `orders` sheet is the shared layout of `benchmark/common/orders_sheet.py`: one row per (order variant, `k`, `w`, `L` = k+w-1) with `mod_t` of that combination, `expected_density_factor_mean/std` for the row's own k/w/t, then per-dataset `<dataset>_particular_density_factor_mean/std` filled on the rows of the dataset's own (k, w, L, mod_t) and export fails if an applicable cell is empty; `smoke` excluded throughout; every sheet with `k` and `w` has `L` = k+w-1 right after `w`: datasets, summary, runs, orders; column order of `summary`/`runs`: dataset, n_reads (, n_reads_in_repeat), k, w, L, order, gm_k, projection, gm_table, mod_t, then app parameters `mid_occ_threshold` (, `threads`), then repeats (summary `n_repeats` / runs `seed, seed_idx`), then metrics banded density, space (`index_size_gb`), total, build, map (`map_load_s`/`map_run_s`), accuracy, and on runs only the build phase `build_collect_s` / `build_collect_pct` (from build.log), the mapping profile `map_prof_{sketch,seed,chain,align}_pct` (after `experiments/profile_run.sh`) and `telemetry` (`build_pinned_cpus`, `map_pinned_cpus`); `analysis_*.csv` carry the hidden block key `mid_occ_cap` right after `mod_t`) |
| `lib/sketch_kmers.py` | Counts the k-mers minimap2 sketches per genome, HPC or raw (density_factor's denominator), cached as `genome.fa.{hpc,raw}_k<k>.json` |
| `lib/expected_density.py` | Theoretical expected `density_factor` per order (Monte Carlo on i.i.d. random DNA; thin adapter over the shared `benchmark/common/expected_density_driver.py`) for a dataset's k, w, t and HPC setting; feeds the `orders` sheet's `expected_density_factor_mean/std` (one call per (k, w, t) combination), cached as `csvs/minimap2_density_expected_k<k>_w<w>_t<mod_t>_{hpc,nohpc}.csv` (columns `order, gm_k, projection, gm_table, mod_t, k, w, L, ...`; one row per order variant: greedymini once per projection) |
| `experiments/sweep_orders.yaml` | Documents the default sweep's shape (see that file's own header for why it isn't wired through the shared YAML engine yet) |
| `tests/test_order_sanity.py` | Offline pytest: every order builds/maps, densities are non-degenerate, mapped-read sets agree across orders |
| `environment.yml`, `setup.sh` | Conda env for `pbsim2`/`samtools`, plus PBSIM 1.0.3 built from source |
| `analysis/` | Scheme analysis (writes `results/csvs/analysis_*.csv`, per-read tables under `results/csvs/analysis/<ds>/{seeds,timing}/`, sampled reads under `results/analysis/<ds>/`): paired per-read accuracy tests, per-read seed/anchor/time stats, timing replicates, statistics (`accuracy_paired.py`, `seed_stats.py`, `run_seed_experiments.sh`, `stats.py`) |

## Usage

```bash
cd benchmark/minimap2
./minimap2_bench.sh make-dataset           # regenerate data/smoke/ (built-in simulator)
./setup.sh && ./minimap2_bench.sh pbsim-reads 5   # swap in PBSIM 1.0.3 CLR reads (paper's D1 setup) at 5x depth
./minimap2_bench.sh list                   # orders the default sweep covers
./minimap2_bench.sh sweep                  # build+map+score every default order (skips runs that already have meta.json)
FORCE=1 ./minimap2_bench.sh sweep          # rerun them anyway
# MOD_T=... overrides mod-minimizer's t (default 4 + (k-4) % w); t is not in the run
# directory name, so a complete run made with another t is refused, never reused or overwritten
./minimap2_bench.sh run xor-pcg-hash 0x123  # just one order/seed
MID_OCC=1244 ./minimap2_bench.sh sweep    # same sweep, map-time repeat cap fixed (-f 1244) for every order
./minimap2_bench.sh export                 # -> $MINIMAP2_RESULTS_DIR/csvs/minimap2_{overview,datasets,summary,orders,runs}.csv + minimap2_results.xlsx

# offline sanity tests (needs a python3 new enough for these type hints --
# system python3 here is 3.6; use the conda env's or /usr/bin/python3.11)
/usr/bin/python3.11 -m pytest tests/ -q
```

## Measurement protocol

Index build and mapping are timed under the protocol in [`../README.md`](../README.md): node-quiet wait,
pinning (`$THREADS` cores + 2 spare, one NUMA node, one thread per core; lists in `cpus_build.txt`
and `cpus_map.txt`, exported as `build_pinned_cpus` / `map_pinned_cpus`; with its memory bound to the slot's node and a step record (`step_*.json`, exported as `build|map_node`, `_huge_fraction`, `_foreign_busy_cores`, ... -- see the protocol); the
per-read seed experiments use `-t 1`; `sweep` runs seed-major: all orders for seed 0, then seed 1, ...) and per-step wall and CPU seconds (`build_cpu_s`, `map_cpu_s`), RSS, major page faults and reads
(runs sheet). The mapping time splits into `map_load_s` (minimap2's own "loaded/built the index" stamp) and
`map_run_s`. `sweep` runs every order for every seed, `native` included (its repeats build identical indexes
and give the baseline replicates). Every `index.mmi` is kept (the sweep needs about 620 GB);
`experiments/full_run.sh` runs both repeat-filter conditions, the seed experiments, the statistics and the export.

## Known limitations

**mod-minimizer is only relevant when `w < k`**: it improves on a classic
minimizer's density bound specifically in the regime where the sampled
substring `k` is large *relative to* the window `w`. This harness's default
sweep (`k=19` or `k=15`, `w=10`) is in that regime; don't expect the same benefit (or
run it at all) for `w >= k`.

**mod-minimizer-xor is strand-symmetric**: plain mod-sampling ("reduce the
winning t-mer's position mod w") is not reverse-complement symmetric, so
minimap2's mod-minimizer-xor uses a canonical variant shared with kraken2 and
sshash (`common/minimizer-orders/mod_minimizer_rule.h`): with t == k (mod w)
a unique minimal t-mer samples offset x mod w, and a tie at the minimum takes
the offset closest to the window centre, two equally close ones by the window
<= its reverse complement (the centre rule of Cologni & Pibiri, as in
SSHash). A read and its reverse complement then select the same loci, as for
the per-locus canonical orders, and the scheme is forward.
`analysis/canonical_mod.py` compares plain mod-sampling with the earlier
G/T-strand variant and simd-minimizers' canonical minimizer (density, strand
agreement).

`minimap2/order.cpp`'s `sketch_mod_minimizer` uses kraken2's single
monotone deque over all live t-mers; `minimap2/validate_mod_minimizer.py`
(run by `make -C minimap2 test`) checks the compiled sketch against a brute
force, including N gaps and palindromic k-mers.
