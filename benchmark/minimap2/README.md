# minimap2 benchmark harness

Compares minimap2's minimizer orders (native, xor-hash, xor-pcg-hash,
greedymini, and mod-minimizer-{xor-hash,pcg-xor-hash} -- no
mod-minimizer-greedymini here, unlike kraken2/sshash; see
[`orders/README.md`](../../orders/README.md)) the way Jain et al. 2020
("Weighted minimizer sampling improves long read mapping" -- Winnowmap,
Bioinformatics 36(Supp_1):i111-i118) evaluated Winnowmap's weighted
minimizers against minimap2's stock (uniform-hash) ones: minimizer density,
index build/mapping time and memory, and mapping accuracy specifically on
reads originating inside a long, near-identical repeat.

## What is measured

Per (dataset, order) run, from `minimap2 -d` (index build) and a plain
mapping pass:

- **Minimizer density**: `distinct minimizers` (parsed from minimap2's own
  `-d` stderr) and the derived `density_factor = distinct_minimizers /
  genome_len * (w+1)` (the same quantity `benchmark/common/density.py`
  computes for kmc/sshash; a uniformly random order has expected
  density_factor -> 2 as w grows, GreedyMini/mod-minimizer push it below
  that -- see that module's docstring).
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
  numbers.

## Datasets

**`smoke`** (default, `$DATASET` unset): a synthetic ~300kb genome with a
20kb region duplicated verbatim at a distant location (`lib/make_dataset.py`,
no download dependency, fast to regenerate) plus simulated long reads with
known ground truth (`truth.tsv`), 30% deliberately sampled from inside the
repeat. Good for fast iteration and CI-scale sanity checks; too small and
too artificial for real accuracy conclusions.

**`chm13_chrX`**: the real dataset -- see "Results" below. To regenerate:

```bash
mkdir -p "$MINIMAP2_DATA_DIR/chm13_chrX" && cd "$MINIMAP2_DATA_DIR/chm13_chrX"
curl -sL -o genome.fa.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.chrX_v0.7.fasta.gz
gunzip genome.fa.gz  # verify md5 89b3dd61db66177dd830527b920956fa first
printf 'chrX_v0.7\t57828561\t60934693\tcentromere_dxz1_array\n' > repeats.bed  # Winnowmap Table 3's own coordinates
cd - && ./minimap2_bench.sh pbsim-reads P6C4.model 5   # PacBio model, 5x depth (matches Winnowmap's D1)
```

Any other dataset works the same way: point `minimap2_bench.sh run`/`sweep`
at `$DATASET` and provide `data/<dataset>/{genome.fa,repeats.bed}` yourself,
then run `pbsim-reads` (or `make_dataset.py`'s own simulator, genome/reads
only, no real download) to fill in `reads.fa`/`truth.tsv`.

`make_dataset.py`'s own read simulation is a lightweight, dependency-free
stand-in for a real long-read simulator. `environment.yml`/`setup.sh` set up
a conda env (`minimap2_bench`) with `pbsim2` (a real ONT/PacBio error-model
simulator); once set up, `./minimap2_bench.sh pbsim-reads` replaces
`reads.fa`/`truth.tsv` with pbsim2's own simulated reads (still sampled from
the same genome+repeat, ground truth parsed straight from pbsim2's own
`.maf` output) -- `lib/pbsim_reads.py`'s docstring has the details. This is
what actually produced the numbers below.

## `.env` setup

```bash
cp benchmark/minimap2/.env.example benchmark/minimap2/.env   # then edit if needed
```

| Variable | Default | Purpose |
|---|---|---|
| `MINIMAP2_BENCH_ROOT` | *(required)* | Root for `data/` and `results/` (SSD, off-repo -- e.g. `/localdata/<user>/minimap2`) |
| `MINIMAP2_DATA_DIR` | `$MINIMAP2_BENCH_ROOT/data` | Datasets |
| `MINIMAP2_RESULTS_DIR` | `$MINIMAP2_BENCH_ROOT/results` | Per-run outputs + aggregated CSV |

CLI/env overrides take precedence; a variable already set in your shell is
never clobbered by `.env`.

## Layout

| Path | Role |
|---|---|
| `minimap2_bench.sh` | Entry point: `make-dataset`, `run <order> [seed]`, `sweep [orders...]`, `export`, `list` |
| `lib/env.sh` | Loads `.env`, requires `MINIMAP2_BENCH_ROOT` |
| `lib/constants.sh` | Defaults (`k=15 w=10 mod_t=5`, `$PYTHON` resolution) + `config.yaml` seed loading |
| `lib/orders.sh` | Order name -> CLI flags (`ALL_ORDERS` == `DEFAULT_ORDERS`; no `mod-minimizer-greedymini` -- see README below) |
| `lib/make_dataset.py` | Synthetic genome + repeat + simulated reads generator |
| `lib/run.sh` | Single (order, seed) run: build, map, score, `meta.json` |
| `lib/score_accuracy.py` | PAF vs. `truth.tsv` accuracy/identity scoring |
| `lib/write_meta.py` | Assembles one run's `meta.json` from `/usr/bin/time -v` logs + accuracy JSON |
| `lib/export.py` | Aggregates every `runs/**/meta.json` into `minimap2_runs.csv` |
| `experiments/sweep_orders.yaml` | Documents the default sweep's shape (see that file's own header for why it isn't wired through the shared YAML engine yet) |
| `tests/test_order_sanity.py` | Offline pytest: every order builds/maps, densities are non-degenerate, mapped-read sets agree across orders |
| `environment.yml`, `setup.sh` | Conda env for `pbsim2`/`samtools` (see "Default dataset" above) |

## Usage

```bash
cd benchmark/minimap2
./minimap2_bench.sh make-dataset           # regenerate data/smoke/ (built-in simulator)
./setup.sh && ./minimap2_bench.sh pbsim-reads R94.model 8  # swap in real pbsim2 ONT reads at 8x depth
./minimap2_bench.sh list                   # orders the default sweep covers
./minimap2_bench.sh sweep                  # build+map+score every default order
./minimap2_bench.sh run xor-pcg-hash 0x123  # just one order/seed
./minimap2_bench.sh export                 # -> $MINIMAP2_RESULTS_DIR/minimap2_runs.csv

# offline sanity tests (needs a python3 new enough for these type hints --
# system python3 here is 3.6; use the conda env's or /usr/bin/python3.11)
/usr/bin/python3.11 -m pytest tests/ -q
```

## Results

### Real chromosome (T2T CHM13 chrX v0.7, 154.3 Mbp)

`data/chm13_chrX/` -- the exact assembly and repeat-region coordinates
Winnowmap's own D1/D3 evaluation used (Jain et al. 2020): downloaded from
[github.com/marbl/CHM13](https://github.com/marbl/CHM13) (md5-verified), the
centromeric repeat span (57,828,561-60,934,693) taken directly from that
paper's Table 3 caption. 85,324 real pbsim2 PacBio-model reads (P6C4, 5x
depth). Swept at `k=15` (minimap2's own default) and both `w=10` (also
minimap2's own default) and `w=50` -- Winnowmap's own default, and the
window size where its headline centromere result (3.6%->0%, 34.7%->0%
incorrectly-mapped at w=10 vs w=50) is most dramatic, so it's the window
size where custom orders' effect on this same comparison should show up
most clearly too:

| w | order | density_factor | frac_mapped | frac_correct | frac_correct_in_repeat |
|---|---|---|---|---|---|
| 10 | native | 1.193 | 0.999 | 0.991 | 0.716 |
| 10 | xor-pcg-hash | 1.157 | 0.999 | 0.991 | 0.731 |
| 10 | xor-hash | 1.080 | 0.998 | 0.992 | **0.786** |
| 10 | mod-minimizer-xor-hash | 0.850 | 0.998 | 0.989 | 0.745 |
| 10 | mod-minimizer-pcg-xor-hash | 0.859 | 0.998 | 0.989 | 0.742 |
| 10 | **greedymini** | **0.827** | 0.998 | **0.991** | 0.745 |
| 50 | native | 1.238 | 0.985 | 0.965 | 0.276 |
| 50 | xor-pcg-hash | 1.186 | 0.985 | 0.966 | 0.296 |
| 50 | xor-hash | 1.174 | 0.985 | 0.968 | 0.310 |
| 50 | **greedymini** | **1.267** | **0.986** | **0.968** | **0.337** |
| 50 | mod-minimizer-xor-hash | 1.783 | 0.686 | 0.568 | 0.195 |
| 50 | mod-minimizer-pcg-xor-hash | 1.767 | 0.685 | 0.565 | 0.195 |

At `w=10`: `greedymini` has the lowest density of all six orders with no
accuracy cost; mod-minimizer orders are similarly low-density.

At `w=50`: `greedymini` again gives the best density *and* the best overall
and in-repeat accuracy -- consistent with (a small-scale echo of) Winnowmap's
own finding that a better order matters more as sampling gets sparser.
**mod-minimizer orders collapse here** (`frac_mapped` 0.69 vs ~0.99 for
everything else) -- see "mod-minimizer's known limitation" below before
drawing any conclusion from their numbers at any `w`.

`greedymini` uses `orders/GreedyMini/w10_k15_canonical.gm` (gm_k=15, matching
`-k 15` exactly, generated with the
[GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `-mode
expected -w 10 -k 15`, per `orders/README.md`'s "Regenerating a `.gm`
table"). There is no `mod-minimizer-greedymini` for minimap2, unlike
kraken2/sshash -- it would need a second, differently-sized GreedyMini
table (gm_k == the mod-minimizer `t`, not `k`) just for that one variant,
so it was dropped instead of adding that second table.

### mod-minimizer's known limitation

While investigating the `w=50` collapse, two real bugs were found and fixed
in `minimap2/order.cpp`'s `sketch_mod_minimizer` (see that file's git
history): a tie-break that compared `(residue - window_start) % w` instead
of the true (unreduced) t-mer offset -- wrong on any rank tie, and ties are
common with xor-hash's weak avalanche -- and a per-residue deque capacity
formula that could silently overflow its fixed-size backing array for large
`w`. Both are fixed and verified against a brute-force reference (see the
commit for the exact Python cross-check) across several `(w,k,t,seed)`
combinations.

After those fixes, a **third, unfixed, fundamental issue** remains: the
mod-minimizer scheme's "reduce the winning t-mer's position mod w" rule is
provably not reverse-complement symmetric (the two t-mer offsets under
strand reversal sum to a fixed constant, but reducing each mod `w`
independently doesn't generally land on mirrored l-mer slots). Measured
directly (build an index from a 3kb real-genome sequence, query with its own
reverse complement): mod-minimizer-xor-hash shared only 4 of 140 minimizer
keys between the two strands, versus ~90%+ for xor-hash/native/greedymini,
which use per-locus canonical parsing and don't have this issue. This is the
same property `orders/README.md`'s "mod-minimizer" section already documents
sshash working around with a dual-probe/"regular parsing" fallback --
**not yet ported to minimap2**. It has little visible effect at `w=10`
(there's enough redundant minimizer coverage that mapping still mostly
succeeds) but becomes fatal at `w=50` (sparse sampling leaves no redundancy
to fall back on), which is why the `w=10` numbers above look reasonable
while `w=50`'s don't. **Treat every mod-minimizer-\* result here as
provisional** until that fallback is ported.
