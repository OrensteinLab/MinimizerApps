# minimap2 benchmark harness

Compares minimap2's minimizer orders (native, xor-hash, xor-pcg-hash,
greedymini, and mod-minimizer-{xor-hash,pcg-xor-hash,greedymini}; see
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
  alignment overlaps their known true origin by at least 50% of the read's
  true span (`lib/score_accuracy.py`), reported both overall
  (`frac_correct`) and restricted to reads originating inside the injected
  repeat (`frac_correct_in_repeat` -- the paper's own headline comparison),
  plus mean per-base identity from PAF's `dv:f:`/`de:f:` tag.

## Default dataset

A synthetic ~300kb genome with a 20kb region duplicated verbatim at a
distant location (`lib/make_dataset.py`, no download dependency, fast to
regenerate) plus ~400 simulated long reads (2-10kb, ONT-like ~6%
substitution + 3%+3% indel error, 30% deliberately sampled from inside the
repeat) with known ground truth (`truth.tsv`). Swap in a real genome later
by pointing `minimap2_bench.sh run`/`sweep` at a different `$DATASET` whose
`data/<dataset>/{genome.fa,reads.fa,truth.tsv}` you provide the same way
(e.g. a human chromosome with real segmental duplications, closer to the
paper's own scale) -- see that script's own docstring for the format.

Read simulation is a lightweight, dependency-free stand-in for a real
long-read simulator; `environment.yml`/`setup.sh` set up a conda env with
`pbsim2` (a real ONT/PacBio error-model simulator) for when that matters
more than fast iteration.

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
| `lib/orders.sh` | Order name -> CLI flags (`ALL_ORDERS`, `DEFAULT_ORDERS` -- greedymini variants excluded until their `.gm` tables are supplied) |
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
./minimap2_bench.sh make-dataset          # regenerate data/smoke/
./minimap2_bench.sh list                  # orders the default sweep covers
./minimap2_bench.sh sweep                  # build+map+score every default order
./minimap2_bench.sh run xor-pcg-hash 0x123 # just one order/seed
./minimap2_bench.sh export                 # -> $MINIMAP2_RESULTS_DIR/minimap2_runs.csv

# offline sanity tests (needs a python3 new enough for these type hints --
# system python3 here is 3.6; use the conda env's or /usr/bin/python3.11)
/usr/bin/python3.11 -m pytest tests/ -q
```

`greedymini` and `mod-minimizer-greedymini` are excluded from `list`/`sweep`
until `orders/GreedyMini/w10_k15_canonical.gm` (k=15,w=10, matching this
harness's fixed k/w) and a `gm_k=5` table for the mod-minimizer case are
supplied; `lib/orders.sh`'s `order_gm_file()` documents the expected
filenames. Run them directly once available:

```bash
./minimap2_bench.sh run greedymini
./minimap2_bench.sh run mod-minimizer-greedymini
```
