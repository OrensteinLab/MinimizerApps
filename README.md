# MinimizerApps

A monorepo folding two forked bioinformatics tools together so their shared
minimizer/signature-order logic and benchmark harnesses live in one place,
while each tool's own upstream-tracked source stays isolated from that
shared/harness code.

## Layout

| Path | What it is |
|---|---|
| `kraken2/` | Fork of [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) — taxonomic sequence classification. Contains only upstream-adjacent changes (the order-ranking integration); see [`kraken2/README.md`](kraken2/README.md). |
| `kmc/` | Fork of [refresh-bio/KMC](https://github.com/refresh-bio/KMC) — k-mer counting. Same rule: only order-implementation changes relative to upstream; see [`kmc/README.md`](kmc/README.md). |
| `common/minimizer-orders/` | The minimizer/signature-order implementation shared by both apps (GreedyMini LUT, PCG scramble, reverse-complement, `.gm` table tools) — see [`orders/README.md`](orders/README.md). |
| `orders/` | The canonical `.gm`/`.txt` order tables both apps load, and the order options/CLI-flag reference — [`orders/README.md`](orders/README.md). |
| `benchmark/{kraken2,kmc,common}/` | Both apps' benchmark/scoring harnesses (100% lab-authored, not part of either upstream) plus the shared YAML multi-config experiment-sweep engine they're both built on. |
| `UPSTREAMS.md` | How each app's folder tracks and re-pulls from its real upstream — see [`UPSTREAMS.md`](UPSTREAMS.md). |

## Build and test

```bash
# kraken2
make -C kraken2 all
make -C kraken2 test

# kmc
make -C kmc all
make -C benchmark/kmc venv   # once: creates kmc/.venv, installs pytest into it
make -C benchmark/kmc test
```

## Running a benchmark sweep

Both apps' benchmark harnesses read a YAML file of `defaults:` + `configs:`
(see `benchmark/common/experiments_lib.py`) and run each config as one
sweep entry.

```bash
# kraken2 (needs KRAKEN2_BENCHMARK_DIR or benchmark/kraken2/.env; see benchmark/kraken2/README.md)
benchmark/kraken2/run_config_sweep.sh list
benchmark/kraken2/run_config_sweep.sh sweep_schemes

# kmc (needs KMC3_BENCH_ROOT; see benchmark/kmc/README.md)
benchmark/kmc/kmc3_bench.sh experiment list
benchmark/kmc/kmc3_bench.sh experiment sweep_config
```

## More docs

- [`kraken2/README.md`](kraken2/README.md), [`kmc/README.md`](kmc/README.md) — each app's own documentation.
- [`orders/README.md`](orders/README.md) — the minimizer/signature order options, their shared implementation, and the CLI-flag mapping between the two apps.
- [`benchmark/kraken2/README.md`](benchmark/kraken2/README.md), [`benchmark/kmc/README.md`](benchmark/kmc/README.md) — each benchmark harness's setup and usage.
- [`UPSTREAMS.md`](UPSTREAMS.md) — tracking and re-pulling from each app's real upstream.
