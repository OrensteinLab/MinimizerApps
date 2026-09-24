# MinimizerApps

A monorepo folding three forked bioinformatics tools together so their shared
minimizer/signature-order logic and benchmark harnesses live in one place,
while each tool's own upstream-tracked source stays isolated from that
shared/harness code.

## Layout

| Path | What it is |
|---|---|
| `kraken2/` | Fork of [DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) — taxonomic sequence classification. Contains only upstream-adjacent changes (the order-ranking integration); see [`kraken2/README.md`](kraken2/README.md). |
| `kmc/` | Fork of [refresh-bio/KMC](https://github.com/refresh-bio/KMC) — k-mer counting. Same rule: only order-implementation changes relative to upstream; see [`kmc/README.md`](kmc/README.md). |
| `sshash/` | Fork of [jermp/sshash](https://github.com/jermp/sshash) — compact k-mer dictionary index. Same rule: only order-implementation changes relative to upstream; see [`sshash/README.md`](sshash/README.md). |
| `common/minimizer-orders/` | The minimizer/signature-order implementation shared by all three apps (GreedyMini LUT, PCG scramble, reverse-complement, `.gm` table tools) — see [`orders/README.md`](orders/README.md). |
| `orders/` | The canonical `.gm`/`.txt` order tables all three apps load, and the order options/CLI-flag reference — [`orders/README.md`](orders/README.md). |
| `benchmark/{kraken2,kmc,sshash,common}/` | Each app's benchmark/scoring harness (100% lab-authored, not part of any upstream) plus the shared YAML multi-config experiment-sweep engine they're all built on. |
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

# sshash
make -C sshash all
```

## Running a benchmark sweep

Each app's benchmark harness reads a YAML file of `defaults:` + `configs:`
(see `benchmark/common/experiments_lib.py`) and runs each config as one
sweep entry.

```bash
# kraken2 (needs KRAKEN2_BENCHMARK_DIR or benchmark/kraken2/.env; see benchmark/kraken2/README.md)
benchmark/kraken2/run_config_sweep.sh list
benchmark/kraken2/run_config_sweep.sh sweep_schemes

# kmc (needs KMC3_BENCH_ROOT; see benchmark/kmc/README.md)
benchmark/kmc/kmc3_bench.sh experiment list
benchmark/kmc/kmc3_bench.sh experiment sweep_config

# sshash (needs SSHASH_BENCH_ROOT or benchmark/sshash/.env; see benchmark/sshash/README.md)
benchmark/sshash/sshash_bench.sh experiment list
benchmark/sshash/sshash_bench.sh experiment sweep_orders
```

## More docs

- [`kraken2/README.md`](kraken2/README.md), [`kmc/README.md`](kmc/README.md), [`sshash/README.md`](sshash/README.md) — each app's own documentation.
- [`orders/README.md`](orders/README.md) — the minimizer/signature order options, their shared implementation, and the CLI-flag mapping between the three apps.
- [`benchmark/kraken2/README.md`](benchmark/kraken2/README.md), [`benchmark/kmc/README.md`](benchmark/kmc/README.md), [`benchmark/sshash/README.md`](benchmark/sshash/README.md) — each benchmark harness's setup and usage.
- [`UPSTREAMS.md`](UPSTREAMS.md) — tracking and re-pulling from each app's real upstream.
