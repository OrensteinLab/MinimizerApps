# MinimizerApps

Canonical sketching schemes with low-density orders (`xor`, `xor-pcg-hash`,
the canonical GreedyMini, and the canonical mod-minimizer on an `xor` or a
GreedyMini order) integrated into four high-throughput sequencing methods, and
the harnesses that benchmark them end to end:

- **Kraken2**: taxonomic classification.
- **KMC3**: k-mer counting.
- **SSHash**: k-mer indexing.
- **Minimap2**: read mapping.

Each method is a fork whose upstream-tracked source stays apart from the shared
order code and the benchmark harnesses.

## Layout

| Path | What it is | Docs |
|---|---|---|
| `kraken2/`, `kmc/`, `sshash/`, `minimap2/` | The four forks, each with only the changes needed to plug in the orders | [`kraken2`](kraken2/README.md), [`kmc`](kmc/README.md), [`sshash`](sshash/README.md), [`minimap2`](minimap2/README.fork.md) |
| `orders/` | The order tables (`.txt` rank lists; the `.gm` tables are built from them), the order names and options, and the CLI flags of each method | [`orders/README.md`](orders/README.md) |
| `common/minimizer-orders/` | The order code shared by all four methods: the GreedyMini tables, the PCG scramble and the mod-minimizer rule with its tie handling | [`common/minimizer-orders/README.md`](common/minimizer-orders/README.md) |
| `benchmark/` | The shared measurement protocol, the full run of all four methods (`full_run_all.sh`), the exporters and the analyses | [`benchmark/README.md`](benchmark/README.md) |
| `benchmark/<method>/` | Each method's harness: data, setup (`.env`), sweeps, outputs | [`kraken2`](benchmark/kraken2/README.md), [`kmc`](benchmark/kmc/README.md), [`sshash`](benchmark/sshash/README.md), [`minimap2`](benchmark/minimap2/README.md) |
| `UPSTREAMS.md` | How each fork tracks and re-pulls its upstream | [`UPSTREAMS.md`](UPSTREAMS.md) |

## Getting started

```bash
make -C kraken2 all && make -C kraken2 test
make -C kmc all
make -C sshash all && make -C sshash test
make -C minimap2 all && make -C minimap2 test
make -C common/minimizer-orders all test
```

Then see each method's README for its own options, and each harness README for
its data and how to run it.

## Where to read about

- **The orders, their names and each method's flags:** [`orders/README.md`](orders/README.md).
- **How the orders and the mod-minimizer are implemented, and how each method breaks ties:** [`common/minimizer-orders/README.md`](common/minimizer-orders/README.md).
- **The measurement protocol, the full benchmark run and the workbook columns:** [`benchmark/README.md`](benchmark/README.md).
- **The analyses:** also in [`benchmark/README.md`](benchmark/README.md): significance tests of order differences, how closely each metric follows the density, and the run-to-run stability of two result sets.
- **Running one method's benchmark, its datasets and its outputs:** `benchmark/<method>/README.md`.
