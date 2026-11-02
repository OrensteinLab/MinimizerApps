# Expected density tools

The expected density of a scheme is the limit of its particular density on a
uniformly random string. One tool per app estimates it on random DNA, by the
app's own selection code, and counts what the app counts for its particular
density, so the two are on the same scale and the tie-breaking rule is the
app's by construction:

| Tool | App code it runs | Selected k-mers | k-mers |
|---|---|---|---|
| `kraken2_expected_density` | `MinimizerScanner` (`kraken2/src/mmscanner.cc`) | steps with a new minimizer (as `benchmark_minimizers`) | non-ambiguous steps |
| `minimap2_expected_density` | `mm_sketch` / `mm_sketch_scheme` (`sketch.c`, `scheme.cpp`) | entries of the index | bases after homopolymer compression (`--hpc`), minus k-1 |
| `kmc3_expected_density` | `CSplitter::ProcessReadsT` with KMC's `CMmer` and signature orders (`kmc_api`) | superwindows (`PutExtendedKmer` calls) | k-mers of the superwindows |
| `sshash_expected_density` | `compute_minimizer_tuples` with SSHash's `kmer_iterator` and `scheme_minimizer_iterator` | super-k-mers (a new minimizer position) | k-mers |

`driver.h` holds what the tools share: the random string (`--n-bases`, default
10^9, cut into independent sequences of `--chunk` bases, default 10^7, which biases
the density by less than 10^-6 relative), a string of its own per seed
(`--seeds`; each seed also configures the order), threads (`--threads`), and the
output, one CSV line per seed:

```
seed,n_bases,chunk_bases,n_kmers,n_selected,density_factor
```

With `--fasta FILE` a tool counts the sequences of FILE instead. On random
FASTA/FASTQ files, every tool gave exactly the counts of its app's own binary
(`benchmark_minimizers -d`, the occurrences in a `minimap2 -d` index, `kmc`'s
`#Total_super-k-mers`, `sshash build --verbose`'s `num_minimizer_positions`)
for every scheme, both polarities of the GreedyMini projection and several seeds.

The apps' adapters (`benchmark/*/expected_density.py`) run the tools through
`benchmark/common/expected_density_driver.py` (`tool_density_factor`), which
builds them on first use.

```bash
make                                    # all four, into build/
make SSH_EXT=/path/to/sshash/external   # when sshash/external/pthash is not checked out here
build/kmc3_expected_density --k 28 --p 9 --w 20 --scheme default --n-bases 1000000000
```

Options per tool are listed at the top of each `*_backend.c*` file.
