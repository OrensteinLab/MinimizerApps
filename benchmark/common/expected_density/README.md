# Expected density tools

The expected density of a scheme is the limit of its particular density on a
uniformly random string. One tool per method estimates it on random DNA, by the
method's own selection code, and counts what the method counts for its
particular density, so the two are on the same scale and the tie-breaking rule
is the method's by construction:

| Tool | Method code it runs | Selected k-mers | k-mers |
|---|---|---|---|
| `kraken2_expected_density` | `MinimizerScanner` (`kraken2/src/mmscanner.cc`) | steps with a new minimizer (as `benchmark_minimizers`) | non-ambiguous steps |
| `minimap2_expected_density` | `mm_sketch` / `mm_sketch_scheme` (`sketch.c`, `scheme.cpp`) | entries of the index | bases after homopolymer compression (`--hpc`), minus k-1 |
| `kmc3_expected_density` | `CSplitter::ProcessReadsT` with KMC's `CMmer` and signature orders (`kmc_api`) | superwindows (`PutExtendedKmer` calls) | k-mers of the superwindows |
| `sshash_expected_density` | `compute_minimizer_tuples` with SSHash's `kmer_iterator` and `scheme_minimizer_iterator` | super-k-mers (a new minimizer position) | k-mers |

`driver.h` holds what the tools share, including the shared options:

- `--w W` (required): the number of k-mers per window, L - k + 1; it is the
  w + 1 of `density_factor = n_selected / n_kmers * (w + 1)`. For KMC3 it is
  KMC's k - p + 1 (`--w 20` for k = 28, p = 9 below).
- `--n-bases N` (default 10^9) of random string, cut into independent sequences
  of `--chunk C` bases (default 10^7, which biases the density by less than
  10^-6 relative).
- `--seeds s1,s2,...`: a string of its own per seed; each seed also configures
  the order.
- `--threads T`, and `--fasta FILE` to count the sequences of FILE instead.

The method-specific options (k, scheme, GreedyMini table, ...) are listed at the
top of each `*_backend.c*` file. The output is one CSV line per seed:

```
seed,n_bases,chunk_bases,n_kmers,n_selected,density_factor
```

On random FASTA/FASTQ files, every tool gave exactly the counts of its method's
own binary (`benchmark_minimizers -d`, the occurrences in a `minimap2 -d` index,
`kmc`'s `#Total_super-k-mers`, `sshash build --verbose`'s
`num_minimizer_positions`) for every scheme, both polarities of Φ and several
seeds. To re-check, run a tool with `--fasta FILE` and the method's binary on
the same FILE with the same scheme and seed, and compare `n_selected` with the
binary's count.

The methods' adapters run the tools through
`benchmark/common/expected_density_driver.py` (`tool_density_factor`), which
builds them on first use: `benchmark/kraken2/expected_density.py`,
`benchmark/minimap2/lib/expected_density.py`,
`benchmark/kmc/scripts/expected_density.py` and
`benchmark/sshash/scripts/expected_density.py`.

From `benchmark/common/expected_density/`:

```bash
make                                    # all four, into build/
make SSH_EXT=/path/to/sshash/external   # when sshash/external/pthash is not checked out here
build/kmc3_expected_density --k 28 --p 9 --w 20 --scheme default --n-bases 1000000000
```
