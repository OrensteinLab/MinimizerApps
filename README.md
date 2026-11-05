# MinimizerApps

High-throughput sequencing analysis methods such as Kraken2 (taxonomic
classification), Minimap2 (read mapping), KMC3 (*k*-mer counting) and SSHash
(*k*-mer indexing) sample *k*-mers from DNA sequences as an important step of input
processing. All of them select the sampled *k*-mers with minimizer schemes, which
provide a window guarantee: each window of a given length *L* contains a *k*-mer
from the sample. Almost all of them use pseudo-random minimizers, although
recent tools such as GreedyMini produce minimizers that create smaller samples
with the same window guarantee, which may shorten runtimes and reduce space
usage. Moreover, selection schemes are defined on single-stranded DNA, while
these methods are applied to double-stranded DNA and thus require selection
schemes that are canonical, i.e., make the same selection on both strands, up
to reverse complements.

We bridge both gaps. We modified GreedyMini to produce canonical minimizers,
designed a canonical version of another promising scheme, the mod-minimizer,
integrated the new canonical schemes into the four methods, and evaluated the
methods end to end; this repository holds the code and the benchmarks. The new
schemes create samples up to 32% smaller than those of pseudo-random minimizers,
which led to reductions of up to 36% in runtime on certain datasets and to
reductions in space usage in proportion to the sample size. In Kraken2 and
Minimap2, whose outputs depend on the sampled *k*-mers, performance accuracy was
negligibly affected, and KMC3 and SSHash produced identical outputs.

Each method is evaluated with the schemes below, under the names used in its
CLI flags, run folders and result tables:

| Scheme | Name | Kraken2 | Minimap2 | KMC3 | SSHash |
|---|---|:-:|:-:|:-:|:-:|
| The method's default | its own | `xor` | `native` | `default` | `native` |
| Pseudo-random minimizers | `xor`, `xor-pcg-hash` | ✓ | ✓ | ✓ | ✓ |
| Canonical GreedyMini minimizer | `greedymini` | ✓ | ✓ | ✓ | ✓ |
| Canonical mod-minimizer, `xor` order on *t*-mers | `mod-minimizer-xor` | ✓ | ✓ | | ✓ |
| Canonical mod-minimizer, GreedyMini order on *t*-mers | `mod-minimizer-greedymini` | ✓ | *k* = 19 | | ✓ |

## Layout

Each method is a fork whose upstream-tracked source stays apart from the shared
scheme code and the benchmark scripts.

| Path | What it is | Docs |
|---|---|---|
| `kraken2/`, `minimap2/`, `kmc/`, `sshash/` | The four forks, each with only the changes needed to plug in the schemes | [`kraken2`](kraken2/README.md), [`minimap2`](minimap2/README.fork.md), [`kmc`](kmc/README.md), [`sshash`](sshash/README.md) |
| `orders/` | The GreedyMini order tables (`.txt` rank lists; the `.gm` tables are built from them), the scheme names and options, and the CLI flags of each method | [`orders/README.md`](orders/README.md) |
| `common/minimizer-schemes/` | The scheme code shared by all four methods: the GreedyMini tables, the PCG scramble and the mod-minimizer rule with its tie handling | [`common/minimizer-schemes/README.md`](common/minimizer-schemes/README.md) |
| `benchmark/` | The shared measurement protocol, the full run of all four methods (`full_run_all.sh`), the exporters and the analyses | [`benchmark/README.md`](benchmark/README.md) |
| `benchmark/<method>/` | Each method's benchmark scripts: data, setup (`.env`), experiments, outputs | [`kraken2`](benchmark/kraken2/README.md), [`minimap2`](benchmark/minimap2/README.md), [`kmc`](benchmark/kmc/README.md), [`sshash`](benchmark/sshash/README.md) |
| `paper-export/` | The scripts that compute the paper's results tables, statistics and figures from the result sets, and draw the figures of its methods | [`paper-export/README.md`](paper-export/README.md) |
| `supplementary_tables/` | The paper's supplementary tables, as CSV files | [`supplementary_tables/README.md`](supplementary_tables/README.md) |
| `UPSTREAMS.md` | How each fork tracks and re-pulls its upstream | [`UPSTREAMS.md`](UPSTREAMS.md) |

## Requirements

- **Build:** Linux, `g++` with C++17 and OpenMP, `make`, CMake >= 3.5 (SSHash)
  and the zlib headers.
- **Benchmarks:** Bash, Python >= 3.9 and GNU `time`. The Kraken2 and Minimap2
  benchmarks set up their own tools (`benchmark/<method>/setup.sh`); each
  method's benchmark README lists what it needs.

## Reproducing the results

1. **Build** the shared scheme code, the GreedyMini tables (only their `.txt`
   rank lists are tracked; the methods read the `.gm` tables built from them)
   and the four methods:

   ```bash
   git submodule update --init --recursive    # KMC3's zlib and SSHash's PTHash
   make -C common/minimizer-schemes all test
   for f in orders/GreedyMini/*.txt; do
     b=$(basename "$f" .txt); k=${b##*_k}; k=${k%%_*}
     common/minimizer-schemes/gm_serializer --canonical -k "$k" -o "${f%.txt}.gm" "$f"
   done
   make -C kraken2 all && make -C kraken2 test
   make -C minimap2 all && make -C minimap2 test
   make -C kmc all
   make -C sshash all && make -C sshash test
   ```

2. **Set up each method's benchmark**: copy `benchmark/<method>/.env.example`
   to `.env` and set the data and results directory, outside the repository
   (`benchmark/<method>/README.md`).
3. **Download or generate the datasets** with the commands in [Datasets](#datasets).
4. **Run the paper's experiments**: `benchmark/<method>/experiments/full_run.sh`
   for one method, or `benchmark/full_run_all.sh` for all four. The
   measurement protocol and the outputs are described in
   [`benchmark/README.md`](benchmark/README.md).
5. **Regenerate the paper's tables and figures** from the result sets:
   [`paper-export/README.md`](paper-export/README.md).

## Datasets

The benchmarks use the datasets below. Each method's benchmark README provides
their sources, their checksums where available, and the commands that
download them or, for simulated reads, generate them. Disk space is the size of
each dataset after download and preprocessing.

| Method | Dataset | Disk space (GB) | Download or generation |
|---|---|---|---|
| Kraken2 | Strain-exclusion library and Mason reads (Zenodo) | 43 | `prepare_strex_data.sh` (library and taxonomy), then `download_data.sh` (reads); see [`benchmark/kraken2`](benchmark/kraken2/README.md#quick-start) |
|  | FDA-ARGOS, 25 runs (SRA) | 0.2 | `download_data.sh`; see [`benchmark/kraken2`](benchmark/kraken2/README.md#quick-start) |
|  | CAMI II marine, 250,000 read pairs | 0.2 | `prepare_cami_data.sh` (~6 GB download, deleted after sampling); see [`benchmark/kraken2`](benchmark/kraken2/README.md#quick-start) |
| Minimap2 | CHM13 chrX v0.7, PBSIM CLR reads | 1.8 | `curl` for the genome, then `lib/pbsim_reads.py`; see [`benchmark/minimap2`](benchmark/minimap2/README.md#datasets) |
|  | GRCh38 (UCSC hg38), PBSIM CLR reads | 5.8 | `curl` for the genome, then `lib/pbsim_reads.py`; see [`benchmark/minimap2`](benchmark/minimap2/README.md#datasets) |
|  | CHM13 rel3 ONT, 1,000 reads (mapped to GRCh38) | 0.1 | `curl` streams the 355 GB FASTQ (not stored) into `lib/sample_reads.py`; see [`benchmark/minimap2`](benchmark/minimap2/README.md#datasets) |
| KMC3 | F. vesca | 3.4 | `benchmark/kmc/kmc3_bench.sh download -d f_vesca` |
|  | G. gallus | 24 | `benchmark/kmc/kmc3_bench.sh download -d g_gallus` |
|  | M. balbisiana | 51 | `benchmark/kmc/kmc3_bench.sh download -d m_balbisiana` |
|  | H. sapiens 2 | 106 | `benchmark/kmc/kmc3_bench.sh download -d h_sapiens_2` |
|  | H. sapiens 3 | 614 | `benchmark/kmc/kmc3_bench.sh download -d h_sapiens_3` |
| SSHash | cod | 0.2 | `benchmark/sshash/sshash_bench.sh download -d cod` |
|  | kestrel | 0.8 | `benchmark/sshash/sshash_bench.sh download -d kestrel` |
|  | human | 2.6 | `benchmark/sshash/sshash_bench.sh download -d human` |
|  | NCBI viruses | 0.1 | `benchmark/sshash/sshash_bench.sh download -d ncbi-virus` |
|  | S. enterica | 1.0 | `benchmark/sshash/sshash_bench.sh download -d se` |
|  | HPRC | 1.9 | `benchmark/sshash/sshash_bench.sh download -d hprc` (reuses human's query reads) |

## License

The four methods keep their upstream licenses; everything else is under the MIT
License.

| Path | License | File |
|---|---|---|
| `kraken2/` | MIT | [`kraken2/LICENSE`](kraken2/LICENSE) |
| `minimap2/` | MIT | [`minimap2/LICENSE.txt`](minimap2/LICENSE.txt) |
| `kmc/` | GNU GPL v3 (changes listed in [`kmc/MODIFICATIONS.md`](kmc/MODIFICATIONS.md)) | [`kmc/LICENSE`](kmc/LICENSE) |
| `sshash/` | MIT | [`sshash/LICENSE`](sshash/LICENSE) |
| Everything else: `common/`, `benchmark/`, `paper-export/`, and the order tables in `orders/` (generated by [GreedyMini](https://github.com/OrensteinLab/GreedyMini), MIT) | MIT | [`LICENSE`](LICENSE) |
