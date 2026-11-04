# Minimap2 (fork)

Fork of [lh3/minimap2](https://github.com/lh3/minimap2) (long-read
mapper/aligner), pinned to release `v2.31`. See
[`UPSTREAMS.md`](../UPSTREAMS.md) for how this folder tracks and re-pulls
from the real upstream project.

Minimap2 is C (the other three forks, `kraken2/`, `kmc/` and `sshash/`, are
C++), so the schemes are a small, self-contained addition:

| File | Role |
|---|---|
| `scheme.h` | Plain C-callable declarations (extern-"C" guarded, like `minimap.h`) |
| `scheme.cpp` | The only C++ translation unit: wraps `common/minimizer-schemes/` **unmodified** behind that flat API |
| `sketch_gm.c` | The `greedymini` sketch, `mm_sketch_gm()`: `mm_sketch()`'s loop in plain C with the GreedyMini order compared rank first (see below) |
| `sketch.c`, `index.c`, `map.c`, `minimap.h`, `main.c` | Small, targeted edits -- see below |
| `Makefile`, `setup.py` | Build `scheme.cpp` with CXX and link `-lstdc++`, and `sketch_gm.c` (Makefile also: `sketch_dump`, `make test`); `setup.py` builds both into `mappy` |
| `sketch_dump.c`, `test_scheme.sh` | Test driver (one sketch to stdout; `-x` prints every minimizer's x and y, `-H` HPC) and the `make test` suite |
| `.gitignore` | Ignores `sketch_dump` |

Everything else in this folder is unmodified upstream (`git diff --stat
f8d2f59e -- minimap2`, the import commit, lists the changes).

## Custom minimizer scheme

Same schemes as the other three methods (see
[`orders/README.md`](../orders/README.md) for the shared implementation and
the CLI flags of every method): `xor`, `xor-pcg-hash`, `greedymini`
(binary alphabet with on-the-fly projection only, no DNA-alphabet table;
a `gm_k` table is center-extended to `-k = gm_k + 2i`, as in Kraken2 and
SSHash; either polarity of Φ, `--gm-projection`),
`mod-minimizer-xor` and `mod-minimizer-greedymini` (canonical
mod-sampling with `xor`, respectively a GreedyMini table of t-mers, as the
t-mer order; see below).

```bash
./minimap2 -a ref.fa reads.fq                                                    # native, unchanged
./minimap2 -a --scheme xor --order-seed 42 ref.fa reads.fq
./minimap2 -a --scheme xor-pcg-hash --order-seed 42 ref.fa reads.fq
./minimap2 -a --scheme greedymini \
    --greedymini-file ../orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq   # gm_k=15 table: -k 15, 17, 19, 21 or 23
./minimap2 -a --scheme greedymini --gm-projection at1-cg0 \
    --greedymini-file ../orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq   # A,T -> 1 / C,G -> 0 instead of the default at0-cg1

# canonical mod-minimizer, xor as O_t
./minimap2 -a --scheme mod-minimizer-xor --order-seed 42 --mod-t 5 ref.fa reads.fq
# canonical mod-minimizer, a GreedyMini table of t-mers as O_t (t = gm_k + 2i;
# the w20_k9 table at -k 19 -w 10, t = 9)
./minimap2 -a --scheme mod-minimizer-greedymini --order-seed 42 --mod-t 9 -k 19 -w 10 \
    --greedymini-file ../orders/GreedyMini/w20_k9_canonical.gm ref.fa reads.fq
```

All scheme flags are long (`--scheme`, `--order-seed`, `--greedymini-file`,
`--gm-projection`, `--mod-t`): `-O` (gap-open penalty) and `-s` (minimal DP
score) are taken by Minimap2.

The chosen scheme is embedded in a dumped index (`-d`, including split `-I`
multi-part indexes): a `.mmi` reload never needs `--greedymini-file`
re-supplied, and an index built by upstream minimap2 (or this fork with the
default scheme) still loads -- the format is versioned by an explicit
magic-byte bump (`MMI\3` vs. upstream's `MMI\2`), not inferred from a short
read. See `scheme.cpp`'s top-of-file comment for the full design rationale,
in particular why the order only ever decides *which* k-mer wins a window
and never becomes the value actually stored in the index.

An index whose stored scheme cannot be restored (an unknown scheme kind or
mod-minimizer tie rule, a truncated or corrupt scheme section, a table that does not match the
index's `-k`) is refused rather than mapped with the native scheme. Scheme
flags given together with a prebuilt index are ignored, as `-k`/`-w`/`-H` are:
the index's own scheme is used, with a warning when the flags describe a
different scheme. Custom schemes take `-k` up to 28, as native minimap2 does.

### Design notes / known limitations

- The native scheme's own sketching path (`sketch.c`'s `mm_sketch()`) is
  byte-for-byte untouched; only the new `mm_sketch_scheme()` (`scheme.cpp`)
  is added, dispatched on `mi->scheme != NULL`. Zero added cost for the
  default case. `mm_sketch_scheme()` runs `mm_sketch_gm()` (`sketch_gm.c`)
  for `greedymini`, `sketch_mod_minimizer<Kind>` for the mod-minimizer
  schemes and `sketch_single_window<Kind>` (a port of `mm_sketch()`'s loop)
  for `xor` and `xor-pcg-hash`, in index building and in mapping alike.
- `mm_idx_str()` (the in-memory index-from-strings API) and `mappy`'s
  index-from-FASTA path do not take a custom scheme -- only the CLI path
  through `mm_idx_reader_open`/`mm_idx_gen` does. Indexes built this way are
  always native-scheme; `mappy` does map against a custom-scheme `.mmi`.
- `greedymini` uses the shared center extension
  (`common/minimizer-schemes/greedymini_center.h`, injective packing): with
  `-k = gm_k + 2i` the table ranks the middle `gm_k` bases' AT/CG projection,
  then the k-mer's Sigma and Gamma planes break ties (min over both strands),
  a strand-invariant total order on canonical k-mers; `--order-seed` XORs
  its low k bits into Gamma (as kraken2's and sshash's greedymini seeds).
  `-k - gm_k` must be
  even and `-k <= 23` (the packed value must stay below UINT64_MAX).
  `mm_sketch_gm()` (`sketch_gm.c`) is `mm_sketch()`'s loop (HPC, palindrome
  skip, N reset, duplicate flush, `>=` tie rule) in plain C, comparing the
  order rank first, as Kraken2's scanner does: per base it rolls one Sigma
  bit of the forward strand and reads the k-mer's rank (the center's table
  entry, 0xFFFF outside the UHS) directly from the `uint16_t` table; only
  when two candidates' ranks are equal does it compute their tie-break
  fields (min over both strands of Sigma and Gamma XOR the toggle) from the
  stored canonical k-mers, which gives the same comparisons as the full
  value. The stored `x` is `mm_hash64()` of each selected k-mer, computed
  when it is output.
  `benchmark/minimap2/tests/test_greedymini_sketch.py` checks the compiled rank and sketch against two
  independent Python versions; `test_scheme.sh` pins the x and y of both
  greedymini schemes' sketches (checksums from the reference implementation
  that computes the full value per k-mer).
- `mod-minimizer-xor` / `mod-minimizer-greedymini` are a canonical variant of mod-sampling (Groot
  Koerkamp & Pibiri, WABI 2024), with the rule shared with kraken2 and sshash
  in `common/minimizer-schemes/mod_minimizer_rule.h`: `--mod-t` must satisfy
  t == k (mod w) (t = 9 at k = 19, t = 5 at k = 15, w = 10); a unique minimal
  t-mer at `x` samples offset `x mod w`, and a tie takes the offset closest to
  the window's center `(w-1)/2`, between two equally close ones the left if the
  window of L = w+k-1 bases is lexicographically <= its reverse complement (the
  center rule of Cologni & Pibiri, as in SSHash). A read and its reverse
  complement thus select the same loci, and the selection never moves back.
  The order on t-mers is a template parameter of the one scan
  (`sketch_mod_minimizer<Kind>`, `TmerRank`): xor ranks the canonical t-mer
  XOR the seed; the GreedyMini order is `TmerRank<MM_SCHEME_MOD_GREEDYMINI>`,
  the injective center-extension value at length t (rank << 2t | min over
  both strands of Sigma << t | Gamma XOR the toggle, with the same projection
  and `--order-seed` Gamma toggle as `greedymini`) from the two strands'
  words rolled one base per step and one table load, its table and length
  stored in the `.mmi`.
  The scan keeps the window minimum and its leftmost/rightmost t-mers by
  conditional moves and rescans only when the minimum leaves the window (as
  native `mm_sketch()` does), and emits a selected k-mer when it leaves the
  window (output stays in position order); `benchmark/minimap2/tests/test_mod_minimizer_sketch.py`
  checks it, with both orders on t-mers, against a brute force of the definition. A window whose sampled
  k-mer is its own reverse complement (even `-k` only), or spans 256 or more
  bases under `-H`, emits no minimizer; native minimap2 instead skips such a
  k-mer and lets the window pick another.

## Build and test

```bash
make all     # CC compiles the C sources; CXX (g++) compiles scheme.cpp; -lstdc++ added to LIBS
make test    # builds the sketch_dump test driver, then ./test_scheme.sh -- sanity checks for every scheme,
             # including test_mod_minimizer_sketch.py (brute force vs. the compiled mod-minimizer sketch)
             # and test_greedymini_sketch.py (both in benchmark/minimap2/tests/)
```
