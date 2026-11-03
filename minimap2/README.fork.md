# minimap2 (fork)

Fork of [lh3/minimap2](https://github.com/lh3/minimap2) (long-read
mapper/aligner), pinned to release `v2.31`. See
[`UPSTREAMS.md`](../UPSTREAMS.md) for how this folder tracks and re-pulls
from the real upstream project.

Unlike the other three forks in this monorepo (`kraken2/`, `kmc/`,
`sshash/`, all C++), minimap2 is pure C, so the custom-minimizer-scheme
integration doesn't take the shape of scattered upstream-adjacent edits --
it's a small, self-contained addition:

| File | Role |
|---|---|
| `scheme.h` | Plain C-callable declarations (extern-"C" guarded, like `minimap.h` already is) |
| `scheme.cpp` | The only C++ translation unit here: wraps `common/minimizer-schemes/` **unmodified** behind that flat API |
| `sketch.c`, `index.c`, `map.c`, `minimap.h`, `main.c` | Small, targeted edits -- see below |

Everything else in this folder is an unmodified import of upstream.

## Custom minimizer scheme

Same scheme menu as the other three apps (see
[`orders/README.md`](../orders/README.md) for the shared implementation and
the cross-app CLI-flag table): `xor`, `xor-pcg-hash`, `greedymini`
(binary alphabet + on-the-fly projection only, no DNA-alphabet backend;
a `gm_k` table is centre-extended to `-k = gm_k + 2i`, as in kraken2 and
sshash; either AT/CG projection polarity, `--gm-projection`, as in sshash),
`mod-minimizer-xor` and `mod-minimizer-greedymini` (strand-symmetric
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

# strand-symmetric mod-minimizer, xor as O_t
./minimap2 -a --scheme mod-minimizer-xor --order-seed 42 --mod-t 5 ref.fa reads.fq
# strand-symmetric mod-minimizer, a GreedyMini table of t-mers as O_t (t = gm_k + 2i;
# the w20_k9 table at -k 19 -w 10, t = 9)
./minimap2 -a --scheme mod-minimizer-greedymini --order-seed 42 --mod-t 9 -k 19 -w 10 \
    --greedymini-file ../orders/GreedyMini/w20_k9_canonical.gm ref.fa reads.fq
```

Long-only flags (`--scheme`, `--order-seed`, `--greedymini-file`,
`--gm-projection`, `--mod-t`): unlike kraken2/kmc/sshash's short flags, `-O`
(gap-open penalty) and `-s` (min DP score) were already taken by minimap2
itself.

The chosen scheme is embedded in a dumped index (`-d`, including split `-I`
multi-part indexes): a `.mmi` reload never needs `--greedymini-file`
re-supplied, and an index built by upstream minimap2 (or this fork with the
default scheme) still loads fine -- the format is versioned by an explicit
magic-byte bump (`MMI\3` vs. upstream's `MMI\2`), not inferred from a short
read. See `scheme.cpp`'s top-of-file comment for the full design rationale,
in particular why the order only ever decides *which* k-mer wins a window
and never becomes the value actually stored in the index.

An index whose stored scheme cannot be restored (an unknown or retired scheme
kind, a truncated or corrupt scheme section, a table that does not match the
index's `-k`) is refused rather than mapped with the native scheme. Scheme
flags given together with a prebuilt index are ignored, as `-k`/`-w`/`-H` are:
the index's own scheme is used, with a warning when the flags describe a
different scheme. Custom schemes take `-k` up to 28, as native minimap2 does.

### Design notes / known limitations

- The native scheme's own sketching path (`sketch.c`'s `mm_sketch()`) is
  byte-for-byte untouched; only the new `mm_sketch_scheme()` (`scheme.cpp`)
  is added, dispatched on `mi->scheme != NULL`. Zero added cost for the
  default case.
- `mm_idx_str()` (the in-memory index-from-strings API) and `mappy`'s
  index-from-FASTA path do not take a custom scheme -- only the CLI path
  through `mm_idx_reader_open`/`mm_idx_gen` does. Indexes built this way are
  always native-scheme; `mappy` does map against a custom-scheme `.mmi`.
- `greedymini` uses the shared centre extension
  (`common/minimizer-schemes/greedymini_centre.h`, injective packing): with
  `-k = gm_k + 2i` the table ranks the middle `gm_k` bases' AT/CG projection,
  then the k-mer's Sigma and Gamma planes break ties (min over both strands),
  a strand-invariant total order on canonical k-mers; `--order-seed` XORs
  its low k bits into Gamma (as kraken2's and sshash's greedymini seeds).
  `-k - gm_k` must be
  even and `-k <= 23` (the packed value must stay below UINT64_MAX).
  `benchmark/minimap2/tests/test_greedymini_sketch.py` checks the compiled rank and sketch against two
  independent Python versions.
- `mod-minimizer-xor` / `mod-minimizer-greedymini` are a strand-symmetric variant of mod-sampling (Groot
  Koerkamp & Pibiri, WABI 2024), with the rule shared with kraken2 and sshash
  in `common/minimizer-schemes/mod_minimizer_rule.h`: `--mod-t` must satisfy
  t == k (mod w) (t = 9 at k = 19, t = 5 at k = 15, w = 10); a unique minimal
  t-mer at `x` samples offset `x mod w`, and a tie takes the offset closest to
  the window centre `(w-1)/2`, between two equally close ones the left if the
  window of `w+k-1` bases is lexicographically <= its reverse complement (the
  centre rule of Cologni & Pibiri, as in SSHash). A read and its reverse
  complement thus select the same loci, and the selection never moves back.
  The order on t-mers is a template parameter of the one scan
  (`sketch_mod_minimizer<Kind>`, `TmerRank`): xor ranks the canonical t-mer
  XOR the seed; the GreedyMini order is `RankState<MM_SCHEME_GREEDYMINI>` at
  length t (the shared centre-extension order, built for t, with the same
  projection and `--order-seed` Gamma toggle as `greedymini`), stored in the
  `.mmi` with its table and length.
  The scan keeps the window minimum and its leftmost/rightmost t-mers by
  conditional moves and rescans only when the minimum leaves the window (as
  native `mm_sketch()` does), and emits a selected k-mer when it leaves the
  window (output stays in position order); an index built with the earlier
  G/T-strand tie rule is refused (rebuild it); `benchmark/minimap2/tests/test_mod_minimizer_sketch.py`
  checks it, with both orders on t-mers, against a brute force of the definition. A window whose sampled
  k-mer is its own reverse complement (even `-k` only), or spans 256 or more
  bases under `-H`, emits no minimizer; native minimap2 instead skips such a
  k-mer and lets the window pick another. Indexes built with the earlier
  forward-only mod-minimizer (scheme kind 4) are rejected: rebuild them.

## Build and test

```bash
make all     # CC compiles the C sources; CXX (g++) compiles scheme.cpp; -lstdc++ added to LIBS
make test    # builds the sketch_dump test driver, then ./test_scheme.sh -- sanity checks for every scheme,
             # including test_mod_minimizer_sketch.py (brute force vs. the compiled mod-minimizer sketch)
             # and test_greedymini_sketch.py (both in benchmark/minimap2/tests/)
```
