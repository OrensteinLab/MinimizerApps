# minimap2 (fork)

Fork of [lh3/minimap2](https://github.com/lh3/minimap2) (long-read
mapper/aligner), pinned to release `v2.31`. See
[`UPSTREAMS.md`](../UPSTREAMS.md) for how this folder tracks and re-pulls
from the real upstream project.

Unlike the other three forks in this monorepo (`kraken2/`, `kmc/`,
`sshash/`, all C++), minimap2 is pure C, so the custom-minimizer-order
integration doesn't take the shape of scattered upstream-adjacent edits --
it's a small, self-contained addition:

| File | Role |
|---|---|
| `order.h` | Plain C-callable declarations (extern-"C" guarded, like `minimap.h` already is) |
| `order.cpp` | The only C++ translation unit here: wraps `common/minimizer-orders/` **unmodified** behind that flat API |
| `sketch.c`, `index.c`, `map.c`, `minimap.h`, `main.c`, `options.c` | Small, targeted edits -- see below |

Everything else in this folder is an unmodified import of upstream.

## Custom minimizer order

Same order menu as the other three apps (see
[`orders/README.md`](../orders/README.md) for the shared implementation and
the cross-app CLI-flag table): `xor`, `xor-pcg-hash`, `greedymini`
(binary alphabet + on-the-fly projection only, no DNA-alphabet backend;
a `gm_k` table is centre-extended to `-k = gm_k + 2i`, as in kraken2 and
sshash; either AT/CG projection polarity, `--gm-projection`, as in sshash), and
`mod-minimizer-xor` (strand-symmetric mod-sampling with `xor` as the t-mer
order; see below).

```bash
./minimap2 -a ref.fa reads.fq                                                    # native, unchanged
./minimap2 -a --minimizer-order xor --order-seed 42 ref.fa reads.fq
./minimap2 -a --minimizer-order xor-pcg-hash --order-seed 42 ref.fa reads.fq
./minimap2 -a --minimizer-order greedymini \
    --greedymini-file ../orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq   # gm_k=15 table: -k 15, 17, 19, 21 or 23
./minimap2 -a --minimizer-order greedymini --gm-projection at1-cg0 \
    --greedymini-file ../orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq   # A,T -> 1 / C,G -> 0 instead of the default at0-cg1

# strand-symmetric mod-minimizer, xor as O_t
./minimap2 -a --minimizer-order mod-minimizer-xor --order-seed 42 --mod-t 5 ref.fa reads.fq
```

Long-only flags (`--minimizer-order`, `--order-seed`, `--greedymini-file`,
`--gm-projection`, `--mod-t`): unlike kraken2/kmc/sshash's short flags, `-O`
(gap-open penalty) and `-s` (min DP score) were already taken by minimap2
itself.

The chosen order is embedded in a dumped index (`-d`, including split `-I`
multi-part indexes): a `.mmi` reload never needs `--greedymini-file`
re-supplied, and an index built by upstream minimap2 (or this fork with the
default order) still loads fine -- the format is versioned by an explicit
magic-byte bump (`MMI\3` vs. upstream's `MMI\2`), not inferred from a short
read. See `order.cpp`'s top-of-file comment for the full design rationale,
in particular why the order only ever decides *which* k-mer wins a window
and never becomes the value actually stored in the index.

### Design notes / known limitations

- The native order's own sketching path (`sketch.c`'s `mm_sketch()`) is
  byte-for-byte untouched; only the new `mm_sketch_ordered()` (`order.cpp`)
  is added, dispatched on `mi->order != NULL`. Zero added cost for the
  default case.
- `mm_idx_str()` (the in-memory index-from-strings API, e.g. used by the
  `mappy` Python bindings) does not support a custom order -- only the CLI
  path through `mm_idx_reader_open`/`mm_idx_gen` does. Indexes built this
  way are always native-order.
- `greedymini` uses the shared centre extension
  (`common/minimizer-orders/greedymini_centre.h`, injective packing): with
  `-k = gm_k + 2i` the table ranks the middle `gm_k` bases' AT/CG projection,
  then the k-mer's Sigma and Gamma planes break ties (min over both strands),
  a strand-invariant total order on canonical k-mers. `-k - gm_k` must be
  even and `-k <= 23` (the packed value must stay below UINT64_MAX).
  `validate_greedymini.py` checks the compiled rank and sketch against two
  independent Python versions.
- `mod-minimizer-xor` is a strand-symmetric variant of mod-sampling (Groot
  Koerkamp & Pibiri, WABI 2024), not a published scheme: each window of
  `w+k-1` bases prefers the strand with more G/T (tie: the lexicographically
  smaller strand); a forward-preferred window samples offset `x_L mod w`
  (`x_L` = leftmost minimal t-mer, the paper's rule), a reverse-preferred one
  `(w-1) - (((T-1) - x_R) mod w)` (`x_R` = rightmost minimal t-mer, `T =
  w+k-t`). A read and its reverse complement thus select the same loci,
  which plain mod-sampling does not guarantee. The scan uses kraken2's single
  monotone deque over all live t-mers and emits a selected k-mer when it
  leaves the window (output stays in position order); `validate_mod_minimizer.py`
  checks it against a brute force of the definition. Indexes built with the
  earlier forward-only mod-minimizer (order kind 4) are rejected: rebuild them.

## Build and test

```bash
make all     # CC compiles the C sources; CXX (g++) compiles order.cpp; -lstdc++ added to LIBS
make test    # builds the sketch_dump test driver, then ./test_order.sh -- sanity checks for every order,
             # including validate_mod_minimizer.py (brute force vs. the compiled mod-minimizer sketch)
             # and validate_greedymini.py
```
