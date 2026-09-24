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
the cross-app CLI-flag table): `xor-hash`, `xor-pcg-hash`, `greedymini`
(binary alphabet + on-the-fly projection only -- no DNA-alphabet backend, no
Extend-k: the table's `gm_k` must equal `-k` exactly), and
`mod-minimizer-{xor-hash,pcg-xor-hash,greedymini}`.

```bash
./minimap2 -a ref.fa reads.fq                                                    # native, unchanged
./minimap2 -a --minimizer-order xor-hash --order-seed 42 ref.fa reads.fq
./minimap2 -a --minimizer-order xor-pcg-hash --order-seed 42 ref.fa reads.fq
./minimap2 -a --minimizer-order greedymini \
    --greedymini-file ../orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq   # needs -k 15 -w 10 (or matching gm_k)

# mod-minimizer (Groot Koerkamp & Pibiri, WABI 2024): any base order as O_t
./minimap2 -a --minimizer-order mod-minimizer-xor-hash --order-seed 42 --mod-t 5 ref.fa reads.fq
./minimap2 -a --minimizer-order mod-minimizer-greedymini \
    --greedymini-file ../orders/GreedyMini/w10_k5_canonical.gm --mod-t 5 ref.fa reads.fq  # table's gm_k must equal --mod-t
```

Long-only flags (`--minimizer-order`, `--order-seed`, `--greedymini-file`,
`--gm-alphabet`, `--mod-t`): unlike kraken2/kmc/sshash's short flags, `-O`
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
- `greedymini`/`mod-minimizer-greedymini` require a table whose `gm_k`
  matches exactly (`-k` for the former, `--mod-t` for the latter) --
  Extend-k (a table built for a shorter length than it's used at) isn't
  supported here, unlike kraken2/kmc/sshash.

## Build and test

```bash
make all     # CC compiles the C sources; CXX (g++) compiles order.cpp; -lstdc++ added to LIBS
make test    # ./test_order.sh -- sanity checks for every order except greedymini (needs a table; see above)
```
