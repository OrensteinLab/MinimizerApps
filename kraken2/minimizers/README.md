# Minimizer rank tables

GreedyMini (and related) UHS rank tables used by Kraken 2 when
`--order GreedyMini` is selected. Binary tables are σ-extended to DNA at
rank time: Kraken packs bases as ACGT=`00,01,10,11`, remaps to ATCG bitplanes
so the key bit is **AT vs CG**, and breaks ties with the residual bit
(lexicographic τ). Orientation picks the DNA strand with the better
gm_k UHS primary, then applies Extend-k / Extend-σ in that frame.

Tables now live under [`../../orders/GreedyMini/`](../../orders/GreedyMini/),
shared with kmc's own GreedyMini order (see `orders/README.md`). Only the
`.txt` rank-list sources are tracked in git; `.gm` binaries are generated
from them on demand (`../benchmark/kraken2/benchmark_interface.sh`'s
`bench_ensure_gm_file`). Defaults used by the benchmark scripts:

| File | Role |
|------|------|
| `w5_k15_canonical.gm` | Minimizer + binary GreedyMini (`-g` / `GM_FILE`; default). Extend-k + AT/CG σ-extension keeps dens factor ~1.47 for ℓ ≥ gm_k=15. |
| `w5_k6_canonical.gm` | Mod + binary GreedyMini (`-M` / `MOD_GM_FILE`; default). gm_k=t=6 with Kraken ℓ=31,k=35 (w=5). |

Kraken 2 only ever uses the binary-alphabet GreedyMini backend (its gm_k
values make a `4^gm_k` DNA table impractical); a few other tables that used
to live here (including a DNA-alphabet one that turned out to be
non-canonical) are unused and have been moved to `~/unused-orders/kraken2/`
outside the repo. See `orders/README.md` for the full picture.

Classic `.gm` encoding (see `src/greedymini_rank.h`):

```
[uint64 length N][N × uint64 ranks]
```

- `N` is `2^k` (binary) or `4^k` (DNA alphabet)
- UHS members store a finite rank `0 .. m-1`
- Non-UHS keys store sentinel `N`

## Ranked list → `.gm`

[`common/minimizer-orders/gm_serializer.cpp`](../../common/minimizer-orders/gm_serializer.cpp)
(built via that directory's `Makefile`, shared with kmc) turns a plain-text
list of integer k-mers (rank order) into a `.gm` table. It's usually built
and run automatically (`../benchmark/kraken2/benchmark_interface.sh`'s
`bench_ensure_gm_file`), but can be run directly:

- **Default:** each listed k-mer alone gets that rank.
- **`--canonical`:** each entry is the smaller of `{x, reverse(x)}`; both
  share the rank (palindromes once).

```bash
make -C ../common/minimizer-orders

# reverse-collapsed (canonical) UHS
../common/minimizer-orders/gm_serializer --canonical \
    ../orders/GreedyMini/w5_k15_canonical.txt \
    -o ../orders/GreedyMini/w5_k15_canonical.gm

../common/minimizer-orders/gm_serializer --canonical \
    ../orders/GreedyMini/w5_k6_canonical.txt \
    -o ../orders/GreedyMini/w5_k6_canonical.gm -k 6

# plain ranked list (no reverse pairing)
../common/minimizer-orders/gm_serializer ranked_kmers.txt -o ranked_kmers.gm -k 15
```

Options:

- `--canonical` — assign the same rank to each k-mer and its binary reverse
- `-k INT` — binary length (default: inferred from the largest k-mer)
- `-o PATH` — output path (default: same basename with `.gm`)

σ-extension for DNA is **not** written into the file: with `--gm-alphabet
binary`, Kraken 2 remaps ACGT→ATCG (AT/CG key, residual lex tie-break) when
ranking. See `src/greedymini_rank.h`.
