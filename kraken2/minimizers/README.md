# GreedyMini order in Kraken2

How Kraken2 ranks *L*-mers with a GreedyMini table (`--order greedymini`), and
the density it reaches. Notation as in the paper: *k* = length of the selected
*k*-mer, *w* = number of *k*-mers in a window, *L* = *k* + *w* − 1, *t* = *t*-mer length of
the mod-minimizer. Kraken2's own names differ: its minimizer length ℓ (`-l`)
is *k*, and its *k*-mer length (`-k`) is *L*, so `-l 31 -k 35` is *k* = 31, *w* = 5.
Paths below are from the repo root.

GreedyMini tables are binary (σ = 2). Kraken2 packs bases as ACGT = `00, 01,
10, 11` and splits them into two bit planes: Σ, the A,T vs C,G bit (the
projection Φ; `--gm-projection at1-cg0` flips its polarity, default
`at0-cg1`), and Γ, the residual bit.

**Center extension (*k* > gm_k).** A canonical (mirror-symmetric) table is
extended to *k* = gm_k + 2i by viewing the *k*-mer as `x y z`, `|x| = |z| = i`,
`|y| = gm_k`: the rank is `table[Σ(y)]` (one lookup; strand invariant because
the reverse complement has the mirror of Σ(y) at its center); ties are broken
by the strand with the smaller Σ(xyz), then by the DNA (Γ), so the order is
total on canonical *k*-mers. It keeps the table's designed density at every
valid *k* (see [Density](#density)). The algorithm is the shared
`common/minimizer-schemes/greedymini_center.h` (the same extension SSHash and
Minimap2 use, including the plane extraction `TwoBitPlanes`);
`kraken2/src/minimizer_rankers.h` only supplies Kraken2's ACGT encoding. For
*k* = gm_k (and for the mod-minimizer with a GreedyMini order, *t* = gm_k) the
table ranks the *k*-mer itself and Γ breaks ties. GreedyMini is DNA only
(protein databases refuse it).

**Scanning (rank first).** The scanner compares each *k*-mer's table rank (one
Σ `pext` and one table load, on the unmasked Σ, also with a spaced seed) and
computes the full packed value only when a finite rank ties with the back of
the sliding-window queue (about 0.3% of the comparisons on E. coli unitigs). It
selects exactly the *k*-mers of the full-`Rank()` path (`make -C kraken2 test`
checks both, plus a rolled Σ plane, `MinimizerScanner::SetRankFirst`;
`benchmark_minimizers --rank-path full|first|roll` times them).

- **Valid *k*:** `k − gm_k` even and *k* ≤ 31, e.g. `w5_k15_canonical`:
  *k* = 15, 17, …, 31 (default 31). Other *k* are rejected by
  `kraken2-build` / `build_db` / the scanner.
- **Table:** binary alphabet, mirror-symmetric, each rank shared only by a key
  and its mirror (checked at load; all `orders/GreedyMini/*_canonical.gm`
  qualify).
- **Spaced seed (`--minimizer-spaces`, default 7):** the rank lookup always
  sees the unmasked center; the mask (ANDed into each strand's word) only
  coarsens the tie-break fields. Masking y would cost the table's density
  (1.77 instead of 1.47 at *k* = 31, *s* = 7). `s ≤ k/4`, so *s* = 7 needs *k* ≥ 28.
- **Databases:** `opts.k2d` stores the scheme (`db_version` 2, 80 bytes, magic
  `K2SK`: `uint8 sketching_scheme; uint8 order; uint16 gm_k; uint32 magic;
  uint64 gm_num_entries` after the original fields). `classify` and
  `dump_table` refuse any other database, including upstream ones, with a
  "rebuild" error, and `k2 build` does not count one as already built.

## Tables

The tables live in [`orders/GreedyMini/`](../../orders/GreedyMini/), shared
with the other methods (see [`orders/README.md`](../../orders/README.md)).
Only the `.txt` rank lists are tracked; the `.gm` tables are built from them
(see the root README, or `bench_ensure_gm_file` in
`benchmark/kraken2/benchmark_interface.sh`). Kraken2 supports only
binary-alphabet tables. The benchmark uses:

| File | Role |
|------|------|
| `w5_k15_canonical.gm` | `greedymini` (`-g` / `GM_FILE`). Density factor 1.466 at every valid *k* = 15, 17, …, 31 |
| `w5_k6_canonical.gm` | `mod-minimizer-greedymini` (`-M` / `MOD_GM_FILE`): gm_k = *t* = 6 for *k* = 31, *w* = 5 |

`.gm` format (`kraken2/src/greedymini_rank.h`):

```
[uint64 length N][N × uint64 ranks]
```

- `N` = `2^gm_k` (binary alphabet)
- UHS members store a finite rank `0 .. m-1`
- other keys store the sentinel `N`

## Ranked list → `.gm`

[`common/minimizer-schemes/gm_serializer.cpp`](../../common/minimizer-schemes/gm_serializer.cpp)
turns a plain-text list of integer *k*-mers (in rank order) into a `.gm` table:

- **Default:** each listed *k*-mer alone gets that rank.
- **`--canonical`:** each entry is the smaller of `{x, reverse(x)}`; both
  share the rank (palindromes once).

```bash
make -C common/minimizer-schemes

common/minimizer-schemes/gm_serializer --canonical \
    orders/GreedyMini/w5_k15_canonical.txt \
    -o orders/GreedyMini/w5_k15_canonical.gm

common/minimizer-schemes/gm_serializer --canonical \
    orders/GreedyMini/w5_k6_canonical.txt \
    -o orders/GreedyMini/w5_k6_canonical.gm -k 6

# plain ranked list (no reverse pairing)
common/minimizer-schemes/gm_serializer ranked_kmers.txt -o ranked_kmers.gm -k 15
```

Options: `--canonical` (same rank for a *k*-mer and its binary reverse), `-k INT`
(binary length; default inferred from the largest *k*-mer), `-o PATH` (default:
same basename with `.gm`). The DNA projection is not stored in the file;
Kraken2 applies it when ranking.

## Density

Density factor = (selected *k*-mers / windows) × (*w* + 1), *w* = 5; a random order
gives ≈ 2, `w5_k15_canonical` is designed for 1.466. "Selected" is what
`benchmark_minimizers` counts: successive windows whose selected *k*-mer
changes. Numbers are mean ± std (sample) over 6 independent i.i.d. random DNA
streams of 3 Mb (stream i: `--seed i` and its own toggle mask), spaced seed off
(*s* = 0) unless stated. At *k* = 15 = gm_k the order is the table itself (the
design density, 1.4662 ± 0.0003 measured). `xor` is Kraken2's default order
(xor with the toggle mask; its density varies with the toggle), `xor-pcg-hash`
a proper random hash.

**Expected density, `w5_k15_canonical`**
([`benchmark/kraken2/tables/gm_center_density_expected.csv`](../../benchmark/kraken2/tables/gm_center_density_expected.csv))

| *k* | greedymini (center) | xor | xor-pcg-hash |
|---|---|---|---|
| 15 (= gm_k, design) | 1.4662 ± 0.0003 | 2.019 ± 0.064 | 1.9993 ± 0.0011 |
| 17 | 1.4659 ± 0.0003 | 2.034 ± 0.058 | 1.9997 ± 0.0005 |
| 19 | 1.4667 ± 0.0003 | 2.012 ± 0.067 | 2.0006 ± 0.0009 |
| 21 | 1.4674 ± 0.0003 | 2.047 ± 0.048 | 1.9997 ± 0.0008 |
| 23 | 1.4668 ± 0.0003 | 2.027 ± 0.062 | 1.9997 ± 0.0009 |
| 25 | 1.4661 ± 0.0004 | 2.001 ± 0.075 | 1.9997 ± 0.0008 |
| 27 | 1.4661 ± 0.0003 | 2.022 ± 0.053 | 1.9996 ± 0.0005 |
| 29 | 1.4662 ± 0.0003 | 2.034 ± 0.046 | 1.9999 ± 0.0011 |
| **31 (default)** | **1.4661 ± 0.0003** | 2.024 ± 0.074 | 2.0001 ± 0.0005 |

The center extension stays at the table's value for every valid *k*. With the
spaced seed on (*s* = 7, *k* ≥ 28) it is unchanged to the 4th decimal (*k* = 29:
1.4662; *k* = 31: 1.4661).

Other tables at their own design window: `w11_k13` (*w* = 11): *k* = 13
1.5913 ± 0.0004, *k* = 21 1.5895 ± 0.0005, *k* = 31 1.5914 ± 0.0005. `w11_k15`:
*k* = 15 1.5236 ± 0.0002, *k* = 21 1.5248, *k* = 31 1.5231. `w20_k9` (*w* = 20): *k* = 9
1.7180 ± 0.0007, *k* = 21 1.7200, *k* = 31 1.7179. `w5_k6` is the mod-minimizer
table (*k* = gm_k = *t*).

**Particular density** (real sequence): the library FASTAs of the benchmark's
`standard_p0.5` library (`benchmark_minimizers.sh -z`; archaea, bacteria,
human, UniVec_Core, viral; ≈ 1.42 Gbp, scanned in full, one pass per scheme),
`w5_k15_canonical`, one run per cell
([`benchmark/kraken2/tables/gm_center_density_particular.csv`](../../benchmark/kraken2/tables/gm_center_density_particular.csv)):

| *k* | greedymini (center) | xor | xor-pcg-hash |
|---|---|---|---|
| 15 (= gm_k) | 1.4730 | 1.954 | 2.000 |
| 17 | 1.4727 | 2.099 | 1.9998 |
| 21 | 1.4749 | 2.051 | 2.000 |
| 25 | 1.4726 | 2.020 | 1.9998 |
| 29 | 1.4731 | 2.055 | 1.9999 |
| **31 (default)** | **1.4731** | 2.100 | 2.000 |

The same holds with the default spaced seed (*s* = 7) at *k* = 29 / 31: 1.4731 /
1.4731. Real data sits ≈ 0.007 above the random-DNA value (0.4–0.5 %).
`NextMinimizer` costs ≈ 93–95 ns per selected *k*-mer (xor ≈ 91 ns, on a loaded
machine).

**Reproduce** (from the repo root; `nice -n 19`, ≤ 4 threads):

```bash
make -C kraken2/src benchmark_minimizers
S=benchmark/kraken2/density_sweep.sh; G=orders/GreedyMini/w5_k15_canonical.gm
# expected (6 streams x 3 Mb; -s = minimizer spaces; add -L "29 31" with -s 7; -R adds the random baselines)
$S -b kraken2/src/benchmark_minimizers -g $G -v center -o e.csv -s 0 -n 3000000 -S 6 -R
python3 benchmark/kraken2/density_summary.py e.csv     # mean +- std per k
# particular (one stream is enough; -d = the library)
DB=$KRAKEN2_BENCH_ROOT/db/standard_p0.5
for l in 31 25 21 17; do
  $S -b kraken2/src/benchmark_minimizers -g $G -v center -o p.csv -L $l -S 1 -n 100000 -d $DB -j 3 -R
done
python3 benchmark/kraken2/density_summary.py p.csv --col fac_particular
```

`-v` is a free label written to the CSV's `variant` column, `-L` takes
Kraken2's ℓ (= *k*), and `-P` the polarities of Φ (default both). Per (*k*,
polarity, stream i), `density_sweep.sh` calls `nice -n 19 BIN -S greedymini
[xor xor-pcg-hash] -g $G --gm-projection P -l k -k k+4 -n 3000000 -j 1
-s SPACES -T TOGGLE_i --seed i --seed-idx i -d DB -o out.csv`. Unit tests:
`make -C kraken2 test` (per-base reference over ≥ 1e7 random *k*-mers, strand
invariance, injectivity, rejections, database guard, equality with SSHash's
Python vectors).
