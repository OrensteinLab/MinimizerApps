# Minimizer rank tables

GreedyMini (and related) UHS rank tables used by Kraken 2 when
`--order GreedyMini` is selected. Binary tables are σ-extended to DNA at
rank time: Kraken packs bases as ACGT=`00,01,10,11`, remaps to ATCG bitplanes
so the key bit is **AT vs CG** (Σ), and the residual bit is Γ.

**Centre extension (ℓ > gm_k).** A canonical (mirror-symmetric) table is
extended to ℓ = gm_k + 2i by viewing the ℓ-mer as `x y z`, `|x| = |z| = i`,
`|y| = gm_k`: the rank is `table[Σ(y)]` (one lookup; strand invariant because
the reverse complement has the mirror of Σ(y) at its centre), ties are broken
by the strand frame with the smaller Σ(xyz), then by the DNA (Γ), so the order
is total on canonical ℓ-mers. It keeps the table's designed density at every
valid ℓ (below). It replaces the earlier *prefix* extension (first gm_k bases
of the better strand), whose density depends erratically on ℓ because the two
strands' prefixes are different windows of the ℓ-mer. The algorithm is the
shared `common/minimizer-orders/greedymini_centre.h` (the same extension
sshash uses, including the Σ/Γ plane extraction: `TwoBitPlanes`);
`src/minimizer_rankers.h` only picks Kraken's ACGT encoding. For ℓ = gm_k
(and Mod+GM, t = gm_k) nothing changed.

**Scanning (rank first).** The scanner compares each ℓ-mer's table rank (one
Σ pext + one table load, the unmasked Σ, also with a spaced seed) and computes
the full packed value only when a finite rank ties with the back of the
sliding-window queue (about 0.3% of the comparisons on E. coli unitigs); it
selects exactly the minimizers of the full-`Rank()` path (`make test` checks
both, plus a rolled Σ plane, `MinimizerScanner::SetRankFirst`;
`benchmark_minimizers --rank-path full|first|roll` times them).

- **Valid ℓ:** `ℓ - gm_k` must be even and ℓ ≤ 31, e.g. `w5_k15_canonical`:
  ℓ = 15, 17, …, 31 (default 31). Other ℓ (16, 30, …) are rejected with a
  clear error by `kraken2-build`/`build_db`/the scanner.
- **Table:** binary alphabet, mirror-symmetric, each rank shared by a key and
  its mirror only (checked at load; all shipped `*_canonical.gm` qualify).
- **Spaced seed (`--minimizer-spaces`, default 7):** the rank lookup always
  sees the *unmasked* centre; the mask (ANDed into each strand's word, never
  RC(x & mask)) only coarsens the tie-break fields. A mask applied to y would
  cost the table's density (measured 1.77 instead of 1.47 at ℓ = 31, s = 7).
  Note `s ≤ ℓ/4`, so s = 7 needs ℓ ≥ 28.
- **Databases:** the new order selects different minimizers, so
  `db_version` is 2. `classify` and `dump_table` accept only `db_version` 2
  and refuse any other database (built by an older/different Kraken 2,
  including upstream) with a "rebuild" error; `k2 build` does not reuse it as
  "already built". Rebuild such databases (`k2 build ...`/`kraken2-build ...
  --order GreedyMini`).

Tables now live under [`../../orders/GreedyMini/`](../../orders/GreedyMini/),
shared with kmc's own GreedyMini order (see `orders/README.md`). Only the
`.txt` rank-list sources are tracked in git; `.gm` binaries are generated
from them on demand (`../benchmark/kraken2/benchmark_interface.sh`'s
`bench_ensure_gm_file`). Defaults used by the benchmark scripts:

| File | Role |
|------|------|
| `w5_k15_canonical.gm` | Minimizer + binary GreedyMini (`-g` / `GM_FILE`; default). Centre extension keeps density factor 1.466 (the table's value at ℓ = gm_k = 15) for every valid ℓ = 15, 17, …, 31; see [Density](#density). |
| `w5_k6_canonical.gm` | Mod + binary GreedyMini (`-M` / `MOD_GM_FILE`; default). gm_k=t=6 with Kraken ℓ=31,k=35 (w=5). |

Kraken 2 supports only binary-alphabet GreedyMini tables; a few other tables
that used to live here (including a DNA-alphabet one that turned out to be
non-canonical) are unused and have been moved to `~/unused-orders/kraken2/`
outside the repo. See `orders/README.md` for the full picture.

Classic `.gm` encoding (see `src/greedymini_rank.h`):

```
[uint64 length N][N × uint64 ranks]
```

- `N` is `2^k` (binary alphabet)
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

σ-extension for DNA is **not** written into the file: Kraken 2 remaps
ACGT→ATCG (AT/CG key, residual lex tie-break) when ranking. See `src/greedymini_rank.h`.

## Density

Density factor = (selected minimizers / windows) × (w+1), w = k − ℓ + 1 = 5
(k = ℓ + 4); random order ≈ 2, `w5_k15_canonical` designed 1.466. "Selected"
is what `benchmark_minimizers` counts: successive windows whose minimizer key
changes. Numbers are mean ± std (sample) over 6 independent i.i.d. random DNA
streams of 3 Mb (stream i: `--seed i` and a different tie-break/hash toggle),
spaced seed off (s = 0) unless stated; "old" = the prefix extension of
`575b123`, "new" = centre. At ℓ = 15 = gm_k both are the table itself, i.e. the
design density (1.4662 ± 0.0003 measured). RandomXor is Kraken's default order (xor with the toggle mask; its
density varies with the toggle), RandomPcg a proper random hash.

**Expected density, `w5_k15_canonical`** ([`../../benchmark/kraken2/tables/gm_centre_density_expected.csv`](../../benchmark/kraken2/tables/gm_centre_density_expected.csv))

| ℓ | old (prefix) | new (centre) | RandomXor | RandomPcg |
|---|---|---|---|---|
| 15 (= gm_k, design) | 1.4662 ± 0.0003 | 1.4662 ± 0.0003 | 2.019 ± 0.064 | 1.9993 ± 0.0011 |
| 17 | 1.9853 ± 0.0006 | 1.4659 ± 0.0003 | 2.034 ± 0.058 | 1.9997 ± 0.0005 |
| 19 | 1.7463 ± 0.0009 | 1.4667 ± 0.0003 | 2.012 ± 0.067 | 2.0006 ± 0.0009 |
| 21 | 1.6505 ± 0.0004 | 1.4674 ± 0.0003 | 2.047 ± 0.048 | 1.9997 ± 0.0008 |
| 23 | 1.5147 ± 0.0003 | 1.4668 ± 0.0003 | 2.027 ± 0.062 | 1.9997 ± 0.0009 |
| 25 | **1.4373 ± 0.0005** | 1.4661 ± 0.0004 | 2.001 ± 0.075 | 1.9997 ± 0.0008 |
| 27 | 1.5213 ± 0.0002 | 1.4661 ± 0.0003 | 2.022 ± 0.053 | 1.9996 ± 0.0005 |
| 29 | **1.4510 ± 0.0003** | 1.4662 ± 0.0003 | 2.034 ± 0.046 | 1.9999 ± 0.0011 |
| **31 (default)** | 1.5294 ± 0.0004 | **1.4661 ± 0.0003** | 2.024 ± 0.074 | 2.0001 ± 0.0005 |

The centre order is flat at the table's value for every valid ℓ. The prefix
order is worse than the table at most ℓ (up to +0.52 at ℓ = 17; +0.063, 4.3 %
at the default ℓ = 31) but happens to be *better* than the design value at
ℓ = 25 and 29 (−0.029, −0.015); that is an accident of the offset ℓ − gm_k
between the two strands' prefix windows, not a property to rely on. The old
numbers with the spaced seed on (s = 7, valid for ℓ ≥ 28) are unchanged to
the 4th decimal (ℓ = 29: 1.4512; ℓ = 31: 1.5294) and so are the new ones
(1.4662, 1.4661).

Other shipped tables at their own design window (k = ℓ + w − 1; old cannot
run where `ext + ℓ > 48`): `w11_k13` (w = 11): ℓ = 13 1.5913 ± 0.0004, centre
ℓ = 21 1.5895 ± 0.0005 (old 1.7276 ± 0.0007), ℓ = 31 1.5914 ± 0.0005 (old
unsupported). `w11_k15`: ℓ = 15 1.5236 ± 0.0002; ℓ = 21 centre 1.5248
(old 1.8863), ℓ = 31 centre 1.5231 (old 1.6023). `w20_k9` (w = 20): ℓ = 9
1.7180 ± 0.0007; ℓ = 21 centre 1.7200 (old 1.8016), ℓ = 31 centre 1.7179 (old
unsupported). `w5_k6` is the Mod table (ℓ = gm_k = t, no extension).

**Particular density** (real sequence): the library FASTAs of the benchmark's
`standard_p0.5` Kraken library (`benchmark_minimizers.sh -z`; archaea, bacteria,
human, UniVec_Core, viral; ≈1.42 Gbp, scanned in full, one pass per
scheme), `w5_k15_canonical`, one run per cell ([`../../benchmark/kraken2/tables/gm_centre_density_particular.csv`](../../benchmark/kraken2/tables/gm_centre_density_particular.csv)):

| ℓ | old (prefix) | new (centre) | RandomXor | RandomPcg |
|---|---|---|---|---|
| 15 (= gm_k) | 1.4730 | 1.4730 | 1.954 | 2.000 |
| 17 | 1.9882 | 1.4727 | 2.099 | 1.9998 |
| 21 | 1.6575 | 1.4749 | 2.051 | 2.000 |
| 25 | 1.4513 | 1.4726 | 2.020 | 1.9998 |
| 29 | 1.4598 | 1.4731 | 2.055 | 1.9999 |
| **31 (default)** | 1.5357 | **1.4731** | 2.100 | 2.000 |

(The same holds with the default spaced seed, s = 7, at ℓ = 29 / 31: old
1.4601 / 1.5357, new 1.4731 / 1.4731.) Real data sits ≈ 0.007 above the
random-DNA value for the centre order (0.4–0.5 %) and ≈ 0.006 above it for the
prefix order, so the ranking and the sizes of the gaps carry over from the
random model: at ℓ = 31 the centre order selects 4.1 % fewer minimizers
(1.4731 vs 1.5357), at ℓ = 17 / 21 26 % / 11 % fewer, and 1.4 % / 0.9 % *more*
at ℓ = 25 / 29. The benchmark does not measure database size from density, so
no size change is claimed; fewer selected minimizers means fewer distinct
minimizer lookups per Mbp. `ns/minimizer` of `NextMinimizer` is unchanged
(≈ 93–95 ns centre vs ≈ 100 ns prefix, RandomXor ≈ 91 ns on a loaded machine).

**Reproduce** (run from the repo root; `nice -n 19`, ≤ 4 threads):

```bash
make -C kraken2/src benchmark_minimizers                 # new (centre)
# old (prefix) binary: build the same target from the parent commit
git archive 575b123 kraken2 common | tar -x -C /tmp/old && make -C /tmp/old/kraken2/src benchmark_minimizers
S=benchmark/kraken2/density_sweep.sh; G=orders/GreedyMini/w5_k15_canonical.gm
# expected (6 streams x 3 Mb; -s = minimizer spaces; add -L "29 31" with -s 7)
$S -b /tmp/old/kraken2/src/benchmark_minimizers -g $G -v prefix -o e.csv -s 0 -n 3000000 -S 6 -R
$S -b kraken2/src/benchmark_minimizers          -g $G -v centre -o e.csv -s 0 -n 3000000 -S 6
python3.9 benchmark/kraken2/density_summary.py e.csv     # mean +- std per l
# particular (one stream is enough; -d = the library, -R adds the random baselines)
DB=$KRAKEN2_BENCHMARK_DIR/db/standard_p0.5
for l in 31 25 21 17; do
  $S -b /tmp/old/kraken2/src/benchmark_minimizers -g $G -v prefix -o p.csv -L $l -S 1 -n 100000 -d $DB -j 3 -R
  $S -b kraken2/src/benchmark_minimizers          -g $G -v centre -o p.csv -L $l -S 1 -n 100000 -d $DB -j 3
done
python3.9 benchmark/kraken2/density_summary.py p.csv --col fac_particular
```

`density_sweep.sh` calls, per (ℓ, stream i): `nice -n 19 BIN -S GreedyMini
[Random RandomPcg] -g $G -l ℓ -k ℓ+4 -n 3000000 -j 1 -s SPACES -T TOGGLE_i
--seed i --seed-idx i -d DB -o out.csv`. Unit tests: `make -C kraken2/src test`
(per-base reference over ≥ 1e7 random ℓ-mers, strand invariance, injectivity,
rejections, DB guard, equality with sshash's Python vectors).

