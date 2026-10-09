# Minimizer/signature orders

All four apps pick a *substring order* — the ranking used to choose which
k-mer substring (kraken2's minimizer, kmc's signature, sshash's/minimap2's
minimizer) represents a window — from the same small set of options, and
share the C++ implementation of every non-trivial one instead of each
maintaining its own copy. minimap2 is pure C, so it wraps that shared C++
implementation behind a small extern-"C" bridge instead of using it directly
— see [`minimap2/README.md`](../minimap2/README.md) and `minimap2/order.cpp`.

## Order modes

| Order | What it does |
|---|---|
| default / native | Each app's own default order, not shared code: kraken2 — xor (the canonical ℓ-mer XOR its toggle mask, i.e. the xor order below); kmc — its signature order (lexicographically smallest canonical signature, skipping those that start with AAA or ACA, contain AA or end in TG); sshash — `mixer_64` (multiplicative hash with the FxHash constant ⌊2^64/π⌋ + 1, XOR a seed-derived constant); minimap2 — `hash64()`, Thomas Wang's 64-bit integer hash. |
| xor (`xor`, `hash`, `xor-hash`) | A random order: XOR the substring with a seed, then rank canonically (with a tie-break for revcomp pairs). |
| xor-pcg (`xor-pcg`, `pcg`, `xor-pcg-hash`) | The xor order's XOR step, then scrambled through the PCG RXS-M-XS 64-bit output permutation (O'Neill 2014) before ranking. |
| GreedyMini / `greedymini` | Ranks substrings by a precomputed lookup table (a "GreedyMini" order, built by the [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool, which greedily constructs low-density minimizer orders) instead of a formula. |

## Shared implementation

`common/minimizer-orders/`:

- `pcg_rxs_m_xs.h` — the PCG RXS-M-XS scramble (`pcg_rxs_m_xs64`), used by all four apps' xor-pcg order.
- `greedymini_rank.h` — the full GreedyMini implementation: the bit-manipulation helpers for extracting a DNA k-mer key, the `GreedyMiniTable` LUT class (binary-alphabet tables only; kmc's DNA-alphabet tables are its own, see below), and the `.gm` binary file I/O (`load_gm_file_raw`/`write_gm_file_raw`/`load_gm_table`). All at global C++ scope, so kraken2 (via unqualified lookup from inside `namespace kraken2 {}`), kmc, and sshash (via unqualified lookup from inside `namespace sshash {}`) can all use it directly without a wrapper.
- `greedymini_centre.h` — additive header on top of `greedymini_rank.h`: the centre extension of a canonical table (`gm_centre::Order`, compressed packing for kraken2, injective packing for sshash and minimap2), the encoding-agnostic plane extraction (`TwoBitPlanes`) and rolling Sigma helpers, a slow per-base `Reference` and a density helper. Used by kraken2's `GMRanker`, sshash's `custom_order_hasher` and minimap2's `order.cpp`; kmc does not use it.
- `mod_minimizer_rule.h` — the strand-symmetric mod-minimizer rule (preferred strand by G/T count and its tie-break, the two offset rules, the residue-class mirror slot, sshash's centre tie rule), shared by kraken2's `ModSampler`, minimap2's `sketch_mod_minimizer` and sshash's `compute_mod_minimizer` / `mod_minimizer_iterator`; each app keeps its own scan.
- `dna_revcomp.h` — reverse-complement of a packed 2-bit-per-base k-mer. Kept in its own header (not `greedymini_rank.h`) because kraken2 already has its own `dna_reverse_complement` in `dna_kmer.h`; combining the two into one header caused an ambiguous-overload compile error under kraken2's `using namespace kraken2;`. kmc includes this header directly (sshash does not: it runs on its own nucleotide encoding, see `sshash/README.md`); kraken2 keeps using its own (verified bit-identical).
- `gm_serializer.cpp` (+ `Makefile`) — builds a `.gm` binary table from a GreedyMini tool's `.txt` rank-list output. C++, not Python, so all four apps share one implementation instead of independent re-serializers.

`benchmark/kmc/tools/gm_sigma_extend.cpp` (+ `Makefile`) — kmc-only C++ tool that expands a binary-alphabet `.gm` into a DNA-alphabet `_dna.gm` table, reusing the shared `gm_extract_dna_key`/`dna_reverse_complement` functions rather than reimplementing the bit math.

`orders/GreedyMini/` — the canonical rank tables all four apps load: each table's `.txt` source (the GreedyMini tool's raw output) alongside its serialized `.gm` (and, where present, `_dna.gm`).

sshash's `custom_order_hasher` (`sshash/include/custom_order_hasher.hpp`) additionally folds each GreedyMini key into its packed rank as a secondary tie-break below the raw UHS rank. A canonical table deliberately assigns the same rank to a mirror pair of keys (and to every non-UHS key), which kraken2's/kmc's own consumers tolerate but which breaks an invariant sshash's `minimizer_iterator` relies on (a hash tie must imply a value tie) — see that file's comments for the full reasoning. This is sshash-specific packing, not a change to the shared `GreedyMiniTable`/rank semantics itself.

## GreedyMini's two backends

A GreedyMini table can be looked up two ways:

- **Binary-alphabet LUT (default for all four apps)**: a `2^gm_k`-entry table over an AT/CG-collapsed alphabet, with an on-the-fly projection from the real DNA k-mer at lookup time. Compact even for large `gm_k`.
- **DNA-alphabet LUT**: a native `4^gm_k`-entry table, keyed directly by the DNA k-mer — no projection step, but the table is 2^gm_k times larger.

kraken2, sshash and minimap2 support only the binary backend (the shared `GreedyMiniTable` is binary-only; kraken2's `gm_k` values would make a `4^gm_k` DNA-alphabet table impractically large anyway). kmc supports both, with its own DNA-alphabet code, with binary as the default backend (flipped from an earlier DNA-expanded default); `_dna.gm` tables and `-ode` opt into the DNA-alphabet backend for a binary-alphabet source. sshash's minimizer length `m` has no minimizer-length-specific compile-time cap (the practical bound is `Kmer::max_m=32`, `sshash/include/kmer.hpp`, for the 2-bit DNA alphabet — `sshash/include/constants.hpp`'s `min_l`/`max_l` are unrelated, they size the skew index's bucket-size partitioning) — a `4^gm_k` DNA-alphabet table stays practical at that scale, but binary is still the default.

## The centre extension (kraken2, sshash, minimap2): a table built for a shorter length, used at a longer one

kmc uses a table only at its own length: the signature length `-p` must equal the table's `gm_k` (`CCustomOrder::check_k`; `-p 9` with the `w20_k9` table), so it never extends one. kraken2, sshash and minimap2 may rank a length `l > gm_k` with a `gm_k` table: the centre `gm_k` bases decide the order (looked up in the table) and the rest only break ties, so `l - gm_k` must be even. kraken2's binary DNA order for `l > gm_k`, sshash's `-O greedymini` and minimap2 (every `-k`, in the injective packing) use this centre extension, implemented once in `common/minimizer-orders/greedymini_centre.h` (the order, the Sigma/Gamma plane extraction for the apps' base encodings, the rolling Sigma plane, and a per-base reference for tests); sshash's variant is described below. The GreedyMini order is DNA only; kraken2 refuses it for protein databases.

## DNA -> binary projection polarity (all four apps)

The binary backend classifies every base into {A,T} vs {C,G} (the Sigma plane) and uses the Gamma plane as tie-break. The polarity of Sigma is a free choice: `at0-cg1` (default, the historical behaviour, bit-identical hashes) maps A,T->0 and C,G->1; `at1-cg0` flips every Sigma bit (Gamma untouched). The two are the same order up to a relabelling of table keys, so on i.i.d. uniform DNA the expected density is the same, but a table is not symmetric under key complement, so the *particular* density on a real (AT/CG-skewed) collection can differ. sshash: `sshash build --gm-projection {at0-cg1,at1-cg0}` (stored in the index, `sshash query --verbose` prints it), implemented in `sshash/include/custom_order_hasher.hpp`. minimap2: `--gm-projection {at0-cg1,at1-cg0}` with `--minimizer-order greedymini` (stored in the `.mmi`), implemented in `minimap2/order.cpp` through the shared header's `TwoBitPlanes` flip. kraken2: `k2 build --gm-projection {at0-cg1,at1-cg0}` (`build_db`/`estimate_capacity` `-P`; stored in `opts.k2d`, so `classify` follows the database), for `--order greedymini` with either sketching, implemented in `src/minimizer_rankers.h` (the `TwoBitPlanes` flip in the centre extension, a complemented table key when `l == gm_k`). kmc: `-oat1` with `-co` and a binary `.gm` (kmc ranks by the binary key and breaks ties on the DNA, so the flip is a one-time permutation `rank'[key] = rank[~key]` of the binary table, which the sidecar stores; a DNA `_dna.gm` table is refused). `sshash/test/test_gm_centre.cpp` checks that `hash()` stays injective in the m-mer and reverse-complement invariant for both polarities; kraken2's `test_scheme_rankers` checks that the flipped order of x equals the default order of x with the bases relabelled A<->C, G<->T (which complements Sigma and keeps Gamma).

### sshash: centre extension of the table to the m-mer

sshash's `-O greedymini` (binary-alphabet tables only) extends an explicit canonical (mirror-symmetric) table of `gm_k`-mers to `m = gm_k + 2i`: write the m-mer as `xyz` with `|x| = |z| = i` and `|y| = gm_k`, rank the **binary projection of `y`** by the table (one lookup; strand invariant since `rc` maps `y` to its mirror and `rank[key] == rank[reverse(key)]`), then break ties by the lexicographic minimum of the projection of `xyz` and its reversal, and finally by the Gamma plane (so the hash stays injective on canonical m-mers). It needs `m - gm_k` even (at `m = 21` only the odd-`gm_k` tables `w11_k13` and `w11_k15`; an even `gm_k` is rejected with a clear error), the binary alphabet and a mirror-symmetric table (checked at configure time). On i.i.d. DNA at `m = 21, w = 11` it gives the tables' designed densities (1.592 / 1.524). The harness names the orders `greedymini-k13`, `greedymini-k15` (+ `-at1`). sshash's `-O mod-greedymini` ranks the k-mer's t-mers with a table of t-mers by the same hash (`t = gm_k + 2i`, in practice `t = gm_k`, so an even `gm_k` such as `w22_k10` is fine there).

minimap2's `--minimizer-order greedymini` uses the same centre order on the k-mer (`-k = gm_k + 2i`, `-k <= 23`), in the injective packing, on minimap2's ACGT = 0123 encoding (`TwoBitPlanes<true>`); the `w10_k15` table at minimap2's `-k 19` keeps its designed density (1.52 on i.i.d. DNA; ~1.55 under homopolymer compression, `map-pb`'s `-H`, which makes the projected bit sequence non-i.i.d.).

## Usage

**kraken2** (`k2 build` / `kraken2-build`, or the `build_db`/`estimate_capacity` binaries directly):

```bash
k2 build --sketching minimizer --order xor ...                                      # default
k2 build --sketching minimizer --order xor-pcg ...
k2 build --sketching minimizer --order greedymini \
         --greedymini-file orders/GreedyMini/w5_k15_canonical.gm
k2 build --sketching mod-minimizer --order greedymini --greedymini-file ... ...     # mod-minimizer + any order above
```

Only exact lowercase names are accepted (no aliases). The scheme label used by
the benchmark harness is `xor`, `xor-pcg`, `greedymini`, `mod-minimizer-xor`
or `mod-minimizer-greedymini` (`--sketching` + `--order`). The kraken2 flag is
`--greedymini-file`; the sshash flag is still `--greedy-mini-file` and kmc has
its own `-co` / `-or*` options, and minimap2 `--greedymini-file`: the apps' flag names differ.

**kmc**:

```bash
kmc -orhash -ol<label> [-os<seed>] ...                        # xor (-ol is mandatory with -or)
kmc -orpcg -ol<label> [-os<seed>] ...                         # xor-pcg
kmc -co<file>.gm -ol<label> [-os<seed>] [-ode] ...             # GreedyMini (binary .gm, or _dna.gm/-ode for kmc's DNA-alphabet backend)
kmc ...                                                        # no -or/-co: kmc's own lex order (not shared code)
```

**sshash**:

```bash
./sshash build -i FILE -k 31 -m 21 -O native ...                                  # default, unchanged
./sshash build -i FILE -k 31 -m 21 -O xor-hash -s SEED ...
./sshash build -i FILE -k 31 -m 21 -O xor-pcg-hash -s SEED ...
./sshash build -i FILE -k 31 -m 21 -O greedymini \
    --greedy-mini-file orders/GreedyMini/w11_k15_canonical.gm ...
```

**minimap2** (long-only flags: `-O`/`-s` were already taken by minimap2 itself, for gap-open penalty and min DP score):

```bash
./minimap2 -a ref.fa reads.fq                                                     # native, unchanged
./minimap2 -a --minimizer-order xor --order-seed SEED ref.fa reads.fq
./minimap2 -a --minimizer-order xor-pcg-hash --order-seed SEED ref.fa reads.fq
./minimap2 -a --minimizer-order greedymini \
    --greedymini-file orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq       # binary alphabet only; centre-extended to -k (15, 17, ..., 23)
./minimap2 -a --minimizer-order greedymini --gm-projection at1-cg0 \
    --greedymini-file orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq       # flipped AT/CG projection (see below)

# mod-minimizer: O_t is xor + --mod-t (minimap2's own -k/-w already match
# the paper's m/w notation directly).
./minimap2 -a --minimizer-order mod-minimizer-xor --order-seed SEED --mod-t T ref.fa reads.fq
```

### CLI flag mapping

| Order | kraken2 | kmc | sshash | minimap2 |
|---|---|---|---|---|
| default/native | `--order xor` (default: XOR toggle mask) | no `-or`/`-co` (default: KMC's signature order) | no `-O` (default: `mixer_64`, multiplicative hash with the FxHash constant ⌊2^64/π⌋ + 1) | no `--minimizer-order` (default: `hash64()`, Wang's integer hash) |
| xor | `--order xor` (default) | `-orhash -ol<label> [-os<seed>]` | `-O xor-hash -s<seed>` | `--minimizer-order xor --order-seed SEED` |
| xor-pcg | `--order xor-pcg` | `-orpcg -ol<label> [-os<seed>]` | `-O xor-pcg-hash -s<seed>` | `--minimizer-order xor-pcg-hash --order-seed SEED` |
| greedymini | `--order greedymini --greedymini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` | `-co FILE.gm -ol<label> [-os<seed>] [-ode] [-oat1]` | `-O greedymini --greedy-mini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` (binary alphabet only) | `--minimizer-order greedymini --greedymini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` (binary alphabet only; centre extension to `-k = gm_k + 2i`) |
| mod-minimizer sketching (strand-symmetric, see below) | `--sketching mod-minimizer` (combine with any `--order` above) | (n/a — kmc has no analogous two-stage scheme) | `-O mod-xor-hash --mod-t T -s<seed>` or `-O mod-greedymini --greedy-mini-file FILE.gm --mod-t T` (a table of t-mers) | `--minimizer-order mod-minimizer-xor --order-seed SEED --mod-t T` (O_t = xor only) |

The mod-minimizer of all three apps is a strand-symmetric variant of
mod-sampling (Groot Koerkamp & Pibiri, WABI 2024). Its canonical definition is
simd-minimizers' (Groot Koerkamp & Martayan, "SimdMinimizers: Computing random
minimizers, fast", [bioRxiv 2025.01.27.634998](https://www.biorxiv.org/content/10.1101/2025.01.27.634998v1.full.pdf)):
each window prefers the strand with more G/T bases; applying it to
mod-sampling and the tie rules (window <= its reverse complement on equal G/T
counts, ties at the minimum below) are this repo's. A forward-preferred window samples the k-mer at
offset `x_L mod w` (`x_L` = leftmost minimal t-mer), a reverse-preferred one at
`(w-1) - (((T-1) - x_R) mod w)` (`x_R` = rightmost minimal t-mer, `T` t-mers per
window), so a read and its reverse complement select the same k-mers.
The rule is `common/minimizer-orders/mod_minimizer_rule.h`; the scans are
`kraken2/src/mod_minimizer.h`, `minimap2/order.cpp` and
`sshash/include/util.hpp` (`compute_mod_minimizer`; sshash's window is the k-mer, `w = k-m+1`,
`T = k-t+1`). With the usual `t = 4 + (m - 4) mod w` (`(T-1) ≡ w-1 mod w`) both rules reduce
to `x mod w` and the strands differ only on ties, so the density is mod-sampling's. Ties
(equal minimal canonical t-mers): minimap2 and kraken2 take the leftmost tied t-mer on a
forward-preferred window and the rightmost otherwise; sshash, whose super-k-mers need a
forward scheme, takes the tied t-mer whose sampled offset is closest to the window centre
`(w-1)/2` (a mirror pair by strand), which keeps the scheme forward (checked by its builder).
For any other `t` the strand rules disagree and the density degrades (sshash, k=31, m=21:
a t=11 GreedyMini table gave ~2.2 against 1.46 at t=10 with `w22_k10`), so only t = m mod w
is used.

## Regenerating a `.gm` table

The [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `.txt` rank-list output is the source of truth; the `.gm` binary is a derived, serialized form:

```bash
cd common/minimizer-orders && make               # builds gm_serializer
./gm_serializer [--canonical] -k N -o OUT.gm IN.txt

# kmc only, to also produce a DNA-alphabet table from a binary one:
cd benchmark/kmc/tools && make                    # builds gm_sigma_extend
./gm_sigma_extend --seed N -o OUT_dna.gm IN.gm
```
