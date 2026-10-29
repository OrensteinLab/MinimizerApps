# Minimizer schemes and their orders

All four apps pick a *substring order* — the ranking used to choose which
k-mer substring (kraken2's minimizer, kmc's signature, sshash's/minimap2's
minimizer) represents a window — from the same small set of options, and
share the C++ implementation of every non-trivial one instead of each
maintaining its own copy. minimap2 is pure C, so it wraps that shared C++
implementation behind a small extern-"C" bridge instead of using it directly
— see [`minimap2/README.fork.md`](../minimap2/README.fork.md) and `minimap2/scheme.cpp`.

## Schemes

| Scheme | What it does |
|---|---|
| default / native | Each app's own default order, not shared code: kraken2 — xor (the canonical ℓ-mer XOR its toggle mask, i.e. the xor order below); kmc — its signature order (lexicographically smallest canonical signature, skipping those that start with AAA or ACA, contain AA or end in TG); sshash — `mixer_64` (multiplicative hash with the FxHash constant ⌊2^64/π⌋ + 1, XOR a seed-derived constant); minimap2 — `hash64()`, Thomas Wang's 64-bit integer hash. |
| `xor` | A random order: XOR the substring with a seed, then rank canonically (with a tie-break for revcomp pairs). |
| `xor-pcg-hash` | The xor order's XOR step, then scrambled through the PCG RXS-M-XS 64-bit output permutation (O'Neill 2014) before ranking. |
| `greedymini` | Ranks substrings by a precomputed lookup table (a "GreedyMini" order, built by the [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool, which greedily constructs low-density minimizer orders) instead of a formula. The table (`gm_k`, `gm_table`) and the DNA->binary projection (`at0-cg1` / `at1-cg0`) are parameters of the order, not part of its name. |
| `mod-minimizer-xor`, `mod-minimizer-greedymini` | The strand-symmetric mod-minimizer (see below) with xor, respectively a GreedyMini table of t-mers, as the order `O_t` on t-mers; `t` is the extra parameter `mod_t`. |

These five names are the only spellings used anywhere (CLI flags, run directories, `meta.json`, the `scheme` column of the CSV/XLSX tables); the old aliases (`hash`, `xor-hash`, `xor-pcg`, `pcg`, `greedymini-k13`, `greedymini-at1`, `mod-xor-hash`, ...) are gone.

## Shared implementation

`common/minimizer-schemes/` (how each order and the mod-minimizer rule are implemented, per app: [`common/minimizer-schemes/README.md`](../common/minimizer-schemes/README.md)):

- `pcg_rxs_m_xs.h` — the PCG RXS-M-XS scramble (`pcg_rxs_m_xs64`), used by all four apps' `xor-pcg-hash` order.
- `greedymini_rank.h` — the full GreedyMini implementation: the bit-manipulation helpers for extracting a DNA k-mer key, the `GreedyMiniTable` LUT class (binary-alphabet tables only; kmc's DNA-alphabet tables are its own, see below), and the `.gm` binary file I/O (`load_gm_file_raw`/`write_gm_file_raw`/`load_gm_table`). All at global C++ scope, so kraken2 (via unqualified lookup from inside `namespace kraken2 {}`), kmc, and sshash (via unqualified lookup from inside `namespace sshash {}`) can all use it directly without a wrapper.
- `greedymini_centre.h` — additive header on top of `greedymini_rank.h`: the centre extension of a canonical table (`gm_centre::Order`, compressed packing for kraken2, injective packing for sshash and minimap2), the encoding-agnostic plane extraction (`TwoBitPlanes`) and rolling Sigma helpers, a slow per-base `Reference` and a density helper. Used by kraken2's `GMRanker`, sshash's `custom_scheme_hasher` and minimap2's `scheme.cpp`; kmc does not use it.
- `mod_minimizer_rule.h` — the strand-symmetric mod-minimizer rule (offset `x mod w` with t == k mod w, the centre tie rule `CentrePick` / `centre_tie_offset`, the O(1) window-vs-reverse-complement comparison `le_revcomp_planes`), shared by kraken2's `ModSampler`, minimap2's `sketch_mod_minimizer` and sshash's `compute_mod_minimizer` / `mod_minimizer_iterator`; each app keeps its own scan.
- `dna_revcomp.h` — reverse-complement of a packed 2-bit-per-base k-mer. Kept in its own header (not `greedymini_rank.h`) because kraken2 already has its own `dna_reverse_complement` in `dna_kmer.h`; combining the two into one header caused an ambiguous-overload compile error under kraken2's `using namespace kraken2;`. kmc includes this header directly (sshash does not: it runs on its own nucleotide encoding, see `sshash/README.md`); kraken2 keeps using its own (verified bit-identical).
- `gm_serializer.cpp` (+ `Makefile`) — builds a `.gm` binary table from a GreedyMini tool's `.txt` rank-list output. C++, not Python, so all four apps share one implementation instead of independent re-serializers.

`benchmark/kmc/tools/gm_sigma_extend.cpp` (+ `Makefile`) — kmc-only C++ tool that expands a binary-alphabet `.gm` into a DNA-alphabet `_dna.gm` table, reusing the shared `gm_extract_dna_key`/`dna_reverse_complement` functions rather than reimplementing the bit math.

`orders/GreedyMini/` — the canonical rank tables all four apps load: each table's `.txt` source (the GreedyMini tool's raw output) alongside its serialized `.gm` (and, where present, `_dna.gm`).

sshash's `custom_scheme_hasher` (`sshash/include/custom_scheme_hasher.hpp`) additionally folds each GreedyMini key into its packed rank as a secondary tie-break below the raw UHS rank. A canonical table deliberately assigns the same rank to a mirror pair of keys (and to every non-UHS key), which kraken2's/kmc's own consumers tolerate but which breaks an invariant sshash's `minimizer_iterator` relies on (a hash tie must imply a value tie) — see that file's comments for the full reasoning. This is sshash-specific packing, not a change to the shared `GreedyMiniTable`/rank semantics itself.

## GreedyMini's two backends

A GreedyMini table can be looked up two ways:

- **Binary-alphabet LUT (default for all four apps)**: a `2^gm_k`-entry table over an AT/CG-collapsed alphabet, with an on-the-fly projection from the real DNA k-mer at lookup time. Compact even for large `gm_k`.
- **DNA-alphabet LUT**: a native `4^gm_k`-entry table, keyed directly by the DNA k-mer — no projection step, but the table is 2^gm_k times larger.

kraken2, sshash and minimap2 support only the binary backend (the shared `GreedyMiniTable` is binary-only; kraken2's `gm_k` values would make a `4^gm_k` DNA-alphabet table impractically large anyway). kmc supports both, with its own DNA-alphabet code, with binary as the default backend (flipped from an earlier DNA-expanded default); `_dna.gm` tables and `-ode` opt into the DNA-alphabet backend for a binary-alphabet source. sshash's minimizer length `m` has no minimizer-length-specific compile-time cap (the practical bound is `Kmer::max_m=32`, `sshash/include/kmer.hpp`, for the 2-bit DNA alphabet — `sshash/include/constants.hpp`'s `min_l`/`max_l` are unrelated, they size the skew index's bucket-size partitioning) — a `4^gm_k` DNA-alphabet table stays practical at that scale, but binary is still the default.

## The centre extension (kraken2, sshash, minimap2): a table built for a shorter length, used at a longer one

kmc uses a table only at its own length: the signature length `-p` must equal the table's `gm_k` (`CCustomOrder::check_k`; `-p 9` with the `w20_k9` table), so it never extends one. kraken2, sshash and minimap2 may rank a length `l > gm_k` with a `gm_k` table: the centre `gm_k` bases decide the order (looked up in the table) and the rest only break ties, so `l - gm_k` must be even. kraken2's binary DNA order for `l > gm_k`, sshash's `--scheme greedymini` and minimap2 (every `-k`, in the injective packing) use this centre extension, implemented once in `common/minimizer-schemes/greedymini_centre.h` (the order, the Sigma/Gamma plane extraction for the apps' base encodings, the rolling Sigma plane, and a per-base reference for tests); sshash's variant is described below. The GreedyMini order is DNA only; kraken2 refuses it for protein databases.

## DNA -> binary projection polarity (all four apps)

The binary backend classifies every base into {A,T} vs {C,G} (the Sigma plane) and uses the Gamma plane as tie-break. The polarity of Sigma is a free choice: `at0-cg1` (default, the historical behaviour, bit-identical hashes) maps A,T->0 and C,G->1; `at1-cg0` flips every Sigma bit (Gamma untouched). The two are the same order up to a relabelling of table keys, so on i.i.d. uniform DNA the expected density is the same, but a table is not symmetric under key complement, so the *particular* density on a real (AT/CG-skewed) collection can differ. sshash: `sshash build --gm-projection {at0-cg1,at1-cg0}` (stored in the index, `sshash query --verbose` prints it), implemented in `sshash/include/custom_scheme_hasher.hpp`. minimap2: `--gm-projection {at0-cg1,at1-cg0}` with `--scheme greedymini` (stored in the `.mmi`), implemented in `minimap2/scheme.cpp` through the shared header's `TwoBitPlanes` flip. kraken2: `k2 build --gm-projection {at0-cg1,at1-cg0}` (`build_db`/`estimate_capacity` `-P`; stored in `opts.k2d`, so `classify` follows the database), for `--order greedymini` with either sketching, implemented in `src/minimizer_rankers.h` (the `TwoBitPlanes` flip in the centre extension, a complemented table key when `l == gm_k`). kmc: `-oat1` with `-co` and a binary `.gm` (kmc ranks by the binary key and breaks ties on the DNA, so the flip is a one-time permutation `rank'[key] = rank[~key]` of the binary table, which the sidecar stores; a DNA `_dna.gm` table is refused). `sshash/test/test_gm_centre.cpp` checks that `hash()` stays injective in the m-mer and reverse-complement invariant for both polarities; kraken2's `test_scheme_rankers` checks that the flipped order of x equals the default order of x with the bases relabelled A<->C, G<->T (which complements Sigma and keeps Gamma).

### sshash: centre extension of the table to the m-mer

sshash's `--scheme greedymini` (binary-alphabet tables only) extends an explicit canonical (mirror-symmetric) table of `gm_k`-mers to `m = gm_k + 2i`: write the m-mer as `xyz` with `|x| = |z| = i` and `|y| = gm_k`, rank the **binary projection of `y`** by the table (one lookup; strand invariant since `rc` maps `y` to its mirror and `rank[key] == rank[reverse(key)]`), then break ties by the lexicographic minimum of the projection of `xyz` and its reversal, and finally by the Gamma plane (so the hash stays injective on canonical m-mers). It needs `m - gm_k` even (at `m = 21` only the odd-`gm_k` tables `w11_k13` and `w11_k15`; an even `gm_k` is rejected with a clear error), the binary alphabet and a mirror-symmetric table (checked at configure time). On i.i.d. DNA at `m = 21, w = 11` it gives the tables' designed densities (1.592 / 1.524). The harness runs it as `greedymini` once per (table, projection) variant, e.g. `w11_k13` / `w11_k15` x `at0-cg1` / `at1-cg0` (`gm_k` = 13 / 15 in the run directory, not in the order name). sshash's `--scheme mod-minimizer-greedymini` ranks the k-mer's t-mers with a table of t-mers by the same hash (`t = gm_k + 2i`, in practice `t = gm_k`, so an even `gm_k` such as `w22_k10` is fine there).

minimap2's `--scheme greedymini` uses the same centre order on the k-mer (`-k = gm_k + 2i`, `-k <= 23`), in the injective packing, on minimap2's ACGT = 0123 encoding (`TwoBitPlanes<true>`), with `--order-seed` toggling the Gamma tie-break as kraken2's and sshash's greedymini seeds do (it reorders only k-mers with equal rank and Sigma plane: ~1% of the minimizer positions of a real genome at `-k 19`); the `w10_k15` table at minimap2's `-k 19` keeps its designed density (1.52 on i.i.d. DNA; ~1.55 under homopolymer compression, `map-pb`'s `-H`, which makes the projected bit sequence non-i.i.d.).

## Usage

**kraken2** (`k2 build` / `kraken2-build`, or the `build_db`/`estimate_capacity` binaries directly):

```bash
k2 build --sketching minimizer --order xor ...                                      # default
k2 build --sketching minimizer --order xor-pcg-hash ...
k2 build --sketching minimizer --order greedymini \
         --greedymini-file orders/GreedyMini/w5_k15_canonical.gm [--gm-projection at1-cg0]
k2 build --sketching mod-minimizer --order xor ...                                  # mod-minimizer-xor
k2 build --sketching mod-minimizer --order greedymini --greedymini-file ... ...     # mod-minimizer-greedymini
```

Only exact lowercase canonical names are accepted (no aliases). The benchmark
harness's scheme (`-S`) is one of `xor`, `xor-pcg-hash`, `greedymini`,
`mod-minimizer-xor`, `mod-minimizer-greedymini`, i.e. `--sketching` + `--order`;
`-P at0-cg1|at1-cg0` (env `GM_PROJECTION`) sets `--gm-projection` for both
greedymini schemes. The flag names differ per app (kraken2 `--greedymini-file`,
sshash `--greedy-mini-file`, kmc `-co`, minimap2 `--greedymini-file`).

**kmc**:

```bash
kmc -or xor -ol<label> [-os<seed>] ...                         # xor (-ol is mandatory with -or)
kmc -or xor-pcg-hash -ol<label> [-os<seed>] ...                # xor-pcg-hash
kmc -co<file>.gm -ol<label> [-os<seed>] [-ode] [-oat1] ...     # greedymini (binary .gm, or _dna.gm / -ode for kmc's DNA-alphabet backend; -oat1 = flipped projection, binary .gm only)
kmc ...                                                        # no -or/-co: kmc's own lex order (not shared code)
```

The harness (`kmc3_bench.sh -S default|lex|xor|xor-pcg-hash|greedymini`) adds
`--gm-lut dna|binary` (the table backend; same order, different storage) and
`GM_PROJECTION=at1-cg0` (= `-oat1`, needs `gm_lut=binary`). `lex` is kmc's true
canonical lex baseline, run as `-or xor -os 0`.

**sshash**:

```bash
./sshash build -i FILE -k 31 -m 21 --scheme native ...                                  # default, unchanged
./sshash build -i FILE -k 31 -m 21 --scheme xor -s SEED ...
./sshash build -i FILE -k 31 -m 21 --scheme xor-pcg-hash -s SEED ...
./sshash build -i FILE -k 31 -m 21 --scheme greedymini \
    --greedy-mini-file orders/GreedyMini/w11_k15_canonical.gm [--gm-projection at1-cg0] ...
./sshash build -i FILE -k 31 -m 21 --scheme mod-minimizer-xor --mod-t T -s SEED ...
./sshash build -i FILE -k 31 -m 21 --scheme mod-minimizer-greedymini \
    --greedy-mini-file orders/GreedyMini/w5_k10_canonical.gm --mod-t 10 ...
```

The harness takes the scheme list with `-S` (`native`, `xor`, ...) and sweeps the
GreedyMini (table, projection) variants (`--gm-tables`, `--gm-projections`;
default: config.yaml's).

**minimap2** (long-only flags: `-O`/`-s` are taken by minimap2 itself, for gap-open penalty and min DP score):

```bash
./minimap2 -a ref.fa reads.fq                                                     # native, unchanged
./minimap2 -a --scheme xor --order-seed SEED ref.fa reads.fq
./minimap2 -a --scheme xor-pcg-hash --order-seed SEED ref.fa reads.fq
./minimap2 -a --scheme greedymini \
    --greedymini-file orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq       # binary alphabet only; centre-extended to -k (15, 17, ..., 23)
./minimap2 -a --scheme greedymini --gm-projection at1-cg0 \
    --greedymini-file orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq       # flipped AT/CG projection (see above)

# mod-minimizer: O_t is xor + --mod-t (minimap2's own -k/-w already match
# the paper's m/w notation directly).
./minimap2 -a --scheme mod-minimizer-xor --order-seed SEED --mod-t T ref.fa reads.fq
```

### CLI flag mapping

| Scheme | kraken2 | kmc | sshash | minimap2 |
|---|---|---|---|---|
| default/native | `--order xor` (default: XOR toggle mask) | no `-or`/`-co` (default: KMC's signature order) | `--scheme native` / no `--scheme` (default: `mixer_64`, multiplicative hash with the FxHash constant ⌊2^64/π⌋ + 1) | no `--scheme` (default: `hash64()`, Wang's integer hash) |
| xor | `--sketching minimizer --order xor` (default) | `-or xor -ol<label> [-os<seed>]` | `--scheme xor -s<seed>` | `--scheme xor --order-seed SEED` |
| xor-pcg-hash | `--sketching minimizer --order xor-pcg-hash` | `-or xor-pcg-hash -ol<label> [-os<seed>]` | `--scheme xor-pcg-hash -s<seed>` | `--scheme xor-pcg-hash --order-seed SEED` |
| greedymini | `--sketching minimizer --order greedymini --greedymini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` | `-co FILE.gm -ol<label> [-os<seed>] [-ode] [-oat1]` | `--scheme greedymini --greedy-mini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` (binary alphabet only) | `--scheme greedymini --greedymini-file FILE.gm [--order-seed SEED] [--gm-projection at0-cg1\|at1-cg0]` (binary alphabet only; centre extension to `-k = gm_k + 2i`; the seed toggles the Gamma tie-break) |
| mod-minimizer-xor (strand-symmetric, see below) | `--sketching mod-minimizer --order xor` | (n/a — kmc has no analogous two-stage scheme) | `--scheme mod-minimizer-xor --mod-t T -s<seed>` | `--scheme mod-minimizer-xor --order-seed SEED --mod-t T` (O_t = xor only) |
| mod-minimizer-greedymini | `--sketching mod-minimizer --order greedymini --greedymini-file FILE.gm [--gm-projection ...]` | (n/a) | `--scheme mod-minimizer-greedymini --greedy-mini-file FILE.gm --mod-t T [--gm-projection ...]` (a table of t-mers) | (n/a) |

### Run directories and result columns

A scheme is identified by five separate fields, never packed into one string
(`benchmark/common/scheme_names.py`): `scheme` (one of the five names above, plus
an app baseline: `native` for sshash/minimap2, `default`/`lex` for kmc),
`gm_k` (GreedyMini table k), `projection` (`at0-cg1` | `at1-cg0`), `gm_table`
(table stem, e.g. `w5_k15_canonical`) and `mod_t` (mod-minimizer t-mer length).
They are blank where they do not apply. kmc adds `gm_lut` (`dna` | `binary`,
the table backend). Every harness writes them as the columns
`scheme, gm_k, projection, gm_table, mod_t` (kmc: `scheme, gm_k, projection,
gm_table, gm_lut`) of its CSV/XLSX tables, in `meta.json` and in kraken2's run
manifest.

The run directory of a scheme is `scheme_names.scheme_dir_name(scheme, gm_k,
projection, mod_t)`: `<scheme>[/k<gm_k>_<projection>_t<mod_t>]` (the parts that
apply, in that order), followed by the app's own segments:

| Scheme | Run-dir segment |
|---|---|
| `xor`, `xor-pcg-hash`, `native`, `default`, `lex` | the name itself |
| `greedymini` | `greedymini/k15_at0-cg1`, `greedymini/k15_at1-cg0` |
| `mod-minimizer-xor` | `mod-minimizer-xor/t6` |
| `mod-minimizer-greedymini` | `mod-minimizer-greedymini/k6_at0-cg1_t6` |

kmc appends `/lut_<gm_lut>` for greedymini (`greedymini/k9_at1-cg0/lut_binary`,
then `/<seed_idx>`).

The mod-minimizer of all three apps is a strand-symmetric variant of
mod-sampling (Groot Koerkamp & Pibiri, WABI 2024); the scans are
`kraken2/src/mod_minimizer.h`, `minimap2/scheme.cpp` and `sshash/include/util.hpp`
(`compute_mod_minimizer`; sshash's window is the k-mer, `w = k-m+1`, `T = k-t+1`). With
t == m (mod w) its density is mod-sampling's; any other `t` is rejected (for such t the
candidate offsets of a window and its reverse complement are not mirror images; sshash, k=31,
m=21: a t=11 GreedyMini table gave ~2.2 against 1.46 at t=10 with `w22_k10`).

**Tie rules, plain minimizers.** Each application keeps its native rule for breaking ties
between equal-rank positions, because each rule serves how that application uses the
minimizer. In every order we use, the rank is a total order on canonical k-mers (a hash of
kappa(x) = min(x, rc(x)), or GreedyMini's injective packing; in kraken2 the packing is not
injective on its own and the strand-invariant key `dna_hash_key ^ toggle` completes it), so
equal ranks occur only when the same canonical k-mer appears more than once in a window.
kraken2 keeps the leftmost tied position (its monotone deque); it uses only the minimizer's
value, as the hash-table key whose LCA it looks up, so the rule does not change its output.
kmc moves the signature to the newest equal occurrence: a super-k-mer ends when its
signature's position leaves the k-mer, so keeping the newest occurrence keeps the signature
longest, giving the longest super-k-mers and the least bin I/O; bins are chosen by the
canonical value, so counts do not depend on the tie position (at the level of values this is
robust winnowing, Schleimer, Wilkerson & Aiken, SIGMOD 2003, Definition 3; the KMC 2 and
KMC 3 papers do not discuss ties). minimap2 emits every tied position (standard winnowing;
Winnowmap, Jain et al. 2020, replaces it by robust winnowing) and skips k-mers equal to their
reverse complement, which is strand-symmetric without a tie rule; it is also the upstream
baseline. sshash stores super-k-mers by minimizer and offset, so the position matters: the
scheme must be canonical (a k-mer and its reverse complement find the same minimizer
occurrence) and forward (its builder checks it), which no leftmost/rightmost-by-strand rule
achieves; it uses the centre rule of Cologni and Pibiri (SSHash v6.0.0, PR 93): the tied
position closest to the centre of the k-mer, and between two equally close ones the left if
v <= RC(v), else the right. GreedyMini changes only the rank, so every app applies its native
rule to it.

**Tie rule, mod-minimizer** (`common/minimizer-schemes/mod_minimizer_rule.h`, shared by
kraken2 `ModSampler`, minimap2 `sketch_mod_minimizer` and sshash `compute_mod_minimizer` /
`mod_minimizer_iterator`). Here the tie rule changes the output in every app: the window
minimum is over t-mers, and the k-mer sampled for the minimal t-mer at x is the one at offset
p = x mod w, so equal-rank t-mers in different residue classes select different k-mers. All
apps require t == k (mod w) (kraken2 `t = 4 + (l-4) mod w`; minimap2 `--mod-t` and sshash
`--mod-t` are checked), so the mirror of x in RC(v) lies in class w-1-p and the candidate
offsets are mirror-symmetric. The centre rule is applied to them: the candidate closest to the
window centre (w-1)/2, and between the two equally close mirror candidates p and w-1-p, the
left if v <= RC(v), else the right. RC(v) then samples the mirror k-mer (canonical), and the
chosen offset never moves backwards as the window slides (forward). This replaced the G/T-strand
rule of SimdMinimizers (Groot Koerkamp & Martayan, bioRxiv 2025.01.27.634998: the leftmost
minimum on the strand with more G/T bases), which is canonical but not forward and needs a
strand decision in every window; the centre rule does extra work only on a tied minimum, and
compares v with RC(v) only on an exact two-way tie. The native plain rules do not carry over:
leftmost (kraken2) is not canonical for the mod-minimizer, keeping every tied k-mer (minimap2)
inflates density at small t, and kmc has no mod-minimizer. v <= RC(v) is lexicographic in
kraken2 and minimap2 (`le_revcomp_planes`, O(1) on two bit planes) and sshash's packed-word
comparison `!(kmer_rc < kmer)` in sshash, the one its plain rule uses; both are total orders,
which is all canonicity needs. Indexes built with the G/T rule are refused (kraken2
`opts.k2d` `mod_tie_rule`, minimap2 `.mmi` order section, sshash index version 15).

## Regenerating a `.gm` table

The [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `.txt` rank-list output is the source of truth; the `.gm` binary is a derived, serialized form:

```bash
cd common/minimizer-schemes && make               # builds gm_serializer
./gm_serializer [--canonical] -k N -o OUT.gm IN.txt

# kmc only, to also produce a DNA-alphabet table from a binary one:
cd benchmark/kmc/tools && make                    # builds gm_sigma_extend
./gm_sigma_extend --seed N -o OUT_dna.gm IN.gm
```
