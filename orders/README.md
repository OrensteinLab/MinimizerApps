# Minimizer/signature orders

All three apps pick a *substring order* — the ranking used to choose which
k-mer substring (kraken2's minimizer, kmc's signature, sshash's minimizer)
represents a window — from the same small set of options, and share the C++
implementation of every non-trivial one instead of each maintaining its own
copy.

## Order modes

| Order | What it does |
|---|---|
| default / lex / native | Each app's own baseline order (kraken2's original hash-based minimizer order; kmc's canonical lexicographic signature order; sshash's own `mixer_64`-based order). Not shared code — this is each app's pre-existing behavior. |
| RandomXor / `hash` / `xor-hash` | A random order: XOR the substring with a seed, then rank canonically (with a tie-break for revcomp pairs). |
| RandomPcg / `pcg` / `xor-pcg-hash` | RandomXor's XOR step, then scrambled through the PCG RXS-M-XS 64-bit output permutation (O'Neill 2014) before ranking. |
| GreedyMini / `greedymini` | Ranks substrings by a precomputed lookup table (a "GreedyMini" order, built by the [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool, which greedily constructs low-density minimizer orders) instead of a formula. |

## Shared implementation

`common/minimizer-orders/`:

- `pcg_rxs_m_xs.h` — the PCG RXS-M-XS scramble (`pcg_rxs_m_xs64`), used by all three apps' RandomPcg order.
- `greedymini_rank.h` — the full GreedyMini implementation: the bit-manipulation helpers for extracting a DNA k-mer key, the `GreedyMiniTable` LUT class (both alphabet backends, see below), and the `.gm` binary file I/O (`load_gm_file_raw`/`write_gm_file_raw`/`load_gm_table`). All at global C++ scope, so kraken2 (via unqualified lookup from inside `namespace kraken2 {}`), kmc, and sshash (via unqualified lookup from inside `namespace sshash {}`) can all use it directly without a wrapper.
- `dna_revcomp.h` — reverse-complement of a packed 2-bit-per-base k-mer. Kept in its own header (not `greedymini_rank.h`) because kraken2 already has its own `dna_reverse_complement` in `dna_kmer.h`; combining the two into one header caused an ambiguous-overload compile error under kraken2's `using namespace kraken2;`. kmc and sshash include this header directly; kraken2 keeps using its own (verified bit-identical).
- `gm_serializer.cpp` (+ `Makefile`) — builds a `.gm` binary table from a GreedyMini tool's `.txt` rank-list output. C++, not Python, so all three apps share one implementation instead of independent re-serializers.

`benchmark/kmc/tools/gm_sigma_extend.cpp` (+ `Makefile`) — kmc-only C++ tool that expands a binary-alphabet `.gm` into a DNA-alphabet `_dna.gm` table, reusing the shared `gm_extract_dna_key`/`dna_reverse_complement` functions rather than reimplementing the bit math.

`orders/GreedyMini/` — the canonical rank tables all three apps load: each table's `.txt` source (the GreedyMini tool's raw output) alongside its serialized `.gm` (and, where present, `_dna.gm`).

sshash's `custom_order_hasher` (`sshash/include/custom_order_hasher.hpp`) additionally folds each GreedyMini key into its packed rank as a secondary tie-break below the raw UHS rank. A canonical table deliberately assigns the same rank to a mirror pair of keys (and to every non-UHS key), which kraken2's/kmc's own consumers tolerate but which breaks an invariant sshash's `minimizer_iterator` relies on (a hash tie must imply a value tie) — see that file's comments for the full reasoning. This is sshash-specific packing, not a change to the shared `GreedyMiniTable`/rank semantics itself.

## GreedyMini's two backends

A GreedyMini table can be looked up two ways:

- **Binary-alphabet LUT (default for all three apps)**: a `2^gm_k`-entry table over an AT/CG-collapsed alphabet, with an on-the-fly projection from the real DNA k-mer at lookup time. Compact even for large `gm_k`.
- **DNA-alphabet LUT**: a native `4^gm_k`-entry table, keyed directly by the DNA k-mer — no projection step, but the table is 2^gm_k times larger.

kraken2 only ever uses the binary backend in practice: its `gm_k` values make a `4^gm_k` DNA-alphabet table impractically large. kmc supports both, with binary as the default backend (flipped from an earlier DNA-expanded default); `_dna.gm` tables and `-ode` opt into the DNA-alphabet backend for a binary-alphabet source. sshash's minimizer length `m` has no minimizer-length-specific compile-time cap (the practical bound is `Kmer::max_m=32`, `sshash/include/kmer.hpp`, for the 2-bit DNA alphabet this fork always builds with — `sshash/include/constants.hpp`'s `min_l`/`max_l` are unrelated, they size the skew index's bucket-size partitioning) — a `4^gm_k` DNA-alphabet table stays practical at that scale, but binary is still the default.

## Extend-k: a table built for a shorter length, used at a longer one

A GreedyMini table's own `gm_k` need not equal the length `l` being ranked: when `l > gm_k`, `gm_k` bases decide the order (looked up in the table) and the remaining `l - gm_k` bases become a secondary tie-break folded in below the table's rank (`ext_bits`/`gamma_bits` in `common/minimizer-orders/greedymini_rank.h` and each app's own ranker — `kraken2::GMRanker`, `custom_order_hasher`). This is how kraken2 has always used tables shorter than its own minimizer length, and — since this feature — how sshash's `greedymini` order uses `orders/GreedyMini/w11_k15_canonical.gm` (`gm_k=15`) at `-m 21` (`l=21`, 6 bases of Extend-k).

## mod-minimizer

sshash additionally supports the **mod-minimizer** (Groot Koerkamp & Pibiri, ["The mod-minimizer: a simple and efficient sampling algorithm for long k-mers"](https://doi.org/10.4230/LIPIcs.WABI.2024.11), WABI 2024) — a two-step sampling scheme with provably lower density than a plain (random) minimizer when the minimizer length `m` is large relative to the window `w = k-m+1`. Instead of directly ranking the `w` candidate m-mers, it ranks the more numerous **t-mers** (`t < m`) across the *whole* k-mer with a separate order **O_t** (any of the 4 base orders above), and reduces the winning t-mer's position mod `w` to pick which m-mer is sampled. This is the same algorithm as kraken2's own "Mod" sketching scheme below (both implement Groot Koerkamp & Pibiri's mod-sampling, and — since the incremental-iterator work — even share the `common/minimizer-orders/mod_residue_deque.h` primitive); sshash just exposes it as its own `mod-minimizer-*` orders rather than a separate scheme flag.

The 4 base orders above each get a `mod-minimizer-*` counterpart, reusing that order as O_t:

| Order | O_t | sshash CLI |
|---|---|---|
| mod-minimizer-native-hash | native | `-O mod-minimizer-native-hash --t-length T -s SEED` |
| mod-minimizer-xor-hash | xor-hash | `-O mod-minimizer-xor-hash --t-length T -s SEED` |
| mod-minimizer-pcg-xor-hash | xor-pcg-hash | `-O mod-minimizer-pcg-xor-hash --t-length T -s SEED` |
| mod-minimizer-greedymini | greedymini | `-O mod-minimizer-greedymini --greedy-mini-file FILE.gm --t-length T -s SEED` |

`t` has no single canonical value: the paper's own formula (`t = r + ((k-r) mod w)` for a small lower bound `r`, Definition 16) gives `t=10` for sshash's own tested regime (`k=31, m=21, w=11`) and is what reproduces the paper's published numbers; `t=11` instead matches `orders/GreedyMini/w21_k11_canonical.gm`'s own length exactly (no Extend-k needed for `mod-minimizer-greedymini`), so the benchmark harness (`benchmark/sshash/`) sweeps both for the 3 hash-based orders and fixes `t=11` for `mod-minimizer-greedymini`.

**Canonical (double-stranded) k-mers**: mod-sampling's position-selection rule isn't mirror-equivariant the way the classical per-locus-canonicalized scheme is (see `sshash/include/util.hpp`'s `compute_mod_minimizer_single_strand` doc comment for why, including a worked counterexample), so sshash falls back to its *original* mechanism for this instead — "regular parsing" from the 2022 SSHash paper (also still live in `jermp/sshash-lite`): construction keeps each k-mer's minimizer as whichever single strand the input was walked in, and Lookup "will first search for g and — only if not found — will also search for ĝ [effectively doubling the query time in the worst case]" (2022 paper, "Double Strandedness", quoted verbatim). This applies only to `mod-minimizer-*` orders — the other 3 orders keep sshash's newer single-probe "canonical parsing" scheme unchanged.

## Usage

**kraken2** (`k2 build` / `kraken2-build`, or the `build_db`/`estimate_capacity` binaries directly):

```bash
k2 build --sketching-scheme Minimizer --order RandomXor ...                         # default
k2 build --sketching-scheme Minimizer --order RandomPcg ...
k2 build --sketching-scheme Minimizer --order GreedyMini \
         --greedy-mini-file orders/GreedyMini/w5_k15_canonical.gm [--gm-alphabet binary|dna] ...
k2 build --sketching-scheme Mod --order GreedyMini --greedy-mini-file ... ...        # Mod sketching + any order above
```

**kmc**:

```bash
kmc -orhash -ol<label> [-os<seed>] ...                        # RandomXor (-ol is mandatory with -or)
kmc -orpcg -ol<label> [-os<seed>] ...                         # RandomPcg
kmc -co<file>.gm -ol<label> [-os<seed>] [-ode] ...             # GreedyMini (binary .gm, or _dna.gm/-ode for DNA-alphabet)
kmc ...                                                        # no -or/-co: kmc's own lex order (not shared code)
```

**sshash**:

```bash
./sshash build -i FILE -k 31 -m 21 -O native ...                                  # default, unchanged
./sshash build -i FILE -k 31 -m 21 -O xor-hash -s SEED ...
./sshash build -i FILE -k 31 -m 21 -O xor-pcg-hash -s SEED ...
./sshash build -i FILE -k 31 -m 21 -O greedymini \
    --greedy-mini-file orders/GreedyMini/w11_k15_canonical.gm [--gm-alphabet binary|dna] ...

# mod-minimizer (see "mod-minimizer" above): O_t is any base order + --t-length
./sshash build -i FILE -k 31 -m 21 -O mod-minimizer-xor-hash --t-length 10 -s SEED ...
./sshash build -i FILE -k 31 -m 21 -O mod-minimizer-greedymini \
    --greedy-mini-file orders/GreedyMini/w21_k11_canonical.gm --t-length 11 -s SEED ...
```

### CLI flag mapping

| Order | kraken2 | kmc | sshash |
|---|---|---|---|
| default/lex/native | (n/a — kraken2 has no lex order) | no `-or`/`-co` (kmc's own order) | no `-O` (sshash's own order) |
| RandomXor | `--order RandomXor` (default) | `-orhash -ol<label> [-os<seed>]` | `-O xor-hash -s<seed>` |
| RandomPcg | `--order RandomPcg` | `-orpcg -ol<label> [-os<seed>]` | `-O xor-pcg-hash -s<seed>` |
| GreedyMini | `--order GreedyMini --greedy-mini-file FILE.gm [--gm-alphabet binary\|dna]` | `-co FILE.gm -ol<label> [-os<seed>] [-ode]` | `-O greedymini --greedy-mini-file FILE.gm [--gm-alphabet binary\|dna]` |
| mod-minimizer / Mod sketching scheme | `--sketching-scheme Mod` (combine with any `--order` above) | (n/a — kmc has no analogous two-stage scheme) | `-O mod-minimizer-{native-hash,xor-hash,pcg-xor-hash,greedymini} --t-length T [--greedy-mini-file FILE.gm]` — see "mod-minimizer" above |

## Regenerating a `.gm` table

The [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `.txt` rank-list output is the source of truth; the `.gm` binary is a derived, serialized form:

```bash
cd common/minimizer-orders && make               # builds gm_serializer
./gm_serializer [--canonical] -k N -o OUT.gm IN.txt

# kmc only, to also produce a DNA-alphabet table from a binary one:
cd benchmark/kmc/tools && make                    # builds gm_sigma_extend
./gm_sigma_extend --seed N -o OUT_dna.gm IN.gm
```
