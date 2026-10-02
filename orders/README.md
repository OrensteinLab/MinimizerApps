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

## Extend-k: a table built for a shorter length, used at a longer one (kraken2, kmc)

A GreedyMini table's own `gm_k` need not equal the length `l` being ranked: when `l > gm_k`, `gm_k` bases decide the order (looked up in the table) and the remaining `l - gm_k` bases become a secondary tie-break folded in below the table's rank (`ext_bits`/`gamma_bits` in `common/minimizer-orders/greedymini_rank.h` and each app's own ranker — `kraken2::GMRanker`). sshash does not use Extend-k; it uses the centre extension below.

## DNA -> binary projection polarity (sshash `--gm-projection`)

The binary backend classifies every base into {A,T} vs {C,G} (the Sigma plane) and uses the Gamma plane as tie-break. The polarity of Sigma is a free choice: `at0-cg1` (default, the historical behaviour, bit-identical hashes) maps A,T->0 and C,G->1; `at1-cg0` flips every Sigma bit (Gamma untouched). The two are the same order up to a relabelling of table keys, so on i.i.d. uniform DNA the expected density is the same, but a table is not symmetric under key complement, so the *particular* density on a real (AT/CG-skewed) collection can differ. sshash only: `sshash build --gm-projection {at0-cg1,at1-cg0}` (stored in the index, `sshash query --verbose` prints it). Implemented in `sshash/include/custom_order_hasher.hpp`, not in the shared header. `sshash/test/test_gm_centre.cpp` checks that `hash()` stays injective in the m-mer and reverse-complement invariant for both polarities.

### sshash: centre extension of the table to the m-mer

sshash's `-O greedymini` (binary-alphabet tables only) extends an explicit canonical (mirror-symmetric) table of `gm_k`-mers to `m = gm_k + 2i`: write the m-mer as `xyz` with `|x| = |z| = i` and `|y| = gm_k`, rank the **binary projection of `y`** by the table (one lookup; strand invariant since `rc` maps `y` to its mirror and `rank[key] == rank[reverse(key)]`), then break ties by the lexicographic minimum of the projection of `xyz` and its reversal, and finally by the Gamma plane (so the hash stays injective on canonical m-mers). It needs `m - gm_k` even (at `m = 21` only the odd-`gm_k` tables `w11_k13` and `w11_k15`; an even `gm_k` is rejected with a clear error), the binary alphabet and a mirror-symmetric table (checked at configure time). On i.i.d. DNA at `m = 21, w = 11` it gives the tables' designed densities (1.592 / 1.524). The harness names the orders `greedymini-k13`, `greedymini-k15` (+ `-at1`).

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
    --greedy-mini-file orders/GreedyMini/w11_k15_canonical.gm ...
```

### CLI flag mapping

| Order | kraken2 | kmc | sshash |
|---|---|---|---|
| default/lex/native | (n/a — kraken2 has no lex order) | no `-or`/`-co` (kmc's own order) | no `-O` (sshash's own order) |
| RandomXor | `--order RandomXor` (default) | `-orhash -ol<label> [-os<seed>]` | `-O xor-hash -s<seed>` |
| RandomPcg | `--order RandomPcg` | `-orpcg -ol<label> [-os<seed>]` | `-O xor-pcg-hash -s<seed>` |
| GreedyMini | `--order GreedyMini --greedy-mini-file FILE.gm [--gm-alphabet binary\|dna]` | `-co FILE.gm -ol<label> [-os<seed>] [-ode]` | `-O greedymini --greedy-mini-file FILE.gm` (binary alphabet only) |
| Mod sketching scheme | `--sketching-scheme Mod` (combine with any `--order` above) | (n/a — kmc has no analogous two-stage scheme) | (n/a) |

## Regenerating a `.gm` table

The [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `.txt` rank-list output is the source of truth; the `.gm` binary is a derived, serialized form:

```bash
cd common/minimizer-orders && make               # builds gm_serializer
./gm_serializer [--canonical] -k N -o OUT.gm IN.txt

# kmc only, to also produce a DNA-alphabet table from a binary one:
cd benchmark/kmc/tools && make                    # builds gm_sigma_extend
./gm_sigma_extend --seed N -o OUT_dna.gm IN.gm
```
