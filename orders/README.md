# Minimizer/signature orders

Both apps pick a *substring order* — the ranking used to choose which k-mer
substring (kraken2's minimizer, kmc's signature) represents a window — from
the same small set of options, and share the C++ implementation of every
non-trivial one instead of each maintaining its own copy.

## Order modes

| Order | What it does |
|---|---|
| default / lex | Each app's own baseline order (kraken2's original hash-based minimizer order; kmc's canonical lexicographic signature order). Not shared code — this is each app's pre-existing behavior. |
| RandomXor / `hash` | A random order: XOR the substring with a seed, then rank canonically (with a tie-break for revcomp pairs). |
| RandomPcg / `pcg` | RandomXor's XOR step, then scrambled through the PCG RXS-M-XS 64-bit output permutation (O'Neill 2014) before ranking. |
| GreedyMini | Ranks substrings by a precomputed lookup table (a "GreedyMini" order, built by the [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool, which greedily constructs low-density minimizer orders) instead of a formula. |

## Shared implementation

`common/minimizer-orders/`:

- `pcg_rxs_m_xs.h` — the PCG RXS-M-XS scramble (`pcg_rxs_m_xs64`), used by both apps' RandomPcg order.
- `greedymini_rank.h` — the full GreedyMini implementation: the bit-manipulation helpers for extracting a DNA k-mer key, the `GreedyMiniTable` LUT class (both alphabet backends, see below), and the `.gm` binary file I/O (`load_gm_file_raw`/`write_gm_file_raw`/`load_gm_table`). All at global C++ scope, so both kraken2 (via unqualified lookup from inside `namespace kraken2 {}`) and kmc can use it directly without a wrapper.
- `dna_revcomp.h` — reverse-complement of a packed 2-bit-per-base k-mer. Kept in its own header (not `greedymini_rank.h`) because kraken2 already has its own `dna_reverse_complement` in `dna_kmer.h`; combining the two into one header caused an ambiguous-overload compile error under kraken2's `using namespace kraken2;`. kmc includes this header directly; kraken2 keeps using its own (verified bit-identical).
- `gm_serializer.cpp` (+ `Makefile`) — builds a `.gm` binary table from a GreedyMini tool's `.txt` rank-list output. C++, not Python, so both apps share one implementation instead of two independent re-serializers.

`benchmark/kmc/tools/gm_sigma_extend.cpp` (+ `Makefile`) — kmc-only C++ tool that expands a binary-alphabet `.gm` into a DNA-alphabet `_dna.gm` table, reusing the shared `gm_extract_dna_key`/`dna_reverse_complement` functions rather than reimplementing the bit math.

`orders/GreedyMini/` — the canonical rank tables both apps load: each table's `.txt` source (the GreedyMini tool's raw output) alongside its serialized `.gm` (and, where present, `_dna.gm`).

## GreedyMini's two backends

A GreedyMini table can be looked up two ways:

- **Binary-alphabet LUT (default for both apps)**: a `2^gm_k`-entry table over an AT/CG-collapsed alphabet, with an on-the-fly projection from the real DNA k-mer at lookup time. Compact even for large `gm_k`.
- **DNA-alphabet LUT**: a native `4^gm_k`-entry table, keyed directly by the DNA k-mer — no projection step, but the table is 2^gm_k times larger.

kraken2 only ever uses the binary backend in practice: its `gm_k` values make a `4^gm_k` DNA-alphabet table impractically large. kmc supports both, with binary as the default backend (flipped from an earlier DNA-expanded default); `_dna.gm` tables and `-ode` opt into the DNA-alphabet backend for a binary-alphabet source.

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
kmc -or hash -os<seed> ...                                    # RandomXor
kmc -or pcg -os<seed> ...                                     # RandomPcg
kmc -co<file>.gm -ol<label> [-os<seed>] [-ode] ...             # GreedyMini (binary .gm, or _dna.gm/-ode for DNA-alphabet)
kmc ...                                                        # no -or/-co: kmc's own lex order (not shared code)
```

### CLI flag mapping

| Order | kraken2 | kmc |
|---|---|---|
| default/lex | (n/a — kraken2 has no lex order) | no `-or`/`-co` (kmc's own order) |
| RandomXor | `--order RandomXor` (default) | `-or hash -os<seed>` |
| RandomPcg | `--order RandomPcg` | `-or pcg -os<seed>` |
| GreedyMini | `--order GreedyMini --greedy-mini-file FILE.gm [--gm-alphabet binary\|dna]` | `-co FILE.gm -ol<label> [-os<seed>] [-ode]` |
| Mod sketching scheme | `--sketching-scheme Mod` (combine with any `--order` above) | (n/a — kmc has no analogous two-stage scheme) |

## Regenerating a `.gm` table

The [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `.txt` rank-list output is the source of truth; the `.gm` binary is a derived, serialized form:

```bash
cd common/minimizer-orders && make               # builds gm_serializer
./gm_serializer [--canonical] -k N -o OUT.gm IN.txt

# kmc only, to also produce a DNA-alphabet table from a binary one:
cd benchmark/kmc/tools && make                    # builds gm_sigma_extend
./gm_sigma_extend --seed N -o OUT_dna.gm IN.gm
```
