# Minimizer schemes and their orders

All four methods choose the selected k-mer of each window (Kraken2's minimizer,
KMC3's signature, SSHash's and Minimap2's minimizer) by an *order* on k-mers,
from the same small set of schemes, and share the C++ implementation of every
non-trivial one instead of each keeping its own copy. Minimap2 is pure C, so it
wraps that shared C++ implementation behind a small extern-"C" bridge — see
[`minimap2/README.fork.md`](../minimap2/README.fork.md) and `minimap2/scheme.cpp`.

**Notation** (as in the paper): k = length of the selected k-mer, w = number of
k-mers per window, L = k + w - 1 (the "L-mer", the method's own k-mer, e.g.
SSHash's `-k 31` with `-m 21` is L = 31, k = 21), t = t-mer length of the
mod-minimizer. Code and flag names map to these terms as follows: the
*projection* (`at0-cg1` / `at1-cg0`) is the polarity of Φ, the DNA -> binary
map whose image is the Sigma plane; `gm_k` is the k of the GreedyMini table;
`mod_t` is t; the code's `l` (Kraken2), `p` (KMC3) and `m` (SSHash) are k; a
super-k-mer (KMC3) is a superwindow.

## Schemes

| Scheme | What it does |
|---|---|
| default / native | Each method's own default order, not shared code: Kraken2 — xor (the canonical k-mer XOR its toggle mask, i.e. the xor order below); KMC3 — its signature order (lexicographically smallest canonical signature, skipping those that start with AAA or ACA, contain AA or end in TG); SSHash — `mixer_64` (multiplicative hash with the FxHash constant ⌊2^64/π⌋ + 1, XOR a seed-derived constant); Minimap2 — `hash64()`, Thomas Wang's 64-bit integer hash. |
| `xor` | A random order: XOR the k-mer with a seed, then rank canonically (with a tie-break for reverse-complement pairs). |
| `xor-pcg-hash` | The xor order's XOR step, then scrambled through the PCG RXS-M-XS 64-bit output permutation (O'Neill 2014) before ranking. |
| `greedymini` | Ranks k-mers by a precomputed lookup table (a GreedyMini order, built by the [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool, which greedily constructs low-density minimizer orders) instead of a formula. The table (`gm_k`, `gm_table`) and the polarity of Φ (`at0-cg1` / `at1-cg0`) are parameters of the scheme, not part of its name. |
| `mod-minimizer-xor`, `mod-minimizer-greedymini` | The canonical mod-minimizer (see below) with xor, respectively a GreedyMini table of t-mers, as the order `O_t` on t-mers; `t` is the extra parameter `mod_t`. |

These names are the only spellings accepted anywhere (CLI flags, run directories, `meta.json`, the `scheme` column of the CSV/XLSX tables).

## Shared implementation

`common/minimizer-schemes/` (how each order and the mod-minimizer rule are implemented, per method: [`common/minimizer-schemes/README.md`](../common/minimizer-schemes/README.md)):

- `pcg_rxs_m_xs.h` — the PCG RXS-M-XS scramble (`pcg_rxs_m_xs64`), used by the `xor-pcg-hash` order of all four methods.
- `greedymini_rank.h` — the GreedyMini table: the bit-manipulation helpers for extracting a DNA k-mer key, the `GreedyMiniTable` LUT class (binary-alphabet tables only; KMC3's DNA-alphabet tables are its own, see below), and the `.gm` binary file I/O (`load_gm_file_raw`/`write_gm_file_raw`/`load_gm_table`). All at global C++ scope, so Kraken2 (via unqualified lookup from inside `namespace kraken2 {}`), KMC3 and SSHash (via unqualified lookup from inside `namespace sshash {}`) use it directly without a wrapper.
- `greedymini_center.h` — on top of `greedymini_rank.h`: the center extension of a canonical table (`gm_center::Order`, compressed packing for Kraken2, injective packing for SSHash and Minimap2), the encoding-agnostic plane extraction (`TwoBitPlanes`) and rolling Sigma helpers, a slow per-base `Reference` and a density helper. Used by Kraken2's `GMRanker`, SSHash's `custom_scheme_hasher` and Minimap2's `scheme.cpp`; KMC3 does not use it.
- `mod_minimizer_rule.h` — the canonical mod-minimizer rule (offset `x mod w` with t == k mod w, the center rule `CenterPick` / `center_tie_offset`, the O(1) window-vs-reverse-complement comparison `le_revcomp_planes`), shared by Kraken2's `ModSampler`, Minimap2's `sketch_mod_minimizer` and SSHash's `compute_mod_minimizer` / `mod_minimizer_iterator`; each method keeps its own scan.
- `dna_revcomp.h` — reverse complement of a packed 2-bit-per-base k-mer. Kept in its own header because Kraken2 has its own `dna_reverse_complement` in `dna_kmer.h`, and one header with both is an ambiguous overload under Kraken2's `using namespace kraken2;`. KMC3 includes this header (SSHash does not: it runs on its own nucleotide encoding, see `sshash/README.md`); Kraken2 uses its own (bit-identical).
- `gm_serializer.cpp` (+ `Makefile`) — builds a `.gm` binary table from a GreedyMini tool's `.txt` rank list, one implementation for all four methods.

`benchmark/kmc/tools/gm_sigma_extend.cpp` (+ `Makefile`) — KMC3-only tool that expands a binary-alphabet `.gm` into a DNA-alphabet `_dna.gm` table, reusing the shared `gm_extract_dna_key`/`dna_reverse_complement` functions.

`orders/GreedyMini/` — the canonical rank tables the four methods load (generated by [GreedyMini](https://github.com/OrensteinLab/GreedyMini), MIT License). Only each table's `.txt` source (the GreedyMini tool's output) is tracked; the `.gm` the methods read is built from it (see [Regenerating a `.gm` table](#regenerating-a-gm-table) and the root [`README.md`](../README.md#getting-started)).

SSHash's `custom_scheme_hasher` (`sshash/include/custom_scheme_hasher.hpp`) additionally folds each GreedyMini key into its packed rank as a secondary tie-break below the table rank. A canonical table assigns the same rank to a mirror pair of keys (and to every key outside the UHS), which Kraken2 and KMC3 tolerate but which breaks an invariant SSHash's `minimizer_iterator` relies on (a hash tie must imply a value tie) — see that file's comments. This is SSHash-specific packing, not a change to the shared `GreedyMiniTable` rank semantics.

## GreedyMini's two backends

A GreedyMini table can be looked up two ways:

- **Binary-alphabet LUT**: a `2^gm_k`-entry table over the {A,T}/{C,G} alphabet, with the DNA k-mer mapped by Φ at lookup time. Compact even for large `gm_k`. All four methods use it, and it is the one the paper uses.
- **DNA-alphabet LUT**: a `4^gm_k`-entry table keyed directly by the DNA k-mer — no projection step, but 2^gm_k times larger. KMC3 only.

Kraken2, SSHash and Minimap2 support only the binary backend (the shared `GreedyMiniTable` is binary-only). KMC3 supports both, with its own code: a binary `.gm` uses the binary backend by default; a `_dna.gm` table (from `gm_sigma_extend`) or `-ode` selects the DNA-alphabet backend.

## The center extension (Kraken2, SSHash, Minimap2): a table built for a shorter length, used at a longer one

KMC3 uses a table only at its own length: the signature length `-p` must equal the table's `gm_k` (`CCustomOrder::check_k`; `-p 9` with the `w20_k9` table). Kraken2, SSHash and Minimap2 may rank a k-mer of length `k > gm_k` with a `gm_k` table: the center `gm_k` bases decide the order (looked up in the table) and the rest only break ties, so `k - gm_k` must be even. Kraken2's GreedyMini order for `k > gm_k`, SSHash's `--scheme greedymini` and Minimap2 (every `-k`, in the injective packing) use this center extension, implemented once in `common/minimizer-schemes/greedymini_center.h` (the order, the Sigma/Gamma plane extraction for each method's base encoding, the rolling Sigma plane, and a per-base reference for tests); SSHash's variant is described below. The GreedyMini order is DNA only; Kraken2 refuses it for protein databases.

## Polarity of Φ (all four methods)

The binary backend maps every base by Φ to {A,T} vs {C,G} (the Sigma plane) and uses the Gamma plane as tie-break. The polarity of Φ is a free choice: `at0-cg1` (the default) maps A,T->0 and C,G->1; `at1-cg0` flips every Sigma bit (Gamma untouched). The two are the same order up to a relabelling of table keys, so on i.i.d. uniform DNA the expected density is the same, but a table is not symmetric under key complement, so the particular density on a real (AT/CG-skewed) dataset can differ. The paper reports `at0-cg1`.

- SSHash: `sshash build --gm-projection {at0-cg1,at1-cg0}` (stored in the index, `sshash query --verbose` prints it), in `sshash/include/custom_scheme_hasher.hpp`.
- Minimap2: `--gm-projection {at0-cg1,at1-cg0}` with `--scheme greedymini` or `mod-minimizer-greedymini` (stored in the `.mmi`), in `minimap2/scheme.cpp` through the shared header's `TwoBitPlanes` flip.
- Kraken2: `k2 build --gm-projection {at0-cg1,at1-cg0}` (`build_db`/`estimate_capacity` `-P`; stored in `opts.k2d`, so `classify` follows the database), for `--order greedymini` with either sketching, in `src/minimizer_rankers.h` (the `TwoBitPlanes` flip in the center extension, a complemented table key when `k == gm_k`).
- KMC3: `-oat1` with `-co` and a binary `.gm` (KMC3 ranks by the binary key and breaks ties on the DNA, so the flip is a one-time permutation `rank'[key] = rank[~key]` of the binary table; a `_dna.gm` table is refused).

`sshash/test/test_gm_center.cpp` checks that `hash()` stays injective on the k-mer and reverse-complement invariant for both polarities; Kraken2's `test_scheme_rankers` checks that the flipped order of x equals the default order of x with the bases relabelled A<->C, G<->T (which complements Sigma and keeps Gamma).

### SSHash: center extension of the table to the k-mer

SSHash's `--scheme greedymini` (binary-alphabet tables only) extends a canonical (mirror-symmetric) table of `gm_k`-mers to `k = gm_k + 2i` (SSHash's `-m`): write the k-mer as `xyz` with `|x| = |z| = i` and `|y| = gm_k`, rank the binary projection of `y` by the table (one lookup; strand invariant since the reverse complement maps `y` to its mirror and `rank[key] == rank[reverse(key)]`), then break ties by the lexicographic minimum of the projection of `xyz` and its reversal, and finally by the Gamma plane (so the hash stays injective on canonical k-mers). It needs `k - gm_k` even (at `k = 21` the odd-`gm_k` tables `w11_k13` and `w11_k15`; an even `gm_k` is rejected), the binary alphabet and a mirror-symmetric table (checked at configure time). On i.i.d. DNA at `k = 21, w = 11` it gives the tables' designed densities (1.592 / 1.524). The SSHash benchmark runs it as `greedymini` once per (table, polarity) variant from `benchmark/sshash/config.yaml` (`w11_k13` and `w11_k15`, `gm_k` = 13 / 15 in the run directory, not in the scheme name); the paper reports `w11_k15`. SSHash's `--scheme mod-minimizer-greedymini` ranks the t-mers with a table of t-mers by the same hash (`t = gm_k + 2i`, in practice `t = gm_k`, so an even `gm_k` such as `w22_k10` is fine there).

Minimap2's `--scheme greedymini` uses the same center order on the k-mer (`-k = gm_k + 2i`, `-k <= 23`), in the injective packing, on Minimap2's ACGT = 0123 encoding (`TwoBitPlanes<true>`), with `--order-seed` toggling the Gamma tie-break as Kraken2's and SSHash's GreedyMini seeds do (it reorders only k-mers with equal rank and Sigma plane: ~1% of the selected positions of a real genome at `-k 19`); the `w10_k15` table at `-k 19` keeps its designed density (1.52 on i.i.d. DNA; ~1.55 under homopolymer compression, `map-pb`'s `-H`, which makes the projected bit sequence non-i.i.d.). Minimap2's `--scheme mod-minimizer-greedymini` ranks the t-mers by the same order at length t (`t = gm_k + 2i`; the benchmark uses the `w20_k9` table at `-k 19 -w 10`, t = 9, so the window holds w + k - t = 20 t-mers, the table's own w).

## Usage

**Kraken2** (`k2 build` / `kraken2-build`, or the `build_db`/`estimate_capacity` binaries directly):

```bash
k2 build --sketching minimizer --order xor ...                                      # default
k2 build --sketching minimizer --order xor-pcg-hash ...
k2 build --sketching minimizer --order greedymini \
         --greedymini-file orders/GreedyMini/w5_k15_canonical.gm [--gm-projection at1-cg0]
k2 build --sketching mod-minimizer --order xor ...                                  # mod-minimizer-xor
k2 build --sketching mod-minimizer --order greedymini --greedymini-file ... ...     # mod-minimizer-greedymini
```

The Kraken2 benchmark's scheme (`-S`) is one of `xor`, `xor-pcg-hash`, `greedymini`,
`mod-minimizer-xor`, `mod-minimizer-greedymini`, i.e. `--sketching` + `--order`;
`-P at0-cg1|at1-cg0` (env `GM_PROJECTION`) sets `--gm-projection` for both
GreedyMini schemes. The table flag differs per method (Kraken2 and Minimap2
`--greedymini-file`, SSHash `--greedy-mini-file`, KMC3 `-co`).

**KMC3**:

```bash
kmc -or xor -ol<label> [-os<seed>] ...                         # xor (-ol is mandatory with -or)
kmc -or xor-pcg-hash -ol<label> [-os<seed>] ...                # xor-pcg-hash
kmc -co<file>.gm -ol<label> [-os<seed>] [-ode] [-oat1] ...     # greedymini (binary .gm; _dna.gm or -ode: DNA-alphabet backend; -oat1: at1-cg0, binary .gm only)
kmc ...                                                        # no -or/-co: KMC3's own signature order (not shared code)
```

The KMC3 benchmark (`kmc3_bench.sh -S default|lex|xor|xor-pcg-hash|greedymini`) adds
`--gm-lut binary|dna` (the table backend, binary by default; same order,
different storage) and `GM_PROJECTION=at1-cg0` (= `-oat1`, needs the binary
backend). `lex` is KMC3's canonical lexicographic baseline, run as `-or xor -os 0`.

**SSHash** (`-k` is L, `-m` is k):

```bash
./sshash build -i FILE -k 31 -m 21 --scheme native ...                                  # default
./sshash build -i FILE -k 31 -m 21 --scheme xor -s SEED ...
./sshash build -i FILE -k 31 -m 21 --scheme xor-pcg-hash -s SEED ...
./sshash build -i FILE -k 31 -m 21 --scheme greedymini \
    --greedy-mini-file orders/GreedyMini/w11_k15_canonical.gm [--gm-projection at1-cg0] ...
./sshash build -i FILE -k 31 -m 21 --scheme mod-minimizer-xor --mod-t T -s SEED ...
./sshash build -i FILE -k 31 -m 21 --scheme mod-minimizer-greedymini \
    --greedy-mini-file orders/GreedyMini/w22_k10_canonical.gm --mod-t 10 ...
```

The SSHash benchmark takes the scheme list with `-S` (`native`, `xor`, ...) and sweeps the
GreedyMini (table, polarity) variants (`--gm-tables`, `--gm-projections`;
default: `config.yaml`'s).

**Minimap2** (`-O` and `-s` are taken by Minimap2, so all scheme flags are long):

```bash
./minimap2 -a ref.fa reads.fq                                                     # native
./minimap2 -a --scheme xor --order-seed SEED ref.fa reads.fq
./minimap2 -a --scheme xor-pcg-hash --order-seed SEED ref.fa reads.fq
./minimap2 -a --scheme greedymini \
    --greedymini-file orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq       # binary alphabet only; center-extended to -k (15, 17, ..., 23)
./minimap2 -a --scheme greedymini --gm-projection at1-cg0 \
    --greedymini-file orders/GreedyMini/w10_k15_canonical.gm ref.fa reads.fq       # at1-cg0 polarity

# mod-minimizer: O_t is xor, or a GreedyMini table of t-mers (t = gm_k + 2i;
# the benchmark uses w20_k9 at -k 19 -w 10, t = 9), plus --mod-t (Minimap2's
# -k/-w are the paper's k/w).
./minimap2 -a --scheme mod-minimizer-xor --order-seed SEED --mod-t T ref.fa reads.fq
./minimap2 -a --scheme mod-minimizer-greedymini --greedymini-file orders/GreedyMini/w20_k9_canonical.gm \
    [--gm-projection at1-cg0] [--order-seed SEED] -k 19 -w 10 --mod-t 9 ref.fa reads.fq
```

### CLI flag mapping

| Scheme | Kraken2 | KMC3 | SSHash | Minimap2 |
|---|---|---|---|---|
| default/native | `--order xor` (default: XOR toggle mask) | no `-or`/`-co` (default: KMC3's signature order) | `--scheme native` / no `--scheme` (default: `mixer_64`, multiplicative hash with the FxHash constant ⌊2^64/π⌋ + 1) | no `--scheme` (default: `hash64()`, Wang's integer hash) |
| xor | `--sketching minimizer --order xor` (default) | `-or xor -ol<label> [-os<seed>]` | `--scheme xor -s<seed>` | `--scheme xor --order-seed SEED` |
| xor-pcg-hash | `--sketching minimizer --order xor-pcg-hash` | `-or xor-pcg-hash -ol<label> [-os<seed>]` | `--scheme xor-pcg-hash -s<seed>` | `--scheme xor-pcg-hash --order-seed SEED` |
| greedymini | `--sketching minimizer --order greedymini --greedymini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` | `-co FILE.gm -ol<label> [-os<seed>] [-ode] [-oat1]` | `--scheme greedymini --greedy-mini-file FILE.gm [--gm-projection at0-cg1\|at1-cg0]` (binary alphabet only) | `--scheme greedymini --greedymini-file FILE.gm [--order-seed SEED] [--gm-projection at0-cg1\|at1-cg0]` (binary alphabet only; center extension to `-k = gm_k + 2i`; the seed toggles the Gamma tie-break) |
| mod-minimizer-xor | `--sketching mod-minimizer --order xor` | (n/a — KMC3 has no mod-minimizer) | `--scheme mod-minimizer-xor --mod-t T -s<seed>` | `--scheme mod-minimizer-xor --order-seed SEED --mod-t T` |
| mod-minimizer-greedymini | `--sketching mod-minimizer --order greedymini --greedymini-file FILE.gm [--gm-projection ...]` | (n/a) | `--scheme mod-minimizer-greedymini --greedy-mini-file FILE.gm --mod-t T [--gm-projection ...]` (a table of t-mers) | `--scheme mod-minimizer-greedymini --greedymini-file FILE.gm --mod-t T [--order-seed SEED] [--gm-projection ...]` (a table of t-mers, t = gm_k + 2i) |

### Run directories and result columns

A scheme is identified by five separate fields, never packed into one string
(`benchmark/common/scheme_names.py`): `scheme` (one of the names above, plus
a method's baseline: `native` for SSHash/Minimap2, `default`/`lex` for KMC3),
`gm_k` (k of the GreedyMini table), `projection` (polarity of Φ: `at0-cg1` | `at1-cg0`), `gm_table`
(table stem, e.g. `w5_k15_canonical`) and `mod_t` (t). They are blank where
they do not apply. KMC3 adds `gm_lut` (`dna` | `binary`, the table backend).
Every benchmark writes them as the columns `scheme, gm_k, projection, gm_table,
mod_t` (KMC3: `scheme, gm_k, projection, gm_table, gm_lut`) of its CSV/XLSX
tables, in `meta.json` and in Kraken2's run manifest.

The run directory of a scheme is `scheme_names.scheme_dir_name(scheme, gm_k,
projection, mod_t)`: `<scheme>[/k<gm_k>_<projection>_t<mod_t>]` (the parts that
apply, in that order), followed by the method's own segments:

| Scheme | Run-dir segment |
|---|---|
| `xor`, `xor-pcg-hash`, `native`, `default`, `lex` | the name itself |
| `greedymini` | `greedymini/k15_at0-cg1`, `greedymini/k15_at1-cg0` |
| `mod-minimizer-xor` | `mod-minimizer-xor/t6` |
| `mod-minimizer-greedymini` | `mod-minimizer-greedymini/k6_at0-cg1_t6` |

KMC3 appends `/lut_<gm_lut>` for greedymini (`greedymini/k9_at1-cg0/lut_binary`,
then `/<seed_idx>`).

The mod-minimizer of Kraken2, Minimap2 and SSHash is a canonical (strand-symmetric) variant of
mod-sampling (Groot Koerkamp & Pibiri, WABI 2024); the scans are
`kraken2/src/mod_minimizer.h`, `minimap2/scheme.cpp` and `sshash/include/util.hpp`
(`compute_mod_minimizer`; SSHash's window is its L-mer, `w = L - k + 1`). With
t == k (mod w) its density is mod-sampling's; any other `t` is rejected (for such t the
candidate offsets of a window and its reverse complement are not mirror images; SSHash, L = 31,
k = 21: a t = 11 GreedyMini table gave ~2.2 against 1.46 at t = 10 with `w22_k10`).

The tie rules (each method's native rule for plain minimizers, the center rule for the
mod-minimizer) are described next to their code, in
[`common/minimizer-schemes/README.md`](../common/minimizer-schemes/README.md).

## Regenerating a `.gm` table

The [GreedyMini](https://github.com/OrensteinLab/GreedyMini) tool's `.txt` rank list is the source of truth; the `.gm` binary is derived from it. From the repository root, `-k` is the `gm_k` in the table's name (`w<w>_k<gm_k>_canonical.txt`):

```bash
make -C common/minimizer-schemes all              # builds gm_serializer
common/minimizer-schemes/gm_serializer --canonical -k 15 \
    -o orders/GreedyMini/w11_k15_canonical.gm orders/GreedyMini/w11_k15_canonical.txt

# all tables at once (as in the root README):
for f in orders/GreedyMini/*.txt; do
  b=$(basename "$f" .txt); k=${b##*_k}; k=${k%%_*}
  common/minimizer-schemes/gm_serializer --canonical -k "$k" -o "${f%.txt}.gm" "$f"
done

# KMC3 only, a DNA-alphabet table from a binary one:
make -C benchmark/kmc/tools                       # builds gm_sigma_extend
benchmark/kmc/tools/gm_sigma_extend --seed N -o OUT_dna.gm IN.gm
```
