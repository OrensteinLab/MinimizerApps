# Minimizer schemes: how each order and the mod-minimizer rule are implemented

This README documents the implementation of the shared scheme library in this
directory and of each method's consumer of it. For what the schemes are, CLI flags,
the polarity of Φ (the code's projection), run-directory naming and the center-extension rationale,
see [`orders/README.md`](../../orders/README.md); nothing is repeated
here beyond what is needed to follow the code. Scheme names (`xor`,
`xor-pcg-hash`, `greedymini`, `mod-minimizer-xor`, `mod-minimizer-greedymini`)
and parameters (`gm_k`, `projection`, `gm_table`, `mod_t`) are as there.

## Notation

As in the paper (glossary in [`orders/README.md`](../../orders/README.md)): **k** is
the length of the selected k-mer (the code calls it `l` in Kraken2, `m` in SSHash,
`p` in KMC3 and `k` in Minimap2), **w** the number of k-mers per window and
**L = k + w - 1** the window span in bases (the L-mer). `gm_k` is the length of the GreedyMini table's keys
(`gm_k <= k`, `k - gm_k` even), `t` the t-mer length of the mod-minimizers
(`mod_t`). Code identifiers are quoted as written, so `lmer`, `mmer` or `l_`
inside backticks all mean the k-mer.

## 1. Layout and which method uses what

| File | Role | Kraken2 | SSHash | Minimap2 | KMC3 |
|---|---|---|---|---|---|
| `pcg_rxs_m_xs.h` | `pcg_rxs_m_xs64(x)`: `count = x>>59; x ^= x>>(5+count); x *= 12605985483714917081; return x ^ (x>>43)` | `PcgRanker` | `plain_hasher` | `RankState<XOR_PCG_HASH>` | `build_gm_rank_by_dna` |
| `greedymini_rank.h` | key/tie-break extraction (`gm_acgt_to_atcg`, `gm_extract_dna_key/_tb`, `pext`/SWAR), `GreedyMiniTable` (uint16 ranks, `kInfRank = 0xFFFF`), `.gm` I/O (`load_gm_file_raw`, `write_gm_file_raw`, `load_gm_table`) | yes (via `kraken2/src/greedymini_rank.h`) | yes | yes | yes (`CMmer::DnaToBinaryKey`, `load_gm_file_raw`) |
| `greedymini_center.h` | namespace `gm_center`: `Order` (center extension, two packings), `TwoBitPlanes<SigmaXor>`, `roll_in_low/high`, `validate_table`, `Reference`, `density_factor` | `GMRanker` | `custom_scheme_hasher` | `scheme.cpp` (`Order` validates the table and length; the sketches roll their own planes) | not used |
| `mod_minimizer_rule.h` | namespace `mod_rule`: `t_is_canonical`, `center_dist2`, `CenterPick`, `center_tie_offset`, `le_revcomp_planes`, `window_le_revcomp` (tested by `test_mod_minimizer_rule.cpp`, `make test`) | `ModSampler` | `util.hpp` | `sketch_mod_minimizer` | not used |
| `dna_revcomp.h` | `dna_reverse_complement` on packed 2-bit k-mers | own copy in `dna_kmer.h` | own encoding | own loops | `CMmer::DnaReverseComplement` |
| `gm_serializer.cpp` (+ `Makefile`) | builds a `.gm` from a rank-list `.txt`; `--canonical` gives each listed k-mer and its reverse the same rank | tool | tool | tool | tool |

`.gm` format: `[size_t N = 2^gm_k][N x uint64]`, rank `0..n-1` for the `n` UHS keys and
the sentinel `N` for non-UHS keys (`greedymini_rank.h` header comment).
`GreedyMiniTable::set_raw` maps that to a dense uint16 table with `kInfRank`
(`0xFFFF`) for "infinity" and rejects ranks that do not fit.

## 2. Per-order implementation

### Shared definitions

- Planes. `TwoBitPlanes<SigmaXor>` extracts an n-bit Sigma and Gamma plane from
  a word of 2-bit bases, first base in the MSB, by `pext` (BMI2) or the SWAR
  fallback `compress_even_bits`. `SigmaXor = true` (Kraken2, KMC3; Minimap2's
  own rolled planes follow it too: ACGT = 0..3) uses Sigma = bit0 ^ bit1 (A,T -> 0; C,G -> 1), Gamma = bit1;
  `SigmaXor = false` (SSHash: A=00 C=01 G=11 T=10) uses Sigma = bit0,
  Gamma = bit1. The `flip` argument complements Sigma (`at1-cg0`); Gamma is
  untouched.
- Rolling. `roll_in_low` (newest base lowest, Kraken2; Minimap2's forward strand
  rolls the same way in C) and `roll_in_high` (newest base highest, SSHash;
  Minimap2's rc strand for `mod-minimizer-greedymini`)
  slide a Sigma plane by one base instead of a `pext` per window.
- Mirror-symmetric canonical table: `rank[key] == rank[reverse_bits(key)]`
  (`mirror_key`); `validate_table` checks it and, when `require_rank_pairs`,
  that every finite rank is held by one key and its mirror only (SSHash and the
  injective packing skip the pair check).

### `xor`

Rank = `canonical ^ seed`, canonical = numerically smaller of the k-mer and its
reverse complement. There is no shared code beyond `dna_revcomp`.

- Kraken2: `RandomRanker::Rank` returns `canonical_lmer ^ toggle_` (code name `lmer` = the k-mer)
  (`minimizer_rankers.h`); the default scan `NextMinimizerXor` in
  `mmscanner.cc` inlines it (`canonical_representation(lmer_) [& spaced mask]`
  then `^ toggle_mask_`). Equal candidates keep the leftmost (the pop test is
  strict `>`).
- SSHash: `plain_hasher<XOR_HASH>::hash(v) = v ^ seed` on the numerically
  canonical k-mer (`canonical_mmer_at`, `util.hpp`). Ties are broken by the
  center-closest rule of `util::compute_minimizer` (`resolve_tie`, and
  `break_tie` in `minimizer_iterator`), mirror images by `kmer <= rc(kmer)`.
- Minimap2: `RankState<MM_SCHEME_XOR>::rank = kmer ^ seed` with `seed` masked to
  2k bits so a rank never equals `UINT64_MAX` (the "no candidate" value of
  `sketch_single_window`); `kmer` is the lex-min strand `kmer[z]`.
- KMC3: no per-window formula; `install_random` (`custom_order.cpp`) fills
  `g_gm_rank_by_dna[i] = i ^ seed` for all `4^k` k-mers and runs
  `build_gm_dna_luts` (below), so lookup is a DNA LUT.

### `xor-pcg-hash`

The same as `xor` with `pcg_rxs_m_xs64(canonical ^ seed)` as the rank:
`PcgRanker` (Kraken2), `plain_hasher` (SSHash), `RankState<MM_SCHEME_XOR_PCG_HASH>`
(Minimap2, which moves the one input hashing to `UINT64_MAX` to `UINT64_MAX-1`),
`build_gm_rank_by_dna` with `Mode::PCG` (KMC3). Kraken2's `ModSampler` tabulates
the t-mer rank in `rank_lut_` when `t <= kMaxLutT = 8`.

### `greedymini`

Table: `2^gm_k` uint16 ranks over the Sigma (AT/CG) projection of a `gm_k`-mer,
`kInfRank` for non-UHS keys.

For `k == gm_k` (all methods except where a longer k is configured) the table ranks
the k-mer itself. For `k = gm_k + 2i` the center extension of
`gm_center::Order` is used (Kraken2, SSHash, Minimap2; not KMC3, see 4). The
k-mer is `x y z` with `|x| = |z| = i`, `|y| = gm_k`; the rank is one load
`table[Sigma(y)]` (`center_key`, `rank_of_sigma`, `rank_at`), strand invariant
because the reverse complement's center is `mirror(y)` and the table is
mirror-symmetric. `check_length` requires `k - gm_k` even.

Packing of the order value (`Order::pack`, constructor flag `injective`):

| | compressed (default; Kraken2) | injective (SSHash, Minimap2) |
|---|---|---|
| value | `rank << sh_rank \| Sx \| ybit \| Sz \| (Gamma^toggle)` | `rank_or_inf << 2k \| (Sigma << k \| Gamma^toggle)` |
| width | `rank_bits + 2i + 1 + k` | `16 + 2k`, so `k <= 24` (`k <= 23` in Minimap2, `kGmMaxK`; SSHash `16+2k <= 63`) |
| infinite center | `UINT64_MAX` (all such windows tie) | rank field `0xFFFF`, still ordered by the planes |
| `ybit` | `key > mirror_key(key)`, precomputed in `ybit_` | not present (full Sigma is in the value) |
| `validate_table` pair check | yes | skipped |

Compressed vs injective, in short. Both are monotone encodings of the same
idea (order by table rank first, then break ties by the strand frame with the
smaller Sigma plane, then by the DNA) into one 64-bit integer, so that comparing
integers compares orders. They differ in what is kept:

- Compressed keeps the minimum of bits: since a canonical table gives each
  finite rank to exactly one key and its mirror, `(rank, ybit)` already
  determines the center's Sigma, so only the flanks `Sx`, `Sz` and Gamma are
  appended. This fits `k` up to 31 (Kraken2), but an infinite center (not in the
  UHS) collapses to `UINT64_MAX`, so all such windows tie and the order is not
  total; the consumer must break those ties itself (Kraken2 uses the
  strand-invariant xor key). It also requires the table's rank-pair property
  (`validate_table` pair check).
- Injective keeps the whole Sigma and Gamma planes of the winning strand below a
  fixed 16-bit rank field, so the value identifies the canonical k-mer: hash ties
  imply equal k-mers (what SSHash's iterators rely on), infinite windows are still
  ordered, and no rank-pair property is needed. The price is `16 + 2k` bits, so
  `k <= 24` (23 in Minimap2, 23 in SSHash by its `<= 63` check).

The order of a strand is that strand's packed value; the order of the k-mer is
the min over the two strands (`value_of_planes`: one lookup, rc planes by
`reverse_low_bits`, rc Gamma = `reverse ^ (2^k-1)`; `value_of_both_strands`:
planes of both strands supplied by a rolling caller; `value_of_two_strands[_masked]`:
strands supplied separately, used for spaced seeds). The toggle (seed) is XORed
into Gamma only. `Reference` is the slow per-base specification used by tests.

Projection `at0-cg1`/`at1-cg0`: Sigma polarity. Applied as the `flip` flag of
`TwoBitPlanes` in Kraken2 (`GMRanker::flip_`, also `key_flip_` for the
`k == gm_k` path), SSHash (`derive`: `TwoBitPlanes<false>(k, projection==AT1_CG0)`)
and Minimap2 (`build_gm_order`: `mm_gm_params_t::flip`, XORed into each base's Sigma bit). KMC3 has no planes: `install_table` permutes the
table once, `flipped[key] = ranks[key ^ key_mask]`.

Per method:

- Kraken2 (`minimizer_rankers.h`, `GMRanker`). `k == gm_k`: `probe_strand` takes
  the key by `pext` of the ATCG word for each strand (fwd and rc);
  `pick_better_primary` keeps the finite/smaller primary, ties by
  `(primary << k) | (Gamma ^ toggle)` (`pack`). `k > gm_k`: `init_center` builds
  `gm_center::Order(table, k, toggle)` (compressed) and
  `TwoBitPlanes<true>(k, flip)`; `rank_first()` is true for revcom version 1.
  The scanner compares only the table rank (`rank_of_sigma`, with Sigma either
  `pext`ed per k-mer or rolled by `roll_in_low` when `roll_sigma_`, default off:
  `kDefaultRollSigma = false`) and calls `value()` only when two finite ranks
  tie (`PopDominatedRankFirst`). Spaced seed: the rank lookup sees the unmasked
  center; each strand's word is ANDed with the mask independently (never
  `rc(x & mask)`) and only the tie-break fields x, z, Gamma come from the masked
  planes (`planes()`, `value_of_two_strands_masked`); with `k == gm_k` the table
  key is taken from the masked word. Because the compressed value is not
  injective, equal finite values and all-infinite windows are broken by the
  strand-invariant `dna_hash_key ^ toggle` (`KeyKeepsBack`,
  `InfiniteWindowMinimizer`).
- SSHash (`custom_scheme_hasher.hpp`). `configure` loads the `.gm`, validates
  (`k - gm_k` even, `16 + 2k <= 63`, mirror symmetry via
  `validate_table(..., false)`), and `derive` builds
  `gm_center::Order(table, k, seed, injective=true)`. `hash_gm(x) =
  hash_from_planes(planes_of_mmer(x))` (code name `mmer` = the k-mer): two `pext`s and one 64-bit bit reversal
  give the rc planes, then `value_of_both_strands`. The hash is strand
  invariant (`strand_invariant()`), so callers hash the raw forward k-mer and
  canonicalize only the winner's value. The injective packing means the Sigma
  plane is folded into the value below the rank, so hash ties are value ties.
  The table is stored in the `.sshash` file (`m_gm_ranks`).
- Minimap2 (`sketch_gm.c`, `mm_sketch_gm`; parameters from `scheme.cpp`'s
  `build_gm_order`, which builds an `Order` to validate the table and length).
  Injective packing, seed toggles Gamma, plain C (no `gm_center` code at sketch
  time). Rank first, like Kraken2's scanner: per base one forward Sigma bit is
  rolled (newest base lowest, as `roll_in_low`), the rank is one `uint16_t`
  table load of the center, candidates compare ranks, and only on equal ranks
  is the tie-break (min over both strands of `Sigma << k | Gamma ^ toggle`)
  computed from the stored canonical k-mers, so the comparisons are those of
  the full injective value. Equal values are equal canonical k-mers, which
  matters because the loop (`mm_sketch()`'s) flushes every tied candidate;
  `UINT64_MAX` is "no candidate". `mod-minimizer-greedymini`
  (`TmerRank<MM_SCHEME_MOD_GREEDYMINI>`, `scheme.cpp`) computes the full value
  of each t-mer from the two strands' `Sigma << t | Gamma` words, rolled one base
  per step (forward newest base lowest, rc newest base highest, rc Gamma `g ^ 1`),
  and one table load; never `UINT64_MAX` (t <= 23).
- KMC3: see 4.

### `mod-minimizer-xor` / `mod-minimizer-greedymini`

Mod-sampling: a window of `L = k + w - 1` bases holds `w` k-mers and
`T = L - t + 1 = w + k - t` t-mers (code names differ: Kraken2 `k, l` = L, k;
SSHash `k, m`; Minimap2 `w, k`); the minimal t-mer at locus `x` selects the
k-mer at offset `x mod w`.
Canonical (strand-symmetric) rule (`mod_minimizer_rule.h`), with t == k (mod w) (`t_is_canonical`):

- a unique minimal t-mer at `x` samples offset `x % w`;
- a tie at the minimum: the center rule, `CenterPick` (`add` each tied offset, `pick(le_rc)`)
  or `center_tie_offset(x_l, x_r, w, is_min, le_rc)`: the offset closest to `(w-1)/2`
  (doubled distance `center_dist2`), between two equally close ones the left if `le_rc()`;
- `le_rc` = v <= RC(v), evaluated only on an exact two-way tie: `le_revcomp_planes(hi, lo, len)`
  (lexicographic, O(1) on two bit planes, `len <= 64`; Kraken2, Minimap2) or
  `!(kmer_rc < kmer)` on packed words (SSHash); `window_le_revcomp` is the per-base version
  (tests).

**Tie rule, mod-minimizer** (`common/minimizer-schemes/mod_minimizer_rule.h`, shared by
Kraken2 `ModSampler`, Minimap2 `sketch_mod_minimizer` and SSHash `compute_mod_minimizer` /
`mod_minimizer_iterator`). Here the tie rule changes the output in every method: the window
minimum is over t-mers, and the k-mer sampled for the minimal t-mer at x is the one at offset
p = x mod w, so equal-rank t-mers at positions with different x mod w select different k-mers. All
methods require t == k (mod w) (Kraken2 `t = 4 + (l-4) mod w`; Minimap2 `--mod-t` and SSHash
`--mod-t` are checked), so the mirror of x in RC(v) lies in class w-1-p and the candidate
offsets are mirror-symmetric. The center rule is applied to them: the candidate closest to the
window center (w-1)/2, and between the two equally close mirror candidates p and w-1-p, the
left if v <= RC(v), else the right. RC(v) then samples the mirror k-mer (canonical), and the
chosen offset never moves backwards as the window slides (forward). The center rule does extra
work only on a tied minimum, and compares v with RC(v) only on an exact two-way tie. The native plain rules do not carry over:
leftmost (Kraken2) is not canonical for the mod-minimizer, keeping every tied k-mer (Minimap2)
inflates density at small t, and KMC3 has no mod-minimizer. v <= RC(v) is lexicographic in
Kraken2 and Minimap2 (`le_revcomp_planes`, O(1) on two bit planes) and SSHash's packed-word
comparison `!(kmer_rc < kmer)`, the one its plain rule uses; both are total orders,
which is all canonicity needs.

**Tie rules, plain minimizers.** Each method keeps its native rule for breaking ties
between equal-rank positions, because each rule serves how that method uses the
minimizer. In every order we use, the rank is a total order on canonical k-mers (a hash of
kappa(x) = min(x, rc(x)), or GreedyMini's injective packing; in Kraken2 the packing is not
injective on its own and the strand-invariant key `dna_hash_key ^ toggle` completes it), so
equal ranks occur only when the same canonical k-mer appears more than once in a window.
Kraken2 keeps the leftmost tied position (its monotone deque); it uses only the minimizer's
value, as the hash-table key whose LCA it looks up, so the rule does not change its output.
KMC3 moves the signature to the newest equal occurrence: a superwindow ends when its
signature's position leaves the k-mer, so keeping the newest occurrence keeps the signature
longest, giving the longest superwindows and the least bin I/O; bins are chosen by the
canonical value, so counts do not depend on the tie position (at the level of values this is
robust winnowing, Schleimer, Wilkerson & Aiken, SIGMOD 2003, Definition 3; the KMC 2 and
KMC 3 papers do not discuss ties). Minimap2 emits every tied position (standard winnowing;
Winnowmap, Jain et al. 2020, replaces it by robust winnowing) and skips k-mers equal to their
reverse complement, which is strand-symmetric without a tie rule; it is also the upstream
baseline. SSHash stores superwindows by minimizer and offset, so the position matters: the
scheme must be canonical (a k-mer and its reverse complement find the same minimizer
occurrence) and forward (its builder checks it), which no leftmost/rightmost-by-strand rule
achieves; it uses the center rule of Cologni and Pibiri (SSHash v6.0.0, PR 93): the tied
position closest to the center of the k-mer, and between two equally close ones the left if
v <= RC(v), else the right. GreedyMini changes only the rank, so every method applies its native
rule to it.

How `t` is set: Kraken2 fixes `t = 4 + (k-4) % w` (`ModSampler` constructor, `r = 4`,
needs `2 <= w <= 16`, and a GreedyMini table must have `gm_k == t`); SSHash and
Minimap2 take `--mod-t` (`mod_t`). For `t = r + (k-r) % w`, `(T-1) = (w-1) mod w` and
both offset rules reduce to `x mod w`. `mod-minimizer-greedymini` ranks t-mers
with the `greedymini` machinery at length t (Kraken2: `GMRanker` with `k = t`;
SSHash: `custom_scheme_hasher` with `m_l = mod_t`, `gm_k <= t` allowed by center
extension; `test_mod_minimizer.cpp` tests each table at `t = gm_k` and `t = gm_k + 2`; Minimap2: `TmerRank<MM_SCHEME_MOD_GREEDYMINI>`, the `greedymini` order at `t = gm_k + 2i`, used in the benchmark with the `w20_k9` table at `-k 19 -w 10 --mod-t 9`).

Sliding structures (read from the code, they differ):

| Method | Structure | Update | Where |
|---|---|---|---|
| Kraken2 | one monotone deque of t-mers, power-of-two ring `dq_[64]` (`kDqCap`), entries `{rank, pos, cls}`, equal ranks kept; plus a k-mer ring `lmer_ring_[w]` indexed by residue class | `push_tmer` pops strictly larger back; `finish_window` drops front `pos < win`; a single front entry is the selection, a front run of equal ranks goes to the center rule (`tie_slot`) | `ModSampler` in `mod_minimizer.h`; driven by `NextMinimizerMod` in `mmscanner.cc` |
| SSHash | ring of the T t-mer hashes (`m_ring[64]`, slot `(m_base + x) & 63`), tracked `m_min_hash`, `m_left`, `m_right`; for GreedyMini a ring of 16-bit ranks plus rolled Sigma plane (`m_rank`, `m_sigma`) | hash only the new t-mer; full re-scan only when the leftmost minimum leaves (`m_left == 0`); GM hashes in full only on rank ties | `mod_minimizer_iterator` in `minimizer_iterator.hpp`; definition `util::compute_mod_minimizer` |
| Minimap2 | not a deque: ring of t-mer ranks and residue classes (`trank`, `tcls`, power-of-two), running `min_rank`, `x_l`, `x_r` updated by conditional moves (`take`); k-mer ring `lring[w]` with a per-slot `pending` flag; the window's two bit planes for `le_revcomp_planes` (`L <= 54`) | rescan `[win, ot]` only when `x_l < win`; `x_l == x_r` selects `tcls[x_l]`, otherwise the center rule over `[x_l, x_r]`; the selected k-mer is marked pending in its slot and emitted when it leaves the window, so output stays in position order | `sketch_mod_minimizer` in `scheme.cpp` (comment explains why not a deque) |

Kraken2 additionally replaces per-t-mer ranking by a LUT `rank_lut_[4^t]` for
`xor-pcg-hash`/GM when `t <= 8`, on the raw forward t-mer for GM and the
canonical t-mer for PCG; xor needs none (`select_kind`: `K_XOR`, `K_LUT_RAW`,
`K_LUT_CAN`, `K_RANK_RAW`, `K_RANK_CAN`).

### Sliding-window minimum for the plain (non-mod) orders

| Method | Machinery | Tie handling |
|---|---|---|
| Kraken2 | `std::deque` of `{rank, candidate, pos}` (`queue_`; the default xor path uses a 16-byte `xq_`, `NextMinimizerXor`). New k-mer pops the back while worse (`PopDominatedRankFirst`/`PopDominatedGM`, plain `>` for xor/pcg); front expires when `pos < queue_pos_ - (L - k)` (code: `k` = L, `l` = k) | strict `>` keeps the leftmost on equal ranks; GM equal finite values: `KeyKeepsBack` (hash key ^ toggle) |
| SSHash | no deque. `minimizer_iterator<Kmer, Mode>`: keeps the current min, re-scans all `w` loci (`rescan`) only when it leaves the window (`m_min_pos_in_kmer == 0`), otherwise one new hash and a compare; counts ties in `m_num_mins`. `strand_invariant_minimizer_iterator` (GreedyMini): ring of 16-bit ranks (`m_ring[136]`: 64 slots, every rank stored twice so a window's ranks are contiguous; SSE min 8 at a time), rolled Sigma plane (`roll_sigma`), full hash only on rank ties | `break_tie`: tied locus closest to the center, mirror images by `kmer <= rc(kmer)` |
| Minimap2 | native `mm_sketch()` structure ported in `sketch_single_window<Kind>` (xor, pcg): a circular `buf[w]` of `Cand{rank, kmer, span, y}`, `min`/`min_pos`; new min on `info.rank <= min.rank`, rescan of `buf` when the old min leaves; HPC span via `tiny_queue_t`. GM: the same loop in C, `mm_sketch_gm` (`sketch_gm.c`), `buf[w]` of `{rank, k-mer << 8 \| span, y}`, compared rank first (tie-break from the k-mer only on equal ranks; the rescan finds the minimal rank and its count first) | ties flush as extra minimizers when `y` differs (native behaviour kept); stored `x` is always `mm_hash64(selected k-mer)`, not the rank (`cand_to_mm128`, `gm_out`) |
| KMC3 | no queue. `CSplitter` (`kmc_core/splitter.cpp`, `CalcStatsT`/`ProcessReadsT`, templated on `CMmer::LutMode()`) holds `current_signature` and `end_mmer`: each base `end_mmer.insert_sym<Mode>`; `<` replaces, `==` also moves the start, and when the current signature leaves the k-mer (`signature_start_pos + kmer_len - 1 < i`) the k-mer is rescanned with `<=` (rightmost on ties) | equal `current_val` moves to the newer position; for GM the value is `(rank << 32) \| min(str^seed, rc^seed)` |

## 3. Per-method summary

| Method | Order evaluated on | Alphabet / encoding | Table backend | Sliding structure | Code |
|---|---|---|---|---|---|
| Kraken2 | canonical k-mer (xor/pcg); raw k-mer with strand handling inside `GMRanker` (greedymini); canonical t-mers (mod xor/pcg), raw t-mers (mod GM) | 2-bit, ACGT = 0..3, `lookup_table_` | binary only (`GreedyMiniTable`); center extension, compressed packing | deque (minimizer); monotone ring deque (mod) | `kraken2/src/{mmscanner.cc,minimizer_rankers.h,mod_minimizer.h}` |
| SSHash | canonical k-mer (xor/pcg), raw k-mer with strand-invariant hash (GM), t-mers for mod | 2-bit, A=00 C=01 G=11 T=10 (traditional encoding not supported for GM) | binary only; center extension, injective packing, table embedded in the index | re-scan iterator; rank ring + rolled Sigma (GM); hash ring (mod) | `sshash/include/{custom_scheme_hasher.hpp,minimizer_iterator.hpp,util.hpp}` |
| Minimap2 | lex-min strand k-mer, k <= 28 (GM: k <= 23) | 2-bit, ACGT = 0..3 (`seq_nt4_table`), pure C caller via extern "C" | binary only; center extension, injective packing; table stored in `.mmi` | circular `buf[w]` with rescan (GM: rank first); ring of t-mer ranks (mod) | `minimap2/scheme.cpp`, `sketch_gm.c`, `scheme.h` |
| KMC3 | signature (a k-mer, k = 5..11) | 2-bit, ACGT = 0..3 | DNA LUT (`4^k`) for xor/pcg and `_dna.gm`; binary LUT (`2^k`, default for `.gm`) | `current_signature`/`end_mmer` rescan | `kmc/kmc_api/{custom_order.cpp,mmer.h,mmer.cpp}`, `kmc/kmc_core/splitter.cpp` |

## 4. KMC3's two backends and why k == gm_k

`CCustomOrder::LutBackend` (`custom_order.h`): `BinaryRuntime` (default) and
`DnaExpanded`. `CMmer::LutMode()` selects one of three templated paths
(`LUT_DEFAULT`, `LUT_DNA`, `LUT_BINARY`).

- DNA LUT (`SetGmDnaLut`): `g_gm_dna_sort[4^k]` is the sort position of the
  canonical k-mer, `g_gm_dna_canonical[4^k]` the canonical k-mer, built by
  `build_gm_dna_luts`. `gm_pick_canonical_dna` picks the smaller rank of `i`
  and its rc, ties by `i ^ seed` vs `rc ^ seed`. `current_val =
  gm_dna_sort_lut[str]`, one load. Built from a binary `.gm` through
  `build_gm_dna_luts_from_binary` (`g_gm_rank_by_dna[i] = ranks[DnaToBinaryKey(i,p)]`),
  loaded from a prebuilt `_dna.gm` (magic `KMC_DNA_GM_MAGIC`, `gm_sigma_extend`
  makes these), and used for `xor`/`xor-pcg-hash` too.
- Binary LUT (`SetGmBinaryLut`): a `2^k` uint32 rank table; `insert_sym<LUT_BINARY>`
  rolls `bin_key` (one Sigma bit per base via `DnaBaseToBinaryBit`) and `str_rc`,
  then `pack_binary`: `current_val = rank << 32 | min(str ^ seed, str_rc ^ seed)`.
  `get()` returns `current_val ^ seed`, the canonical k-mer used as the bin key.
  `-oat1` is applied by permuting the table (above), not at lookup.

Why `k == gm_k` (KMC3's `-p`): `CCustomOrder::check_k` throws unless the table's `gm_k` equals the
signature length `k`, and requires `5 <= k <= 11`. In the code, the binary key
is the whole k-mer (`DnaToBinaryKey(str, p)` = `gm_extract_dna_key(dna, p, p)`,
`bin_key_mask = 2^p - 1`, code name `p` = k) and the DNA tables have exactly `4^k` entries
(`load_dna_gm_file` checks `n == 1 << 2k`); there is no center-key extraction or
Sigma/Gamma tie-break beyond `str ^ seed`.

Measured cost of the backends: see `benchmark/kmc/tools/bench_signature_impls.cpp`.

## 5. Tests and adding an order

Tests that pin the implementations:

| Subject | File |
|---|---|
| Kraken2 GreedyMini ranker / center extension / rank-first vs `Rank()` | `kraken2/src/test_scheme_rankers.cc` |
| Kraken2 `ModSampler` and `MinimizerScanner` (mod) | `kraken2/src/test_mod_minimizer.cc` |
| SSHash center hash vs `hash_reference`, Python vectors | `sshash/test/test_gm_center.cpp`, `sshash/test/gm_center_vectors.inc` |
| SSHash rolling planes, `strand_invariant_minimizer_iterator` | `sshash/test/test_gm_rolling.cpp` |
| SSHash minimizer iterators / schemes | `sshash/test/test_minimizer.cpp`, `sshash/test/test_streaming_scheme.cpp` |
| SSHash mod-minimizer vs string-level reference | `sshash/test/test_mod_minimizer.cpp` |
| Minimap2 schemes end to end | `minimap2/test_scheme.sh` |
| benchmark-level scheme sanity | `benchmark/{sshash,minimap2,kmc}/tests/test_scheme_sanity.py` (KMC3: `test_custom_order_sanity.py`), `benchmark/common/tests/test_scheme_names.py` |

Adding a new order:

1. Put the order's pure math (rank from a word or from Sigma/Gamma planes) in a
   header here, with a slow per-base reference next to it (as `Reference`).
2. Add the name to each method's parser (`sketching_schemes.h` in Kraken2,
   `custom_scheme_mode` in `custom_scheme_hasher.hpp`, `mm_scheme_kind_t` /
   `mm_scheme_parse_kind` in Minimap2, `CCustomOrder::ParseMode` in KMC3) and to
   `benchmark/common/scheme_names.py`.
3. Wire the method consumer: Kraken2 `make_order_ranker`; SSHash `dispatch_scheme`
   and `plain_hasher`/`hash_gm`; Minimap2 `RankState<Kind>` plus the
   `mm_sketch_scheme` dispatch; KMC3 `Install`.
4. Make it strand invariant (a read and its rc must select the same minimizer)
   and injective on canonical k-mers, or define a strand-invariant tie-break.
5. Add the tests above (reference equality, rc symmetry, resets at `N`), and
   record any new flag in `orders/README.md`.
