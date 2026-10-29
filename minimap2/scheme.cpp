// Custom minimizer order integration for minimap2 (lab fork).
//
// minimap2 is pure C; the shared minimizer/signature-order library
// (common/minimizer-schemes/) is unavoidably C++ (GreedyMiniTable class,
// templates, in-class member initializers -- none of it is includable from
// a .c file). Per the existing precedent in this
// monorepo ("each app keeps its own thin wrapper around it, since the
// call-site shapes differ" -- orders/README.md), this file is that thin
// wrapper: it #includes the shared headers UNMODIFIED and exposes a flat
// extern-"C" API (scheme.h) that minimap2's C code calls into. It is the only
// C++ translation unit in minimap2/.
//
// -- Design: order selection vs. index key (read this before touching the
//    sketch functions below) --
//
// In native mm_sketch() (sketch.c), mm_hash64(kmer, mask) plays two roles at
// once: it is the argmin criterion used to pick a window's winning k-mer,
// AND it becomes `x` (packed with the k-mer span into the low 8 bits), the
// literal key stored in the index and compared against query minimizers for
// seed matching -- a match is "same x", minimap2 never re-checks the
// underlying sequence. Bucket routing in index.c also relies on x's bits
// above the low 8 (span) being well-mixed (mm_idx_get: `mi->B[minier&mask]`).
//
// A GreedyMini rank (or a plain xor) is a poor choice to reuse directly
// as that key: GreedyMini collapses every non-UHS k-mer towards one packed
// region, and a plain XOR doesn't redistribute compositional bias the way a
// real avalanche hash does -- either would corrupt seed matching (or badly
// skew index buckets) if stored as `x` directly.
//
// So: **the order decides only which k-mer wins a window; `x` for the
// winner is always derived the same way native does** -- `mm_hash64()` of
// the actual selected/emitted k-mer, computed once per *emitted* minimizer,
// never from the order's own rank value. This costs one extra mm_hash64()
// call per output minimizer (not per candidate scanned), which is cheap
// (amortized ~once every `w` bases) and exactly mirrors how kraken2's own
// GMRanker/ModSampler are used elsewhere in this monorepo: Rank() only
// picks a window's winner, a separate mechanism derives the stored key.
//
// Every rank function below must also be *effectively injective* over the
// space of k-mers/t-mers it ranks: mm_sketch()'s own winnowing loop doesn't
// just tie-break on exact rank equality, it explicitly *flushes every tied
// candidate as an additional minimizer* (its "identical k-mers" duplicate
// logic) -- correct for native's near-injective hash64, but if reused
// naively for GreedyMini's rank (which collapses many distinct non-UHS
// k-mers to the same value) it would pathologically over-select right where
// GreedyMini is supposed to reduce density. RankState<MM_SCHEME_GREEDYMINI> below uses the
// shared centre-extension order in its injective packing (UHS rank, then the
// Sigma and Gamma planes of the k-mer), which is a total order on canonical
// k-mers, so ties only occur where the original algorithm intends duplicate
// detection to fire (genuine repeats), not as an artifact of the ranking's
// own coarseness.

#include "scheme.h"
#include "kvec.h" // kv_resize/kv_push (also pulls in kalloc.h for kmalloc/krealloc)

#include "../common/minimizer-schemes/pcg_rxs_m_xs.h"
#include "../common/minimizer-schemes/greedymini_rank.h"
#include "../common/minimizer-schemes/greedymini_centre.h"
#include "../common/minimizer-schemes/mod_minimizer_rule.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>
#include <sys/types.h> // ssize_t

// From sketch.c: seq_nt4_table is an ordinary (non-static) global there.
extern "C" unsigned char seq_nt4_table[256];

// ---------------------------------------------------------------------------
// mm_scheme_t
// ---------------------------------------------------------------------------

struct mm_scheme_s {
	mm_scheme_kind_t kind;
	uint64_t seed;
	int t;
	std::shared_ptr<const GreedyMiniTable> gm; // GREEDYMINI only
	// GREEDYMINI only, built once at init/deserialize for the index's k
	// (read-only afterwards, so shared by all sketching threads and by
	// mm_scheme_clone's copies): the centre-extension order of the gm_k table
	// at length k (k = gm_k + 2i), and the Sigma/Gamma planes of minimap2's
	// ACGT=0123 k-mers. --gm-projection at1-cg0 complements Sigma (sshash's
	// custom_scheme_hasher.hpp does the same on its own encoding).
	int gm_l = 0;
	mm_gm_projection_t gm_projection = MM_GM_PROJ_AT0_CG1;
	std::shared_ptr<const gm_centre::Order> gm_order;
	gm_centre::TwoBitPlanes<true> gm_planes;
};

namespace {

void set_err(char *errbuf, size_t n, const char *msg) {
	if (!errbuf || n == 0) return;
	std::strncpy(errbuf, msg, n - 1);
	errbuf[n - 1] = '\0';
}

bool is_mod_kind(mm_scheme_kind_t k) {
	return k == MM_SCHEME_MOD_XOR;
}

bool valid_kind(uint32_t k) {
	return k == MM_SCHEME_NATIVE || k == MM_SCHEME_XOR || k == MM_SCHEME_XOR_PCG_HASH ||
	       k == MM_SCHEME_GREEDYMINI || k == MM_SCHEME_MOD_XOR;
}

bool needs_gm_kind(mm_scheme_kind_t k) {
	return k == MM_SCHEME_GREEDYMINI;
}

// GreedyMini rank of a k-mer: the shared centre extension
// (common/minimizer-schemes/greedymini_centre.h), as kraken2's GMRanker and
// sshash's custom_scheme_hasher use it. With k = gm_k + 2i the k-mer is
// x y z (|x| = |z| = i, |y| = gm_k): rank the binary (AT/CG) projection of
// the centre y by the table, then break ties by the smaller of the two
// strands' (Sigma, Gamma ^ toggle), toggle = the low k bits of --order-seed
// (as kraken2's and sshash's greedymini seeds) -- the injective packing, a total
// order on canonical k-mers that is strand invariant (the reverse
// complement's centre is the mirror of y, and the table is mirror
// symmetric). A centre outside the UHS has rank 0xFFFF in the top field, so
// the value is never UINT64_MAX (k <= 23, checked at init) and such k-mers
// still compete, ordered by their planes. i = 0 (k == gm_k) is the plain
// table order.
//
// Per-sketch rank state of the single-window path, one specialisation per
// order kind (sketch_single_window is a template on it, so the order is
// resolved at compile time, not per k-mer): push(c) sees every base of the
// current N-free run, in order, before the palindrome skip; rank(kmer)
// ranks the k-mer ending at the last pushed base. `kmer` is the lex-min
// canonical strand (native's z = kmer[0]<kmer[1]?0:1); the GreedyMini order
// is strand invariant and ignores it.
template <int Kind> struct RankState;

// xor: the seed's bits above the 2k k-mer bits add the same constant to every
// rank, so masking them keeps the order and keeps a rank off UINT64_MAX
// (sketch_single_window's "no candidate"; seed ~0 would otherwise drop
// poly-A).
template <> struct RankState<MM_SCHEME_XOR> {
	uint64_t seed;
	RankState(const mm_scheme_t *ord, int k) : seed(ord->seed & ((1ULL << 2 * k) - 1)) {}
	void push(int) {}
	uint64_t rank(uint64_t kmer) const { return kmer ^ seed; }
};

// xor-pcg-hash: the one input hashing to UINT64_MAX is moved to UINT64_MAX - 1.
template <> struct RankState<MM_SCHEME_XOR_PCG_HASH> {
	uint64_t seed;
	RankState(const mm_scheme_t *ord, int) : seed(ord->seed) {}
	void push(int) {}
	uint64_t rank(uint64_t kmer) const {
		const uint64_t h = pcg_rxs_m_xs64(kmer ^ seed);
		return h - (h == UINT64_MAX);
	}
};

// GreedyMini: the Sigma (AT/CG) and Gamma planes of both strands are rolled
// one base per step, as kraken2's GMRanker and sshash's iterator do with the
// same header's helpers, instead of extracted (two pext) and bit-reversed per
// k-mer; Order::value_of_both_strands then costs one table load. Planes:
// bit p = base p of the strand word, first base = most significant bit. The
// forward word's newest base is its lowest (roll_in_low); the reverse
// complement's newest base is its highest (roll_in_high), and complementing
// a base keeps its AT/CG class (same Sigma bit) and flips Gamma.
template <> struct RankState<MM_SCHEME_GREEDYMINI> {
	const gm_centre::Order &order;
	const gm_centre::TwoBitPlanes<true> &planes;
	int k;
	uint64_t mask, sigma_f = 0, gamma_f = 0, sigma_r = 0, gamma_r = 0;
	RankState(const mm_scheme_t *ord, int k_)
		: order(*ord->gm_order), planes(ord->gm_planes), k(k_), mask((1ULL << k_) - 1) {}
	void push(int c) {
		const uint64_t s = planes.base_sigma((uint64_t)c);
		const uint64_t g = (uint64_t)c >> 1;
		sigma_f = gm_centre::roll_in_low(sigma_f, s, mask);
		gamma_f = ((gamma_f << 1) | g) & mask;
		sigma_r = gm_centre::roll_in_high(sigma_r, s, k);
		gamma_r = (gamma_r >> 1) | ((g ^ 1ULL) << (k - 1));
	}
	uint64_t rank(uint64_t) const { return order.value_of_both_strands(sigma_f, gamma_f, sigma_r, gamma_r); }
};

// ---------------------------------------------------------------------------
// HPC span bookkeeping (copied from sketch.c's own private tiny_queue_t --
// that one is `static inline`/file-local, not shared library code, so this
// is minimap2's own helper duplicated the same pragmatic way common/minimizer-schemes/dna_revcomp.h's
// own doc comment describes for consumers that don't already have a copy in
// scope; not part of common/minimizer-schemes/'s order-algorithm surface).
// ---------------------------------------------------------------------------

struct tiny_queue_t { int front, count; int a[32]; };

inline void tq_push(tiny_queue_t *q, int x) { q->a[((q->count++) + q->front) & 0x1f] = x; }

inline int tq_shift(tiny_queue_t *q) {
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

// ---------------------------------------------------------------------------
// Single-window ordered sketch (xor / xor-pcg-hash / greedymini).
//
// Mechanical port of mm_sketch()'s own robust-winnowing loop (sketch.c),
// with the per-candidate scoring function swapped for RankState::rank() and the
// stored `x` reconstructed from the *selected* k-mer at each push site
// (Cand::kmer) rather than carried through as the comparison value itself
// (Cand::rank) -- see the file header for why. Everything else (HPC span
// tracking, palindrome skipping, the duplicate-tie flush, the tie-break
// rules) is preserved exactly so the ordered path degrades to native's own
// behaviour in every respect except which candidate wins.
// ---------------------------------------------------------------------------

struct Cand { uint64_t rank; uint64_t kmer; uint32_t span; uint64_t y; };

inline mm128_t cand_to_mm128(const Cand &c, uint64_t mask) {
	mm128_t r;
	if (c.rank == UINT64_MAX) { r.x = UINT64_MAX; r.y = UINT64_MAX; return r; }
	r.x = (mm_hash64(c.kmer, mask) << 8) | c.span;
	r.y = c.y;
	return r;
}

template <int Kind>
void sketch_single_window(const mm_scheme_t *ord, void *km, const char *str, int len,
                           int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	RankState<Kind> rs(ord, k);
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	std::vector<Cand> buf(w);
	Cand min = { UINT64_MAX, 0, 0, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28));
	for (i = 0; i < w; ++i) buf[i] = min; // {UINT64_MAX, ...}
	std::memset(&tq, 0, sizeof(tq));
	kv_resize(mm128_t, km, *p, p->n + len / w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		Cand info = { UINT64_MAX, 0, 0, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1;
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | (uint64_t)c) & mask;
			kmer[1] = (kmer[1] >> 2) | ((3ULL ^ (uint64_t)c) << shift1);
			rs.push(c);
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers", same as native
			z = kmer[0] < kmer[1] ? 0 : 1;
			++l;
			if (l >= k && kmer_span < 256) {
				info.rank = rs.rank(kmer[z]);
				info.kmer = kmer[z];
				info.span = (uint32_t)kmer_span;
				info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | (uint64_t)z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info;
		if (l == w + k - 1 && min.rank != UINT64_MAX) {
			for (j = buf_pos + 1; j < w; ++j)
				if (min.rank == buf[j].rank && buf[j].y != min.y) kv_push(mm128_t, km, *p, cand_to_mm128(buf[j], mask));
			for (j = 0; j < buf_pos; ++j)
				if (min.rank == buf[j].rank && buf[j].y != min.y) kv_push(mm128_t, km, *p, cand_to_mm128(buf[j], mask));
		}
		if (info.rank <= min.rank) { // a new minimum; then write the old min
			if (l >= w + k && min.rank != UINT64_MAX) kv_push(mm128_t, km, *p, cand_to_mm128(min, mask));
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.rank != UINT64_MAX) kv_push(mm128_t, km, *p, cand_to_mm128(min, mask));
			min.rank = UINT64_MAX;
			for (j = buf_pos + 1; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.rank >= buf[j].rank) min = buf[j], min_pos = j;
			for (j = 0; j <= buf_pos; ++j)
				if (min.rank >= buf[j].rank) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.rank != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j)
					if (min.rank == buf[j].rank && min.y != buf[j].y) kv_push(mm128_t, km, *p, cand_to_mm128(buf[j], mask));
				for (j = 0; j <= buf_pos; ++j)
					if (min.rank == buf[j].rank && min.y != buf[j].y) kv_push(mm128_t, km, *p, cand_to_mm128(buf[j], mask));
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.rank != UINT64_MAX)
		kv_push(mm128_t, km, *p, cand_to_mm128(min, mask));
}

// ---------------------------------------------------------------------------
// mod-minimizer-xor: strand-symmetric (canonical) mod-sampling with O_t = xor.
//
// mod-sampling (Groot Koerkamp & Pibiri, WABI 2024) in minimap2's notation: k
// is the sampled k-mer length, w the window's k-mer count, the window spans
// L = w + k - 1 bases (mm_sketch()'s own first-output condition) and holds
// T = w + k - t t-mers, ranked as canonical t-mers XOR the seed. Plain
// mod-sampling picks the k-mer at offset x_L mod w (x_L = leftmost minimal
// t-mer); that is not reverse-complement symmetric (reducing the winner's
// position mod w does not commute with mirroring the window), so a read and
// its reverse complement select different loci (~5% at t = 5 on human DNA).
//
// This variant (the rule of common/minimizer-schemes/mod_minimizer_rule.h,
// shared with the kraken2 and sshash forks) requires t == k (mod w), so the
// offset p = x mod w of a minimal t-mer x and the mirror offset w-1-p in the
// reverse complement correspond. A unique minimum samples its offset; a tie
// takes the offset closest to the window centre (w-1)/2, and between two
// equally close ones the left if W <= revcomp(W) lexicographically, the right
// otherwise (the centre rule of Cologni & Pibiri, as in SSHash). The reverse
// complement of a window then samples the mirror k-mer, so both strands select
// the same loci, and the selection never moves backwards (forward).
// validate_mod_minimizer.py checks the code against a brute force of this
// definition.
//
// Streaming: the live t-mers' ranks and residue classes in a power-of-two
// ring indexed by ordinal, with the window minimum and its leftmost and
// rightmost ordinals updated per push by conditional moves; the window is
// rescanned (also branch-free) only when the leftmost minimum leaves it, as
// native mm_sketch() does. (A monotone deque, as kraken2::ModSampler uses,
// mispredicts about one branch per base on random ranks: its pop count is
// data dependent; this costs ~2 rank compares per base and few mispredicts.)
// Residue classes by wrap-around counters. A selected k-mer is marked pending
// in its residue slot and emitted when it leaves the window (no later window
// can select it): output stays in position order, each k-mer once, no sort.
// ---------------------------------------------------------------------------
void sketch_mod_minimizer(const mm_scheme_t *ord, void *km, const char *str, int len,
                          int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	const int t = ord->t;
	assert(t > 0 && t < k);
	assert(w > 0 && w < 256); // same bound native mm_sketch() itself asserts
	const int tmer_span = w + k - t;
	const uint64_t t_mask = (1ULL << (2 * t)) - 1, t_shift1 = (uint64_t)(2 * (t - 1));
	const uint64_t k_mask = (1ULL << (2 * k)) - 1, k_shift1 = (uint64_t)(2 * (k - 1));
	const uint64_t seed = ord->seed;

	// The window's tmer_span t-mers, by ordinal: rank and residue class.
	uint32_t cap = 1;
	while (cap < (uint32_t)tmer_span) cap <<= 1;
	const uint32_t rmask = cap - 1;
	std::vector<uint64_t> trank(cap);
	std::vector<uint8_t> tcls(cap); // ordinal mod w
	uint64_t min_rank = 0; // minimum over the t-mers since the last rescan
	ssize_t x_l = 0, x_r = 0; // its leftmost / rightmost ordinal
	// Equal ranks keep the leftmost x_l and move the rightmost x_r.
	auto take = [&](uint64_t rank, ssize_t x) {
		const bool lt = rank < min_rank, le = rank <= min_rank;
		x_l = lt ? x : x_l;
		x_r = le ? x : x_r;
		min_rank = lt ? rank : min_rank;
	};

	struct LEntry { uint64_t kmer; uint32_t span; uint64_t y; int valid; };
	std::vector<LEntry> lring(w); // indexed by k-mer ordinal mod w

	uint64_t tmer_fwd = 0, tmer_rc = 0, kmer[2] = { 0, 0 };
	ssize_t run_pos = -1;   // bases seen since the last reset (0-indexed once >= 0)
	int t_cls = 0;          // (next t-mer ordinal) mod w
	int l_cls = 0;          // (next k-mer ordinal) mod w
	int kmer_span = 0;
	tiny_queue_t tq;
	std::memset(&tq, 0, sizeof(tq));
	for (int a = 0; a < w; ++a) lring[a].valid = 0;

	// The window's bases as two bit planes, newest base lowest, read only on
	// an exact two-way tie (v <= RC(v)). t == k (mod w) with 0 < t < k gives
	// w <= k - t, so L <= 2k - 2 <= 54 bases fit. Per residue slot: whether
	// its k-mer is selected but not yet emitted.
	const int L = w + k - 1;
	assert(L <= 64);
	std::vector<uint8_t> pending((size_t)w, 0);
	uint64_t win_hi = 0, win_lo = 0;
	auto emit = [&](const LEntry &e) {
		mm128_t out;
		out.x = (mm_hash64(e.kmer, k_mask) << 8) | e.span;
		out.y = e.y;
		kv_push(mm128_t, km, *p, out);
	};
	// Emit the still-pending k-mers of the current run in position order (the
	// oldest live ordinal sits in slot l_cls).
	auto flush_pending = [&]() {
		for (int j = 0, a = l_cls; j < w; ++j, a = a + 1 == w ? 0 : a + 1)
			if (pending[(size_t)a]) { emit(lring[(size_t)a]); pending[(size_t)a] = 0; }
	};

	kv_resize(mm128_t, km, *p, p->n + len / w);

	for (int i = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		if (c >= 4) { // ambiguous base: reset everything, same as native
			flush_pending();
			win_hi = win_lo = 0;
			for (int a = 0; a < w; ++a) lring[a].valid = 0;
			run_pos = -1; t_cls = l_cls = 0; kmer_span = 0;
			tq.count = tq.front = 0;
			continue;
		}
		if (is_hpc) {
			int skip_len = 1;
			if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
				for (skip_len = 2; i + skip_len < len; ++skip_len)
					if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
						break;
				i += skip_len - 1;
			}
			tq_push(&tq, skip_len);
			kmer_span += skip_len;
			if (tq.count > k) kmer_span -= tq_shift(&tq);
		} else {
			kmer_span = kmer_span + 1 < k ? kmer_span + 1 : k;
		}
		++run_pos;

		tmer_fwd = ((tmer_fwd << 2) | (uint64_t)c) & t_mask;
		tmer_rc = (tmer_rc >> 2) | (((uint64_t)(3 ^ c)) << t_shift1);
		kmer[0] = ((kmer[0] << 2) | (uint64_t)c) & k_mask;
		kmer[1] = (kmer[1] >> 2) | (((uint64_t)(3 ^ c)) << k_shift1);
		win_hi = (win_hi << 1) | (uint64_t)(c >> 1);
		win_lo = (win_lo << 1) | (uint64_t)(c & 1);

		const ssize_t ot = run_pos - t + 1; // ordinal of the t-mer ending here
		if (ot >= 0) {
			const uint64_t rank = (tmer_fwd < tmer_rc ? tmer_fwd : tmer_rc) ^ seed;
			trank[(size_t)(ot & rmask)] = rank;
			tcls[(size_t)(ot & rmask)] = (uint8_t)t_cls;
			if (++t_cls == w) t_cls = 0;
			if (ot == 0) min_rank = rank, x_l = x_r = 0; // first t-mer of the run
			else take(rank, ot);
		}

		const ssize_t ol = run_pos - k + 1; // ordinal of the k-mer ending here
		if (ol < 0) continue;
		const int slot = l_cls;
		if (++l_cls == w) l_cls = 0;
		LEntry &le = lring[(size_t)slot];
		if (pending[(size_t)slot]) { // k-mer ol - w leaves every future window
			emit(le);
			pending[(size_t)slot] = 0;
		}
		if (kmer[0] == kmer[1] || kmer_span >= 256) {
			le.valid = 0; // palindrome or oversized span: no valid k-mer here
		} else {
			const int z = kmer[0] < kmer[1] ? 0 : 1;
			le.kmer = kmer[z];
			le.span = (uint32_t)kmer_span;
			le.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | (uint64_t)z;
			le.valid = 1;
		}

		// The window whose last k-mer is `ol` starts at k-mer ordinal win; its
		// t-mers are exactly ordinals [win, ot] (ot == win + tmer_span - 1).
		const ssize_t win = ol - w + 1;
		if (win < 0) continue;

		if (x_l < win) { // the leftmost minimum left the window: rescan [win, ot]
			min_rank = trank[(size_t)(win & rmask)], x_l = x_r = win;
			for (ssize_t x = win + 1; x <= ot; ++x) take(trank[(size_t)(x & rmask)], x);
		}
		// The k-mer of a minimal t-mer x sits in the residue slot of x
		// ((win + (x - win) mod w) mod w == x mod w): a unique minimum is that
		// slot. A tie takes the centre rule over the tied t-mers' offsets
		// p = (tcls[x] - l_cls) mod w (l_cls == win mod w here); the chosen
		// offset p is the k-mer in slot (l_cls + p) mod w.
		int sel_slot = tcls[(size_t)(x_l & rmask)];
		if (x_l != x_r) {
			mod_rule::CentrePick cp;
			for (ssize_t x = x_l; x <= x_r; ++x) {
				if (trank[(size_t)(x & rmask)] != min_rank) continue;
				int p_x = (int)tcls[(size_t)(x & rmask)] - l_cls;
				if (p_x < 0) p_x += w;
				cp.add((uint64_t)p_x, (uint64_t)w);
			}
			const int p_sel = (int)cp.pick([&]() { return mod_rule::le_revcomp_planes(win_hi, win_lo, L); });
			sel_slot = l_cls + p_sel;
			if (sel_slot >= w) sel_slot -= w;
		}
		pending[(size_t)sel_slot] |= (uint8_t)lring[(size_t)sel_slot].valid;
	}
	flush_pending();
}

} // namespace

// ---------------------------------------------------------------------------
// extern "C" API (scheme.h)
// ---------------------------------------------------------------------------

namespace {

// Max k for greedymini: the injective packing takes 16 + 2k bits, and a
// value must stay below UINT64_MAX (sketch_single_window's "no candidate").
const int kGmMaxK = 23;

// Builds ord->gm_order for minimizer length l (ord->gm and ord->gm_projection
// already set). Returns false with a message in errbuf on an invalid length
// (l < gm_k, l - gm_k odd, l > kGmMaxK) or a table the centre extension
// rejects (not mirror symmetric).
bool build_gm_order(mm_scheme_t *ord, int l, char *errbuf, size_t errbuf_len) {
	if (!ord->gm) return true;
	if (l > kGmMaxK) {
		char msg[160];
		std::snprintf(msg, sizeof msg, "greedymini supports -k up to %d (got %d)", kGmMaxK, l);
		set_err(errbuf, errbuf_len, msg);
		return false;
	}
	try {
		// the seed toggles the Gamma tie-break (its low l bits), as kraken2's
		// and sshash's greedymini seeds do; seed 0 is the untoggled order
		ord->gm_order = std::make_shared<const gm_centre::Order>(ord->gm, l, ord->seed, /*injective=*/true);
	} catch (const std::exception &e) {
		set_err(errbuf, errbuf_len, e.what());
		return false;
	}
	ord->gm_l = l;
	ord->gm_planes = gm_centre::TwoBitPlanes<true>(l, ord->gm_projection == MM_GM_PROJ_AT1_CG0);
	return true;
}

} // namespace

extern "C" mm_scheme_t *mm_scheme_init(const mm_scheme_conf_t *conf, char *errbuf, size_t errbuf_len)
{
	if (errbuf && errbuf_len) errbuf[0] = '\0';
	if (!conf) { set_err(errbuf, errbuf_len, "mm_scheme_init: null conf"); return nullptr; }
	if (conf->kind == MM_SCHEME_NATIVE) {
		set_err(errbuf, errbuf_len, "mm_scheme_init: MM_SCHEME_NATIVE does not need a scheme context");
		return nullptr;
	}
	const bool mod = is_mod_kind(conf->kind);
	if (mod && (conf->t <= 0 || conf->t >= conf->k)) {
		set_err(errbuf, errbuf_len, "mod-minimizer scheme requires 0 < --mod-t < -k");
		return nullptr;
	}
	if (mod && !mod_rule::t_is_canonical(conf->k, conf->w, conf->t)) {
		set_err(errbuf, errbuf_len, "mod-minimizer scheme requires --mod-t == -k (mod -w) (the centre tie rule)");
		return nullptr;
	}
	if (!mod && conf->t != 0) {
		set_err(errbuf, errbuf_len, "--mod-t given but the chosen scheme is not a mod-minimizer scheme");
		return nullptr;
	}
	std::shared_ptr<const GreedyMiniTable> gm;
	if (needs_gm_kind(conf->kind)) {
		if (!conf->greedymini_file || !*conf->greedymini_file) {
			set_err(errbuf, errbuf_len, "greedymini scheme requires --greedymini-file");
			return nullptr;
		}
		try {
			gm = load_gm_table(conf->greedymini_file);
		} catch (const std::exception &e) {
			set_err(errbuf, errbuf_len, e.what());
			return nullptr;
		}
		if (!gm) {
			set_err(errbuf, errbuf_len, "failed to load GreedyMini table (missing or corrupt file)");
			return nullptr;
		}
	} else if (conf->greedymini_file && *conf->greedymini_file) {
		set_err(errbuf, errbuf_len, "--greedymini-file given but the chosen scheme does not use a GreedyMini table");
		return nullptr;
	}
	if (conf->gm_projection != MM_GM_PROJ_AT0_CG1 && !needs_gm_kind(conf->kind)) {
		set_err(errbuf, errbuf_len, "--gm-projection only applies to the greedymini scheme");
		return nullptr;
	}
	mm_scheme_t *ord = new (std::nothrow) mm_scheme_t();
	if (!ord) { set_err(errbuf, errbuf_len, "out of memory"); return nullptr; }
	ord->kind = conf->kind;
	ord->seed = conf->seed;
	ord->t = conf->t;
	ord->gm = gm;
	ord->gm_projection = conf->gm_projection;
	if (!build_gm_order(ord, conf->k, errbuf, errbuf_len) || mm_scheme_validate(ord, conf->k, conf->w, errbuf, errbuf_len) != 0) {
		delete ord;
		return nullptr;
	}
	return ord;
}

extern "C" int mm_scheme_validate(const mm_scheme_t *ord, int k, int w, char *errbuf, size_t errbuf_len)
{
	char msg[160];
	if (!ord) { set_err(errbuf, errbuf_len, "null scheme"); return -1; }
	if (k <= 0 || k > 28) { // native mm_sketch()'s own bound
		std::snprintf(msg, sizeof msg, "custom minimizer schemes support -k up to 28 (got %d)", k);
		set_err(errbuf, errbuf_len, msg);
		return -1;
	}
	if (is_mod_kind(ord->kind) && (ord->t <= 0 || ord->t >= k || !mod_rule::t_is_canonical(k, w, ord->t))) {
		std::snprintf(msg, sizeof msg, "mod-minimizer scheme requires 0 < t < k and t == k (mod w) (t = %d, k = %d, w = %d)", ord->t, k, w);
		set_err(errbuf, errbuf_len, msg);
		return -1;
	}
	if (needs_gm_kind(ord->kind) != (ord->gm != nullptr)) {
		set_err(errbuf, errbuf_len, needs_gm_kind(ord->kind) ? "greedymini scheme without a GreedyMini table"
		                                                     : "GreedyMini table given for a scheme that does not use one");
		return -1;
	}
	if (ord->gm && (!ord->gm_order || ord->gm_l != k)) {
		std::snprintf(msg, sizeof msg, "GreedyMini order built for k = %d, index k = %d", ord->gm_l, k);
		set_err(errbuf, errbuf_len, msg);
		return -1;
	}
	return 0;
}

extern "C" int mm_scheme_matches_conf(const mm_scheme_t *ord, const mm_scheme_conf_t *conf)
{
	if (!ord || !conf || ord->kind != conf->kind) return 0;
	if (needs_gm_kind(ord->kind)) {
		if (ord->gm_projection != conf->gm_projection || ord->seed != conf->seed) return 0;
		if (conf->greedymini_file && *conf->greedymini_file) {
			try {
				const std::shared_ptr<const GreedyMiniTable> gm = load_gm_table(conf->greedymini_file);
				if (!gm || gm->size() != ord->gm->size()
				    || std::memcmp(gm->data(), ord->gm->data(), gm->size() * sizeof(uint16_t)) != 0) return 0;
			} catch (const std::exception &) {
				return 0;
			}
		}
		return 1;
	}
	return ord->seed == conf->seed && (!is_mod_kind(ord->kind) || ord->t == conf->t);
}

extern "C" void mm_scheme_destroy(mm_scheme_t *ord) { delete ord; }

extern "C" mm_scheme_t *mm_scheme_clone(const mm_scheme_t *ord)
{
	if (!ord) return nullptr;
	return new (std::nothrow) mm_scheme_t(*ord); // shares the GreedyMiniTable via shared_ptr copy
}

extern "C" mm_scheme_kind_t mm_scheme_kind(const mm_scheme_t *ord) { return ord ? ord->kind : MM_SCHEME_NATIVE; }

extern "C" void mm_sketch_scheme(const mm_scheme_t *ord, void *km, const char *str, int len,
                                   int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	switch (ord->kind) {
		case MM_SCHEME_XOR:
			sketch_single_window<MM_SCHEME_XOR>(ord, km, str, len, w, k, rid, is_hpc, p);
			return;
		case MM_SCHEME_XOR_PCG_HASH:
			sketch_single_window<MM_SCHEME_XOR_PCG_HASH>(ord, km, str, len, w, k, rid, is_hpc, p);
			return;
		case MM_SCHEME_GREEDYMINI:
			sketch_single_window<MM_SCHEME_GREEDYMINI>(ord, km, str, len, w, k, rid, is_hpc, p);
			return;
		case MM_SCHEME_MOD_XOR:
			sketch_mod_minimizer(ord, km, str, len, w, k, rid, is_hpc, p);
			return;
		default:
			assert(0 && "mm_sketch_scheme: MM_SCHEME_NATIVE should call mm_sketch() directly");
	}
}

// .mmi order section: [u64 len][blob], blob = u32 kind, u64 seed, i32 t,
// u8 has_gm, (u64 n, n x u64 table entries in the .gm encoding when has_gm),
// u8 projection, i32 gm_l (the two trailing fields were appended later; a
// blob without them loads as at0-cg1 / a pre-centre-extension index). Both
// directions stream the table in chunks: a full-size temporary buffer, once
// freed, would raise glibc's mmap threshold for the rest of the run and
// inflate the peak memory of an index build or of mapping.

namespace {

const size_t kChunk = 4096; // table entries per fwrite/fread
// Mod-minimizer tie rule stored in an .mmi order section (after gm_l): the
// centre rule; an index without it was built with the retired G/T-strand rule.
constexpr uint8_t kModTieRuleCentre = 1;

uint64_t blob_len(const mm_scheme_t *ord) {
	uint64_t len = sizeof(uint32_t) + sizeof(uint64_t) + sizeof(int32_t) + sizeof(uint8_t);
	if (ord->gm) len += sizeof(uint64_t) + (uint64_t)ord->gm->size() * sizeof(uint64_t);
	return len + sizeof(uint8_t) + sizeof(int32_t) + sizeof(uint8_t);
}

// Reads exactly the section's `left` bytes, field by field.
struct SectionReader {
	FILE *fp;
	uint64_t left;
	bool get(void *dst, size_t n) {
		if (n > left || std::fread(dst, 1, n, fp) != n) return false;
		left -= n;
		return true;
	}
	bool skip_rest() { // unknown trailing fields of a newer writer
		uint8_t tmp[256];
		while (left > 0) if (!get(tmp, left < sizeof tmp ? (size_t)left : sizeof tmp)) return false;
		return true;
	}
};

} // namespace

extern "C" int mm_scheme_dump(const mm_scheme_t *ord, FILE *fp)
{
	if (!ord || !fp) return -1;
	const uint64_t len = blob_len(ord);
	const uint32_t kind32 = (uint32_t)ord->kind;
	const int32_t t32 = (int32_t)ord->t;
	const uint8_t has_gm = ord->gm ? 1 : 0;
	bool ok = std::fwrite(&len, sizeof len, 1, fp) == 1 && std::fwrite(&kind32, sizeof kind32, 1, fp) == 1
		&& std::fwrite(&ord->seed, sizeof ord->seed, 1, fp) == 1 && std::fwrite(&t32, sizeof t32, 1, fp) == 1
		&& std::fwrite(&has_gm, sizeof has_gm, 1, fp) == 1;
	if (ok && has_gm) {
		const uint64_t n = (uint64_t)ord->gm->size(), inf = n;
		const uint16_t *r = ord->gm->data();
		ok = std::fwrite(&n, sizeof n, 1, fp) == 1;
		uint64_t chunk[kChunk];
		for (uint64_t i = 0; ok && i < n;) {
			const size_t c = n - i < kChunk ? (size_t)(n - i) : kChunk;
			for (size_t j = 0; j < c; ++j)
				chunk[j] = r[i + j] == GreedyMiniTable::kInfRank ? inf : r[i + j];
			ok = std::fwrite(chunk, sizeof(uint64_t), c, fp) == c;
			i += c;
		}
	}
	const uint8_t proj = (uint8_t)ord->gm_projection;
	const int32_t gm_l = (int32_t)ord->gm_l;
	const uint8_t mod_rule_id = kModTieRuleCentre;
	ok = ok && std::fwrite(&proj, sizeof proj, 1, fp) == 1 && std::fwrite(&gm_l, sizeof gm_l, 1, fp) == 1
		&& std::fwrite(&mod_rule_id, sizeof mod_rule_id, 1, fp) == 1;
	return ok ? 0 : -1;
}

extern "C" mm_scheme_t *mm_scheme_load(FILE *fp, char *errbuf, size_t errbuf_len)
{
	if (errbuf && errbuf_len) errbuf[0] = '\0';
	uint64_t len = 0;
	if (std::fread(&len, sizeof len, 1, fp) != 1 || len == 0) {
		set_err(errbuf, errbuf_len, "truncated scheme metadata");
		return nullptr;
	}
	SectionReader in{fp, len};
	uint32_t kind32; uint64_t seed; int32_t t32; uint8_t has_gm;
	if (!in.get(&kind32, sizeof kind32) || !in.get(&seed, sizeof seed) || !in.get(&t32, sizeof t32) || !in.get(&has_gm, sizeof has_gm)) {
		set_err(errbuf, errbuf_len, "truncated scheme metadata");
		return nullptr;
	}
	if (!valid_kind(kind32) || kind32 == MM_SCHEME_NATIVE) {
		set_err(errbuf, errbuf_len, "index was built with a minimizer scheme this binary no longer supports; rebuild it");
		return nullptr;
	}
	if (has_gm != (needs_gm_kind((mm_scheme_kind_t)kind32) ? 1 : 0)) {
		set_err(errbuf, errbuf_len, "corrupt scheme metadata (GreedyMini table flag does not match the scheme)");
		return nullptr;
	}

	std::shared_ptr<const GreedyMiniTable> gm;
	if (has_gm) {
		uint64_t n;
		if (!in.get(&n, sizeof n)) { set_err(errbuf, errbuf_len, "truncated GreedyMini table length"); return nullptr; }
		if (n > in.left / sizeof(uint64_t)) { set_err(errbuf, errbuf_len, "truncated GreedyMini table bytes"); return nullptr; }
		try {
			auto table = std::make_shared<GreedyMiniTable>();
			table->reset((size_t)n);
			uint64_t chunk[kChunk];
			for (uint64_t i = 0; i < n;) {
				const size_t c = n - i < kChunk ? (size_t)(n - i) : kChunk;
				if (!in.get(chunk, c * sizeof(uint64_t))) { set_err(errbuf, errbuf_len, "truncated GreedyMini table bytes"); return nullptr; }
				for (size_t j = 0; j < c; ++j) table->set_raw((size_t)(i + j), chunk[j]);
				i += c;
			}
			gm = table;
		} catch (const std::exception &e) {
			set_err(errbuf, errbuf_len, e.what());
			return nullptr;
		}
	}
	uint8_t proj = MM_GM_PROJ_AT0_CG1;
	if (in.left >= sizeof proj) in.get(&proj, sizeof proj);
	if (proj > MM_GM_PROJ_AT1_CG0) { set_err(errbuf, errbuf_len, "invalid GreedyMini projection in index"); return nullptr; }
	int32_t gm_l = 0;
	if (in.left >= sizeof gm_l) in.get(&gm_l, sizeof gm_l);
	uint8_t mod_rule_id = 0; // absent: written before the centre tie rule
	if (in.left >= sizeof mod_rule_id) in.get(&mod_rule_id, sizeof mod_rule_id);
	if (!in.skip_rest()) { set_err(errbuf, errbuf_len, "truncated scheme metadata"); return nullptr; }
	if (is_mod_kind((mm_scheme_kind_t)kind32) && mod_rule_id != kModTieRuleCentre) {
		set_err(errbuf, errbuf_len, "mod-minimizer index built with the retired G/T-strand tie rule; rebuild it");
		return nullptr;
	}
	if (gm && gm_l <= 0) {
		set_err(errbuf, errbuf_len, "greedymini index built before the GreedyMini centre extension; rebuild it");
		return nullptr;
	}
	mm_scheme_t *ord = new (std::nothrow) mm_scheme_t();
	if (!ord) { set_err(errbuf, errbuf_len, "out of memory"); return nullptr; }
	ord->kind = (mm_scheme_kind_t)kind32;
	ord->seed = seed;
	ord->t = (int)t32;
	ord->gm = gm;
	ord->gm_projection = (mm_gm_projection_t)proj;
	if (!build_gm_order(ord, gm_l, errbuf, errbuf_len)) { delete ord; return nullptr; }
	return ord;
}

extern "C" const char *mm_scheme_kind_name(mm_scheme_kind_t kind)
{
	switch (kind) {
		case MM_SCHEME_NATIVE:            return "native";
		case MM_SCHEME_XOR:               return "xor";
		case MM_SCHEME_XOR_PCG_HASH:      return "xor-pcg-hash";
		case MM_SCHEME_GREEDYMINI:        return "greedymini";
		case MM_SCHEME_MOD_XOR:           return "mod-minimizer-xor";
	}
	return "unknown";
}

extern "C" int mm_scheme_parse_kind(const char *s, mm_scheme_kind_t *out)
{
	if (!s || !out) return -1;
	struct Entry { const char *name; mm_scheme_kind_t kind; };
	static const Entry table[] = {
		{ "native",                     MM_SCHEME_NATIVE },
		{ "xor",                        MM_SCHEME_XOR },
		{ "xor-pcg-hash",               MM_SCHEME_XOR_PCG_HASH },
		{ "greedymini",                 MM_SCHEME_GREEDYMINI },
		{ "mod-minimizer-xor",          MM_SCHEME_MOD_XOR },
	};
	for (const Entry &e : table)
		if (std::strcmp(s, e.name) == 0) { *out = e.kind; return 0; }
	return -1;
}
