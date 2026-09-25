// Custom minimizer order integration for minimap2 (lab fork).
//
// minimap2 is pure C; the shared minimizer/signature-order library
// (common/minimizer-orders/) is unavoidably C++ (GreedyMiniTable class,
// `enum class GmAlphabet`, templates, in-class member initializers -- none
// of it is includable from a .c file). Per the existing precedent in this
// monorepo ("each app keeps its own thin wrapper around it, since the
// call-site shapes differ" -- orders/README.md), this file is that thin
// wrapper: it #includes the shared headers UNMODIFIED and exposes a flat
// extern-"C" API (order.h) that minimap2's C code calls into. It is the only
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
// A GreedyMini rank (or a plain xor-hash) is a poor choice to reuse directly
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
// GreedyMini is supposed to reduce density. greedymini_rank() below fixes
// this by packing the UHS primary above a PCG-scrambled tiebreak of the full
// k-mer, so ties only occur where the original algorithm intends duplicate
// detection to fire (genuine repeats), not as an artifact of the ranking's
// own coarseness.

#include "order.h"
#include "kvec.h" // kv_resize/kv_push (also pulls in kalloc.h for kmalloc/krealloc)

#include "../common/minimizer-orders/pcg_rxs_m_xs.h"
#include "../common/minimizer-orders/dna_revcomp.h"
#include "../common/minimizer-orders/greedymini_rank.h"
#include "../common/minimizer-orders/mod_residue_deque.h"

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
// mm_order_t
// ---------------------------------------------------------------------------

struct mm_order_s {
	mm_order_kind_t kind;
	uint64_t seed;
	int t;
	std::shared_ptr<const GreedyMiniTable> gm; // GREEDYMINI only
};

namespace {

void set_err(char *errbuf, size_t n, const char *msg) {
	if (!errbuf || n == 0) return;
	std::strncpy(errbuf, msg, n - 1);
	errbuf[n - 1] = '\0';
}

bool is_mod_kind(mm_order_kind_t k) {
	return k == MM_ORDER_MOD_XOR_HASH || k == MM_ORDER_MOD_PCG_XOR_HASH;
}

bool needs_gm_kind(mm_order_kind_t k) {
	return k == MM_ORDER_GREEDYMINI;
}

// Dual-orientation GreedyMini rank (binary alphabet + on-the-fly projection,
// per the user's requirement -- no DNA-alphabet backend, no Extend-k: the
// table's own gm_k must equal the ranked length exactly, checked in
// mm_order_init). Ported in spirit from kraken2::GMRanker (see
// kraken2/src/minimizer_rankers.h) but simplified: minimap2 never persists
// this rank as an index key (see file header), so there is no need for
// GMRanker's Extend-k/gamma tiebreak *packing* -- only for *injectivity*,
// which the PCG-scrambled tiebreak below provides directly from the shared
// pcg_rxs_m_xs.h primitive instead of hand-rolled bit packing.
uint64_t greedymini_rank(uint64_t lmer, int l, const GreedyMiniTable &tbl) {
	const uint64_t atcg = gm_acgt_to_atcg(lmer);
	const uint32_t key_fwd = gm_extract_dna_key_atcg(atcg, l, tbl.gm_k());
	uint16_t best = tbl.primary(key_fwd);
	const uint64_t rc = dna_reverse_complement(lmer, (uint8_t)l);
	if (rc != lmer) {
		const uint64_t atcg_rc = gm_acgt_to_atcg(rc);
		const uint32_t key_rc = gm_extract_dna_key_atcg(atcg_rc, l, tbl.gm_k());
		const uint16_t p_rc = tbl.primary(key_rc);
		if (p_rc < best) best = p_rc;
	}
	// best < kInfRank (0xFFFF) for any UHS member; kInfRank for non-UHS. The
	// tiebreak is unmasked/full-width so every distinct k-mer is (with
	// overwhelming probability) individually ranked, never bucketed together.
	const uint64_t tiebreak = pcg_rxs_m_xs64(lmer) & 0xFFFFFFFFFFFFULL; // low 48 bits
	return ((uint64_t)best << 48) | tiebreak;
}

// Selection rank for the plain single-window path (xor-hash / xor-pcg-hash /
// greedymini). `kmer` is already the lex-min-canonicalized strand (matching
// native's own z = kmer[0]<kmer[1]?0:1), same convention for all three kinds.
uint64_t order_rank(const mm_order_t *ord, uint64_t kmer, int k) {
	switch (ord->kind) {
		case MM_ORDER_XOR_HASH:     return kmer ^ ord->seed;
		case MM_ORDER_XOR_PCG_HASH: return pcg_rxs_m_xs64(kmer ^ ord->seed);
		case MM_ORDER_GREEDYMINI:   return greedymini_rank(kmer, k, *ord->gm);
		default:                    assert(0 && "order_rank: unexpected kind"); return kmer;
	}
}

// Selection rank for a canonical t-mer, used as O_t by the mod-minimizer
// path. GreedyMini's own dual-orientation lookup makes feeding it the
// canonical (already-min'd) t-mer equivalent to feeding it either raw
// strand -- unlike kraken2::ModSampler::rank_tmer, no raw/canonical split is
// needed here (see greedymini_rank's own fwd/rc probe above).
uint64_t order_rank_tmer(const mm_order_t *ord, uint64_t canon_t, int t) {
	switch (ord->kind) {
		case MM_ORDER_MOD_XOR_HASH:     return canon_t ^ ord->seed;
		case MM_ORDER_MOD_PCG_XOR_HASH: return pcg_rxs_m_xs64(canon_t ^ ord->seed);
		default:                        assert(0 && "order_rank_tmer: unexpected kind"); return canon_t;
	}
}

// ---------------------------------------------------------------------------
// HPC span bookkeeping (copied from sketch.c's own private tiny_queue_t --
// that one is `static inline`/file-local, not shared library code, so this
// is minimap2's own helper duplicated the same pragmatic way dna_revcomp.h's
// own doc comment describes for consumers that don't already have a copy in
// scope; not part of common/minimizer-orders/'s order-algorithm surface).
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
// Single-window ordered sketch (xor-hash / xor-pcg-hash / greedymini).
//
// Mechanical port of mm_sketch()'s own robust-winnowing loop (sketch.c),
// with the per-candidate scoring function swapped for order_rank() and the
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

void sketch_single_window(const mm_order_t *ord, void *km, const char *str, int len,
                           int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
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
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers", same as native
			z = kmer[0] < kmer[1] ? 0 : 1;
			++l;
			if (l >= k && kmer_span < 256) {
				info.rank = order_rank(ord, kmer[z], k);
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
// mod-minimizer ordered sketch (Groot Koerkamp & Pibiri, WABI 2024).
//
// minimap2's own (w, k) already match the paper's (and kraken2's/sshash's)
// notation directly: k is the sampled substring length ("m"/"l" elsewhere),
// w is the window's candidate count, and the *paper's* window span is
// exactly minimap2's own w+k-1 bases (matching mm_sketch()'s own
// `l >= w+k-1` first-output condition). t-mers are ranked over that whole
// (w+k-1)-base span, i.e. w+k-t candidate t-mer positions per window --
// `tmer_span` below, matching kraken2::ModSampler's own formula exactly.
//
// Single canonical (min(fwd,rc)) pass, one shared base-offset counter
// (`run_pos`) driving both the t-mer and k-mer rolling windows -- ported in
// spirit from kraken2::ModSampler's dna2-canonical path (mod_minimizer.h),
// reusing ModResidueDeque (common/minimizer-orders/mod_residue_deque.h)
// unmodified for the per-residue-class monotone tracking. Deliberately not
// porting kraken2's amortized O(1) incremental-expire fast path (its
// `win == win_start_+1` special case): every window here re-scans all `w`
// residues' fronts, an O(w) op per output position -- the same complexity
// class as native mm_sketch()'s own rescan branch, not a new source of
// redundant work.
//
// KNOWN LIMITATION, not fixed here: this canonical single-pass approach
// (ranking t-mers via min(fwd,rc), same as kraken2::ModSampler) is NOT
// reverse-complement symmetric, because "reduce the winning t-mer's position
// mod w" does not commute with mirroring a window end-to-end (verified: for
// a window mapped under strand reversal, the two i-offsets sum to a fixed
// constant `tmer_span-1`, but reducing each mod w independently does not
// generally land on mirrored l-mer slots unless w happens to divide that
// constant a particular way). Concretely: building an index from one strand
// and querying with the reverse complement of the *same* physical sequence
// can select largely disjoint sets of minimizers -- measured on a 3kb
// real-genome test, a mod-minimizer-xor-hash mapping between forward and
// reverse-complement copies of the same 3kb sequence shared only 4 of 140
// minimizer keys, versus ~90% overlap for xor-hash/native/greedymini's
// per-locus canonical-parsing designs, which don't have this issue. This is
// the same property sshash's own mod-minimizer implementation works around
// with a dual-probe/"regular parsing" fallback (orders/README.md's
// "mod-minimizer" section, sshash/include/util.hpp's
// compute_mod_minimizer_single_strand doc comment) -- not yet ported here;
// mod-minimizer-* order results should be treated as unreliable for
// cross-strand mapping until it is.
// ---------------------------------------------------------------------------

void sketch_mod_minimizer(const mm_order_t *ord, void *km, const char *str, int len,
                           int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	const int t = ord->t;
	assert(t > 0 && t < k);
	assert(w > 0 && w < 256); // same bound native mm_sketch() itself asserts
	const int tmer_span = w + k - t;
	// Each residue-class deque's `cap` bounds ModResidueDeque's own
	// fixed-size backing array (data[kMaxCap], mod_residue_deque.h).
	// kraken2::ModSampler (which this was ported from) sizes it as
	// tmer_span/w+4 -- but ModResidueDeque::push_back has NO bounds check
	// against `cap`: once size reaches cap, the next push_back silently
	// wraps its write index onto data[head], corrupting front() (the
	// selection result) instead of erroring. A run of push_monotone calls
	// with strictly increasing ranks (nothing to pop) longer than that
	// formula's estimate hits exactly this, and empirically does happen --
	// found via a w=50 benchmark run collapsing to ~71% reads mapped
	// (mod-minimizer-xor-hash on real chrX data) where tmer_span/w+4 = 5
	// was too tight. Since the backing array is a fixed 16 slots regardless
	// of `cap` (kraken2 itself never exercises large-w configurations where
	// its own formula would be this tight), just always use the full
	// kMaxCap instead of a formula that under-provisions here.
	const uint64_t t_mask = (1ULL << (2 * t)) - 1, t_shift1 = (uint64_t)(2 * (t - 1));
	const uint64_t k_mask = (1ULL << (2 * k)) - 1, k_shift1 = (uint64_t)(2 * (k - 1));

	std::vector<ModResidueDeque> deques(w);
	for (int a = 0; a < w; ++a) { deques[a].clear(); deques[a].cap = ModResidueDeque::kMaxCap; }

	struct LEntry { uint64_t kmer; uint32_t span; uint64_t y; int valid; };
	std::vector<LEntry> lring(w);
	for (int a = 0; a < w; ++a) lring[a].valid = 0;

	uint64_t tmer_fwd = 0, tmer_rc = 0, kmer[2] = { 0, 0 };
	ssize_t run_pos = -1; // bases seen since the last reset (0-indexed once >=0)
	int kmer_span = 0;
	ssize_t win_start = -1;
	ssize_t last_emitted = -1; // absolute l-mer ordinal last pushed, so a winner that
	                           // stays the same across consecutive slides (the common
	                           // case) is emitted once, same as native's own dedup
	int lmers_stored = 0;
	tiny_queue_t tq;
	std::memset(&tq, 0, sizeof(tq));

	kv_resize(mm128_t, km, *p, p->n + len / w);

	for (int i = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		if (c >= 4) { // ambiguous base: reset everything, same as native
			for (int a = 0; a < w; ++a) { deques[a].clear(); lring[a].valid = 0; }
			run_pos = -1; kmer_span = 0; win_start = -1; lmers_stored = 0;
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

		const ssize_t ot = run_pos - t + 1;
		if (ot >= 0) {
			const uint64_t canon_t = tmer_fwd < tmer_rc ? tmer_fwd : tmer_rc;
			const uint64_t rank = order_rank_tmer(ord, canon_t, t);
			deques[(size_t)(ot % w)].push_monotone(ModTmerData{ rank, ot });
		}

		const ssize_t ol = run_pos - k + 1;
		if (ol < 0) continue;
		const size_t slot = (size_t)(ol % w);
		if (kmer[0] == kmer[1] || kmer_span >= 256) {
			lring[slot].valid = 0; // palindrome or oversized span: no valid k-mer here
		} else {
			const int z = kmer[0] < kmer[1] ? 0 : 1;
			lring[slot].kmer = kmer[z];
			lring[slot].span = (uint32_t)kmer_span;
			lring[slot].y = (uint64_t)rid << 32 | (uint32_t)i << 1 | (uint64_t)z;
			lring[slot].valid = 1;
		}
		if (lmers_stored < w) ++lmers_stored;

		// The window ending at l-mer ordinal `ol` (its w-th and last member)
		// starts at ordinal `ol - w + 1` -- *that* is the window ready to
		// fire this step, not `ol` itself (which is merely its right edge).
		const ssize_t win = ol - w + 1;
		if (win < 0 || lmers_stored < w || ot < win + tmer_span - 1) continue;

		const int wm = (int)(win % w);
		// Tie-break by the TRUE absolute t-mer position (smallest wins --
		// the paper's "prefer smallest i", i = pos-win), not `(a-wm)%w`:
		// that's only equal to the true offset when pos-win < w. Since
		// tmer_span = w+k-t is always > w (k>t always), a residue class's
		// front can legitimately be more than one window-width ahead, and
		// comparing the mod-w-reduced offset instead of the real one picks
		// the wrong winner on a tie -- and ties are common with xor-hash's
		// weak avalanche. Worse, the wrong comparison depends on `win`
		// itself (via wm), i.e. on *where the current scan started*, not
		// just on window content -- so the same physical repeat could tie-
		// break differently between an index build (scanning a whole
		// contig from position 0) and a query (scanning a read from its
		// own position 0), silently breaking the seed match. Found via a
		// w=50 benchmark run on real chrX data collapsing to ~71% reads
		// mapped; verified fixed against a brute-force reference for
		// several (w,k,t,seed) combinations, including this exact case.
		uint64_t best_rank = UINT64_MAX; int best_a = -1; ssize_t best_pos = -1;
		for (int a = 0; a < w; ++a) {
			ModResidueDeque &dq = deques[(size_t)a];
			while (!dq.empty() && dq.front().pos < win) dq.pop_front();
			if (dq.empty()) continue;
			const uint64_t r = dq.front().rank;
			const ssize_t pos = dq.front().pos;
			if (best_a < 0 || r < best_rank || (r == best_rank && pos < best_pos)) {
				best_rank = r; best_pos = pos; best_a = a;
			}
		}
		if (best_a >= 0) {
			const int best_p = ((best_a - wm) % w + w) % w;
			const ssize_t abs_ord = win + best_p;
			if (abs_ord != last_emitted) {
				const size_t sel = (size_t)best_a;
				if (lring[sel].valid) {
					mm128_t out;
					out.x = (mm_hash64(lring[sel].kmer, k_mask) << 8) | lring[sel].span;
					out.y = lring[sel].y;
					kv_push(mm128_t, km, *p, out);
				}
				last_emitted = abs_ord;
			}
			win_start = win;
		}
	}
	(void)win_start; // retained for readability/parity with the design comment; not re-read
}

} // namespace

// ---------------------------------------------------------------------------
// extern "C" API (order.h)
// ---------------------------------------------------------------------------

extern "C" mm_order_t *mm_order_init(const mm_order_conf_t *conf, char *errbuf, size_t errbuf_len)
{
	if (errbuf && errbuf_len) errbuf[0] = '\0';
	if (!conf) { set_err(errbuf, errbuf_len, "mm_order_init: null conf"); return nullptr; }
	if (conf->kind == MM_ORDER_NATIVE) {
		set_err(errbuf, errbuf_len, "mm_order_init: MM_ORDER_NATIVE does not need an order context");
		return nullptr;
	}
	const bool mod = is_mod_kind(conf->kind);
	if (mod && (conf->t <= 0 || conf->t >= conf->k)) {
		set_err(errbuf, errbuf_len, "mod-minimizer order requires 0 < --mod-t < -k");
		return nullptr;
	}
	if (!mod && conf->t != 0) {
		set_err(errbuf, errbuf_len, "--mod-t given but the chosen order is not a mod-minimizer order");
		return nullptr;
	}
	std::shared_ptr<const GreedyMiniTable> gm;
	if (needs_gm_kind(conf->kind)) {
		if (!conf->greedymini_file || !*conf->greedymini_file) {
			set_err(errbuf, errbuf_len, "greedymini order requires --greedymini-file");
			return nullptr;
		}
		try {
			gm = load_gm_table(conf->greedymini_file, GmAlphabet::BINARY);
		} catch (const std::exception &e) {
			set_err(errbuf, errbuf_len, e.what());
			return nullptr;
		}
		if (!gm) {
			set_err(errbuf, errbuf_len, "failed to load GreedyMini table (missing or corrupt file)");
			return nullptr;
		}
		const int need_gm_k = mod ? conf->t : conf->k;
		if (gm->gm_k() != need_gm_k) {
			char msg[256];
			std::snprintf(msg, sizeof msg,
			               "GreedyMini table gm_k=%d does not match required length %d (%s); "
			               "Extend-k is not supported for minimap2, the table must match exactly",
			               gm->gm_k(), need_gm_k, mod ? "--mod-t" : "-k");
			set_err(errbuf, errbuf_len, msg);
			return nullptr;
		}
	} else if (conf->greedymini_file && *conf->greedymini_file) {
		set_err(errbuf, errbuf_len, "--greedymini-file given but the chosen order does not use a GreedyMini table");
		return nullptr;
	}
	mm_order_t *ord = new (std::nothrow) mm_order_t();
	if (!ord) { set_err(errbuf, errbuf_len, "out of memory"); return nullptr; }
	ord->kind = conf->kind;
	ord->seed = conf->seed;
	ord->t = conf->t;
	ord->gm = gm;
	return ord;
}

extern "C" void mm_order_destroy(mm_order_t *ord) { delete ord; }

extern "C" mm_order_t *mm_order_clone(const mm_order_t *ord)
{
	if (!ord) return nullptr;
	return new (std::nothrow) mm_order_t(*ord); // shares the GreedyMiniTable via shared_ptr copy
}

extern "C" mm_order_kind_t mm_order_kind(const mm_order_t *ord) { return ord ? ord->kind : MM_ORDER_NATIVE; }

extern "C" void mm_sketch_ordered(const mm_order_t *ord, void *km, const char *str, int len,
                                   int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	switch (ord->kind) {
		case MM_ORDER_XOR_HASH:
		case MM_ORDER_XOR_PCG_HASH:
		case MM_ORDER_GREEDYMINI:
			sketch_single_window(ord, km, str, len, w, k, rid, is_hpc, p);
			return;
		case MM_ORDER_MOD_XOR_HASH:
		case MM_ORDER_MOD_PCG_XOR_HASH:
			sketch_mod_minimizer(ord, km, str, len, w, k, rid, is_hpc, p);
			return;
		default:
			assert(0 && "mm_sketch_ordered: MM_ORDER_NATIVE should call mm_sketch() directly");
	}
}

namespace {

void append_bytes(std::vector<uint8_t> &buf, const void *data, size_t n) {
	const uint8_t *b = (const uint8_t *)data;
	buf.insert(buf.end(), b, b + n);
}
template <typename T> void append_pod(std::vector<uint8_t> &buf, const T &v) { append_bytes(buf, &v, sizeof(v)); }

} // namespace

extern "C" int mm_order_serialize(const mm_order_t *ord, uint8_t **out, uint64_t *out_len)
{
	if (!ord || !out || !out_len) return -1;
	std::vector<uint8_t> buf;
	const uint32_t kind32 = (uint32_t)ord->kind;
	append_pod(buf, kind32);
	append_pod(buf, ord->seed);
	const int32_t t32 = (int32_t)ord->t;
	append_pod(buf, t32);
	const uint8_t has_gm = ord->gm ? 1 : 0;
	append_pod(buf, has_gm);
	if (has_gm) {
		const std::vector<uint64_t> raw = gm_table_to_raw(*ord->gm);
		const uint64_t n = (uint64_t)raw.size();
		append_pod(buf, n);
		if (n) append_bytes(buf, raw.data(), raw.size() * sizeof(uint64_t));
	}
	uint8_t *mem = buf.empty() ? nullptr : (uint8_t *)std::malloc(buf.size());
	if (!buf.empty() && !mem) return -1;
	if (!buf.empty()) std::memcpy(mem, buf.data(), buf.size());
	*out = mem;
	*out_len = (uint64_t)buf.size();
	return 0;
}

extern "C" mm_order_t *mm_order_deserialize(const uint8_t *buf, uint64_t len, char *errbuf, size_t errbuf_len)
{
	if (errbuf && errbuf_len) errbuf[0] = '\0';
	size_t off = 0;
	auto need = [&](size_t n) { return off + n <= (size_t)len; };
	if (!need(sizeof(uint32_t) + sizeof(uint64_t) + sizeof(int32_t) + sizeof(uint8_t))) {
		set_err(errbuf, errbuf_len, "truncated order metadata");
		return nullptr;
	}
	uint32_t kind32; std::memcpy(&kind32, buf + off, sizeof(kind32)); off += sizeof(kind32);
	uint64_t seed;   std::memcpy(&seed, buf + off, sizeof(seed));     off += sizeof(seed);
	int32_t t32;     std::memcpy(&t32, buf + off, sizeof(t32));       off += sizeof(t32);
	uint8_t has_gm;  std::memcpy(&has_gm, buf + off, sizeof(has_gm)); off += sizeof(has_gm);

	std::shared_ptr<const GreedyMiniTable> gm;
	if (has_gm) {
		if (!need(sizeof(uint64_t))) { set_err(errbuf, errbuf_len, "truncated GreedyMini table length"); return nullptr; }
		uint64_t n; std::memcpy(&n, buf + off, sizeof(n)); off += sizeof(n);
		if (!need((size_t)n * sizeof(uint64_t))) { set_err(errbuf, errbuf_len, "truncated GreedyMini table bytes"); return nullptr; }
		std::vector<uint64_t> raw(n);
		if (n) std::memcpy(raw.data(), buf + off, (size_t)n * sizeof(uint64_t));
		off += (size_t)n * sizeof(uint64_t);
		try {
			gm = std::make_shared<const GreedyMiniTable>(raw, GmAlphabet::BINARY);
		} catch (const std::exception &e) {
			set_err(errbuf, errbuf_len, e.what());
			return nullptr;
		}
	}
	mm_order_t *ord = new (std::nothrow) mm_order_t();
	if (!ord) { set_err(errbuf, errbuf_len, "out of memory"); return nullptr; }
	ord->kind = (mm_order_kind_t)kind32;
	ord->seed = seed;
	ord->t = (int)t32;
	ord->gm = gm;
	return ord;
}

extern "C" const char *mm_order_kind_name(mm_order_kind_t kind)
{
	switch (kind) {
		case MM_ORDER_NATIVE:            return "native";
		case MM_ORDER_XOR_HASH:          return "xor-hash";
		case MM_ORDER_XOR_PCG_HASH:      return "xor-pcg-hash";
		case MM_ORDER_GREEDYMINI:        return "greedymini";
		case MM_ORDER_MOD_XOR_HASH:      return "mod-minimizer-xor-hash";
		case MM_ORDER_MOD_PCG_XOR_HASH:  return "mod-minimizer-pcg-xor-hash";
	}
	return "unknown";
}

extern "C" int mm_order_parse_kind(const char *s, mm_order_kind_t *out)
{
	if (!s || !out) return -1;
	struct Entry { const char *name; mm_order_kind_t kind; };
	static const Entry table[] = {
		{ "native",                     MM_ORDER_NATIVE },
		{ "xor-hash",                   MM_ORDER_XOR_HASH },
		{ "xor-pcg-hash",               MM_ORDER_XOR_PCG_HASH },
		{ "greedymini",                 MM_ORDER_GREEDYMINI },
		{ "mod-minimizer-xor-hash",     MM_ORDER_MOD_XOR_HASH },
		{ "mod-minimizer-pcg-xor-hash", MM_ORDER_MOD_PCG_XOR_HASH },
	};
	for (const Entry &e : table)
		if (std::strcmp(s, e.name) == 0) { *out = e.kind; return 0; }
	return -1;
}
