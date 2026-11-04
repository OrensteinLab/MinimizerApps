// The greedymini sketch (lab fork): mm_sketch() (sketch.c) with the hash
// replaced by the GreedyMini order, in plain C and in the same loop: the same
// HPC span tracking, palindrome skip, N reset, duplicate flush and >= tie rule,
// so it selects as native does in every respect except which k-mer wins a
// window. scheme.cpp's mm_sketch_scheme() calls it for MM_SCHEME_GREEDYMINI.
//
// The order is the shared center extension in its injective packing
// (common/minimizer-schemes/greedymini_center.h): with k = gm_k + 2i, the
// k-mer x y z (|x| = |z| = i, |y| = gm_k) has the value
//   rank << 2k | min(Sigma_f << k | Gamma_f ^ toggle, Sigma_r << k | Gamma_r ^ toggle)
// where rank = table[Sigma(y)] of the forward strand (0xFFFF for a center
// outside the UHS: such k-mers still compete, ordered by their planes), the
// same for the reverse complement (its center is the mirror of y and the table
// is mirror symmetric), Sigma = the AT/CG projection (A,T -> 0, C,G -> 1;
// complemented for at1-cg0), Gamma = A,C -> 0, G,T -> 1 (complemented on the
// reverse complement strand), first base in the most significant bit, toggle =
// the low k bits of --order-seed. The value is a total order on canonical
// k-mers, so equal values mean equal canonical k-mers (genuine repeats, which
// the duplicate flush is for).
//
// Rank first (as kraken2's scanner): per base only the forward Sigma plane is
// rolled, and a k-mer's rank is one load of the uint16 table. Windows compare
// ranks; the tie-break field (the min(...) above) is computed from the stored
// canonical k-mer only when two candidates' ranks are equal, which keeps the
// comparison of the full value. Equal tie-breaks <=> equal canonical k-mers, so
// equality tests compare k-mers. As for every scheme, the stored x is
// mm_hash64() of the selected k-mer (scheme.cpp's file header), computed only
// when a k-mer is output.
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kvec.h"
#include "mmpriv.h"
#include "scheme.h"

typedef struct { // a simplified version of kdq (sketch.c's tiny_queue_t, file-local there)
	int front, count;
	int a[32];
} gm_tiny_queue_t;

static inline void gm_tq_push(gm_tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int gm_tq_shift(gm_tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

// A candidate k-mer: r = rank (UINT64_MAX = no candidate), ks = canonical
// k-mer << 8 | span, y as native.
typedef struct { uint64_t r, ks, y; } gm_cand_t;

static inline mm128_t gm_out(const gm_cand_t *c, uint64_t mask)
{
	mm128_t m;
	m.x = mm_hash64(c->ks >> 8, mask) << 8 | (c->ks & 0xff);
	m.y = c->y;
	return m;
}

// Even bits of x packed into the low half (bit 2j -> bit j).
static inline uint64_t gm_compress_even(uint64_t x)
{
	x &= 0x5555555555555555ULL;
	x = (x | (x >> 1)) & 0x3333333333333333ULL;
	x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0fULL;
	x = (x | (x >> 4)) & 0x00ff00ff00ff00ffULL;
	x = (x | (x >> 8)) & 0x0000ffff0000ffffULL;
	x = (x | (x >> 16)) & 0x00000000ffffffffULL;
	return x;
}

// Sigma << k | Gamma ^ toggle of one strand of a k-mer (2-bit ACGT word, first
// base highest; Sigma = bit 0 ^ bit 1 of a base, Gamma = bit 1).
static inline uint64_t gm_planes_word(uint64_t x, int k, uint64_t flip_k, uint64_t tb)
{
	return (gm_compress_even(x ^ x >> 1) ^ flip_k) << k | (gm_compress_even(x >> 1) ^ tb);
}

// Reverse complement of a k-mer word.
static inline uint64_t gm_revcomp(uint64_t x, int k)
{
	x = ~x;
	x = (x >> 2 & 0x3333333333333333ULL) | (x & 0x3333333333333333ULL) << 2;
	x = (x >> 4 & 0x0f0f0f0f0f0f0f0fULL) | (x & 0x0f0f0f0f0f0f0f0fULL) << 4;
	x = __builtin_bswap64(x);
	return x >> (64 - 2 * k);
}

// The tie-break field of a k-mer's value (the bits below the rank): the
// smaller of its two strands' (Sigma << k | Gamma ^ toggle).
static inline uint64_t gm_tie(uint64_t kmer, int k, uint64_t flip_k, uint64_t tb)
{
	uint64_t a = gm_planes_word(kmer, k, flip_k, tb), b = gm_planes_word(gm_revcomp(kmer, k), k, flip_k, tb);
	return a < b? a : b;
}

void mm_sketch_gm(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, const mm_gm_params_t *gp, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	const uint64_t key_mask = (1ULL << gp->gm_k) - 1, smask = (1ULL << k) - 1, tb = gp->toggle, flip = gp->flip;
	const uint64_t flip_k = flip? smask : 0;
	const uint16_t *ranks = gp->ranks;
	uint64_t sigma = 0; // forward Sigma plane of the last k bases, first base highest
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	gm_cand_t buf[256], min = { UINT64_MAX, 0, UINT64_MAX };
	gm_tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28) && gp->l == k);
	for (j = 0; j < w; ++j) buf[j] = min;
	memset(&tq, 0, sizeof(gm_tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

// a and b (both candidates, equal ranks) hold the same canonical k-mer
#define GM_SAME(a, b) ((a).r == (b).r && (a).ks >> 8 == (b).ks >> 8)
// a <= b in the order (UINT64_MAX ranks, no candidate, are all equal)
#define GM_LE(a, b) ((a).r < (b).r || ((a).r == (b).r && ((a).r == UINT64_MAX || gm_tie((a).ks >> 8, k, flip_k, tb) <= gm_tie((b).ks >> 8, k, flip_k, tb))))

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		gm_cand_t info = { UINT64_MAX, 0, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				gm_tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= gm_tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			sigma = (sigma << 1 | ((uint64_t)((c ^ c >> 1) & 1) ^ flip)) & smask;
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.r = ranks[(sigma >> gp->i) & key_mask]; // the center y's rank
				info.ks = kmer[z] << 8 | kmer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.r != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (GM_SAME(min, buf[j]) && buf[j].y != min.y) kv_push(mm128_t, km, *p, gm_out(&buf[j], mask));
			for (j = 0; j < buf_pos; ++j)
				if (GM_SAME(min, buf[j]) && buf[j].y != min.y) kv_push(mm128_t, km, *p, gm_out(&buf[j], mask));
		}
		if (GM_LE(info, min)) { // a new minimum; then write the old min
			if (l >= w + k && min.r != UINT64_MAX) kv_push(mm128_t, km, *p, gm_out(&min, mask));
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			uint64_t r = UINT64_MAX, min_tie = 0;
			int n_r = 0;
			if (l >= w + k - 1 && min.r != UINT64_MAX) kv_push(mm128_t, km, *p, gm_out(&min, mask));
			// rank first: the minimal rank r and its count n_r, then, if several
			// candidates have it, the smallest tie-break among them; of equal
			// ones the last in window order wins (native's >=, s.t. min is
			// always the closest k-mer)
			for (j = 0; j < w; ++j) {
				if (buf[j].r < r) r = buf[j].r, n_r = 1;
				else n_r += buf[j].r == r;
			}
			min.r = UINT64_MAX;
#define GM_TAKE(j) do { \
				if (buf[j].r == r) { \
					if (r == UINT64_MAX || n_r == 1) min = buf[j], min_pos = (j); \
					else { \
						uint64_t t = gm_tie(buf[j].ks >> 8, k, flip_k, tb); \
						if (min.r == UINT64_MAX || t <= min_tie) min = buf[j], min_pos = (j), min_tie = t; \
					} \
				} \
			} while (0)
			for (j = buf_pos + 1; j < w; ++j) GM_TAKE(j); // the two loops are necessary when there are identical k-mers
			for (j = 0; j <= buf_pos; ++j) GM_TAKE(j);
#undef GM_TAKE
			if (l >= w + k - 1 && min.r != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (GM_SAME(min, buf[j]) && min.y != buf[j].y) kv_push(mm128_t, km, *p, gm_out(&buf[j], mask));
				for (j = 0; j <= buf_pos; ++j)
					if (GM_SAME(min, buf[j]) && min.y != buf[j].y) kv_push(mm128_t, km, *p, gm_out(&buf[j], mask));
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.r != UINT64_MAX)
		kv_push(mm128_t, km, *p, gm_out(&min, mask));
#undef GM_LE
#undef GM_SAME
}
