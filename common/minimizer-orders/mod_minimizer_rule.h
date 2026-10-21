/*
 * The strand-symmetric (canonical) mod-minimizer rule, shared by the kraken2,
 * minimap2 and sshash forks. Header only, C++11; each app keeps its own scan
 * (deque or ring) and calls these to decide which k-mer a window samples.
 *
 * Mod-sampling (Groot Koerkamp & Pibiri, WABI 2024): a window v of w k-mers
 * (L = w + k - 1 bases) holds T = w + k - t t-mers; with x the locus of its
 * minimal t-mer, mod-sampling picks the k-mer at offset p = x mod w.
 *
 * We require t == k (mod w) (t_is_canonical; the usual t = r + (k - r) mod w
 * satisfies it). Then T - 1 == w - 1 (mod w): the mirror locus T-1-x of x in
 * RC(v) has offset w-1-p, so with strand-invariant t-mer ranks the minimal
 * t-mers of RC(v) offer exactly the mirror offsets of those of v.
 *
 * Ties at the minimum (several minimal t-mers in different residue classes,
 * which select different k-mers) are broken by the centre rule of Cologni and
 * Pibiri ("Canonical Schemes and Minimizers", Proposition 24; SSHash v6.0.0,
 * PR 93, which applies it to plain minimizers):
 *   1. the offset closest to the window centre (w-1)/2, i.e. the smallest
 *      doubled distance |2p - (w-1)|;
 *   2. between the two equally close mirror offsets p and w-1-p, the left one
 *      if v <= RC(v), the right one otherwise.
 * RC(v) then samples the mirror k-mer (canonical), and since a tied offset
 * stays put while the window slides and new candidates enter on the right,
 * the selection never moves backwards (forward). A unique minimum needs no
 * strand computation at all; v <= RC(v) is evaluated only on an exact
 * two-way tie. The comparison is any fixed total order on windows: kraken2
 * and minimap2 compare lexicographically (le_revcomp_planes,
 * window_le_revcomp), sshash compares its packed words as its plain rule does.
 *
 * This replaces the G/T-strand rule of SimdMinimizers (Groot Koerkamp &
 * Martayan, bioRxiv 2025.01.27.634998: leftmost minimum on the strand with
 * more G/T bases), which is canonical but not forward.
 */

#ifndef MINIMIZER_ORDERS_MOD_MINIMIZER_RULE_H_
#define MINIMIZER_ORDERS_MOD_MINIMIZER_RULE_H_

#include <cstdint>

namespace mod_rule {

// Whether the centre rule applies: t == k (mod w).
inline bool t_is_canonical(int k, int w, int t) {
  return w > 0 && t > 0 && t % w == k % w;
}

// Doubled distance of offset p from the centre (w-1)/2 of w offsets.
inline uint64_t centre_dist2(uint64_t p, uint64_t w) {
  return 2 * p > w - 1 ? 2 * p - (w - 1) : (w - 1) - 2 * p;
}

// The centre rule over candidate offsets added one by one (any order).
struct CentrePick {
  uint64_t best = ~uint64_t(0), lo = 0, hi = 0;

  void add(uint64_t p, uint64_t w) {
    const uint64_t d = centre_dist2(p, w);
    if (d < best) {
      best = d;
      lo = hi = p;
    } else if (d == best) {  // the mirror offset, or p itself again
      lo = p < lo ? p : lo;
      hi = p > hi ? p : hi;
    }
  }

  // le_rc() = v <= RC(v), called only when two offsets remain.
  template <typename LeRc>
  uint64_t pick(LeRc le_rc) const {
    return lo == hi || le_rc() ? lo : hi;
  }
};

// The offset over the loci x in [x_l, x_r] with is_min(x).
template <typename IsMin, typename LeRc>
inline uint64_t centre_tie_offset(uint64_t x_l, uint64_t x_r, uint64_t w,
                                  IsMin is_min, LeRc le_rc) {
  CentrePick c;
  for (uint64_t x = x_l; x <= x_r; ++x)
    if (is_min(x))
      c.add(x % w, w);
  return c.pick(le_rc);
}

// Reverse the low n bits (n in 1..64) of v.
inline uint64_t reverse_low_bits(uint64_t v, int n) {
  v = __builtin_bswap64(v);
  v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
  v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
  v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
  return v >> (64 - n);
}

// window <= revcomp(window) lexicographically, in O(1), for len <= 64 bases
// with 2-bit codes A,C,G,T = 0,1,2,3 (complement 3 - b) held as two bit
// planes: bit i of hi / lo is bit 1 / bit 0 of the base len-1-i, i.e. the
// first base is the highest bit and the newest the lowest (a plane rolled by
// (plane << 1) | bit). Bits at and above len are ignored.
inline bool le_revcomp_planes(uint64_t hi, uint64_t lo, int len) {
  const uint64_t mask = len >= 64 ? ~0ULL : (1ULL << len) - 1;
  hi &= mask;
  lo &= mask;
  // RC base j = 3 - base len-1-j: complement both planes, reverse their order
  const uint64_t rhi = reverse_low_bits(~hi, len);
  const uint64_t rlo = reverse_low_bits(~lo, len);
  const uint64_t diff = ((hi ^ rhi) | (lo ^ rlo)) & mask;
  if (diff == 0)  // its own reverse complement
    return true;
  const int b = 63 - __builtin_clzll(diff);  // the first differing base
  const unsigned a = (unsigned)(((hi >> b) & 1) << 1 | ((lo >> b) & 1));
  const unsigned r = (unsigned)(((rhi >> b) & 1) << 1 | ((rlo >> b) & 1));
  return a < r;
}

// The same for any length, base_at(j) = the window's j-th base (0..3).
template <typename BaseAt>
inline bool window_le_revcomp(int len, BaseAt base_at) {
  for (int j = 0; j < len; ++j) {
    const int a = base_at(j);
    const int b = 3 - base_at(len - 1 - j);
    if (a != b)
      return a < b;
  }
  return true;
}

}  // namespace mod_rule

#endif  // MINIMIZER_ORDERS_MOD_MINIMIZER_RULE_H_
