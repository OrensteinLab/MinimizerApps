/*
 * The strand-symmetric (canonical) mod-minimizer rule, shared by the kraken2,
 * minimap2 and sshash forks. Header only, C++11; each app keeps its own scan
 * (deque or ring) and calls these to decide which k-mer a window samples.
 *
 * Mod-sampling (Groot Koerkamp & Pibiri, WABI 2024): a window of w k-mers
 * (L = w + k - 1 bases) holds T = w + k - t t-mers; with x the locus of its
 * minimal t-mer, plain mod-sampling picks the k-mer at offset x mod w. That
 * is not reverse-complement symmetric. For the canonical variant we use the
 * definition of canonical sampling in simd-minimizers (Groot Koerkamp &
 * Martayan, "SimdMinimizers: Computing random minimizers, fast", bioRxiv
 * 2025.01.27.634998,
 * https://www.biorxiv.org/content/10.1101/2025.01.27.634998v1.full.pdf):
 * every window gets a preferred strand by its G/T count, and the scheme's
 * rule is applied from that strand's side. Applying it to mod-sampling, and
 * the tie rules below (including even L), are these forks' own:
 *
 *  - prefer_forward: the strand with more G/T bases (the reverse complement
 *    has as many G/T as the window has A/C); a tie, possible only for even L,
 *    goes to the forward strand when window <= revcomp(window), so exactly
 *    one of a window and its reverse complement is forward-preferred unless
 *    it is its own reverse complement. kraken2 and minimap2 compare
 *    lexicographically (window_le_revcomp); sshash compares its packed words.
 *  - offset: a forward-preferred window samples x mod w, a reverse-preferred
 *    one applies the same rule from the other end, (w-1) - (((T-1) - x) mod w).
 *    The reverse complement of a window then samples the mirror k-mer.
 *  - Ties at the minimum (several minimal t-mers):
 *      kraken2, minimap2: x = the leftmost minimal t-mer when forward-
 *        preferred, the rightmost otherwise (in a residue-class ring the
 *        rightmost one's slot is mirror_slot of its class);
 *      sshash: centre_tie_offset, the offset closest to the centre (w-1)/2
 *        over all minimal t-mers, which keeps the scheme forward when
 *        t == m mod w (sshash's builder relies on that).
 *  With t = 4 + (m - 4) mod w (m = minimizer length), T-1 == w-1 (mod w) and
 *  both offset rules reduce to x mod w.
 */

#ifndef MINIMIZER_ORDERS_MOD_MINIMIZER_RULE_H_
#define MINIMIZER_ORDERS_MOD_MINIMIZER_RULE_H_

#include <cstdint>

namespace mod_rule {

// Whether the window prefers its forward strand. two_gt_minus_len = 2 * (G/T
// count) - L; tie() (called only when it is 0) breaks the tie.
template <typename Tie>
inline bool prefer_forward(int two_gt_minus_len, Tie tie) {
#if defined(__GNUC__) || defined(__clang__)
  if (__builtin_expect(two_gt_minus_len == 0, 0))
#else
  if (two_gt_minus_len == 0)
#endif
    return tie();
  return two_gt_minus_len > 0;
}

// The same, branch-free, when the tie-break is cheap to evaluate eagerly
// (le_revcomp = window <= revcomp(window)).
inline bool prefer_forward_eager(int two_gt_minus_len, bool le_revcomp) {
  return (two_gt_minus_len > 0) | ((two_gt_minus_len == 0) & le_revcomp);
}

// window <= revcomp(window) lexicographically, for 2-bit bases with
// complement 3 - b (A,C,G,T = 0,1,2,3): base_at(j) is the window's j-th base.
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

// The sampled k-mer's offset in the window for the minimal t-mer at locus x.
inline uint64_t offset(uint64_t x, uint64_t w, uint64_t T, bool forward) {
  return forward ? x % w : (w - 1) - (((T - 1) - x) % w);
}

// Residue-class rings (slot = k-mer position mod w): a reverse-preferred
// window whose rightmost minimal t-mer is in class cls samples the k-mer in
// slot (cls + (w-1) - (T-1)) mod w, independent of the window's position.
inline int mirror_slot(int cls, int w, int T) {
  return (cls + (w - 1) - (T - 1) % w) % w;
}

// Branch-free choice between the forward slot and the mirror slot (the
// preference is close to a fair coin per window).
inline int select_slot(int forward_slot, int mirror, bool forward) {
  const int fmask = -(int)forward;
  return mirror ^ ((forward_slot ^ mirror) & fmask);
}

// sshash's tie rule: over the loci x in [x_l, x_r] with is_min(x), the offset
// closest to the centre (w-1)/2 (doubled distance |2p - (w-1)|), and between
// two equally close ones (mirror images p and w-1-p) the smaller on a
// forward-preferred window, the larger otherwise.
template <typename IsMin>
inline uint64_t centre_tie_offset(uint64_t x_l, uint64_t x_r, uint64_t w,
                                  uint64_t T, bool forward, IsMin is_min) {
  uint64_t best = ~uint64_t(0), lo = 0, hi = 0;
  for (uint64_t x = x_l; x <= x_r; ++x) {
    if (!is_min(x))
      continue;
    const uint64_t p = offset(x, w, T, forward);
    const uint64_t d = 2 * p > w - 1 ? 2 * p - (w - 1) : (w - 1) - 2 * p;
    if (d < best) {
      best = d;
      lo = p;
      hi = p;
    } else if (d == best) {
      lo = p < lo ? p : lo;
      hi = p > hi ? p : hi;
    }
  }
  return forward ? lo : hi;
}

}  // namespace mod_rule

#endif  // MINIMIZER_ORDERS_MOD_MINIMIZER_RULE_H_
