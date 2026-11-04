/*
 * GreedyMini CENTER extension of a canonical binary table (additive header;
 * includes greedymini_rank.h and changes nothing in it).
 *
 * A canonical (mirror-symmetric) binary table ranks the 2^gm_k binary strings
 * of length gm_k. To use it for l-mers with l = gm_k + 2i, view an l-mer as
 * x y z with |x| = |z| = i and |y| = gm_k, and let Sigma (AT/CG = 0/1 or the
 * reverse polarity, the caller's choice) and Gamma (the residual bit that
 * tells A from T and C from G) be its two l-bit planes, first base in the most
 * significant bit.
 *
 *   rank(strand)   = table[Sigma(y)]
 *   packed(strand) = rank | Sx | ybit | Sz | (Gamma ^ toggle)        (MSB..LSB)
 *   ybit           = [Sigma(y) > mirror(Sigma(y))]
 *   order(l-mer)   = min(packed(fwd strand), packed(reverse-complement strand))
 *
 * (an infinite rank is UINT64_MAX and therefore never the minimum unless both
 * strands are infinite). Why this is the right order:
 *  - The reverse complement of x y z has the mirror of Sigma(y) at its center
 *    and a mirror-symmetric table gives the same rank: ONE lookup, strand
 *    invariant by construction.
 *  - Every finite rank is shared by exactly a key and its mirror, so
 *    (rank, ybit) determine Sigma(y) and the packed value is lexicographic in
 *    (rank, Sigma(xyz), Gamma ^ toggle) while needing only
 *    rank_bits + 2i + 1 + l bits. Ties on rank are thus broken by the strand
 *    frame with the smaller Sigma, then by the DNA.
 *  - The winning frame's (Sigma, Gamma) determine the DNA, so the order is
 *    injective on canonical l-mers (hence a total order on them).
 *  - It keeps the density the table was designed for for every l = gm_k + 2i.
 *
 * Encoding-agnostic: everything operates on the Sigma/Gamma planes; the caller
 * extracts them from its own k-mer encoding. For an exact reverse complement
 * Sigma(rc) = reverse(Sigma) and Gamma(rc) = reverse(Gamma) ^ (2^l - 1)
 * (complementing flips A<->T and C<->G, i.e. the Gamma bit, and keeps Sigma);
 * value_of_planes() uses this. With a spaced seed (or any per-strand
 * transform) pass both strands' planes to value_of_two_strands() instead.
 *
 * Two packings (constructor flag `injective`):
 *  - compressed (default; kraken2, l up to 31): the layout above, an infinite
 *    center gives UINT64_MAX. Windows with an infinite center therefore all
 *    tie (the order is not injective on them); a consumer that needs a total
 *    order there must break ties itself (kraken2's scanner breaks equal
 *    values, and a window of infinite l-mers, by the strand-invariant
 *    xor key, so its selection stays strand symmetric).
 *  - injective (sshash; needs 16 + 2l <= 64, i.e. l <= 24): the order is the
 *    lexicographic (rank_or_inf, min over the two frames of
 *    (Sigma << l | Gamma ^ toggle)), packed as
 *        value = rank_or_inf << 2l | min(Sigma << l | Gamma ^ t,
 *                                       Sigma_rc << l | Gamma_rc ^ t)
 *    with rank_or_inf = the table's infinity (0xFFFF) for an infinite center,
 *    so infinite windows are ordered too and the value is injective on
 *    canonical l-mers (the winning frame's (Sigma, Gamma) determine the DNA).
 *    This is bit for bit sshash's custom_scheme_hasher hash.
 * Rank-first comparison works in both modes: finite < infinite, and the rank
 * alone is rank_of_sigma() / rank_at() (one table load).
 *
 * Spaced seed / per-strand masking: the table lookup (rank, ybit) must see the
 * UNMASKED center (masking bits of y destroys the density the table was designed
 * for: kraken2's default 7-space mask lands in the center and costs 0.3 in
 * density factor at l = 31); the mask only coarsens the tie-break fields. The
 * *_masked entry points take the unmasked Sigma for the lookup and the masked
 * Sigma/Gamma planes of each strand for Sx, Sz and Gamma.
 *
 * Also here: a slow per-base reference (tests) and a small density helper.
 */

#ifndef MINIMIZER_ORDERS_GREEDYMINI_CENTER_H_
#define MINIMIZER_ORDERS_GREEDYMINI_CENTER_H_

#include "greedymini_rank.h"

#include <cstdint>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace gm_center {

inline uint64_t reverse_bits64(uint64_t v) {
  v = __builtin_bswap64(v);
  v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
  v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
  v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
  return v;
}

// Reverse the low n bits (n in 1..64) of v.
inline uint64_t reverse_low_bits(uint64_t v, int n) {
  return reverse_bits64(v) >> (64 - n);
}

// Mirror (bit reversal) of a gm_k-bit key.
inline uint32_t mirror_key(uint32_t key, int gm_k) {
  return (uint32_t)reverse_low_bits(key, gm_k);
}

// Key bits [i, i + gm_k) of an l-bit Sigma plane counted from the MSB, i.e.
// the center y of x y z.
inline uint32_t center_key(uint64_t sigma, int i, int gm_k) {
  return (uint32_t)(sigma >> i) & (uint32_t)((1ULL << gm_k) - 1ULL);
}

// ---------------------------------------------------------------------------
// Sigma / Gamma plane extraction from a word of 2-bit bases (pext on BMI2, a
// SWAR fallback otherwise: gm_pext_even / gm_pext_odd of greedymini_rank.h).
// Bit p of a plane describes the base held in word bits 2p, 2p+1, so a word
// that keeps its first base in the highest bits yields planes with the first
// base in the most significant bit, the convention of Order. The two bases
// encodings in use differ only in where Sigma sits:
//   SigmaXor = false: Sigma = bit 0, Gamma = bit 1 (sshash: A=00 C=01 G=11
//                     T=10)
//   SigmaXor = true:  Sigma = bit 0 ^ bit 1, Gamma = bit 1 (kraken2 and kmc
//                     ACGT = 0..3: A,T -> 0 and C,G -> 1; same planes as the
//                     ATCG remap of gm_acgt_to_atcg, without forming it)
// `flip` complements Sigma (the other polarity). The rolling helpers keep the
// Sigma plane of a sliding window one base per step instead of one pext.
// ---------------------------------------------------------------------------
template <bool SigmaXor>
class TwoBitPlanes {
 public:
  TwoBitPlanes() = default;
  // n bases (1..32 for sigma/gamma/planes; sigma_wide takes up to 64).
  explicit TwoBitPlanes(int n, bool flip = false)
      : n_(n), flip_bit_(flip ? 1ULL : 0ULL) {
    const uint64_t evens = 0x5555555555555555ULL;
    even_ = (n >= 32) ? evens : (evens & ((1ULL << (2 * n)) - 1ULL));
    odd_ = even_ << 1;
    n_mask_ = (n >= 64) ? ~0ULL : ((1ULL << n) - 1ULL);
    flip_n_ = flip_bit_ ? n_mask_ : 0ULL;
  }

  int n() const { return n_; }
  uint64_t mask() const { return n_mask_; }

  // The n-bit planes of a word of n bases (word bits above 2n are ignored).
  uint64_t sigma(uint64_t w) const {
    return gm_pext_even(SigmaXor ? (w ^ (w >> 1)) : w, even_) ^ flip_n_;
  }
  uint64_t gamma(uint64_t w) const { return gm_pext_odd(w, odd_); }

  // Sigma of a k-mer held in two words (low / high 64 bits, 32 bases each;
  // bits above 2k zero), k <= 64.
  uint64_t sigma_wide(uint64_t lo, uint64_t hi, int k) const {
    const uint64_t evens = 0x5555555555555555ULL;
    const uint64_t kmask = (k >= 64) ? ~0ULL : ((1ULL << k) - 1ULL);
    const uint64_t a = SigmaXor ? (lo ^ (lo >> 1)) : lo;
    const uint64_t b = SigmaXor ? (hi ^ (hi >> 1)) : hi;
    const uint64_t s = (gm_pext_even(a, evens) | (gm_pext_even(b, evens) << 32)) & kmask;
    return flip_bit_ ? s ^ kmask : s;
  }

  // Sigma bit (0/1) of one base given by its 2-bit code.
  uint64_t base_sigma(uint64_t code) const {
    return (SigmaXor ? ((code ^ (code >> 1)) & 1ULL) : (code & 1ULL)) ^ flip_bit_;
  }

 private:
  int n_ = 0;
  uint64_t even_ = 0, odd_ = 0, n_mask_ = 0, flip_n_ = 0, flip_bit_ = 0;
};

// Slide an n-bit Sigma plane by one base. The new base's Sigma bit `b`
// enters at the highest position (the lowest leaves) when bit p = base p of
// a word whose newest base is the highest one (sshash), or at the lowest
// position (the oldest drops off the top, `mask` = 2^n - 1) when the newest
// base is the lowest one and the first base is the MSB (kraken2's l-mer).
inline uint64_t roll_in_high(uint64_t sigma, uint64_t b, int n) {
  return (sigma >> 1) | (b << (n - 1));
}
inline uint64_t roll_in_low(uint64_t sigma, uint64_t b, uint64_t mask) {
  return ((sigma << 1) | b) & mask;
}

// Valid minimizer lengths: l = gm_k + 2i, l <= 31 (so 2l bits fit a word).
inline void check_length(int l, int gm_k) {
  if (gm_k > l)
    throw std::invalid_argument("GreedyMini gm_k exceeds minimizer length l");
  if (((l - gm_k) & 1) != 0)
    throw std::invalid_argument(
        "GreedyMini center extension needs l - gm_k even (l = " +
        std::to_string(l) + ", gm_k = " + std::to_string(gm_k) +
        "); valid l: gm_k, gm_k+2, ..., 31");
}

// Checks that `table` may be extended to length l: l - gm_k
// even, mirror symmetry rank[key] == rank[mirror(key)], every finite rank
// shared by one key and its mirror only (skipped when !require_rank_pairs: the
// injective packing does not need it). Returns the largest finite rank
// (0 if none). Throws std::invalid_argument with a descriptive message.
inline uint32_t validate_table(const GreedyMiniTable &table, int l,
                               bool require_rank_pairs = true) {
  const int gm_k = table.gm_k();
  check_length(l, gm_k);
  const size_t n = table.size();
  const uint16_t *r = table.data();
  // Keys per finite rank (saturating at 3), for the rank-pair check: with
  // the table mirror symmetric, every rank is held by exactly a key and its
  // mirror iff its count is 1 for a palindromic key and 2 otherwise. Only
  // allocated when checked (one byte per key, so no large temporary).
  std::vector<uint8_t> keys_of_rank;
  if (require_rank_pairs)
    keys_of_rank.assign(GreedyMiniTable::kInfRank, 0);
  uint32_t max_rank = 0;
  for (size_t key = 0; key < n; key++) {
    const size_t mir = mirror_key((uint32_t)key, gm_k);
    if (r[key] != r[mir])
      throw std::invalid_argument(
          "GreedyMini center extension needs a mirror-symmetric (canonical) "
          "table: rank[" + std::to_string(key) + "] != rank[" +
          std::to_string(mir) + "]");
    if (r[key] != GreedyMiniTable::kInfRank && require_rank_pairs &&
        keys_of_rank[r[key]] < 3)
      keys_of_rank[r[key]]++;
  }
  if (require_rank_pairs)
    for (size_t key = 0; key < n; key++) {
      if (r[key] == GreedyMiniTable::kInfRank)
        continue;
      const size_t mir = mirror_key((uint32_t)key, gm_k);
      if (keys_of_rank[r[key]] != (key == mir ? 1 : 2))
        throw std::invalid_argument(
            "GreedyMini center extension needs every rank shared by a key and "
            "its mirror only (rank " + std::to_string(r[key]) + ")");
    }
  for (size_t key = 0; key < n; key++)
    if (r[key] != GreedyMiniTable::kInfRank && r[key] > max_rank)
      max_rank = r[key];
  return max_rank;
}

class Order {
 public:
  Order() = default;

  // `toggle` is XORed into the Gamma plane (low l bits are used). `injective`
  // selects the total-order packing (see the header comment; needs
  // 16 + 2l <= 64). Throws like validate_table(), and if the value does not
  // fit 64 bits.
  Order(std::shared_ptr<const GreedyMiniTable> table, int l, uint64_t toggle,
        bool injective = false)
      : table_(std::move(table)), l_(l), injective_(injective) {
    const uint32_t max_rank = validate_table(*table_, l, !injective);
    gm_k_ = table_->gm_k();
    i_ = (l_ - gm_k_) / 2;
    int rank_bits = 1;
    while ((max_rank >> rank_bits) != 0)
      rank_bits++;
    rank_bits_ = rank_bits;
    sh_z_ = l_;
    sh_y_ = l_ + i_;
    sh_x_ = sh_y_ + 1;
    sh_rank_ = sh_x_ + i_;
    if (injective_) {
      sh_rank_ = 2 * l_;
      rank_bits_ = 16;
    }
    if (sh_rank_ + rank_bits_ > 64)
      throw std::invalid_argument(
          std::string("GreedyMini center extension: l = ") +
          std::to_string(l_) + " too large to pack the " +
          (injective_ ? "injective" : "compressed") + " order");
    z_mask_ = (1ULL << i_) - 1ULL;
    l_mask_ = (l_ >= 64) ? ~0ULL : ((1ULL << l_) - 1ULL);
    tb_ = toggle & l_mask_;
    ranks_ = table_->data();
    ybit_.assign((table_->size() + 63) / 64, 0);
    if (!injective_)
      for (size_t key = 0; key < table_->size(); key++)
        if (key > mirror_key((uint32_t)key, gm_k_))
          ybit_[key >> 6] |= 1ULL << (key & 63);
  }

  int l() const { return l_; }
  int gm_k() const { return gm_k_; }
  int center_offset() const { return i_; }  // i = (l - gm_k) / 2
  bool injective() const { return injective_; }
  // Bits of a value: compressed rank_bits + 2i + 1 + l, injective 16 + 2l.
  int value_bits() const { return sh_rank_ + rank_bits_; }
  // Compressed mode: a window whose center is not in the UHS has the value
  // UINT64_MAX. Injective mode: its rank field is kInfRank and the low bits
  // are a real tie-break.
  enum : uint16_t { kInfRank = GreedyMiniTable::kInfRank };

  // Rank (kInfRank for a center outside the UHS) of the l-bit Sigma plane of
  // an l-mer, or of the window at `shift` of a longer Sigma plane: the l-mer
  // whose own plane is (plane >> shift) & (2^l - 1) has its center key at
  // bits [shift + i, shift + i + gm_k) of the plane. One table load; the
  // rank-first comparison of two windows equals the comparison of their
  // values' top fields in both modes.
  uint16_t rank_of_sigma(uint64_t sigma) const {
    return ranks_[center_key(sigma, i_, gm_k_)];
  }
  uint16_t rank_at(uint64_t plane, int shift) const {
    return ranks_[(uint32_t)(plane >> (shift + i_)) &
                  (uint32_t)((1ULL << gm_k_) - 1ULL)];
  }
  uint16_t rank_of_key(uint32_t key) const { return ranks_[key]; }

  // Packed value of one strand (UINT64_MAX when the center is not in the UHS).
  uint64_t value_of_strand(uint64_t sigma, uint64_t gamma) const {
    const uint32_t key = center_key(sigma, i_, gm_k_);
    return pack(sigma, gamma, ranks_[key], key);
  }

  // Order value from the planes of one strand when the other strand is its
  // exact reverse complement and the table is read unmodified: one lookup
  // (the rc strand has the mirror center, hence the same rank), the rc planes
  // from bit reversals, the tie-break a branchless min. UINT64_MAX when the
  // center is not in the UHS.
  uint64_t value_of_planes(uint64_t sigma, uint64_t gamma) const {
    const uint32_t key = center_key(sigma, i_, gm_k_);
    const uint16_t r = ranks_[key];
    if (!injective_ && r == GreedyMiniTable::kInfRank)  // see value_of_both_strands
      return UINT64_MAX;
    const uint64_t sigma_rc = reverse_low_bits(sigma, l_);
    const uint64_t gamma_rc = reverse_low_bits(gamma, l_) ^ l_mask_;
    const uint64_t pa = pack(sigma, gamma, r, key);
    const uint64_t pb =
        pack(sigma_rc, gamma_rc, r, center_key(sigma_rc, i_, gm_k_));
    return pa < pb ? pa : pb;
  }

  // Order value from the planes of BOTH strands of an l-mer that the caller
  // already has (e.g. a rolling iterator): one lookup (the table is mirror
  // symmetric, so the rc strand has the same rank), no bit reversal. Same
  // result as value_of_planes(sigma, gamma) when (sigma_rc, gamma_rc) are the
  // exact reverse-complement planes.
  uint64_t value_of_both_strands(uint64_t sigma, uint64_t gamma,
                                 uint64_t sigma_rc, uint64_t gamma_rc) const {
    const uint32_t key = center_key(sigma, i_, gm_k_);
    const uint16_t r = ranks_[key];
    // injective_ first: constant per Order, so an injective scan does not
    // branch on whether the center is in the UHS (unpredictable per l-mer)
    if (!injective_ && r == GreedyMiniTable::kInfRank)
      return UINT64_MAX;
    const uint64_t pa = pack(sigma, gamma, r, key);
    const uint64_t pb = pack(sigma_rc, gamma_rc, r, mirror_key(key, gm_k_));
    return pa < pb ? pa : pb;
  }

  // Order value when the two strands' planes are given separately: the
  // smaller finite packed value, UINT64_MAX if both are infinite.
  uint64_t value_of_two_strands(uint64_t sigma_a, uint64_t gamma_a,
                                uint64_t sigma_b, uint64_t gamma_b) const {
    const uint64_t pa = value_of_strand(sigma_a, gamma_a);
    const uint64_t pb = value_of_strand(sigma_b, gamma_b);
    return pa < pb ? pa : pb;
  }

  // Masked variants (see the header comment): `sigma` is the unmasked Sigma
  // (rank and ybit come from its center), `sigma_m` / `gamma_m` the planes of
  // the masked strand word (Sx, Sz, Gamma come from them).
  uint64_t value_of_strand_masked(uint64_t sigma, uint64_t sigma_m,
                                  uint64_t gamma_m) const {
    const uint32_t key = center_key(sigma, i_, gm_k_);
    return pack(sigma_m, gamma_m, ranks_[key], key);
  }

  uint64_t value_of_two_strands_masked(uint64_t sigma_a, uint64_t sigma_ma,
                                       uint64_t gamma_ma, uint64_t sigma_b,
                                       uint64_t sigma_mb,
                                       uint64_t gamma_mb) const {
    const uint64_t pa = value_of_strand_masked(sigma_a, sigma_ma, gamma_ma);
    const uint64_t pb = value_of_strand_masked(sigma_b, sigma_mb, gamma_mb);
    return pa < pb ? pa : pb;
  }

 private:
  uint64_t pack(uint64_t sigma, uint64_t gamma, uint16_t rank,
                uint32_t key) const {
    if (injective_)
      return ((uint64_t)rank << sh_rank_) | (sigma << l_) | (gamma ^ tb_);
    if (rank == GreedyMiniTable::kInfRank)
      return UINT64_MAX;
    const uint64_t yb = (ybit_[key >> 6] >> (key & 63)) & 1ULL;
    return ((uint64_t)rank << sh_rank_) | ((sigma >> (l_ - i_)) << sh_x_) |
           (yb << sh_y_) | ((sigma & z_mask_) << sh_z_) | (gamma ^ tb_);
  }

  std::shared_ptr<const GreedyMiniTable> table_;
  const uint16_t *ranks_ = nullptr;
  int l_ = 0;
  bool injective_ = false;
  int gm_k_ = 0, i_ = 0;
  int sh_z_ = 0, sh_y_ = 0, sh_x_ = 0, sh_rank_ = 0, rank_bits_ = 0;
  uint64_t z_mask_ = 0, l_mask_ = 0, tb_ = 0;
  std::vector<uint64_t> ybit_;
};

// ---------------------------------------------------------------------------
// Slow per-base specification (for tests): works on per-position plane bits,
// position 0 = first base, one int (0/1) per base, builds the packed value one
// bit at a time and offers the lexicographic tuple the packed value must order
// like. Independent of the bit tricks of Order.
// ---------------------------------------------------------------------------
struct Reference {
  const GreedyMiniTable &table;
  int l;
  uint64_t toggle;  // low l bits, bit (l-1-p) belongs to position p
  bool injective;  // the injective packing (see the header comment)

  Reference(const GreedyMiniTable &t, int l_, uint64_t toggle_,
            bool injective_ = false)
      : table(t), l(l_), toggle(toggle_), injective(injective_) {}

  struct Strand {
    int rank;                // -1 = infinity
    uint64_t packed;         // UINT64_MAX when infinite
    std::vector<int> tuple;  // rank, Sigma(xyz) bits, Gamma^toggle bits
  };

  // sigma: unmasked Sigma (rank / ybit); sigma_tie: Sigma used for Sx, Sz and
  // the tuple (default: sigma); gamma: Gamma (masked if a mask applies).
  Strand strand(const std::vector<int> &sigma, const std::vector<int> &gamma,
                const std::vector<int> *sigma_tie = nullptr) const {
    const std::vector<int> &st = sigma_tie ? *sigma_tie : sigma;
    const int gk = table.gm_k(), i = (l - gk) / 2;
    uint64_t key = 0, mirror = 0;
    for (int p = i; p < i + gk; p++)
      key = (key << 1) | (uint64_t)sigma[(size_t)p];
    for (int p = i + gk - 1; p >= i; p--)
      mirror = (mirror << 1) | (uint64_t)sigma[(size_t)p];
    Strand s;
    const uint16_t r = table.primary((uint32_t)key);
    const bool inf = r == GreedyMiniTable::kInfRank;
    s.rank = (inf && !injective) ? -1 : (int)r;
    if (s.rank < 0) {
      s.packed = UINT64_MAX;
      return s;
    }
    uint64_t v = r;
    if (injective) {
      for (int p = 0; p < l; p++)
        v = (v << 1) | (uint64_t)st[(size_t)p];
    } else {
      for (int p = 0; p < i; p++)
        v = (v << 1) | (uint64_t)st[(size_t)p];
      v = (v << 1) | (uint64_t)(key > mirror);
      for (int p = i + gk; p < l; p++)
        v = (v << 1) | (uint64_t)st[(size_t)p];
    }
    s.tuple.push_back(s.rank);
    for (int p = 0; p < l; p++)
      s.tuple.push_back(p >= i && p < i + gk ? sigma[(size_t)p]
                                             : st[(size_t)p]);
    for (int p = 0; p < l; p++) {
      const int g = gamma[(size_t)p] ^ (int)((toggle >> (l - 1 - p)) & 1ULL);
      v = (v << 1) | (uint64_t)g;
      s.tuple.push_back(g);
    }
    s.packed = v;
    return s;
  }

  // Order value of an l-mer from its two strands' planes. The optional
  // sigma_tie_* are the masked Sigma planes (see Strand()).
  uint64_t value(const std::vector<int> &sa, const std::vector<int> &ga,
                 const std::vector<int> &sb, const std::vector<int> &gb,
                 const std::vector<int> *sta = nullptr,
                 const std::vector<int> *stb = nullptr) const {
    const uint64_t a = strand(sa, ga, sta).packed;
    const uint64_t b = strand(sb, gb, stb).packed;
    return a < b ? a : b;
  }

  // The lexicographically smaller finite tuple (rank, Sigma, Gamma^toggle) of
  // the two strands; false when both are infinite.
  bool tuple(const std::vector<int> &sa, const std::vector<int> &ga,
             const std::vector<int> &sb, const std::vector<int> &gb,
             std::vector<int> *out, const std::vector<int> *sta = nullptr,
             const std::vector<int> *stb = nullptr) const {
    const Strand a = strand(sa, ga, sta), b = strand(sb, gb, stb);
    if (a.rank < 0 && b.rank < 0)
      return false;
    if (a.rank < 0)
      *out = b.tuple;
    else if (b.rank < 0)
      *out = a.tuple;
    else
      *out = a.tuple < b.tuple ? a.tuple : b.tuple;
    return true;
  }
};

// ---------------------------------------------------------------------------
// Density helper. ranks[j] = order value of the l-mer starting at position j
// (UINT64_MAX loses to everything finite); a window is w consecutive l-mers
// and its minimizer is the leftmost minimum. Returns
//   (selected positions / windows) * (w + 1)
// (2 for a random order, 2/(w+1) per position; the table's design value for
// its own window w when l == gm_k). 0 if there is no full window.
// ---------------------------------------------------------------------------
inline double density_factor(const std::vector<uint64_t> &ranks, int w) {
  if (w < 1 || ranks.size() < (size_t)w)
    return 0.0;
  const size_t windows = ranks.size() - (size_t)w + 1;
  std::deque<size_t> q;  // increasing ranks, leftmost minimum at the front
  uint64_t selected = 0;
  long long last = -1;
  for (size_t j = 0; j < ranks.size(); j++) {
    while (!q.empty() && ranks[q.back()] > ranks[j])
      q.pop_back();
    q.push_back(j);
    if (j + 1 < (size_t)w)
      continue;
    if (q.front() + (size_t)w <= j)
      q.pop_front();
    if ((long long)q.front() != last)
      selected++;
    last = (long long)q.front();
  }
  return (double)selected / (double)windows * (double)(w + 1);
}

}  // namespace gm_center

#endif  // MINIMIZER_ORDERS_GREEDYMINI_CENTER_H_
