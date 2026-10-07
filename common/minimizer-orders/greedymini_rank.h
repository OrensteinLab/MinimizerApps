/*
 * GreedyMini rank calculation (Golan et al. 2025), shared by kraken2 and kmc.
 *
 * Binary UHS on {0,1}^gm_k sigma-extended to DNA, for an l-mer with l == gm_k
 * (a longer l uses the centre extension of greedymini_centre.h: kraken2's
 * binary order, sshash, minimap2; kmc only uses a table at its own length):
 *   DNA packing is ACGT = 0,1,2,3 (A=00, C=01, G=10, T=11).
 *   1. Remap each base ACGT -> ATCG bitplanes (A=00, T=01, C=10, G=11)
 *      so MSB = AT vs CG (Sigma) and LSB = within-pair residual (Gamma).
 *   2. key        = Sigma of the gm_k bases -> UHS primary
 *   3. tiebreaker = Gamma of the l bases (Extend-sigma lex tau)
 *   4. Orient by the better UHS primary across DNA fwd/RC, then pack the
 *      Gamma tiebreaker in that frame ("canonical gm_k, then sigma-extend").
 *      Spaced-seed AND is applied to each strand's word independently
 *      (never RC of a masked l-mer).
 *   Double-stranded use orients by the key so reverse-symmetric tables stay
 *   valid. DNA hot path uses BMI2 pext when available. DNA only (no protein).
 *
 * Paper / .gm encoding:
 *   - UHS members have contiguous finite ranks 0 .. m-1
 *   - Non-UHS keys store sentinel = table size ("infinity")
 *
 * In-memory table stores dense uint16 ranks with kInfRank for non-UHS.
 *
 * Both apps use this canonically (double-stranded orientation via
 * dna_reverse_complement below); there is no non-canonical code path.
 * See orders/README.md.
 */

#ifndef MINIMIZER_ORDERS_GREEDYMINI_RANK_H_
#define MINIMIZER_ORDERS_GREEDYMINI_RANK_H_

#include <cctype>
#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__BMI2__)
#include <immintrin.h>
#endif

// Reverse complement: kraken2 already has its own copy of this in
// src/dna_kmer.h (kraken2::dna_reverse_complement) for its general
// (non-order-specific) canonical k-mer handling, verified bit-identical to
// the one in dna_revcomp.h -- consumers that already have kraken2's copy in
// scope should keep using it; dna_revcomp.h is for consumers (like kmc) that
// don't. Not included from this file to avoid a same-name collision when
// both end up visible together (e.g. via `using namespace kraken2`).

// ---------------------------------------------------------------------------
// Bit helpers (closed-form masks; sigma-extend path)
// ---------------------------------------------------------------------------

// Even bit positions 0,2,4,... = LSBs of DNA bases (ATCG residual after remap).
inline uint64_t gm_dna_tb_pext_mask(int l) {
  const uint64_t evens = 0x5555555555555555ULL;
  if (l >= 32)
    return evens;
  return evens & ((1ULL << (2 * l)) - 1ULL);
}

// Odd bits in the high gm_k-base window = MSBs of the first gm_k bases
// (AT/CG after ACGT->ATCG remap).
inline uint64_t gm_dna_key_pext_mask(int l, int gm_k) {
  const uint64_t odds = 0xAAAAAAAAAAAAAAAAULL;
  const int lo = 2 * (l - gm_k);
  if (gm_k >= 32)
    return (lo >= 64) ? 0ULL : (odds << lo);
  const uint64_t window = ((1ULL << (2 * gm_k)) - 1ULL) << lo;
  return odds & window;
}

// Pack bits in even positions into a contiguous low field (SWAR).
inline uint64_t compress_even_bits(uint64_t x) {
  x &= 0x5555555555555555ULL;
  x = (x | (x >> 1)) & 0x3333333333333333ULL;
  x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0fULL;
  x = (x | (x >> 4)) & 0x00ff00ff00ff00ffULL;
  x = (x | (x >> 8)) & 0x0000ffff0000ffffULL;
  x = (x | (x >> 16)) & 0x00000000ffffffffULL;
  return x;
}

inline uint64_t compress_odd_bits(uint64_t x) {
  return compress_even_bits((x >> 1) & 0x5555555555555555ULL);
}

// pext even-bit plane (or SWAR fallback).
inline uint64_t gm_pext_even(uint64_t x, uint64_t even_mask) {
#if defined(__BMI2__)
  return _pext_u64(x, even_mask);
#else
  return compress_even_bits(x & even_mask);
#endif
}

// pext odd-bit plane (or SWAR fallback).
inline uint64_t gm_pext_odd(uint64_t x, uint64_t odd_mask) {
#if defined(__BMI2__)
  return _pext_u64(x, odd_mask);
#else
  return compress_even_bits((x & odd_mask) >> 1);
#endif
}

// ACGT (A=00,C=01,G=10,T=11) -> ATCG (A=00,T=01,C=10,G=11).
// After remap: MSB = AT vs CG, LSB = within-pair residual for lex tau.
// Verified bit-identical to kmc's own CMmer::DnaBaseToBinaryBit-based key
// derivation (exhaustively for gm_k<=6, sampled for p=9).
inline uint64_t gm_acgt_to_atcg(uint64_t lmer) {
  const uint64_t evens = 0x5555555555555555ULL;
  const uint64_t msb = (lmer >> 1) & evens;
  const uint64_t lsb = lmer & evens;
  return ((msb ^ lsb) << 1) | msb;
}

inline uint32_t gm_extract_dna_key_atcg(uint64_t atcg, int l, int gm_k) {
#if defined(__BMI2__)
  return (uint32_t)_pext_u64(atcg, gm_dna_key_pext_mask(l, gm_k));
#else
  uint64_t top = atcg >> (2 * (l - gm_k));
  if (gm_k < 32)
    top &= (1ULL << (2 * gm_k)) - 1ULL;
  return (uint32_t)compress_odd_bits(top);
#endif
}

inline uint64_t gm_extract_dna_tb_atcg(uint64_t atcg, int l) {
#if defined(__BMI2__)
  return _pext_u64(atcg, gm_dna_tb_pext_mask(l));
#else
  const uint64_t lmer_mask = (l >= 32) ? ~0ULL : ((1ULL << (2 * l)) - 1ULL);
  return compress_even_bits(atcg & lmer_mask);
#endif
}

inline uint32_t gm_extract_dna_key(uint64_t lmer, int l, int gm_k) {
  return gm_extract_dna_key_atcg(gm_acgt_to_atcg(lmer), l, gm_k);
}

inline uint64_t gm_extract_dna_tb(uint64_t lmer, int l) {
  return gm_extract_dna_tb_atcg(gm_acgt_to_atcg(lmer), l);
}

inline void gm_extract_dna(uint64_t lmer, int l, int gm_k,
                           uint32_t &key, uint64_t &tiebreaker) {
  const uint64_t atcg = gm_acgt_to_atcg(lmer);
  key = gm_extract_dna_key_atcg(atcg, l, gm_k);
  tiebreaker = gm_extract_dna_tb_atcg(atcg, l);
}

inline int gm_k_from_size(size_t table_size) {
  int gm_k = 0;
  size_t sz = table_size;
  while (sz > 1) {
    sz >>= 1;
    gm_k++;
  }
  return gm_k;
}

// ---------------------------------------------------------------------------
// Compact in-memory UHS rank table
// ---------------------------------------------------------------------------

class GreedyMiniTable {
 public:
  // Sentinel for non-UHS keys (paper "infinity"). Fits in uint16.
  // enum avoids C++11 ODR issues with static constexpr used by reference.
  enum : uint16_t { kInfRank = 0xFFFFu };

  GreedyMiniTable() : gm_k_(0), uhs_size_(0) {}

  explicit GreedyMiniTable(const std::vector<uint64_t> &raw_gm) {
    load_from_raw(raw_gm);
  }

  void load_from_raw(const std::vector<uint64_t> &raw_gm) {
    if (raw_gm.empty())
      throw std::invalid_argument("empty GreedyMini table");
    const size_t n = raw_gm.size();
    const int log2n = gm_k_from_size(n);
    if (((size_t)1 << log2n) != n)
      throw std::invalid_argument("GreedyMini table size is not a power of 2");

    gm_k_ = log2n;  // binary bits

    // Sentinel = table length (2^gm_k).
    const uint64_t inf_sentinel = (uint64_t)n;
    ranks_.assign(n, kInfRank);
    uhs_size_ = 0;
    for (size_t i = 0; i < n; i++) {
      const uint64_t r = raw_gm[i];
      if (r < inf_sentinel) {
        if (r >= kInfRank)
          throw std::invalid_argument("UHS rank does not fit in uint16");
        ranks_[i] = (uint16_t)r;
        uhs_size_++;
      } else if (r != inf_sentinel) {
        throw std::invalid_argument(
            "GreedyMini entry is neither UHS rank nor infinity sentinel");
      }
    }
  }

  bool empty() const { return ranks_.empty(); }
  int gm_k() const { return gm_k_; }
  // Bits used to index the table.
  int key_bits() const { return gm_k_; }
  uint32_t uhs_size() const { return uhs_size_; }
  size_t size() const { return ranks_.size(); }
  const uint16_t *data() const { return ranks_.data(); }

  bool in_uhs(uint32_t key) const { return ranks_[key] != kInfRank; }
  uint16_t primary(uint32_t key) const { return ranks_[key]; }

 private:
  int gm_k_;
  uint32_t uhs_size_;
  std::vector<uint16_t> ranks_;
};

// Classic .gm / gm.k2d I/O: [size_t len][len x uint64_t]
inline std::vector<uint64_t> load_gm_file_raw(const std::string &filename) {
  std::vector<uint64_t> vec;
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open())
    return vec;
  size_t size = 0;
  file.read(reinterpret_cast<char *>(&size), sizeof(size));
  if (!file)
    return {};
  vec.resize(size);
  file.read(reinterpret_cast<char *>(vec.data()),
            static_cast<std::streamsize>(size * sizeof(uint64_t)));
  if (!file)
    return {};
  return vec;
}

inline void write_gm_file_raw(const std::string &filename,
                              const std::vector<uint64_t> &ranks) {
  std::ofstream file(filename, std::ios::binary);
  if (!file)
    throw std::runtime_error("Unable to write GreedyMini table: " + filename);
  const size_t gm_len = ranks.size();
  file.write(reinterpret_cast<const char *>(&gm_len), sizeof(gm_len));
  file.write(reinterpret_cast<const char *>(ranks.data()),
             static_cast<std::streamsize>(gm_len * sizeof(uint64_t)));
  if (!file)
    throw std::runtime_error("Error writing GreedyMini table: " + filename);
}

inline std::shared_ptr<const GreedyMiniTable>
load_gm_table(const std::string &filename) {
  auto raw = load_gm_file_raw(filename);
  if (raw.empty())
    return nullptr;
  return std::make_shared<const GreedyMiniTable>(raw);
}

// Export compact table back to classic .gm encoding (UHS ranks / n sentinel).
inline std::vector<uint64_t>
gm_table_to_raw(const GreedyMiniTable &table) {
  const uint64_t inf = (uint64_t)table.size();
  std::vector<uint64_t> raw(table.size(), inf);
  const uint16_t *p = table.data();
  for (size_t i = 0; i < table.size(); i++) {
    if (p[i] != GreedyMiniTable::kInfRank)
      raw[i] = p[i];
  }
  return raw;
}

#endif  // MINIMIZER_ORDERS_GREEDYMINI_RANK_H_
