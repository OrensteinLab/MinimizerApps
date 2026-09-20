/*
 * GreedyMini rank calculation (Golan et al. 2025), shared by kraken2 and kmc.
 *
 * Two alphabet modes (GmAlphabet):
 *
 * BINARY (default): binary UHS on {0,1}^gm_k sigma-extended to DNA, with
 * Extend-k when l > gm_k:
 *   DNA packing is ACGT = 0,1,2,3 (A=00, C=01, G=10, T=11).
 *   1. Remap each base ACGT -> ATCG bitplanes (A=00, T=01, C=10, G=11)
 *      so MSB = AT vs CG (Sigma) and LSB = within-pair residual (Gamma).
 *   2. key         = Sigma of the first gm_k bases -> UHS primary
 *   3. suffix_sigma = Sigma of the remaining (l - gm_k) bases (Extend-k lex)
 *   4. tiebreaker  = Gamma of all l bases (Extend-sigma lex tau)
 *   5. Orient by better gm_k UHS primary across DNA fwd/RC, then pack
 *      Extend-k / Gamma in that frame ("canonical gm_k, then sigma-extend").
 *      Spaced-seed AND is applied to each strand's word independently
 *      (never RC of a masked l-mer).
 *   Extend-k guarantees density does not increase when l > gm_k (SS);
 *   double-stranded use orients by gm_k so reverse-symmetric tables stay valid.
 *   DNA hot path uses BMI2 pext when available.
 *
 * DNA: native DNA UHS on {A,C,G,T}^gm_k (ACGT=0123), table size 4^gm_k:
 *   key        = top 2*gm_k bits (first gm_k bases)
 *   tiebreaker = low 2*(l - gm_k) bits (remaining bases)
 *   Orient by better gm_k UHS primary across DNA fwd/RC, then pack
 *   (primary << tb_shift) | (tiebreaker ^ tb_mask) in that frame.
 *
 * Paper / .gm encoding:
 *   - UHS members have contiguous finite ranks 0 .. m-1
 *   - Non-UHS keys store sentinel = table size ("infinity")
 *
 * In-memory table stores dense uint16 ranks with kInfRank for non-UHS.
 *
 * Both apps use this canonically (double-stranded orientation via
 * dna_reverse_complement below); there is no non-canonical code path.
 * kraken2 only ever instantiates the BINARY alphabet in practice (its gm_k
 * values are too large for a 4^gm_k DNA table); kmc uses both, at its own
 * much smaller p. See orders/README.md.
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
// Alphabet mode (persisted in IndexOptions.gm_alphabet; 0 = legacy binary)
// ---------------------------------------------------------------------------

enum class GmAlphabet : uint8_t {
  BINARY = 0,  // sigma-extend + tiebreaker
  DNA    = 1   // contiguous DNA key + suffix tiebreaker
};

inline const char *gm_alphabet_name(GmAlphabet a) {
  return a == GmAlphabet::DNA ? "dna" : "binary";
}

inline GmAlphabet parse_gm_alphabet(const std::string &s) {
  std::string lower = s;
  for (char &c : lower)
    c = (char)std::tolower((unsigned char)c);
  if (lower == "binary")
    return GmAlphabet::BINARY;
  if (lower == "dna")
    return GmAlphabet::DNA;
  throw std::invalid_argument(
      "unknown gm-alphabet '" + s + "' (choices: binary, dna)");
}

// ---------------------------------------------------------------------------
// Bit helpers (closed-form masks; binary / sigma-extend path)
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

// Odd bits in the low (l - gm_k)-base window = Extend-k suffix Sigma.
inline uint64_t gm_dna_suffix_sigma_pext_mask(int l, int gm_k) {
  const int suf = l - gm_k;
  if (suf <= 0)
    return 0;
  const uint64_t odds = 0xAAAAAAAAAAAAAAAAULL;
  if (suf >= 32)
    return odds;
  return odds & ((1ULL << (2 * suf)) - 1ULL);
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

inline uint64_t gm_extract_suffix_sigma_atcg(uint64_t atcg, int l, int gm_k) {
  const int suf = l - gm_k;
  if (suf <= 0)
    return 0;
#if defined(__BMI2__)
  return _pext_u64(atcg, gm_dna_suffix_sigma_pext_mask(l, gm_k));
#else
  const uint64_t lo =
      (suf >= 32) ? atcg : (atcg & ((1ULL << (2 * suf)) - 1ULL));
  return compress_odd_bits(lo);
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

// Binary / protein: Sigma = MSB of each base, Gamma = remaining bits.
// key = Sigma[0..gm_k), suffix_sigma = Sigma[gm_k..l), tb = Gamma[0..l).
inline void gm_extract_generic(uint64_t lmer, int l, int gm_k, int bpc,
                               uint32_t &key, uint64_t &suffix_sigma,
                               uint64_t &tiebreaker) {
  const uint64_t lsb_mask = (1u << (bpc - 1)) - 1;
  uint64_t tmp = lmer << (64 - bpc * l);
  key = 0;
  suffix_sigma = 0;
  tiebreaker = 0;
  for (int i = 0; i < l; i++, tmp <<= bpc) {
    const uint64_t base = tmp >> (64 - bpc);
    const uint64_t msb = base >> (bpc - 1);
    if (i < gm_k)
      key = (key << 1) | (uint32_t)msb;
    else
      suffix_sigma = (suffix_sigma << 1) | msb;
    tiebreaker = (tiebreaker << (bpc - 1)) | (base & lsb_mask);
  }
}

// DNA-alphabet key: contiguous top 2*gm_k bits (first gm_k ACGT bases).
inline uint32_t gm_extract_dna_alphabet_key(uint64_t lmer, int l, int gm_k) {
  const int key_bits = 2 * gm_k;
  const uint64_t mask =
      (key_bits >= 64) ? ~0ULL : ((1ULL << key_bits) - 1ULL);
  return (uint32_t)((lmer >> (2 * (l - gm_k))) & mask);
}

// DNA-alphabet tiebreaker: contiguous low 2*(l - gm_k) bits (suffix bases).
inline uint64_t gm_extract_dna_alphabet_tb(uint64_t lmer, int l, int gm_k) {
  const int tb_bits = 2 * (l - gm_k);
  if (tb_bits <= 0)
    return 0;
  if (tb_bits >= 64)
    return lmer;
  return lmer & ((1ULL << tb_bits) - 1ULL);
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

  GreedyMiniTable()
      : alphabet_(GmAlphabet::BINARY), gm_k_(0), uhs_size_(0) {}

  explicit GreedyMiniTable(const std::vector<uint64_t> &raw_gm,
                           GmAlphabet alphabet = GmAlphabet::BINARY) {
    load_from_raw(raw_gm, alphabet);
  }

  void load_from_raw(const std::vector<uint64_t> &raw_gm,
                     GmAlphabet alphabet = GmAlphabet::BINARY) {
    if (raw_gm.empty())
      throw std::invalid_argument("empty GreedyMini table");
    alphabet_ = alphabet;
    const size_t n = raw_gm.size();
    const int log2n = gm_k_from_size(n);
    if (((size_t)1 << log2n) != n)
      throw std::invalid_argument("GreedyMini table size is not a power of 2");

    if (alphabet_ == GmAlphabet::DNA) {
      if ((log2n % 2) != 0)
        throw std::invalid_argument(
            "DNA GreedyMini table size must be 4^k (even log2)");
      gm_k_ = log2n / 2;  // bases
    } else {
      gm_k_ = log2n;  // binary bits
    }

    // Sentinel = table length (2^gm_k binary, or 4^gm_k DNA).
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
  GmAlphabet alphabet() const { return alphabet_; }
  int gm_k() const { return gm_k_; }
  // Bits used to index the table (binary: gm_k; DNA: 2*gm_k).
  int key_bits() const {
    return alphabet_ == GmAlphabet::DNA ? 2 * gm_k_ : gm_k_;
  }
  uint32_t uhs_size() const { return uhs_size_; }
  size_t size() const { return ranks_.size(); }
  const uint16_t *data() const { return ranks_.data(); }

  bool in_uhs(uint32_t key) const { return ranks_[key] != kInfRank; }
  uint16_t primary(uint32_t key) const { return ranks_[key]; }

 private:
  GmAlphabet alphabet_;
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
load_gm_table(const std::string &filename,
              GmAlphabet alphabet = GmAlphabet::BINARY) {
  auto raw = load_gm_file_raw(filename);
  if (raw.empty())
    return nullptr;
  return std::make_shared<const GreedyMiniTable>(raw, alphabet);
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
