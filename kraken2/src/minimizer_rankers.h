/*
 * Rankers used by MinimizerScanner: RandomXor, RandomPcg, GreedyMini.
 *
 * GreedyMini Rank() takes a raw ℓ-mer (scanner DNA-lex is only the hash key).
 *
 * Binary DNA:
 *   1. Remap ACGT→ATCG (MSB = AT/CG Σ, LSB = lex residual Γ)
 *   2. For each DNA strand (fwd and RC):
 *        primary ← table[Σ of first gm_k bases]
 *   3. Orient by the better gm_k UHS primary (finite ≺ ∞, then smaller rank).
 *      On primary tie, break with Extend-k suffix Σ then Γ.
 *   4. Pack ((primary << (ℓ−gm_k) | suffix_σ) << ℓ) | (Γ ^ toggle)
 *      from that strand. Spaced-seed AND is applied independently to each
 *      strand’s word (fwd and RC(raw)); never RC(x & mask).
 *
 *   This is “canonical gm_k, then Extend-k / Extend-σ” — not DNA-lex on the
 *   full ℓ-mer (which breaks reverse-symmetric tables when ℓ > gm_k).
 *
 *   DNA alphabet (4^gm_k table): same gm_k orientation; pack primary of the
 *   first gm_k bases with the suffix DNA bits in that frame.
 *   BINARY + protein: MSB Σ / residual Γ with same Extend-k pack (no RC).
 */

#ifndef KRAKEN2_MINIMIZER_RANKERS_H_
#define KRAKEN2_MINIMIZER_RANKERS_H_

#include "greedymini_rank.h"
#include "dna_kmer.h"
#include "sketching_schemes.h"
#include "../../common/minimizer-orders/pcg_rxs_m_xs.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

namespace kraken2 {

struct RankHitStats {
  uint64_t uhs = 0;      // Rank() returned a finite packed rank
  uint64_t non_uhs = 0;  // Rank() returned UINT64_MAX (infinity)
  uint64_t total() const { return uhs + non_uhs; }
};

class MinimizerRanker {
 public:
  virtual ~MinimizerRanker() = default;
  virtual void Reset() = 0;
  virtual uint64_t Rank(uint64_t lmer) = 0;
  virtual void ClearRankStats() {}
  virtual RankHitStats GetRankStats() const { return {}; }
  virtual bool has_uhs_stats() const { return false; }
};

class RandomRanker final : public MinimizerRanker {
 public:
  explicit RandomRanker(uint64_t toggle_mask) : toggle_(toggle_mask) {}
  void Reset() override {}
  uint64_t Rank(uint64_t canonical_lmer) override {
    return canonical_lmer ^ toggle_;
  }

 private:
  uint64_t toggle_;
};

class PcgRanker final : public MinimizerRanker {
 public:
  explicit PcgRanker(uint64_t toggle_mask) : toggle_(toggle_mask) {}
  void Reset() override {}
  uint64_t Rank(uint64_t canonical_lmer) override {
    return pcg_rxs_m_xs64(canonical_lmer ^ toggle_);
  }

 private:
  uint64_t toggle_;
};

class GMRanker final : public MinimizerRanker {
 public:
  GMRanker(int l, int bpc, uint64_t toggle_mask, uint64_t lmer_mask,
           bool dna, std::shared_ptr<const GreedyMiniTable> table,
           uint64_t spaced_seed_mask = 0,
           int revcom_version = DNA_REVCOM_VERSION)
      : l_(l), bpc_(bpc), gm_k_(0), ext_bits_(0), gamma_bits_(0),
        tb_bits_mask_(0), tb_mask_(0), spaced_mask_(0), key_pext_mask_(0),
        suffix_pext_mask_(0), tb_pext_mask_(0), alphabet_(GmAlphabet::BINARY),
        dna_(dna), skip_tb_xor_(true), skip_space_(true),
        revcom_version_(revcom_version), table_(std::move(table)),
        ranks_ptr_(nullptr) {
    if (!table_ || table_->empty())
      throw std::invalid_argument("GreedyMini ranker requires a non-empty table");
    alphabet_ = table_->alphabet();
    gm_k_ = table_->gm_k();
    if (gm_k_ > l_)
      throw std::invalid_argument("GreedyMini gm_k exceeds minimizer length l");
    if (alphabet_ == GmAlphabet::DNA && bpc_ != 2)
      throw std::invalid_argument(
          "DNA GreedyMini alphabet requires DNA (2-bit) sequence");
    ranks_ptr_ = table_->data();

    // Low-field width below the UHS primary:
    //   DNA alphabet: 2*(ℓ−gm_k) suffix DNA bits
    //   Binary:       Extend-k suffix Σ (ℓ−gm_k) + Γ (ℓ·(bpc−1))
    if (alphabet_ == GmAlphabet::DNA) {
      ext_bits_ = 0;
      gamma_bits_ = 2 * (l_ - gm_k_);
    } else {
      ext_bits_ = l_ - gm_k_;
      gamma_bits_ = l_ * (bpc_ - 1);
    }
    // uint16 primary << (ext+gamma) must fit in uint64.
    if (ext_bits_ + gamma_bits_ > 48)
      throw std::invalid_argument(
          "GreedyMini tiebreaker width leaves no room to pack primary");

    tb_bits_mask_ =
        (gamma_bits_ >= 64) ? ~0ULL : ((1ULL << gamma_bits_) - 1ULL);
    tb_mask_ = toggle_mask & tb_bits_mask_;
    skip_tb_xor_ = (tb_mask_ == 0);

    spaced_mask_ = spaced_seed_mask & lmer_mask;
    skip_space_ = (spaced_mask_ == 0);

    key_pext_mask_ = gm_dna_key_pext_mask(l_, gm_k_);
    suffix_pext_mask_ = gm_dna_suffix_sigma_pext_mask(l_, gm_k_);
    tb_pext_mask_ = gm_dna_tb_pext_mask(l_);
  }

  void Reset() override {}

  void ClearRankStats() override { stats_ = {}; }
  RankHitStats GetRankStats() const override { return stats_; }
  bool has_uhs_stats() const override { return true; }

  uint64_t Rank(uint64_t lmer) override {
    if (alphabet_ == GmAlphabet::DNA)
      return rank_dna_alphabet(lmer);
    if (!dna_)
      return rank_binary_protein(lmer);
    return rank_binary_dna(lmer);
  }

  int gm_k() const { return gm_k_; }
  GmAlphabet alphabet() const { return alphabet_; }
  uint32_t uhs_size() const { return table_->uhs_size(); }

 private:
  // Finite primary ≺ ∞, then smaller primary; equal primary → packed rest.
  template <typename FinishA, typename FinishB>
  uint64_t pick_better_primary(uint16_t pa, uint16_t pb, FinishA finish_a,
                               FinishB finish_b) const {
    const bool a_inf = (pa == GreedyMiniTable::kInfRank);
    const bool b_inf = (pb == GreedyMiniTable::kInfRank);
    if (a_inf && b_inf)
      return UINT64_MAX;
    if (a_inf)
      return finish_b();
    if (b_inf)
      return finish_a();
    if (pa < pb)
      return finish_a();
    if (pb < pa)
      return finish_b();
    const uint64_t ra = finish_a();
    const uint64_t rb = finish_b();
    return ra < rb ? ra : rb;
  }

  uint64_t record_oriented(uint64_t r) {
    if (r == UINT64_MAX)
      stats_.non_uhs++;
    else
      stats_.uhs++;
    return r;
  }

  // One DNA orientation: gm_k UHS key (+ ATCG word for later Extend-k pack).
  struct StrandView {
    uint16_t primary;
    uint64_t atcg;
  };

  uint64_t space(uint64_t word) const {
    return skip_space_ ? word : (word & spaced_mask_);
  }

  uint64_t pack_dna_alphabet(uint16_t primary, uint64_t word) const {
    if (primary == GreedyMiniTable::kInfRank)
      return UINT64_MAX;
#if defined(__BMI2__)
    const uint64_t tb = _bzhi_u64(word, (unsigned)gamma_bits_);
#else
    const uint64_t tb = word & tb_bits_mask_;
#endif
    return ((uint64_t)primary << gamma_bits_) |
           (skip_tb_xor_ ? tb : (tb ^ tb_mask_));
  }

  uint64_t rank_dna_alphabet(uint64_t lmer) {
    const uint64_t a = space(lmer);
    const uint16_t pa = ranks_ptr_[(uint32_t)(a >> gamma_bits_)];
    const uint64_t rc =
        dna_reverse_complement(lmer, (uint8_t)l_, revcom_version_);
    if (rc == lmer)
      return record_oriented(pack_dna_alphabet(pa, a));
    const uint64_t b = space(rc);
    const uint16_t pb = ranks_ptr_[(uint32_t)(b >> gamma_bits_)];
    return record_oriented(pick_better_primary(
        pa, pb, [&] { return pack_dna_alphabet(pa, a); },
        [&] { return pack_dna_alphabet(pb, b); }));
  }

  uint64_t rank_binary_dna(uint64_t lmer) {
    const StrandView a = probe_strand(lmer);
    const uint64_t rc =
        dna_reverse_complement(lmer, (uint8_t)l_, revcom_version_);
    if (rc == lmer)
      return record_oriented(finish_strand(a));
    const StrandView b = probe_strand(rc);
    return record_oriented(pick_better_primary(
        a.primary, b.primary, [&] { return finish_strand(a); },
        [&] { return finish_strand(b); }));
  }

  StrandView probe_strand(uint64_t lmer) const {
    const uint64_t atcg = space(gm_acgt_to_atcg(lmer));
#if defined(__BMI2__)
    const uint32_t key = (uint32_t)_pext_u64(atcg, key_pext_mask_);
#else
    const uint32_t key = gm_extract_dna_key_atcg(atcg, l_, gm_k_);
#endif
    return StrandView{ranks_ptr_[key], atcg};
  }

  uint64_t finish_strand(const StrandView &s) const {
    if (s.primary == GreedyMiniTable::kInfRank)
      return UINT64_MAX;
#if defined(__BMI2__)
    const uint64_t suffix = _pext_u64(s.atcg, suffix_pext_mask_);
    const uint64_t gamma = _pext_u64(s.atcg, tb_pext_mask_);
#else
    const uint64_t suffix = gm_extract_suffix_sigma_atcg(s.atcg, l_, gm_k_);
    const uint64_t gamma = gm_extract_dna_tb_atcg(s.atcg, l_);
#endif
    return pack_extend_k(s.primary, suffix, gamma);
  }

  uint64_t rank_binary_protein(uint64_t lmer) {
    uint32_t key = 0;
    uint64_t suffix = 0, gamma = 0;
    gm_extract_generic(space(lmer), l_, gm_k_, bpc_, key, suffix, gamma);
    const uint16_t primary = ranks_ptr_[key];
    if (primary == GreedyMiniTable::kInfRank) {
      stats_.non_uhs++;
      return UINT64_MAX;
    }
    stats_.uhs++;
    return pack_extend_k(primary, suffix, gamma);
  }

  // ((primary << ext_bits) | suffix_σ) << gamma_bits | (Γ ^ toggle)
  uint64_t pack_extend_k(uint16_t primary, uint64_t suffix,
                         uint64_t gamma) const {
    uint64_t mid = (uint64_t)primary << ext_bits_;
    if (ext_bits_ > 0)
      mid |= suffix;
    return (mid << gamma_bits_) |
           (skip_tb_xor_ ? gamma : (gamma ^ tb_mask_));
  }

  int l_;
  int bpc_;
  int gm_k_;
  int ext_bits_;    // Extend-k Σ suffix width (0 for DNA-alphabet mode)
  int gamma_bits_;  // Γ / DNA-alphabet suffix width
  uint64_t tb_bits_mask_;
  uint64_t tb_mask_;
  uint64_t spaced_mask_;
  uint64_t key_pext_mask_;
  uint64_t suffix_pext_mask_;
  uint64_t tb_pext_mask_;
  GmAlphabet alphabet_;
  bool dna_;
  bool skip_tb_xor_;
  bool skip_space_;
  int revcom_version_;
  std::shared_ptr<const GreedyMiniTable> table_;
  const uint16_t *ranks_ptr_;
  RankHitStats stats_;
};

struct RankerConfig {
  Order order = Order::RANDOM_XOR;
  int l = 0;
  int bpc = 2;
  uint64_t toggle_mask = 0;
  uint64_t lmer_mask = 0;
  bool dna = true;
  int revcom_version = DNA_REVCOM_VERSION;
  // MINIMIZER+GM: AND each strand independently inside GMRanker.
  // MOD t-mer rankers ignore this; ModSampler ANDs stored/returned ℓ-mers only.
  uint64_t spaced_seed_mask = 0;
  std::string greedy_mini_filename;
  std::shared_ptr<const GreedyMiniTable> gm_table;
  GmAlphabet gm_alphabet = GmAlphabet::BINARY;
};

inline std::unique_ptr<MinimizerRanker>
make_order_ranker(const RankerConfig &cfg) {
  if (cfg.order == Order::RANDOM_XOR)
    return std::unique_ptr<MinimizerRanker>(new RandomRanker(cfg.toggle_mask));
  if (cfg.order == Order::RANDOM_PCG)
    return std::unique_ptr<MinimizerRanker>(new PcgRanker(cfg.toggle_mask));
  if (cfg.order == Order::GM) {
    if (!cfg.gm_table || cfg.gm_table->empty())
      throw std::invalid_argument("GreedyMini order requires a loaded rank table");
    return std::unique_ptr<MinimizerRanker>(
        new GMRanker(cfg.l, cfg.bpc, cfg.toggle_mask, cfg.lmer_mask, cfg.dna,
                     cfg.gm_table, cfg.spaced_seed_mask, cfg.revcom_version));
  }
  throw std::invalid_argument("unknown order");
}

}  // namespace kraken2

#endif
