/*
 * Rankers used by MinimizerScanner: RandomXor, RandomPcg, GreedyMini.
 *
 * GreedyMini Rank() takes a raw l-mer (scanner DNA-lex is only the hash key).
 *
 * Binary DNA, l == gm_k (this is also what MOD+GM uses, with l = t = gm_k):
 *   1. Remap ACGT->ATCG (MSB = AT/CG Sigma, LSB = lex residual Gamma)
 *   2. primary <- table[Sigma of the l-mer] for each DNA strand (fwd and RC)
 *   3. Orient by the better primary (finite < infinity, then smaller rank).
 *      On primary tie, break with Gamma.
 *   4. Pack (primary << l) | (Gamma ^ toggle) from that strand.
 *
 * Binary DNA, l > gm_k: CENTRE extension (l - gm_k even, mirror-symmetric
 * table; valid l = gm_k, gm_k+2, ..., 31). The l-mer is x y z with
 * |x| = |z| = (l - gm_k)/2, |y| = gm_k; the order is table[Sigma(y)], then the
 * smaller of the two strand frames' Sigma(xyz), then the DNA, packed into one
 * uint64. Strand invariant, injective on canonical l-mers (no spaced seed),
 * and it keeps the table's designed density at every valid l. The algorithm
 * and its proofs/notes live in common/minimizer-schemes/greedymini_centre.h;
 * this file only extracts the Sigma/Gamma planes from kraken2's ACGT encoding
 * (TwoBitPlanes<true>). The scanner compares table ranks first (rank_first():
 * Sigma plane + one table load per l-mer) and asks for the full packed value
 * (value(), what Rank() returns) only when two finite ranks tie.
 *
 * Spaced seed (if any), centre path (l > gm_k): the rank lookup always sees
 * the UNMASKED centre (a mask that reaches into y would cost the table's
 * density; kraken2's default 7 spaces does at l = 31); each strand's word is
 * ANDed with the mask (never RC(x & mask)) only for the tie-break fields x, z,
 * Gamma, and the order is the smaller of the two strands' packed values, so
 * the mask never changes which table entry is looked up. With l == gm_k
 * (probe_strand) the table key is taken from the spaced, masked word, so
 * there the mask does change the entry looked up.
 *
 * With l == gm_k the table ranks the l-mer itself (Gamma breaks ties). The
 * GreedyMini order is DNA only: protein databases are refused.
 */

#ifndef KRAKEN2_MINIMIZER_RANKERS_H_
#define KRAKEN2_MINIMIZER_RANKERS_H_

#include "greedymini_rank.h"
#include "../../common/minimizer-schemes/greedymini_centre.h"
#include "dna_kmer.h"
#include "sketching_schemes.h"
#include "../../common/minimizer-schemes/pcg_rxs_m_xs.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

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
           int revcom_version = DNA_REVCOM_VERSION, bool flip = false)
      : l_(l), gm_k_(0), gamma_bits_(0),
        tb_bits_mask_(0), tb_mask_(0), spaced_mask_(0), key_pext_mask_(0),
        tb_pext_mask_(0),
        dna_(dna), flip_(flip), centre_(false), sym_(false), skip_tb_xor_(true),
        skip_space_(true),
        revcom_version_(revcom_version), table_(std::move(table)),
        ranks_ptr_(nullptr) {
    if (!table_ || table_->empty())
      throw std::invalid_argument("GreedyMini ranker requires a non-empty table");
    if (!dna_ || bpc != 2)
      throw std::invalid_argument(
          "GreedyMini is supported for DNA only (not for protein databases)");
    gm_k_ = table_->gm_k();
    if (gm_k_ > l_)
      throw std::invalid_argument("GreedyMini gm_k exceeds minimizer length l");
    ranks_ptr_ = table_->data();

    // l == gm_k: the table ranks the l-mer itself, Γ (l bits) breaks ties;
    // l > gm_k: the centre extension (greedymini_centre.h).
    gamma_bits_ = l_;
    centre_ = l_ > gm_k_;
    // uint16 primary << gamma must fit in uint64.
    if (!centre_ && gamma_bits_ > 48)
      throw std::invalid_argument(
          "GreedyMini tiebreaker width leaves no room to pack primary");

    tb_bits_mask_ =
        (gamma_bits_ >= 64) ? ~0ULL : ((1ULL << gamma_bits_) - 1ULL);
    tb_mask_ = toggle_mask & tb_bits_mask_;
    skip_tb_xor_ = (tb_mask_ == 0);

    spaced_mask_ = spaced_seed_mask & lmer_mask;
    skip_space_ = (spaced_mask_ == 0);

    key_pext_mask_ = gm_dna_key_pext_mask(l_, gm_k_);
    // --gm-projection at1-cg0 complements Sigma: the table key here, the
    // planes (key and Sigma tie-break) in the centre extension
    key_flip_ = flip_ ? (uint32_t)((1ULL << gm_k_) - 1ULL) : 0u;
    tb_pext_mask_ = gm_dna_tb_pext_mask(l_);
    odd_pext_mask_ = tb_pext_mask_ << 1;

    if (centre_)
      init_centre();
  }

  void Reset() override {}

  void ClearRankStats() override { stats_ = {}; }
  RankHitStats GetRankStats() const override { return stats_; }
  bool has_uhs_stats() const override { return true; }

  uint64_t Rank(uint64_t lmer) override { return rank_binary_dna(lmer); }

  // ---- rank-first interface (centre extension, exact reverse complement) ----
  // The order is (table rank, tie-break fields); the rank is one table load
  // from the Sigma plane of the UNMASKED l-mer (also with a spaced seed, and
  // strand invariant), so a scanner can compare ranks and ask for the full
  // value() only when two ranks are equal and finite. Infinite ranks all tie
  // (their value is UINT64_MAX): no value needed.
  bool rank_first() const { return centre_ && revcom_version_ == 1; }
  // Sigma plane of an ACGT l-mer word; the scanner may instead roll it
  // (sigma_bit / roll_in_low) one base per step.
  uint64_t sigma_of(uint64_t lmer) const { return planes_.sigma(lmer); }
  uint64_t sigma_bit(uint64_t base_code) const {
    return planes_.base_sigma(base_code);
  }
  uint64_t sigma_mask() const { return planes_.mask(); }
  enum : uint64_t { kInfRank = GreedyMiniTable::kInfRank };
  uint64_t rank_of_sigma(uint64_t sigma) {
    const uint64_t r = centre_order_.rank_of_sigma(sigma);
    if (r == kInfRank)
      stats_.non_uhs++;
    else
      stats_.uhs++;
    return r;
  }
  // The full packed value of an l-mer (what Rank() returns), no statistics.
  uint64_t value(uint64_t lmer) const { return rank_centre(lmer); }

  int gm_k() const { return gm_k_; }
  bool centre_extension() const { return centre_; }
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

  // ---- centre extension (binary DNA, l > gm_k): see greedymini_centre.h ----

  void init_centre() {
    // Validates l - gm_k even, mirror symmetry, rank sharing and packed width.
    centre_order_ = gm_centre::Order(table_, l_, tb_mask_);
    planes_ = gm_centre::TwoBitPlanes<true>(l_, flip_);
    sigma_flip_ = flip_ ? ((l_ >= 64) ? ~0ULL : ((1ULL << l_) - 1ULL)) : 0ULL;
    // Unmasked + exact reverse complement: the rc strand's planes come from
    // bit reversals of this strand's and share its rank.
    sym_ = skip_space_ && revcom_version_ == 1;
    // A mask that keeps or drops whole bases (every expanded kraken2 mask)
    // acts on the ACGT word directly; any other mask acts on the ATCG remap.
    const uint64_t evens = 0x5555555555555555ULL;
    mask_whole_bases_ = ((spaced_mask_ & evens) << 1) == (spaced_mask_ & ~evens);
  }

  // Planes of one strand word, l bits each, first base in the highest bit:
  // Sigma (AT/CG) of the unmasked word (rank lookup), and Sigma / Gamma of the
  // word with the spaced seed applied (tie-break fields).
  struct Planes {
    uint64_t S, Sm, Gm;
  };
  Planes planes(uint64_t acgt) const {
    Planes p;
    p.S = planes_.sigma(acgt);
    if (skip_space_) {
      p.Sm = p.S;
      p.Gm = planes_.gamma(acgt);
    } else if (mask_whole_bases_) {
      const uint64_t m = acgt & spaced_mask_;
      p.Sm = planes_.sigma(m);
      p.Gm = planes_.gamma(m);
    } else {
      // ATCG: Sigma in the odd bits, Gamma in the even bits
      const uint64_t atcg_m = gm_acgt_to_atcg(acgt) & spaced_mask_;
      p.Sm = gm_pext_odd(atcg_m, odd_pext_mask_) ^ sigma_flip_;
      p.Gm = gm_pext_even(atcg_m, tb_pext_mask_);
    }
    return p;
  }

  uint64_t rank_centre(uint64_t lmer) const {
    const Planes a = planes(lmer);
    if (sym_)
      return centre_order_.value_of_planes(a.S, a.Gm);
    const Planes b =
        planes(dna_reverse_complement(lmer, (uint8_t)l_, revcom_version_));
    return centre_order_.value_of_two_strands_masked(a.S, a.Sm, a.Gm, b.S,
                                                      b.Sm, b.Gm);
  }

  // One DNA orientation (l == gm_k): UHS key (+ ATCG word for the Γ pack).
  struct StrandView {
    uint16_t primary;
    uint64_t atcg;
  };

  uint64_t space(uint64_t word) const {
    return skip_space_ ? word : (word & spaced_mask_);
  }

  uint64_t rank_binary_dna(uint64_t lmer) {
    if (centre_)
      return record_oriented(rank_centre(lmer));
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
    const uint32_t key = (uint32_t)_pext_u64(atcg, key_pext_mask_) ^ key_flip_;
#else
    const uint32_t key = gm_extract_dna_key_atcg(atcg, l_, gm_k_) ^ key_flip_;
#endif
    return StrandView{ranks_ptr_[key], atcg};
  }

  uint64_t finish_strand(const StrandView &s) const {
    if (s.primary == GreedyMiniTable::kInfRank)
      return UINT64_MAX;
#if defined(__BMI2__)
    const uint64_t gamma = _pext_u64(s.atcg, tb_pext_mask_);
#else
    const uint64_t gamma = gm_extract_dna_tb_atcg(s.atcg, l_);
#endif
    return pack(s.primary, gamma);
  }

  // (primary << gamma_bits) | (Γ ^ toggle)
  uint64_t pack(uint16_t primary, uint64_t gamma) const {
    return ((uint64_t)primary << gamma_bits_) |
           (skip_tb_xor_ ? gamma : (gamma ^ tb_mask_));
  }

  int l_;
  int gm_k_;
  int gamma_bits_;  // Γ width
  uint64_t tb_bits_mask_;
  uint64_t tb_mask_;
  uint64_t spaced_mask_;
  uint64_t key_pext_mask_;
  uint64_t tb_pext_mask_;
  bool dna_;
  bool flip_;           // --gm-projection at1-cg0
  uint32_t key_flip_ = 0;
  uint64_t sigma_flip_ = 0;
  bool centre_;
  bool sym_;
  gm_centre::Order centre_order_;
  gm_centre::TwoBitPlanes<true> planes_;
  uint64_t odd_pext_mask_ = 0;
  bool mask_whole_bases_ = true;
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
  bool gm_flip = false;  // GreedyMini --gm-projection at1-cg0
  std::string greedy_mini_filename;
  std::shared_ptr<const GreedyMiniTable> gm_table;
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
                     cfg.gm_table, cfg.spaced_seed_mask, cfg.revcom_version,
                     cfg.gm_flip));
  }
  throw std::invalid_argument("unknown order");
}

}  // namespace kraken2

#endif
