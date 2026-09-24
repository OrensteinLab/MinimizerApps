/*
 * Streaming mod-sampling (Groot Koerkamp & Pibiri 2024).
 *
 * For Kraken (k,l) with w=k-l+1 and t = r+((l-r)%w), r=4:
 *   ranks successive canonical t-mers with Order (RandomXor, RandomPcg, or GM),
 *   maintains w residue deques, returns the ℓ-mer at position
 *   p = argmin_j front_rank(j) (ties → smallest j).
 *
 * Spaced-seed mask is applied to stored/returned canonical ℓ-mers only.
 * t-mer ranking is unmasked by design (selection order vs index key).
 * MOD+GM ranks raw t-mers with GMRanker (ACGT→ATCG σ-extend), same as
 * MINIMIZER+GM on ℓ-mers.
 */

#ifndef KRAKEN2_MOD_MINIMIZER_H_
#define KRAKEN2_MOD_MINIMIZER_H_

#include "minimizer_rankers.h"
#include "sketching_schemes.h"
#include "dna_kmer.h"
#include "../../common/minimizer-orders/mod_residue_deque.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#define MOD_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define MOD_ALWAYS_INLINE inline
#endif

namespace kraken2 {

class ModSampler {
 public:
  static constexpr int kDefaultR = 4;
  static constexpr int kMaxResidues = 16;

  ModSampler(int k, int l, int bpc, uint64_t toggle_mask, bool dna,
             Order order, int revcom_version = 1,
             uint64_t spaced_seed_mask = 0,
             std::shared_ptr<const GreedyMiniTable> gm_table = nullptr)
      : k_(k), l_(l), bpc_(bpc), dna_bpc2_(dna && bpc == 2), order_(order),
        revcom_version_(revcom_version),
        w_(k - l + 1), t_(0), tmer_span_(0),
        t_mask_(0), toggle_t_(0), lmer_mask_(0), spaced_seed_mask_(0),
        t_rc_shift_(0), l_rc_shift_(0),
        tmer_fwd_(0), tmer_rc_(0), loaded_t_(0),
        lmer_fwd_(0), lmer_rc_(0), loaded_l_(0),
        abs_base_pos_(0), abs_tmer_pos_(0), win_start_(-1), win_mod_(0),
        residue_a_(0), lmer_ring_ix_(0), lmers_stored_(0),
        ranked_tmer_log_(nullptr), ranked_lmer_log_(nullptr) {
    if (w_ < 2 || w_ > kMaxResidues)
      throw std::invalid_argument("ModSampler requires 2 <= w <= kMaxResidues");
    if (l < kDefaultR)
      throw std::invalid_argument("ModSampler requires l >= r=4");
    t_ = kDefaultR + ((l - kDefaultR) % w_);
    if (t_ < 1 || t_ > l)
      throw std::invalid_argument("ModSampler: invalid t");
    tmer_span_ = w_ + l - t_;
    const int t_bits = t_ * bpc_;
    t_mask_ = (t_bits >= 64) ? ~0ULL : ((1ULL << t_bits) - 1ULL);
    toggle_t_ = toggle_mask & t_mask_;
    const int l_bits = l_ * bpc_;
    lmer_mask_ = (l_bits >= 64) ? ~0ULL : ((1ULL << l_bits) - 1ULL);
    // Mask canonical ℓ-mers before store/return. t-mer ranks stay unmasked.
    spaced_seed_mask_ = spaced_seed_mask & lmer_mask_;
    if (dna_bpc2_) {
      t_rc_shift_ = 2 * (t_ - 1);
      l_rc_shift_ = 2 * (l_ - 1);
    }

    const int per_residue = tmer_span_ / w_ + 4;
    if (per_residue > ModResidueDeque::kMaxCap)
      throw std::invalid_argument("ModSampler residue deque capacity exceeded");
    for (int i = 0; i < w_; i++) {
      deques_[(size_t)i].cap = per_residue;
      front_rank_[(size_t)i] = UINT64_MAX;
    }

    for (int wm = 0; wm < w_; wm++) {
      for (int a = 0; a < w_; a++) {
        int p = a - wm;
        if (p < 0)
          p += w_;
        residue_pos_lut_[(size_t)wm][(size_t)a] = (uint8_t)p;
#if !defined(MOD_DISABLE_W5_FASTPATH)
        if (w_ == 5)
          ring_ix_lut_[(size_t)wm][(size_t)p] = (uint8_t)((wm + p) % 5);
#endif
      }
    }

    RankerConfig cfg;
    cfg.order = order;
    cfg.l = t_;
    cfg.bpc = bpc;
    cfg.toggle_mask = toggle_t_;
    cfg.lmer_mask = t_mask_;
    cfg.dna = dna;
    cfg.revcom_version = revcom_version_;
    // t-mer rankers do not AND the spaced mask (ℓ-mer store path does).
    cfg.spaced_seed_mask = 0;
    cfg.gm_table = gm_table;
    if (order == Order::GM) {
      if (!gm_table || gm_table->empty())
        throw std::invalid_argument("MOD+GM requires a GM table");
      if (gm_table->alphabet() != GmAlphabet::BINARY)
        throw std::invalid_argument(
            "MOD+GM requires BINARY alphabet (w5_k6_canonical)");
      if (gm_table->gm_k() != t_)
        throw std::invalid_argument("MOD+GM table gm_k must equal t");
      cfg.gm_alphabet = GmAlphabet::BINARY;
    }
    t_ranker_ = make_order_ranker(cfg);
  }

  void Reset() { ResetAt(0); }

  // Restart as if the non-ambiguous stream index of the next base is
  // `origin` (e.g. LoadSequence start). Seeds residue class and ℓ-mer ring
  // slot so chunked scans match a full-contig pass. After 'N', use Reset()/0.
  void ResetAt(ssize_t origin) {
    if (origin < 0)
      origin = 0;
    clear_deques();
    tmer_fwd_ = 0;
    tmer_rc_ = 0;
    loaded_t_ = 0;
    lmer_fwd_ = 0;
    lmer_rc_ = 0;
    loaded_l_ = 0;
    abs_base_pos_ = origin;
    abs_tmer_pos_ = origin;
    win_start_ = -1;
    win_mod_ = 0;
    residue_a_ = (int)(origin % w_);
    lmer_ring_ix_ = (int)(origin % w_);
    lmers_stored_ = 0;
    if (t_ranker_)
      t_ranker_->Reset();
  }

  void ClearRankStats() {
    if (t_ranker_)
      t_ranker_->ClearRankStats();
  }
  RankHitStats GetRankStats() const {
    return t_ranker_ ? t_ranker_->GetRankStats() : RankHitStats{};
  }
  bool has_uhs_stats() const {
    return order_ == Order::GM && t_ranker_ && t_ranker_->has_uhs_stats();
  }

  void SetRankedTmerLog(std::vector<uint64_t> *log) { ranked_tmer_log_ = log; }
  void SetRankedLmerLog(std::vector<uint64_t> *log) { ranked_lmer_log_ = log; }

  int t() const { return t_; }
  int w() const { return w_; }
  Order order() const { return order_; }

  MOD_ALWAYS_INLINE bool PushAndFinish(uint8_t base_code, ssize_t win,
                                       uint64_t *out_lmer) {
    if (dna_bpc2_)
      push_base_dna2(base_code);
    else
      push_base_generic(base_code);
    if (win < 0)
      return false;
    return finish_window(win, out_lmer);
  }

  bool PushBase(uint8_t base_code, uint64_t * /*unused*/) {
    if (dna_bpc2_)
      push_base_dna2(base_code);
    else
      push_base_generic(base_code);
    return false;
  }

  bool FinishWindow(ssize_t win, uint64_t *out_lmer) {
    return finish_window(win, out_lmer);
  }

  void OnAmbiguous() { Reset(); }

 private:
  inline uint64_t rank_tmer(uint64_t canonical_t, uint64_t raw_fwd_t) {
    if (order_ == Order::RANDOM_XOR)
      return canonical_t ^ toggle_t_;
    if (order_ == Order::GM)
      return t_ranker_->Rank(raw_fwd_t);
    return t_ranker_->Rank(canonical_t);
  }

  // Version-1 incremental RC matches dna_canonical_representation; older
  // revcom_version=0 uses the packed-word helper (pre-K2.0.8 encoding).
  MOD_ALWAYS_INLINE uint64_t canonical_dna_word(uint64_t fwd, uint64_t rc_inc,
                                               uint8_t n) const {
    if (revcom_version_ == DNA_REVCOM_VERSION)
      return fwd < rc_inc ? fwd : rc_inc;
    return dna_canonical_representation(fwd, n, revcom_version_);
  }

  MOD_ALWAYS_INLINE void push_base_dna2(uint8_t base_code) {
    abs_base_pos_++;
    const uint64_t base = (uint64_t)base_code;
    const uint64_t rev_base = base ^ 0x3ULL;

    tmer_fwd_ = ((tmer_fwd_ << 2) | base) & t_mask_;
    tmer_rc_ = (tmer_rc_ >> 2) | (rev_base << t_rc_shift_);
    tmer_rc_ &= t_mask_;
    if (loaded_t_ < t_)
      loaded_t_++;
    if (loaded_t_ >= t_) {
      const uint64_t canonical_t =
          canonical_dna_word(tmer_fwd_, tmer_rc_, (uint8_t)t_);
      if (ranked_tmer_log_)
        ranked_tmer_log_->push_back(order_ == Order::GM ? tmer_fwd_
                                                        : canonical_t);
      const uint64_t rank = rank_tmer(canonical_t, tmer_fwd_);
      const int a = residue_a_;
      ModResidueDeque &dq = deques_[(size_t)a];
      if (dq.empty()) {
        dq.push_monotone(ModTmerData{rank, abs_tmer_pos_});
        front_rank_[(size_t)a] = rank;
      } else {
        dq.push_monotone(ModTmerData{rank, abs_tmer_pos_});
      }
      residue_a_++;
      if (residue_a_ >= w_)
        residue_a_ = 0;
      abs_tmer_pos_++;
    }

    lmer_fwd_ = ((lmer_fwd_ << 2) | base) & lmer_mask_;
    lmer_rc_ = (lmer_rc_ >> 2) | (rev_base << l_rc_shift_);
    lmer_rc_ &= lmer_mask_;
    if (loaded_l_ < l_)
      loaded_l_++;
    if (loaded_l_ >= l_) {
      uint64_t canonical_l =
          canonical_dna_word(lmer_fwd_, lmer_rc_, (uint8_t)l_);
      if (spaced_seed_mask_)
        canonical_l &= spaced_seed_mask_;
      if (ranked_lmer_log_)
        ranked_lmer_log_->push_back(canonical_l);
      lmer_ring_[(size_t)lmer_ring_ix_] = canonical_l;
      lmer_ring_ix_++;
      if (lmer_ring_ix_ >= w_)
        lmer_ring_ix_ = 0;
      if (lmers_stored_ < w_)
        lmers_stored_++;
    }
  }

  void push_base_generic(uint8_t base_code) {
    abs_base_pos_++;
    const uint64_t base = (uint64_t)base_code;

    tmer_fwd_ = ((tmer_fwd_ << bpc_) | base) & t_mask_;
    if (loaded_t_ < t_)
      loaded_t_++;
    if (loaded_t_ >= t_) {
      if (ranked_tmer_log_)
        ranked_tmer_log_->push_back(tmer_fwd_);
      const uint64_t rank = rank_tmer(tmer_fwd_, tmer_fwd_);
      const int a = residue_a_;
      ModResidueDeque &dq = deques_[(size_t)a];
      if (dq.empty()) {
        dq.push_monotone(ModTmerData{rank, abs_tmer_pos_});
        front_rank_[(size_t)a] = rank;
      } else {
        dq.push_monotone(ModTmerData{rank, abs_tmer_pos_});
      }
      residue_a_++;
      if (residue_a_ >= w_)
        residue_a_ = 0;
      abs_tmer_pos_++;
    }

    lmer_fwd_ = ((lmer_fwd_ << bpc_) | base) & lmer_mask_;
    if (loaded_l_ < l_)
      loaded_l_++;
    if (loaded_l_ >= l_) {
      uint64_t lmer = lmer_fwd_;
      if (spaced_seed_mask_)
        lmer &= spaced_seed_mask_;
      if (ranked_lmer_log_)
        ranked_lmer_log_->push_back(lmer);
      lmer_ring_[(size_t)lmer_ring_ix_] = lmer;
      lmer_ring_ix_++;
      if (lmer_ring_ix_ >= w_)
        lmer_ring_ix_ = 0;
      if (lmers_stored_ < w_)
        lmers_stored_++;
    }
  }

  static inline void sync_front_rank(uint64_t *front_rank, int a,
                                     ModResidueDeque &dq) {
    front_rank[a] = dq.empty() ? UINT64_MAX : dq.front().rank;
  }

  // Update best (rank, p) with residue a; UINT64_MAX means empty residue.
  MOD_ALWAYS_INLINE static void consider_residue(uint64_t r, uint8_t p,
                                                 uint64_t *best_rank,
                                                 int *best_p, bool *any) {
    if (r == UINT64_MAX)
      return;
    if (!*any || r < *best_rank || (r == *best_rank && p < *best_p)) {
      *best_rank = r;
      *best_p = p;
      *any = true;
    }
  }

#if !defined(MOD_DISABLE_W5_FASTPATH)
  // Fully unrolled finish for default Kraken (k,l)=(35,31) → w=5.
  MOD_ALWAYS_INLINE bool finish_window_w5(ssize_t win, uint64_t *out_lmer) {
    if (win < 0 || lmers_stored_ < 5)
      return false;
    if (abs_tmer_pos_ < win + tmer_span_)
      return false;

    uint64_t best_rank = UINT64_MAX;
    int best_p = 0;
    bool any = false;
    const int wm = (int)(win % 5);
    const uint8_t *pos = residue_pos_lut_[(size_t)wm];

    if (win_start_ >= 0 && win == win_start_ + 1) {
      const int expire_a = (int)((win - 1) % 5);
      ModResidueDeque &dq = deques_[expire_a];
      while (!dq.empty() && dq.front().pos < win)
        dq.pop_front();
      sync_front_rank(front_rank_, expire_a, dq);
      consider_residue(front_rank_[0], pos[0], &best_rank, &best_p, &any);
      consider_residue(front_rank_[1], pos[1], &best_rank, &best_p, &any);
      consider_residue(front_rank_[2], pos[2], &best_rank, &best_p, &any);
      consider_residue(front_rank_[3], pos[3], &best_rank, &best_p, &any);
      consider_residue(front_rank_[4], pos[4], &best_rank, &best_p, &any);
    } else {
      for (int a = 0; a < 5; a++) {
        ModResidueDeque &dq = deques_[(size_t)a];
        while (!dq.empty() && dq.front().pos < win)
          dq.pop_front();
        front_rank_[(size_t)a] =
            dq.empty() ? UINT64_MAX : dq.front().rank;
      }
      consider_residue(front_rank_[0], pos[0], &best_rank, &best_p, &any);
      consider_residue(front_rank_[1], pos[1], &best_rank, &best_p, &any);
      consider_residue(front_rank_[2], pos[2], &best_rank, &best_p, &any);
      consider_residue(front_rank_[3], pos[3], &best_rank, &best_p, &any);
      consider_residue(front_rank_[4], pos[4], &best_rank, &best_p, &any);
    }

    if (!any)
      return false;

    *out_lmer = lmer_ring_[ring_ix_lut_[(size_t)wm][(size_t)best_p]];
    win_start_ = win;
    win_mod_ = wm;
    return true;
  }
#endif

  MOD_ALWAYS_INLINE bool finish_window(ssize_t win, uint64_t *out_lmer) {
#if !defined(MOD_DISABLE_W5_FASTPATH)
    if (w_ == 5)
      return finish_window_w5(win, out_lmer);
#endif
    if (win < 0 || lmers_stored_ < w_)
      return false;
    if (abs_tmer_pos_ < win + tmer_span_)
      return false;

    uint64_t best_rank = UINT64_MAX;
    int best_p = 0;
    const int wm = (int)(win % w_);
    const uint8_t *pos = residue_pos_lut_[(size_t)wm];
    bool any = false;

    if (win_start_ >= 0 && win == win_start_ + 1) {
      const int expire_a = (int)((win - 1) % w_);
      ModResidueDeque &dq = deques_[expire_a];
      while (!dq.empty() && dq.front().pos < win)
        dq.pop_front();
      sync_front_rank(front_rank_, expire_a, dq);
      for (int a = 0; a < w_; a++) {
        consider_residue(front_rank_[(size_t)a], pos[a], &best_rank, &best_p,
                         &any);
      }
    } else {
      for (int a = 0; a < w_; a++) {
        ModResidueDeque &dq = deques_[(size_t)a];
        while (!dq.empty() && dq.front().pos < win)
          dq.pop_front();
        const uint64_t r = dq.empty() ? UINT64_MAX : dq.front().rank;
        front_rank_[(size_t)a] = r;
        consider_residue(r, pos[a], &best_rank, &best_p, &any);
      }
    }

    if (!any)
      return false;

    *out_lmer = lmer_ring_[(size_t)((wm + best_p) % w_)];

    win_start_ = win;
    win_mod_ = wm;
    return true;
  }

  void clear_deques() {
    for (int i = 0; i < w_; i++) {
      deques_[(size_t)i].clear();
      front_rank_[(size_t)i] = UINT64_MAX;
    }
  }

  int k_, l_, bpc_;
  bool dna_bpc2_;
  Order order_;
  int revcom_version_;
  int w_, t_, tmer_span_;
  uint64_t t_mask_, toggle_t_, lmer_mask_, spaced_seed_mask_;
  int t_rc_shift_, l_rc_shift_;
  uint64_t tmer_fwd_, tmer_rc_;
  int loaded_t_;
  uint64_t lmer_fwd_, lmer_rc_;
  int loaded_l_;
  ssize_t abs_base_pos_;
  ssize_t abs_tmer_pos_;
  ssize_t win_start_;
  int win_mod_;
  int residue_a_;
  int lmer_ring_ix_;
  int lmers_stored_;
  uint8_t residue_pos_lut_[kMaxResidues][kMaxResidues];
  uint8_t ring_ix_lut_[kMaxResidues][kMaxResidues];
  uint64_t front_rank_[kMaxResidues];
  uint64_t lmer_ring_[kMaxResidues];
  ModResidueDeque deques_[kMaxResidues];
  std::unique_ptr<MinimizerRanker> t_ranker_;
  std::vector<uint64_t> *ranked_tmer_log_;
  std::vector<uint64_t> *ranked_lmer_log_;
};

}  // namespace kraken2

#endif
