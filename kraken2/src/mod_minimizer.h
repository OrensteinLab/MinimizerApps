/*
 * Streaming strand-symmetric mod-sampling (mod-sampling: Groot Koerkamp &
 * Pibiri 2024).
 *
 * For Kraken (k,l) with w=k-l+1 and t = r+((l-r)%w), r=4:
 *   ranks successive canonical t-mers with Order (RandomXor, RandomPcg, or GM)
 *   and maintains one monotone deque over all residue classes. Plain
 *   mod-sampling returns the ℓ-mer at offset p = x mod w over the window's
 *   minimal t-mers x (ties → smallest p).
 *
 * Spaced-seed mask is applied to stored/returned canonical ℓ-mers only.
 * t-mer ranking is unmasked by design (selection order vs index key).
 * MOD+GM ranks raw t-mers with GMRanker (ACGT→ATCG σ-extend), same as
 * MINIMIZER+GM on ℓ-mers.
 *
 * Plain mod-sampling is not strand-symmetric, so on DNA this sampler uses the
 * canonical variant of common/minimizer-orders/mod_minimizer_rule.h (the same
 * rule as the minimap2 and sshash forks' mod-minimizers): t == l (mod w), so
 * the offset p = x mod w of a minimal t-mer x and the mirror offset w-1-p in
 * the reverse complement correspond; a unique minimum samples its offset, and
 * a tie at the minimum takes the offset closest to the window centre
 * (w-1)/2, between two equally close ones the left if the window (the k bases
 * of the k-mer) is <= its reverse complement lexicographically, the right
 * otherwise (the centre rule of Cologni & Pibiri, as in SSHash). A read and
 * its reverse complement then give the same minimizers (the t-mer ranks are
 * strand-invariant for every order here), and the selection never moves
 * backwards. Protein (no strands) keeps the plain rule above.
 */

#ifndef KRAKEN2_MOD_MINIMIZER_H_
#define KRAKEN2_MOD_MINIMIZER_H_

#include "minimizer_rankers.h"
#include "sketching_schemes.h"
#include "dna_kmer.h"
#include "../../common/minimizer-orders/mod_minimizer_rule.h"

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
             std::shared_ptr<const GreedyMiniTable> gm_table = nullptr,
             bool gm_flip = false)
      : k_(k), l_(l), bpc_(bpc), dna_bpc2_(dna && bpc == 2), order_(order),
        revcom_version_(revcom_version),
        w_(k - l + 1), t_(0), tmer_span_(0),
        t_mask_(0), toggle_t_(0), lmer_mask_(0), spaced_seed_mask_(0),
        t_rc_shift_(0), l_rc_shift_(0),
        tmer_fwd_(0), tmer_rc_(0), loaded_t_(0),
        lmer_fwd_(0), lmer_rc_(0), loaded_l_(0),
        abs_tmer_pos_(0), win_start_(-1), win_mod_(0),
        residue_a_(0), lmer_ring_ix_(0), lmers_stored_(0),
        dq_head_(0), dq_tail_(0), n_ranked_(0), n_inf_(0),
        ranked_tmer_log_(nullptr), ranked_lmer_log_(nullptr) {
    if (w_ < 2 || w_ > kMaxResidues)
      throw std::invalid_argument("ModSampler requires 2 <= w <= kMaxResidues");
    // the old revcom_version 0 reverse complement is not exact, so its t-mer
    // ranks are not strand-invariant
    if (dna_bpc2_ && revcom_version != DNA_REVCOM_VERSION)
      throw std::invalid_argument("mod-sampling on DNA needs the current revcom version");
    if (dna_bpc2_ && k_ > 64)  // one bit per window base in a 64-bit word
      throw std::invalid_argument("mod-sampling on DNA needs k <= 64");
    if (l < kDefaultR)
      throw std::invalid_argument("ModSampler requires l >= r=4");
    t_ = kDefaultR + ((l - kDefaultR) % w_);
    if (t_ < 1 || t_ > l)
      throw std::invalid_argument("ModSampler: invalid t");
    if (dna_bpc2_ && !mod_rule::t_is_canonical(l_, w_, t_))
      throw std::invalid_argument("ModSampler: the centre rule needs t == l (mod w)");
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

    // One monotone deque over all t-mers of the window (<= tmer_span_ + 1
    // live entries); power-of-two ring, no division on the hot path.
    if ((uint32_t)tmer_span_ + 2 > kDqCap)
      throw std::invalid_argument("ModSampler window exceeds deque capacity");

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
    cfg.gm_flip = gm_flip;
    if (order == Order::GM) {
      if (!gm_table || gm_table->empty())
        throw std::invalid_argument("MOD+GM requires a GM table");
      if (gm_table->gm_k() != t_)
        throw std::invalid_argument("MOD+GM table gm_k must equal t");
    }
    t_ranker_ = make_order_ranker(cfg);

    // t-mer rank is a pure function of the (raw or canonical) t-mer word:
    // tabulate it for the table-driven orders (GM / Pcg) when 4^t is small.
    if (dna_bpc2_ && order_ != Order::RANDOM_XOR && t_ <= kMaxLutT) {
      const size_t n = (size_t)1 << (2 * t_);
      rank_lut_.resize(n);
      for (size_t x = 0; x < n; x++)
        rank_lut_[x] = t_ranker_->Rank((uint64_t)x);
      t_ranker_->ClearRankStats();
      use_lut_ = true;
    }
    select_kind();
  }

  void Reset() { ResetAt(0); }

  // Restart as if the non-ambiguous stream index of the next base is
  // `origin` (e.g. LoadSequence start). Seeds residue class and ℓ-mer ring
  // slot so chunked scans match a full-contig pass. After 'N', use Reset()/0.
  void ResetAt(ssize_t origin) {
    if (origin < 0)
      origin = 0;

    tmer_fwd_ = 0;
    tmer_rc_ = 0;
    loaded_t_ = 0;
    lmer_fwd_ = 0;
    lmer_rc_ = 0;
    loaded_l_ = 0;
    abs_tmer_pos_ = origin;
    win_start_ = -1;
    win_mod_ = 0;
    residue_a_ = (int)(origin % w_);
    lmer_ring_ix_ = (int)(origin % w_);
    lmers_stored_ = 0;
    dq_head_ = dq_tail_ = 0;
    whi_ = wlo_ = 0;
    if (t_ranker_)
      t_ranker_->Reset();
  }

  void ClearRankStats() {
    n_ranked_ = n_inf_ = 0;
    if (t_ranker_)
      t_ranker_->ClearRankStats();
  }
  RankHitStats GetRankStats() const {
    RankHitStats hs = t_ranker_ ? t_ranker_->GetRankStats() : RankHitStats{};
    if (use_lut_) {  // LUT bypasses the ranker's own counters
      hs.uhs += n_ranked_ - n_inf_;
      hs.non_uhs += n_inf_;
    }
    return hs;
  }
  bool has_uhs_stats() const {
    return order_ == Order::GM && t_ranker_ && t_ranker_->has_uhs_stats();
  }

  void SetRankedTmerLog(std::vector<uint64_t> *log) {
    ranked_tmer_log_ = log;
    select_kind();
  }
  void SetRankedLmerLog(std::vector<uint64_t> *log) {
    ranked_lmer_log_ = log;
    select_kind();
  }

  int t() const { return t_; }
  int w() const { return w_; }
  Order order() const { return order_; }

  MOD_ALWAYS_INLINE bool PushAndFinish(uint8_t base_code, ssize_t win,
                                       uint64_t *out_lmer) {
    push_dispatch(base_code);
    if (win < 0)
      return false;
    return finish_window(win, out_lmer);
  }

  bool PushBase(uint8_t base_code, uint64_t * /*unused*/) {
    push_dispatch(base_code);
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
    if (use_lut_) {
      const uint64_t r = rank_lut_[order_ == Order::GM ? raw_fwd_t : canonical_t];
      n_ranked_++;
      n_inf_ += (r == UINT64_MAX);
      return r;
    }
    if (order_ == Order::GM)
      return t_ranker_->Rank(raw_fwd_t);
    return t_ranker_->Rank(canonical_t);
  }

  // Monotone push (equal ranks are kept: they are ties between residues).
  MOD_ALWAYS_INLINE void push_tmer(uint64_t rank) {
    while (dq_tail_ != dq_head_ && dq_[(dq_tail_ - 1) & kDqMask].rank > rank)
      --dq_tail_;
    Ent &e = dq_[dq_tail_ & kDqMask];
    e.rank = rank;
    e.pos = abs_tmer_pos_;
    e.cls = residue_a_;
    ++dq_tail_;
    if (++residue_a_ >= w_)
      residue_a_ = 0;
    abs_tmer_pos_++;
  }

  // Version-1 incremental RC matches dna_canonical_representation; older
  // revcom_version=0 uses the packed-word helper (pre-K2.0.8 encoding).
  MOD_ALWAYS_INLINE uint64_t canonical_dna_word(uint64_t fwd, uint64_t rc_inc,
                                               uint8_t n) const {
    if (revcom_version_ == DNA_REVCOM_VERSION)
      return fwd < rc_inc ? fwd : rc_inc;
    return dna_canonical_representation(fwd, n, revcom_version_);
  }

  // Per-base push specialised at compile time on how the t-mer rank is taken
  // (the order and the LUT choice are fixed at construction). The debug-log
  // and generic-alphabet cases keep the runtime-branching code below.
  enum Kind {
    K_XOR = 0,       // rank = canonical t-mer ^ toggle (no table)
    K_LUT_RAW,       // LUT on the raw forward t-mer (GM)
    K_LUT_CAN,       // LUT on the canonical t-mer (Pcg)
    K_RANK_RAW,      // ranker on the raw forward t-mer (GM, 4^t too large)
    K_RANK_CAN,      // ranker on the canonical t-mer
    K_SLOW,          // logs enabled: generic runtime-branching path
    K_GENERIC        // non-2-bit alphabet
  };

  void select_kind() {
    if (!dna_bpc2_)
      kind_ = K_GENERIC;
    else if (ranked_tmer_log_ || ranked_lmer_log_)
      kind_ = K_SLOW;
    else if (order_ == Order::RANDOM_XOR)
      kind_ = K_XOR;
    else if (order_ == Order::GM)
      kind_ = use_lut_ ? K_LUT_RAW : K_RANK_RAW;
    else
      kind_ = use_lut_ ? K_LUT_CAN : K_RANK_CAN;
  }

  MOD_ALWAYS_INLINE void push_dispatch(uint8_t base_code) {
    if (dna_bpc2_) { // the window's k bases, for the centre rule's v <= RC(v)
      // bit planes of the recent bases, newest lowest
      whi_ = (whi_ << 1) | (uint64_t)(base_code >> 1);
      wlo_ = (wlo_ << 1) | (uint64_t)(base_code & 1);
    }
    switch (kind_) {
      case K_XOR:      push_base_k<K_XOR>(base_code); break;
      case K_LUT_RAW:  push_base_k<K_LUT_RAW>(base_code); break;
      case K_LUT_CAN:  push_base_k<K_LUT_CAN>(base_code); break;
      case K_RANK_RAW: push_base_k<K_RANK_RAW>(base_code); break;
      case K_RANK_CAN: push_base_k<K_RANK_CAN>(base_code); break;
      case K_SLOW:     push_base_dna2(base_code); break;
      default:         push_base_generic(base_code); break;
    }
  }

  template <int K>
  MOD_ALWAYS_INLINE void push_base_k(uint8_t base_code) {
    const bool kCanonical = (K == K_XOR || K == K_LUT_CAN || K == K_RANK_CAN);
    const uint64_t base = (uint64_t)base_code;
    const uint64_t rev_base = base ^ 0x3ULL;

    tmer_fwd_ = ((tmer_fwd_ << 2) | base) & t_mask_;
    if (kCanonical)
      tmer_rc_ = ((tmer_rc_ >> 2) | (rev_base << t_rc_shift_)) & t_mask_;
    if (loaded_t_ < t_)
      loaded_t_++;
    if (loaded_t_ >= t_) {
      uint64_t rank;
      if (K == K_XOR) {
        rank = canonical_dna_word(tmer_fwd_, tmer_rc_, (uint8_t)t_) ^ toggle_t_;
      } else if (K == K_LUT_RAW) {
        rank = rank_lut_[tmer_fwd_];
        n_ranked_++;
        n_inf_ += (rank == UINT64_MAX);
      } else if (K == K_LUT_CAN) {
        rank = rank_lut_[canonical_dna_word(tmer_fwd_, tmer_rc_, (uint8_t)t_)];
        n_ranked_++;
        n_inf_ += (rank == UINT64_MAX);
      } else if (K == K_RANK_RAW) {
        rank = t_ranker_->Rank(tmer_fwd_);
      } else {
        rank = t_ranker_->Rank(
            canonical_dna_word(tmer_fwd_, tmer_rc_, (uint8_t)t_));
      }
      push_tmer(rank);
    }
    store_lmer_dna2(base, rev_base);
  }

  MOD_ALWAYS_INLINE void store_lmer_dna2(uint64_t base, uint64_t rev_base) {
    lmer_fwd_ = ((lmer_fwd_ << 2) | base) & lmer_mask_;
    lmer_rc_ = ((lmer_rc_ >> 2) | (rev_base << l_rc_shift_)) & lmer_mask_;
    if (loaded_l_ < l_)
      loaded_l_++;
    if (loaded_l_ >= l_) {
      uint64_t canonical_l =
          canonical_dna_word(lmer_fwd_, lmer_rc_, (uint8_t)l_);
      if (spaced_seed_mask_)
        canonical_l &= spaced_seed_mask_;
      lmer_ring_[(size_t)lmer_ring_ix_] = canonical_l;
      lmer_ring_ix_++;
      if (lmer_ring_ix_ >= w_)
        lmer_ring_ix_ = 0;
      if (lmers_stored_ < w_)
        lmers_stored_++;
    }
  }

  void push_base_dna2(uint8_t base_code) {
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
      push_tmer(rank);
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
    const uint64_t base = (uint64_t)base_code;

    tmer_fwd_ = ((tmer_fwd_ << bpc_) | base) & t_mask_;
    if (loaded_t_ < t_)
      loaded_t_++;
    if (loaded_t_ >= t_) {
      if (ranked_tmer_log_)
        ranked_tmer_log_->push_back(tmer_fwd_);
      const uint64_t rank = rank_tmer(tmer_fwd_, tmer_fwd_);
      push_tmer(rank);
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

  // Window minimum over the single monotone deque. The entries equal to the
  // front rank are exactly the window t-mers of minimal rank. DNA picks among
  // them with centre_slot (the centre rule). Protein uses plain mod-sampling: the smallest
  // offset p = (pos - win) mod w (ties between residues -> smallest p),
  // identical to combining w residue deques. Non-finite ranks (UINT64_MAX,
  // non-UHS t-mers) are never candidates.
  MOD_ALWAYS_INLINE bool finish_window(ssize_t win, uint64_t *out_lmer) {
    if (win < 0 || lmers_stored_ < w_)
      return false;
    if (abs_tmer_pos_ < win + tmer_span_)
      return false;

    int wm;
    if (win_start_ >= 0 && win == win_start_ + 1) {
      wm = win_mod_ + 1;
      if (wm >= w_)
        wm = 0;
    } else {
      wm = (int)(win % w_);
    }
    win_start_ = win;
    win_mod_ = wm;

    while (dq_head_ != dq_tail_ && dq_[dq_head_ & kDqMask].pos < win)
      ++dq_head_;
    if (dq_head_ == dq_tail_)
      return false;
    const uint64_t best_rank = dq_[dq_head_ & kDqMask].rank;
    if (best_rank == UINT64_MAX)
      return false;
    if (dna_bpc2_) {
      *out_lmer = lmer_ring_[(size_t)centre_slot(best_rank, wm)];
      return true;
    }
    // protein: plain mod-sampling, smallest p over the minimal t-mers
    int best_a = dq_[dq_head_ & kDqMask].cls;
    int best_p = best_a - wm;
    if (best_p < 0)
      best_p += w_;
    for (uint32_t i = dq_head_ + 1; i != dq_tail_; ++i) {
      const Ent &e = dq_[i & kDqMask];
      if (e.rank != best_rank)
        break;
      int p = e.cls - wm;
      if (p < 0)
        p += w_;
      if (p < best_p) {
        best_p = p;
        best_a = e.cls;
      }
    }
    // The l-mer ring is indexed by residue class: slot (win + p) % w == cls.
    *out_lmer = lmer_ring_[(size_t)best_a];
    return true;
  }

  // DNA: l-mer ring slot of the window's selection. Equal minimal ranks
  // sit at the deque front in position order (a push pops only strictly
  // larger ranks). The l-mer ring is indexed by residue class, so a minimal
  // t-mer's l-mer is in the slot of its class; its offset in the window is
  // (cls - wm) mod w. A unique minimum (the common case) needs nothing else.
  MOD_ALWAYS_INLINE int centre_slot(uint64_t best_rank, int wm) const {
    const Ent &front = dq_[dq_head_ & kDqMask];
    uint32_t r = dq_head_ + 1;
    if (r == dq_tail_ || dq_[r & kDqMask].rank != best_rank)
      return front.cls;
    return tie_slot(best_rank, wm);
  }

  // A tie at the minimum: the centre rule over the tied entries' offsets.
  int tie_slot(uint64_t best_rank, int wm) const {
    mod_rule::CentrePick c;
    for (uint32_t r = dq_head_;
         r != dq_tail_ && dq_[r & kDqMask].rank == best_rank; ++r) {
      int p = dq_[r & kDqMask].cls - wm;
      if (p < 0)
        p += w_;
      c.add((uint64_t)p, (uint64_t)w_);
    }
    // the window is the last k pushed bases, newest in bit 0 of the planes
    const int p = (int)c.pick([this]() {
      return mod_rule::le_revcomp_planes(whi_, wlo_, k_);
    });
    const int slot = wm + p;
    return slot >= w_ ? slot - w_ : slot;
  }

  int k_, l_, bpc_;
  bool dna_bpc2_;
  Order order_;
  uint64_t whi_ = 0, wlo_ = 0;  // DNA: bits 1 / 0 of the last bases
  int revcom_version_;
  int w_, t_, tmer_span_;
  uint64_t t_mask_, toggle_t_, lmer_mask_, spaced_seed_mask_;
  int t_rc_shift_, l_rc_shift_;
  uint64_t tmer_fwd_, tmer_rc_;
  int loaded_t_;
  uint64_t lmer_fwd_, lmer_rc_;
  int loaded_l_;
  ssize_t abs_tmer_pos_;
  ssize_t win_start_;
  int win_mod_;
  int residue_a_;
  int lmer_ring_ix_;
  int lmers_stored_;
  static constexpr uint32_t kDqCap = 64;  // power of two
  static constexpr uint32_t kDqMask = kDqCap - 1;
  static constexpr int kMaxLutT = 8;      // 4^8 = 64K LUT entries (512 KB max)
  struct Ent {
    uint64_t rank;
    ssize_t pos;  // absolute t-mer start index since last reset
    int cls;      // pos mod w
  };
  Ent dq_[kDqCap];
  uint32_t dq_head_, dq_tail_;  // free-running; slot = index & kDqMask
  uint64_t lmer_ring_[kMaxResidues];
  std::vector<uint64_t> rank_lut_;
  bool use_lut_ = false;
  int kind_ = K_GENERIC;
  uint64_t n_ranked_, n_inf_;
  std::unique_ptr<MinimizerRanker> t_ranker_;
  std::vector<uint64_t> *ranked_tmer_log_;
  std::vector<uint64_t> *ranked_lmer_log_;
};

}  // namespace kraken2

#endif
