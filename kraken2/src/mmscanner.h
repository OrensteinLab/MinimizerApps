/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_MMSCANNER_H_
#define KRAKEN2_MMSCANNER_H_

#include "kraken2_headers.h"
#include "kv_store.h"
#include "minimizer_rankers.h"
#include "mod_minimizer.h"
#include "sketching_schemes.h"
#include "dna_kmer.h"

#include <memory>

namespace kraken2 {

const uint64_t DEFAULT_TOGGLE_MASK = 0xe37e28c4271b5a2dULL;
const uint64_t DEFAULT_SPACED_SEED_MASK = 0;
const int BITS_PER_CHAR_DNA = 2;
const int BITS_PER_CHAR_PRO = 4;
const int CURRENT_REVCOM_VERSION = 1;

struct MinimizerData {
  uint64_t rank;
  uint64_t candidate;
  ssize_t pos;
};

// Cycle accumulators for NextMinimizer phase split (bench builds only).
struct NextMinPhaseStats {
  uint64_t steps = 0;  // completed encode→canonical→rank→deque windows
  uint64_t encode_cycles = 0;
  uint64_t canonical_cycles = 0;
  uint64_t rank_cycles = 0;
  uint64_t deque_cycles = 0;

  uint64_t total_cycles() const {
    return encode_cycles + canonical_cycles + rank_cycles + deque_cycles;
  }
};

class MinimizerScanner {
  public:
  // Create scanner for seq over interval [start, finish)
  // Will report length l minimizers for all k-mers in the interval
  // Minimizers ordered by XORing candidates w/ toggle_mask
  //   toggle_mask == 0 implies lexicographical ordering
  MinimizerScanner(ssize_t k, ssize_t l,
                   uint64_t spaced_seed_mask = DEFAULT_SPACED_SEED_MASK,
                   bool dna_sequence = true,
                   uint64_t toggle_mask = DEFAULT_TOGGLE_MASK,
                   int revcom_version = CURRENT_REVCOM_VERSION);

  MinimizerScanner(ssize_t k, ssize_t l,
                   uint64_t spaced_seed_mask,
                   bool dna_sequence,
                   uint64_t toggle_mask,
                   int revcom_version,
                   SketchingScheme sketching,
                   Order order,
                   std::string greedy_mini_filename = "",
                   std::shared_ptr<const GreedyMiniTable> gm_table = nullptr);

  static std::vector<uint64_t> load_gm_file(const std::string &filename);
  static void write_gm_file(const std::string &filename,
                            const std::vector<uint64_t> &ranks);
  static std::string gm_k2d_path(const std::string &hashtable_filename);
  static std::shared_ptr<const GreedyMiniTable>
      shared_gm_table(const std::string &filename,
                      GmAlphabet alphabet = GmAlphabet::BINARY);
  static void register_shared_gm_table(
      const std::string &filename,
      std::shared_ptr<const GreedyMiniTable> table,
      GmAlphabet alphabet = GmAlphabet::BINARY);
  static void embed_gm_table(const std::string &source_path,
                             const std::string &hashtable_filename,
                             std::shared_ptr<const GreedyMiniTable> table);
  static int gm_k_from_size(size_t table_size);

  void LoadSequence(const std::string &seq, size_t start = 0,
      size_t finish = SIZE_MAX);
  // Chunked Mod builds: pass precomputed stream origin at `start` (non-ambig
  // count since last N / contig start).
  void LoadSequence(const std::string &seq, size_t start, size_t finish,
                    ssize_t mod_stream_origin);

  // origin_at[i] = Mod stream index when the next base to read is seq[i].
  void ComputeModStreamOrigins(const std::string &seq,
                               std::vector<uint32_t> *origin_at) const;

  uint64_t *NextMinimizer();
  // Return last minimizer, only valid if NextMinimizer last returned non-NULL
  uint64_t last_minimizer() const { return last_minimizer_; }
  ssize_t k() const { return k_; }
  ssize_t l() const { return l_; }
  bool is_dna() const { return dna_; }
  bool is_ambiguous() const {
    if (sketching_ == SketchingScheme::MOD)
      return (seq_pos_ < k_) || (!! last_ambig_);
    return (queue_pos_ < k_ - l_) || (!! last_ambig_);
  }

  SketchingScheme sketching_scheme() const { return sketching_; }
  Order order() const { return order_; }

  void ClearRankStats() {
    if (mod_sampler_)
      mod_sampler_->ClearRankStats();
    else if (ranker_)
      ranker_->ClearRankStats();
  }
  RankHitStats GetRankStats() const {
    if (mod_sampler_)
      return mod_sampler_->GetRankStats();
    return ranker_ ? ranker_->GetRankStats() : RankHitStats{};
  }
  bool has_uhs_stats() const {
    if (order_ != Order::GM)
      return false;
    if (mod_sampler_)
      return mod_sampler_->has_uhs_stats();
    return ranker_ && ranker_->has_uhs_stats();
  }

  void ClearPhaseStats() {
    phase_stats_ = {};
  }
  NextMinPhaseStats GetPhaseStats() const {
    return phase_stats_;
  }

  // When non-null, every ℓ-mer passed to Rank() is appended (bench).
  void SetRankedLmerLog(std::vector<uint64_t> *log) {
    ranked_lmer_log_ = log;
    if (mod_sampler_)
      mod_sampler_->SetRankedLmerLog(log);
  }
  void SetRankedTmerLog(std::vector<uint64_t> *log) {
    if (mod_sampler_)
      mod_sampler_->SetRankedTmerLog(log);
  }

  private:
  uint64_t reverse_complement(uint64_t kmer, uint8_t n);
  uint64_t canonical_representation(uint64_t kmer, uint8_t n);
  uint64_t dna_hash_key(uint64_t lmer);
  void set_lookup_table_character(char ch, uint8_t val);
  void init_ranker(const std::string &greedy_mini_filename,
                   std::shared_ptr<const GreedyMiniTable> gm_table);
  uint64_t *NextMinimizerMod();

  const std::string *str_;  // pointer to sequence
  ssize_t k_;
  ssize_t l_;
  size_t str_pos_, start_, finish_;
  uint64_t spaced_seed_mask_;
  bool dna_;
  uint64_t toggle_mask_;
  uint64_t lmer_;
  uint64_t lmer_mask_;
  uint64_t last_minimizer_;
  ssize_t loaded_ch_;
  std::deque<MinimizerData> queue_;
  ssize_t queue_pos_;
  uint64_t last_ambig_;
  uint8_t lookup_table_[UINT8_MAX + 1];
  const int revcom_version_;

  SketchingScheme sketching_;
  Order order_;
  std::unique_ptr<MinimizerRanker> ranker_;
  std::unique_ptr<ModSampler> mod_sampler_;
  ssize_t seq_pos_;  // non-ambiguous bases since last reset (MOD + shared)
  // LoadSequence start (or 0 after 'N'); absolute k-mer/window index for Mod.
  ssize_t stream_origin_ = 0;
  // Lazy cache for LoadSequence without an origin hint (tests / classify).
  const std::string *mod_origin_cache_seq_ = nullptr;
  std::vector<uint32_t> mod_origin_cache_;
  std::vector<uint64_t> *ranked_lmer_log_ = nullptr;
  // Always present so BENCH vs non-BENCH TUs share one MinimizerScanner layout
  // under -flto (RDTSC updates remain gated in mmscanner.cc).
  NextMinPhaseStats phase_stats_;
};

}

#endif
