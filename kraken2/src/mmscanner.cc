/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "mmscanner.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <limits.h>
#include <map>
#include <mutex>
#include <sys/stat.h>
#include <unistd.h>

#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
#include <x86intrin.h>
#endif

using std::string;
using std::vector;

namespace kraken2 {

namespace {

#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
inline uint64_t nm_rdtsc() {
  unsigned aux;
  return __rdtscp(&aux);
}
#endif

std::mutex gm_table_cache_mu;
std::map<std::string, std::weak_ptr<const GreedyMiniTable>> gm_table_cache;

bool same_gm_file(const std::string &a, const std::string &b) {
  if (a == b)
    return true;
  char resolved_a[PATH_MAX];
  char resolved_b[PATH_MAX];
  if (::realpath(a.c_str(), resolved_a) && ::realpath(b.c_str(), resolved_b))
    return std::strcmp(resolved_a, resolved_b) == 0;
  return false;
}

}  // namespace

void MinimizerScanner::set_lookup_table_character(char c, uint8_t val) {
  lookup_table_[(int) c] = val;
  lookup_table_[tolower(c)] = val;
}

int MinimizerScanner::gm_k_from_size(size_t table_size) {
  return ::gm_k_from_size(table_size);
}

std::vector<uint64_t> MinimizerScanner::load_gm_file(const std::string &filename) {
  auto vec = load_gm_file_raw(filename);
  if (vec.empty())
    std::cerr << "Unable to open file for reading: " << filename << std::endl;
  return vec;
}

void MinimizerScanner::write_gm_file(const std::string &filename,
                                     const std::vector<uint64_t> &ranks) {
  try {
    write_gm_file_raw(filename, ranks);
  } catch (const std::exception &ex) {
    errx(EX_OSERR, "%s", ex.what());
  }
}

std::string MinimizerScanner::gm_k2d_path(const std::string &hashtable_filename) {
  auto slash = hashtable_filename.rfind('/');
  if (slash != string::npos)
    return hashtable_filename.substr(0, slash + 1) + "gm.k2d";
  return "gm.k2d";
}

void MinimizerScanner::register_shared_gm_table(
    const std::string &filename,
    std::shared_ptr<const GreedyMiniTable> table) {
  if (filename.empty() || !table || table->empty())
    return;
  std::lock_guard<std::mutex> lock(gm_table_cache_mu);
  gm_table_cache[filename] = table;
}

std::shared_ptr<const GreedyMiniTable>
MinimizerScanner::shared_gm_table(const std::string &filename) {
  const std::string &key = filename;
  std::lock_guard<std::mutex> lock(gm_table_cache_mu);
  auto it = gm_table_cache.find(key);
  if (it != gm_table_cache.end()) {
    if (auto existing = it->second.lock())
      return existing;
  }
  std::shared_ptr<const GreedyMiniTable> loaded;
  try {
    loaded = load_gm_table(filename);
  } catch (const std::exception &ex) {
    errx(EX_DATAERR, "GreedyMini table '%s': %s",
         filename.c_str(), ex.what());
  }
  if (!loaded || loaded->empty())
    return nullptr;
  gm_table_cache[key] = loaded;
  return loaded;
}

void MinimizerScanner::embed_gm_table(
    const std::string &source_path,
    const std::string &hashtable_filename,
    std::shared_ptr<const GreedyMiniTable> table) {
  if (!table || table->empty())
    errx(EX_DATAERR, "cannot embed empty GreedyMini rank table");
  const string dest_path = gm_k2d_path(hashtable_filename);
  register_shared_gm_table(source_path, table);
  if (!same_gm_file(source_path, dest_path)) {
    write_gm_file(dest_path, gm_table_to_raw(*table));
    register_shared_gm_table(dest_path, table);
    std::cerr << "GreedyMini rank table written to " << dest_path
              << " (gm_k=" << table->gm_k()
              << ", UHS=" << table->uhs_size()
              << "/" << table->size() << ")" << std::endl;
  }
  else {
    register_shared_gm_table(dest_path, table);
    std::cerr << "GreedyMini rank table already at " << dest_path << std::endl;
  }
}

void MinimizerScanner::init_ranker(
    const std::string &greedy_mini_filename,
    std::shared_ptr<const GreedyMiniTable> gm_table) {
  RankerConfig cfg;
  cfg.order = order_;
  cfg.l = (int)l_;
  cfg.bpc = dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO;
  cfg.toggle_mask = toggle_mask_;
  cfg.lmer_mask = lmer_mask_;
  cfg.dna = dna_;
  cfg.spaced_seed_mask = spaced_seed_mask_;
  cfg.revcom_version = revcom_version_;
  cfg.greedy_mini_filename = greedy_mini_filename;
  cfg.gm_table = gm_table;
  cfg.gm_flip = gm_flip_;

  if (order_ == Order::GM && !cfg.gm_table) {
    if (greedy_mini_filename.empty())
      errx(EX_USAGE, "GreedyMini order requires a rank file (.gm)");
    cfg.gm_table = shared_gm_table(greedy_mini_filename);
    if (!cfg.gm_table || cfg.gm_table->empty())
      errx(EX_NOINPUT, "failed to load GreedyMini rank file: %s",
           greedy_mini_filename.c_str());
  }

  if (sketching_ == SketchingScheme::MOD) {
    try {
      mod_sampler_.reset(new ModSampler(
          (int)k_, (int)l_, cfg.bpc, toggle_mask_, dna_, order_,
          revcom_version_, spaced_seed_mask_, cfg.gm_table, gm_flip_));
    } catch (const std::exception &ex) {
      errx(EX_DATAERR, "ModSampler: %s", ex.what());
    }
    ranker_.reset();
    return;
  }

  try {
    ranker_ = make_order_ranker(cfg);
  } catch (const std::exception &ex) {
    errx(EX_DATAERR, "%s", ex.what());
  }
  gm_rf_ = nullptr;
  gm_ = nullptr;
  if (order_ == Order::GM && dna_) {
    gm_ = dynamic_cast<GMRanker *>(ranker_.get());
    if (gm_ && gm_->rank_first()) {
      sigma_mask_ = gm_->sigma_mask();
      SetRankFirst(true, kDefaultRollSigma);
    }
  }
}

// Test / benchmark hook: rank_first = false takes the Rank() path (full packed
// value of every l-mer, the reference behaviour); roll_sigma = true keeps the
// Sigma plane by one shift per base instead of one pext per l-mer. All
// combinations select the same minimizers. No-op unless the ranker supports
// rank-first comparison (GM centre order, revcom version 1).
void MinimizerScanner::SetRankFirst(bool rank_first, bool roll_sigma) {
  gm_rf_ = (gm_ && gm_->rank_first() && rank_first) ? gm_ : nullptr;
  roll_sigma_ = gm_rf_ && roll_sigma;
}

MinimizerScanner::MinimizerScanner(ssize_t k, ssize_t l,
    uint64_t spaced_seed_mask, bool dna_sequence, uint64_t toggle_mask,
    int revcom_version)
    : MinimizerScanner(k, l, spaced_seed_mask, dna_sequence, toggle_mask,
                       revcom_version, SketchingScheme::MINIMIZER,
                       Order::RANDOM_XOR, "")
{}

MinimizerScanner::MinimizerScanner(ssize_t k, ssize_t l,
    uint64_t spaced_seed_mask, bool dna_sequence, uint64_t toggle_mask,
    int revcom_version, SketchingScheme sketching, Order order,
    std::string greedy_mini_filename,
    std::shared_ptr<const GreedyMiniTable> gm_table, bool gm_flip)
    : gm_flip_(gm_flip), str_(nullptr), k_(k), l_(l), str_pos_(0), start_(0), finish_(0),
      spaced_seed_mask_(spaced_seed_mask), dna_(dna_sequence),
      toggle_mask_(toggle_mask), loaded_ch_(0),
      last_ambig_(0), revcom_version_(revcom_version),
      sketching_(sketching), order_(order), seq_pos_(0)
{
  if (l_ > (ssize_t) ((sizeof(uint64_t) * 8 - 1) / (dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO)))
    errx(EX_SOFTWARE, "l exceeds size limits for minimizer %s scanner",
         dna_ ? "nucleotide" : "protein");
  lmer_mask_ = 1;
  lmer_mask_ <<= (l_ * (dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO));
  lmer_mask_--;
  toggle_mask_ &= lmer_mask_;
  if (finish_ == SIZE_MAX)
    finish_ = str_->size();
  if ((ssize_t) (finish_ - start_) + 1 < l_)  // Invalidate scanner if interval < 1 l-mer
    str_pos_ = finish_;
  for (int i = 0; i < UINT8_MAX + 1; i++)
    lookup_table_[i] = UINT8_MAX;
  if (dna_) {
    set_lookup_table_character('A', 0x00);
    set_lookup_table_character('C', 0x01);
    set_lookup_table_character('G', 0x02);
    set_lookup_table_character('T', 0x03);
  }
  else {
    // Reduced alphabet uses 15-letter alphabet from AD Solis (2015),
    // in Proteins (doi:10.1002/prot.24936)
    // This means that the ambiguous codes B (N/D) and Z (E/Q) aren't
    // usable here because they span multiple groups.  Although J (I/L)
    // could be used because I & L are in the same group, we treat it as
    // we do B, Z, and X for consistency and to ease future changes to
    // the code base.

    // stop codons/rare amino acids
    set_lookup_table_character('*', 0x00);
    set_lookup_table_character('U', 0x00);
    set_lookup_table_character('O', 0x00);
    // alanine
    set_lookup_table_character('A', 0x01);
    // asparagine, glutamine, serine
    set_lookup_table_character('N', 0x02);
    set_lookup_table_character('Q', 0x02);
    set_lookup_table_character('S', 0x02);
    // cysteine
    set_lookup_table_character('C', 0x03);
    // aspartic acid, glutamic acid
    set_lookup_table_character('D', 0x04);
    set_lookup_table_character('E', 0x04);
    // phenylalanine
    set_lookup_table_character('F', 0x05);
    // glycine
    set_lookup_table_character('G', 0x06);
    // histidine
    set_lookup_table_character('H', 0x07);
    // isoleucine, leucine
    set_lookup_table_character('I', 0x08);
    set_lookup_table_character('L', 0x08);
    // lysine
    set_lookup_table_character('K', 0x09);
    // proline
    set_lookup_table_character('P', 0x0a);
    // arginine
    set_lookup_table_character('R', 0x0b);
    // methionine, valine
    set_lookup_table_character('M', 0x0c);
    set_lookup_table_character('V', 0x0c);
    // threonine
    set_lookup_table_character('T', 0x0d);
    // tryptophan
    set_lookup_table_character('W', 0x0e);
    // tyrosine
    set_lookup_table_character('Y', 0x0f);
  }

  init_ranker(greedy_mini_filename, gm_table);
  if (sketching_ == SketchingScheme::MOD && mod_sampler_)
    scan_kind_ = kScanMod;
#if !defined(KRAKEN2_MINIMIZER_BENCHMARK)
  else if (sketching_ == SketchingScheme::MINIMIZER &&
           order_ == Order::RANDOM_XOR && dna_)
    scan_kind_ = kScanXor;
#endif
}

void MinimizerScanner::ComputeModStreamOrigins(
    const string &seq, std::vector<uint32_t> *origin_at) const {
  origin_at->resize(seq.size() + 1);
  uint32_t o = 0;
  for (size_t i = 0; i < seq.size(); i++) {
    (*origin_at)[i] = o;
    if (lookup_table_[(unsigned char)seq[i]] == UINT8_MAX)
      o = 0;
    else
      ++o;
  }
  (*origin_at)[seq.size()] = o;
}

void MinimizerScanner::LoadSequence(const string &seq, size_t start,
                                    size_t finish) {
  LoadSequence(seq, start, finish, (ssize_t)-1);
}

void MinimizerScanner::LoadSequence(const string &seq, size_t start,
                                    size_t finish, ssize_t mod_stream_origin) {
  str_ = &seq;
  start_ = start;
  finish_ = finish;
  str_pos_ = start_;
  if (finish_ > str_->size())
    finish_ = str_->size();
  if ((ssize_t) (finish_ - start_) + 1 < l_)  // Invalidate scanner if interval < 1 l-mer
    str_pos_ = finish_;
  queue_.clear();
  xq_.clear();
  queue_pos_ = 0;
  loaded_ch_ = 0;
  last_minimizer_ = ~0;
  last_ambig_ = 0;
  lmer_ = 0;
  sigma_ = 0;
  seq_pos_ = 0;
  stream_origin_ = 0;
  if (ranker_)
    ranker_->Reset();
  if (mod_sampler_) {
    if (mod_stream_origin >= 0) {
      stream_origin_ = mod_stream_origin;
    } else {
      // One O(n) table per distinct sequence string; O(1) per later chunk.
      if (mod_origin_cache_seq_ != str_ ||
          mod_origin_cache_.size() != str_->size() + 1) {
        ComputeModStreamOrigins(*str_, &mod_origin_cache_);
        mod_origin_cache_seq_ = str_;
      }
      stream_origin_ = (ssize_t)mod_origin_cache_[start_];
    }
    mod_sampler_->ResetAt(stream_origin_);
  }
}

__attribute__((noinline)) uint64_t *MinimizerScanner::NextMinimizerMod() {
  if (str_pos_ >= finish_)
    return nullptr;

  const auto bits_per_char = dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO;
  const auto ambig_code = (1u << bits_per_char) - 1;

  // Batch warm-up like classic Minimizer: do not return until a full k-mer
  // window is ready (seq_pos_ >= k_). Ambiguous bases still return immediately
  // so classify/build see TOK_AMBIG / is_ambiguous() at the N.
  while (true) {
    if (str_pos_ >= finish_)
      return nullptr;

#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_enc0 = nm_rdtsc();
#endif
    const auto lookup_code = lookup_table_[(int)(*str_)[str_pos_++]];
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_enc1 = nm_rdtsc();
#endif

    // Shift ambig register every base (same as Minimizer path).
    last_ambig_ <<= bits_per_char;
    if (lookup_code == UINT8_MAX) {
      last_ambig_ |= ambig_code;
      last_ambig_ &= lmer_mask_;
      seq_pos_ = 0;
      stream_origin_ = 0;  // post-N stream is 0-relative (matches full scan)
      mod_sampler_->OnAmbiguous();
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
      phase_stats_.steps++;
      phase_stats_.encode_cycles += (t_enc1 - t_enc0);
#endif
      return &last_minimizer_;
    }
    last_ambig_ &= lmer_mask_;

    seq_pos_++;

    uint64_t sampled = last_minimizer_;
    // Absolute window index so residue/ring state matches a contig-long scan.
    const ssize_t win =
        seq_pos_ >= k_ ? stream_origin_ + seq_pos_ - k_ : (ssize_t)-1;
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    mod_sampler_->PushBase(lookup_code, nullptr);
    const uint64_t t_rank1 = nm_rdtsc();
    if (win >= 0 && mod_sampler_->FinishWindow(win, &sampled))
      last_minimizer_ = sampled;
    const uint64_t t_deq1 = nm_rdtsc();
    phase_stats_.steps++;
    phase_stats_.encode_cycles += (t_enc1 - t_enc0);
    phase_stats_.rank_cycles += (t_rank1 - t_enc1);
    phase_stats_.deque_cycles += (t_deq1 - t_rank1);
#else
    if (mod_sampler_->PushAndFinish(lookup_code, win, &sampled))
      last_minimizer_ = sampled;
#endif
    if (seq_pos_ >= k_)
      return &last_minimizer_;
  }
}

// Pop the queue back while it is worse than the l-mer lmer_ with table rank
// `rank` (the leftmost minimum stays on ties). Equal finite ranks compare the
// full packed values, which is what Rank() would have returned; infinite
// ranks are all equal.
void MinimizerScanner::PopDominatedRankFirst(uint64_t rank) {
  uint64_t value = 0;
  bool have_value = false;
  while (! queue_.empty()) {
    const uint64_t back = queue_.back().rank;
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    rank_cmp_++;
#endif
    if (back < rank)
      break;
    if (back == rank) {
      if (rank == GMRanker::kInfRank)
        break;
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
      tie_cmp_++;
#endif
      if (! have_value) {
        value = gm_rf_->value(lmer_);
        have_value = true;
      }
      if (gm_rf_->value(queue_.back().candidate) <= value)
        break;
    }
    queue_.pop_back();
  }
}

// Default path (Minimizer sketching, RandomXor order, DNA): the original
// Kraken 2 scan loop (16-byte queue elements, no virtual ranker call, no GM
// or Mod bookkeeping). Returns exactly what the general loop below returns.
__attribute__((noinline)) uint64_t *MinimizerScanner::NextMinimizerXor() {
  if (str_pos_ >= finish_)
    return nullptr;
  bool changed_minimizer = false;
  while (! changed_minimizer) {
    if (loaded_ch_ == l_)
      loaded_ch_--;
    while (loaded_ch_ < l_ && str_pos_ < finish_) {
      loaded_ch_++;
      lmer_ <<= BITS_PER_CHAR_DNA;
      last_ambig_ <<= BITS_PER_CHAR_DNA;
      auto lookup_code = lookup_table_[ (int) (*str_)[str_pos_++] ];
      if (lookup_code == UINT8_MAX) {
        xq_.clear();
        queue_pos_ = 0;
        lmer_ = 0;
        loaded_ch_ = 0;
        last_ambig_ |= (1u << BITS_PER_CHAR_DNA) - 1;
      }
      else
        lmer_ |= lookup_code;
      lmer_ &= lmer_mask_;
      last_ambig_ &= lmer_mask_;
      if ((str_pos_ - start_) >= (size_t) k_ && loaded_ch_ < l_)
        return &last_minimizer_;
    }
    if (loaded_ch_ < l_)
      return nullptr;
    uint64_t canonical_lmer = canonical_representation(lmer_, l_);
    if (spaced_seed_mask_)
      canonical_lmer &= spaced_seed_mask_;
    uint64_t candidate_lmer = canonical_lmer ^ toggle_mask_;
    if (k_ == l_) {
      last_minimizer_ = canonical_lmer;
      return &last_minimizer_;
    }
    while (! xq_.empty() && xq_.back().candidate > candidate_lmer)
      xq_.pop_back();
    XorMinimizerData data = { candidate_lmer, queue_pos_ };
    if (xq_.empty() && queue_pos_ >= k_ - l_)
      changed_minimizer = true;
    xq_.push_back(data);
    if (xq_.front().pos < queue_pos_ - k_ + l_) {
      xq_.pop_front();
      changed_minimizer = true;
    }
    if (queue_pos_ == k_ - l_)
      changed_minimizer = true;
    queue_pos_++;
    if (str_pos_ >= (size_t) k_)
      break;
  }
  last_minimizer_ = xq_.front().candidate ^ toggle_mask_;
  return &last_minimizer_;
}

__attribute__((noinline)) uint64_t *MinimizerScanner::NextMinimizerGeneral() {
  if (sketching_ == SketchingScheme::MOD && mod_sampler_)
    return NextMinimizerMod();
  if (str_pos_ >= finish_)  // Abort if we've exhausted string interval
    return nullptr;
  bool changed_minimizer = false;
  auto bits_per_char = dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO;
  auto ambig_code = (1u << bits_per_char) - 1;
  while (! changed_minimizer) {
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_enc0 = nm_rdtsc();
#endif
    // Incorporate next character (and more if needed to fill l-mer)
    if (loaded_ch_ == l_)
      loaded_ch_--;
    while (loaded_ch_ < l_ && str_pos_ < finish_) {  // char loading loop
      loaded_ch_++;
      lmer_ <<= bits_per_char;
      last_ambig_ <<= bits_per_char;
      auto lookup_code = lookup_table_[ (int) (*str_)[str_pos_++] ];
      if (lookup_code == UINT8_MAX) {
        queue_.clear();
        queue_pos_ = 0;
        lmer_ = 0;
        sigma_ = 0;
        loaded_ch_ = 0;
        last_ambig_ |= ambig_code;
        seq_pos_ = 0;
        if (ranker_)
          ranker_->Reset();
      }
      else {
        lmer_ |= lookup_code;
        if (roll_sigma_)
          sigma_ = gm_centre::roll_in_low(sigma_, gm_rf_->sigma_bit(lookup_code),
                                          sigma_mask_);
        seq_pos_++;
      }
      lmer_ &= lmer_mask_;
      last_ambig_ &= lmer_mask_;
      // If we haven't filled up first k-mer, don't return
      // Otherwise, if l-mer is incomplete, still return
      // If l-mer is complete, go through queue management
      // If we don't have a full l-mer, there's no point in
      // queue management for minimizer calculation
      if ((str_pos_ - start_) >= (size_t) k_ && loaded_ch_ < l_)
        return &last_minimizer_;
    }  // end character loading loop
    if (loaded_ch_ < l_)  // Abort if we've exhausted string interval w/o
      return nullptr;     //   filling the l-mer
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_enc1 = nm_rdtsc();
    const uint64_t t_can0 = t_enc1;
#endif
    if (k_ == l_) {  // Short-circuit queue work
      last_minimizer_ = dna_hash_key(lmer_);
      return &last_minimizer_;
    }

    // RandomXor/Pcg: DNA-lex every ℓ-mer (rank = canonical ^ toggle).
    // GM: Rank the raw ℓ-mer (gm_k orientation is in the ranker); DNA-lex
    //     only the window winner for the hash key.
    uint64_t rank_arg = lmer_;
    uint64_t candidate = lmer_;
    if (order_ != Order::GM) {
      candidate = dna_hash_key(lmer_);
      rank_arg = candidate;
    }
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_can1 = nm_rdtsc();
#endif

    // Sliding window minimum calculation (MINIMIZER sketching)
    if (ranked_lmer_log_)
      ranked_lmer_log_->push_back(rank_arg);
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_rank0 = nm_rdtsc();
#endif
    // GM centre order, rank first: the rank is one table load of the Sigma
    // plane (rolled, or one pext of lmer_); the full packed value is only
    // computed when a finite rank ties with the back of the queue.
    const uint64_t rank =
        gm_rf_ ? gm_rf_->rank_of_sigma(roll_sigma_ ? sigma_
                                                    : gm_rf_->sigma_of(lmer_))
               : ranker_->Rank(rank_arg);
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_rank1 = nm_rdtsc();
    const uint64_t t_deq0 = t_rank1;
#endif
    if (gm_rf_)
      PopDominatedRankFirst(rank);
    else
      while (! queue_.empty() && queue_.back().rank > rank)
        queue_.pop_back();
    MinimizerData data = { rank, candidate, queue_pos_ };
    if (queue_.empty() && queue_pos_ >= k_ - l_) {
      // Empty queue means front will change
      // Minimizer will change iff we've processed enough l-mers
      changed_minimizer = true;
    }
    queue_.push_back(data);
    // expire an l-mer not in the current window
    if (queue_.front().pos < queue_pos_ - k_ + l_) {
      // queue_.erase(queue_.begin());
      queue_.pop_front();
      // Change in front means minimizer changed
      changed_minimizer = true;
    }

    // Change from no minimizer (beginning of sequence/near ambig. char)
    if (queue_pos_ == k_ - l_)
      changed_minimizer = true;
    queue_pos_++;
#if defined(KRAKEN2_MINIMIZER_BENCHMARK)
    const uint64_t t_deq1 = nm_rdtsc();
    phase_stats_.steps++;
    phase_stats_.encode_cycles += (t_enc1 - t_enc0);
    phase_stats_.canonical_cycles += (t_can1 - t_can0);
    phase_stats_.rank_cycles += (t_rank1 - t_rank0);
    phase_stats_.deque_cycles += (t_deq1 - t_deq0);
#endif

    // Return only if we've read in at least one k-mer's worth of chars
    if (str_pos_ >= (size_t) k_) {
      break;
    }
  }  // end while ! changed_minimizer
  assert(! queue_.empty());
  if (order_ == Order::GM && dna_) {
    if (changed_minimizer)
      last_minimizer_ = dna_hash_key(queue_.front().candidate);
  } else {
    last_minimizer_ = queue_.front().candidate;
  }
  return &last_minimizer_;
}

// Adapted for 64-bit DNA use from public domain code at:
//   https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
uint64_t MinimizerScanner::reverse_complement(uint64_t kmer, uint8_t n) {
  return dna_reverse_complement(kmer, n, revcom_version_);
}

uint64_t MinimizerScanner::canonical_representation(uint64_t kmer, uint8_t n) {
  return dna_canonical_representation(kmer, n, revcom_version_);
}

uint64_t MinimizerScanner::dna_hash_key(uint64_t lmer) {
  // DNA-lex (scanner revcom_version), then spaced AND. Used as RandomXor/Pcg
  // rank input and as the returned hash key (including GM window winners).
  uint64_t key = dna_ ? canonical_representation(lmer, (uint8_t)l_) : lmer;
  if (spaced_seed_mask_)
    key &= spaced_seed_mask_;
  return key;
}

} // end namespace

// using namespace kraken2;

// int main(void) {
//   MinimizerScanner mm(4, 3, DEFAULT_SPACED_SEED_MASK, true, 0);
//   // std::string sequence =
//   // "ACCAGCAGGGCGCACGCATCGCGCTGCCCGAGCCCTCGACCCCCTTTCTCTGTGAGCAGCAGCCCGATTGGCCGTTGAAGGCACCACAACAACCCGACTGGCCGCTGCCCAGGCCGGCGGCGCCAACCCGGCCTGCGGCGCGGCCCCAGAGCCTGCCTTATTTCTCCCCCTACCGACCATCTCCGCAGCATGCCAGCCGTCCTCTCGCGATAGGGGACCGCCACAATAGCAACGCCAATCACCTTCCCACCCCAGCCACCAGCAGCAGCAGCACGCTCTCAGATGACGCCCTCTGGCAGGCCCTCGGCGGCGGGCGCCGCCCCTCGTCGTCATCGGCGGCGACGAACAACAACAGCAGCAACAACAACAGGCTCAGCTACCACCCAGCAAACGGGGAAGCTCGCCGCCTGCATCCAGCTACTGCGACGGCCACGGCGAAGGCGAAGGCGACGGCAGAACACCAGGCCCATCCCGACTACCCGCATATGAGCCCCTATGACGACGACAACGCCGACGCCGACGCCGACGCCGATAACAACGGCGATGTGGGTGGCGGGACTGGTGTGGTGGGCCTGGGCGTCGTGATGGCGCCGCGCCCGGCGAGCATGTTGCCGGCGCGGTTGCGGACGCGTGATCGGGCTCCGCCGGCGCTGCCGCTGCTGCCTGCGACGACGCTGGCGTCGGTGTCGTTGGGGCCGACGTTGAGTGGCCCGTTTGTTGCGCAGTTTGAGTGAGCGGTGGTTTGTCTTGNNNNNNNNNNNNNNNNCTTGCTTGTCTAGAGGGGGTGTGGATTGGGAAGTGGTGGCGGGGAAGAAGGGCGGGCGTTGGTTGATGGGACAGGACGAAGGGCTGGGTTTCATTATGAACCAATTGTCCACGGTTGGGTAAGATATAGAGGAGCAAAGCGCATTGCATATGATACCCCTGGATAAATGAGAGCAGCAATGATCGACTGAATCCTTCTGATGTCGGTGTCCGGAGAAAGGGATCGGTTCTGAAC";
//   std::string sequence = "GTCATGCACGTTCAC";
//   mm.LoadSequence(sequence);
//   size_t mm_count = 0;
//   while (mm.NextMinimizer()) {
//     mm_count += 1;
//   }

//   std::cout << "minimizers found = " << mm_count << std::endl;

//   return 0;
// }
