/*
 * Correctness tests for MOD streaming sampler (ModSampler + MinimizerScanner).
 */

#include "dna_kmer.h"
#include "mmscanner.h"
#include "mod_minimizer.h"
#include "utilities.h"

#include <deque>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

using namespace kraken2;

static int failures = 0;

static void check(bool cond, const char *msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << "\n";
    failures++;
  }
}

static uint64_t xor_rank(uint64_t canonical, uint64_t toggle_mask) {
  return canonical ^ toggle_mask;
}

// Slow reference MOD stream (no fast-path LUTs/counters).
class NaiveModReference {
 public:
  NaiveModReference(int k, int l, uint64_t toggle_mask,
                    uint64_t spaced_seed_mask = 0,
                    int revcom_version = CURRENT_REVCOM_VERSION,
                    std::shared_ptr<const GreedyMiniTable> gm_table = nullptr,
                    bool mask_tmer_rank = false)
      : k_(k), l_(l), w_(k - l + 1),
        t_(ModSampler::kDefaultR + ((l - ModSampler::kDefaultR) % w_)),
        tmer_span_(w_ + l - t_), toggle_(toggle_mask),
        revcom_version_(revcom_version), mask_tmer_rank_(mask_tmer_rank) {
    const int t_bits = 2 * t_;
    t_mask_ = (t_bits >= 64) ? ~0ULL : ((1ULL << t_bits) - 1ULL);
    toggle_t_ = toggle_ & t_mask_;
    const int l_bits = 2 * l_;
    lmer_mask_ = (l_bits >= 64) ? ~0ULL : ((1ULL << l_bits) - 1ULL);
    spaced_seed_mask_ = spaced_seed_mask & lmer_mask_;
    deques_.assign((size_t)w_, {});
    lmer_ring_.assign((size_t)w_, 0);
    if (gm_table) {
      const int t_bits = 2 * t_;
      const uint64_t t_mask =
          (t_bits >= 64) ? ~0ULL : ((1ULL << t_bits) - 1ULL);
      t_ranker_.reset(new GMRanker(t_, 2, toggle_t_, t_mask, true, gm_table, 0,
                                   revcom_version));
    }
    reset();
  }

  void reset() {
    for (auto &dq : deques_)
      dq.clear();
    t_bases_.clear();
    l_bases_.clear();
    abs_base_ = 0;
    abs_tmer_ = 0;
    lmers_stored_ = 0;
    win_start_ = -1;
  }

  void push_base(uint8_t code) {
    abs_base_++;
    t_bases_.push_back(code);
    if ((int)t_bases_.size() > t_)
      t_bases_.pop_front();
    if ((int)t_bases_.size() == t_) {
      const uint64_t fwd = encode(t_bases_, t_);
      uint64_t canonical =
          dna_canonical_representation(fwd, (uint8_t)t_, revcom_version_);
      if (mask_tmer_rank_ && spaced_seed_mask_)
        canonical &= spaced_seed_mask_;
      const uint64_t rank =
          t_ranker_ ? t_ranker_->Rank((mask_tmer_rank_ && spaced_seed_mask_)
                                          ? (fwd & spaced_seed_mask_)
                                          : fwd)
                    : xor_rank(canonical, toggle_t_);
      const int a = (int)(abs_tmer_ % w_);
      push_monotone(deques_[(size_t)a], rank, abs_tmer_);
      abs_tmer_++;
    }

    l_bases_.push_back(code);
    if ((int)l_bases_.size() > l_)
      l_bases_.pop_front();
    if ((int)l_bases_.size() == l_) {
      const uint64_t fwd = encode(l_bases_, l_);
      uint64_t canonical =
          dna_canonical_representation(fwd, (uint8_t)l_, revcom_version_);
      if (spaced_seed_mask_)
        canonical &= spaced_seed_mask_;
      lmer_ring_[(size_t)((abs_base_ - l_) % w_)] = canonical;
      if (lmers_stored_ < w_)
        lmers_stored_++;
    }
  }

  bool finish_window(ssize_t win, uint64_t *out) {
    if (win < 0 || lmers_stored_ < w_)
      return false;
    if (abs_tmer_ < win + tmer_span_)
      return false;

    for (int a = 0; a < w_; a++)
      expire_before(deques_[(size_t)a], win);

    const int wm = (int)(win % w_);
    uint64_t best_rank = UINT64_MAX;
    int best_p = 0;
    bool any = false;
    for (int a = 0; a < w_; a++) {
      if (deques_[(size_t)a].empty())
        continue;
      const uint64_t r = deques_[(size_t)a].front().rank;
      // Match ModSampler: non-UHS GM ranks (UINT64_MAX) are not candidates.
      if (r == UINT64_MAX)
        continue;
      const int p = residue_pos(a, wm);
      if (!any || r < best_rank || (r == best_rank && p < best_p)) {
        best_rank = r;
        best_p = p;
        any = true;
      }
    }
    if (!any)
      return false;
    *out = lmer_ring_[(size_t)((wm + best_p) % w_)];
    win_start_ = win;
    return true;
  }

 private:
  struct Entry {
    uint64_t rank;
    ssize_t pos;
  };

  static uint64_t encode(const std::deque<uint8_t> &bases, int n) {
    uint64_t v = 0;
    for (int i = 0; i < n; i++)
      v = (v << 2) | (uint64_t)bases[(size_t)i];
    return v;
  }

  static void push_monotone(std::deque<Entry> &dq, uint64_t rank, ssize_t pos) {
    while (!dq.empty() && dq.back().rank > rank)
      dq.pop_back();
    dq.push_back(Entry{rank, pos});
  }

  static void expire_before(std::deque<Entry> &dq, ssize_t win) {
    while (!dq.empty() && dq.front().pos < win)
      dq.pop_front();
  }

  int residue_pos(int a, int win_mod) const {
    int p = a - win_mod;
    return p < 0 ? p + w_ : p;
  }

  int k_, l_, w_, t_, tmer_span_;
  uint64_t toggle_, t_mask_, toggle_t_, lmer_mask_, spaced_seed_mask_;
  int revcom_version_;
  bool mask_tmer_rank_;
  std::unique_ptr<MinimizerRanker> t_ranker_;
  std::vector<std::deque<Entry>> deques_;
  std::deque<uint8_t> t_bases_;
  std::deque<uint8_t> l_bases_;
  std::vector<uint64_t> lmer_ring_;
  ssize_t abs_base_;
  ssize_t abs_tmer_;
  int lmers_stored_;
  ssize_t win_start_;
};

static uint8_t dna_lut(char ch) {
  switch (ch) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      return 255;
  }
}
static std::vector<uint64_t> scan_mod_minimizers(
    const std::string &dna, int k, int l, uint64_t toggle,
    uint64_t spaced = 0, int revcom_version = CURRENT_REVCOM_VERSION,
    Order order = Order::RANDOM_XOR,
    std::shared_ptr<const GreedyMiniTable> gm_table = nullptr) {
  MinimizerScanner sc(k, l, spaced, true, toggle, revcom_version,
                      SketchingScheme::MOD, order, "", gm_table);
  sc.LoadSequence(dna);
  std::vector<uint64_t> out;
  // After warm-up batching: ~L-k+1 returns on clean DNA (plus ambig returns).
  out.reserve(dna.size() >= (size_t)k ? dna.size() - (size_t)k + 1 : 0);
  while (uint64_t *mm = sc.NextMinimizer())
    out.push_back(*mm);
  return out;
}

static std::vector<uint64_t> naive_mod_stream(
    const std::string &dna, int k, int l, uint64_t toggle,
    uint64_t spaced = 0, int revcom_version = CURRENT_REVCOM_VERSION,
    std::shared_ptr<const GreedyMiniTable> gm_table = nullptr,
    bool mask_tmer_rank = false) {
  NaiveModReference ref(k, l, toggle, spaced, revcom_version, gm_table,
                        mask_tmer_rank);
  std::vector<uint64_t> out;
  out.reserve(dna.size() >= (size_t)k ? dna.size() - (size_t)k + 1 : 0);
  ssize_t seq_pos = 0;
  uint64_t last = UINT64_MAX;
  for (char ch : dna) {
    const uint8_t code = dna_lut(ch);
    if (code > 3) {
      ref.reset();
      seq_pos = 0;
      // Match NextMinimizerMod: one return at the ambiguous base.
      out.push_back(last);
      continue;
    }
    ref.push_base(code);
    seq_pos++;
    uint64_t sampled = last;
    const ssize_t win = seq_pos >= k ? seq_pos - k : (ssize_t)-1;
    if (win >= 0 && ref.finish_window(win, &sampled))
      last = sampled;
    // Match warm-up batching: emit only once a full k-mer window is ready.
    if (seq_pos >= k)
      out.push_back(last);
  }
  return out;
}

static void test_naive_matches_scanner(
    int k, int l, const std::string &dna, uint64_t toggle, const char *label,
    uint64_t spaced = 0, int revcom_version = CURRENT_REVCOM_VERSION,
    Order order = Order::RANDOM_XOR,
    std::shared_ptr<const GreedyMiniTable> gm_table = nullptr) {
  const auto fast = scan_mod_minimizers(dna, k, l, toggle, spaced,
                                        revcom_version, order, gm_table);
  const auto naive = naive_mod_stream(dna, k, l, toggle, spaced,
                                      revcom_version, gm_table);
  check(fast.size() == naive.size(),
        (std::string(label) + " size").c_str());
  if (fast.size() == naive.size()) {
    for (size_t i = 0; i < fast.size(); i++) {
      if (fast[i] != naive[i]) {
        std::cerr << "FAIL: " << label << " mismatch at " << i << " fast=0x"
                  << std::hex << fast[i] << " naive=0x" << naive[i] << std::dec
                  << "\n";
        failures++;
        break;
      }
    }
  }
}

// Default NT spaced seed for l=31, s=7 after ExpandSpacedSeedMask(×2).
static constexpr uint64_t kDefaultNtSpacedMask = 0x3ffffffff3333333ULL;

static void test_mod_spaced_seed_masking() {
  const int k = 35, l = 31;
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  const uint64_t spaced = kDefaultNtSpacedMask;
  std::mt19937_64 rng(123);
  std::string dna(5000, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];

  test_naive_matches_scanner(k, l, dna, toggle, "spaced naive", spaced);

  MinimizerScanner sc(k, l, spaced, true, toggle, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  sc.LoadSequence(dna);
  size_t checked = 0;
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (sc.is_ambiguous())
      continue;
    check((*mm & spaced) == *mm, "returned minimizer respects spaced mask");
    checked++;
  }
  check(checked > 100, "spaced seed produced samples");

  // Masked stream must differ from unmasked on random DNA.
  const auto plain = scan_mod_minimizers(dna, k, l, toggle, 0);
  const auto masked = scan_mod_minimizers(dna, k, l, toggle, spaced);
  check(plain.size() == masked.size(), "plain vs masked size");
  bool differ = false;
  for (size_t i = 0; i < plain.size() && i < masked.size(); i++) {
    if (plain[i] != masked[i]) {
      differ = true;
      break;
    }
  }
  check(differ, "spaced mask changes some minimizer values");
}

static void test_mod_sampler_direct(int k, int l, const std::string &dna,
                                    uint64_t toggle) {
  ModSampler sampler(k, l, 2, toggle, true, Order::RANDOM_XOR,
                     CURRENT_REVCOM_VERSION);

  std::vector<uint64_t> from_sampler;
  ssize_t seq_pos = 0;
  uint64_t last_direct = UINT64_MAX;
  for (char ch : dna) {
    const uint8_t code = dna_lut(ch);
    if (code > 3) {
      sampler.OnAmbiguous();
      seq_pos = 0;
      continue;
    }
    seq_pos++;
    const ssize_t win = seq_pos >= k ? seq_pos - k : (ssize_t)-1;
    uint64_t sampled = last_direct;
    if (sampler.PushAndFinish(code, win, &sampled)) {
      last_direct = sampled;
      from_sampler.push_back(sampled);
    }
  }

  MinimizerScanner sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  sc.LoadSequence(dna);
  std::vector<uint64_t> from_scanner;
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (!sc.is_ambiguous())
      from_scanner.push_back(*mm);
  }

  check(from_sampler.size() == from_scanner.size(),
        "ModSampler vs scanner size");
  if (from_sampler.size() == from_scanner.size()) {
    for (size_t i = 0; i < from_sampler.size(); i++) {
      if (from_sampler[i] != from_scanner[i]) {
        std::cerr << "FAIL: ModSampler vs scanner mismatch at " << i << "\n";
        failures++;
        break;
      }
    }
  }
}

static void test_incremental_rc_equals_canonical(int l) {
  std::mt19937_64 rng(99);
  uint64_t fwd = 0, rc = 0;
  const uint64_t mask =
      (2 * l >= 64) ? ~0ULL : ((1ULL << (2 * l)) - 1ULL);
  const int rc_shift = 2 * (l - 1);
  for (int step = 0; step < 5000; step++) {
    const uint64_t base = rng() & 3ULL;
    const uint64_t rev = base ^ 3ULL;
    fwd = ((fwd << 2) | base) & mask;
    rc = (rc >> 2) | (rev << rc_shift);
    rc &= mask;
    if (step + 1 < l)
      continue;  // ModSampler compares canonical only once l bases are loaded.
    const uint64_t inc = fwd < rc ? fwd : rc;
    const uint64_t full =
        dna_canonical_representation(fwd, (uint8_t)l, CURRENT_REVCOM_VERSION);
    if (inc != full) {
      std::cerr << "FAIL: incremental RC l=" << l << " step=" << step
                << " inc=0x" << std::hex << inc << " full=0x" << full
                << std::dec << "\n";
      failures++;
      break;
    }
  }
}

static void test_lut_ring_index(int w) {
  for (int wm = 0; wm < w; wm++) {
    for (int p = 0; p < w; p++) {
      const int expect = (wm + p) % w;
      int got = p;
      if (w == 5) {
        uint8_t lut[16][16] = {};
        for (int wi = 0; wi < w; wi++) {
          for (int a = 0; a < w; a++) {
            int pos = a - wi;
            if (pos < 0)
              pos += w;
            lut[wi][pos] = (uint8_t)((wi + pos) % w);
          }
        }
        got = lut[wm][p];
      } else {
        got = (wm + p) % w;
      }
      if (got != expect) {
        std::cerr << "FAIL: ring_ix_lut wm=" << wm << " p=" << p << "\n";
        failures++;
      }
    }
  }
}

static size_t count_mod_unique(const std::string &dna, int k, int l,
                               uint64_t toggle, size_t chunk = 0) {
  std::unordered_set<uint64_t> uniq;
  MinimizerScanner sc(k, l, DEFAULT_SPACED_SEED_MASK, true, toggle,
                      CURRENT_REVCOM_VERSION, SketchingScheme::MOD,
                      Order::RANDOM_XOR, "", nullptr);
  auto collect = [&](ssize_t start, ssize_t finish) {
    sc.LoadSequence(dna, start, finish);
    uint64_t *m;
    while ((m = sc.NextMinimizer())) {
      if (!sc.is_ambiguous())
        uniq.insert(*m);
    }
  };
  if (!chunk || chunk >= dna.size())
    collect(0, (ssize_t)dna.size());
  else {
    for (size_t j = 0; j < dna.size(); j += chunk) {
      size_t fin = j + chunk + (size_t)k - 1;
      if (fin > dna.size())
        fin = dna.size();
      collect((ssize_t)j, (ssize_t)fin);
    }
  }
  return uniq.size();
}

static void test_mod_full_scan_matches_chunked() {
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  const int k = 35, l = 31;
  std::mt19937_64 rng(42);
  std::string dna(200003, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];

  const size_t full = count_mod_unique(dna, k, l, toggle);
  const size_t chunked = count_mod_unique(dna, k, l, toggle, 16384);
  if (full != chunked) {
    std::cerr << "FAIL: Mod full-scan unique=" << full
              << " chunked=" << chunked << "\n";
    failures++;
  }
}

// Position-aware ResetAt: each chunk's first returns must match full-scan
// windows at the same genomic starts (build_db subblock overlap = k-1).
static void test_mod_chunked_windows_match_full() {
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  const int k = 35, l = 31;
  const size_t chunk = 4096;
  std::mt19937_64 rng(99);
  std::string dna(50000, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];

  const auto full = scan_mod_minimizers(dna, k, l, toggle);
  check(full.size() == dna.size() - (size_t)k + 1, "full stream L-k+1");

  MinimizerScanner sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  for (size_t i = 0; i < dna.size(); i += chunk) {
    size_t fin = i + chunk + (size_t)k - 1;
    if (fin > dna.size())
      fin = dna.size();
    if (fin <= i + (size_t)k - 1 && i > 0)
      break;
    sc.LoadSequence(dna, i, fin);
    size_t local = 0;
    while (uint64_t *mm = sc.NextMinimizer()) {
      if (sc.is_ambiguous())
        continue;
      const size_t win = i + local;
      if (win >= full.size()) {
        check(false, "chunk produced past end of full stream");
        break;
      }
      if (*mm != full[win]) {
        std::cerr << "FAIL: chunked window mismatch i=" << i
                  << " local=" << local << " win=" << win << "\n";
        failures++;
        return;
      }
      local++;
      // Owned prefix is `chunk` windows; overlap beyond that is optional.
      if (local >= chunk)
        break;
    }
    check(local > 0 || i + (size_t)k > dna.size(),
          "chunk should yield windows when long enough");
  }
}

// Non-zero origin must not change samples vs scanning the suffix in isolation
// only when origin==0 on that suffix — here we check offset origins against full.
static void test_mod_reset_at_nonzero_origin() {
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  const int k = 35, l = 31;
  std::string dna(500, 'A');
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = "ACGT"[i % 4];

  const auto full = scan_mod_minimizers(dna, k, l, toggle);
  const size_t origin = 17;  // not divisible by w=5
  MinimizerScanner sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  sc.LoadSequence(dna, origin, dna.size());
  size_t n = 0;
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (sc.is_ambiguous())
      continue;
    check(origin + n < full.size(), "origin stream within full");
    if (origin + n < full.size() && *mm != full[origin + n]) {
      std::cerr << "FAIL: ResetAt origin=" << origin << " n=" << n << "\n";
      failures++;
      return;
    }
    n++;
  }
  check(n == full.size() - origin, "suffix window count");
}

// Chunk starting after an N must use segment-relative origin, not raw offset.
static void test_mod_chunk_after_n_matches_full() {
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  const int k = 35, l = 31;
  std::string dna(400, 'A');
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = "ACGT"[i % 4];
  dna[80] = 'N';

  // Full non-ambig stream after the N (warm-up batched): windows of the
  // suffix starting at genomic 81.
  MinimizerScanner full_sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                           SketchingScheme::MOD, Order::RANDOM_XOR, "",
                           nullptr);
  full_sc.LoadSequence(dna);
  std::vector<uint64_t> after_n;
  bool past_n = false;
  while (uint64_t *mm = full_sc.NextMinimizer()) {
    if (full_sc.is_ambiguous()) {
      past_n = true;
      after_n.clear();  // restart collection after ambig return
      continue;
    }
    if (past_n)
      after_n.push_back(*mm);
  }
  check(!after_n.empty(), "expected windows after N");

  const size_t chunk_start = 120;  // well after N + k
  MinimizerScanner sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  sc.LoadSequence(dna, chunk_start, dna.size());
  // Stream origin should be chunk_start - 81 = 39 non-ambig bases after N.
  const size_t expect_origin = chunk_start - 81;
  size_t n = 0;
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (sc.is_ambiguous())
      continue;
    const size_t win = expect_origin + n;
    if (win >= after_n.size() || *mm != after_n[win]) {
      std::cerr << "FAIL: post-N chunk win=" << win << " n=" << n << "\n";
      failures++;
      return;
    }
    n++;
  }
  check(n == after_n.size() - expect_origin, "post-N suffix length");
}

static void test_deterministic_checksum() {
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  const int k = 35, l = 31;
  std::mt19937_64 rng(1);
  std::string dna(1000035, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];

  uint64_t checksum = 0;
  const auto vals = scan_mod_minimizers(dna, k, l, toggle);
  // Warm-up batching: one return per window after first k-mer.
  check(vals.size() == dna.size() - (size_t)k + 1,
        "checksum stream length is L-k+1");
  for (uint64_t v : vals)
    checksum ^= v;

  // Regression anchor: if sampled ℓ-mer stream changes, update this.
  // Warm-up batching does not change window samples (only drops pre-window
  // returns of ~0, which XOR out of the old L-length checksum).
  const uint64_t expect = 0x11c7d2e0bfbefd85ULL;
  if (checksum != expect) {
    std::cerr << "FAIL: checksum regression got=0x" << std::hex << checksum
              << " expect=0x" << expect << std::dec << "\n";
    failures++;
  } else {
    std::cout << "checksum regression 0x" << std::hex << checksum << std::dec
              << " OK\n";
  }
}

static void test_split_push_finish_matches_fused(const std::string &dna, int k,
                                                int l, uint64_t toggle) {
  ModSampler fused(k, l, 2, toggle, true, Order::RANDOM_XOR,
                   CURRENT_REVCOM_VERSION);
  ModSampler split(k, l, 2, toggle, true, Order::RANDOM_XOR,
                   CURRENT_REVCOM_VERSION);

  ssize_t seq_pos = 0;
  uint64_t last_f = UINT64_MAX, last_s = UINT64_MAX;
  for (char ch : dna) {
    const uint8_t code = dna_lut(ch);
    if (code > 3) {
      fused.OnAmbiguous();
      split.OnAmbiguous();
      seq_pos = 0;
      continue;
    }
    seq_pos++;
    const ssize_t win = seq_pos >= k ? seq_pos - k : (ssize_t)-1;
    uint64_t f = last_f, s = last_s;
    if (fused.PushAndFinish(code, win, &f))
      last_f = f;
    split.PushBase(code, nullptr);
    if (win >= 0 && split.FinishWindow(win, &s))
      last_s = s;
    check(last_f == last_s, "PushAndFinish vs PushBase+FinishWindow");
  }
}

static void test_mod_recovers_after_n() {
  // One N must not leave is_ambiguous sticky for the rest of the read.
  const int k = 35, l = 31;
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  std::string dna(120, 'A');
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = "ACGT"[i % 4];
  dna[40] = 'N';

  MinimizerScanner sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  sc.LoadSequence(dna);

  int ambig_returns = 0;
  int non_ambig_after = 0;
  bool saw_ambig = false;
  while (uint64_t *mm = sc.NextMinimizer()) {
    (void)mm;
    if (sc.is_ambiguous()) {
      ambig_returns++;
      saw_ambig = true;
      non_ambig_after = 0;
    } else if (saw_ambig) {
      non_ambig_after++;
    }
  }
  check(ambig_returns >= 1, "expected at least one ambiguous return at N");
  check(non_ambig_after > l, "Mod should resume non-ambiguous after N");
}

static void test_mod_warmup_batches_like_minimizer() {
  // Clean DNA: Mod should return once per window after first k-mer
  // (L-k+1), matching classic Minimizer cadence — not once per base.
  const int k = 35, l = 31;
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  std::string dna(150, 'A');
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = "ACGT"[i % 4];

  MinimizerScanner mod_sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                          SketchingScheme::MOD, Order::RANDOM_XOR, "", nullptr);
  mod_sc.LoadSequence(dna);
  size_t mod_calls = 0;
  size_t mod_ambig = 0;
  while (uint64_t *mm = mod_sc.NextMinimizer()) {
    (void)mm;
    mod_calls++;
    if (mod_sc.is_ambiguous())
      mod_ambig++;
  }

  MinimizerScanner min_sc(k, l, 0, true, toggle, CURRENT_REVCOM_VERSION,
                          SketchingScheme::MINIMIZER, Order::RANDOM_XOR, "",
                          nullptr);
  min_sc.LoadSequence(dna);
  size_t min_calls = 0;
  while (min_sc.NextMinimizer())
    min_calls++;

  check(mod_ambig == 0, "clean DNA should have no ambig returns");
  check(mod_calls == dna.size() - (size_t)k + 1,
        "Mod returns L-k+1 on clean DNA");
  check(mod_calls == min_calls,
        "Mod return count matches classic Minimizer on clean DNA");
}

static uint64_t pack_dna_lmer(const std::string &s, size_t pos, int l) {
  uint64_t x = 0;
  for (int i = 0; i < l; i++) {
    const char c = s[pos + (size_t)i];
    const uint64_t b = (c == 'C') ? 1ULL : (c == 'G') ? 2ULL
                      : (c == 'T') ? 3ULL : 0ULL;
    x = (x << 2) | b;
  }
  return x;
}

// GM ranks the raw ℓ-mer (gm_k orientation); the returned hash key is DNA-lex.
static void test_gm_scanner_rank_raw_key_canonical() {
  std::vector<uint64_t> raw(16, 16);
  // mirror-symmetric (canonical): ranks shared by a key and its mirror
  raw[0] = 0;
  raw[1] = raw[8] = 1;
  raw[2] = raw[4] = 2;
  raw[5] = raw[10] = 3;
  auto table = std::make_shared<const GreedyMiniTable>(raw);
  const int k = 10, l = 6;
  const uint64_t lmer_mask = (1ULL << (2 * l)) - 1ULL;
  MinimizerScanner sc(k, l, 0, true, 0, CURRENT_REVCOM_VERSION,
                      SketchingScheme::MINIMIZER, Order::GM, "", table);
  const std::string seq = "ACGTAAAACGTAAAACGTAA";
  sc.LoadSequence(seq);
  std::vector<uint64_t> got;
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (!sc.is_ambiguous())
      got.push_back(*mm);
  }

  GMRanker ranker(l, 2, 0, lmer_mask, true, table);
  const int w = k - l + 1;
  std::vector<uint64_t> expect;
  for (size_t start = 0; start + (size_t)k <= seq.size(); start++) {
    const uint64_t mer0 = pack_dna_lmer(seq, start, l);
    uint64_t best_rank = ranker.Rank(mer0);
    uint64_t best_mer = mer0;
    for (int j = 1; j < w; j++) {
      const uint64_t mer = pack_dna_lmer(seq, start + (size_t)j, l);
      const uint64_t r = ranker.Rank(mer);
      if (r < best_rank) {
        best_rank = r;
        best_mer = mer;
      }
    }
    const uint64_t key =
        dna_canonical_representation(best_mer, (uint8_t)l,
                                     CURRENT_REVCOM_VERSION);
    check(key <= dna_reverse_complement(key, (uint8_t)l, CURRENT_REVCOM_VERSION),
          "expected GM key is DNA-lex");
    expect.push_back(key);
  }
  check(got.size() == expect.size(), "GM scanner window count");
  check(got == expect, "GM scanner: Rank(raw) winner, DNA-lex hash key");
}

static std::shared_ptr<const GreedyMiniTable> synthetic_gm_k4() {
  std::vector<uint64_t> raw(16, 16);
  // mirror-symmetric (canonical): ranks shared by a key and its mirror
  raw[0] = 0;
  raw[1] = raw[8] = 1;
  raw[2] = raw[4] = 2;
  raw[5] = raw[10] = 3;
  return std::make_shared<const GreedyMiniTable>(raw);
}

static std::vector<uint64_t> scan_scheme_keys(
    const std::string &seq, int k, int l, uint64_t spaced, int revcom_version,
    SketchingScheme sketching, Order order,
    std::shared_ptr<const GreedyMiniTable> table) {
  MinimizerScanner sc(k, l, spaced, true, 0, revcom_version, sketching, order,
                      "", table);
  sc.LoadSequence(seq);
  std::vector<uint64_t> got;
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (!sc.is_ambiguous())
      got.push_back(*mm);
  }
  return got;
}

static std::vector<uint64_t> expect_gm_minimizer_keys(
    const std::string &seq, int k, int l, uint64_t spaced, int revcom_version,
    std::shared_ptr<const GreedyMiniTable> table) {
  const uint64_t lmer_mask =
      (l >= 32) ? ~0ULL : ((1ULL << (2 * l)) - 1ULL);
  GMRanker ranker(l, 2, 0, lmer_mask, true, table, spaced, revcom_version);
  const int w = k - l + 1;
  std::vector<uint64_t> expect;
  for (size_t start = 0; start + (size_t)k <= seq.size(); start++) {
    const uint64_t mer0 = pack_dna_lmer(seq, start, l);
    uint64_t best_rank = ranker.Rank(mer0);
    uint64_t best_mer = mer0;
    for (int j = 1; j < w; j++) {
      const uint64_t mer = pack_dna_lmer(seq, start + (size_t)j, l);
      const uint64_t r = ranker.Rank(mer);
      if (r < best_rank) {
        best_rank = r;
        best_mer = mer;
      }
    }
    uint64_t key =
        dna_canonical_representation(best_mer, (uint8_t)l, revcom_version);
    if (spaced)
      key &= spaced;
    expect.push_back(key);
  }
  return expect;
}

// MINIMIZER+GM: Rank ANDs each strand; returned key is DNA-lex(winner) AND mask.
static void test_gm_scanner_spaced_seed() {
  auto table = synthetic_gm_k4();
  const int k = 10, l = 6;
  const uint64_t lmer_mask = (1ULL << (2 * l)) - 1ULL;
  const uint64_t mask = lmer_mask & ~0xFULL;
  const std::string seq = "ACGTAAAACGTAAAACGTAAACGTAA";
  const auto got = scan_scheme_keys(seq, k, l, mask, CURRENT_REVCOM_VERSION,
                                    SketchingScheme::MINIMIZER, Order::GM,
                                    table);
  const auto expect = expect_gm_minimizer_keys(
      seq, k, l, mask, CURRENT_REVCOM_VERSION, table);
  check(got.size() == expect.size(), "GM spaced scanner window count");
  check(got == expect, "GM spaced scanner matches per-strand Rank + DNA-lex key");
  check(!got.empty(), "GM spaced produced keys");
  for (uint64_t key : got)
    check((key & mask) == key, "GM spaced hash key respects mask");

  const auto unmasked = scan_scheme_keys(
      seq, k, l, 0, CURRENT_REVCOM_VERSION, SketchingScheme::MINIMIZER,
      Order::GM, table);
  check(unmasked.size() == got.size(), "GM masked vs plain size");
  bool differ = false;
  for (size_t i = 0; i < got.size() && i < unmasked.size(); i++) {
    if (got[i] != unmasked[i]) {
      differ = true;
      break;
    }
  }
  check(differ, "GM spaced mask changes some hash keys");
}

// Rank-first scanning (compare table ranks, full packed values only on equal
// finite ranks; Sigma plane pext'd or rolled) must select exactly the
// minimizers of the Rank() path. Random mirror-symmetric tables with many
// equal and infinite ranks, centre extension l = gm_k + 2i, with and without
// (also non-base-aligned) spaced seeds, sequences with N and low complexity.
static std::shared_ptr<const GreedyMiniTable> random_symmetric_table(
    int gm_k, std::mt19937_64 &rng) {
  const size_t n = (size_t)1 << gm_k;
  std::vector<uint64_t> raw(n, n);
  uint64_t next_rank = 0;
  for (size_t key = 0; key < n; key++) {
    const size_t mir = (size_t)gm_centre::mirror_key((uint32_t)key, gm_k);
    if (mir < key || raw[key] != n || (rng() & 3ULL) == 0)
      continue;
    raw[key] = raw[mir] = next_rank++;
  }
  return std::make_shared<const GreedyMiniTable>(raw);
}

static std::vector<uint64_t> scan_all(MinimizerScanner &sc,
                                      const std::vector<std::string> &seqs) {
  std::vector<uint64_t> out;
  for (const std::string &seq : seqs) {
    sc.LoadSequence(seq);
    while (uint64_t *mm = sc.NextMinimizer())
      out.push_back(sc.is_ambiguous() ? ~0ULL : *mm);
    sc.LoadSequence(seq, 3, seq.size() > 20 ? seq.size() - 7 : seq.size());
    while (uint64_t *mm = sc.NextMinimizer())
      out.push_back(sc.is_ambiguous() ? ~0ULL : *mm);
  }
  return out;
}

static void test_gm_rank_first_matches_rank_path() {
  std::mt19937_64 rng(99);
  std::vector<std::string> seqs;
  const char *bases = "ACGT";
  for (int i = 0; i < 30; i++) {
    std::string s;
    const size_t n = (i < 4) ? (size_t)(rng() % 40) : 100 + rng() % 600;
    for (size_t j = 0; j < n; j++) {
      uint64_t r = rng();
      char c = bases[r & 3ULL];
      if (i % 3 == 1 && (r >> 8) % 40 == 0)
        c = 'N';
      if (i % 3 == 2)
        c = bases[(r & 3ULL) & ((r >> 5) % 3 ? 1ULL : 3ULL)];
      s.push_back(c);
    }
    seqs.push_back(s);
  }
  seqs.push_back(std::string(200, 'A'));
  seqs.push_back("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

  for (int gm_k : {4, 6, 9}) {
    for (int trial = 0; trial < 3; trial++) {
      auto table = random_symmetric_table(gm_k, rng);
      for (int l = gm_k; l <= 31; l += 2) {
        const int k = l + 1 + (int)(rng() % 6);
        if (k > 31 + 8)
          continue;
        const uint64_t lmer_mask = (l >= 32) ? ~0ULL : ((1ULL << (2 * l)) - 1);
        for (int mode = 0; mode < 3; mode++) {
          uint64_t mask = 0;
          if (mode == 1)
            mask = MakeSpacedSeedMask(l, l / 4, 2);
          else if (mode == 2)
            mask = rng() & lmer_mask;  // not aligned to bases
          const uint64_t toggle = (trial == 0) ? 0 : rng();
          for (int revcom : {1, 0}) {
            MinimizerScanner sc(k, l, mask, true, toggle, revcom,
                                SketchingScheme::MINIMIZER, Order::GM, "",
                                table);
            sc.SetRankFirst(false);
            const auto ref = scan_all(sc, seqs);
            sc.SetRankFirst(true, false);
            const bool same_pext = (scan_all(sc, seqs) == ref);
            sc.SetRankFirst(true, true);
            const bool same_roll = (scan_all(sc, seqs) == ref);
            check(same_pext && same_roll && !ref.empty(),
                  "GM rank-first scanner == Rank() scanner");
            if (!(same_pext && same_roll))
              std::cerr << "  gm_k=" << gm_k << " l=" << l << " k=" << k
                        << " mode=" << mode << " revcom=" << revcom << "\n";
          }
        }
      }
    }
  }
}

static void test_gm_scanner_revcom_v0() {
  auto table = synthetic_gm_k4();
  const int k = 10, l = 6;
  std::mt19937_64 rng(11);
  std::string seq(80, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < seq.size(); i++)
    seq[i] = bases[rng() & 3ULL];

  const auto got0 = scan_scheme_keys(seq, k, l, 0, /*revcom*/0,
                                     SketchingScheme::MINIMIZER, Order::GM,
                                     table);
  const auto exp0 =
      expect_gm_minimizer_keys(seq, k, l, 0, /*revcom*/0, table);
  check(got0 == exp0, "GM scanner revcom_version=0 matches Rank+key");

  const auto got1 = scan_scheme_keys(seq, k, l, 0, CURRENT_REVCOM_VERSION,
                                     SketchingScheme::MINIMIZER, Order::GM,
                                     table);
  check(got0.size() == got1.size(), "GM v0 vs v1 window count");
  bool differ = false;
  for (size_t i = 0; i < got0.size() && i < got1.size(); i++) {
    if (got0[i] != got1[i]) {
      differ = true;
      break;
    }
  }
  check(differ, "GM scanner revcom_version=0 differs from v1");
}

// MOD+GM: GMRanker on raw t-mers (ACGT→ATCG σ-extend); ℓ-mer store is masked.
static void test_mod_gm_matches_naive() {
  auto table = synthetic_gm_k4();  // gm_k=4 == t for (k,l)=(6,4)
  const int k = 6, l = 4;
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  std::mt19937_64 rng(21);
  std::string dna(400, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];

  test_naive_matches_scanner(k, l, dna, toggle, "MOD+GM naive", 0,
                             CURRENT_REVCOM_VERSION, Order::GM, table);

  const uint64_t lmer_mask = (1ULL << (2 * l)) - 1ULL;
  const uint64_t mask = lmer_mask & ~0x3ULL;
  test_naive_matches_scanner(k, l, dna, toggle, "MOD+GM spaced naive", mask,
                             CURRENT_REVCOM_VERSION, Order::GM, table);

  const auto masked = scan_mod_minimizers(dna, k, l, toggle, mask,
                                          CURRENT_REVCOM_VERSION, Order::GM,
                                          table);
  check(!masked.empty(), "MOD+GM spaced produced samples");
  MinimizerScanner sc_mask(k, l, mask, true, toggle, CURRENT_REVCOM_VERSION,
                           SketchingScheme::MOD, Order::GM, "", table);
  sc_mask.LoadSequence(dna);
  size_t mask_checked = 0;
  size_t valid_keys = 0;
  while (uint64_t *mm = sc_mask.NextMinimizer()) {
    if (sc_mask.is_ambiguous())
      continue;
    if (*mm == UINT64_MAX)
      continue;
    valid_keys++;
    check((*mm & mask) == *mm, "MOD+GM returned key respects spaced mask");
    mask_checked++;
  }
  check(valid_keys > 50, "MOD+GM spaced produced valid keys");
  check(mask_checked == valid_keys, "MOD+GM spaced checked all valid keys");

  // By design: ranking masked t-mers would disagree with ℓ-mer-only masking.
  const auto designed =
      naive_mod_stream(dna, k, l, toggle, mask, CURRENT_REVCOM_VERSION, table,
                       /*mask_tmer_rank=*/false);
  const auto wrong =
      naive_mod_stream(dna, k, l, toggle, mask, CURRENT_REVCOM_VERSION, table,
                       /*mask_tmer_rank=*/true);
  check(designed.size() == wrong.size(), "MOD spaced t-mer-rank size");
  bool differ = false;
  for (size_t i = 0; i < designed.size() && i < wrong.size(); i++) {
    if (designed[i] != wrong[i]) {
      differ = true;
      break;
    }
  }
  check(differ, "MOD t-mer rank is unmasked (differs from AND-then-rank)");
  check(masked == designed, "MOD+GM scanner matches unmasked t-mer rank");
}

static void test_mod_revcom_v0() {
  const int k = 35, l = 31;
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  std::mt19937_64 rng(13);
  std::string dna(800, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];

  test_naive_matches_scanner(k, l, dna, toggle, "MOD revcom v0", 0, 0);

  const auto v0 = scan_mod_minimizers(dna, k, l, toggle, 0, /*revcom*/0);
  const auto v1 = scan_mod_minimizers(dna, k, l, toggle, 0,
                                      CURRENT_REVCOM_VERSION);
  check(v0.size() == v1.size(), "MOD v0 vs v1 size");
  bool differ = false;
  for (size_t i = 0; i < v0.size() && i < v1.size(); i++) {
    if (v0[i] != v1[i]) {
      differ = true;
      break;
    }
  }
  check(differ, "MOD revcom_version=0 differs from v1");
}

// The per-base push is specialised on the order / LUT choice; the debug-log
// path uses the generic runtime-branching code. Both must give the same
// stream, and the logs must hold the canonical (GM: raw) t-mers and the
// canonical, spaced-seed-masked l-mers.
static void test_mod_kinds_and_logs() {
  std::mt19937_64 rng(77);
  std::string dna(3000, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna.size(); i++)
    dna[i] = bases[rng() & 3ULL];
  dna[1000] = 'N';
  struct Cfg { int k, l; Order order; bool gm; uint64_t spaced; int revcom; };
  const uint64_t smask = 0x3ffffffff3333333ULL;
  const Cfg cfgs[] = {
      {35, 31, Order::RANDOM_XOR, false, 0, 1},
      {35, 31, Order::RANDOM_XOR, false, smask, 1},
      {35, 31, Order::RANDOM_XOR, false, smask, 0},
      {35, 31, Order::RANDOM_PCG, false, smask, 1},   // t=6: canonical LUT
      {40, 31, Order::RANDOM_PCG, false, smask, 1},   // t=11: ranker, no LUT
      {40, 31, Order::RANDOM_XOR, false, 0, 0},       // t=11 xor, revcom v0
      {6, 4, Order::GM, true, 0, 1},                  // t=4: raw LUT
      {6, 4, Order::GM, true, 0x3fULL << 2, 1},
  };
  for (const Cfg &c : cfgs) {
    std::shared_ptr<const GreedyMiniTable> table;
    if (c.gm)
      table = synthetic_gm_k4();
    const uint64_t toggle = DEFAULT_TOGGLE_MASK;
    ModSampler plain(c.k, c.l, 2, toggle, true, c.order, c.revcom, c.spaced,
                     table);
    ModSampler logged(c.k, c.l, 2, toggle, true, c.order, c.revcom, c.spaced,
                      table);
    std::vector<uint64_t> tlog, llog;
    logged.SetRankedTmerLog(&tlog);
    logged.SetRankedLmerLog(&llog);
    const int t = plain.t();
    const uint64_t lmask = (1ULL << (2 * c.l)) - 1ULL;
    std::vector<uint64_t> a, b, exp_t, exp_l;
    ssize_t seq_pos = 0;
    uint64_t fwd_t = 0, fwd_l = 0;
    int loaded = 0;
    uint64_t la = UINT64_MAX, lb = UINT64_MAX;
    for (char ch : dna) {
      const uint8_t code = dna_lut(ch);
      if (code > 3) {
        plain.OnAmbiguous();
        logged.OnAmbiguous();
        seq_pos = 0;
        loaded = 0;
        fwd_t = fwd_l = 0;
        continue;
      }
      seq_pos++;
      loaded++;
      fwd_t = ((fwd_t << 2) | code) & ((1ULL << (2 * t)) - 1ULL);
      fwd_l = ((fwd_l << 2) | code) & lmask;
      if (loaded >= t)
        exp_t.push_back(c.gm ? fwd_t
                             : dna_canonical_representation(fwd_t, (uint8_t)t,
                                                            c.revcom));
      if (loaded >= c.l) {
        uint64_t cl =
            dna_canonical_representation(fwd_l, (uint8_t)c.l, c.revcom);
        exp_l.push_back(c.spaced ? (cl & c.spaced & lmask) : cl);
      }
      const ssize_t win = seq_pos >= c.k ? seq_pos - c.k : (ssize_t)-1;
      if (plain.PushAndFinish(code, win, &la))
        a.push_back(la);
      if (logged.PushAndFinish(code, win, &lb))
        b.push_back(lb);
    }
    check(!a.empty() && a == b, "kinds: specialised stream == logged stream");
    check(tlog == exp_t, "kinds: t-mer log");
    check(llog == exp_l, "kinds: l-mer log");
  }
}

int main() {
  const uint64_t toggle = DEFAULT_TOGGLE_MASK;
  std::mt19937_64 rng(7);

  for (int l : {15, 31}) {
    test_incremental_rc_equals_canonical(l);
  }
  test_lut_ring_index(5);
  test_lut_ring_index(7);

  for (int trial = 0; trial < 8; trial++) {
    const int l = 31;
    const int k = (trial % 2 == 0) ? 35 : 37;  // w=5 and w=7
    std::string dna((size_t)(2000 + trial * 131), 'A');
    const char *bases = "ACGT";
    for (size_t i = 0; i < dna.size(); i++)
      dna[i] = bases[rng() & 3ULL];
    test_naive_matches_scanner(k, l, dna, toggle,
                               (std::string("naive k=") + std::to_string(k))
                                   .c_str());
  }

  {
    std::string dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    test_mod_sampler_direct(35, 31, dna, toggle);
    test_split_push_finish_matches_fused(dna, 35, 31, toggle);
  }

  test_mod_full_scan_matches_chunked();
  test_mod_chunked_windows_match_full();
  test_mod_reset_at_nonzero_origin();
  test_mod_chunk_after_n_matches_full();
  test_deterministic_checksum();
  test_mod_recovers_after_n();
  test_mod_warmup_batches_like_minimizer();
  test_mod_spaced_seed_masking();
  test_gm_scanner_rank_raw_key_canonical();
  test_gm_scanner_spaced_seed();
  test_gm_scanner_revcom_v0();
  test_gm_rank_first_matches_rank_path();
  test_mod_gm_matches_naive();
  test_mod_revcom_v0();
  test_mod_kinds_and_logs();

  if (failures) {
    std::cerr << failures << " failure(s)\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
