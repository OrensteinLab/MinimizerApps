/*
 * Microbenchmark: RandomXor / RandomPcg / GreedyMini via
 * MinimizerScanner::NextMinimizer.
 * Reports mean throughput, batch latency distribution, and (GreedyMini only)
 * Rank() UHS stats and isolated Rank() replay by UHS vs non-UHS path.
 */

#include "mmscanner.h"
#include "utilities.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <err.h>
#include <fstream>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <random>
#include <string>
#include <sysexits.h>
#include <vector>
#include <x86intrin.h>

using namespace kraken2;
using Clock = std::chrono::steady_clock;

struct DistStats {
  double mean_ns;
  double std_ns;
  double min_ns;
  double p50_ns;
  double p90_ns;
  double p99_ns;
  double p999_ns;
  double max_ns;
  size_t n_samples;
};

struct BenchResult {
  double ns_per_op;
  double ops_per_sec;
  uint64_t n_ops;
  uint64_t checksum;
  DistStats dist;
  bool has_uhs_stats = false;
  RankHitStats uhs_stats;
  NextMinPhaseStats phases;
};

struct PathDurationStats {
  uint64_t n_calls = 0;
  double mean_ns = 0;       // from total cycles (all calls)
  DistStats sample_dist;    // from subsampled per-call timings
};

static double percentile_sorted(const std::vector<double> &sorted, double p) {
  if (sorted.empty())
    return 0;
  const double idx = p * (sorted.size() - 1);
  const size_t lo = (size_t)idx;
  const size_t hi = std::min(lo + 1, sorted.size() - 1);
  const double frac = idx - lo;
  return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

static DistStats summarize(std::vector<double> samples) {
  DistStats d = {};
  d.n_samples = samples.size();
  if (samples.empty())
    return d;
  std::sort(samples.begin(), samples.end());
  d.min_ns = samples.front();
  d.max_ns = samples.back();
  d.p50_ns = percentile_sorted(samples, 0.50);
  d.p90_ns = percentile_sorted(samples, 0.90);
  d.p99_ns = percentile_sorted(samples, 0.99);
  d.p999_ns = percentile_sorted(samples, 0.999);
  double sum = 0;
  for (double x : samples)
    sum += x;
  d.mean_ns = sum / samples.size();
  double var = 0;
  for (double x : samples) {
    const double dlt = x - d.mean_ns;
    var += dlt * dlt;
  }
  d.std_ns = std::sqrt(var / samples.size());
  return d;
}

static inline uint64_t rdtsc_start() {
  unsigned aux;
  _mm_mfence();
  return __rdtscp(&aux);
}

static inline uint64_t rdtsc_end() {
  unsigned aux;
  const uint64_t t = __rdtscp(&aux);
  _mm_lfence();
  return t;
}

// Calibrate TSC → ns using steady_clock over a short spin.
static double ns_per_cycle() {
  const auto wall0 = Clock::now();
  const uint64_t c0 = rdtsc_start();
  // Busy-wait ~20ms so calibration is stable.
  while (std::chrono::duration<double, std::milli>(Clock::now() - wall0).count() <
         20.0) {
  }
  const uint64_t c1 = rdtsc_end();
  const auto wall1 = Clock::now();
  const double ns =
      std::chrono::duration<double, std::nano>(wall1 - wall0).count();
  const uint64_t cycles = c1 - c0;
  return cycles ? (ns / (double)cycles) : 0.0;
}

// Time NextMinimizer in batches so clock overhead does not dominate (~ns calls).
static BenchResult scan_once(MinimizerScanner &sc, const std::string &dna,
                             size_t batch) {
  sc.LoadSequence(dna);
  while (sc.NextMinimizer()) {
  }

  sc.LoadSequence(dna);
  sc.ClearRankStats();
  sc.ClearPhaseStats();
  std::vector<double> batch_ns;
  batch_ns.reserve(dna.size() / batch + 8);
  uint64_t checksum = 0;
  uint64_t n_ops = 0;
  uint64_t in_batch = 0;

  auto batch_t0 = Clock::now();
  auto t0_all = batch_t0;
  while (uint64_t *mm = sc.NextMinimizer()) {
    checksum ^= *mm;
    n_ops++;
    in_batch++;
    if (in_batch == batch) {
      auto t1 = Clock::now();
      double ns = std::chrono::duration<double, std::nano>(t1 - batch_t0).count();
      batch_ns.push_back(ns / batch);
      batch_t0 = t1;
      in_batch = 0;
    }
  }
  if (in_batch > 0) {
    auto t1 = Clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1 - batch_t0).count();
    batch_ns.push_back(ns / in_batch);
  }
  auto t1_all = Clock::now();
  double sec = std::chrono::duration<double>(t1_all - t0_all).count();

  BenchResult r;
  r.n_ops = n_ops ? n_ops : 1;
  r.ns_per_op = (sec * 1e9) / r.n_ops;
  r.ops_per_sec = r.n_ops / sec;
  r.checksum = checksum;
  r.dist = summarize(std::move(batch_ns));
  r.has_uhs_stats = sc.has_uhs_stats();
  r.uhs_stats = sc.GetRankStats();
  r.phases = sc.GetPhaseStats();
  return r;
}

static void print_dist(const char *name, const DistStats &d, const char *unit) {
  std::printf(
      "%-14s  mean=%6.2f  std=%6.2f  min=%6.2f  p50=%6.2f  p90=%6.2f  "
      "p99=%6.2f  p999=%6.2f  max=%6.2f  (%s, %zu samples)\n",
      name, d.mean_ns, d.std_ns, d.min_ns, d.p50_ns, d.p90_ns, d.p99_ns,
      d.p999_ns, d.max_ns, unit, d.n_samples);
}

static void print_phase_breakdown(const char *scheme, const NextMinPhaseStats &p,
                                  double nspc, double wall_ns_per_op) {
  if (!p.steps)
    return;
  const double enc = (p.encode_cycles * nspc) / p.steps;
  const double can = (p.canonical_cycles * nspc) / p.steps;
  const double rank = (p.rank_cycles * nspc) / p.steps;
  const double deq = (p.deque_cycles * nspc) / p.steps;
  const double sum = enc + can + rank + deq;
  const double tot = p.total_cycles() ? (double)p.total_cycles() : 1.0;
  std::printf(
      "\n%s NextMinimizer phase split (RDTSC in-scanner; %llu steps):\n"
      "  Phase        ns/step    %% of phases   note\n"
      "  ----------  ---------  ------------  ----\n"
      "  encode      %9.3f  %11.1f%%  base→ℓ-mer\n"
      "  canonical   %9.3f  %11.1f%%  DNA-lex ℓ-mer (skipped for GM)\n"
      "  rank        %9.3f  %11.1f%%  Rank()\n"
      "  deque       %9.3f  %11.1f%%  window min\n"
      "  sum(phases) %9.3f               (incl. rdtscp boundaries)\n"
      "  wall/op     %9.3f               (steady_clock batch mean)\n",
      scheme, (unsigned long long)p.steps, enc, 100.0 * p.encode_cycles / tot,
      can, 100.0 * p.canonical_cycles / tot, rank,
      100.0 * p.rank_cycles / tot, deq, 100.0 * p.deque_cycles / tot, sum,
      wall_ns_per_op);
}

// Stateful Rank() replay timed with RDTSC, bucketed by UHS vs non-UHS return.
static void time_rank_paths(const std::vector<uint64_t> &mers,
                            MinimizerRanker &ranker, size_t sample_stride,
                            size_t max_samples, double nspc,
                            PathDurationStats &uhs, PathDurationStats &non_uhs) {
  uhs = {};
  non_uhs = {};
  std::vector<double> uhs_samples;
  std::vector<double> non_samples;
  uhs_samples.reserve(std::min(max_samples, mers.size() / sample_stride + 8));
  non_samples.reserve(std::min(max_samples, mers.size() / sample_stride + 8));

  uint64_t uhs_cycles = 0, non_cycles = 0;

  ranker.Reset();
  ranker.ClearRankStats();
  for (size_t i = 0; i < mers.size(); i++) {
    const uint64_t t0 = rdtsc_start();
    const uint64_t rank = ranker.Rank(mers[i]);
    const uint64_t t1 = rdtsc_end();
    const uint64_t cyc = t1 - t0;
    const bool is_uhs = (rank != UINT64_MAX);
    if (is_uhs) {
      uhs.n_calls++;
      uhs_cycles += cyc;
      if ((uhs.n_calls % sample_stride) == 0 && uhs_samples.size() < max_samples)
        uhs_samples.push_back(cyc * nspc);
    } else {
      non_uhs.n_calls++;
      non_cycles += cyc;
      if ((non_uhs.n_calls % sample_stride) == 0 &&
          non_samples.size() < max_samples)
        non_samples.push_back(cyc * nspc);
    }
  }

  uhs.mean_ns = uhs.n_calls ? (uhs_cycles * nspc) / uhs.n_calls : 0;
  non_uhs.mean_ns =
      non_uhs.n_calls ? (non_cycles * nspc) / non_uhs.n_calls : 0;
  uhs.sample_dist = summarize(std::move(uhs_samples));
  non_uhs.sample_dist = summarize(std::move(non_samples));
}

// Density = fraction of successive windows whose selected minimizer value
// changes (TOK_LOOKUP / scan-step proxy). Also track super-mer run lengths
// (maximal runs of identical non-ambiguous keys).
struct DensityStats {
  double density = 0;
  double sm_mean = 0;
  uint64_t sm_max = 0;
  uint64_t n_steps = 0;
  uint64_t n_supermers = 0;
};

static void accumulate_scan_stats(MinimizerScanner &sc, const std::string &dna,
                                  uint64_t *distinct, uint64_t *total,
                                  uint64_t *sm_sum, uint64_t *sm_max,
                                  uint64_t *sm_count) {
  if (dna.size() < 2)
    return;
  sc.LoadSequence(dna);
  uint64_t last = UINT64_MAX;
  uint64_t run = 0;
  bool have = false;
  auto close_run = [&]() {
    if (!have)
      return;
    (*sm_sum) += run;
    (*sm_count)++;
    if (run > *sm_max)
      *sm_max = run;
    have = false;
    run = 0;
    last = UINT64_MAX;
  };
  while (uint64_t *mm = sc.NextMinimizer()) {
    if (sc.is_ambiguous()) {
      // Ambiguous / no-selection breaks the sequence; gap bound is per segment.
      close_run();
      continue;
    }
    (*total)++;
    if (!have || *mm != last) {
      close_run();
      (*distinct)++;
      last = *mm;
      run = 1;
      have = true;
    } else {
      run++;
    }
  }
  close_run();
}

static DensityStats density_stats_from_counts(uint64_t distinct, uint64_t total,
                                              uint64_t sm_sum, uint64_t sm_max,
                                              uint64_t sm_count) {
  DensityStats st;
  st.n_steps = total;
  st.n_supermers = sm_count;
  st.density = total ? (double)distinct / (double)total : 0.0;
  st.sm_mean = sm_count ? (double)sm_sum / (double)sm_count : 0.0;
  st.sm_max = sm_max;
  return st;
}

static DensityStats density_proxy(MinimizerScanner &sc,
                                  const std::vector<std::string> &seqs) {
  uint64_t distinct = 0, total = 0, sm_sum = 0, sm_max = 0, sm_count = 0;
  for (const auto &dna : seqs)
    accumulate_scan_stats(sc, dna, &distinct, &total, &sm_sum, &sm_max,
                          &sm_count);
  return density_stats_from_counts(distinct, total, sm_sum, sm_max, sm_count);
}

static DensityStats density_proxy(MinimizerScanner &sc, const std::string &dna) {
  std::vector<std::string> one{dna};
  return density_proxy(sc, one);
}

// Random minimizer: 2/(w+1). Mod: 1/w.
static double density_theory(SketchingScheme sketching, int k, int l) {
  const int w = k - l + 1;
  if (w < 1)
    return 0.0;
  if (sketching == SketchingScheme::MOD)
    return 1.0 / (double)w;
  return 2.0 / (double)(w + 1);
}

static bool ends_with_ci(const std::string &name, const char *suf) {
  const size_t n = std::strlen(suf);
  if (name.size() < n)
    return false;
  for (size_t i = 0; i < n; i++) {
    char a = name[name.size() - n + i];
    char b = suf[i];
    if (a >= 'A' && a <= 'Z')
      a = (char)(a - 'A' + 'a');
    if (b >= 'A' && b <= 'Z')
      b = (char)(b - 'A' + 'a');
    if (a != b)
      return false;
  }
  return true;
}

static bool is_library_fasta_name(const std::string &name) {
  return ends_with_ci(name, ".fa") || ends_with_ci(name, ".fna") ||
         ends_with_ci(name, ".fasta") || ends_with_ci(name, ".faa");
}

// Density factor normalizes by the random-minimizer scale: dens × (w+1).
// Random minimizer theory ⇒ factor 2; Mod 1/w ⇒ factor (w+1)/w.
static double density_factor(double dens, int w) {
  return dens * (double)(w + 1);
}

static std::string short_fasta_label(const std::string &path) {
  // Prefer library/<name>/library.fna → <name>
  const std::string marker = "/library/";
  const size_t pos = path.rfind(marker);
  if (pos != std::string::npos) {
    std::string rest = path.substr(pos + marker.size());
    const size_t slash = rest.find('/');
    if (slash != std::string::npos)
      return rest.substr(0, slash);
    return rest;
  }
  const size_t slash = path.find_last_of('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

// Thread-safe particular-density progress (byte-weighted across scheme×file).
struct DensityProgress {
  std::mutex mu;
  uint64_t total_bytes = 0;
  uint64_t done_bytes = 0;
  int total_tasks = 0;
  int done_tasks = 0;
  std::vector<std::string> active;  // scheme:file currently in flight
  std::chrono::steady_clock::time_point t0 =
      std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point last_print =
      std::chrono::steady_clock::now();

  static std::string label(const char *scheme, const std::string &path) {
    return std::string(scheme ? scheme : "?") + ":" + short_fasta_label(path);
  }

  void set_active(const char *scheme, const std::string &path) {
    std::lock_guard<std::mutex> lock(mu);
    active.push_back(label(scheme, path));
    print_unlocked(true);
  }

  void add_bytes(uint64_t n) {
    if (!n)
      return;
    std::lock_guard<std::mutex> lock(mu);
    done_bytes += n;
    if (done_bytes > total_bytes)
      done_bytes = total_bytes;
    print_unlocked(false);
  }

  void finish_task(const char *scheme, const std::string &path) {
    std::lock_guard<std::mutex> lock(mu);
    done_tasks++;
    const std::string lab = label(scheme, path);
    for (auto it = active.begin(); it != active.end(); ++it) {
      if (*it == lab) {
        active.erase(it);
        break;
      }
    }
    print_unlocked(true);
  }

  void finish_all() {
    std::lock_guard<std::mutex> lock(mu);
    done_bytes = total_bytes;
    done_tasks = total_tasks;
    active.clear();
    print_unlocked(true);
    std::cerr << "\n";
  }

  void print_unlocked(bool force) {
    const auto now = std::chrono::steady_clock::now();
    if (!force &&
        std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print)
                .count() < 250)
      return;
    last_print = now;
    const double pct =
        total_bytes ? (100.0 * (double)done_bytes / (double)total_bytes) : 0.0;
    const int filled = (int)(pct / 5.0 + 0.5);  // 20-char bar
    char bar[21];
    for (int i = 0; i < 20; i++)
      bar[i] = i < filled ? '#' : '-';
    bar[20] = '\0';
    const double elapsed =
        std::chrono::duration<double>(now - t0).count();
    double eta = 0.0;
    if (done_bytes > 0 && done_bytes < total_bytes && elapsed > 0.5)
      eta = elapsed * ((double)(total_bytes - done_bytes) / (double)done_bytes);
    std::cerr << "\r[particular] " << (int)(pct + 0.5) << "% |" << bar << "| "
              << done_tasks << "/" << total_tasks << " files  "
              << (done_bytes / 1000000) << "/" << (total_bytes / 1000000)
              << " Mbp  ";
    if (!active.empty()) {
      std::cerr << "active:";
      const size_t nshow = active.size() < 3 ? active.size() : 3;
      for (size_t i = 0; i < nshow; i++)
        std::cerr << (i ? "," : "") << active[i];
      if (active.size() > nshow)
        std::cerr << ",+" << (active.size() - nshow);
      std::cerr << "  ";
    }
    std::cerr << elapsed << "s";
    if (eta > 0.0)
      std::cerr << "  ETA " << (int)(eta + 0.5) << "s";
    std::cerr << "    " << std::flush;
  }
};

// Byte ranges [begin, end) of each FASTA sequence body (excludes '>' headers).
static std::vector<std::pair<uint64_t, uint64_t>>
index_fasta_seq_ranges(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    errx(EX_NOINPUT, "cannot open %s", path.c_str());
  std::vector<std::pair<uint64_t, uint64_t>> ranges;
  std::string line;
  uint64_t offset = 0;
  bool in_seq = false;
  while (std::getline(in, line)) {
    const uint64_t line_off = offset;
    std::streamoff after_off = in.tellg();
    uint64_t after =
        after_off >= 0 ? (uint64_t)after_off : offset + line.size() + 1;
    if (!line.empty() && line[0] == '>') {
      if (in_seq)
        ranges.back().second = line_off;
      in_seq = true;
      ranges.push_back({after, 0});
    }
    offset = after;
  }
  if (in_seq)
    ranges.back().second = offset;
  return ranges;
}

static void read_fasta_seq_range(const std::string &path, uint64_t begin,
                                 uint64_t end, std::string *dna) {
  dna->clear();
  if (end <= begin)
    return;
  std::ifstream in(path, std::ios::binary);
  if (!in)
    errx(EX_NOINPUT, "cannot open %s", path.c_str());
  in.seekg((std::streamoff)begin);
  dna->reserve((size_t)std::min(end - begin, (uint64_t)(64ull << 20)));
  std::string line;
  while (true) {
    const std::streamoff pos = in.tellg();
    if (pos < 0 || (uint64_t)pos >= end)
      break;
    if (!std::getline(in, line))
      break;
    if (!line.empty() && line[0] == '>')
      break;
    for (char c : line) {
      if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' ||
          c == 'a' || c == 'c' || c == 'g' || c == 't' || c == 'n')
        dna->push_back(c);
    }
  }
}

// Parallel over sequences inside one library FASTA (each record independent).
struct ScannerSpec {
  int k = 0;
  int l = 0;
  uint64_t toggle = 0;
  uint64_t spaced_seed_mask = 0;
  SketchingScheme sketching = SketchingScheme::MINIMIZER;
  Order order = Order::RANDOM_XOR;
  std::string gm_path;
  std::shared_ptr<const GreedyMiniTable> gm_tbl;
};

static MinimizerScanner make_scanner(const ScannerSpec &spec) {
  return MinimizerScanner(spec.k, spec.l, spec.spaced_seed_mask, true,
                          spec.toggle, CURRENT_REVCOM_VERSION, spec.sketching,
                          spec.order, spec.gm_path, spec.gm_tbl);
}

static void density_proxy_one_file_mt(const ScannerSpec &spec,
                                      const std::string &path, size_t file_size_est,
                                      int n_threads, uint64_t *distinct,
                                      uint64_t *total, uint64_t *sm_sum,
                                      uint64_t *sm_max, uint64_t *sm_count,
                                      DensityProgress *prog, const char *scheme) {
  if (prog)
    prog->set_active(scheme, path);

  const std::vector<std::pair<uint64_t, uint64_t>> ranges =
      index_fasta_seq_ranges(path);
  if (ranges.empty()) {
    if (prog) {
      if (file_size_est)
        prog->add_bytes(file_size_est);
      prog->finish_task(scheme, path);
    }
    return;
  }

  if (n_threads < 1)
    n_threads = 1;
  // Cap threads by sequence count (tiny files).
  if (n_threads > (int)ranges.size())
    n_threads = (int)ranges.size();

  if (n_threads == 1) {
    MinimizerScanner sc = make_scanner(spec);
    std::string dna;
    uint64_t span_sum = 0;
    for (size_t i = 0; i < ranges.size(); i++) {
      read_fasta_seq_range(path, ranges[i].first, ranges[i].second, &dna);
      accumulate_scan_stats(sc, dna, distinct, total, sm_sum, sm_max, sm_count);
      const uint64_t span = ranges[i].second > ranges[i].first
                                ? ranges[i].second - ranges[i].first
                                : 0;
      span_sum += span;
      if (prog && span)
        prog->add_bytes(span);
    }
    if (prog) {
      if (file_size_est > span_sum)
        prog->add_bytes(file_size_est - span_sum);
      prog->finish_task(scheme, path);
    }
    return;
  }

  uint64_t ld = 0, lt = 0, ls = 0, lm = 0, lc = 0, span_sum = 0;
#pragma omp parallel num_threads(n_threads) reduction(+ : ld, lt, ls, lc, span_sum)
  {
    MinimizerScanner sc = make_scanner(spec);
    std::string dna;
    uint64_t local_max = 0;
    uint64_t local_reported = 0;
#pragma omp for schedule(dynamic, 8)
    for (int i = 0; i < (int)ranges.size(); i++) {
      uint64_t d = 0, t = 0, s = 0, m = 0, c = 0;
      read_fasta_seq_range(path, ranges[(size_t)i].first,
                           ranges[(size_t)i].second, &dna);
      accumulate_scan_stats(sc, dna, &d, &t, &s, &m, &c);
      ld += d;
      lt += t;
      ls += s;
      lc += c;
      if (m > local_max)
        local_max = m;
      const uint64_t span =
          ranges[(size_t)i].second > ranges[(size_t)i].first
              ? ranges[(size_t)i].second - ranges[(size_t)i].first
              : 0;
      span_sum += span;
      local_reported += span;
      if (prog && local_reported >= (16ull << 20)) {
        prog->add_bytes(local_reported);
        local_reported = 0;
      }
    }
    if (prog && local_reported)
      prog->add_bytes(local_reported);
#pragma omp critical
    {
      if (local_max > lm)
        lm = local_max;
    }
  }
  *distinct += ld;
  *total += lt;
  *sm_sum += ls;
  *sm_count += lc;
  if (lm > *sm_max)
    *sm_max = lm;

  if (prog) {
    if (file_size_est > span_sum)
      prog->add_bytes(file_size_est - span_sum);
    prog->finish_task(scheme, path);
  }
}

static DensityStats density_proxy_files_mt(
    const ScannerSpec &spec, const std::vector<std::string> &paths,
    const std::vector<size_t> &sizes, int n_threads, DensityProgress *prog,
    const char *scheme) {
  if (paths.empty())
    return density_stats_from_counts(0, 0, 0, 0, 0);
  if (n_threads < 1)
    n_threads = 1;

  // Process one library FASTA at a time; parallelize across its sequences so a
  // huge bacteria/library.fna can use all threads (files are size-skewed).
  uint64_t distinct = 0, total = 0, sm_sum = 0, sm_max = 0, sm_count = 0;
  for (size_t i = 0; i < paths.size(); i++) {
    const size_t sz = i < sizes.size() ? sizes[i] : 0;
    density_proxy_one_file_mt(spec, paths[i], sz, n_threads, &distinct, &total,
                              &sm_sum, &sm_max, &sm_count, prog, scheme);
  }
  return density_stats_from_counts(distinct, total, sm_sum, sm_max, sm_count);
}

// Recursively collect library FASTA paths under root; estimate bases ≈ file size.
static void collect_library_fastas(const std::string &root,
                                   std::vector<std::string> *paths,
                                   std::vector<size_t> *sizes,
                                   size_t *n_bases_est) {
  DIR *dp = opendir(root.c_str());
  if (!dp)
    return;
  std::vector<std::string> subdirs;
  std::vector<std::string> files;
  while (dirent *ent = readdir(dp)) {
    if (ent->d_name[0] == '.')
      continue;
    const std::string name = ent->d_name;
    const std::string path = root + "/" + name;
    if (ent->d_type == DT_DIR) {
      subdirs.push_back(path);
      continue;
    }
    if (ent->d_type == DT_UNKNOWN || ent->d_type == DT_LNK) {
      DIR *probe = opendir(path.c_str());
      if (probe) {
        closedir(probe);
        subdirs.push_back(path);
        continue;
      }
    }
    if (is_library_fasta_name(name))
      files.push_back(path);
  }
  closedir(dp);
  std::sort(files.begin(), files.end());
  for (const auto &path : files) {
    paths->push_back(path);
    size_t sz = 0;
    std::ifstream probe(path, std::ios::binary | std::ios::ate);
    if (probe)
      sz = (size_t)probe.tellg();
    sizes->push_back(sz);
    *n_bases_est += sz;
  }
  std::sort(subdirs.begin(), subdirs.end());
  for (const auto &sub : subdirs)
    collect_library_fastas(sub, paths, sizes, n_bases_est);
}

// Particular density uses the Kraken DB library FASTAs (not query reads).
// Accepts a DB dir with library/, or a directory of .fna/.fa/.fasta/.faa.
static std::vector<std::string> list_particular_files(
    const std::string &dir, size_t *n_bases_est,
    std::vector<size_t> *sizes_out) {
  *n_bases_est = 0;
  std::vector<std::string> paths;
  std::vector<size_t> sizes;
  const std::string library_dir = dir + "/library";
  DIR *lib_probe = opendir(library_dir.c_str());
  if (lib_probe) {
    closedir(lib_probe);
    collect_library_fastas(library_dir, &paths, &sizes, n_bases_est);
  } else {
    collect_library_fastas(dir, &paths, &sizes, n_bases_est);
  }
  if (paths.empty())
    errx(EX_NOINPUT,
         "no library FASTA (.fna/.fa/.fasta/.faa) under %s "
         "(expected Kraken DB with library/, or a dir of reference FASTAs) — "
         "run ../benchmark/kraken2/download_data.sh",
         dir.c_str());
  std::cerr << "[*] particular DNA (DB library): " << paths.size()
            << " FASTA file(s) (~" << (*n_bases_est / 1000000) << " Mbp est) from "
            << dir << "\n";
  for (size_t i = 0; i < paths.size(); i++)
    std::cerr << "      " << paths[i] << "  (~" << (sizes[i] / 1000000)
              << " Mbp)\n";
  if (sizes_out)
    *sizes_out = std::move(sizes);
  return paths;
}

struct SchemeLabel {
  SketchingScheme sketching;
  Order order;
  const char *name;       // CLI / display (matches run_benchmark.sh)
  bool needs_gm_bin;      // Minimizer + GreedyMini binary (-g)
  bool needs_gm_dna;      // Minimizer + GreedyMini DNA (-G)
  bool needs_mod_gm;      // Mod + GreedyMini (-M)
};

// Same labels as ../benchmark/kraken2/run_benchmark.sh -S.
static const SchemeLabel kAllSchemes[] = {
    {SketchingScheme::MINIMIZER, Order::RANDOM_XOR, "Random", false, false,
     false},
    {SketchingScheme::MINIMIZER, Order::GM, "GreedyMini", true, false, false},
    {SketchingScheme::MINIMIZER, Order::GM, "GreedyMiniDNA", false, true,
     false},
    {SketchingScheme::MOD, Order::RANDOM_XOR, "ModRandom", false, false, false},
    {SketchingScheme::MOD, Order::GM, "ModGreedyMini", false, false, true},
    {SketchingScheme::MINIMIZER, Order::RANDOM_PCG, "RandomPcg", false, false,
     false},
};

// Default list matches run_benchmark.sh (GreedyMiniDNA is opt-in via -S).
static const char *kDefaultSchemeNames[] = {
    "Random", "GreedyMini", "ModRandom", "ModGreedyMini", "RandomPcg",
};

static std::string normalize_scheme_name(std::string name) {
  if (name == "RandomXor")
    return "Random";
  if (name == "ModRandomXor")
    return "ModRandom";
  if (name == "RandomPCG" || name == "randompcg" || name == "random_pcg")
    return "RandomPcg";
  return name;
}

static const SchemeLabel *find_scheme(const std::string &name) {
  const std::string n = normalize_scheme_name(name);
  for (const auto &s : kAllSchemes) {
    if (n == s.name)
      return &s;
  }
  return nullptr;
}

static void list_schemes_stderr() {
  std::cerr << "  Random | GreedyMini | GreedyMiniDNA | ModRandom | "
               "ModGreedyMini | RandomPcg\n";
}

struct RunRow {
  SchemeLabel cfg;
  BenchResult bench;
  double dens_expected = 0;    // random i.i.d. DNA (Monte Carlo E[density])
  double dens_particular = 0;  // measured on Kraken DB library FASTAs
  double dens_theory = 0;      // 2/(w+1) Minimizer, 1/w Mod
  double sm_mean_exp = 0;
  uint64_t sm_max_exp = 0;
  double sm_mean_part = 0;
  uint64_t sm_max_part = 0;
  bool ran;
};

// Kraken ClassifySequence cost proxy per scan step (one base):
//   nextmin_ns + density×lookup_ns_per_tok + resolve_ns_per_base
// density matches TOK_LOOKUP / (TOK_LOOKUP + TOK_REPEAT) trigger rate.
struct KrakenCostParams {
  double lookup_ns_per_tok = 0.0;
  double resolve_ns_per_base = 0.0;
  bool have_lookup = false;
  std::string source;
};

static bool parse_timing_double(const std::string &text, const char *key,
                                double *out) {
  const std::string needle = std::string(key) + "=";
  size_t pos = 0;
  while ((pos = text.find(needle, pos)) != std::string::npos) {
    pos += needle.size();
    char *end = nullptr;
    const double v = std::strtod(text.c_str() + pos, &end);
    if (end != text.c_str() + pos) {
      *out = v;
      return true;
    }
  }
  return false;
}

static bool load_kraken_cost_params(const std::string &path,
                                    KrakenCostParams *params) {
  FILE *fp = std::fopen(path.c_str(), "r");
  if (!fp)
    return false;
  std::string text;
  char buf[512];
  while (std::fgets(buf, sizeof buf, fp))
    text += buf;
  std::fclose(fp);

  double lookup_s = 0.0, resolve_s = 0.0, lookups = 0.0, bases = 0.0;
  const bool got_lookup_s = parse_timing_double(text, "lookup_s", &lookup_s);
  const bool got_lookups = parse_timing_double(text, "lookups", &lookups);
  parse_timing_double(text, "resolve_s", &resolve_s);
  parse_timing_double(text, "bases", &bases);

  if (got_lookup_s && got_lookups && lookups > 0.0) {
    params->lookup_ns_per_tok = lookup_s * 1e9 / lookups;
    params->have_lookup = true;
    params->source = path;
  }
  if (bases > 0.0)
    params->resolve_ns_per_base = resolve_s * 1e9 / bases;
  return params->have_lookup;
}

static void print_kraken_payoff(const std::vector<RunRow> &rows,
                                const KrakenCostParams &cost) {
  const RunRow *min_xor = nullptr;
  const RunRow *mod_xor = nullptr;
  const RunRow *min_pcg = nullptr;
  const RunRow *mod_pcg = nullptr;
  const RunRow *min_gm = nullptr;
  const RunRow *mod_gm = nullptr;
  for (const auto &row : rows) {
    if (!row.ran)
      continue;
    if (row.cfg.sketching == SketchingScheme::MINIMIZER &&
        row.cfg.order == Order::RANDOM_XOR)
      min_xor = &row;
    if (row.cfg.sketching == SketchingScheme::MOD &&
        row.cfg.order == Order::RANDOM_XOR)
      mod_xor = &row;
    if (row.cfg.sketching == SketchingScheme::MINIMIZER &&
        row.cfg.order == Order::RANDOM_PCG)
      min_pcg = &row;
    if (row.cfg.sketching == SketchingScheme::MOD &&
        row.cfg.order == Order::RANDOM_PCG)
      mod_pcg = &row;
    if (row.cfg.sketching == SketchingScheme::MINIMIZER &&
        row.cfg.order == Order::GM)
      min_gm = &row;
    if (row.cfg.sketching == SketchingScheme::MOD && row.cfg.order == Order::GM)
      mod_gm = &row;
  }

  const double lookup_ns = cost.lookup_ns_per_tok;
  const double resolve_ns = cost.resolve_ns_per_base;

  auto kraken_proxy = [&](const RunRow *row) {
    const double scan = row->bench.ns_per_op;
    const double lookup = row->dens_expected * lookup_ns;
    return scan + lookup + resolve_ns;
  };

  auto print_pair = [&](const char *label, const RunRow *mn,
                        const RunRow *md) {
    if (!mn || !md || mn->dens_expected <= 0.0 || md->dens_expected <= 0.0)
      return;
    const double scan_mn = mn->bench.ns_per_op;
    const double scan_md = md->bench.ns_per_op;
    const double scan_overhead = scan_md / scan_mn;
    const double density_ratio = md->dens_expected / mn->dens_expected;
    const double sparsity_gain = mn->dens_expected / md->dens_expected;
    const double lookup_mn = mn->dens_expected * lookup_ns;
    const double lookup_md = md->dens_expected * lookup_ns;
    const double lookup_savings = lookup_mn - lookup_md;
    const double scan_penalty = scan_md - scan_mn;
    const double total_mn = kraken_proxy(mn);
    const double total_md = kraken_proxy(md);
    const double net_delta = total_md - total_mn;
    const double total_ratio = total_md / total_mn;
    const double lookups_per_mbp_mn = mn->dens_expected * 1e6;
    const double lookups_per_mbp_md = md->dens_expected * 1e6;
    const double dens_gap = mn->dens_expected - md->dens_expected;
    const double breakeven_lookup =
        dens_gap > 0.0 ? scan_penalty / dens_gap : 0.0;
    const bool pays_off = total_md < total_mn;

    std::cout << "\n=== Kraken classify cost proxy: " << label << " ===\n"
              << "  model per base: nextmin_ns + density×lookup_ns + resolve_ns\n"
              << "    (density = expected TOK_LOOKUP rate on random DNA;"
                 " lower = sparser = good)\n";
    if (cost.have_lookup) {
      std::cout << "  lookup_ns/tok: " << lookup_ns;
      if (!cost.source.empty())
        std::cout << "  (from " << cost.source << ")";
      std::cout << "\n";
    } else {
      std::cout << "  lookup_ns/tok: (not set — pass -L or -F from classify "
                   "PHASE_TIMING)\n";
    }
    if (resolve_ns > 0.0)
      std::cout << "  resolve_ns/base: " << resolve_ns << "\n";

    std::cout << "  density expected (TOK_LOOKUP / scan step):\n"
              << "    Minimizer: " << mn->dens_expected << "   Mod: "
              << md->dens_expected << "   (" << (100.0 * (1.0 - density_ratio))
              << "% fewer lookups/Mbp, sparsity " << sparsity_gain << "x)\n"
              << "  density particular (DB library FASTAs):\n"
              << "    Minimizer: " << mn->dens_particular << "   Mod: "
              << md->dens_particular << "\n"
              << "  dens_theory: Minimizer " << mn->dens_theory << "   Mod "
              << md->dens_theory << "\n"
              << "  lookups/Mbp (expected): Minimizer " << lookups_per_mbp_mn
              << "   Mod " << lookups_per_mbp_md << "\n"
              << "  nextmin ns/step: Minimizer " << scan_mn << "   Mod "
              << scan_md << "   (overhead " << scan_overhead << "x)\n";

    if (cost.have_lookup) {
      std::cout << "  component ns/base (scan + lookup + resolve = total):\n"
                << "    Minimizer  " << scan_mn << " + " << lookup_mn << " + "
                << resolve_ns << " = " << total_mn << "\n"
                << "    Mod        " << scan_md << " + " << lookup_md << " + "
                << resolve_ns << " = " << total_md << "\n"
                << "  lookup savings from sparsity: " << lookup_savings
                << " ns/base\n"
                << "  scan penalty from Mod: " << scan_penalty << " ns/base\n"
                << "  net (Mod − Min): " << net_delta << " ns/base  (ratio "
                << total_ratio << "x)\n"
                << "  breakeven lookup_ns/tok: " << breakeven_lookup
                << "  (Mod wins above this)\n";
    } else {
      std::cout << "  scan-only net (Mod − Min): " << scan_penalty
                << " ns/base  (overhead " << scan_overhead << "x)\n"
                << "  breakeven lookup_ns/tok: " << breakeven_lookup
                << "  (Mod wins when lookup cost exceeds this)\n";
    }

    std::cout << "  Verdict: ";
    if (!cost.have_lookup) {
      std::cout << "lookup phase not modeled — supply -L or -F to include "
                   "GetBatch; breakeven above assumes scan+lookup tradeoff.\n";
    } else if (pays_off) {
      std::cout << "MOD pays off — " << lookup_savings
                << " ns/base lookup savings exceed " << scan_penalty
                << " ns/base scan penalty (total " << total_ratio << "x).\n";
    } else {
      std::cout << "MOD does NOT pay off — " << scan_penalty
                << " ns/base scan penalty exceeds " << lookup_savings
                << " ns/base lookup savings (total " << total_ratio << "x).\n";
    }
  };

  print_pair("RandomXor", min_xor, mod_xor);
  print_pair("RandomPcg", min_pcg, mod_pcg);
  print_pair("GreedyMini", min_gm, mod_gm);
}

static void usage() {
  std::cerr
      << "Usage: benchmark_minimizers [options]\n"
      << "  -S NAME   Scheme(s) to run (repeatable / space-separated):\n"
      << "            Random | GreedyMini | GreedyMiniDNA | ModRandom |\n"
      << "            ModGreedyMini | RandomPcg\n"
      << "            (default: Random GreedyMini ModRandom ModGreedyMini "
         "RandomPcg)\n"
      << "  -g FILE   GreedyMini binary ℓ-mer table for GreedyMini\n"
      << "  -G FILE   GreedyMini DNA alphabet table for GreedyMiniDNA\n"
      << "  -M FILE   Mod+GreedyMini binary t-mer table (also --mod-gm)\n"
      << "  -l INT    minimizer length ℓ (default: 31)\n"
      << "  -k INT    k-mer length (default: l+4 ⇒ w=5)\n"
      << "  -n INT    DNA length ≈ n + k for expected density (default: 50000000)\n"
      << "  -b INT    batch size for latency samples (default: 1024)\n"
      << "  -j INT    parallel jobs for density (schemes × files; "
         "default: OMP threads)\n"
      << "  -d DIR    Kraken DB dir for particular density "
         "(library/ FASTAs)\n"
      << "  -o FILE   write density table to FILE (.csv or text summary)\n"
      << "  --seed-idx INT  seed index column for CSV (default: 0)\n"
      << "  -L FLOAT  ns per TOK_LOOKUP (GetBatch); from lookup_s/lookups×1e9\n"
      << "  -F FILE   classify timing file (nextminimizer_timing_*.txt)\n"
      << "  -R FLOAT  resolve phase ns per base (default: from -F or 0)\n"
      << "  -T HEX    toggle mask\n"
      << "  -s INT    minimizer-spaces s (default: 7, same as kraken2-build NT)\n"
      << "  --seed INT  RNG seed (default: 1)\n"
      << "Density: theory (2/(w+1) Minimizer, 1/w Mod), expected, "
         "particular (DB library).\n"
      << "Super-mer: mean/max run length of identical keys (scan steps).\n"
      << "Density factor = density × (w+1)  (random minimizer theory ⇒ 2).\n";
  exit(EX_USAGE);
}

static void print_density_summary(std::ostream &out, int l, int k,
                                  size_t dna_len, size_t particular_bases_est,
                                  const std::string &data_dir,
                                  const std::vector<RunRow> &rows) {
  const int w = k - l + 1;
  out << "=== Minimizer density summary ===\n"
      << "  l=" << l << "  k=" << k << "  w=" << w << "\n"
      << "  theory: Minimizer 2/(w+1)=" << (2.0 / (w + 1))
      << "  Mod 1/w=" << (1.0 / w) << "\n"
      << "  density factor = density × (w+1)  "
         "(random minimizer theory ⇒ 2)\n"
      << "  expected DNA: random i.i.d. len=" << dna_len << "\n"
      << "  particular DNA (DB library): ~" << particular_bases_est
      << " bases (est) from " << data_dir << "\n\n"
      << "Scheme                   dens_theory  dens_expected  dens_particular"
         "  fac_theory  fac_expected  fac_particular  ns/minimizer"
         "  sm_mean_exp  sm_max_exp  sm_mean_part  sm_max_part\n"
      << "-----------------------  -----------  -------------  ---------------"
         "  ----------  ------------  --------------  ------------"
         "  -----------  ----------  ------------  -----------\n";
  for (const auto &row : rows) {
    if (!row.ran)
      continue;
    char buf[480];
    std::snprintf(
        buf, sizeof buf,
        "%-23s  %11.6f  %13.6f  %15.6f  %10.4f  %12.4f  %14.4f  %12.3f"
        "  %11.4f  %10llu  %12.4f  %11llu\n",
        row.cfg.name, row.dens_theory, row.dens_expected, row.dens_particular,
        density_factor(row.dens_theory, w),
        density_factor(row.dens_expected, w),
        density_factor(row.dens_particular, w), row.bench.ns_per_op,
        row.sm_mean_exp, (unsigned long long)row.sm_max_exp, row.sm_mean_part,
        (unsigned long long)row.sm_max_part);
    out << buf;
  }
  out << "\n"
      << "Notes:\n"
      << "  dens_theory     = 2/(w+1) Minimizer; 1/w Mod\n"
      << "  dens_expected   = change-rate on random i.i.d. ACGT (Monte Carlo)\n"
      << "  dens_particular = change-rate on DB library FASTAs under -d\n"
      << "  fac_*           = dens_* × (w+1); lower than 2 ⇒ better than "
         "random minimizer\n"
      << "  sm_mean/max     = super-mer run length (identical keys, scan steps)\n"
      << "  lower density ⇒ fewer distinct lookups / Mbp\n";
}

static bool ends_with_csv(const std::string &path) {
  return path.size() >= 4 && path.compare(path.size() - 4, 4, ".csv") == 0;
}

static void print_density_csv(std::ostream &out, int l, int k, uint64_t toggle,
                              int seed_idx, const std::vector<RunRow> &rows) {
  const int w = k - l + 1;
  out << "scheme,seed_idx,toggle_hex,l,k,w,dens_expected,"
         "dens_particular,fac_expected,fac_particular,"
         "ns_per_minimizer,sm_mean_exp,sm_max_exp,sm_mean_part,sm_max_part\n";
  for (const auto &row : rows) {
    if (!row.ran)
      continue;
    char buf[640];
    std::snprintf(
        buf, sizeof buf,
        "%s,%d,0x%llx,%d,%d,%d,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%llu,%.6f,%llu\n",
        row.cfg.name, seed_idx, (unsigned long long)toggle, l, k, w,
        row.dens_expected, row.dens_particular,
        density_factor(row.dens_expected, w),
        density_factor(row.dens_particular, w), row.bench.ns_per_op,
        row.sm_mean_exp, (unsigned long long)row.sm_max_exp, row.sm_mean_part,
        (unsigned long long)row.sm_max_part);
    out << buf;
  }
}

static void print_uhs_hits(const char *label, const BenchResult &r,
                           const GreedyMiniTable *table) {
  if (!r.has_uhs_stats)
    return;
  const RankHitStats &hs = r.uhs_stats;
  const double tot = hs.total() ? (double)hs.total() : 1.0;
  std::cout << "\n" << label << " Rank() UHS hits (timed NextMinimizer):\n"
            << "  Rank calls     " << hs.total() << "\n"
            << "  UHS (finite)   " << hs.uhs << "  ("
            << (100.0 * hs.uhs / tot) << "%)\n"
            << "  non-UHS (inf)  " << hs.non_uhs << "  ("
            << (100.0 * hs.non_uhs / tot) << "%)\n";
  if (table && table->size())
    std::cout << "  table UHS dens " << table->uhs_size() << "/"
              << table->size() << "  ("
              << (100.0 * table->uhs_size() / table->size()) << "%)\n";
}

static void print_rank_path(const char *label, size_t n_mers, double nspc,
                            const PathDurationStats &uhs,
                            const PathDurationStats &non_uhs) {
  const uint64_t path_tot = uhs.n_calls + non_uhs.n_calls;
  std::cout << "\n" << label << " Rank() duration by path "
               "(stateful replay, RDTSC; ns/Rank):\n"
            << "  mers           " << n_mers << "  calibrated " << nspc
            << " ns/cycle\n"
            << "  UHS calls      " << uhs.n_calls;
  if (path_tot)
    std::cout << "  (" << (100.0 * uhs.n_calls / path_tot) << "%)";
  std::cout << "  mean=" << uhs.mean_ns << " ns\n"
            << "  non-UHS calls  " << non_uhs.n_calls;
  if (path_tot)
    std::cout << "  (" << (100.0 * non_uhs.n_calls / path_tot) << "%)";
  std::cout << "  mean=" << non_uhs.mean_ns << " ns\n";
  if (non_uhs.mean_ns > 0)
    std::cout << "  UHS / non-UHS  " << (uhs.mean_ns / non_uhs.mean_ns)
              << "x\n";
  print_dist("UHS", uhs.sample_dist, "ns/Rank");
  print_dist("non-UHS", non_uhs.sample_dist, "ns/Rank");
}

int main(int argc, char **argv) {
  std::string gm_file;
  std::string gm_dna_file;
  std::string mod_gm_file;
  int l = 31;
  int k = 0;
  size_t n_bases = 50000000;
  size_t batch = 1024;
  uint64_t toggle = DEFAULT_TOGGLE_MASK;
  int spaces = 7;
  uint64_t seed = 1;
  int n_jobs = 0;  // 0 = omp_get_max_threads()
  std::string timing_file;
  std::string data_dir;
  std::string out_file;
  int seed_idx = 0;
  double lookup_ns_cli = -1.0;
  double resolve_ns_cli = -1.0;
  KrakenCostParams kraken_cost;
  std::vector<const SchemeLabel *> selected;
  bool scheme_filter_set = false;

  auto append_scheme_arg = [&](const char *raw) {
    std::string chunk(raw);
    size_t pos = 0;
    while (pos < chunk.size()) {
      while (pos < chunk.size() &&
             (chunk[pos] == ' ' || chunk[pos] == '\t'))
        pos++;
      if (pos >= chunk.size())
        break;
      size_t end = pos;
      while (end < chunk.size() && chunk[end] != ' ' && chunk[end] != '\t')
        end++;
      const std::string name = chunk.substr(pos, end - pos);
      pos = end;
      if (!scheme_filter_set) {
        selected.clear();
        scheme_filter_set = true;
      }
      const SchemeLabel *spec = find_scheme(name);
      if (!spec) {
        std::cerr << "Error: unknown scheme '" << name << "'.\n";
        list_schemes_stderr();
        exit(EX_USAGE);
      }
      selected.push_back(spec);
    }
  };

  for (int i = 1; i < argc; i++) {
    if (!std::strcmp(argv[i], "-S") && i + 1 < argc) {
      append_scheme_arg(argv[++i]);
      while (i + 1 < argc && argv[i + 1][0] != '-')
        append_scheme_arg(argv[++i]);
    } else if (!std::strcmp(argv[i], "-g") && i + 1 < argc)
      gm_file = argv[++i];
    else if (!std::strcmp(argv[i], "-G") && i + 1 < argc)
      gm_dna_file = argv[++i];
    else if ((!std::strcmp(argv[i], "-M") ||
              !std::strcmp(argv[i], "--mod-gm")) &&
             i + 1 < argc)
      mod_gm_file = argv[++i];
    else if (!std::strcmp(argv[i], "-l") && i + 1 < argc)
      l = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "-k") && i + 1 < argc)
      k = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "-n") && i + 1 < argc)
      n_bases = (size_t)std::strtoull(argv[++i], nullptr, 10);
    else if (!std::strcmp(argv[i], "-b") && i + 1 < argc)
      batch = (size_t)std::strtoull(argv[++i], nullptr, 10);
    else if (!std::strcmp(argv[i], "-j") && i + 1 < argc)
      n_jobs = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "-d") && i + 1 < argc)
      data_dir = argv[++i];
    else if (!std::strcmp(argv[i], "-o") && i + 1 < argc)
      out_file = argv[++i];
    else if (!std::strcmp(argv[i], "-L") && i + 1 < argc)
      lookup_ns_cli = std::strtod(argv[++i], nullptr);
    else if (!std::strcmp(argv[i], "-F") && i + 1 < argc)
      timing_file = argv[++i];
    else if (!std::strcmp(argv[i], "-R") && i + 1 < argc)
      resolve_ns_cli = std::strtod(argv[++i], nullptr);
    else if (!std::strcmp(argv[i], "-T") && i + 1 < argc)
      toggle = std::strtoull(argv[++i], nullptr, 0);
    else if (!std::strcmp(argv[i], "-s") && i + 1 < argc)
      spaces = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--seed") && i + 1 < argc)
      seed = std::strtoull(argv[++i], nullptr, 10);
    else if (!std::strcmp(argv[i], "--seed-idx") && i + 1 < argc)
      seed_idx = std::atoi(argv[++i]);
    else
      usage();
  }
  if (l < 1 || l > 31 || batch < 1)
    usage();
  if (n_jobs < 0)
    errx(EX_USAGE, "-j must be >= 0 (0 = auto)");
  if (k <= 0)
    k = l + 4;
  if (k < l)
    errx(EX_USAGE, "k must be >= l");
  if (spaces < 0)
    errx(EX_USAGE, "-s minimizer-spaces must be >= 0");
  if (spaces > 0 && spaces > l / 4)
    errx(EX_USAGE, "-s %d exceeds max ℓ/4 = %d", spaces, l / 4);
  const uint64_t spaced_seed_mask =
      MakeSpacedSeedMask(l, spaces, BITS_PER_CHAR_DNA);
  if (data_dir.empty())
    errx(EX_USAGE,
         "particular density requires -d DIR "
         "(Kraken DB with library/ FASTAs)");
  if (scheme_filter_set && selected.empty())
    errx(EX_USAGE, "-S requires at least one scheme");
  if (!scheme_filter_set) {
    for (const char *name : kDefaultSchemeNames) {
      const SchemeLabel *spec = find_scheme(name);
      if (spec)
        selected.push_back(spec);
    }
  }

  bool need_gm_bin = false, need_gm_dna = false, need_mod_gm = false;
  for (const SchemeLabel *s : selected) {
    need_gm_bin = need_gm_bin || s->needs_gm_bin;
    need_gm_dna = need_gm_dna || s->needs_gm_dna;
    need_mod_gm = need_mod_gm || s->needs_mod_gm;
  }
  if (need_gm_bin && gm_file.empty())
    errx(EX_USAGE, "GreedyMini requires -g FILE");
  if (need_gm_dna && gm_dna_file.empty())
    errx(EX_USAGE, "GreedyMiniDNA requires -G FILE");
  if (need_mod_gm && mod_gm_file.empty())
    errx(EX_USAGE, "ModGreedyMini requires -M FILE");

  if (!timing_file.empty())
    load_kraken_cost_params(timing_file, &kraken_cost);
  if (lookup_ns_cli >= 0.0) {
    kraken_cost.lookup_ns_per_tok = lookup_ns_cli;
    kraken_cost.have_lookup = true;
    if (kraken_cost.source.empty())
      kraken_cost.source = "-L";
  }
  if (resolve_ns_cli >= 0.0)
    kraken_cost.resolve_ns_per_base = resolve_ns_cli;

  std::shared_ptr<const GreedyMiniTable> table_bin;
  if (need_gm_bin) {
    table_bin =
        MinimizerScanner::shared_gm_table(gm_file, GmAlphabet::BINARY);
    if (!table_bin || table_bin->empty())
      errx(EX_NOINPUT, "failed to load %s", gm_file.c_str());
    if (table_bin->gm_k() > l)
      errx(EX_DATAERR, "gm_k=%d > l=%d", table_bin->gm_k(), l);
  }

  std::shared_ptr<const GreedyMiniTable> table_dna;
  if (need_gm_dna) {
    table_dna =
        MinimizerScanner::shared_gm_table(gm_dna_file, GmAlphabet::DNA);
    if (!table_dna || table_dna->empty())
      errx(EX_NOINPUT, "failed to load %s", gm_dna_file.c_str());
    if (table_dna->gm_k() > l)
      errx(EX_DATAERR, "gm_k=%d > l=%d", table_dna->gm_k(), l);
  }

  std::shared_ptr<const GreedyMiniTable> mod_table;
  if (need_mod_gm) {
    mod_table =
        MinimizerScanner::shared_gm_table(mod_gm_file, GmAlphabet::BINARY);
    if (!mod_table || mod_table->empty())
      errx(EX_NOINPUT, "failed to load %s", mod_gm_file.c_str());
  }

  size_t particular_bases_est = 0;
  std::vector<size_t> particular_sizes;
  const std::vector<std::string> particular_paths =
      list_particular_files(data_dir, &particular_bases_est, &particular_sizes);

  std::mt19937_64 rng(seed);
  const size_t dna_len = n_bases + (size_t)k;
  std::string dna(dna_len, 'A');
  const char *bases = "ACGT";
  for (size_t i = 0; i < dna_len; i++)
    dna[i] = bases[rng() & 3ULL];

  std::vector<RunRow> rows(selected.size());
  std::vector<std::string> errors(selected.size());

  std::cerr << "[*] latency microbench (serial)...\n";
  // Latency microbench stays serial so ns/minimizer is not contending.
  for (size_t i = 0; i < selected.size(); i++) {
    const SchemeLabel &cfg = *selected[i];
    RunRow &row = rows[i];
    row.cfg = cfg;
    row.ran = false;
    ScannerSpec spec;
    spec.k = k;
    spec.l = l;
    spec.toggle = toggle;
    spec.spaced_seed_mask = spaced_seed_mask;
    spec.sketching = cfg.sketching;
    spec.order = cfg.order;
    if (cfg.needs_gm_bin) {
      spec.gm_path = gm_file;
      spec.gm_tbl = table_bin;
    } else if (cfg.needs_gm_dna) {
      spec.gm_path = gm_dna_file;
      spec.gm_tbl = table_dna;
    } else if (cfg.needs_mod_gm) {
      spec.gm_path = mod_gm_file;
      spec.gm_tbl = mod_table;
    }
    try {
      std::cerr << "    " << cfg.name << " NextMinimizer timing...\n";
      MinimizerScanner sc = make_scanner(spec);
      row.bench = scan_once(sc, dna, batch);
      row.dens_theory = density_theory(cfg.sketching, k, l);
      row.ran = true;
    } catch (const std::exception &ex) {
      errx(EX_DATAERR, "%s: %s", cfg.name, ex.what());
    }
  }

  if (n_jobs == 0)
    n_jobs = omp_get_max_threads();
  if (n_jobs < 1)
    n_jobs = 1;
  const int n_schemes = (int)selected.size();
  int scheme_threads = n_jobs < n_schemes ? n_jobs : n_schemes;
  if (scheme_threads < 1)
    scheme_threads = 1;
  int file_threads = n_jobs / scheme_threads;
  if (file_threads < 1)
    file_threads = 1;
  if (n_schemes == 1)
    file_threads = n_jobs;
  omp_set_nested(scheme_threads > 1 && file_threads > 1 ? 1 : 0);

  uint64_t bytes_per_scheme = 0;
  for (size_t sz : particular_sizes)
    bytes_per_scheme += sz;
  DensityProgress progress;
  progress.total_bytes = bytes_per_scheme * (uint64_t)n_schemes;
  progress.total_tasks = n_schemes * (int)particular_paths.size();
  progress.t0 = std::chrono::steady_clock::now();

  std::cerr << "[*] density: expected (random DNA) + particular (DB library)\n"
            << "    schemes=" << n_schemes << "  files/scheme="
            << particular_paths.size() << "  jobs=" << n_jobs
            << " (scheme_threads=" << scheme_threads
            << " seq_threads/file=" << file_threads << ")\n";

  // Expected + particular density: parallel over schemes (and files if -j allows).
#pragma omp parallel for num_threads(scheme_threads) schedule(dynamic)
  for (int i = 0; i < n_schemes; i++) {
    if (!rows[(size_t)i].ran)
      continue;
    const SchemeLabel &cfg = *selected[(size_t)i];
    ScannerSpec spec;
    spec.k = k;
    spec.l = l;
    spec.toggle = toggle;
    spec.spaced_seed_mask = spaced_seed_mask;
    spec.sketching = cfg.sketching;
    spec.order = cfg.order;
    if (cfg.needs_gm_bin) {
      spec.gm_path = gm_file;
      spec.gm_tbl = table_bin;
    } else if (cfg.needs_gm_dna) {
      spec.gm_path = gm_dna_file;
      spec.gm_tbl = table_dna;
    } else if (cfg.needs_mod_gm) {
      spec.gm_path = mod_gm_file;
      spec.gm_tbl = mod_table;
    }
    try {
      MinimizerScanner sc = make_scanner(spec);
      {
#pragma omp critical(bench_status)
        std::cerr << "\n[*] " << cfg.name << " expected density (random DNA)...\n";
      }
      const DensityStats exp = density_proxy(sc, dna);
      rows[(size_t)i].dens_expected = exp.density;
      rows[(size_t)i].sm_mean_exp = exp.sm_mean;
      rows[(size_t)i].sm_max_exp = exp.sm_max;
      {
#pragma omp critical(bench_status)
        std::cerr << "[*] " << cfg.name << " particular density (DB library)...\n";
      }
      const DensityStats part = density_proxy_files_mt(
          spec, particular_paths, particular_sizes, file_threads, &progress,
          cfg.name);
      rows[(size_t)i].dens_particular = part.density;
      rows[(size_t)i].sm_mean_part = part.sm_mean;
      rows[(size_t)i].sm_max_part = part.sm_max;
    } catch (const std::exception &ex) {
      errors[(size_t)i] = ex.what();
      rows[(size_t)i].ran = false;
    }
  }
  progress.finish_all();
  for (size_t i = 0; i < errors.size(); i++) {
    if (!errors[i].empty())
      errx(EX_DATAERR, "%s: %s", selected[i]->name, errors[i].c_str());
  }

  const double nspc = ns_per_cycle();
  const uint64_t lmer_mask =
      (l >= 32) ? ~0ULL : ((1ULL << (2 * l)) - 1ULL);
  const int w = k - l + 1;

  std::cout << "=== MinimizerScanner microbenchmark (Sketching × Order) ===\n"
            << "  l=" << l << "  k=" << k << "  w=" << w
            << "  dna_len=" << dna_len << "  batch=" << batch
            << "  toggle=0x" << std::hex << toggle
            << "  spaces=" << std::dec << spaces
            << "  spaced=0x" << std::hex << spaced_seed_mask << std::dec
            << "  jobs=" << n_jobs << " (schemes=" << scheme_threads
            << " files=" << file_threads << ")\n"
            << "  schemes:";
  for (const SchemeLabel *s : selected)
    std::cout << " " << s->name;
  std::cout << "\n"
            << "  dens_theory: Minimizer 2/(w+1)=" << (2.0 / (w + 1))
            << "  Mod 1/w=" << (1.0 / w)
            << "  (factor = dens×(w+1); random ⇒ 2)\n"
            << "  particular (DB library): " << particular_paths.size()
            << " file(s), ~" << particular_bases_est << " bases from "
            << data_dir << "\n";
  if (table_bin)
    std::cout << "  GreedyMini:    " << gm_file
              << "  alphabet=binary"
              << "  gm_k=" << table_bin->gm_k()
              << "  UHS=" << table_bin->uhs_size() << "/" << table_bin->size()
              << "\n";
  if (table_dna)
    std::cout << "  GreedyMiniDNA: " << gm_dna_file
              << "  alphabet=dna"
              << "  gm_k=" << table_dna->gm_k()
              << "  UHS=" << table_dna->uhs_size() << "/" << table_dna->size()
              << "\n";
  if (mod_table)
    std::cout << "  ModGreedyMini: " << mod_gm_file
              << "  alphabet=binary"
              << "  gm_k=" << mod_table->gm_k()
              << "  UHS=" << mod_table->uhs_size() << "/" << mod_table->size()
              << "\n";
  std::cout << "  BMI2 pext: "
#if defined(__BMI2__)
            << "yes"
#else
            << "no"
#endif
            << "\n\n"
            << "Scheme                   ns/minimizer   #calls   dens_theory"
               "  dens_exp  dens_part  fac_exp  fac_part  sm_mean  sm_max"
               "  checksum\n"
            << "-----------------------  -------------  ----------  -----------"
               "  --------  ---------  -------  --------  -------  ------"
               "  ----------------\n";
  for (const auto &row : rows) {
    if (!row.ran)
      continue;
    std::printf(
        "%-23s  %13.3f  %10llu  %11.6f  %8.4f  %9.4f  %7.4f  %8.4f  %7.3f"
        "  %6llu  %016llx\n",
        row.cfg.name, row.bench.ns_per_op, (unsigned long long)row.bench.n_ops,
        row.dens_theory, row.dens_expected, row.dens_particular,
        density_factor(row.dens_expected, w),
        density_factor(row.dens_particular, w), row.sm_mean_exp,
        (unsigned long long)row.sm_max_exp,
        (unsigned long long)row.bench.checksum);
  }

  std::cout << "\n";
  if (!out_file.empty()) {
    std::ofstream ofs(out_file);
    if (!ofs)
      errx(EX_CANTCREAT, "cannot write %s", out_file.c_str());
    if (ends_with_csv(out_file))
      print_density_csv(ofs, l, k, toggle, seed_idx, rows);
    else
      print_density_summary(ofs, l, k, dna_len, particular_bases_est, data_dir,
                            rows);
  }

  print_kraken_payoff(rows, kraken_cost);

  std::cout << "\nLatency distribution (ns/minimizer from batches of " << batch
            << "):\n";
  for (const auto &row : rows) {
    if (row.ran)
      print_dist(row.cfg.name, row.bench.dist, "ns/minimizer");
  }

  for (const auto &row : rows) {
    if (row.ran)
      print_phase_breakdown(row.cfg.name, row.bench.phases, nspc,
                            row.bench.ns_per_op);
  }

  for (const auto &row : rows) {
    if (!row.ran || row.cfg.order != Order::GM)
      continue;
    const GreedyMiniTable *t = nullptr;
    if (row.cfg.needs_gm_bin)
      t = table_bin.get();
    else if (row.cfg.needs_gm_dna)
      t = table_dna.get();
    else if (row.cfg.needs_mod_gm)
      t = mod_table.get();
    print_uhs_hits(row.cfg.name, row.bench, t);
  }

  // Isolated Rank() replay (GreedyMini only): ℓ-mers for MINIMIZER, t-mers for MOD.
  if (table_bin) {
    MinimizerScanner sc_log(k, l, spaced_seed_mask, true, toggle,
                            CURRENT_REVCOM_VERSION, SketchingScheme::MINIMIZER,
                            Order::GM, gm_file, table_bin);
    std::vector<uint64_t> ranked_lmers;
    ranked_lmers.reserve(dna_len);
    sc_log.SetRankedLmerLog(&ranked_lmers);
    sc_log.LoadSequence(dna);
    while (sc_log.NextMinimizer()) {
    }
    sc_log.SetRankedLmerLog(nullptr);

    GMRanker gm_ranker(l, BITS_PER_CHAR_DNA, toggle, lmer_mask, true,
                       table_bin, spaced_seed_mask);
    PathDurationStats uhs, non;
    const size_t stride = ranked_lmers.size() > 2000000 ? 32 : 8;
    time_rank_paths(ranked_lmers, gm_ranker, stride, 200000, nspc, uhs, non);
    print_rank_path("GreedyMini", ranked_lmers.size(), nspc, uhs, non);
  }

  if (table_dna) {
    MinimizerScanner sc_log(k, l, spaced_seed_mask, true, toggle,
                            CURRENT_REVCOM_VERSION, SketchingScheme::MINIMIZER,
                            Order::GM, gm_dna_file, table_dna);
    std::vector<uint64_t> ranked_lmers;
    ranked_lmers.reserve(dna_len);
    sc_log.SetRankedLmerLog(&ranked_lmers);
    sc_log.LoadSequence(dna);
    while (sc_log.NextMinimizer()) {
    }
    sc_log.SetRankedLmerLog(nullptr);

    GMRanker gm_ranker(l, BITS_PER_CHAR_DNA, toggle, lmer_mask, true,
                       table_dna, spaced_seed_mask);
    PathDurationStats uhs, non;
    const size_t stride = ranked_lmers.size() > 2000000 ? 32 : 8;
    time_rank_paths(ranked_lmers, gm_ranker, stride, 200000, nspc, uhs, non);
    print_rank_path("GreedyMiniDNA", ranked_lmers.size(), nspc, uhs, non);
  }

  if (mod_table) {
    MinimizerScanner sc_mod(k, l, spaced_seed_mask, true, toggle,
                            CURRENT_REVCOM_VERSION, SketchingScheme::MOD,
                            Order::GM, "", mod_table);
    std::vector<uint64_t> ranked_tmers;
    ranked_tmers.reserve(dna_len);
    sc_mod.SetRankedTmerLog(&ranked_tmers);
    sc_mod.LoadSequence(dna);
    while (sc_mod.NextMinimizer()) {
    }
    sc_mod.SetRankedTmerLog(nullptr);

    const int t = ModSampler::kDefaultR + ((l - ModSampler::kDefaultR) % w);
    const int t_bits = 2 * t;
    const uint64_t t_mask =
        (t_bits >= 64) ? ~0ULL : ((1ULL << t_bits) - 1ULL);
    const size_t stride = ranked_tmers.size() > 2000000 ? 32 : 8;

    GMRanker t_gm(t, BITS_PER_CHAR_DNA, toggle & t_mask, t_mask, true,
                  mod_table, 0, CURRENT_REVCOM_VERSION);
    PathDurationStats uhs, non;
    time_rank_paths(ranked_tmers, t_gm, stride, 200000, nspc, uhs, non);
    print_rank_path("ModGreedyMini", ranked_tmers.size(), nspc, uhs, non);
  }

  std::cout << "\n";
  print_density_summary(std::cout, l, k, dna_len, particular_bases_est,
                        data_dir, rows);
  if (!out_file.empty())
    std::cout << "\n[*] density summary written to " << out_file << "\n";

  return 0;
}
