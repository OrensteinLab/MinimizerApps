// Microbenchmark: ways of computing a GreedyMini signature per base, as in
// kmc's splitter, plus a realistic consumer (sliding-window minimum over
// W = L-p+1 signature values, counting minimizer changes = super-k-mers).
//
// Variants (all driven by a precomputed 2-bit ACGT stream, one insert per base):
//   1 default      kmc CMmer LUT_DEFAULT: norm[str] (lexicographic canonical)
//   2 dna          kmc CMmer LUT_DNA:     gm_dna_sort_lut[str] (4^p u32 table)
//   3 binary       kmc CMmer LUT_BINARY:  rolling bin_key + str_rc + tie-break
//   4a shared      gm_center::Order as kraken2 uses it: rolled Sigma plane
//                  (roll_in_low), Gamma via TwoBitPlanes::gamma on the 2-bit
//                  word, value_of_planes (p == gm_k -> center offset 0)
//   4b shared-both rolling Sigma/Gamma planes of BOTH strands + value_of_both_strands
//   4c shared-rank rank_of_sigma only (rank-first, ties -> leftmost; NOT a
//                  total order, lower bound of the shared path)
//   5 lean-lb      HAND VARIANT, not an implementation of kmc's order: rolling
//                  bin_key + rank LUT only, no tie-break/rc (lower bound of
//                  variant 3; ties -> leftmost, so not strand symmetric)
//
// Build (from this directory; kmc's own flags are -O3 -fsigned-char, no -march;
// build a second binary with -march=native to see the BMI2 (pext) effect):
//   R=../../..; g++ -O3 -std=c++17 -fsigned-char [-march=native] \
//     -o bench_signature_impls bench_signature_impls.cpp \
//     $R/kmc/kmc_api/mmer.cpp $R/kmc/kmc_api/custom_order.cpp
// Run (pin to one core):
//   taskset -c 5 ./bench_signature_impls P TABLE.gm TABLE_dna.gm DNA_SEED L SRC NBASES REPS
//   SRC = "random" or a FASTQ path (ACGT bases of the sequence lines, N dropped,
//   reads concatenated). TABLE_dna.gm: gm_sigma_extend --seed S -o X TABLE.gm.
//   DNA_SEED = the seed baked into the _dna.gm (also given to the binary LUT so
//   tie-breaks are comparable; kmc's default seed for a binary table is 0).

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <immintrin.h>

#define class struct  // CMmer has no explicit private: -> expose current_val (what kmc compares)
#include "../../../kmc/kmc_api/mmer.h"
#undef class
#include "../../../kmc/kmc_api/custom_order.h"
#include "../../../common/minimizer-schemes/greedymini_center.h"

using clk = std::chrono::steady_clock;
static int P, W;

// ---------------------------------------------------------------- signatures
template <int Mode> struct SigKmc {
  CMmer m;
  explicit SigKmc(int p) : m(p) {}
  inline uint64_t next(uint8_t s) { m.insert_sym<Mode>(s); return m.current_val; }
};

struct SigShared {  // 4a
  gm_center::Order ord; gm_center::TwoBitPlanes<true> pl; uint64_t word = 0, sigma = 0, wmask, smask;
  SigShared(std::shared_ptr<const GreedyMiniTable> t, int p, uint64_t tog)
      : ord(t, p, tog), pl(p, false), wmask((1ULL << (2 * p)) - 1), smask(pl.mask()) {}
  inline uint64_t next(uint8_t s) {
    word = ((word << 2) | s) & wmask;
    sigma = gm_center::roll_in_low(sigma, pl.base_sigma(s), smask);
    return ord.value_of_planes(sigma, pl.gamma(word));
  }
};

struct SigSharedBoth {  // 4b
  gm_center::Order ord; gm_center::TwoBitPlanes<true> pl;
  uint64_t sg = 0, gm = 0, sgr = 0, gmr = 0, mask; int p;
  SigSharedBoth(std::shared_ptr<const GreedyMiniTable> t, int p_, uint64_t tog)
      : ord(t, p_, tog), pl(p_, false), mask(pl.mask()), p(p_) {}
  inline uint64_t next(uint8_t s) {
    sg = gm_center::roll_in_low(sg, pl.base_sigma(s), mask);
    gm = gm_center::roll_in_low(gm, s >> 1, mask);
    sgr = (sgr >> 1) | (pl.base_sigma(s) << (p - 1));        // rc: reversed order, same Sigma
    gmr = (gmr >> 1) | ((uint64_t)((s >> 1) ^ 1) << (p - 1));  // complement flips Gamma
    return ord.value_of_both_strands(sg, gm, sgr, gmr);
  }
};

struct SigSharedRank {  // 4c
  gm_center::Order ord; gm_center::TwoBitPlanes<true> pl; uint64_t sigma = 0, smask;
  SigSharedRank(std::shared_ptr<const GreedyMiniTable> t, int p, uint64_t tog)
      : ord(t, p, tog), pl(p, false), smask(pl.mask()) {}
  inline uint64_t next(uint8_t s) {
    sigma = gm_center::roll_in_low(sigma, pl.base_sigma(s), smask);
    return ord.rank_of_sigma(sigma);
  }
};

struct SigLean {  // 5
  const uint16_t* r; uint32_t key = 0, mask;
  SigLean(const uint16_t* r_, int p) : r(r_), mask((1u << p) - 1) {}
  inline uint64_t next(uint8_t s) {
    key = ((key << 1) | ((s ^ (s >> 1)) & 1u)) & mask;
    return r[key];
  }
};

// ------------------------------------------------------------------ drivers
template <class S> static uint64_t run_insert(S& s, const uint8_t* x, size_t n) {
  uint64_t acc = 0;
  for (size_t i = 0; i < n; i++) acc += s.next(x[i]);
  return acc;
}

// Monotone-deque sliding minimum over W values; counts changes of the selected
// (leftmost-minimum) position; optionally records the selected position/window.
template <class S>
static uint64_t run_window(S& s, const uint8_t* x, size_t n, std::vector<uint32_t>* sel) {
  uint64_t val[32]; uint32_t pos[32];
  uint32_t head = 0, tail = 0; int64_t last = -1; uint64_t changes = 0;
  for (size_t i = 0; i < n; i++) {
    const uint64_t v = s.next(x[i]);
    if ((int)i < P - 1) continue;
    const uint32_t j = (uint32_t)(i - (P - 1));
    while (tail > head && val[(tail - 1) & 31] > v) tail--;
    val[tail & 31] = v; pos[tail & 31] = j; tail++;
    if ((int64_t)pos[head & 31] + W <= (int64_t)j) head++;
    if ((int)j >= W - 1) {
      const uint32_t p0 = pos[head & 31];
      if ((int64_t)p0 != last) { changes++; last = p0; }
      if (sel) sel->push_back(p0);
    }
  }
  return changes;
}

// ------------------------------------------------------------------- input
static std::vector<uint8_t> make_random(size_t n) {
  std::vector<uint8_t> x(n); uint64_t z = 88172645463325252ULL;
  for (size_t i = 0; i < n; i += 32) {
    z += 0x9E3779B97F4A7C15ULL; uint64_t t = z;
    t = (t ^ (t >> 30)) * 0xBF58476D1CE4E5B9ULL; t = (t ^ (t >> 27)) * 0x94D049BB133111EBULL; t ^= t >> 31;
    for (size_t k = 0; k < 32 && i + k < n; k++) x[i + k] = (t >> (2 * k)) & 3;
  }
  return x;
}
static std::vector<uint8_t> read_fastq(const std::string& path, size_t n) {
  std::ifstream f(path); if (!f) throw std::runtime_error("cannot open " + path);
  std::vector<uint8_t> x; x.reserve(n); std::string line; size_t ln = 0;
  while (x.size() < n && std::getline(f, line)) {
    if ((ln++ & 3) != 1) continue;
    for (char c : line) {
      int v = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
      if (v >= 0 && x.size() < n) x.push_back((uint8_t)v);
    }
  }
  return x;
}

static double med(std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; }

struct Row { std::string name; double mi, ni, mw, nw; uint64_t changes; };

int main(int argc, char** argv) {
  if (argc < 9) { fprintf(stderr, "usage: P TABLE.gm TABLE_dna.gm DNA_SEED L SRC NBASES REPS\n"); return 1; }
  P = atoi(argv[1]); const std::string gm = argv[2], dna = argv[3];
  const uint64_t seed = strtoull(argv[4], 0, 0); const int L = atoi(argv[5]); W = L - P + 1;
  const std::string src = argv[6]; const size_t n = (size_t)atof(argv[7]); const int reps = atoi(argv[8]);
  std::vector<uint8_t> x = src == "random" ? make_random(n) : read_fastq(src, n);
  fprintf(stderr, "p=%d W=%d src=%s bases=%zu reps=%d\n", P, W, src.c_str(), x.size(), reps);

  auto table = load_gm_table(gm);
  if (!table || table->gm_k() != P) { fprintf(stderr, "bad table\n"); return 1; }
  const uint64_t tog = 0x2a5;

  auto inst = [&](CCustomOrder::Mode m, CCustomOrder::LutBackend b, const std::string& file, uint64_t sd) {
    CCustomOrder::Config c; c.mode = m; c.lut_backend = b; c.file = file; c.signature_len = P; c.seed = sd;
    CCustomOrder::Install(c);
  };

  // ------------------------------ correctness (first CN bases)
  {
    const size_t cn = std::min<size_t>(x.size(), 20000000);
    std::vector<uint32_t> s2, s3, s3b, s4a, s4b, s4c, s5, s1;
    inst(CCustomOrder::Mode::DEFAULT, CCustomOrder::LutBackend::BinaryRuntime, "", 0);
    { SigKmc<CMmer::LUT_DEFAULT> g(P); run_window(g, x.data(), cn, &s1); }
    inst(CCustomOrder::Mode::TABLE, CCustomOrder::LutBackend::DnaExpanded, dna, 0);
    { SigKmc<CMmer::LUT_DNA> g(P); run_window(g, x.data(), cn, &s2); }
    inst(CCustomOrder::Mode::TABLE, CCustomOrder::LutBackend::BinaryRuntime, gm, seed);
    { SigKmc<CMmer::LUT_BINARY> g(P); run_window(g, x.data(), cn, &s3); }
    inst(CCustomOrder::Mode::TABLE, CCustomOrder::LutBackend::BinaryRuntime, gm, 0);
    { SigKmc<CMmer::LUT_BINARY> g(P); run_window(g, x.data(), cn, &s3b); }
    { SigShared g(table, P, tog); run_window(g, x.data(), cn, &s4a); }
    { SigSharedBoth g(table, P, tog); run_window(g, x.data(), cn, &s4b); }
    { SigSharedRank g(table, P, tog); run_window(g, x.data(), cn, &s4c); }
    { SigLean g(table->data(), P); run_window(g, x.data(), cn, &s5); }
    // independent rank per p-mer (inf -> 0xFFFF; kmc's binary LUT ranks inf as table size, still the max)
    const size_t np = cn - P + 1; std::vector<uint32_t> rk(np); uint32_t key = 0, km = (1u << P) - 1;
    for (size_t i = 0; i < cn; i++) {
      key = ((key << 1) | ((x[i] ^ (x[i] >> 1)) & 1u)) & km;
      if ((int)i >= P - 1) rk[i - (P - 1)] = table->data()[key];
    }
    const size_t nw = s2.size(); std::vector<uint32_t> mn(nw);
    { std::deque<size_t> q; size_t w = 0;
      for (size_t j = 0; j < np; j++) {
        while (!q.empty() && rk[q.back()] > rk[j]) q.pop_back();
        q.push_back(j); if (q.front() + W <= j) q.pop_front();
        if ((int)j >= W - 1) mn[w++] = rk[q.front()];
      } }
    auto rank_ok = [&](const char* nm, const std::vector<uint32_t>& s) {
      size_t bad = 0; for (size_t w = 0; w < nw; w++) if (rk[s[w]] != mn[w]) bad++;
      printf("  rank-of-selected != window-min rank: %-10s %zu / %zu\n", nm, bad, nw);
    };
    auto diff = [&](const char* a, const std::vector<uint32_t>& u, const char* b, const std::vector<uint32_t>& v) {
      size_t d = 0; for (size_t w = 0; w < nw; w++) d += u[w] != v[w];
      printf("  selected position differs: %s vs %s: %zu / %zu\n", a, b, d, nw);
    };
    auto cnt = [&](const std::vector<uint32_t>& s) { size_t c = 0; for (size_t w = 0; w < nw; w++) c += (w == 0 || s[w] != s[w - 1]); return c; };
    printf("CORRECTNESS on first %zu bases (%zu windows)\n", cn, nw);
    rank_ok("2 dna", s2); rank_ok("3 binary", s3); rank_ok("3 bin sd0", s3b); rank_ok("4a", s4a); rank_ok("4b", s4b); rank_ok("4c", s4c); rank_ok("5 lean", s5);
    diff("2 dna", s2, "3 binary(seed)", s3); diff("3 binary(seed)", s3, "4a", s4a); diff("4a", s4a, "4b", s4b);
    diff("3 binary(seed)", s3, "3 binary(seed0)", s3b);
    printf("  super-k-mer counts: 1 %zu | 2 %zu | 3 %zu | 3(seed0) %zu | 4a %zu | 4b %zu | 4c %zu | 5 %zu\n",
           cnt(s1), cnt(s2), cnt(s3), cnt(s3b), cnt(s4a), cnt(s4b), cnt(s4c), cnt(s5));
    fflush(stdout);
  }

  // ------------------------------ timing
  std::vector<Row> rows;
  auto bench = [&](const char* name, auto make) {
    std::vector<double> ti, tw; uint64_t ch = 0, sink = 0;
    for (int r = -1; r < reps; r++) {  // r=-1 warm-up
      { auto s = make(); auto t0 = clk::now(); sink += run_insert(s, x.data(), x.size());
        double d = std::chrono::duration<double>(clk::now() - t0).count() * 1e9 / x.size(); if (r >= 0) ti.push_back(d); }
      { auto s = make(); auto t0 = clk::now(); ch = run_window(s, x.data(), x.size(), nullptr);
        double d = std::chrono::duration<double>(clk::now() - t0).count() * 1e9 / x.size(); if (r >= 0) tw.push_back(d); }
    }
    rows.push_back({name, med(ti), *std::min_element(ti.begin(), ti.end()), med(tw), *std::min_element(tw.begin(), tw.end()), ch});
    fprintf(stderr, "  done %s (sink %llu)\n", name, (unsigned long long)(sink & 1));
  };
  inst(CCustomOrder::Mode::DEFAULT, CCustomOrder::LutBackend::BinaryRuntime, "", 0);
  bench("1 default norm[]", [&] { return SigKmc<CMmer::LUT_DEFAULT>(P); });
  inst(CCustomOrder::Mode::TABLE, CCustomOrder::LutBackend::DnaExpanded, dna, 0);
  bench("2 kmc DNA LUT", [&] { return SigKmc<CMmer::LUT_DNA>(P); });
  inst(CCustomOrder::Mode::TABLE, CCustomOrder::LutBackend::BinaryRuntime, gm, seed);
  bench("3 kmc binary LUT", [&] { return SigKmc<CMmer::LUT_BINARY>(P); });
  bench("4a shared value_of_planes", [&] { return SigShared(table, P, tog); });
  bench("4b shared value_of_both", [&] { return SigSharedBoth(table, P, tog); });
  bench("4c shared rank only (lb)", [&] { return SigSharedRank(table, P, tog); });
  bench("5 lean bin rank (lb)", [&] { return SigLean(table->data(), P); });

  printf("TIMING ns/base (median / min); ratio = median vs variant 1; window changes = super-k-mers per 1000 bases\n");
  printf("%-28s %8s %8s %6s | %8s %8s %6s | %8s\n", "variant", "ins med", "ins min", "x v1", "win med", "win min", "x v1", "skm/kb");
  for (auto& r : rows)
    printf("%-28s %8.3f %8.3f %6.2f | %8.3f %8.3f %6.2f | %8.2f\n", r.name.c_str(), r.mi, r.ni, r.mi / rows[0].mi,
           r.mw, r.nw, r.mw / rows[0].mw, 1000.0 * r.changes / x.size());
  return 0;
}
