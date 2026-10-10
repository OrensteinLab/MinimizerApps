// scan_dump: test tool that exercises ONLY the public MinimizerScanner API
// (constructor, LoadSequence, NextMinimizer, is_ambiguous) so that the same
// source compiles against the fork's mmscanner.h and against upstream
// Kraken 2's (define SCAN_DUMP_UPSTREAM; only the default xor scheme exists).
//
//   scan_dump --reads FILE|random:N:LEN [--scheme S] [--mode M] [--out FILE]
//             [--toggle HEX] [--gm-dir DIR] [--reps N] [--chunk C]
//
//   S (labels as in the harness): xor (default), xor-pcg-hash, greedymini,
//     mod-minimizer-xor, mod-minimizer-greedymini.
//   M:  checksum  print FNV-1a checksum of the key stream (default)
//       dump      write the key stream to --out: one uint64 per NextMinimizer
//                 return (kAmbig for ambiguous returns), kSeqMark after each
//                 sequence; prints the checksum
//       chunk     scan each sequence in overlapping chunks like build_db
//                 (--chunk bases + k-1 overlap, LoadSequence(seq,start,fin));
//                 writes the sorted unique non-ambiguous key set of every
//                 sequence (then kSeqMark) and reports whether it equals the
//                 set of the whole-sequence scan
//       time      time the scan loop over --reps repetitions, ns/base
//                 (pin with taskset; discard rep 0 if the cache is cold)
//
// Reads: FASTQ (first char '@'), FASTA ('>') or random:N:LEN.
// Parameters are those of classify for a default DB: k=35 l=31, 7 spaces.
#include "mmscanner.h"
#include "utilities.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace kraken2;

static const uint64_t kAmbig = 0xFFFFFFFFFFFFFFFEULL;   // ambiguous return
static const uint64_t kSeqMark = 0xFFFFFFFFFFFFFFFFULL;  // end of sequence

// Same construction as MakeSpacedSeedMask(l, spaces, 2) in the fork; written
// out here so the file also compiles against upstream (ExpandSpacedSeedMask
// exists in both).
static uint64_t spaced_mask(int l, int spaces) {
  if (spaces <= 0)
    return 0;
  uint64_t m = 0;
  for (int i = 0; i < l - 2 * spaces; i++)
    m = (m << 1) | 1ULL;
  for (int i = 0; i < spaces; i++)
    m = (m << 2) | 1ULL;
  ExpandSpacedSeedMask(m, 2);
  return m;
}

static void usage() {
  fprintf(stderr,
          "usage: scan_dump --reads FILE|random:N:LEN [--scheme xor|xor-pcg-hash|"
          "greedymini|mod-minimizer-xor|mod-minimizer-greedymini]\n"
          "       [--gm-projection at0-cg1|at1-cg0]\n"
          "       [--mode checksum|dump|chunk|time] [--out FILE] [--toggle HEX]\n"
          "       [--gm-dir DIR] [--reps N] [--chunk C] [--k 35] [--l 31] [--spaces 7]\n");
  exit(2);
}

static std::vector<std::string> read_sequences(const std::string &src) {
  std::vector<std::string> reads;
  if (src.rfind("random:", 0) == 0) {
    int n, len;
    if (sscanf(src.c_str(), "random:%d:%d", &n, &len) != 2) usage();
    std::mt19937_64 rng(1);
    for (int i = 0; i < n; i++) {
      std::string s(len, 'A');
      for (auto &c : s) c = "ACGT"[rng() & 3];
      reads.push_back(s);
    }
    return reads;
  }
  std::ifstream f(src);
  if (!f) { fprintf(stderr, "cannot open %s\n", src.c_str()); exit(1); }
  std::string line;
  if (!std::getline(f, line)) return reads;
  if (!line.empty() && line[0] == '@') {
    std::string seq, plus, qual;
    while (std::getline(f, seq) && std::getline(f, plus) && std::getline(f, qual)) {
      reads.push_back(seq);
      if (!std::getline(f, line)) break;
    }
  } else {
    std::string cur;
    bool have = false;
    do {
      if (!line.empty() && line[0] == '>') {
        if (have) reads.push_back(cur);
        cur.clear();
        have = true;
      } else {
        cur += line;
      }
    } while (std::getline(f, line));
    if (have) reads.push_back(cur);
  }
  return reads;
}

static inline void fnv(uint64_t &h, uint64_t v) {
  h = (h ^ v) * 1099511628211ULL;
}

int main(int argc, char **argv) {
  std::string scheme = "xor", reads_src, mode = "checksum", out, toggle_s,
      gm_projection = "at0-cg1";
  std::string gm_dir = "orders/GreedyMini/";
  int reps = 3, k = 35, l = 31, spaces = 7;
  size_t chunk = 4096;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto next = [&]() -> std::string { if (i + 1 >= argc) usage(); return argv[++i]; };
    if (a == "--scheme") scheme = next();
    else if (a == "--reads") reads_src = next();
    else if (a == "--mode") mode = next();
    else if (a == "--out") out = next();
    else if (a == "--toggle") toggle_s = next();
    else if (a == "--gm-dir") gm_dir = next();
    else if (a == "--gm-projection") gm_projection = next();
    else if (a == "--reps") reps = atoi(next().c_str());
    else if (a == "--chunk") chunk = strtoull(next().c_str(), 0, 10);
    else if (a == "--k") k = atoi(next().c_str());
    else if (a == "--l") l = atoi(next().c_str());
    else if (a == "--spaces") spaces = atoi(next().c_str());
    else usage();
  }
  if (reads_src.empty()) usage();
  if (!gm_dir.empty() && gm_dir.back() != '/') gm_dir += '/';
  uint64_t toggle = toggle_s.empty() ? DEFAULT_TOGGLE_MASK
                                     : strtoull(toggle_s.c_str(), 0, 16);
  const uint64_t sp = spaced_mask(l, spaces);

#ifdef SCAN_DUMP_UPSTREAM
  if (scheme != "xor") {
    fprintf(stderr, "upstream build supports only --scheme xor\n");
    return 2;
  }
  MinimizerScanner sc(k, l, sp, true, toggle, CURRENT_REVCOM_VERSION);
#else
  SketchingScheme sk = SketchingScheme::MINIMIZER;
  Order order = Order::RANDOM_XOR;
  std::string gm_file;
  bool gm_flip = (parse_gm_projection(gm_projection) == GmProjection::AT1_CG0);
  if (scheme == "xor") {
  } else if (scheme == "xor-pcg-hash") {
    order = Order::RANDOM_PCG;
  } else if (scheme == "greedymini") {
    order = Order::GM; gm_file = gm_dir + "w5_k15_canonical.gm";
  } else if (scheme == "mod-minimizer-xor") {
    sk = SketchingScheme::MOD;
  } else if (scheme == "mod-minimizer-greedymini") {
    sk = SketchingScheme::MOD; order = Order::GM;
    gm_file = gm_dir + "w5_k6_canonical.gm";
  } else {
    usage();
  }
  std::shared_ptr<const GreedyMiniTable> tab;
  if (!gm_file.empty()) {
    tab = MinimizerScanner::shared_gm_table(gm_file);
    if (!tab) { fprintf(stderr, "cannot load %s\n", gm_file.c_str()); return 1; }
  }
  MinimizerScanner sc(k, l, sp, true, toggle, CURRENT_REVCOM_VERSION, sk, order,
                      "", tab, gm_flip);
#endif

  std::vector<std::string> reads = read_sequences(reads_src);
  uint64_t bases = 0;
  for (auto &r : reads) bases += r.size();

  std::ofstream dump;
  if (!out.empty()) dump.open(out, std::ios::binary);
  auto put = [&](uint64_t v) {
    if (dump.is_open()) dump.write((const char *)&v, sizeof v);
  };

  if (mode == "checksum" || mode == "dump") {
    if (mode == "dump" && !dump.is_open()) usage();
    uint64_t h = 14695981039346656037ULL, calls = 0, amb = 0;
    for (auto &r : reads) {
      sc.LoadSequence(r);
      uint64_t *m;
      while ((m = sc.NextMinimizer())) {
        calls++;
        if (sc.is_ambiguous()) { amb++; fnv(h, kAmbig); put(kAmbig); continue; }
        fnv(h, *m); put(*m);
      }
      fnv(h, kSeqMark); put(kSeqMark);
    }
    printf("%s seqs=%zu bases=%lu calls=%lu ambiguous=%lu checksum=%016lx\n",
           scheme.c_str(), reads.size(), (unsigned long)bases,
           (unsigned long)calls, (unsigned long)amb, (unsigned long)h);
  } else if (mode == "chunk") {
    if (!dump.is_open()) usage();
    uint64_t h = 14695981039346656037ULL, mismatches = 0, total = 0;
    std::vector<uint64_t> full, part;
    for (auto &r : reads) {
      full.clear(); part.clear();
      uint64_t *m;
      sc.LoadSequence(r);
      while ((m = sc.NextMinimizer()))
        if (!sc.is_ambiguous()) full.push_back(*m);
      for (size_t s = 0; s < r.size(); s += chunk) {
        size_t fin = std::min(r.size(), s + chunk + (size_t)k - 1);
        sc.LoadSequence(r, s, fin);
        while ((m = sc.NextMinimizer()))
          if (!sc.is_ambiguous()) part.push_back(*m);
      }
      std::sort(full.begin(), full.end());
      full.erase(std::unique(full.begin(), full.end()), full.end());
      std::sort(part.begin(), part.end());
      part.erase(std::unique(part.begin(), part.end()), part.end());
      if (full != part) mismatches++;
      total += part.size();
      for (uint64_t v : part) { fnv(h, v); put(v); }
      fnv(h, kSeqMark); put(kSeqMark);
    }
    printf("%s chunk=%zu seqs=%zu unique_keys=%lu seqs_differing_from_full_scan=%lu checksum=%016lx\n",
           scheme.c_str(), chunk, reads.size(), (unsigned long)total,
           (unsigned long)mismatches, (unsigned long)h);
  } else if (mode == "time") {
    for (int rep = 0; rep < reps; rep++) {
      uint64_t chk = 0, calls = 0;
      auto t0 = std::chrono::steady_clock::now();
      for (auto &r : reads) {
        sc.LoadSequence(r);
        uint64_t *m;
        while ((m = sc.NextMinimizer())) {
          calls++;
          chk = chk * 1099511628211ULL + (sc.is_ambiguous() ? 7 : *m);
        }
      }
      auto t1 = std::chrono::steady_clock::now();
      double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
      printf("%s rep=%d bases=%lu calls=%lu ns/base=%.4f chk=%lx\n",
             scheme.c_str(), rep, (unsigned long)bases, (unsigned long)calls,
             ns / bases, (unsigned long)chk);
    }
  } else {
    usage();
  }
  return 0;
}
