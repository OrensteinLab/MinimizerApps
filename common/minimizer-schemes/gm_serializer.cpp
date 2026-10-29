// Serialize a ranked binary k-mer list to a classic GreedyMini .gm file.
// Shared by kraken2 and kmc: both call this one binary instead of keeping
// their own copies of this logic.
//
// Input (plain text): a comma-separated list of non-negative integer k-mers
// in rank order (rank 0 is the first entry), optionally wrapped in [ ].
//
//   --canonical
//       Each listed k-mer is the smaller of {x, reverse(x)}; both receive
//       that rank. Palindromes appear once.
//   (default: not canonical)
//       Each listed k-mer alone receives its rank; reverse is not implied.
//
// Output (.gm): [uint64_t length N=2^k][N x uint64_t ranks] (see
// greedymini_rank.h's load_gm_file_raw/write_gm_file_raw).
//   - ranks[key] = UHS rank, or sentinel N for non-UHS keys
//   - binary alphabet; DNA sigma-extension at rank time remaps ACGT->ATCG
//     (AT/CG key, residual lex tie-break), not baked into the table

#include "greedymini_rank.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

[[noreturn]] void die(const std::string &msg) {
  std::cerr << msg << "\n";
  std::exit(1);
}

std::vector<uint64_t> parse_kmer_list(const std::string &text) {
  std::vector<uint64_t> out;
  std::string tok;
  auto flush = [&]() {
    if (tok.empty())
      return;
    for (char c : tok)
      if (!isdigit((unsigned char)c))
        die("entry is not a non-negative integer: " + tok);
    out.push_back(std::stoull(tok));
    tok.clear();
  };
  for (char c : text) {
    if (c == '[' || c == ']' || c == ',' || isspace((unsigned char)c)) {
      flush();
    } else {
      tok.push_back(c);
    }
  }
  flush();
  if (out.empty())
    die("empty input");
  return out;
}

int infer_k(const std::vector<uint64_t> &kmers, int k) {
  uint64_t mx = 0;
  for (auto x : kmers)
    if (x > mx)
      mx = x;
  if (k > 0) {
    if (mx >= (1ULL << k))
      die("k=" + std::to_string(k) + " too small for max k-mer " +
          std::to_string(mx));
    return k;
  }
  int need = 0;
  while ((1ULL << need) <= mx)
    need++;
  return need < 1 ? 1 : need;
}

uint64_t revbits(uint64_t x, int k) {
  uint64_t r = 0;
  for (int i = 0; i < k; i++) {
    r = (r << 1) | (x & 1ULL);
    x >>= 1;
  }
  return r;
}

std::vector<uint64_t> build_ranks(const std::vector<uint64_t> &kmers, int k,
                                  bool canonical) {
  const uint64_t n = 1ULL << k;
  const uint64_t inf = n;
  std::vector<uint64_t> ranks(n, inf);
  for (size_t rank = 0; rank < kmers.size(); rank++) {
    const uint64_t x = kmers[rank];
    if (x >= n)
      die("k-mer " + std::to_string(x) + " out of range for k=" +
          std::to_string(k));
    const uint64_t keys[2] = {x, canonical ? revbits(x, k) : x};
    const int nkeys = canonical ? 2 : 1;
    for (int i = 0; i < nkeys; i++) {
      const uint64_t key = keys[i];
      const uint64_t prev = ranks[key];
      if (prev != inf && prev != rank)
        die("conflict at key " + std::to_string(key) + ": rank " +
            std::to_string(prev) + " vs " + std::to_string(rank));
      ranks[key] = rank;
    }
  }
  return ranks;
}

std::string read_file(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    die("cannot open input file: " + path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::string default_output(const std::string &input) {
  auto dot = input.rfind('.');
  auto slash = input.find_last_of("/\\");
  if (dot != std::string::npos && (slash == std::string::npos || dot > slash))
    return input.substr(0, dot) + ".gm";
  return input + ".gm";
}

}  // namespace

int main(int argc, char **argv) {
  std::string input, output;
  int k = 0;
  bool canonical = false;

  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "--canonical") {
      canonical = true;
    } else if (arg == "-o" || arg == "--output") {
      if (++i >= argc)
        die("-o requires a path");
      output = argv[i];
    } else if (arg == "-k") {
      if (++i >= argc)
        die("-k requires an integer");
      k = std::atoi(argv[i]);
    } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
      die("unknown option: " + arg);
    } else {
      if (!input.empty())
        die("only one input file may be given");
      input = arg;
    }
  }
  if (input.empty())
    die("usage: gm_serializer [--canonical] [-k N] [-o OUTPUT] INPUT.txt");
  if (output.empty())
    output = default_output(input);

  const auto kmers = parse_kmer_list(read_file(input));
  const int kk = infer_k(kmers, k);
  const auto ranks = build_ranks(kmers, kk, canonical);
  write_gm_file_raw(output, ranks);

  const uint64_t n = ranks.size();
  uint64_t uhs = 0;
  for (auto r : ranks)
    if (r < n)
      uhs++;
  std::cerr << "wrote " << output << "  k=" << kk << "  N=" << n
            << "  listed=" << kmers.size() << "  UHS=" << uhs
            << "  max_rank=" << (kmers.size() - 1)
            << "  canonical=" << (canonical ? "true" : "false") << "\n";
  return 0;
}
