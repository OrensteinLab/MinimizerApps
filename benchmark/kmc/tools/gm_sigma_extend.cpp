// Extend a binary GreedyMini .gm table (2^p ranks) to KMC's own DNA-alphabet
// _dna.gm format: a precomputed 4^p-entry sort-key + canonical-DNA lookup,
// with a fixed tie-break seed baked in (see kmc_api/custom_order.cpp's
// load_dna_gm_file/KMC_DNA_GM_MAGIC for the reader).
//
// This is a KMC-specific artifact, not the same as kraken2's DNA-alphabet
// .gm tables (see common/minimizer-schemes/greedymini_rank.h's header
// comment) -- but it's built from the same shared AT/CG-key and
// reverse-complement primitives kmc's own runtime uses, via
// common/minimizer-schemes/greedymini_rank.h and dna_revcomp.h, so this tool
// and kmc_api/mmer.h can never compute the DNA key / reverse complement
// differently from each other.
//
// (Extend-k here is kmc's prefix extension; kraken2 and sshash extend a
// binary table to longer l with the center extension of
// common/minimizer-schemes/greedymini_center.h instead.)

#include "../../../common/minimizer-schemes/greedymini_rank.h"
#include "../../../common/minimizer-schemes/dna_revcomp.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

constexpr uint64_t KMC_DNA_GM_MAGIC = 0x314D444E414D4B43ull;  // "KMCANDM1"

[[noreturn]] void die(const std::string &msg) {
  std::cerr << msg << "\n";
  std::exit(1);
}

// dna_to_binary_key(dna, p) == gm_extract_dna_key(dna, p, p) (l == gm_k, no
// Extend-k) -- verified bit-identical to kmc's own CMmer::DnaToBinaryKey.
uint32_t dna_to_binary_key(uint64_t dna, int p) {
  return gm_extract_dna_key(dna, p, p);
}

struct DnaTables {
  std::vector<uint32_t> sort_key;
  std::vector<uint32_t> sig_id;
};

DnaTables build_dna_tables(const std::vector<uint64_t> &ranks, int k,
                          uint32_t seed) {
  const uint32_t bin_k = (uint32_t)ranks.size();
  if (bin_k != (1u << k))
    die("rank table size " + std::to_string(bin_k) + " != 2^k=" +
        std::to_string(1u << k));
  const uint64_t dna_n = 1ull << (2 * k);

  std::vector<uint64_t> rank_lut(dna_n);
  for (uint64_t i = 0; i < dna_n; i++)
    rank_lut[i] = ranks[dna_to_binary_key(i, k)];

  std::vector<uint32_t> sig_id(dna_n);
  for (uint64_t i = 0; i < dna_n; i++) {
    const uint64_t rev = dna_reverse_complement(i, (uint8_t)k, DNA_REVCOM_VERSION);
    const uint64_t ri = rank_lut[i], rr = rank_lut[rev];
    if (ri != rr) {
      sig_id[i] = (uint32_t)(ri < rr ? i : rev);
      continue;
    }
    const uint64_t ti = i ^ seed, tr = rev ^ seed;
    sig_id[i] = (uint32_t)(ti != tr ? (ti < tr ? i : rev) : std::min(i, rev));
  }

  std::vector<uint32_t> order(dna_n);
  std::iota(order.begin(), order.end(), 0u);
  std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
    const uint64_t ra = rank_lut[a], rb = rank_lut[b];
    if (ra != rb)
      return ra < rb;
    return (a ^ seed) < (b ^ seed);
  });

  std::vector<uint32_t> packed(dna_n);
  for (uint32_t pos = 0; pos < dna_n; pos++)
    packed[order[pos]] = pos;

  DnaTables out;
  out.sort_key.resize(dna_n);
  out.sig_id = std::move(sig_id);
  for (uint64_t i = 0; i < dna_n; i++)
    out.sort_key[i] = packed[out.sig_id[i]];
  return out;
}

void write_dna_gm(const std::string &path, int k, uint64_t seed,
                  const DnaTables &t) {
  const uint64_t n = t.sort_key.size();
  std::ofstream out(path, std::ios::binary);
  if (!out)
    die("cannot write " + path);
  const uint64_t header[4] = {KMC_DNA_GM_MAGIC, (uint64_t)k, seed, n};
  out.write(reinterpret_cast<const char *>(header), sizeof(header));
  out.write(reinterpret_cast<const char *>(t.sort_key.data()),
            (std::streamsize)(n * sizeof(uint32_t)));
  out.write(reinterpret_cast<const char *>(t.sig_id.data()),
            (std::streamsize)(n * sizeof(uint32_t)));
  if (!out)
    die("error writing " + path);
}

std::string default_output(const std::string &input) {
  const std::string suffix = ".gm";
  if (input.size() > suffix.size() &&
      input.compare(input.size() - suffix.size(), suffix.size(), suffix) == 0)
    return input.substr(0, input.size() - suffix.size()) + "_dna.gm";
  return input + "_dna.gm";
}

}  // namespace

int main(int argc, char **argv) {
  std::string input, output;
  uint64_t seed = 0;

  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "-o" || arg == "--output") {
      if (++i >= argc)
        die("-o requires a path");
      output = argv[i];
    } else if (arg == "--seed") {
      if (++i >= argc)
        die("--seed requires a value");
      seed = std::strtoull(argv[i], nullptr, 0);
    } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
      die("unknown option: " + arg);
    } else {
      if (!input.empty())
        die("only one input file may be given");
      input = arg;
    }
  }
  if (input.empty())
    die("usage: gm_sigma_extend [--seed N] [-o OUTPUT] BINARY.gm");
  if (output.empty())
    output = default_output(input);

  const auto ranks = load_gm_file_raw(input);
  if (ranks.empty())
    die("cannot open or read binary .gm: " + input);
  const int k = gm_k_from_size(ranks.size());
  if ((1ull << k) != ranks.size())
    die("binary .gm length " + std::to_string(ranks.size()) +
        " is not a power of 2");

  const auto tables = build_dna_tables(ranks, k, (uint32_t)seed);
  write_dna_gm(output, k, seed, tables);

  std::cerr << "wrote " << output << " (k=" << k << ", seed=0x" << std::hex
            << seed << std::dec << ")\n";
  return 0;
}
