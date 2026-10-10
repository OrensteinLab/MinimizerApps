/*
 * Canonical DNA k-mer helpers shared by mmscanner and mod_minimizer.
 *
 * reverse_complement adapted from public domain code at:
 *   https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
 */

#ifndef KRAKEN2_DNA_KMER_H_
#define KRAKEN2_DNA_KMER_H_

#include <climits>
#include <cstdint>

namespace kraken2 {

inline uint64_t dna_reverse_complement(uint64_t kmer, uint8_t n,
                                       int revcom_version) {
  kmer = ((kmer & 0xCCCCCCCCCCCCCCCCUL) >> 2)
       | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer & 0xF0F0F0F0F0F0F0F0UL) >> 4)
       | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer & 0xFF00FF00FF00FF00UL) >> 8)
       | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer & 0xFFFF0000FFFF0000UL) >> 16)
       | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = (kmer >> 32) | (kmer << 32);
  if (revcom_version == 0)
    return (~kmer) & ((1ull << (n * 2)) - 1);
  return ((~kmer) >> (sizeof(kmer) * CHAR_BIT - n * 2)) & ((1ull << (n * 2)) - 1);
}

inline uint64_t dna_canonical_representation(uint64_t kmer, uint8_t n,
                                             int revcom_version) {
  const uint64_t revcom = dna_reverse_complement(kmer, n, revcom_version);
  return kmer < revcom ? kmer : revcom;
}

// Matches MinimizerScanner::CURRENT_REVCOM_VERSION (Kraken 2 default).
enum { DNA_REVCOM_VERSION = 1 };

}  // namespace kraken2

#endif
