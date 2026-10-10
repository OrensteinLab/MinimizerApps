// Reverse complement of a 2-bit-per-base packed DNA k-mer (public-domain
// bit-reversal trick, adapted from
// https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel).
//
// Verified bit-identical to kraken2's own kraken2::dna_reverse_complement
// (src/dna_kmer.h) for DNA_REVCOM_VERSION 1, kraken2's default -- kept as a
// separate header (rather than folded into greedymini_rank.h) so kraken2's
// existing consumers, which already have kraken2::dna_reverse_complement in
// scope, don't collide with this one; this header is for consumers (like
// kmc) that don't already have a copy.

#ifndef MINIMIZER_ORDERS_DNA_REVCOMP_H_
#define MINIMIZER_ORDERS_DNA_REVCOMP_H_

#include <cstdint>

enum { DNA_REVCOM_VERSION = 1 };

inline uint64_t dna_reverse_complement(uint64_t kmer, uint8_t n,
                                       int revcom_version = DNA_REVCOM_VERSION) {
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
  return ((~kmer) >> (sizeof(kmer) * 8 - n * 2)) & ((1ull << (n * 2)) - 1);
}

#endif  // MINIMIZER_ORDERS_DNA_REVCOMP_H_
