// PCG RXS-M-XS 64-bit output permutation (O'Neill, "PCG: A Family of Simple
// Fast Space-Efficient Statistically Good Algorithms for Random Number
// Generation", 2014). Used by both kraken2 and kmc as the "random PCG"
// minimizer/signature order: rank(x) = pcg_rxs_m_xs64(x ^ seed).
//
// This header is the single source of the scramble itself; each app keeps
// its own thin wrapper around it (kraken2's PcgRanker in
// src/minimizer_rankers.h, kmc's install_random()/build_gm_rank_by_dna() in
// kmc_api/custom_order.cpp) since the call-site shapes differ.

#ifndef MINIMIZER_ORDERS_PCG_RXS_M_XS_H_
#define MINIMIZER_ORDERS_PCG_RXS_M_XS_H_

#include <cstdint>

inline uint64_t pcg_rxs_m_xs64(uint64_t x) {
  const unsigned count = (unsigned)(x >> 59);
  x ^= x >> (5 + count);
  x *= 12605985483714917081ull;
  return x ^ (x >> 43);
}

#endif  // MINIMIZER_ORDERS_PCG_RXS_M_XS_H_
