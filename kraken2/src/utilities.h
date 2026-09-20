/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_UTILITIES_H_
#define KRAKEN2_UTILITIES_H_

#include "kraken2_headers.h"

// Functions used by 2+ programs that I couldn't think of a better place for.

namespace kraken2 {

// Turns a simple bitstring of 0s and 1s into one where the 1s and 0s are expanded,
// e.g. 010110 expanded by a factor of 2 would become 001100111100
// Allows specification of the spaced seed to represent positions in the sequence
// rather than actual bits in the internal representation
inline void ExpandSpacedSeedMask(uint64_t &spaced_seed_mask,
                                 const int bit_expansion_factor) {
  uint64_t new_mask = 0;
  uint64_t bits = 1ULL << bit_expansion_factor;
  bits--;
  for (int i = 64 / bit_expansion_factor - 1; i >= 0; i--) {
    new_mask <<= bit_expansion_factor;
    if ((spaced_seed_mask >> i) & 1)
      new_mask |= bits;
  }
  spaced_seed_mask = new_mask;
}

// Kraken `--minimizer-spaces`: bit pattern 1*(ℓ−2s)+(01)*s, then expand by bpc.
// s==0 → 0 (no AND), matching DEFAULT_SPACED_SEED_MASK.
// Where the mask is applied (by design):
//   MINIMIZER+RandomXor/Pcg: DNA-lex then AND (rank input and hash key).
//   MINIMIZER+GM: AND each strand inside Rank; hash key is DNA-lex(winner) AND.
//   MOD: AND stored/returned canonical ℓ-mers only; t-mer ranks are unmasked.
inline uint64_t MakeSpacedSeedMask(int l, int spaces, int bpc) {
  if (spaces <= 0 || l <= 0 || bpc <= 0)
    return 0;
  uint64_t m = 0;
  const int kept = l - 2 * spaces;
  for (int i = 0; i < kept; i++)
    m = (m << 1) | 1ULL;
  for (int i = 0; i < spaces; i++)
    m = (m << 2) | 1ULL;
  ExpandSpacedSeedMask(m, bpc);
  return m;
}

std::vector<std::string> SplitString(const std::string &str,
  const std::string &delim = "\t", const size_t max_fields = (size_t) -1);

}

#endif
