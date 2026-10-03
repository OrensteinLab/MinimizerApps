/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_KRAKEN2_DATA_H_
#define KRAKEN2_KRAKEN2_DATA_H_

#include "kraken2_headers.h"
#include "readcounts.h"
#include "sketching_schemes.h"

namespace kraken2 {

// db_version 0 — original Kraken 2 opts.k2d (Minimizer + Random only).
// db_version 1 — sketching_scheme, order, GreedyMini metadata in trailing fields.
// db_version 2 — same layout; binary GreedyMini with l > gm_k ranks by the
//                centre extension (minimizer_rankers.h). A version-1
//                Minimizer+GreedyMini database with l > gm_k was built with the
//                old prefix extension and selects different minimizers.
static const int KRAKEN2_DB_VERSION_ORIGINAL = 0;
static const int KRAKEN2_DB_VERSION_SCHEMES = 1;
static const int KRAKEN2_DB_VERSION_GM_CENTRE = 2;
static const int KRAKEN2_DB_VERSION = KRAKEN2_DB_VERSION_GM_CENTRE;

struct IndexOptions {
  // ---- db_version 0 fields (do not reorder) ----
  size_t k;
  size_t l;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  bool dna_db;
  uint64_t minimum_acceptable_hash_value;
  int revcom_version;   // Fix bug from before K2.0.8
  int db_version;  // To allow for future database structural changes
  int db_type;  // To allow for future use of other data structures

  // ---- db_version >= 1 ----
  uint8_t  sketching_scheme;
  uint8_t  gm_alphabet;
  uint16_t gm_k;
  uint32_t reserved;          // low 8 bits = Order
  uint64_t gm_num_entries;
};

inline void DecodeIndexSchemes(const IndexOptions &opts,
                               SketchingScheme &sketching, Order &order) {
  if (opts.db_version >= KRAKEN2_DB_VERSION_SCHEMES) {
    sketching = (SketchingScheme)opts.sketching_scheme;
    order = (Order)(opts.reserved & 0xffu);
  } else {
    sketching = SketchingScheme::MINIMIZER;
    order = Order::RANDOM_XOR;
  }
}

inline void EncodeIndexSchemes(IndexOptions &opts, SketchingScheme sketching,
                               Order order) {
  opts.db_version = KRAKEN2_DB_VERSION;
  opts.sketching_scheme = (uint8_t)sketching;
  opts.reserved = (uint32_t)(uint8_t)order;
}

// True for a database whose minimizers were chosen with the retired prefix
// Extend-k order (Minimizer + binary GreedyMini, nucleotide, l > gm_k, built
// before db_version 2). Such a database must be rebuilt.
inline bool IndexUsesLegacyGmPrefix(const IndexOptions &opts) {
  if (opts.db_version >= KRAKEN2_DB_VERSION_GM_CENTRE ||
      opts.db_version < KRAKEN2_DB_VERSION_SCHEMES)
    return false;
  SketchingScheme sketching;
  Order order;
  DecodeIndexSchemes(opts, sketching, order);
  return order == Order::GM && sketching == SketchingScheme::MINIMIZER &&
         opts.dna_db && opts.gm_alphabet == 0 && opts.gm_k > 0 &&
         opts.l > opts.gm_k;
}

typedef uint64_t taxid_t;
const taxid_t TAXID_MAX = (taxid_t) -1;

typedef std::unordered_map<taxid_t, uint64_t> taxon_counts_t;

#ifdef EXACT_COUNTING
  // Exact unique k-mer sets for taxon read counters.
  typedef ReadCounts<unordered_set<uint64_t>> READCOUNTER;
#else
  typedef ReadCounts<HyperLogLogPlusMinus<uint64_t>> READCOUNTER;
#endif
  typedef std::map<taxid_t, READCOUNTER> taxon_counters_t;
}

#endif
