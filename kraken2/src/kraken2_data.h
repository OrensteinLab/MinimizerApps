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

// db_version 2: opts.k2d holds the sketching scheme / order / GreedyMini
// metadata below, and Minimizer+GreedyMini with l > gm_k ranks by the centre
// extension (minimizer_rankers.h). Databases with any other db_version (the
// original Kraken 2 layout, or earlier fork versions) are not supported.
static const int KRAKEN2_DB_VERSION = 2;

struct IndexOptions {
  // ---- original Kraken 2 fields (do not reorder) ----
  size_t k;
  size_t l;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  bool dna_db;
  uint64_t minimum_acceptable_hash_value;
  int revcom_version;   // Fix bug from before K2.0.8
  int db_version;  // To allow for future database structural changes
  int db_type;  // To allow for future use of other data structures

  // ---- fork fields ----
  uint8_t  sketching_scheme;
  uint8_t  gm_alphabet;       // reserved byte, always written as 0
  uint16_t gm_k;
  uint32_t reserved;          // low 8 bits = Order
  uint64_t gm_num_entries;
};

inline void DecodeIndexSchemes(const IndexOptions &opts,
                               SketchingScheme &sketching, Order &order) {
  sketching = (SketchingScheme)opts.sketching_scheme;
  order = (Order)(opts.reserved & 0xffu);
}

inline void EncodeIndexSchemes(IndexOptions &opts, SketchingScheme sketching,
                               Order order) {
  opts.db_version = KRAKEN2_DB_VERSION;
  opts.sketching_scheme = (uint8_t)sketching;
  opts.reserved = (uint32_t)(uint8_t)order;
}

// Message for a database whose opts.k2d is not the current format.
inline const char *IndexRebuildMessage() {
  return "database was built by an older/different Kraken 2 build; rebuild it";
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
