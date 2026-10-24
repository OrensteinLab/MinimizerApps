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
// original Kraken 2 layout, or earlier fork versions), the wrong file size or
// a wrong magic are not supported.
static const int KRAKEN2_DB_VERSION = 2;
// "K2SK" little endian. Databases written by earlier fork builds carry the
// Order in the low byte of this field (0, 1 or 2), so they never match.
static const uint32_t KRAKEN2_OPTS_MAGIC = 0x4b32534bu;
// Tie rule of a MOD database's mod-minimizer (mod_minimizer_rule.h): the
// centre rule. MOD databases built with the G/T-strand rule store 0 and are
// rejected; MINIMIZER databases store 0 and are unaffected.
static const uint8_t kModTieRuleCentre = 1;

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
  uint8_t  sketching_scheme;  // SketchingScheme
  uint8_t  order;             // Order
  uint16_t gm_k;
  uint32_t magic;             // KRAKEN2_OPTS_MAGIC
  uint8_t  gm_projection;     // GmProjection (0 in databases built before it)
  uint8_t  mod_tie_rule;      // MOD: kModTieRuleCentre (0 = the retired G/T rule)
  uint8_t  reserved_[2];
  uint64_t gm_num_entries;
};
static_assert(sizeof(IndexOptions) == 80, "opts.k2d layout is 80 bytes");

// True iff `opts` (read from a file of `file_size` bytes) is in the current
// format: exact size, db_version and magic.
inline bool IndexOptionsValid(const IndexOptions &opts, size_t file_size) {
  return file_size == sizeof(IndexOptions) &&
         opts.db_version == KRAKEN2_DB_VERSION &&
         opts.magic == KRAKEN2_OPTS_MAGIC &&
         (opts.sketching_scheme == (uint8_t)SketchingScheme::MINIMIZER ||
          opts.sketching_scheme == (uint8_t)SketchingScheme::MOD) &&
         opts.order <= (uint8_t)Order::RANDOM_PCG &&
         opts.gm_projection <= (uint8_t)GmProjection::AT1_CG0 &&
         (opts.sketching_scheme != (uint8_t)SketchingScheme::MOD ||
          opts.mod_tie_rule == kModTieRuleCentre);
}

inline void DecodeIndexSchemes(const IndexOptions &opts,
                               SketchingScheme &sketching, Order &order) {
  sketching = (SketchingScheme)opts.sketching_scheme;
  order = (Order)opts.order;
}

inline void EncodeIndexSchemes(IndexOptions &opts, SketchingScheme sketching,
                               Order order) {
  opts.db_version = KRAKEN2_DB_VERSION;
  opts.sketching_scheme = (uint8_t)sketching;
  opts.order = (uint8_t)order;
  opts.magic = KRAKEN2_OPTS_MAGIC;
  opts.mod_tie_rule = sketching == SketchingScheme::MOD ? kModTieRuleCentre : 0;
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
