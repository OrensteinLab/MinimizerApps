/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "build_db.h"
#include "kraken2_headers.h"
#include "taxonomy.h"
#include "mmscanner.h"
#include "compact_hash.h"
#include "kraken2_data.h"
#include "sketching_schemes.h"
#include "utilities.h"

using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::vector;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using namespace kraken2;

// These will likely be overridden by wrapper script,
// just keeping sane defaults in case they aren't
#define DEFAULT_BLOCK_SIZE (10 * 1024 * 1024)  // 10 MB
#define DEFAULT_SUBBLOCK_SIZE (1024)

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
void ReadIDToTaxonMap(map<string, taxid_t> &id_map, string &filename);
void GenerateTaxonomy(Options &opts, map<string, taxid_t> &id_map);

struct TaxonSeqPair {
  taxid_t taxon;
  string seq;
};

int main(int argc, char **argv) {
  Options opts;
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
  opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.input_is_protein = false;
  opts.num_threads = 1;
  opts.block_size = DEFAULT_BLOCK_SIZE;
  opts.subblock_size = DEFAULT_SUBBLOCK_SIZE;
  opts.cht_cell_size = 32;
  opts.requested_bits_for_taxid = 0;
  opts.min_clear_hash_value = 0;
  opts.maximum_capacity = 0;
  opts.deterministic_build = true;
  opts.sketching_scheme = SketchingScheme::MINIMIZER;
  opts.order = Order::RANDOM_XOR;
  ParseCommandLine(argc, argv, opts);

  // Load GreedyMini rank table once (cached); embed as gm.k2d for classify.
  int gm_k = 0;
  if (opts.order == Order::GM) {
    if (opts.greedy_mini_filename.empty())
      errx(EX_USAGE, "must supply -G <file> when using GreedyMini order");
    try {
      opts.gm_table = MinimizerScanner::shared_gm_table(
          opts.greedy_mini_filename);
    } catch (const std::exception &ex) {
      errx(EX_DATAERR, "GreedyMini table load failed: %s", ex.what());
    }
    if (!opts.gm_table || opts.gm_table->empty())
      errx(EX_NOINPUT, "failed to load GreedyMini file: %s",
           opts.greedy_mini_filename.c_str());
    gm_k = opts.gm_table->gm_k();
    if (opts.sketching_scheme == SketchingScheme::MOD) {
      const int t = ModSampler::kDefaultR +
                    (((int)opts.l - ModSampler::kDefaultR) %
                     ((int)opts.k - (int)opts.l + 1));
      if (gm_k != t)
        errx(EX_DATAERR, "MOD+GM table gm_k=%d must equal t=%d", gm_k, t);
    } else if (gm_k > (int)opts.l) {
      errx(EX_DATAERR, "GreedyMini k_gm=%d exceeds minimizer l=%d", gm_k,
           (int)opts.l);
    } else if (!opts.input_is_protein && gm_k < (int)opts.l &&
               (((int)opts.l - gm_k) & 1) != 0) {
      errx(EX_DATAERR,
           "GreedyMini centre extension needs l - k_gm even (l=%d, k_gm=%d); "
           "valid l: %d, %d, ... (max 31)",
           (int)opts.l, gm_k, gm_k, gm_k + 2);
    }
    std::cerr << "Loaded GreedyMini rank table: k_gm=" << gm_k
              << ", UHS=" << opts.gm_table->uhs_size()
              << "/" << opts.gm_table->size()
              << "  sketching=" << sketching_scheme_name(opts.sketching_scheme)
              << " order=" << order_name(opts.order) << std::endl;
  }

  omp_set_num_threads( opts.num_threads );

  map<string, taxid_t> ID_to_taxon_map;

  ReadIDToTaxonMap(ID_to_taxon_map, opts.ID_to_taxon_map_filename);
  GenerateTaxonomy(opts, ID_to_taxon_map);

  std::cerr << "Taxonomy parsed and converted." << std::endl;

  Taxonomy taxonomy(opts.taxonomy_filename.c_str(), false);
  taxonomy.GenerateExternalToInternalIDMap();
  size_t bits_needed_for_value = 1;
  while ((1 << bits_needed_for_value) < (ssize_t) taxonomy.node_count())
    bits_needed_for_value++;
  if (opts.requested_bits_for_taxid > 0 &&
      bits_needed_for_value > opts.requested_bits_for_taxid)
    errx(EX_DATAERR, "more bits required for storing taxid");

  size_t bits_for_taxid = bits_needed_for_value;
  if (bits_for_taxid < opts.requested_bits_for_taxid)
    bits_for_taxid = opts.requested_bits_for_taxid;

  auto actual_capacity = opts.capacity;
  if (opts.maximum_capacity) {
    double frac = opts.maximum_capacity * 1.0 / opts.capacity;
    if (frac > 1)
      errx(EX_DATAERR, "maximum capacity larger than requested capacity");
    opts.min_clear_hash_value = (uint64_t) ((1 - frac) * UINT64_MAX);
    actual_capacity = opts.maximum_capacity;
  }

  if (opts.cht_cell_size == 32) {
    build<CompactHashCell>(taxonomy, ID_to_taxon_map, opts, actual_capacity, bits_for_taxid);
  } else if (opts.cht_cell_size == 40) {
    build<CompactHashCell40>(taxonomy, ID_to_taxon_map, opts, actual_capacity, bits_for_taxid);
  } else {
    errx(EX_DATAERR, "Unsupported CHT cell size");
  }

  IndexOptions index_opts = {0};
  index_opts.k = opts.k;
  index_opts.l = opts.l;
  index_opts.spaced_seed_mask = opts.spaced_seed_mask;
  index_opts.toggle_mask = opts.toggle_mask;
  index_opts.dna_db = ! opts.input_is_protein;
  index_opts.minimum_acceptable_hash_value = opts.min_clear_hash_value;
  index_opts.revcom_version = CURRENT_REVCOM_VERSION;
  index_opts.db_type = 0;
  EncodeIndexSchemes(index_opts, opts.sketching_scheme, opts.order);
  index_opts.gm_k = (uint16_t)gm_k;
  index_opts.gm_num_entries = opts.gm_table ? (uint64_t)opts.gm_table->size() : 0;
  ofstream opts_fs(opts.options_filename);
  opts_fs.write((char *) &index_opts, sizeof(index_opts));
  if (! opts_fs.good())
    errx(EX_OSERR, "Unable to write options file %s",
         opts.options_filename.c_str());
  opts_fs.close();

  if (opts.order == Order::GM && opts.gm_table)
    MinimizerScanner::embed_gm_table(opts.greedy_mini_filename,
                                     opts.hashtable_filename,
                                     opts.gm_table);

  std::cerr << " complete." << std::endl;

  return 0;
}

// This function exists to deal with NCBI's use of \x01 characters to denote
// the start of a new FASTA header in the same line (for non-redundant DBs).
// We return all sequence IDs in a header line, not just the first.
vector<string> ExtractNCBISequenceIDs(const string &header) {
  vector<string> list;
  string current_str;
  bool in_id = true;

  for (size_t i = 0; i < header.size(); ++i) {
    if (header[i] == 0x01) {
      // 0x01 starts new ID at next char
      if (! current_str.empty())
        list.push_back(current_str);
      current_str.clear();
      in_id = true;
    }
    else if (in_id && isspace(header[i])) {
      // spaces end ID
      if (! current_str.empty())
        list.push_back(current_str);
      current_str.clear();
      in_id = false;
    }
    else if (in_id) {
      // Build ID string char by char
      current_str.push_back(header[i]);
    }
  }
  if (! current_str.empty())
    list.push_back(current_str);
  return list;
}

void ReadIDToTaxonMap(map<string, taxid_t> &id_map, string &filename) {
  ifstream map_file(filename.c_str());
  if (! map_file.good())
    err(EX_NOINPUT, "unable to read from '%s'", filename.c_str());
  string line;

  while (getline(map_file, line)) {
    string seq_id;
    taxid_t taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;
    if (taxid)
      id_map[seq_id] = taxid;
  }
}

void GenerateTaxonomy(Options &opts, map<string, taxid_t> &id_map) {
  NCBITaxonomy ncbi_taxonomy(
    opts.ncbi_taxonomy_directory + "/nodes.dmp",
    opts.ncbi_taxonomy_directory + "/names.dmp"
  );
  for (auto &kv_pair : id_map) {
    if (kv_pair.second != 0) {
      ncbi_taxonomy.MarkNode(kv_pair.second);
    }
  }
  ncbi_taxonomy.ConvertToKrakenTaxonomy(opts.taxonomy_filename.c_str());
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;
  long long sig;

  while ((opt = getopt(argc, argv, "?hB:b:c:FH:m:n:o:t:k:l:M:p:r:S:T:XC:G:s:O:")) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'B' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "must have positive block size");
        opts.block_size = sig;
        break;
      case 'b' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "must have positive subblock size");
        opts.subblock_size = sig;
        break;
      case 'r' :
        sig = atoll(optarg);
        if (sig < 0)
          errx(EX_USAGE, "can't have negative bit storage");
        if (sig > 31)
          errx(EX_USAGE, "can't have more than 31 bits of storage for taxid");
        opts.requested_bits_for_taxid = sig;
        break;
      case 'p' :
        opts.num_threads = atoll(optarg);
        if (opts.num_threads < 1)
          errx(EX_USAGE, "can't have negative number of threads");
        if (opts.num_threads > omp_get_max_threads())
          errx(EX_USAGE, "OMP only wants you to use %d threads", omp_get_max_threads());
        break;
      case 'H' :
        opts.hashtable_filename = optarg;
        break;
      case 'm' :
        opts.ID_to_taxon_map_filename = optarg;
        break;
      case 'n' :
        opts.ncbi_taxonomy_directory = optarg;
        break;
      case 'o' :
        opts.options_filename = optarg;
        break;
      case 't' :
        opts.taxonomy_filename = optarg;
        break;
      case 'S' :
        // strtoull: 64-bit masks with MSB set overflow signed strtol → LONG_MAX
        opts.spaced_seed_mask = strtoull(optarg, nullptr, 2);
        break;
      case 'T' :
        opts.toggle_mask = strtoull(optarg, nullptr, 2);
        break;
      case 'k' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "k must be positive integer");
        opts.k = sig;
        break;
      case 'l' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "l must be positive integer");
        if (sig > 31)
          errx(EX_USAGE, "l must be no more than 31");
        opts.l = sig;
        break;
      case 'c' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "capacity must be positive integer");
        opts.capacity = sig;
        break;
      case 'M' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "max capacity must be positive integer");
        opts.maximum_capacity = sig;
        break;
      case 'F' :
        opts.deterministic_build = false;
        break;
      case 'X' :
        opts.input_is_protein = true;
        break;
      case 'C':
        sig = atoll(optarg);
        if (sig != 32 && sig != 40)
          errx(EX_USAGE, "CHT cell size should be either 32 or 40 bits");
        opts.cht_cell_size = sig;
        break;
      case 's': {
        try {
          opts.sketching_scheme = parse_sketching_scheme(optarg);
        } catch (const std::exception &ex) {
          errx(EX_USAGE, "%s", ex.what());
        }
        break;
      }
      case 'O': {
        try {
          opts.order = parse_order(optarg);
        } catch (const std::exception &ex) {
          errx(EX_USAGE, "%s", ex.what());
        }
        break;
      }
      case 'G':
        opts.greedy_mini_filename = optarg;
        break;
    }
  }

  if (opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK)
    ExpandSpacedSeedMask(opts.spaced_seed_mask,
      opts.input_is_protein ? BITS_PER_CHAR_PRO : BITS_PER_CHAR_DNA);
  if (opts.hashtable_filename.empty() ||
      opts.ID_to_taxon_map_filename.empty() ||
      opts.ncbi_taxonomy_directory.empty() ||
      opts.options_filename.empty() ||
      opts.taxonomy_filename.empty())
  {
    cerr << "missing mandatory filename parameter" << endl;
    usage();
  }
  if (opts.k == 0 || opts.l == 0 || opts.capacity == 0) {
    cerr << "missing mandatory integer parameter" << endl;
    usage();
  }
  if (opts.k < opts.l) {
    cerr << "k cannot be less than l" << endl;
    usage();
  }
  if (opts.block_size < opts.subblock_size) {
    cerr << "block size cannot be less than subblock size\n";
    usage();
  }
  if (opts.maximum_capacity > opts.capacity) {
    cerr << "maximum capacity option shouldn't specify larger capacity than normal" << endl;
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: build_db <options>\n"
       << "\n"
       << "Options (*mandatory):\n"
       << "* -H FILENAME   Kraken 2 hash table filename\n"
       << "* -m FILENAME   Sequence ID to taxon map filename\n"
       << "* -t FILENAME   Kraken 2 taxonomy filename\n"
       << "* -n DIR        NCBI taxonomy directory name\n"
       << "* -o FILENAME   Kraken 2 options filename\n"
       << "* -k INT        Set length of k-mers\n"
       << "* -l INT        Set length of minimizers\n"
       << "* -c INT        Set capacity of hash table\n"
       << "  -M INT        Set maximum capacity of hash table (MiniKraken)\n"
       << "  -S BITSTRING  Spaced seed mask\n"
       << "  -T BITSTRING  Minimizer toggle mask\n"
       << "  -X            Input seqs. are proteins\n"
       << "  -p INT        Number of threads\n"
       << "  -F            Use fast, nondeterministic building method\n"
       << "  -B INT        Read block size\n"
       << "  -b INT        Read subblock size\n"
       << "  -r INT        Bit storage requested for taxid\n"
       << "  -C INT        CHT cell size (default: 32)\n"
       << "  -s STRING     Sketching: Minimizer (default), Mod\n"
       << "  -O STRING     Order: RandomXor (default), GreedyMini, RandomPcg\n"
       << "  -G FILENAME   GreedyMini rank file (.gm); required for GreedyMini order\n";
  exit(exit_code);
}
