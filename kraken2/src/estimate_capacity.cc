/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include <getopt.h>
#include "kraken2_headers.h"
#include "kv_store.h"
#include "mmscanner.h"
#include "sketching_schemes.h"
#include <memory>
#include "seqreader.h"
#include "utilities.h"
#include <ios>
#include <stdio.h>
#include <unistd.h>
#include "kseq.h"

using std::string;
using std::cout;
using std::queue;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::unordered_set;
using std::vector;
using namespace kraken2;

#define RANGE_SECTIONS 1024  // must be power of 2
#define RANGE_MASK (RANGE_SECTIONS - 1)
#define MAX_N RANGE_SECTIONS
#define DEFAULT_N 4
#define DEFAULT_BLOCK_SIZE (30 * 1024 * 1024)  // yes, 30 MB

struct Options {
  size_t k, l, n;
  bool input_is_protein;
  int threads;
  size_t block_size;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  SketchingScheme sketching_scheme;
  Order order;
  std::string greedy_mini_filename;
  std::shared_ptr<const GreedyMiniTable> gm_table;
  GmProjection gm_projection = GmProjection::AT0_CG1;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
void ProcessSequence(string &seq, Options &opts,
    vector<unordered_set<uint64_t>> &sets);
void ProcessSequences(Options &opts);

int main(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  Options opts;
  opts.k = 0;
  opts.l = 0;
  opts.n = DEFAULT_N;
  opts.threads = 1;
  opts.input_is_protein = false;
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
  opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.block_size = DEFAULT_BLOCK_SIZE;
  opts.sketching_scheme = SketchingScheme::MINIMIZER;
  opts.order = Order::RANDOM_XOR;
  ParseCommandLine(argc, argv, opts);

  if (opts.order == Order::GM) {
    if (opts.greedy_mini_filename.empty())
      errx(EX_USAGE, "--greedymini-file is required when --order greedymini");
    try {
      opts.gm_table = MinimizerScanner::shared_gm_table(
          opts.greedy_mini_filename);
    } catch (const std::exception &ex) {
      errx(EX_DATAERR, "greedymini table load failed: %s", ex.what());
    }
    if (!opts.gm_table || opts.gm_table->empty())
      errx(EX_NOINPUT, "failed to load greedymini file: %s",
           opts.greedy_mini_filename.c_str());
  }
  omp_set_num_threads(opts.threads);
  ProcessSequences(opts);
  return 0;
}

void ProcessSequences(Options& opts)
{
  std::vector<vector<unordered_set<uint64_t>>> vector_of_sets(opts.threads);
  BatchSequenceReader reader;

#pragma omp parallel for
  for (auto i = 0; i < opts.threads; i++) {
    BatchSequenceReader reader_copy(reader);
    vector_of_sets[i].resize(opts.n);
    bool have_work = true;
    Sequence sequence;

    while (have_work) {
#pragma omp critical(batch_reading)
      {
        have_work = reader_copy.NextSequence(sequence);
      }
      if (have_work) {
        ProcessSequence(sequence.seq, opts, vector_of_sets[i]);
      }
    }
  }
  auto& sets = vector_of_sets[0];
  for (auto i = 1; i < vector_of_sets.size(); i++) {
    for (auto j = 0; j < opts.n; j++) {
      sets[j].insert(vector_of_sets[i][j].begin(), vector_of_sets[i][j].end());
    }
  }
  size_t sum_set_sizes = 0;
  for (auto& s : sets) {
    sum_set_sizes += s.size();
  }
  sum_set_sizes++; // ensure non-zero estimate
  cout << (size_t)(sum_set_sizes * RANGE_SECTIONS * 1.0 / opts.n) << endl;
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;
  long long sig;

  static const struct option long_options[] = {
    {"sketching", required_argument, nullptr, 's'},
    {"order", required_argument, nullptr, 'O'},
    {"greedymini-file", required_argument, nullptr, 'G'},
    {nullptr, 0, nullptr, 0}};
  while ((opt = getopt_long(argc, argv, "?hk:l:n:S:T:B:p:XG:s:O:P:",
                            long_options, nullptr)) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'p' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "must have at least 1 thread");
        opts.threads = sig;
        break;
      case 'B' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "block size must be positive");
        opts.block_size = sig;
        break;
      case 'n' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "n must be positive integer");
        if (sig > MAX_N)
          errx(EX_USAGE, "n must be no more than %d", MAX_N);
        opts.n = sig;
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
      case 'X' :
        opts.input_is_protein = true;
        break;
      case 'S' :
        // strtoull: 64-bit masks with MSB set overflow signed strtol → LONG_MAX
        opts.spaced_seed_mask = strtoull(optarg, nullptr, 2);
        break;
      case 'T' :
        opts.toggle_mask = strtoull(optarg, nullptr, 2);
        break;
      case 's':
        try {
          opts.sketching_scheme = parse_sketching_scheme(optarg);
        } catch (const std::exception &ex) {
          errx(EX_USAGE, "%s", ex.what());
        }
        break;
      case 'O':
        try {
          opts.order = parse_order(optarg);
        } catch (const std::exception &ex) {
          errx(EX_USAGE, "%s", ex.what());
        }
        break;
      case 'G':
        opts.greedy_mini_filename = optarg;
        break;
      case 'P':
        try {
          opts.gm_projection = parse_gm_projection(optarg);
        } catch (const std::exception &ex) {
          errx(EX_USAGE, "%s", ex.what());
        }
        break;
    }
  }
  if (opts.gm_projection != GmProjection::AT0_CG1 && opts.order != Order::GM)
    errx(EX_USAGE, "-P (--gm-projection) applies to --order greedymini only");

  if (opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK)
    ExpandSpacedSeedMask(opts.spaced_seed_mask,
      opts.input_is_protein ? BITS_PER_CHAR_PRO : BITS_PER_CHAR_DNA);
  if (opts.k == 0 || opts.l == 0) {
    cerr << "missing mandatory integer parameter" << endl;
    usage();
  }
  if (opts.k < opts.l) {
    cerr << "k cannot be less than l" << endl;
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: estimate_capacity <options>" << endl
       << endl
       << "Options (*mandatory):" << endl
       << "* -k INT        Set length of k-mers" << endl
       << "* -l INT        Set length of minimizers" << endl
       << "  -n INT        Set maximum qualifying hash code" << endl
       << "  -X            Input sequences are proteins" << endl
       << "  -S BITSTRING  Spaced seed mask" << endl
       << "  -T BITSTRING  Minimizer ordering toggle mask" << endl
       << "  -B INT        Read block size" << endl
       << "  -p INT        Number of threads" << endl
       << "  -s STRING     Sketching (--sketching): minimizer (default), mod-minimizer" << endl
       << "  -O STRING     Order (--order): xor (default), xor-pcg, greedymini" << endl
       << "  -G FILENAME   GreedyMini rank file (--greedymini-file, .gm);" << endl
       << "                required for --order greedymini" << endl
       << "  -P STRING     GreedyMini projection (--gm-projection): at0-cg1 (default)"
       << " or at1-cg0" << endl;
  exit(exit_code);
}

void ProcessSequence(string &seq, Options &opts,
    vector<unordered_set<uint64_t>> &sets)
{
  MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                           ! opts.input_is_protein, opts.toggle_mask,
                           CURRENT_REVCOM_VERSION,
                           opts.sketching_scheme, opts.order,
                           opts.greedy_mini_filename,
                           opts.gm_table,
                           opts.gm_projection == GmProjection::AT1_CG0);
  // Add terminator for protein sequences if not already there
  if (opts.input_is_protein && seq.back() != '*')
    seq.push_back('*');
  scanner.LoadSequence(seq);
  uint64_t *minimizer_ptr;
  while ((minimizer_ptr = scanner.NextMinimizer())) {
    if (scanner.is_ambiguous())
      continue;
    uint64_t hash_code = MurmurHash3(*minimizer_ptr);
    if ((hash_code & RANGE_MASK) < opts.n) {
      sets[hash_code & RANGE_MASK].insert(*minimizer_ptr);
    }
  }
}
