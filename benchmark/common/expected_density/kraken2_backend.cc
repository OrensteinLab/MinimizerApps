// Expected density of a Kraken2 scheme (see driver.h). Counts as the particular
// density of benchmark_minimizers.cc (accumulate_scan_stats) does: the k-mers
// are the non-ambiguous NextMinimizer() steps and the selected k-mers are the
// steps whose minimizer differs from the previous step's, with Kraken2's own
// MinimizerScanner (mmscanner.cc, built as for build_db and classify).
//
// Options: --l L (minimizer length) --k K (k-mer length) --spaces S
// --scheme xor|xor-pcg-hash|greedymini|mod-minimizer-xor|mod-minimizer-greedymini
// [--gm FILE] [--projection at0-cg1|at1-cg0]; --seeds are the toggle masks.

#include "driver.h"

#include "mmscanner.h"
#include "utilities.h"

namespace {
using namespace kraken2;

struct Kraken2Backend {
  int k, l;
  uint64_t spaced_seed_mask;
  SketchingScheme sketching;
  Order order;
  std::string gm_path;
  std::shared_ptr<const GreedyMiniTable> gm_tbl;
  bool gm_flip;
  std::unique_ptr<MinimizerScanner> sc;

  explicit Kraken2Backend(const expdens::Args &a)
      : k(a.need_i("k")), l(a.need_i("l")),
        spaced_seed_mask(MakeSpacedSeedMask(l, a.i("spaces", 0), BITS_PER_CHAR_DNA)),
        gm_path(a.str("gm")), gm_flip(a.str("projection", "at0-cg1") == "at1-cg0") {
    const std::string s = a.need("scheme");
    const bool mod = s.compare(0, 14, "mod-minimizer-") == 0;
    const std::string o = mod ? s.substr(14) : s;
    sketching = mod ? SketchingScheme::MOD : SketchingScheme::MINIMIZER;
    if (o == "xor") order = Order::RANDOM_XOR;
    else if (o == "xor-pcg-hash" && !mod) order = Order::RANDOM_PCG;
    else if (o == "greedymini") order = Order::GM;
    else throw std::runtime_error("unknown kraken2 scheme " + s);
    if (order == Order::GM) {
      gm_tbl = MinimizerScanner::shared_gm_table(gm_path);
      if (!gm_tbl || gm_tbl->empty()) throw std::runtime_error("cannot load --gm " + gm_path);
    }
  }

  void set_seed(uint64_t toggle) {
    sc.reset(new MinimizerScanner(k, l, spaced_seed_mask, true, toggle, CURRENT_REVCOM_VERSION,
                                  sketching, order, gm_path, gm_tbl, gm_flip));
  }

  void count(const char *seq, size_t len, uint64_t &n_kmers, uint64_t &n_selected) {
    if (len < 2) return;
    const std::string dna(seq, len);
    sc->LoadSequence(dna);
    uint64_t last = 0;
    bool have = false;
    while (uint64_t *mm = sc->NextMinimizer()) {
      if (sc->is_ambiguous()) {
        have = false;
        continue;
      }
      n_kmers++;
      if (!have || *mm != last) {
        n_selected++;
        last = *mm;
        have = true;
      }
    }
  }
};

}  // namespace

int main(int argc, char **argv) { return expdens::run<Kraken2Backend>(argc, argv); }
