// Expected density of a Minimap2 scheme (see driver.h). Counts as the particular
// density of benchmark/minimap2 does (lib/export.py attach_sketch_density): the
// selected k-mers are the entries the index stores, i.e. those mm_sketch
// (native scheme) or mm_sketch_scheme (any other) writes, as index.c calls them;
// the k-mers are those minimap2 sketches (lib/sketch_kmers.py): the bases after
// homopolymer compression under --hpc (map-pb), minus k-1.
//
// Options: --k K --w W [--hpc]
// --scheme native|xor|xor-pcg-hash|greedymini|mod-minimizer-xor|mod-minimizer-greedymini
// [--t T] [--gm FILE] [--projection at0-cg1|at1-cg0]; --seeds are the order seeds.

#include "driver.h"

extern "C" {
#include "minimap.h"
#include "mmpriv.h"
#include "scheme.h"
}

namespace {

struct Minimap2Backend {
  int k, w, is_hpc;
  mm_scheme_conf_t conf;
  std::string gm_path;
  mm_scheme_t *scheme = nullptr;
  mm128_v p = {0, 0, 0};

  explicit Minimap2Backend(const expdens::Args &a)
      : k(a.need_i("k")), w(a.need_i("w")), is_hpc(a.has("hpc") ? 1 : 0), gm_path(a.str("gm")) {
    std::memset(&conf, 0, sizeof conf);
    const std::string s = a.need("scheme");
    if (s == "native") conf.kind = MM_SCHEME_NATIVE;
    else if (mm_scheme_parse_kind(s.c_str(), &conf.kind) != 0)
      throw std::runtime_error("unknown minimap2 scheme " + s);
    conf.greedymini_file = gm_path.empty() ? nullptr : gm_path.c_str();
    conf.gm_projection = a.str("projection", "at0-cg1") == "at1-cg0" ? MM_GM_PROJ_AT1_CG0 : MM_GM_PROJ_AT0_CG1;
    conf.t = a.i("t", 0);
    conf.k = k;
    conf.w = w;
  }
  ~Minimap2Backend() {
    if (scheme) mm_scheme_destroy(scheme);
    std::free(p.a);
  }

  void set_seed(uint64_t seed) {
    if (conf.kind == MM_SCHEME_NATIVE) return;
    if (scheme) mm_scheme_destroy(scheme);
    conf.seed = seed;
    char err[512];
    scheme = mm_scheme_init(&conf, err, sizeof err);
    if (!scheme) throw std::runtime_error(err);
  }

  void count(const char *seq, size_t len, uint64_t &n_kmers, uint64_t &n_selected) {
    // random chunks have no ambiguous bases: one segment
    size_t h = len;
    if (is_hpc) {
      h = len ? 1 : 0;
      for (size_t i = 1; i < len; i++) h += seq[i] != seq[i - 1];
    }
    if (h >= (size_t)k) n_kmers += h - k + 1;
    p.n = 0;
    if (scheme)
      mm_sketch_scheme(scheme, 0, seq, (int)len, w, k, 0, is_hpc, &p);
    else
      mm_sketch(0, seq, (int)len, w, k, 0, is_hpc, &p);
    n_selected += p.n;
  }
};

}  // namespace

int main(int argc, char **argv) { return expdens::run<Minimap2Backend>(argc, argv); }
