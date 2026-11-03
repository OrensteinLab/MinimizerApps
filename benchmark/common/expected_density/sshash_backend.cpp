// Expected density of an SSHash scheme (see driver.h). Counts as SSHash's own
// num_minimizer_positions and num_kmers do: the loop below is
// dictionary_builder::compute_minimizer_tuples (sshash/src/builder/
// compute_minimizer_tuples.cpp), one super-k-mer per change of the selected
// position, over the strings as encode_strings packs them, with SSHash's own
// kmer_iterator, hasher and scheme_minimizer_iterator (canonical parsing).
//
// Options: --k K (the k-mer length, sshash -k) --m M (the minimizer length, -m)
// --scheme native|xor|xor-pcg-hash|greedymini|mod-minimizer-xor|mod-minimizer-greedymini
// [--gm FILE] [--projection at0-cg1|at1-cg0] [--t T]; --seeds are sshash -s.

#include "driver.h"

#include "include/util.hpp"
#include "include/kmer_iterator.hpp"
#include "include/minimizer_iterator.hpp"

namespace {

using namespace sshash;
typedef default_kmer_t Kmer;

struct SshashBackend {
  uint64_t k, m, mod_t;
  custom_scheme_mode scheme;
  gm_projection proj;
  std::string gm_path;
  hasher_type hasher;

  explicit SshashBackend(const expdens::Args &a)
      : k(a.need_i("k")), m(a.need_i("m")), mod_t(a.i("t", 0)),
        scheme(parse_custom_scheme_mode(a.need("scheme"))),
        proj(parse_gm_projection(a.str("projection", "at0-cg1"))), gm_path(a.str("gm")) {}

  void set_seed(uint64_t seed) {
    hasher = hasher_type(seed);
    hasher.configure(scheme, m, gm_path, proj, mod_t);
  }

  void count(const char *seq, size_t n, uint64_t &n_kmers, uint64_t &n_selected) {
    if (n < k) return;
    // encode_strings: one string, then the sentinel words kmer_iterator reads past the end
    bits::bit_vector::builder strings_builder;
    for (size_t i = 0; i < n; ++i) strings_builder.append_bits(Kmer::char_to_uint(seq[i]), Kmer::bits_per_char);
    for (int dummy_bits = Kmer::uint_kmer_bits; dummy_bits; dummy_bits -= 64) strings_builder.append_bits(0, 64);

    dispatch_scheme(hasher.mode(), [&](auto mode) {
      scheme_minimizer_iterator<Kmer, decltype(mode)::value> minimizer_it(k, m, hasher);
      kmer_iterator<Kmer, bits::bit_vector::builder> kmer_it(strings_builder, k);
      const uint64_t begin = 0;
      kmer_it.at(Kmer::bits_per_char * begin);
      minimizer_it.set_position(begin);
      uint64_t prev_pos = constants::invalid_uint64;
      for (uint64_t j = 0; j != n - k + 1; ++j) {
        minimizer_info mini_info = minimizer_it.next(kmer_it.get());
        if (prev_pos == constants::invalid_uint64) prev_pos = mini_info.pos_in_seq;
        if (mini_info.pos_in_seq != prev_pos) {
          n_selected++;  // save(prev_mini_info, ...)
          prev_pos = mini_info.pos_in_seq;
        }
        kmer_it.next();
      }
      n_selected++;  // the last save
    });
    n_kmers += n - k + 1;
  }
};

}  // namespace

int main(int argc, char **argv) { return expdens::run<SshashBackend>(argc, argv); }
