#pragma once

// Test-only reference for util::compute_minimizer: the minimizer definition
// evaluated literally over all loci (used by test_gm_rolling.cpp).

#include "include/util.hpp"

namespace sshash::util {

/* The definition evaluated literally, see compute_minimizer_all_loci. */
template <typename Kmer>
minimizer_info compute_minimizer_reference(Kmer const& kmer, Kmer const& kmer_rc, const uint64_t k,
                                           const uint64_t m, hasher_type const& hasher)  //
{
    return dispatch_order(hasher.mode(), [&](auto mode) {
        if constexpr (is_mod_mode(decltype(mode)::value)) {
            return compute_mod_minimizer<decltype(mode)::value, Kmer>(kmer, kmer_rc, k, m, hasher);
        } else if constexpr (decltype(mode)::value == custom_order_mode::GREEDYMINI) {
            return compute_minimizer_all_loci<Kmer>(kmer, kmer_rc, k, m, hasher);
        } else {
            return compute_minimizer_all_loci<Kmer>(kmer, kmer_rc, k, m,
                                                    hasher.plain<decltype(mode)::value>());
        }
    });
}

}  // namespace sshash::util
