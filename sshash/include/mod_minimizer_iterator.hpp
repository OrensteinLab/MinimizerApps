#pragma once

#include "kmer.hpp"
#include "util.hpp"
#include "../../common/minimizer-orders/mod_residue_deque.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace sshash {

/*
    Incremental mod-sampling: amortized O(1) per kmer, replacing
    util::compute_mod_minimizer_single_strand's O(k-t+1) brute-force scan
    with a monotone deque per t-mer-locus residue class (mod w) -- the same
    algorithm (Groot Koerkamp & Pibiri 2024, Algorithm 6; same single-strand,
    "prefer the smaller i mod w" semantics), ported from kraken2::ModSampler
    (kraken2/src/mod_minimizer.h), sharing its ModTmerData/ModResidueDeque
    primitive (common/minimizer-orders/mod_residue_deque.h).

    Unlike kraken2 (which streams individual bases and must remember old
    m-mer values in a ring buffer), this iterator receives the whole current
    kmer on every call, so the winning m-mer is always freshly re-derivable
    in O(1) from `kmer` once the winning residue p is known -- no ring
    buffer of stored candidates needed, just the rank bookkeeping.

    Self-checked against the brute-force reference on every call (like
    minimizer_iterator's own assert), so any divergence aborts immediately
    in a debug-assertions build.
*/
template <typename Kmer>
struct mod_minimizer_iterator {
    mod_minimizer_iterator() {}

    mod_minimizer_iterator(uint64_t k, uint64_t m, uint64_t t, hasher_type const& hasher,
                           uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_t(t)
        , m_w(k - m + 1)
        , m_hasher(hasher)  //
    {
        assert(t <= m and m <= k);
        const uint64_t num_t_loci = k - t + 1;
        const uint64_t per_residue = num_t_loci / m_w + 4;
        if (per_residue > uint64_t(ModResidueDeque::kMaxCap)) {
            throw std::invalid_argument("mod_minimizer_iterator: residue deque capacity exceeded");
        }
        m_deques.assign(m_w, ModResidueDeque());
        for (auto& dq : m_deques) dq.cap = int(per_residue);
        m_front_rank.assign(m_w, constants::invalid_uint64);
        set_position(position);
    }

    void set_position(uint64_t position) {
        m_position = position;
        m_abs_tmer_pos = position;
        for (auto& dq : m_deques) dq.clear();
        std::fill(m_front_rank.begin(), m_front_rank.end(), constants::invalid_uint64);
        m_primed = false;
    }

    minimizer_info next(Kmer kmer) {
        auto tau = [&](Kmer const& window) {
            Kmer tmer = window;
            tmer.take_chars(m_t);
            return uint64_t(tmer);
        };

        if (!m_primed) {
            /* First kmer since set_position: prime all w deques with a full
               brute-force pass (mirrors minimizer_iterator's own
               rescan-on-first-call trick). */
            Kmer window = kmer;
            for (uint64_t i = 0; i != m_k - m_t + 1; ++i) {
                uint64_t rank = m_hasher.hash(tau(window));
                uint64_t residue = m_abs_tmer_pos % m_w;
                m_deques[residue].push_monotone(ModTmerData{rank, (ssize_t)m_abs_tmer_pos});
                ++m_abs_tmer_pos;
                window.drop_char();
            }
            for (uint64_t r = 0; r != m_w; ++r) {
                m_front_rank[r] =
                    m_deques[r].empty() ? constants::invalid_uint64 : m_deques[r].front().rank;
            }
            m_primed = true;
        } else {
            m_position += 1;

            /* Push the one new rightmost t-mer locus. */
            Kmer window = kmer;
            window.drop_chars(m_k - m_t);
            uint64_t rank = m_hasher.hash(tau(window));
            uint64_t residue = m_abs_tmer_pos % m_w;
            m_deques[residue].push_monotone(ModTmerData{rank, (ssize_t)m_abs_tmer_pos});
            m_front_rank[residue] = m_deques[residue].front().rank;
            ++m_abs_tmer_pos;

            /* Evict the one residue class that just expired: the previous
               window's own start locus, m_position-1. */
            uint64_t expire_residue = (m_position - 1) % m_w;
            ModResidueDeque& dq = m_deques[expire_residue];
            while (!dq.empty() && dq.front().pos < (ssize_t)m_position) dq.pop_front();
            m_front_rank[expire_residue] =
                dq.empty() ? constants::invalid_uint64 : dq.front().rank;
        }

        /* Combine: smallest rank wins; ties -> smallest p, matching
           compute_mod_minimizer_single_strand's exact tie rule. Residue `r`
           (a fixed absolute-position class, used only for internal
           bookkeeping so each locus stays in the same deque across the
           slide) and output position `p = i mod w` (relative to the
           *current* kmer's own start) are related by a rotation that shifts
           by one step every call: p = (r - (m_position mod w)) mod w. */
        const uint64_t start_residue = m_position % m_w;
        uint64_t best_rank = constants::invalid_uint64;
        uint64_t best_p = 0;
        for (uint64_t r = 0; r != m_w; ++r) {
            if (m_front_rank[r] == constants::invalid_uint64) continue;
            uint64_t p = (r + m_w - start_residue) % m_w;
            if (m_front_rank[r] < best_rank or (m_front_rank[r] == best_rank and p < best_p)) {
                best_rank = m_front_rank[r];
                best_p = p;
            }
        }

        Kmer final_window = kmer;
        final_window.drop_chars(best_p);
        final_window.take_chars(m_m);
        minimizer_info mini_info(uint64_t(final_window), best_p);

        assert(mini_info ==
               util::compute_mod_minimizer_single_strand<Kmer>(kmer, m_k, m_m, m_t, m_hasher));
        return mini_info;
    }

private:
    uint64_t m_k, m_m, m_t, m_w;
    hasher_type m_hasher;
    uint64_t m_position;
    uint64_t m_abs_tmer_pos;
    bool m_primed = false;
    std::vector<ModResidueDeque> m_deques;
    std::vector<uint64_t> m_front_rank;
};

}  // namespace sshash
