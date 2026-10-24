#pragma once

#include "kmer.hpp"
#include "util.hpp"

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

namespace sshash {

/*
    NOTE: this is the iterator for the orders whose hash expects the numerically-canonical
    m-mer (native, xor, xor-pcg-hash), instantiated per order: `Mode` is a compile-time
    parameter so that the hash is inlined (for native, upstream's multiply and xor with the
    magic precomputed) rather than selected at every call. A strand-invariant order (GreedyMini)
    uses `strand_invariant_minimizer_iterator` below, selected by the callers: it is
    deliberately a separate class and not a branch in `next`, because a bigger `next` stops
    being inlined into the streaming lookup and slowed native queries by ~16%.

    "Re-scan" method: computes the minimizer of each kmer of a sequence,
    sliding one character at a time.

    See `util::compute_minimizer` for the definition of the minimizer; this
    iterator computes exactly the same thing incrementally, which the assertion
    at the end of `next` checks.

    The only extra state compared to a plain forward minimizer is `m_num_mins`,
    the number of loci of the current window attaining the minimum hash: a tie
    cannot be broken by position alone without breaking mirror-equivariance, so
    when there is one we have to look at the kmer's own orientation.

    The kmer is reverse-complemented once per call and that single value serves
    every locus, via rc(x_i) = rc(x)_{k-m-i}, so no m-mer is ever
    reverse-complemented on its own.
*/
template <typename Kmer, custom_order_mode Mode>
struct minimizer_iterator {
    static_assert(Mode != custom_order_mode::GREEDYMINI);

    minimizer_iterator() {}

    minimizer_iterator(uint64_t k, uint64_t m, hasher_type const& hasher, uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_min_value(constants::invalid_uint64)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher(hasher.plain<Mode>())
#ifndef NDEBUG
        , m_hasher_p(&hasher)
#endif
        //
    {
        assert(k > 0 and m <= k);
        set_position(position);
    }

    void set_position(uint64_t position) {
        m_position = position;
        reset();
    }

    void reset() {
        m_min_pos_in_kmer = 0;
        m_min_position = m_position - 1;
        m_num_mins = 0;
    }

    __attribute__((always_inline)) inline minimizer_info next(Kmer kmer) {
        Kmer kmer_rc = kmer;
        kmer_rc.reverse_complement_inplace(m_k);

        if (m_min_pos_in_kmer == 0) {
            /* min leaves the window: re-scan to compute the new min */
            m_position = m_min_position + 1;
            rescan(kmer, kmer_rc);
        } else {
            m_position += 1;
            Kmer window = kmer;
            window.drop_chars(m_k - m_m);
            uint64_t value = util::canonical_mmer_at<Kmer>(window, kmer_rc, m_k, m_m, m_k - m_m);
            uint64_t hash = m_hasher.hash(value);
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_min_value = value;
                m_min_position = m_position;
                m_min_pos_in_kmer = m_k - m_m;
                m_num_mins = 1;
            } else {
                /* Only the leftmost minimum can leave the window without a
                   re-scan, so the count stays valid across the slide. */
                assert(m_min_pos_in_kmer > 0);
                m_min_pos_in_kmer -= 1;
                if (hash == m_min_hash) m_num_mins += 1;
            }
        }

        minimizer_info mini_info(m_min_value, m_min_position, m_min_pos_in_kmer);
        if (m_num_mins > 1) break_tie(kmer, kmer_rc, mini_info);

        assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, kmer_rc, m_k, m_m, hasher()));

        return mini_info;
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_pos_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    uint64_t m_num_mins;
    /* Just the word the strand-dependent orders hash with (see plain_hasher), by value: a copy
       of the full hasher (GreedyMini: rank table + derived constants) would bloat -- and shift
       the hot members of -- the streaming query object. */
    plain_hasher<Mode> m_hasher;
#ifndef NDEBUG
    hasher_type const* m_hasher_p = nullptr;  // for the assertions only; must outlive the iterator
    hasher_type const& hasher() const { return *m_hasher_p; }
#endif

    void rescan(Kmer kmer, Kmer const& kmer_rc) {
        const uint64_t begin = m_position;

        /* first locus, peeled off the loop: see `util::compute_minimizer` */
        uint64_t min_value = util::canonical_mmer_at<Kmer>(kmer, kmer_rc, m_k, m_m, 0);
        kmer.drop_char();
        uint64_t min_hash = m_hasher.hash(min_value);
        uint64_t min_pos_in_kmer = 0;
        uint64_t num_mins = 1;

        for (uint64_t i = 1; i != m_k - m_m + 1; ++i) {
            uint64_t value = util::canonical_mmer_at<Kmer>(kmer, kmer_rc, m_k, m_m, i);
            kmer.drop_char();
            uint64_t hash = m_hasher.hash(value);
            /* Branch-free selection (the same as `if (hash < min_hash) {...} else if (hash ==
               min_hash) ++num_mins`): whether a locus beats the minimum so far is as
               unpredictable as the hashes, and with a branch here the compiler's choice (upstream's
               is branch-free) made streaming queries 2-3% slower on some data sets. */
            const uint64_t less = -uint64_t(hash < min_hash);  // all ones: a new (leftmost) minimum
            num_mins = (less & 1) | (~less & (num_mins + (hash == min_hash)));
            min_value = (less & value) | (~less & min_value);
            min_pos_in_kmer = (less & i) | (~less & min_pos_in_kmer);
            min_hash = (less & hash) | (~less & min_hash);
        }

        m_min_value = min_value;
        m_min_hash = min_hash;
        m_min_pos_in_kmer = min_pos_in_kmer;
        m_num_mins = num_mins;
        m_position += m_k - m_m;
        m_min_position = begin + m_min_pos_in_kmer;
    }

    /*
        Two or more loci of the window attain the minimum hash (all carrying the
        same class, so only the sampled position is at stake). Apply the
        centre-closest tie-break of `util::compute_minimizer`: take the tied
        locus closest to the window centre, and between two equally close ones
        -- mirror images i and k-m-i -- the smaller index if kmer <= rc(kmer),
        the larger otherwise. Forward and mirror-equivariant. It runs on ~1e-4
        of the windows, so the rescan below costs nothing overall.

        Tied loci cannot precede the leftmost minimum, so the scan starts there.
    */
    void break_tie(Kmer kmer, Kmer const& kmer_rc, minimizer_info& mini_info) const {
        assert(m_num_mins > 1);

        const uint64_t window_begin = m_min_position - m_min_pos_in_kmer;
        const uint64_t two_c = m_k - m_m;
        uint64_t best_dist = constants::invalid_uint64;
        uint64_t lo = m_min_pos_in_kmer;
        uint64_t hi = m_min_pos_in_kmer;

        Kmer window = kmer;
        window.drop_chars(m_min_pos_in_kmer);
        for (uint64_t i = m_min_pos_in_kmer; i != m_k - m_m + 1; ++i) {
            const uint64_t v = util::canonical_mmer_at<Kmer>(window, kmer_rc, m_k, m_m, i);
            if (m_hasher.hash(v) == m_min_hash) {
                assert(v == m_min_value);
                const uint64_t dist = 2 * i > two_c ? 2 * i - two_c : two_c - 2 * i;
                if (dist < best_dist) {
                    best_dist = dist;
                    lo = i;
                    hi = i;
                } else if (dist == best_dist) {
                    hi = i;
                }
            }
            window.drop_char();
        }

        const uint64_t chosen = (lo == hi or !(kmer_rc < kmer)) ? lo : hi;
        mini_info = minimizer_info(m_min_value, window_begin + chosen, chosen);
    }
};

/*
    Same selection as `minimizer_iterator`, for a strand-invariant order (GreedyMini:
    hash(x) == hash(rc(x))): hash the RAW forward m-mer of each window instead of the
    numerically-canonical one, and build the canonical value (the stored minimizer
    identity) only when a locus actually becomes the minimum. The reverse complement of the
    kmer is likewise materialized on demand (new minimum, re-scan result, tie, assertions)
    -- most calls need none of it. Hash ties are still value ties (hash() is injective on
    canonical m-mers), so the counting of ties is unchanged: the result is bit-identical to
    `minimizer_iterator`'s (the assertion in `next` checks it against util::compute_minimizer).

    Ranks first (GreedyMini): a hash is (rank, tie-break)
    with the rank in its top bits, so two windows with different ranks are ordered by their
    ranks alone. The iterator therefore keeps only the Sigma plane of the kmer, slid by one
    base per call, and the rank of every window (one table load, custom_order_hasher::
    centre_rank), in a ring (m_ring: 64 slots, each rank stored twice so that the k-m+1 window
    ranks of the current kmer are contiguous, window i at slot (m_base + i) mod 64; a re-scan
    takes their minimum 8 ranks at a time). The full hash of a window (pext + bit reversal, hash_gm) is
    computed only when ranks tie: a new window whose rank equals the minimum rank, or the
    windows sharing the minimum rank in a re-scan; the hash of the current minimum is
    computed lazily (m_min_hash_valid). Most windows (~75% on real data) carry the shared
    "infinity" rank of non-UHS keys and cost one load and one compare. If the kmer passed to
    `next` is the first one after reset() / set_position() (a new read or string), the plane
    and the ring are recomputed from the kmer; every other call must pass the one-base slide
    of the previous kmer (asserted in debug builds, where the rolled plane and the ring are
    also compared with fresh ones on every call).
*/
template <typename Kmer>
struct strand_invariant_minimizer_iterator {
    strand_invariant_minimizer_iterator() {}

    strand_invariant_minimizer_iterator(uint64_t k, uint64_t m, hasher_type const& hasher,
                                        uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_min_value(constants::invalid_uint64)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher_p(&hasher)  //
    {
        assert(k > 0 and m <= k);
        assert(k <= 63 and m <= 31);
        set_position(position);
    }

    void set_position(uint64_t position) {
        m_position = position;
        reset();
    }

    void reset() {
        m_min_pos_in_kmer = 0;
        m_min_position = m_position - 1;
        m_num_mins = 0;
        m_sigma_valid = false;
    }

    /* Not inlined: keeps the callers' hot paths (shared with the native orders) small. */
    __attribute__((noinline)) minimizer_info next(Kmer kmer) {
        assert(hasher().strand_invariant());
        update_ranks(kmer);
        Kmer kmer_rc = kmer;
        bool have_rc = false;
        auto rc = [&]() -> Kmer const& {
            if (!have_rc) {
                kmer_rc.reverse_complement_inplace(m_k);
                have_rc = true;
            }
            return kmer_rc;
        };

        if (m_min_pos_in_kmer == 0) {
            /* min leaves the window: re-scan to compute the new min */
            m_position = m_min_position + 1;
            rescan(kmer, rc);
        } else {
            m_position += 1;
            const uint64_t rank = m_ring[last_slot()];
            bool new_min = rank < m_min_rank;
            if (rank == m_min_rank) {
                /* equal ranks: the tie-break decides, as the full hashes would */
                const uint64_t min_hash = current_min_hash(kmer, m_min_pos_in_kmer - 1);
                const uint64_t hash = hash_of(kmer, m_k - m_m);
                new_min = hash < min_hash;
                if (new_min) {
                    m_min_hash = hash;
                    m_min_hash_valid = true;
                } else if (hash == min_hash) {
                    m_num_mins += 1;
                }
            }
            if (new_min) {
                if (rank < m_min_rank) m_min_hash_valid = false;
                set_last_window_as_min(kmer, rank, rc);
            } else {
                assert(m_min_pos_in_kmer > 0);
                m_min_pos_in_kmer -= 1;
            }
        }

        minimizer_info mini_info(m_min_value, m_min_position, m_min_pos_in_kmer);
        if (m_num_mins > 1) break_tie(kmer, rc(), mini_info);

        assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, rc(), m_k, m_m, hasher()));

        return mini_info;
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_pos_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    uint64_t m_num_mins;
    hasher_type const* m_hasher_p = nullptr;  // the hasher must outlive the iterator
    hasher_type const& hasher() const { return *m_hasher_p; }

    /* rank mode (see the class comment) */
    uint64_t m_sigma = 0;     // Sigma plane of the kmer last passed to `next`
#ifndef NDEBUG
    Kmer m_prev_kmer;  // debug-only check of the slide contract
#endif
    uint64_t m_min_rank = 0;  // rank of the current minimum; m_min_hash is valid iff m_min_hash_valid
    uint64_t m_base = 0;
    uint16_t m_ring[136];  // 2 x 64 slots (every rank stored twice) + padding for 8-lane loads
    bool m_sigma_valid = false;
    bool m_min_hash_valid = true;

    static uint64_t sigma_of(hasher_type const& h, Kmer const& kmer, uint64_t k) {
        uint64_t lo = uint64_t(kmer), hi = 0;
        if constexpr (sizeof(kmer.bits) > 8) hi = uint64_t(kmer.bits >> 64);
        return h.sigma_of_kmer(lo, hi, int(k));
    }

    /* window i of the current kmer is at slots (m_base + i) mod 64 and +64, so the w window
       ranks are contiguous from m_ring + m_base */
    uint64_t slot_of(uint64_t i) const { return (m_base + i) & 63; }
    uint64_t last_slot() const { return (m_base + (m_k - m_m)) & 63; }
    void set_rank(uint64_t slot, uint16_t r) {
        m_ring[slot] = r;
        m_ring[slot + 64] = r;
    }

    /* Contract: `kmer` is the one-base slide of the previous kmer passed to `next`, unless
       reset()/set_position() was called since (then the plane and the ring are recomputed).
       Checked in debug builds. */
    void update_ranks(Kmer const& kmer) {
        const uint64_t w = m_k - m_m + 1;
        if (m_sigma_valid) {
#ifndef NDEBUG
            Kmer a = m_prev_kmer;
            a.drop_char();
            Kmer c = kmer;
            c.take_chars(m_k - 1);
            assert(a == c);
#endif
            m_sigma = hasher().roll_sigma(m_sigma, kmer.at(m_k - 1), int(m_k));
            m_base = (m_base + 1) & 63;
            set_rank(last_slot(), hasher().centre_rank(m_sigma, int(w - 1)));
        } else {
            m_sigma = sigma_of(hasher(), kmer, m_k);
            m_base = 0;
            for (uint64_t i = 0; i != w; ++i) set_rank(i, hasher().centre_rank(m_sigma, int(i)));
        }
        m_sigma_valid = true;
#ifndef NDEBUG
        m_prev_kmer = kmer;
        const uint64_t kmask = (uint64_t(1) << m_k) - 1;
        assert((m_sigma & kmask) == sigma_of(hasher(), kmer, m_k));
        for (uint64_t i = 0; i != w; ++i) {
            assert(m_ring[slot_of(i)] == hasher().centre_rank(m_sigma, int(i)));
            assert(m_ring[slot_of(i) + 64] == m_ring[slot_of(i)]);
        }
#endif
    }

    /* full hash of the m-mer at positions i .. i+m-1 of the kmer */
    inline uint64_t hash_of(Kmer kmer, uint64_t i) const {
        kmer.drop_chars(i);
        kmer.take_chars(m_m);
        return hasher().hash_gm(uint64_t(kmer));
    }

    /* hash of the current minimum, computed on demand (the minimum is the window `pos`) */
    inline uint64_t current_min_hash(Kmer const& kmer, uint64_t pos) {
        if (!m_min_hash_valid) {
            m_min_hash = hash_of(kmer, pos);
            m_min_hash_valid = true;
        }
        return m_min_hash;
    }

    template <typename GetRc>
    inline void set_last_window_as_min(Kmer const& kmer, uint64_t rank, GetRc& rc) {
        Kmer last = kmer;
        last.drop_chars(m_k - m_m);
        m_min_value = util::canonical_mmer_at<Kmer>(last, rc(), m_k, m_m, m_k - m_m);
        m_min_rank = rank;
        m_min_position = m_position;
        m_min_pos_in_kmer = m_k - m_m;
        m_num_mins = 1;
    }

    template <typename GetRc>
    void finish_rescan(Kmer const& kmer, uint64_t begin, GetRc& rc) {
        Kmer window = kmer;
        window.drop_chars(m_min_pos_in_kmer);
        m_min_value = util::canonical_mmer_at<Kmer>(window, rc(), m_k, m_m, m_min_pos_in_kmer);
        m_min_position = begin + m_min_pos_in_kmer;
        m_position = begin + m_k - m_m;
    }

    /* Re-scan: the minimum rank over the window ranks; the windows attaining it are ranked by
       their full hashes only if there are several. */
    template <typename GetRc>
    void rescan(Kmer const& kmer, GetRc& rc) {
        const uint64_t begin = m_position;
        const uint64_t w = m_k - m_m + 1;
        const uint16_t* r = m_ring + m_base;  // the w window ranks, contiguous
        uint64_t min_rank = 0xFFFF, first = 0, count = 0;
#ifdef __SSE4_1__
        /* minimum rank, its leftmost window and the number of windows attaining it, 8 ranks at
           a time (lanes past the last window are masked) */
        const __m128i lane = _mm_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7);
        for (uint64_t c = 0; c < w; c += 8) {
            __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(r + c));
            if (w - c < 8) {
                v = _mm_or_si128(v, _mm_cmpgt_epi16(lane, _mm_set1_epi16(short(w - c - 1))));
            }
            const __m128i mp = _mm_minpos_epu16(v);
            const uint64_t mn = uint16_t(_mm_extract_epi16(mp, 0));
            if (mn < min_rank) {
                min_rank = mn;
                first = c + uint16_t(_mm_extract_epi16(mp, 1));
            }
        }
        const __m128i vm = _mm_set1_epi16(short(min_rank));
        for (uint64_t c = 0; c < w; c += 8) {
            const __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(r + c));
            uint32_t mask = uint32_t(_mm_movemask_epi8(_mm_cmpeq_epi16(v, vm)));
            if (w - c < 8) mask &= (1u << (2 * (w - c))) - 1;
            count += uint64_t(__builtin_popcount(mask)) / 2;
        }
#else
        min_rank = r[0];
        count = 1;
        for (uint64_t i = 1; i != w; ++i) {
            if (r[i] < min_rank) {
                min_rank = r[i];
                first = i;
                count = 1;
            } else if (r[i] == min_rank) {
                count += 1;
            }
        }
#endif
        m_min_rank = min_rank;
        m_min_pos_in_kmer = first;
        m_num_mins = 1;
        if (count == 1) {
            m_min_hash_valid = false;
        } else {
            /* several windows share the minimum rank: full hashes decide, leftmost wins */
            m_min_hash = hash_of(kmer, first);
            m_min_hash_valid = true;
            for (uint64_t i = first + 1; i != w; ++i) {
                if (r[i] != min_rank) continue;
                const uint64_t hash = hash_of(kmer, i);
                if (hash < m_min_hash) {
                    m_min_hash = hash;
                    m_min_pos_in_kmer = i;
                    m_num_mins = 1;
                } else if (hash == m_min_hash) {
                    m_num_mins += 1;
                }
            }
        }
        finish_rescan(kmer, begin, rc);
    }

    /* Centre-closest tie-break, see `minimizer_iterator::break_tie`. */
    void break_tie(Kmer kmer, Kmer const& kmer_rc, minimizer_info& mini_info) const {
        assert(m_num_mins > 1);

        const uint64_t window_begin = m_min_position - m_min_pos_in_kmer;
        const uint64_t two_c = m_k - m_m;
        uint64_t best_dist = constants::invalid_uint64;
        uint64_t lo = m_min_pos_in_kmer;
        uint64_t hi = m_min_pos_in_kmer;

        Kmer window = kmer;
        window.drop_chars(m_min_pos_in_kmer);
        for (uint64_t i = m_min_pos_in_kmer; i != m_k - m_m + 1; ++i) {
            Kmer mmer = window;
            mmer.take_chars(m_m);
            if (hasher().hash_gm(uint64_t(mmer)) == m_min_hash) {
                assert(util::canonical_mmer_at<Kmer>(window, kmer_rc, m_k, m_m, i) == m_min_value);
                const uint64_t dist = 2 * i > two_c ? 2 * i - two_c : two_c - 2 * i;
                if (dist < best_dist) {
                    best_dist = dist;
                    lo = i;
                    hi = i;
                } else if (dist == best_dist) {
                    hi = i;
                }
            }
            window.drop_char();
        }

        const uint64_t chosen = (lo == hi or !(kmer_rc < kmer)) ? lo : hi;
        mini_info = minimizer_info(m_min_value, window_begin + chosen, chosen);
    }
};

/*
    The strand-symmetric mod-minimizer (util::compute_mod_minimizer), incrementally (GreedyMini:
    rank first, see next_gm): the hashes of the kmer's T = k - t + 1 t-mers live in a ring (64 slots, locus x of the current kmer at
    slot (m_base + x) & 63), with the minimum hash and its leftmost and rightmost loci (a tie
    between them is resolved from the ring by the centre rule, util::mod_tie_offset). A slide
    hashes the one new t-mer; the whole ring is re-scanned only when the leftmost minimum leaves
    the kmer (the "re-scan" method of minimizer_iterator). A unique minimum samples offset
    x mod w in O(1). If the kmer passed to `next` is the first one after
    reset() / set_position(), the ring is filled from the kmer; every other call must pass the
    one-base slide of the previous kmer (the assertion against util::compute_minimizer checks
    the result in debug builds).
*/
template <typename Kmer, custom_order_mode Mode>
struct mod_minimizer_iterator {
    static_assert(is_mod_mode(Mode));

    mod_minimizer_iterator() {}

    mod_minimizer_iterator(uint64_t k, uint64_t m, hasher_type const& hasher,
                           uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_t(hasher.mod_t())
        , m_w(k - m + 1)
        , m_T(k - hasher.mod_t() + 1)
        , m_hasher_p(&hasher)  //
    {
        assert(k > 0 and m <= k and k <= 63);
        assert(m_t > 0 and m_t < m);
        set_position(position);
    }

    void set_position(uint64_t position) {
        m_position = position;
        reset();
    }

    void reset() { m_valid = false; }

    __attribute__((noinline)) minimizer_info next(Kmer kmer) {
        Kmer kmer_rc = kmer;
        kmer_rc.reverse_complement_inplace(m_k);
        if constexpr (is_gm_mode(Mode)) return next_gm(kmer, kmer_rc);

        if (!m_valid) {
            m_base = 0;
            Kmer window = kmer;
            for (uint64_t x = 0; x != m_T; ++x, window.drop_char()) {
                m_ring[x] = hash_at(window, kmer_rc, x);
            }
            rescan();
            m_valid = true;
        } else {
            m_base = (m_base + 1) & 63;
            Kmer window = kmer;
            window.drop_chars(m_T - 1);
            const uint64_t hash = hash_at(window, kmer_rc, m_T - 1);
            m_ring[slot(m_T - 1)] = hash;
            if (m_left == 0) {  // the leftmost minimum left the kmer
                rescan();
            } else {
                m_left -= 1;
                m_right -= 1;
                if (hash < m_min_hash) {
                    m_min_hash = hash;
                    m_left = m_T - 1;
                    m_right = m_T - 1;
                } else if (hash == m_min_hash) {
                    m_right = m_T - 1;
                }
            }
        }

        const uint64_t p =
            m_left == m_right
                ? m_left % m_w
                : util::mod_tie_offset(m_ring, m_base, m_w, m_min_hash, m_left, m_right,
                                       [&]() { return util::mod_le_rc(kmer, kmer_rc); });
        Kmer mwin = kmer;
        mwin.drop_chars(p);
        minimizer_info mini_info(util::canonical_mmer_at<Kmer>(mwin, kmer_rc, m_k, m_m, p),
                                 m_position + p, p);
        m_position += 1;

        assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, kmer_rc, m_k, m_m, hasher()));
        return mini_info;
    }

private:
    /*
        GreedyMini, rank first (as strand_invariant_minimizer_iterator): the ring keeps the
        16-bit table rank of every t-mer (one load from the Sigma plane of the kmer, rolled one
        base per call); a full hash is computed only for a t-mer whose rank equals the minimum
        rank (the hash of the current leftmost minimum lazily), and in a re-scan only for the
        t-mers sharing the minimum rank. A hash is (rank, tie-break) with the rank in its top
        bits, so the result is the one of comparing full hashes.
    */
    minimizer_info next_gm(Kmer const& kmer, Kmer const& kmer_rc) {
        if (!m_valid) {
            m_base = 0;
            m_sigma = sigma_of(kmer);
            for (uint64_t x = 0; x != m_T; ++x) {
                m_rank[x] = hasher().centre_rank(m_sigma, int(x));
            }
            rescan_gm(kmer);
            m_valid = true;
        } else {
            m_base = (m_base + 1) & 63;
            m_sigma = hasher().roll_sigma(m_sigma, kmer.at(m_k - 1), int(m_k));
            const uint16_t rank = hasher().centre_rank(m_sigma, int(m_T - 1));
            m_rank[slot(m_T - 1)] = rank;
            if (m_left == 0) {  // the leftmost minimum left the kmer
                rescan_gm(kmer);
            } else {
                m_left -= 1;
                m_right -= 1;
                if (rank < m_min_rank) {
                    m_min_rank = rank;
                    m_left = m_T - 1;
                    m_right = m_T - 1;
                    m_min_hash_valid = false;
                } else if (rank == m_min_rank) {
                    const uint64_t hash = hash_of(kmer, m_T - 1);
                    const uint64_t min_hash = current_min_hash(kmer);
                    if (hash < min_hash) {
                        m_min_hash = hash;
                        m_left = m_T - 1;
                        m_right = m_T - 1;
                    } else if (hash == min_hash) {
                        m_right = m_T - 1;
                    }
                }
            }
        }
#ifndef NDEBUG
        {
            const uint64_t kmask = m_k >= 64 ? ~uint64_t(0) : (uint64_t(1) << m_k) - 1;
            assert((m_sigma & kmask) == sigma_of(kmer));
        }
#endif

        uint64_t p;
        if (m_left == m_right) {
            p = m_left % m_w;
        } else {  // a tie at the minimum hash: the hashes of the minimal-rank loci, on demand
            const uint64_t min_hash = current_min_hash(kmer);
            uint64_t hashes[64];
            for (uint64_t x = m_left; x <= m_right; ++x) {
                hashes[x] = m_rank[slot(x)] == m_min_rank ? hash_of(kmer, x)
                                                          : constants::invalid_uint64;
            }
            p = util::mod_tie_offset(hashes, 0, m_w, min_hash, m_left, m_right,
                                     [&]() { return util::mod_le_rc(kmer, kmer_rc); });
        }
        Kmer mwin = kmer;
        mwin.drop_chars(p);
        minimizer_info mini_info(util::canonical_mmer_at<Kmer>(mwin, kmer_rc, m_k, m_m, p),
                                 m_position + p, p);
        m_position += 1;

        assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, kmer_rc, m_k, m_m, hasher()));
        return mini_info;
    }

    uint64_t sigma_of(Kmer const& kmer) const {
        uint64_t lo = uint64_t(kmer), hi = 0;
        if constexpr (sizeof(kmer.bits) > 8) hi = uint64_t(kmer.bits >> 64);
        return hasher().sigma_of_kmer(lo, hi, int(m_k));
    }

    /* full GreedyMini hash of the t-mer at locus x of the kmer */
    uint64_t hash_of(Kmer kmer, uint64_t x) const {
        kmer.drop_chars(x);
        kmer.take_chars(m_t);
        return hasher().hash_gm(uint64_t(kmer));
    }

    uint64_t current_min_hash(Kmer const& kmer) {
        if (!m_min_hash_valid) {
            m_min_hash = hash_of(kmer, m_left);
            m_min_hash_valid = true;
        }
        return m_min_hash;
    }

    /* minimum rank over the ring; full hashes only if several t-mers share it */
    void rescan_gm(Kmer const& kmer) {
        uint16_t min_rank = m_rank[slot(0)];
        uint64_t first = 0, count = 1;
        for (uint64_t x = 1; x != m_T; ++x) {
            const uint16_t r = m_rank[slot(x)];
            if (r < min_rank) {
                min_rank = r;
                first = x;
                count = 1;
            } else if (r == min_rank) {
                count += 1;
            }
        }
        m_min_rank = min_rank;
        m_left = first;
        m_right = first;
        if (count == 1) {
            m_min_hash_valid = false;
            return;
        }
        m_min_hash = hash_of(kmer, first);
        m_min_hash_valid = true;
        for (uint64_t x = first + 1; x != m_T; ++x) {
            if (m_rank[slot(x)] != min_rank) continue;
            const uint64_t hash = hash_of(kmer, x);
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_left = x;
                m_right = x;
            } else if (hash == m_min_hash) {
                m_right = x;
            }
        }
    }

    uint64_t m_k = 0, m_m = 0, m_t = 0, m_w = 0, m_T = 0;
    uint64_t m_position = 0;  // position of the next kmer in the sequence
    uint64_t m_base = 0;
    uint64_t m_min_hash = 0, m_left = 0, m_right = 0;  // loci in the current kmer
    bool m_valid = false;
    hasher_type const* m_hasher_p = nullptr;  // the hasher must outlive the iterator
    uint64_t m_ring[64];  // xor: the t-mer hashes
    /* GreedyMini (rank first): the t-mer ranks, the rolled Sigma plane, the minimum rank and
       the (lazily computed) full hash of the leftmost minimum */
    uint16_t m_rank[64];
    uint64_t m_sigma = 0;
    uint16_t m_min_rank = 0;
    bool m_min_hash_valid = false;

    hasher_type const& hasher() const { return *m_hasher_p; }
    uint64_t slot(uint64_t x) const { return (m_base + x) & 63; }

    /* hash of the t-mer at locus x; `window` is the kmer shifted so that x is at its low end */
    inline uint64_t hash_at(Kmer window, Kmer const& kmer_rc, uint64_t x) const {
        if constexpr (is_gm_mode(Mode)) {
            window.take_chars(m_t);
            return hasher().hash_gm(uint64_t(window));
        } else {
            return hasher().template plain<Mode>().hash(
                util::canonical_mmer_at<Kmer>(window, kmer_rc, m_k, m_t, x));
        }
    }

    void rescan() {
        m_min_hash = m_ring[slot(0)];
        m_left = 0;
        m_right = 0;
        for (uint64_t x = 1; x != m_T; ++x) {
            const uint64_t h = m_ring[slot(x)];
            if (h < m_min_hash) {
                m_min_hash = h;
                m_left = x;
                m_right = x;
            } else if (h == m_min_hash) {
                m_right = x;
            }
        }
    }
};

/* The iterator of the classical path for the order `Mode`. */
template <typename Kmer, custom_order_mode Mode>
using order_minimizer_iterator = std::conditional_t<
    is_mod_mode(Mode), mod_minimizer_iterator<Kmer, Mode>,
    std::conditional_t<Mode == custom_order_mode::GREEDYMINI,
                       strand_invariant_minimizer_iterator<Kmer>, minimizer_iterator<Kmer, Mode>>>;

}  // namespace sshash
