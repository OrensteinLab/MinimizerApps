#pragma once

#include <array>
#include <cassert>
#include <fstream>
#include <thread>
#include <cmath>  // for std::ceil on linux

#include "hash_util.hpp"

namespace sshash {

enum bucket_t : int {
    SINGLETON = 0,  // minimizer appears only once
    MIDLOAD = 1,    // minimizer appears > 1 but < 2^l times
    HEAVYLOAD = 3   // minimizer appears >= 2^l times
};

enum input_file_t { fasta, cf_seg };

struct streaming_query_report {
    streaming_query_report()
        : num_kmers(0)
        , num_positive_kmers(0)
        , num_negative_kmers(0)
        , num_invalid_kmers(0)
        , num_searches(0)
        , num_extensions(0)
        , num_skipped_singleton_lookups(0)
        , num_bucket_cache_hits(0) {}

    uint64_t num_kmers;
    uint64_t num_positive_kmers;
    uint64_t num_negative_kmers;
    uint64_t num_invalid_kmers;
    uint64_t num_searches;
    uint64_t num_extensions;

    /* seeds answered from the state of the last seed, saving a full lookup:
       negatives implied by a singleton bucket (skipped entirely), and seeds
       verified from the cached locate set (see `streaming_query::seed`) */
    uint64_t num_skipped_singleton_lookups;
    uint64_t num_bucket_cache_hits;
};

struct lookup_result {
    lookup_result(bool mf = true)
        : kmer_id(constants::invalid_uint64)
        , kmer_id_in_string(constants::invalid_uint64)
        , kmer_offset(constants::invalid_uint64)
        , kmer_orientation(constants::forward_orientation)

        , string_id(constants::invalid_uint64)
        , string_begin(constants::invalid_uint64)
        , string_end(constants::invalid_uint64)

        , minimizer_found(mf) {}

    uint64_t kmer_id;            // "absolute" kmer-id
    uint64_t kmer_id_in_string;  // "relative" kmer-id: 0 <= kmer_id_in_string < string_size,
                                 // where string_size = string_end - string_begin - k + 1
    uint64_t kmer_offset;
    int64_t kmer_orientation;

    uint64_t string_id;
    uint64_t string_begin;
    uint64_t string_end;

    bool minimizer_found;
};

inline std::ostream& operator<<(std::ostream& os, lookup_result const& res) {
    os << "  == kmer_id = " << res.kmer_id << '\n';
    os << "  == kmer_id_in_string = " << res.kmer_id_in_string << '\n';
    os << "  == kmer_offset = " << res.kmer_offset << '\n';
    os << "  == kmer_orientation = " << res.kmer_orientation << '\n';
    os << "  == string_id = " << res.string_id << '\n';
    os << "  == string_begin = " << res.string_begin << '\n';
    os << "  == string_end = " << res.string_end << '\n';
    os << "  == string_length = " << (res.string_end - res.string_begin) << '\n';
    os << "  == minimizer_found = " << (res.minimizer_found ? "true" : "false") << '\n';
    return os;
}

template <typename Kmer>
struct neighbourhood {
    std::array<lookup_result, Kmer::alphabet_size> forward;
    std::array<lookup_result, Kmer::alphabet_size> backward;
};

struct minimizer_info {
    minimizer_info()
        : minimizer(constants::invalid_uint64)
        , pos_in_seq(constants::invalid_uint64)
        , pos_in_kmer(constants::invalid_uint64) {}

    minimizer_info(uint64_t mm, uint64_t pk)
        : minimizer(mm), pos_in_seq(constants::invalid_uint64), pos_in_kmer(pk) {}

    minimizer_info(uint64_t mm, uint64_t ps, uint64_t pk)
        : minimizer(mm), pos_in_seq(ps), pos_in_kmer(pk) {}

    uint64_t minimizer;
    uint64_t pos_in_seq;
    uint64_t pos_in_kmer;

    bool operator==(minimizer_info rhs) const {
        return minimizer == rhs.minimizer and    //
               pos_in_seq == rhs.pos_in_seq and  //
               pos_in_kmer == rhs.pos_in_kmer;   //
    }
    bool operator!=(minimizer_info rhs) const { return !(*this == rhs); }
};

[[maybe_unused]] static bool equal_lookup_result(lookup_result expected, lookup_result got) {
    bool good = true;
    if (expected.kmer_id != got.kmer_id) {
        std::cout << "expected kmer_id " << expected.kmer_id << " but got " << got.kmer_id
                  << std::endl;
        good = false;
    }
    if (expected.kmer_id_in_string != got.kmer_id_in_string) {
        std::cout << "expected kmer_id_in_string " << expected.kmer_id_in_string << " but got "
                  << got.kmer_id_in_string << std::endl;
        good = false;
    }
    if (got.kmer_id != constants::invalid_uint64 and
        expected.kmer_orientation != got.kmer_orientation) {
        std::cout << "expected kmer_orientation " << expected.kmer_orientation << " but got "
                  << got.kmer_orientation << std::endl;
        good = false;
    }
    if (expected.string_id != got.string_id) {
        std::cout << "expected string_id " << expected.string_id << " but got " << got.string_id
                  << std::endl;
        good = false;
    }
    if (expected.string_begin != got.string_begin) {
        std::cout << "expected string_begin " << expected.string_begin << " but got "
                  << got.string_begin << std::endl;
        good = false;
    }
    if (expected.string_end != got.string_end) {
        std::cout << "expected string_end " << expected.string_end << " but got " << got.string_end
                  << std::endl;
        good = false;
    }
    return good;
}

struct build_configuration {
    build_configuration()
        : k(31)
        , m(20)
        , seed(constants::seed)
        , num_threads(1)
        , ram_limit_in_GiB(constants::default_ram_limit_in_GiB)

        , lambda(constants::lambda)

        , weighted(false)
        , verbose(true)

        , tmp_dirname(constants::default_tmp_dirname)

        , order(custom_order_mode::NATIVE)
        , gm_alphabet(GmAlphabet::BINARY)
        , mod_minimizer(false)
        , t(0)

    {}

    uint64_t k;  // kmer length
    uint64_t m;  // minimizer length
    uint64_t seed;
    uint64_t num_threads;
    uint64_t ram_limit_in_GiB;

    double lambda;  // drive PTHash trade-off

    bool weighted;
    bool verbose;

    std::string tmp_dirname;

    /* Custom minimizer order (see orders/README.md); NATIVE reproduces
       upstream's own mixer_64-based order unchanged. */
    custom_order_mode order;
    std::string greedy_mini_filename;  // required when order == GREEDYMINI
    GmAlphabet gm_alphabet;

    /* Mod-minimizer (Groot Koerkamp & Pibiri, WABI 2024): when true, the
       minimizer of a kmer is not the smallest m-mer directly, but the m-mer
       at the position mod-sampling selects using `order` as the order O_t on
       t-mers (t < m). See util::compute_mod_minimizer_single_strand. */
    bool mod_minimizer;
    uint64_t t;  // t-mer length for O_t; only meaningful when mod_minimizer

    void print() const {
        std::cout << "k = " << k                                              //
                  << ", m = " << m                                            //
                  << ", seed = " << seed                                      //
                  << ", num_threads = " << num_threads                        //
                  << ", ram_limit_in_GiB = " << ram_limit_in_GiB              //
                  << ", lambda = " << lambda                                  //
                  << ", weighted = " << (weighted ? "true" : "false")         //
                  << ", verbose = " << (verbose ? "true" : "false")           //
                  << ", tmp_dirname = '" << tmp_dirname << "'"                //
                  << ", order = '" << order_name(order, mod_minimizer) << "'"              //
                  << (mod_minimizer ? ", t = " + std::to_string(t) : "") << std::endl;     //
    }
};

namespace util {

static inline void check_version_number(essentials::version_number const& vnum) {
    if (vnum.x != constants::current_version_number::x) {
        throw std::runtime_error("MAJOR index version mismatch: SSHash index needs rebuilding");
    }
}

static inline uint64_t get_seed_for_hash_function(build_configuration const& build_config) {
    static const uint64_t my_favourite_seed = 1234567890;
    return build_config.seed != my_favourite_seed ? my_favourite_seed : ~my_favourite_seed;
}

[[maybe_unused]] static bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

template <typename Kmer>
[[maybe_unused]] static Kmer string_to_uint_kmer(char const* str, uint64_t k) {
    assert(k <= Kmer::max_k);
    Kmer x = 0;
    for (uint64_t i = 0; i != k; ++i) x.set(i, Kmer::char_to_uint(str[i]));
    return x;
}

template <typename Kmer>
static void uint_kmer_to_string(Kmer x, char* str, uint64_t k) {
    assert(k <= Kmer::max_k);
    for (uint64_t i = 0; i != k; ++i) str[i] = Kmer::uint64_to_char(x.pop_char());
}

template <typename Kmer>
[[maybe_unused]] static std::string uint_kmer_to_string(Kmer x, uint64_t k) {
    assert(k <= Kmer::max_k);
    std::string str;
    str.resize(k);
    uint_kmer_to_string(x, str.data(), k);
    return str;
}

template <typename Kmer>
[[maybe_unused]] static std::string uint_minimizer_to_string(uint64_t minimizer, uint64_t m) {
    assert(m <= Kmer::max_m);
    std::string str;
    str.resize(m);
    Kmer x = minimizer;
    uint_kmer_to_string(x, str.data(), m);
    return str;
}

template <typename Kmer>
[[maybe_unused]] static bool is_valid(char const* str, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        if (!Kmer::is_valid(str[i])) return false;
    }
    return true;
}

template <typename Kmer>
static Kmer read_kmer_at(bits::bit_vector const& bv, const uint64_t k, const uint64_t pos) {
    static_assert(Kmer::uint_kmer_bits % 64 == 0);
    Kmer kmer = 0;
    for (int i = Kmer::uint_kmer_bits - 64; i >= 0; i -= 64) {
        if (pos + i < bv.num_bits()) kmer.append64(bv.get_word64(pos + i));
    }
    kmer.take(Kmer::bits_per_char * k);
    return kmer;
}

/*
    The canonical m-mer: the smaller of the m-mer and its reverse
    complement, under the numeric order on the packed encoding. For alphabets
    that have no reverse complement (e.g. amino acids) this is the identity.
*/
template <typename Kmer>
inline uint64_t canonical_mmer(const uint64_t mmer, const uint64_t m) {
    if constexpr (!Kmer::has_reverse_complement) {
        (void)m;
        return mmer;
    } else {
        return std::min(mmer, Kmer::reverse_complement_mmer(mmer, m));
    }
}

/*
    The canonical m-mer at locus i of a kmer x, given rc(x).

    Uses rc(x_i) = rc(x)_{k-m-i}: the reverse complements of all k-m+1 loci are
    windows of the single reverse complement of the whole kmer, so the kmer is
    reverse-complemented once instead of once per locus. `window` is x shifted so
    that locus i sits at its low end -- the forward sliding window the callers
    already carry along.

    When the whole kmer fits in a word -- always so for the 64-bit kmer type, and
    for k <= 32 with the wider one -- the reverse window is extracted with a plain
    64-bit shift; otherwise with a shift of the wide type, which costs a few
    instructions more but still beats reverse-complementing every locus.
*/
template <typename Kmer>
inline uint64_t canonical_mmer_at(Kmer window, [[maybe_unused]] Kmer const& kmer_rc,
                                  [[maybe_unused]] const uint64_t k, const uint64_t m,
                                  [[maybe_unused]] const uint64_t i)  //
{
    window.take_chars(m);
    if constexpr (!Kmer::has_reverse_complement) {
        return uint64_t(window);
    } else {
        constexpr uint64_t b = Kmer::bits_per_char;
        assert(i + m <= k);
        uint64_t rc;
        if (k * b <= 64) {
            const uint64_t mask = m * b == 64 ? ~uint64_t(0) : (uint64_t(1) << (m * b)) - 1;
            rc = (uint64_t(kmer_rc) >> (b * (k - m - i))) & mask;
        } else {
            Kmer mmer_rc = kmer_rc;
            mmer_rc.drop_chars(k - m - i);
            mmer_rc.take_chars(m);
            rc = uint64_t(mmer_rc);
        }
        assert(rc == Kmer::reverse_complement_mmer(uint64_t(window), m));
        return std::min(uint64_t(window), rc);
    }
}

/*
    True if `kmer` is the canonical one of the pair (kmer, rc(kmer)). A kmer that
    equals its own reverse complement (possible only for even k) is deemed
    canonical, so that the answer is always well defined.
*/
template <typename Kmer>
inline bool is_canonical(Kmer kmer, const uint64_t k) {
    if constexpr (!Kmer::has_reverse_complement) {
        (void)k;
        return true;
    } else {
        Kmer kmer_rc = kmer;
        kmer_rc.reverse_complement_inplace(k);
        return !(kmer_rc < kmer);
    }
}

/*
    This implements the random minimizer of a kmer x: the locus i* in [0, k-m]
    minimizing h(kappa(i)), where kappa(i) := min(x_i, rc(x_i)) is the canonical
    m-mer at locus i. The minimizer's value is kappa(i*).

    Because kappa(i) is invariant under reverse-complementation and the sequence
    (kappa(0), ..., kappa(k-m)) of rc(x) is the reversal of that of x, the locus
    selected for rc(x) is the mirror image of the one selected for x, and the
    value is literally the same m-mer. So x and rc(x) always land in the same
    bucket -- a lookup costs a single probe -- while the minimizer remains a
    minimizer of the window under a single random order, which is what keeps the
    density (and hence the number of super-kmers) at that of the plain forward
    minimizer, 2/(k-m+2).

    When the alphabet has no reverse complement, kappa is the identity and this
    reduces exactly to the plain forward minimizer.

    Ties -- h(kappa(i)) == h(kappa(j)) for i != j, which happens when x_i = x_j
    or x_i = rc(x_j) -- must be broken in a mirror-equivariant way, or x and
    rc(x) would be sent to different buckets, which is a correctness failure and
    not merely a density one. No leftmost rule can be equivariant: the mirror
    reverses the order of the tied loci, so "leftmost" names different loci on
    the two strands. We use the tie-break of [Cologni and Pibiri, "Canonical
    Schemes and Minimizers", Proposition 24]: among the tied loci, take the one
    closest to the centre (k-m)/2 of the window -- the one reference point both
    strands agree on -- and when two are equally close (they are then mirror
    images i and k-m-i), take the smaller index if x <= rc(x) and the larger
    otherwise. This rule is mirror-equivariant (for odd k, where x != rc(x)
    always) and, unlike the previous frame-of-the-canonical-kmer rule, it is
    also *forward*: the distance of a fixed locus to the centre is monotone as
    the window slides right, so a locus at least as close as another stays so,
    and the anchor never moves backwards. Forwardness makes the number of
    super-kmers equal the number of sampled positions -- nothing is ever
    abandoned and re-opened.

    Since rc(x_i) = rc(x)_{k-m-i}, the reverse complements of all k-m+1 loci are
    just windows of the single reverse complement of the whole kmer: `kmer_rc` is
    reverse-complemented once and each rc(x_i) is read out of it with a shift,
    instead of reverse-complementing every locus separately. The caller usually
    has rc(x) already -- a lookup needs it anyway to recognise which orientation
    it found -- in which case the reverse complementation is free.
*/
/*
    A tie at the minimum (rare, rate ~ sigma^{-m/2}): all tied loci carry the
    same class, so only the sampled position is at stake, not the value. Take
    the tied locus closest to the centre of the window; the doubled distance
    |2i - (k-m)| avoids the half-integer centre of an odd window.

    Deliberately not inlined: this runs on ~1e-4 of the windows, and letting its
    loop inline into the lookup path costs measurably more than the call.
*/
template <typename Kmer>
__attribute__((noinline)) minimizer_info resolve_tie(Kmer const& kmer, Kmer const& kmer_rc,
                                                     const uint64_t k, const uint64_t m,
                                                     hasher_type const& hasher,
                                                     const uint64_t min_hash,
                                                     const uint64_t minimizer,
                                                     const uint64_t leftmost,
                                                     const uint64_t rightmost)  //
{
    const uint64_t two_c = k - m;
    uint64_t best_dist = constants::invalid_uint64;
    uint64_t lo = leftmost;
    uint64_t hi = leftmost;
    Kmer window = kmer;
    window.drop_chars(leftmost);
    for (uint64_t i = leftmost; i <= rightmost; ++i) {
        if (hasher.hash(canonical_mmer_at<Kmer>(window, kmer_rc, k, m, i)) == min_hash) {
            assert(canonical_mmer_at<Kmer>(window, kmer_rc, k, m, i) == minimizer);
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
    /* lo and hi are the two equally-closest loci (mirror images), or one locus */
    const uint64_t chosen = (lo == hi or !(kmer_rc < kmer)) ? lo : hi;
    return {minimizer, chosen};
}

template <typename Kmer>
minimizer_info compute_minimizer(Kmer kmer, Kmer const& kmer_rc, const uint64_t k, const uint64_t m,
                                 hasher_type const& hasher)  //
{
    assert(m <= Kmer::max_m);
    assert(m <= k);

    auto kappa = [&](Kmer const& window, const uint64_t i) {
        return canonical_mmer_at<Kmer>(window, kmer_rc, k, m, i);
    };

    /* The first locus is peeled off the loop so that `min_hash` starts out at a
       real hash value: initializing it to invalid_uint64 would make an actual
       hash of invalid_uint64 register as a tie rather than as the minimum. */
    Kmer window = kmer;
    uint64_t minimizer = kappa(window, 0);
    uint64_t min_hash = hasher.hash(minimizer);
    uint64_t leftmost = 0;
    uint64_t rightmost = 0;
    bool tie = false;
    window.drop_char();

    for (uint64_t i = 1; i != k - m + 1; ++i) {
        uint64_t value = kappa(window, i);
        uint64_t hash = hasher.hash(value);
        if (hash < min_hash) {  // leftmost
            min_hash = hash;
            minimizer = value;
            leftmost = i;
            rightmost = i;
            tie = false;
        } else if (hash == min_hash) {
            rightmost = i;
            tie = true;
        }
        window.drop_char();
    }

    if (__builtin_expect(tie, false)) {
        return resolve_tie(kmer, kmer_rc, k, m, hasher, min_hash, minimizer, leftmost, rightmost);
    }
    return {minimizer, leftmost};
}

template <typename Kmer>
minimizer_info compute_minimizer(Kmer kmer, const uint64_t k, const uint64_t m,
                                 hasher_type const& hasher)  //
{
    Kmer kmer_rc = kmer;
    kmer_rc.reverse_complement_inplace(k);
    return compute_minimizer(kmer, kmer_rc, k, m, hasher);
}

/*
    Mod-sampling (Groot Koerkamp & Pibiri, "The mod-minimizer: a simple and
    efficient sampling algorithm for long k-mers", WABI 2024, Algorithm 6 /
    Definition 7), specialised to sshash's own (k, m) model: the paper's
    generic "k-mer being sampled" is sshash's m-mer, and its window count w is
    sshash's derived w = k-m+1. Instead of directly ranking the w candidate
    m-mers, mod-sampling ranks the (more numerous) t-mers (t < m) across the
    *whole* kmer by a separate order `hasher` = O_t, and reduces the winning
    t-mer's locus mod w to pick which m-mer is sampled:

        o_min = +inf; p = +inf
        for i = 0; i < k-t+1; i++:
            o = O_t(kmer[i..i+t))
            if o == o_min: p = min(p, i mod w)
            elif o < o_min: o_min = o; p = i mod w
        return p

    (`k-t+1`, the number of t-mer loci in the whole kmer, is the paper's
    `w+k-t` after substituting w=k-m+1.) The sampled minimizer is the m-mer at
    position p -- a raw, single-strand substring, not separately canonicalized
    per locus.

    Mod-sampling as stated operates on one fixed forward strand only; its
    "prefer the smaller i mod w on a tie" rule is not symmetric under
    reversal, so its per-locus values do NOT satisfy the mirror relationship
    (kappa(i) in kmer's frame == kappa(k-m-i) in kmer_rc's frame) that
    compute_minimizer's canonical_mmer_at/resolve_tie machinery relies on --
    concretely, running this routine independently on kmer and on kmer_rc can
    select genuinely different m-mers, for *any* input, not just at ties.

    Worse, combining the two independently (taking the smaller of the two
    resulting m-mers, as sshash's own *canonical parsing* modality does for
    the classical scheme -- see the 2022 SSHash paper, "Double Strandedness")
    is not just a mirror-equivariance problem here: it can make the SAME
    physical kmer get claimed by two different super-kmers during
    construction (a "backward jump", Section 4.1 of the mod-minimizer paper),
    corrupting the index (duplicate kmer ids) rather than merely picking a
    different-but-consistent bucket. Verified empirically: even t=10 (which
    Lemma 13 of the mod-minimizer paper proves gives a true, forward
    minimizer scheme on a single strand) still hits this once combined with
    its own reverse-complement's independent computation.

    So mod-minimizer instead uses sshash's *regular parsing* modality (the
    2022 paper's other, simpler modality, also still real code today in
    jermp/sshash-lite's streaming_query_regular_parsing): construction keeps
    a kmer's minimizer as whichever single strand it was actually walked in
    (this function, called once, no combination -- see
    dictionary_builder::compute_minimizer_tuples), and Lookup instead "will
    first search for g and -- only if not found -- will also search for ĝ.
    This effectively doubles the query time for Lookup in the worst case"
    (ibid., quoted verbatim) -- see dictionary::lookup's two-probe fallback.

    Deliberately a straightforward O(k-t+1) computation per kmer (mirroring
    compute_minimizer's own role as the brute-force reference), not an
    amortized-O(1) incremental scheme like `minimizer_iterator`; see the
    mod-minimizer integration notes for why this is an acceptable starting
    point.
*/
template <typename Kmer>
minimizer_info compute_mod_minimizer_single_strand(Kmer kmer, const uint64_t k, const uint64_t m,
                                                   const uint64_t t,
                                                   hasher_type const& hasher)  //
{
    assert(t <= m and m <= k);
    assert(m <= Kmer::max_m);

    const uint64_t w = k - m + 1;
    const uint64_t num_t_loci = k - t + 1;

    auto tau = [&](Kmer const& window) {
        Kmer tmer = window;
        tmer.take_chars(t);
        return uint64_t(tmer);
    };

    /* first locus peeled off the loop, same reasoning as compute_minimizer:
       o_min must start out as a real hash value, not a sentinel that a real
       hash could legitimately collide with. */
    Kmer window = kmer;
    uint64_t o_min = hasher.hash(tau(window));
    uint64_t p = 0;  // 0 % w == 0
    window.drop_char();

    for (uint64_t i = 1; i != num_t_loci; ++i) {
        uint64_t o = hasher.hash(tau(window));
        uint64_t pos_mod_w = i % w;
        if (o < o_min) {
            o_min = o;
            p = pos_mod_w;
        } else if (o == o_min and pos_mod_w < p) {
            p = pos_mod_w;
        }
        window.drop_char();
    }

    Kmer final_window = kmer;
    final_window.drop_chars(p);
    final_window.take_chars(m);
    return {uint64_t(final_window), p};
}

}  // namespace util

struct buffered_lines_iterator {
    static const uint64_t BUFFER_SIZE = 1024;

    buffered_lines_iterator(std::istream& is, uint64_t buffer_size = BUFFER_SIZE)
        : m_is(is), m_buffer_size(buffer_size), m_read_chars(0) {}

    bool fill_buffer(std::string& buffer) {
        if (buffer.size() >= m_buffer_size) return false;

        bool empty_line_was_read = false;
        uint64_t size = buffer.size();
        buffer.resize(m_buffer_size);

        char* ptr = buffer.data() + size;
        while (size < m_buffer_size) {
            // read until '\n' or rest is filled
            uint64_t rest = m_buffer_size - size;
            m_is.getline(ptr, rest, '\n');
            uint64_t read_chars = m_is.gcount();
            m_read_chars += read_chars;

            if (!m_is) {
                if (read_chars + 1 == rest) {  // '\n' not found
                    m_is.clear();
                    size += read_chars;
                    break;
                }
            } else if (!eof()) {
                assert(read_chars > 0);
                --read_chars;  // discard the delimiter
            }

            if (read_chars == 0) {  // empty line was read
                empty_line_was_read = true;
                break;
            }

            size += read_chars;
            ptr += read_chars;
        }

        buffer.resize(size);
        return empty_line_was_read;
    }

    bool eof() const { return m_is.eof(); }

    uint64_t read_chars() const { return m_read_chars; }

private:
    std::istream& m_is;
    uint64_t m_buffer_size;
    uint64_t m_read_chars;
};

}  // namespace sshash
