#pragma once

#include <array>
#include <cassert>
#include <fstream>
#include <thread>
#include <cmath>  // for std::ceil on linux

#include "hash_util.hpp"
#include "../../common/minimizer-orders/mod_minimizer_rule.h"

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
        , bucket_stats(false)

        , tmp_dirname(constants::default_tmp_dirname)

        , order(custom_order_mode::NATIVE)
        , gm_proj(gm_projection::AT0_CG1)
        , mod_t(0)
        , regular_parsing(false)

    {}

    uint64_t k;  // kmer length
    uint64_t m;  // minimizer length
    uint64_t seed;
    uint64_t num_threads;
    uint64_t ram_limit_in_GiB;

    double lambda;  // drive PTHash trade-off

    bool weighted;
    bool verbose;
    bool bucket_stats;  // also print the bucket-size and super-k-mer histograms

    std::string tmp_dirname;

    /* Custom minimizer order (see orders/README.md); NATIVE reproduces
       upstream's own mixer_64-based order unchanged. */
    custom_order_mode order;
    std::string greedy_mini_filename;  // required when order == GREEDYMINI
    gm_projection gm_proj;  // DNA->binary Sigma polarity
    uint64_t mod_t;         // t-mer length of the mod-minimizer orders (0 otherwise)

    /* Forces the two-probe "regular parsing" modality (see
       util::compute_minimizer_single_strand) -- an explicit comparison mode,
       not a production default. */
    bool regular_parsing;

    void print() const {
        std::cout << "k = " << k                                              //
                  << ", m = " << m                                            //
                  << ", seed = " << seed                                      //
                  << ", num_threads = " << num_threads                        //
                  << ", ram_limit_in_GiB = " << ram_limit_in_GiB              //
                  << ", lambda = " << lambda                                  //
                  << ", weighted = " << (weighted ? "true" : "false")         //
                  << ", verbose = " << (verbose ? "true" : "false")           //
                  << ", bucket_stats = " << (bucket_stats ? "true" : "false") //
                  << ", tmp_dirname = '" << tmp_dirname << "'"                //
                  << ", order = '" << custom_order_mode_name(order) << "'"              //
                  << (is_gm_mode(order) ? std::string(", gm_projection = ") + gm_projection_name(gm_proj) : std::string()) //
                  << (is_mod_mode(order) ? ", mod_t = " + std::to_string(mod_t) : std::string()) //
                  << (regular_parsing ? ", regular_parsing = true" : "") << std::endl;      //
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

    Always inlined: it is the heart of the minimizer iterators' per-k-mer step, and whether the
    compiler's size heuristics inline it depends on everything else in the translation unit
    (out of line, it made the native streaming query 7% slower than upstream's).
*/
template <typename Kmer>
__attribute__((always_inline)) inline uint64_t canonical_mmer_at(Kmer window, [[maybe_unused]] Kmer const& kmer_rc,
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
    Hash of the canonical m-mer at locus i. For the strand-invariant order (GreedyMini, the hasher
    itself: hash(x) == hash(rc(x))) the raw forward window is hashed directly and the
    canonicalization -- the min with the reverse-complement window -- is skipped: it only matters
    for the VALUE of a locus that turns out to be the minimum, which callers compute lazily with
    canonical_mmer_at. Every other order (a plain_hasher<Mode>) hashes canonical_mmer_at(...).
*/
template <typename Kmer, typename Hasher>
inline uint64_t hash_mmer_at(Hasher const& hasher, Kmer window, Kmer const& kmer_rc,
                             const uint64_t k, const uint64_t m, const uint64_t i)  //
{
    if constexpr (strand_invariant_hasher_v<Hasher>) {
        window.take_chars(m);
        return hasher.hash_gm(uint64_t(window));
    } else {
        return hasher.hash(canonical_mmer_at<Kmer>(window, kmer_rc, k, m, i));
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
template <typename Kmer, typename Hasher>
__attribute__((noinline)) minimizer_info resolve_tie(Kmer const& kmer, Kmer const& kmer_rc,
                                                     const uint64_t k, const uint64_t m,
                                                     Hasher const& hasher,
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
        if (hash_mmer_at<Kmer>(hasher, window, kmer_rc, k, m, i) == min_hash) {
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

/* The definition above, evaluated literally: the full hash of every locus, with `hasher` a
   plain_hasher<Mode> or the GreedyMini hasher. Slow for GreedyMini (see compute_minimizer_ranked,
   which skips the tie-break of the loci that cannot win); the reference the tests compare
   compute_minimizer and the iterators against. */
template <typename Kmer, typename Hasher>
__attribute__((always_inline)) inline minimizer_info compute_minimizer_all_loci(
    Kmer kmer, Kmer const& kmer_rc, const uint64_t k, const uint64_t m,
    Hasher const& hasher)  //
{
    assert(m <= Kmer::max_m);
    assert(m <= k);

    auto kappa = [&](Kmer const& window, const uint64_t i) {
        return canonical_mmer_at<Kmer>(window, kmer_rc, k, m, i);
    };
    /* the strand-invariant order hashes raw windows and canonicalizes the winner only */
    constexpr bool raw = strand_invariant_hasher_v<Hasher>;

    /* The first locus is peeled off the loop so that `min_hash` starts out at a
       real hash value: initializing it to invalid_uint64 would make an actual
       hash of invalid_uint64 register as a tie rather than as the minimum. */
    Kmer window = kmer;
    uint64_t minimizer = 0;
    uint64_t min_hash;
    if constexpr (raw) {
        min_hash = hash_mmer_at<Kmer>(hasher, window, kmer_rc, k, m, 0);
    } else {
        minimizer = kappa(window, 0);
        min_hash = hasher.hash(minimizer);
    }
    uint64_t leftmost = 0;
    uint64_t rightmost = 0;
    bool tie = false;
    window.drop_char();

    for (uint64_t i = 1; i != k - m + 1; ++i) {
        uint64_t value = 0;
        uint64_t hash;
        if constexpr (raw) {
            hash = hash_mmer_at<Kmer>(hasher, window, kmer_rc, k, m, i);
        } else {
            value = kappa(window, i);
            hash = hasher.hash(value);
        }
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
    if constexpr (raw) {
        Kmer w = kmer;
        w.drop_chars(leftmost);
        minimizer = kappa(w, leftmost);
    }

    if (__builtin_expect(tie, false)) {
        return resolve_tie(kmer, kmer_rc, k, m, hasher, min_hash, minimizer, leftmost, rightmost);
    }
    return {minimizer, leftmost};
}

/*
    compute_minimizer_all_loci for the GreedyMini order, where the loci are first ranked by the
    rank field of their hash (one table load each, from the Sigma plane of the whole kmer) and
    only the loci sharing the minimum rank get a full hash: a hash is (rank, tie-break) with the
    rank in its top bits, so the other loci cannot be the minimum. Same result, same tie handling.
*/
template <typename Kmer>
minimizer_info compute_minimizer_ranked(Kmer kmer, Kmer const& kmer_rc, const uint64_t k,
                                        const uint64_t m, hasher_type const& hasher)  //
{
    assert(m <= Kmer::max_m);
    assert(m <= k);
    assert(k <= 63);

    uint64_t lo = uint64_t(kmer), hi = 0;
    if constexpr (sizeof(kmer.bits) > 8) hi = uint64_t(kmer.bits >> 64);
    const uint64_t sigma = hasher.sigma_of_kmer(lo, hi, int(k));
    const uint64_t num_loci = k - m + 1;

    uint16_t ranks[64];
    uint16_t min_rank = 0xFFFF;
    for (uint64_t i = 0; i != num_loci; ++i) {
        ranks[i] = hasher.centre_rank(sigma, int(i));
        min_rank = ranks[i] < min_rank ? ranks[i] : min_rank;
    }

    uint64_t min_hash = constants::invalid_uint64;
    uint64_t leftmost = 0;
    uint64_t rightmost = 0;
    bool tie = false;
    bool first = true;
    Kmer window = kmer;
    for (uint64_t i = 0; i != num_loci; ++i, window.drop_char()) {
        if (ranks[i] != min_rank) continue;
        const uint64_t hash = hash_mmer_at<Kmer>(hasher, window, kmer_rc, k, m, i);
        if (first or hash < min_hash) {  // leftmost
            first = false;
            min_hash = hash;
            leftmost = i;
            rightmost = i;
            tie = false;
        } else if (hash == min_hash) {
            rightmost = i;
            tie = true;
        }
    }
    Kmer w = kmer;
    w.drop_chars(leftmost);
    const uint64_t minimizer = canonical_mmer_at<Kmer>(w, kmer_rc, k, m, leftmost);

    if (__builtin_expect(tie, false)) {
        return resolve_tie(kmer, kmer_rc, k, m, hasher, min_hash, minimizer, leftmost, rightmost);
    }
    return {minimizer, leftmost};
}

/*
    v <= RC(v) for the centre tie rule of the mod-minimizer: the same comparison of packed words
    as resolve_tie's `!(kmer_rc < kmer)`; true without strands (plain mod-sampling).
*/
template <typename Kmer>
inline bool mod_le_rc(Kmer const& kmer, Kmer const& kmer_rc) {
    if constexpr (!Kmer::has_reverse_complement) {
        (void)kmer, (void)kmer_rc;
        return true;
    } else {
        return !(kmer_rc < kmer);
    }
}

/* Tie at the minimum of the mod-minimizer: every locus in [x_l, x_r] whose hash
   (hashes[(base + x) & 63]) is min_hash offers the offset x mod w; the centre rule
   (mod_rule::centre_tie_offset) picks one, calling le_rc() only on an exact two-way tie. */
template <typename LeRc>
inline uint64_t mod_tie_offset(uint64_t const* hashes, uint64_t base, uint64_t w, uint64_t min_hash,
                               uint64_t x_l, uint64_t x_r, LeRc le_rc) {
    return mod_rule::centre_tie_offset(
        x_l, x_r, w, [&](uint64_t x) { return hashes[(base + x) & 63] == min_hash; }, le_rc);
}

/*
    Strand-symmetric mod-minimizer (mod-sampling, Groot Koerkamp and Pibiri, WABI 2024, with the
    centre tie rule of compute_minimizer; shared definition with the minimap2 and kraken2 forks
    in common/minimizer-orders/mod_minimizer_rule.h). The kmer's T = k - t + 1 t-mers are ranked
    by the order of `hasher` (canonical t-mers; GreedyMini hashes the raw t-mer, its hash being
    strand invariant); a minimal t-mer at locus x offers the m-mer at offset x mod w,
    w = k - m + 1. The builder requires t == m (mod w), so T-1 == w-1 (mod w) and the mirror
    locus T-1-x in rc(kmer) offers the mirror offset w-1-p, whose m-mer has the same canonical
    value: a kmer and its reverse complement share the bucket and the mirrored position.

    A unique minimum samples its offset. A tie at the minimum (equal canonical t-mers, ~1e-3 of
    the kmers on real data) is broken as compute_minimizer's own (resolve_tie), on the offered
    offsets: the one closest to the centre (w-1)/2 of the kmer's m-mers, and between two
    equally close ones (mirror images) the smaller if kmer <= rc(kmer) (packed words), the
    larger otherwise. The scheme keeps mod-sampling's density and is forward -- a tied t-mer's
    offset stays put while the window slides, new candidates enter on the right, and the
    candidate closest to a centre moving right never lies left of the previous choice (the
    builder checks it, see minimizers_tuples).

    The definition evaluated directly (O(T) hashes); mod_minimizer_iterator computes the same
    incrementally.
*/
template <custom_order_mode Mode, typename Kmer>
minimizer_info compute_mod_minimizer(Kmer const& kmer, Kmer const& kmer_rc, const uint64_t k,
                                     const uint64_t m, hasher_type const& hasher)  //
{
    static_assert(is_mod_mode(Mode));
    const uint64_t t = hasher.mod_t();
    assert(t > 0 and t < m and m <= k);
    const uint64_t w = k - m + 1;
    const uint64_t T = k - t + 1;
    assert(T <= 64);  // k <= 63

    uint64_t min_hash = constants::invalid_uint64;
    uint64_t x_l = 0, x_r = 0;
    uint64_t hashes[64];
    if constexpr (is_gm_mode(Mode)) {
        /* Rank first (as compute_minimizer_ranked): a GreedyMini hash is (rank, tie-break) with
           the rank in its top bits, so the t-mers are ranked from the Sigma plane of the whole
           kmer (one table load each) and only those sharing the minimum rank get a full hash;
           the others keep hashes[x] = invalid (never equal to the minimum). */
        uint64_t lo = uint64_t(kmer), hi = 0;
        if constexpr (sizeof(kmer.bits) > 8) hi = uint64_t(kmer.bits >> 64);
        const uint64_t sigma = hasher.sigma_of_kmer(lo, hi, int(k));
        uint16_t ranks[64];
        uint16_t min_rank = 0xFFFF;
        for (uint64_t x = 0; x != T; ++x) {
            ranks[x] = hasher.centre_rank(sigma, int(x));
            min_rank = ranks[x] < min_rank ? ranks[x] : min_rank;
        }
        bool first = true;
        Kmer window = kmer;
        for (uint64_t x = 0; x != T; ++x, window.drop_char()) {
            hashes[x] = constants::invalid_uint64;
            if (ranks[x] != min_rank) continue;
            Kmer tmer = window;
            tmer.take_chars(t);
            const uint64_t hash = hasher.hash_gm(uint64_t(tmer));
            hashes[x] = hash;
            if (first or hash < min_hash) {
                first = false;
                min_hash = hash;
                x_l = x;
                x_r = x;
            } else if (hash == min_hash) {
                x_r = x;
            }
        }
    } else {
        Kmer window = kmer;
        for (uint64_t x = 0; x != T; ++x, window.drop_char()) {
            const uint64_t hash =
                hasher.plain<Mode>().hash(canonical_mmer_at<Kmer>(window, kmer_rc, k, t, x));
            hashes[x] = hash;
            if (x == 0 or hash < min_hash) {
                min_hash = hash;
                x_l = x;
                x_r = x;
            } else if (hash == min_hash) {
                x_r = x;
            }
        }
    }
    const uint64_t p = x_l == x_r ? x_l % w
                                  : mod_tie_offset(hashes, 0, w, min_hash, x_l, x_r,
                                                   [&]() { return mod_le_rc(kmer, kmer_rc); });
    Kmer mwin = kmer;
    mwin.drop_chars(p);
    return {canonical_mmer_at<Kmer>(mwin, kmer_rc, k, m, p), p};
}

/* The minimizer under the order `Mode`, known at compile time. */
template <custom_order_mode Mode, typename Kmer>
inline minimizer_info compute_minimizer_for(Kmer const& kmer, Kmer const& kmer_rc,
                                            const uint64_t k, const uint64_t m,
                                            hasher_type const& hasher)  //
{
    if constexpr (is_mod_mode(Mode)) {
        return compute_mod_minimizer<Mode, Kmer>(kmer, kmer_rc, k, m, hasher);
    } else if constexpr (Mode == custom_order_mode::GREEDYMINI) {
        return compute_minimizer_ranked<Kmer>(kmer, kmer_rc, k, m, hasher);
    } else {
        return compute_minimizer_all_loci<Kmer>(kmer, kmer_rc, k, m, hasher.plain<Mode>());
    }
}

/* The same, with the order of `hasher` chosen at run time (not inlined: the callers' hot path is
   the native order, see compute_minimizer). */
template <typename Kmer>
__attribute__((noinline)) minimizer_info compute_minimizer_dispatched(
    Kmer const& kmer, Kmer const& kmer_rc, const uint64_t k, const uint64_t m,
    hasher_type const& hasher)  //
{
    return dispatch_order(hasher.mode(), [&](auto mode) {
        return compute_minimizer_for<decltype(mode)::value, Kmer>(kmer, kmer_rc, k, m, hasher);
    });
}

/* The minimizer under the order of `hasher`. Native, the default order, is inlined at its call
   sites (the lookup) exactly as with upstream's single hasher; the others go through a call. */
template <typename Kmer>
inline minimizer_info compute_minimizer(Kmer const& kmer, Kmer const& kmer_rc, const uint64_t k,
                                        const uint64_t m, hasher_type const& hasher)  //
{
    if (__builtin_expect(hasher.mode() == custom_order_mode::NATIVE, true)) {
        return compute_minimizer_for<custom_order_mode::NATIVE, Kmer>(kmer, kmer_rc, k, m, hasher);
    }
    return compute_minimizer_dispatched<Kmer>(kmer, kmer_rc, k, m, hasher);
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
    Plain, single-strand minimizer: argmin over the w = k-m+1 raw (not
    per-locus canonicalized) m-mer loci, leftmost on a tie -- the textbook
    minimizer scheme, and sshash's *regular parsing* modality (the 2022
    SSHash paper's original mechanism). build_configuration::regular_parsing
    forces this (instead of compute_minimizer/minimizer_iterator), paired with
    dictionary::lookup's two-probe fallback: construction keeps a kmer's
    minimizer for whichever single strand it was walked in, and Lookup first
    searches for g and, only if not found, for its reverse complement.

    Forward by construction (plain sliding-window argmin, leftmost tie-break
    never revisits an abandoned locus), so it needs no separate forwardness
    verification.
*/
template <typename Kmer>
minimizer_info compute_minimizer_single_strand(Kmer kmer, const uint64_t k, const uint64_t m,
                                               hasher_type const& hasher)  //
{
    assert(m <= Kmer::max_m);
    assert(m <= k);

    auto mmer_at = [&](Kmer const& window) {
        Kmer mmer = window;
        mmer.take_chars(m);
        return uint64_t(mmer);
    };

    Kmer window = kmer;
    uint64_t minimizer = mmer_at(window);
    uint64_t min_hash = hasher.hash(minimizer);
    uint64_t pos = 0;
    window.drop_char();

    for (uint64_t i = 1; i != k - m + 1; ++i) {
        uint64_t value = mmer_at(window);
        uint64_t hash = hasher.hash(value);
        if (hash < min_hash) {  // leftmost
            min_hash = hash;
            minimizer = value;
            pos = i;
        }
        window.drop_char();
    }

    return {minimizer, pos};
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
