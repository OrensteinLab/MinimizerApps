#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "include/util.hpp"
#include "include/kmer.hpp"
#include "include/minimizer_iterator.hpp"

using namespace sshash;

static uint64_t num_failed = 0;
static uint64_t num_checked = 0;
static uint64_t num_ties = 0;

static void fail(std::string const& what) {
    if (num_failed < 10) std::cerr << "FAILED: " << what << std::endl;
    ++num_failed;
}

/*
    The property the whole design rests on: the minimizer must be
    mirror-equivariant. The m-mer selected for a kmer and for its reverse
    complement must be literally the same string -- otherwise the two would land
    in different buckets and a lookup by one orientation could not find the other
    -- and the selected locus must be the mirror image.
*/
template <typename kmer_t>
static void check_equivariance(kmer_t kmer, uint64_t k, uint64_t m, hasher_type const& hasher) {
    kmer_t kmer_rc = kmer;
    kmer_rc.reverse_complement_inplace(k);

    auto a = util::compute_minimizer(kmer, k, m, hasher);
    auto b = util::compute_minimizer(kmer_rc, k, m, hasher);

    ++num_checked;

    if (a.minimizer != b.minimizer) {
        fail("minimizer of kmer '" + util::uint_kmer_to_string(kmer, k) + "' (m=" +
             std::to_string(m) + ") is " + util::uint_minimizer_to_string<kmer_t>(a.minimizer, m) +
             " but that of its reverse complement is " +
             util::uint_minimizer_to_string<kmer_t>(b.minimizer, m));
        return;
    }
    if (a.pos_in_kmer + b.pos_in_kmer != k - m) {
        fail("locus of kmer '" + util::uint_kmer_to_string(kmer, k) + "' is " +
             std::to_string(a.pos_in_kmer) + " but that of its reverse complement is " +
             std::to_string(b.pos_in_kmer) + " (should mirror to " +
             std::to_string(k - m - a.pos_in_kmer) + ")");
        return;
    }

    /* the reported value must really be the canonical m-mer at the reported locus */
    kmer_t anchored = kmer;
    anchored.drop_chars(a.pos_in_kmer);
    anchored.take_chars(m);
    if (a.minimizer != util::canonical_mmer<kmer_t>(uint64_t(anchored), m)) {
        fail("minimizer value does not match the m-mer at the selected locus");
    }
}

/* Count the loci attaining the minimum hash, to report how often ties fire. */
template <typename kmer_t>
static bool has_tie(kmer_t kmer, uint64_t k, uint64_t m, hasher_type const& hasher) {
    uint64_t min_hash = constants::invalid_uint64;
    uint64_t count = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mmer = kmer;
        mmer.take_chars(m);
        uint64_t hash = hasher.hash(util::canonical_mmer<kmer_t>(uint64_t(mmer), m));
        if (i == 0 or hash < min_hash) {
            min_hash = hash;
            count = 1;
        } else if (hash == min_hash) {
            ++count;
        }
        kmer.drop_char();
    }
    return count > 1;
}

template <typename kmer_t>
static void test_kmers(uint64_t k, uint64_t m, uint64_t num_kmers, std::mt19937_64& gen,
                       custom_order_mode mode) {
    hasher_type hasher(1234);
    hasher.configure(mode, m);
    std::string s(k, 0);
    for (uint64_t n = 0; n != num_kmers; ++n) {
        for (uint64_t i = 0; i != k; ++i) s[i] = "ACGT"[gen() & 3];
        kmer_t kmer = util::string_to_uint_kmer<kmer_t>(s.data(), k);
        if (has_tie(kmer, k, m, hasher)) ++num_ties;
        check_equivariance(kmer, k, m, hasher);
    }
}

/*
    The incremental "re-scan" iterator must agree with the brute-force reference
    on every window of a sequence, including the windows where a tie fires.

    The scheme must also be *forward*: the sampled position must never decrease
    as the window slides, so that the number of super-kmers equals the number of
    sampled positions. The centre-closest tie-break is what guarantees this
    (leftmost-if-canonical / rightmost-otherwise, the previous rule, is
    mirror-equivariant too but moves the anchor backwards on ~1e-5 of the
    windows).
*/
template <typename kmer_t, custom_order_mode Mode>
static void test_iterator(uint64_t k, uint64_t m, uint64_t length, std::mt19937_64& gen) {
    hasher_type hasher(1234);
    hasher.configure(Mode, m);
    std::string s(length, 0);
    for (uint64_t i = 0; i != length; ++i) s[i] = "ACGT"[gen() & 3];

    minimizer_iterator<kmer_t, Mode> it(k, m, hasher);
    it.set_position(0);

    uint64_t prev_pos_in_seq = 0;
    for (uint64_t i = 0; i + k <= length; ++i) {
        kmer_t kmer = util::string_to_uint_kmer<kmer_t>(s.data() + i, k);
        auto got = it.next(kmer);
        auto expected = util::compute_minimizer(kmer, k, m, hasher);
        ++num_checked;
        if (got.minimizer != expected.minimizer or got.pos_in_kmer != expected.pos_in_kmer) {
            fail("iterator disagrees with the reference at position " + std::to_string(i) +
                 " (k=" + std::to_string(k) + ", m=" + std::to_string(m) + ")");
            return;
        }
        /* pos_in_seq must be the absolute position of the selected locus */
        if (got.pos_in_seq != i + got.pos_in_kmer) {
            fail("iterator reports pos_in_seq " + std::to_string(got.pos_in_seq) +
                 " but expected " + std::to_string(i + got.pos_in_kmer));
            return;
        }
        /* forwardness: the sampled position never decreases */
        if (got.pos_in_seq < prev_pos_in_seq) {
            fail("scheme is not forward: sampled position " + std::to_string(got.pos_in_seq) +
                 " after " + std::to_string(prev_pos_in_seq) + " (k=" + std::to_string(k) +
                 ", m=" + std::to_string(m) + ")");
            return;
        }
        prev_pos_in_seq = got.pos_in_seq;
    }
}

/*
    GreedyMini (strand-invariant order): the iterator and util::compute_minimizer hash the RAW
    m-mers and canonicalize lazily. Check both against an independent brute-force reference that
    follows the original definition: hash the numerically-canonical m-mer with the *reference*
    dual-orientation hash, leftmost minimum, ties broken centre-closest (mirror images: smaller
    index if kmer <= rc(kmer)). Sequences are random with injected tandem repeats and
    reverse-complement copies, so that ties (equal hashes at several loci) actually fire.
*/
template <typename kmer_t>
static minimizer_info gm_reference_minimizer(kmer_t kmer, uint64_t k, uint64_t m,
                                             hasher_type const& hasher) {
    kmer_t kmer_rc = kmer;
    kmer_rc.reverse_complement_inplace(k);
    std::vector<uint64_t> h(k - m + 1), val(k - m + 1);
    kmer_t w = kmer;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mm = w;
        mm.take_chars(m);
        val[i] = util::canonical_mmer<kmer_t>(uint64_t(mm), m);
        h[i] = hasher.hash_reference(val[i]);
        w.drop_char();
    }
    uint64_t best = 0;
    for (uint64_t i = 1; i != k - m + 1; ++i) {
        if (h[i] < h[best]) best = i;
    }
    const uint64_t two_c = k - m;
    uint64_t best_dist = constants::invalid_uint64, lo = best, hi = best;
    for (uint64_t i = best; i != k - m + 1; ++i) {
        if (h[i] != h[best]) continue;
        const uint64_t dist = 2 * i > two_c ? 2 * i - two_c : two_c - 2 * i;
        if (dist < best_dist) {
            best_dist = dist;
            lo = hi = i;
        } else if (dist == best_dist) {
            hi = i;
        }
    }
    const uint64_t chosen = (lo == hi or !(kmer_rc < kmer)) ? lo : hi;
    return {val[best], chosen};
}

template <typename kmer_t>
static void test_iterator_gm(uint64_t k, uint64_t m, uint64_t length, std::mt19937_64& gen,
                             hasher_type const& hasher, uint64_t& num_tie_windows) {
    std::string s;
    while (s.size() < length) {
        const uint64_t kind = gen() % 4;
        const uint64_t len = 20 + gen() % 80;
        std::string seg;
        if (kind == 0) {  // random
            for (uint64_t i = 0; i != len; ++i) seg += "ACGT"[gen() & 3];
        } else if (kind == 1) {  // tandem repeat with a short period
            std::string unit;
            for (uint64_t i = 0, p = 1 + gen() % 6; i != p; ++i) unit += "ACGT"[gen() & 3];
            for (uint64_t i = 0; i != len; ++i) seg += unit[i % unit.size()];
        } else if (kind == 2) {  // palindromic hairpin: y + rc(y)
            std::string y;
            for (uint64_t i = 0; i != len / 2; ++i) y += "ACGT"[gen() & 3];
            std::string r(y.rbegin(), y.rend());
            for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
            seg = y + r;
        } else {  // copy of an earlier stretch, or its reverse complement
            if (s.size() > len) {
                seg = s.substr(gen() % (s.size() - len), len);
                if (gen() & 1) {
                    std::string r(seg.rbegin(), seg.rend());
                    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
                    seg = r;
                }
            }
        }
        s += seg;
    }
    s.resize(length);

    strand_invariant_minimizer_iterator<kmer_t> it(k, m, hasher);
    it.set_position(0);
    for (uint64_t i = 0; i + k <= length; ++i) {
        kmer_t kmer = util::string_to_uint_kmer<kmer_t>(s.data() + i, k);
        auto got = it.next(kmer);
        auto direct = util::compute_minimizer(kmer, k, m, hasher);
        auto expected = gm_reference_minimizer<kmer_t>(kmer, k, m, hasher);
        ++num_checked;
        if (has_tie(kmer, k, m, hasher)) ++num_tie_windows;
        if (got.minimizer != expected.minimizer or got.pos_in_kmer != expected.pos_in_kmer or
            direct.minimizer != expected.minimizer or direct.pos_in_kmer != expected.pos_in_kmer) {
            fail("GreedyMini iterator/compute_minimizer disagree with the reference at position " +
                 std::to_string(i) + " (k=" + std::to_string(k) + ", m=" + std::to_string(m) +
                 ")");
            return;
        }
    }
}

int main(int argc, char** argv) {
    using kmer_t = dna_uint_kmer_t<uint64_t>;
    using wide_kmer_t = dna_uint_kmer_t<__uint128_t>;

    std::mt19937_64 gen(42);

    /* the strand-dependent orders: native, xor, xor-pcg-hash */
    const custom_order_mode plain_modes[] = {custom_order_mode::NATIVE, custom_order_mode::XOR_HASH,
                                             custom_order_mode::XOR_PCG_HASH};

    std::cout << "checking mirror-equivariance of the minimizer..." << std::endl;
    for (auto mode : plain_modes) {
        for (uint64_t k : {5, 15, 21, 31}) {
            for (uint64_t m = 1; m <= std::min<uint64_t>(k, kmer_t::max_m); ++m) {
                test_kmers<kmer_t>(k, m, 20000, gen, mode);
            }
        }
        for (uint64_t k : {31, 47, 63}) {
            for (uint64_t m : {1, 2, 3, 7, 13, 21, 31}) {
                if (m > k) continue;
                test_kmers<wide_kmer_t>(k, m, 20000, gen, mode);
            }
        }
    }

    /*
        Small m makes the m-mer universe tiny, so a window very often contains
        two loci carrying the same canonical m-mer: this is what stresses the
        tie-breaking rule, which on real data with m=13 fires on ~1e-5 of the
        windows only.
    */
    std::cout << "  " << num_ties << "/" << num_checked << " of the kmers checked had a tie"
              << std::endl;

    std::cout << "checking the incremental iterator against the reference..." << std::endl;
    auto check_iterators = [&](auto mode) {
        constexpr custom_order_mode Mode = decltype(mode)::value;
        for (uint64_t k : {5, 15, 31}) {
            for (uint64_t m = 1; m <= std::min<uint64_t>(k, kmer_t::max_m); ++m) {
                test_iterator<kmer_t, Mode>(k, m, 20000, gen);
            }
        }
        for (uint64_t k : {47, 63}) {
            for (uint64_t m : {2, 3, 13, 31}) { test_iterator<wide_kmer_t, Mode>(k, m, 20000, gen); }
        }
    };
    check_iterators(order_tag<custom_order_mode::NATIVE>());
    check_iterators(order_tag<custom_order_mode::XOR_HASH>());
    check_iterators(order_tag<custom_order_mode::XOR_PCG_HASH>());

    /* optional: GreedyMini tables given on the command line (strand-invariant raw path) */
    for (int a = 1; a < argc; ++a) {
        for (auto proj : {gm_projection::AT0_CG1, gm_projection::AT1_CG0}) {
            uint64_t before = num_checked, tie_windows = 0;
            for (uint64_t m : {13, 14, 15, 16, 17, 21, 23}) {
                hasher_type hasher(1234);
                try {
                    hasher.configure(custom_order_mode::GREEDYMINI, m, argv[a],
                                     proj);
                } catch (std::exception const&) { continue; }  // gm_k > m, or m - gm_k odd
                for (uint64_t k : {21, 25, 31}) {
                    if (m > k) continue;
                    test_iterator_gm<kmer_t>(k, m, 60000, gen, hasher, tie_windows);
                }
                for (uint64_t k : {31, 47, 63}) {
                    test_iterator_gm<wide_kmer_t>(k, m, 60000, gen, hasher, tie_windows);
                }
            }
            std::cout << "GreedyMini " << argv[a] << " " << gm_projection_name(proj)  << ": "
                      << (num_checked - before) << " windows match the reference (" << tie_windows
                      << " with a hash tie)" << std::endl;
        }
    }

    std::cout << "checked " << num_checked << " kmers" << std::endl;
    if (num_failed != 0) {
        std::cerr << num_failed << " CHECKS FAILED" << std::endl;
        return 1;
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return 0;
}
