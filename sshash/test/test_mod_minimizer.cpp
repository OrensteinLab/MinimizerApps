// Tests of the strand-symmetric mod-minimizer schemes (mod-minimizer-xor, mod-minimizer-greedymini):
// util::compute_mod_minimizer against a string-level reference of the definition, the
// mirror-equivariance a single-probe lookup needs, and mod_minimizer_iterator against
// util::compute_minimizer over whole sequences (with resets).
//
// Usage: test_mod_minimizer [TMER_TABLE.gm ...]  (GreedyMini tables of t-mers, e.g.
//        orders/GreedyMini/w22_k10_canonical.gm; each is tested at t = gm_k and t = gm_k + 2,
//        with the (k, m) below for which t == m (mod k - m + 1))

#include <algorithm>
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
static uint64_t num_reselected = 0;  // sampled positions re-selected after being abandoned

static void fail(std::string const& what) {
    if (num_failed < 10) std::cerr << "FAILED: " << what << std::endl;
    ++num_failed;
}

static std::string revcomp(std::string const& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}

/* random DNA; with `low` set, from {A, C} only, so that equal t-mers (ties) are frequent */
static std::string random_dna(uint64_t n, std::mt19937_64& gen, bool low) {
    static const char bases[] = "ACGT";
    std::string s(n, 'A');
    for (char& c : s) c = bases[low ? gen() % 2 : gen() % 4];
    return s;
}

/*
    The definition, from the k-mer's string (see util::compute_mod_minimizer): hash every t-mer
    (xor: the numerically canonical packed t-mer XOR the seed; GreedyMini: the per-base
    specification hash_reference of the raw t-mer); every minimal locus x offers the offset
    x mod w; the offset closest to (w-1)/2 wins, two equally close ones the smaller if the
    packed kmer <= its packed reverse complement, else the larger; value the canonical packed
    m-mer at the offset.
*/
template <typename kmer_t>
static minimizer_info reference(std::string const& s, uint64_t k, uint64_t m, uint64_t t,
                                hasher_type const& hasher, uint64_t seed) {
    const uint64_t w = k - m + 1, T = k - t + 1;
    auto packed = [](std::string const& x) {
        return uint64_t(util::string_to_uint_kmer<kmer_t>(x.c_str(), x.size()));
    };
    auto canon = [&](std::string const& x) {
        const uint64_t v = packed(x);
        return std::min(v, kmer_t::reverse_complement_mmer(v, x.size()));
    };
    std::vector<uint64_t> h(T);
    for (uint64_t x = 0; x != T; ++x) {
        const std::string tmer = s.substr(x, t);
        h[x] = is_gm_mode(hasher.mode()) ? hasher.hash_reference(packed(tmer)) : canon(tmer) ^ seed;
    }
    const uint64_t mn = *std::min_element(h.begin(), h.end());
    std::vector<uint64_t> tied;
    for (uint64_t x = 0; x != T; ++x) {
        if (h[x] == mn) tied.push_back(x);
    }
    std::vector<uint64_t> offsets;
    for (uint64_t x : tied) offsets.push_back(x % w);
    auto dist = [&](uint64_t p) { return std::max(2 * p, w - 1) - std::min(2 * p, w - 1); };
    uint64_t best = dist(offsets[0]);
    for (uint64_t p : offsets) best = std::min(best, dist(p));
    std::vector<uint64_t> closest;
    for (uint64_t p : offsets) {
        if (dist(p) == best) closest.push_back(p);
    }
    std::sort(closest.begin(), closest.end());
    const kmer_t a = util::string_to_uint_kmer<kmer_t>(s.c_str(), k);
    const kmer_t b = util::string_to_uint_kmer<kmer_t>(revcomp(s).c_str(), k);
    const uint64_t p = !(b < a) ? closest.front() : closest.back();
    return {canon(s.substr(p, m)), p};
}

template <typename kmer_t>
static void test_config(uint64_t k, uint64_t m, hasher_type const& hasher, uint64_t seed,
                        std::mt19937_64& gen, uint64_t num_kmers) {
    const uint64_t t = hasher.mod_t();
    const std::string label = std::string(custom_scheme_mode_name(hasher.mode())) +
                              " k=" + std::to_string(k) + " m=" + std::to_string(m) +
                              " t=" + std::to_string(t);
    uint64_t before = num_failed;
    const uint64_t reselected_before = num_reselected;

    /* 1. compute_minimizer == the definition, and mirror-equivariance */
    for (uint64_t i = 0; i != num_kmers; ++i) {
        const std::string s = random_dna(k, gen, i % 4 == 0);
        const kmer_t kmer = util::string_to_uint_kmer<kmer_t>(s.c_str(), k);
        const kmer_t kmer_rc = util::string_to_uint_kmer<kmer_t>(revcomp(s).c_str(), k);
        const minimizer_info got = util::compute_minimizer(kmer, kmer_rc, k, m, hasher);
        const minimizer_info exp = reference<kmer_t>(s, k, m, t, hasher, seed);
        ++num_checked;
        if (!(got == exp)) {
            fail(label + ": kmer " + s + " samples locus " + std::to_string(got.pos_in_kmer) +
                 ", the definition " + std::to_string(exp.pos_in_kmer));
            break;
        }
        const minimizer_info rc = util::compute_minimizer(kmer_rc, kmer, k, m, hasher);
        if (s == revcomp(s)) continue;  // both strands take the left offset (even k only)
        if (rc.minimizer != got.minimizer or rc.pos_in_kmer + got.pos_in_kmer != k - m) {
            fail(label + ": kmer " + s + " and its reverse complement sample loci " +
                 std::to_string(got.pos_in_kmer) + " and " + std::to_string(rc.pos_in_kmer));
            break;
        }
    }

    /* 2. the incremental iterator == compute_minimizer, over sequences, with resets */
    dispatch_scheme(hasher.mode(), [&](auto mode) {
        constexpr custom_scheme_mode Mode = decltype(mode)::value;
        if constexpr (is_mod_mode(Mode)) {
            for (uint64_t r = 0; r != 4 and num_failed == before; ++r) {
                const std::string seq = random_dna(4000, gen, r == 3);
                const uint64_t base = 1000 * r;
                mod_minimizer_iterator<kmer_t, Mode> it(k, m, hasher, base);
                std::vector<uint64_t> abandoned;  // sampled positions left behind, this run
                uint64_t prev = constants::invalid_uint64;
                for (uint64_t j = 0; j + k <= seq.size(); ++j) {
                    if (j == 1500 or j == 2500) {
                        abandoned.clear();
                        prev = constants::invalid_uint64;
                    }
                    if (j == 1500) it.set_position(base + j);  // a new string
                    if (j == 2500) it.reset();                 // a stream restart
                    const kmer_t kmer = util::string_to_uint_kmer<kmer_t>(seq.c_str() + j, k);
                    const minimizer_info got = it.next(kmer);
                    const minimizer_info exp = util::compute_minimizer(kmer, k, m, hasher);
                    ++num_checked;
                    if (got.minimizer != exp.minimizer or got.pos_in_kmer != exp.pos_in_kmer or
                        got.pos_in_seq != base + j + exp.pos_in_kmer) {
                        fail(label + ": iterator differs at kmer " + std::to_string(j));
                        break;
                    }
                    /* forwardness: a position, once abandoned, is not sampled again */
                    if (got.pos_in_seq != prev) {
                        if (std::find(abandoned.begin(), abandoned.end(), got.pos_in_seq) !=
                            abandoned.end()) {
                            ++num_reselected;
                        }
                        if (prev != constants::invalid_uint64) abandoned.push_back(prev);
                        prev = got.pos_in_seq;
                    }
                }
            }
        }
    });
    const uint64_t reselected = num_reselected - reselected_before;
    std::cout << "  " << label << (num_failed == before ? ": ok" : ": FAILED")
              << ", re-selected positions: " << reselected << std::endl;
    /* t == m (mod w), on i.i.d. and {A,C}-only sequences: forward */
    if (reselected != 0) fail(label + ": a sampled position was re-selected");
}

int main(int argc, char** argv) {
    using kmer_t = dna_uint_kmer_t<uint64_t>;
    using wide_kmer_t = dna_uint_kmer_t<__uint128_t>;
    std::mt19937_64 gen(7);
    const uint64_t seed = 0x2545F4914F6CDD1Dull;

    std::cout << "mod-minimizer-xor:" << std::endl;
    struct cfg { uint64_t k, m, t; };
    /* t == m (mod k - m + 1), which the builder requires (m <= 32: the reference packs m-mers
       into 64 bits) */
    for (cfg c : {cfg{31, 21, 10}, cfg{31, 20, 8}, cfg{30, 20, 9}, cfg{15, 11, 6},
                  cfg{15, 11, 1}, cfg{21, 13, 4}, cfg{12, 8, 3}}) {
        hasher_type hasher(seed);
        hasher.configure(custom_scheme_mode::MOD_XOR_HASH, c.m, "", gm_projection::AT0_CG1, c.t);
        test_config<kmer_t>(c.k, c.m, hasher, seed, gen, 20000);
    }
    for (cfg c : {cfg{47, 31, 14}, cfg{55, 31, 6}, cfg{40, 30, 8}}) {
        hasher_type hasher(seed);
        hasher.configure(custom_scheme_mode::MOD_XOR_HASH, c.m, "", gm_projection::AT0_CG1, c.t);
        test_config<wide_kmer_t>(c.k, c.m, hasher, seed, gen, 20000);
    }

    /* configuration errors */
    {
        hasher_type h(seed);
        auto throws = [&](auto f) {
            try { f(); } catch (std::exception const&) { return true; }
            return false;
        };
        if (!throws([&] { h.configure(custom_scheme_mode::MOD_XOR_HASH, 21, "", gm_projection::AT0_CG1, 0); }))
            fail("mod-minimizer-xor without t accepted");
        if (!throws([&] { h.configure(custom_scheme_mode::MOD_XOR_HASH, 21, "", gm_projection::AT0_CG1, 21); }))
            fail("mod-minimizer-xor with t = m accepted");
        if (!throws([&] { h.configure(custom_scheme_mode::XOR_HASH, 21, "", gm_projection::AT0_CG1, 10); }))
            fail("--mod-t with a plain order accepted");
    }

    for (int a = 1; a < argc; ++a) {
        std::cout << "mod-minimizer-greedymini " << argv[a] << ":" << std::endl;
        for (auto proj : {gm_projection::AT0_CG1, gm_projection::AT1_CG0}) {
            for (uint64_t extra : {0, 2}) {
                const uint64_t t = GreedyMiniTable(load_gm_file_raw(argv[a])).gm_k() + extra;
                for (cfg c : {cfg{31, 21, t}, cfg{29, 20, t}, cfg{25, 19, t}, cfg{31, 22, t}}) {
                    if (c.t >= c.m or c.t % (c.k - c.m + 1) != c.m % (c.k - c.m + 1)) continue;
                    hasher_type hasher(seed);
                    hasher.configure(custom_scheme_mode::MOD_GREEDYMINI, c.m, argv[a], proj, c.t);
                    test_config<kmer_t>(c.k, c.m, hasher, seed, gen, 20000);
                }
            }
        }
    }

    std::cout << num_checked << " checks, " << num_failed << " failed" << std::endl;
    return num_failed ? 1 : 0;
}
