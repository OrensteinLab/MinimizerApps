// Checks the GreedyMini rank/plane API and strand_invariant_minimizer_iterator (rolling Sigma
// plane, ring of window ranks, full hashes only on rank ties) against what it replaces:
//  (a) hash_from_planes(planes_of_mmer(x)) == hash_gm(x) == hash_reference(x) (the centre
//      reference is the independent per-base implementation), over random / biased /
//      palindromic / repeat m-mers;
//  (b) the Sigma plane of a k-mer (sigma_of_kmer, and rolled with roll_sigma along a sequence)
//      gives, via centre_rank, the rank (hash_reference(window) >> 2m) of every m-mer window;
//  (c) the iterator's minimizer_info for every k-mer equals that of the pre-optimization iterator
//      (copied verbatim below as legacy_si_iterator, which hashes each window with hash_gm), and
//      util::compute_minimizer (rank first) equals util::compute_minimizer_reference (every
//      window fully hashed), over random and low-complexity / repeat / palindromic sequences
//      with resets and set_position() jumps, for uint64 and __uint128_t k-mers.
//  (d) tie-heavy synthetic tables (many equal ranks; all ranks equal; no window with a UHS
//      member; all keys share one finite rank) through the same iterator comparison.
// Usage: test_gm_rolling DIR_WITH_w11_k{13,15}_canonical.gm [num_hash_samples_per_config=15000000]
//        [num_iterator_windows_per_config=2500000]
#include <atomic>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

#include "include/util.hpp"
#include "include/kmer.hpp"
#include "include/minimizer_iterator.hpp"

using namespace sshash;

namespace legacy {
template <typename Kmer>
struct legacy_si_iterator {
    legacy_si_iterator() {}

    legacy_si_iterator(uint64_t k, uint64_t m, hasher_type const& hasher,
                                        uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_min_value(constants::invalid_uint64)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher_p(&hasher)  //
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

    /* Not inlined: keeps the callers' hot paths (shared with the native orders) small. */
    __attribute__((noinline)) minimizer_info next(Kmer kmer) {
        assert(hasher().strand_invariant());
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
            Kmer window = kmer;
            window.drop_chars(m_k - m_m);
            window.take_chars(m_m);
            const uint64_t hash = hasher().hash_gm(uint64_t(window));
            if (hash < m_min_hash) {
                m_min_hash = hash;
                Kmer last = kmer;
                last.drop_chars(m_k - m_m);
                m_min_value = util::canonical_mmer_at<Kmer>(last, rc(), m_k, m_m, m_k - m_m);
                m_min_position = m_position;
                m_min_pos_in_kmer = m_k - m_m;
                m_num_mins = 1;
            } else {
                assert(m_min_pos_in_kmer > 0);
                m_min_pos_in_kmer -= 1;
                if (hash == m_min_hash) m_num_mins += 1;
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

    template <typename GetRc>
    void rescan(Kmer kmer, GetRc& rc) {
        const uint64_t begin = m_position;
        const Kmer whole = kmer;

        auto hash_of = [&](Kmer window) {
            window.take_chars(m_m);
            return hasher().hash_gm(uint64_t(window));
        };

        m_min_hash = hash_of(kmer);
        m_min_pos_in_kmer = 0;
        m_num_mins = 1;
        kmer.drop_char();
        for (uint64_t i = 1; i != m_k - m_m + 1; ++i) {
            const uint64_t hash = hash_of(kmer);
            kmer.drop_char();
            if (hash < m_min_hash) {  // leftmost
                m_min_hash = hash;
                m_min_pos_in_kmer = i;
                m_num_mins = 1;
            } else if (hash == m_min_hash) {
                m_num_mins += 1;
            }
        }

        Kmer window = whole;
        window.drop_chars(m_min_pos_in_kmer);
        m_min_value = util::canonical_mmer_at<Kmer>(window, rc(), m_k, m_m, m_min_pos_in_kmer);
        m_min_position = begin + m_min_pos_in_kmer;
        m_position = begin + m_k - m_m;
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


}  // namespace legacy

static std::atomic<uint64_t> failures{0};
static void fail(std::string const& what) {
    if (failures++ < 10) std::cerr << "FAILED: " << what << std::endl;
}

static uint64_t gen_mmer(std::mt19937_64& rng, uint64_t i, uint64_t l) {
    const uint64_t mask = (uint64_t(1) << (2 * l)) - 1;
    uint64_t x = rng() & mask;
    switch (i % 6) {
        case 1: x &= (rng() | rng()); break;
        case 2: x |= (rng() & rng()); break;
        case 3: {
            const uint64_t half = l / 2, hm = (uint64_t(1) << (2 * half)) - 1;
            uint64_t y = rng() & hm;
            x = ((y << (2 * half)) | dna_reverse_complement(y, (uint8_t)half)) & mask;
            break;
        }
        case 4: x = (rng() & 0xFull) * 0x1111111111111111ull & mask; break;
        default: break;
    }
    return x;
}

static std::string gen_sequence(std::mt19937_64& gen, uint64_t length) {
    std::string s;
    while (s.size() < length) {
        const uint64_t len = 1 + gen() % 60;
        const uint64_t kind = gen() % 8;
        std::string seg;
        if (kind <= 2) {
            for (uint64_t i = 0; i != len; ++i) seg += "ACGT"[gen() & 3];
        } else if (kind == 3) {  // AT-rich / CG-rich
            const char* al = (gen() & 1) ? "AT" : "CG";
            for (uint64_t i = 0; i != len; ++i) seg += al[gen() & 1];
        } else if (kind == 4) {  // homopolymer / short repeat
            std::string unit;
            for (uint64_t i = 0, u = 1 + gen() % 4; i != u; ++i) unit += "ACGT"[gen() & 3];
            while (seg.size() < len) seg += unit;
        } else if (kind == 5) {  // reverse-complement palindrome
            std::string y;
            for (uint64_t i = 0; i != len / 2 + 1; ++i) y += "ACGT"[gen() & 3];
            std::string r(y.rbegin(), y.rend());
            for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
            seg = y + r;
        } else if (s.size() > len) {  // copy of an earlier stretch, or its reverse complement
            seg = s.substr(gen() % (s.size() - len), len);
            if (gen() & 1) {
                std::string r(seg.rbegin(), seg.rend());
                for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
                seg = r;
            }
        }
        s += seg;
    }
    s.resize(length);
    return s;
}

template <typename kmer_t>
static void check_planes_and_windows(hasher_type const& h, uint64_t k, uint64_t m, uint64_t n_kmers,
                                     std::mt19937_64& gen, uint64_t& n_windows) {
    const std::string s = gen_sequence(gen, 4000);
    const uint64_t kmask = (uint64_t(1) << k) - 1;
    uint64_t rolled = 0;
    bool have = false;
    for (uint64_t it = 0; it != n_kmers; ++it) {
        const uint64_t pos = gen() % (s.size() - k);
        for (uint64_t j = 0; j != 40 and pos + j + k <= s.size(); ++j) {
            kmer_t kmer = util::string_to_uint_kmer<kmer_t>(s.data() + pos + j, k);
            uint64_t lo = uint64_t(kmer), hi = 0;
            if constexpr (sizeof(kmer.bits) > 8) hi = uint64_t(kmer.bits >> 64);
            const uint64_t fresh = h.sigma_of_kmer(lo, hi, (int)k);
            if (!have or j == 0) {
                rolled = fresh;
                have = true;
            } else {
                rolled = h.roll_sigma(rolled, kmer.at(k - 1), (int)k);
            }
            if ((rolled & kmask) != fresh) fail("rolled sigma != fresh sigma (k=" + std::to_string(k) + ")");
            for (uint64_t i = 0; i != k - m + 1; ++i) {
                kmer_t w = kmer;
                w.drop_chars(i);
                w.take_chars(m);
                const uint64_t ref = h.hash_reference(uint64_t(w));
                ++n_windows;
                if (h.centre_rank(fresh, (int)i) != (ref >> (2 * m)) or
                    h.centre_rank(rolled, (int)i) != (ref >> (2 * m))) {
                    fail("window rank from the Sigma plane != reference rank (k=" + std::to_string(k) +
                         " m=" + std::to_string(m) + " i=" + std::to_string(i) + ")");
                }
            }
        }
    }
}

template <typename kmer_t>
static uint64_t check_iterator(hasher_type const& h, uint64_t k, uint64_t m, uint64_t n_windows,
                               std::mt19937_64& gen, uint64_t& n_ties) {
    strand_invariant_minimizer_iterator<kmer_t> it(k, m, h);
    legacy::legacy_si_iterator<kmer_t> old_it(k, m, h);
    uint64_t checked = 0;
    while (checked < n_windows) {
        const std::string s = gen_sequence(gen, 600 + gen() % 3000);
        uint64_t position = gen() % 1000;
        it.set_position(position);
        old_it.set_position(position);
        for (uint64_t i = 0; i + k <= s.size();) {
            /* events: reset() and set_position() in the middle of a string, followed by a jump
               to another place of the string (the next k-mer is not a slide) */
            const uint64_t ev = gen() % 64;
            if (ev == 0) {
                it.reset();
                old_it.reset();
                i = gen() % (s.size() - k + 1);
            } else if (ev == 1) {
                position = gen() % 100000;
                it.set_position(position);
                old_it.set_position(position);
                i = gen() % (s.size() - k + 1);
            }
            kmer_t kmer = util::string_to_uint_kmer<kmer_t>(s.data() + i, k);
            minimizer_info a = it.next(kmer);
            minimizer_info b = old_it.next(kmer);
            const minimizer_info c = util::compute_minimizer<kmer_t>(kmer, k, m, h);
            kmer_t kmer_rc = kmer;
            kmer_rc.reverse_complement_inplace(k);
            const minimizer_info r = util::compute_minimizer_reference<kmer_t>(kmer, kmer_rc, k, m, h);
            if (c != r) {
                fail("rank-first compute_minimizer != all-windows reference at k=" +
                     std::to_string(k) + " m=" + std::to_string(m));
            }
            ++checked;
            if (a != b) {
                fail("rolling iterator != legacy iterator at k=" + std::to_string(k) + " m=" +
                     std::to_string(m));
            }
            if (a.minimizer != c.minimizer or a.pos_in_kmer != c.pos_in_kmer) {
                fail("rolling iterator != util::compute_minimizer at k=" + std::to_string(k) +
                     " m=" + std::to_string(m) + ": got (" + std::to_string(a.minimizer) + "," +
                     std::to_string(a.pos_in_kmer) + ") legacy (" + std::to_string(b.minimizer) +
                     "," + std::to_string(b.pos_in_kmer) + ") expected (" +
                     std::to_string(c.minimizer) + "," + std::to_string(c.pos_in_kmer) + ")");
            }
            ++i;
        }
        /* strings shorter than k: nothing is fed; the next string starts with set_position() */
    }
    (void)n_ties;
    return checked;
}

/* synthetic mirror-symmetric binary table of 2^gk keys. kind: 0 = 1/3 of the key pairs are UHS
   members with distinct ranks; 1 = every rank is the shared "infinity"; 2 = ~2% UHS members
   (most windows have none); 3 = all keys share the finite rank 0; 4 = UHS members share just
   3 distinct ranks (many ties among finite ranks) */
static void make_table(const char* path, int gk, int kind, uint64_t seed) {
    const uint64_t n = uint64_t(1) << gk;
    std::mt19937_64 rng(seed);
    std::vector<uint64_t> raw(n, n);
    auto mirror = [&](uint64_t key) {
        uint64_t r = 0;
        for (int i = 0; i != gk; ++i) r |= ((key >> i) & 1) << (gk - 1 - i);
        return r;
    };
    uint64_t next_rank = 0;
    for (uint64_t key = 0; key != n; ++key) {
        const uint64_t mir = mirror(key);
        if (mir < key) {
            raw[key] = raw[mir];
            continue;
        }
        switch (kind) {
            case 0: if (rng() % 3 == 0) raw[key] = next_rank++; break;
            case 1: break;
            case 2: if (rng() % 50 == 0) raw[key] = next_rank++; break;
            case 3: raw[key] = 0; break;
            default: if (rng() % 3 == 0) raw[key] = rng() % 3; break;
        }
        raw[mir] = raw[key];
    }
    write_gm_file_raw(path, raw);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " DIR [hash_samples] [iterator_windows]" << std::endl;
        return 1;
    }
    const std::string dir = argv[1];
    const uint64_t n_hash = argc > 2 ? std::stoull(argv[2]) : 15000000ull;
    const uint64_t n_iter = argc > 3 ? std::stoull(argv[3]) : 2500000ull;
    using kmer_t = dna_uint_kmer_t<uint64_t>;
    using wide_t = dna_uint_kmer_t<__uint128_t>;

    uint64_t total_hash = 0, total_windows = 0, total_iter = 0;
    for (const char* table : {"w11_k13_canonical.gm", "w11_k15_canonical.gm"})
            for (auto proj : {gm_projection::AT0_CG1, gm_projection::AT1_CG0}) {
                const std::string path = dir + "/" + table;
                hasher_type h(1234);
                const uint64_t m = 21;
                h.configure(custom_order_mode::GREEDYMINI, m, path, proj);
                const uint64_t before = failures;

                /* (a) */
                std::vector<std::thread> th;
                const unsigned nt = 16;
                for (unsigned t = 0; t != nt; ++t)
                    th.emplace_back([&, t]() {
                        std::mt19937_64 rng(1000 + t);
                        for (uint64_t i = 0; i != n_hash / nt; ++i) {
                            const uint64_t x = gen_mmer(rng, i, m);
                            const uint64_t a = h.hash_from_planes(h.planes_of_mmer(x));
                            const uint64_t b = h.hash_gm(x);
                            const uint64_t c = h.hash_reference(x);
                            if (a != b or b != c) fail("hash_from_planes/hash_gm/reference differ");
                        }
                    });
                for (auto& x : th) x.join();
                total_hash += (n_hash / nt) * nt;

                /* (b) + (c): several (k, Kmer) combinations, in parallel */
                std::atomic<uint64_t> windows{0}, iters{0};
                std::vector<std::thread> th2;
                unsigned id = 0;
                auto spawn = [&](auto kmer_tag, uint64_t k) {
                    using K = decltype(kmer_tag);
                    th2.emplace_back([&, k, seed = ++id]() {
                        std::mt19937_64 gen(77 * seed + 5);
                        uint64_t w = 0, ties = 0;
                        check_planes_and_windows<K>(h, k, m, 40000 / (k - m + 1) + 50, gen, w);
                        windows += w;
                        iters += check_iterator<K>(h, k, m, n_iter / 6, gen, ties);
                    });
                };
                spawn(kmer_t{}, 21); spawn(kmer_t{}, 25); spawn(kmer_t{}, 31);
                spawn(wide_t{}, 31); spawn(wide_t{}, 47); spawn(wide_t{}, 63);
                for (auto& x : th2) x.join();
                total_windows += windows;
                total_iter += iters;
                std::cout << table << " " << gm_projection_name(proj)
                          << ": " << (n_hash / nt) * nt << " m-mers (planes==gm==reference), "
                          << windows << " rank windows, " << iters << " iterator k-mers"
                          << (failures == before ? " ok" : "  FAILED") << std::endl;
            }
    {   /* (d) synthetic tables, m = gk + 4, strings of every k and both kmer widths */
        const char* tmp = "/tmp/test_gm_rolling_synth.gm";
        const char* names[] = {"1/3 UHS", "all INF", "2% UHS", "all rank 0", "3 shared ranks"};
        for (int kind = 0; kind != 5; ++kind)
            for (auto proj : {gm_projection::AT0_CG1, gm_projection::AT1_CG0}) {
                make_table(tmp, 9, kind, 100 + kind);
                hasher_type h(777);
                const uint64_t m = 13;
                h.configure(custom_order_mode::GREEDYMINI, m, tmp, proj);
                const uint64_t before = failures;
                std::atomic<uint64_t> iters{0};
                std::vector<std::thread> th3;
                unsigned id = 0;
                auto spawn = [&](auto kmer_tag, uint64_t k) {
                    using K = decltype(kmer_tag);
                    th3.emplace_back([&, k, seed = ++id]() {
                        std::mt19937_64 gen(991 * seed + kind);
                        uint64_t ties = 0;
                        iters += check_iterator<K>(h, k, m, n_iter / 3, gen, ties);
                    });
                };
                spawn(kmer_t{}, 15); spawn(kmer_t{}, 21); spawn(kmer_t{}, 31);
                spawn(wide_t{}, 31); spawn(wide_t{}, 47); spawn(wide_t{}, 63);
                for (auto& x : th3) x.join();
                total_iter += iters;
                std::cout << "synthetic table '" << names[kind] << "' " << gm_projection_name(proj)
                          << ": " << iters << " iterator k-mers" << (failures == before ? " ok" : "  FAILED")
                          << std::endl;
            }
        std::remove(tmp);
    }
    std::cout << "total: " << total_hash << " m-mers, " << total_windows << " rank windows, "
              << total_iter << " iterator k-mers" << std::endl;
    if (failures) {
        std::cerr << failures << " FAILURES" << std::endl;
        return 1;
    }
    std::cout << "OK" << std::endl;
    return 0;
}
