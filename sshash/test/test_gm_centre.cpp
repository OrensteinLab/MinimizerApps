// Checks the GreedyMini hash (-O greedymini): an explicit canonical
// gm_k table is extended to the m-mer xyz (|x| = |z| = (m - gm_k) / 2, |y| = gm_k) by ranking
// the binary projection of y (one table probe), breaking ties by min(Sigma(xyz),
// reverse(Sigma(xyz))) and then by Gamma.
//  (a) fast hash() == hash_reference() (per-base loop) == an independent std::string
//      implementation, over random / AT- and CG-rich / palindromic / repeat / shared-centre
//      m-mers;
//  (b) hash() == precomputed vectors of the Python reference (gm_centre_vectors.inc);
//  (c) hash(x) == hash(rc(x)); injectivity on canonical m-mers (and hash tie => same m-mer
//      up to reverse complement even for raw m-mers sharing a centre);
//  (d) rejection: m - gm_k odd (a 2^14-entry table at m = 21), asymmetric table, non-greedymini with a projection.
// Usage: test_gm_centre DIR_WITH_GM_TABLES [num_random_samples_per_config=30000000]
#include <atomic>
#include <cstdio>
#include <iostream>
#include <random>
#include <thread>
#include <unordered_map>
#include <vector>

#include "include/util.hpp"
#include "include/kmer.hpp"

using namespace sshash;

static std::atomic<uint64_t> failures{0};

static void report(const char* what, uint64_t x, uint64_t a, uint64_t b) {
    if (failures++ < 10) {
        std::cerr << what << " failed for x=" << x << " got " << a << " expected " << b
                  << std::endl;
    }
}

/* Reverse complement of an m-mer packed in sshash's encoding (A=0,C=1,G=3,T=2: the complement
   of a base is b ^ 2), base by base. */
static uint64_t rc_mmer(uint64_t x, uint64_t l) {
    uint64_t rc = 0;
    for (uint64_t p = 0; p != l; ++p) rc |= (((x >> (2 * p)) & 3) ^ 2) << (2 * (l - 1 - p));
    return rc;
}

struct string_ref {
    std::vector<uint16_t> ranks;
    int l, gk;
    bool flip;
    uint64_t seed;

    /* fully independent: decode to an ACGT string, work on strings */
    uint64_t hash(uint64_t x) const {
        std::string s(l, 'A');
        for (int p = 0; p != l; ++p) s[p] = nucleotides[(x >> (2 * (l - 1 - p))) & 3];
        std::string rc(s.rbegin(), s.rend());
        for (auto& c : rc) c = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C';
        std::string best;
        uint64_t best_rank = 0;
        for (std::string const* v : {&s, &rc}) {
            std::string sig(l, '0'), gam(l, '0');
            for (int p = 0; p != l; ++p) {
                const char c = (*v)[p];
                bool sg = (c == 'C' or c == 'G');
                if (flip) sg = !sg;
                sig[p] = sg ? '1' : '0';
                bool gm = (c == 'G' or c == 'T');
                if ((seed >> (l - 1 - p)) & 1) gm = !gm;
                gam[p] = gm ? '1' : '0';
            }
            const int i = (l - gk) / 2;
            const uint64_t key = std::stoull(sig.substr(i, gk), nullptr, 2);
            const std::string f = sig + gam;  // lexicographic order of the pair (Sigma, Gamma)
            if (best.empty() or f < best) best = f;
            best_rank = ranks[key];  // equal for both orientations (mirror symmetry)
        }
        return (best_rank << (2 * l)) | std::stoull(best, nullptr, 2);
    }
};

static inline void check_one(hasher_type const& h, string_ref const* sr, uint64_t x, uint64_t l) {
    const uint64_t got = h.hash(x);
    const uint64_t ref = h.hash_reference(x);
    if (got != ref) report("hash != hash_reference", x, got, ref);
    if (sr) {
        const uint64_t s = sr->hash(x);
        if (got != s) report("hash != string reference", x, got, s);
    }
    const uint64_t rc = rc_mmer(x, l);
    const uint64_t got_rc = h.hash(rc);
    if (got_rc != got) report("rc-invariance", x, got_rc, got);
}

static uint64_t gen_sample(std::mt19937_64& rng, uint64_t i, uint64_t l, uint64_t mask) {
    uint64_t x = rng() & mask;
    switch (i % 10) {
        case 1: x &= (rng() | rng()); break;  // AT-rich
        case 2: x |= (rng() & rng()); break;  // CG-rich
        case 3: {                             // palindrome (even l) / near palindrome
            const uint64_t half = l / 2, hm = (uint64_t(1) << (2 * half)) - 1;
            uint64_t y = rng() & hm;
            x = ((y << (2 * half)) | rc_mmer(y, half)) & mask;
            if (l & 1) x = (x & ~(uint64_t(3) << (2 * half))) | ((rng() & 3) << (2 * half));
            break;
        }
        case 4: x = (rng() & 0xFull) * 0x1111111111111111ull & mask; break;  // short repeat
        case 5: x = (x & ~(mask >> 10)) | (rc_mmer(x, l) & (mask >> 10)); break;
        case 6: {  // Sigma string is a palindrome (Sigma(x) == reverse Sigma(x)) but x is not
            for (uint64_t p = 0; p != l / 2; ++p) {
                const uint64_t sp = 2 * (l - 1 - p), sq = 2 * p;  // bit offsets of base p and l-1-p
                const uint64_t bp = (x >> sp) & 3, bq = (x >> sq) & 3;
                if ((bp & 1) != (bq & 1)) x ^= (uint64_t(1) << sq);  // bit 0 of a base is Sigma
            }
            break;
        }
        case 7: x = (x & 0x5555555555555555ull) | ((rng() & 1) ? 0 : (x & 0xAAAAAAAAAAAAAAAAull)); x &= mask; break;
        default: break;
    }
    return x;
}

static uint64_t run_checks(hasher_type const& h, string_ref const* sr, uint64_t l,
                           uint64_t num_samples, uint64_t seed) {
    const uint64_t mask = (uint64_t(1) << (2 * l)) - 1;
    const unsigned nt = 16;
    std::vector<std::thread> th;
    for (unsigned t = 0; t != nt; ++t) {
        th.emplace_back([&, t]() {
            std::mt19937_64 rng(seed * 1000 + t);
            const uint64_t n = num_samples / nt;
            for (uint64_t i = 0; i != n; ++i) check_one(h, sr, gen_sample(rng, i, l, mask), l);
        });
    }
    for (auto& t : th) t.join();
    uint64_t total = (num_samples / nt) * nt;
    if (2 * l <= 24) {  // exhaustive
        const uint64_t end = uint64_t(1) << (2 * l);
        std::vector<std::thread> th2;
        for (unsigned t = 0; t != nt; ++t) {
            th2.emplace_back([&, t]() {
                for (uint64_t x = t; x < end; x += nt) check_one(h, sr, x, l);
            });
        }
        for (auto& t : th2) t.join();
        total += end;
    }
    return total;
}

/* hash tie => same m-mer up to reverse complement, over canonical samples AND over raw
   samples that share their centre gm_k-mer (so ranks tie constantly). */
static void run_injectivity(hasher_type const& h, uint64_t l) {
    const uint64_t num_samples = 4000000;
    std::mt19937_64 rng(42);
    const uint64_t mask = (uint64_t(1) << (2 * l)) - 1;
    const int gk = h.gm_k(), i0 = (int)(l - gk) / 2;
    std::unordered_map<uint64_t, uint64_t> seen;
    seen.reserve(num_samples * 2);
    const uint64_t centre_mask = ((uint64_t(1) << (2 * gk)) - 1) << (2 * i0);
    uint64_t fixed_centre = rng() & centre_mask;
    for (uint64_t i = 0; i != num_samples; ++i) {
        uint64_t x = rng() & mask;
        if (i & 1) x &= (rng() | rng());
        if ((i & 3) == 3) x = (x & ~centre_mask) | fixed_centre;  // shared centre
        if ((i & 0xFFFF) == 0xFFFF) fixed_centre = rng() & centre_mask;
        const uint64_t canon = std::min(x, rc_mmer(x, l));
        const uint64_t hv = h.hash(x);
        auto it = seen.find(hv);
        if (it == seen.end()) {
            seen.emplace(hv, canon);
        } else if (it->second != canon) {
            if (failures++ < 10) std::cerr << "hash collision between distinct m-mers" << std::endl;
        }
    }
}

static std::string make_table(const char* path, int gk, bool symmetric, uint64_t seed) {
    const uint64_t n = uint64_t(1) << gk;
    std::mt19937_64 rng(seed);
    std::vector<uint64_t> raw(n, n);
    auto mirror = [&](uint64_t key) -> uint64_t {
        uint64_t r = 0;
        for (int i = 0; i != gk; ++i) r |= ((key >> i) & 1) << (gk - 1 - i);
        return r;
    };
    uint64_t next_rank = 0;
    for (uint64_t key = 0; key != n; ++key) {
        const uint64_t mir = mirror(key);
        if (symmetric and mir < key) {
            raw[key] = raw[mir];
        } else if (rng() % 3 == 0) {
            raw[key] = next_rank++;
            if (symmetric) raw[mir] = raw[key];
        }
    }
    write_gm_file_raw(path, raw);
    return path;
}

template <typename F>
static bool throws(F f, std::string const& what_contains) {
    try {
        f();
    } catch (std::exception const& e) {
        std::cout << "  rejected as expected: " << e.what() << std::endl;
        return std::string(e.what()).find(what_contains) != std::string::npos;
    }
    return false;
}

static const struct { int cfg; uint64_t x, h; } vectors[] = {
#include "gm_centre_vectors.inc"
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " DIR_WITH_w11_k{13,14,15}_canonical.gm [samples]"
                  << std::endl;
        return 1;
    }
    const std::string dir = argv[1];
    const uint64_t num_samples = argc > 2 ? std::stoull(argv[2]) : 30000000ull;
    const uint64_t l = 21;
    uint64_t grand_total = 0;

    struct Cfg { const char* table; gm_projection proj; };
    const Cfg cfgs[4] = {{"w11_k13_canonical.gm", gm_projection::AT0_CG1},
                         {"w11_k13_canonical.gm", gm_projection::AT1_CG0},
                         {"w11_k15_canonical.gm", gm_projection::AT0_CG1},
                         {"w11_k15_canonical.gm", gm_projection::AT1_CG0}};
    hasher_type hs[4];
    for (int c = 0; c != 4; ++c) {
        const std::string path = dir + "/" + cfgs[c].table;
        hs[c].seed(1234);
        hs[c].configure(custom_order_mode::GREEDYMINI, l, path, cfgs[c].proj);
        if (!hs[c].strand_invariant()) {
            std::cerr << "expected a strand-invariant order" << std::endl;
            ++failures;
        }
        string_ref sr;
        {
            auto raw = load_gm_file_raw(path);
            GreedyMiniTable t(raw);
            sr.ranks.assign(t.data(), t.data() + t.size());
            sr.gk = t.gm_k();
        }
        sr.l = (int)l;
        sr.flip = cfgs[c].proj == gm_projection::AT1_CG0;
        sr.seed = 1234;
        const uint64_t before = failures;
        const uint64_t total = run_checks(hs[c], &sr, l, num_samples, 1 + c);
        grand_total += total;
        run_injectivity(hs[c], l);
        std::cout << cfgs[c].table << " " << gm_projection_name(cfgs[c].proj) << " l=" << l << ": "
                  << total << " m-mers: fast == loop reference == string reference, "
                  << "rc-invariant; injective (4M canonical + shared-centre samples)"
                  << (failures == before ? "" : "  FAILED") << std::endl;
    }
    std::cout << "total m-mers compared against both references: " << grand_total << std::endl;

    {  // Python reference vectors
        uint64_t n = 0;
        for (auto const& v : vectors) {
            if (hs[v.cfg].hash(v.x) != v.h) report("python vector", v.x, hs[v.cfg].hash(v.x), v.h);
            ++n;
        }
        std::cout << n << " Python-reference vectors matched" << std::endl;
    }

    {  // rejections
        std::cout << "rejections:" << std::endl;
        hasher_type h;
        const char* tmp = "/tmp/test_gm_centre_synth.gm";
        make_table(tmp, 14, true, 3);  // symmetric, but m - gm_k = 7 is odd
        if (!throws([&] { h.configure(custom_order_mode::GREEDYMINI, 21, tmp,
                                      gm_projection::AT0_CG1); },
                    "even")) {
            std::cerr << "gm_k = 14 at m=21 must be rejected" << std::endl;
            ++failures;
        }
        make_table(tmp, 9, false, 11);  // asymmetric
        if (!throws([&] { h.configure(custom_order_mode::GREEDYMINI, 13, tmp,
                                      gm_projection::AT0_CG1); },
                    "mirror-symmetric")) {
            std::cerr << "asymmetric table must be rejected" << std::endl;
            ++failures;
        }
        if (!throws([&] { h.configure(custom_order_mode::NATIVE, 21, "",
                                      gm_projection::AT1_CG0); },
                    "greedymini")) {
            std::cerr << "a projection with a non-greedymini order must be rejected" << std::endl;
            ++failures;
        }
        std::remove(tmp);
    }

    {  // synthetic symmetric tables: exhaustive sweep (l small), incl. gm_k == l and i = 1..3
        struct Case { int gk; uint64_t l; };
        const char* tmp = "/tmp/test_gm_centre_synth.gm";
        for (Case c : {Case{9, 11}, Case{8, 12}, Case{7, 7}, Case{6, 12}, Case{10, 12}}) {
            make_table(tmp, c.gk, true, 5 + c.gk);
            for (auto proj : {gm_projection::AT0_CG1, gm_projection::AT1_CG0}) {
                hasher_type h(99);
                h.configure(custom_order_mode::GREEDYMINI, c.l, tmp, proj);
                string_ref sr;
                {
                    auto raw = load_gm_file_raw(tmp);
                    GreedyMiniTable t(raw);
                    sr.ranks.assign(t.data(), t.data() + t.size());
                }
                sr.l = (int)c.l; sr.gk = c.gk; sr.flip = proj == gm_projection::AT1_CG0; sr.seed = 99;
                const uint64_t before = failures;
                const uint64_t total = run_checks(h, &sr, c.l, 2000000, 77);
                run_injectivity(h, c.l);
                std::cout << "synthetic symmetric gm_k=" << c.gk << " l=" << c.l << " "
                          << gm_projection_name(proj) << ": " << total << " m-mers ok"
                          << (failures == before ? "" : "  FAILED") << std::endl;
            }
        }
        std::remove(tmp);
    }

    if (failures) {
        std::cerr << failures << " FAILURES" << std::endl;
        return 1;
    }
    std::cout << "OK" << std::endl;
    return 0;
}
