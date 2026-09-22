// Custom minimizer order integration: swaps sshash's own hasher_type
// (upstream: mixer_64, still reproduced unchanged below for the "native"
// order) for one of the three orders implemented once in
// common/minimizer-orders/ and already shared by kraken2 and kmc -- see
// orders/README.md for the order reference and CLI-flag mapping.
//
// sshash pre-canonicalizes the m-mer handed to hasher_type::hash(): see
// util::canonical_mmer_at, which returns min(forward, reverse-complement) as
// a raw packed integer. That is a *different* canonicalization rule than
// GreedyMini's own "orient by whichever strand has the better UHS rank", so
// the GreedyMini path below still needs a full dual-orientation lookup
// (ported from kraken2::GMRanker in kraken2/src/minimizer_rankers.h) rather
// than a single-orientation one: GMRanker::Rank is provably invariant to
// which of {x, rc(x)} it is given (it always recomputes both internally and
// keeps the better), so feeding it either of sshash's canonical picks is
// safe.
//
// Requires building with -D SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING (see
// CMakeLists.txt, default flipped on for this fork), so sshash's 2-bit
// packing (A=0,C=1,G=2,T=3) matches what greedymini_rank.h/dna_revcomp.h
// assume.

#pragma once

#include "hash_util.hpp"

#include "../../common/minimizer-orders/greedymini_rank.h"
#include "../../common/minimizer-orders/pcg_rxs_m_xs.h"
#include "../../common/minimizer-orders/dna_revcomp.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sshash {

enum class custom_order_mode : uint8_t {
    NATIVE = 0,        // sshash's own order (mixer_64), unchanged
    XOR_HASH = 1,      // RandomXor
    XOR_PCG_HASH = 2,  // RandomPcg
    GREEDYMINI = 3     // GreedyMini
};

inline const char* custom_order_mode_name(custom_order_mode mode) {
    switch (mode) {
        case custom_order_mode::XOR_HASH:
            return "xor-hash";
        case custom_order_mode::XOR_PCG_HASH:
            return "xor-pcg-hash";
        case custom_order_mode::GREEDYMINI:
            return "greedymini";
        case custom_order_mode::NATIVE:
        default:
            return "native";
    }
}

inline custom_order_mode parse_custom_order_mode(std::string const& s) {
    if (s == "native") return custom_order_mode::NATIVE;
    if (s == "xor-hash") return custom_order_mode::XOR_HASH;
    if (s == "xor-pcg-hash") return custom_order_mode::XOR_PCG_HASH;
    if (s == "greedymini") return custom_order_mode::GREEDYMINI;
    throw std::invalid_argument("unknown order '" + s +
                                "' (choices: native, xor-hash, xor-pcg-hash, greedymini)");
}

/*
    Full order name -> (base order, is_mod_minimizer). mod-minimizer
    (Groot Koerkamp & Pibiri, WABI 2024) is an orthogonal flag, not a 5th-8th
    custom_order_mode: it changes which locus is sampled (see
    util::compute_mod_minimizer_single_strand), not what O_t/O_m computes for
    a candidate, so the 4 base orders above are reused unchanged as O_t. The 4
    "mod-minimizer-*" names spell their base order slightly differently from
    the bare base-order name (an explicit "-hash" suffix for native, and
    "pcg-xor-hash" rather than "xor-pcg-hash" for the PCG order) purely to
    keep all four mod-minimizer names visually parallel; both spellings are
    accepted here.
*/
inline std::pair<custom_order_mode, bool> parse_order_name(std::string const& s) {
    static const std::string prefix = "mod-minimizer-";
    if (s.rfind(prefix, 0) == 0) {
        std::string base = s.substr(prefix.size());
        if (base == "native-hash") base = "native";
        else if (base == "pcg-xor-hash") base = "xor-pcg-hash";
        return {parse_custom_order_mode(base), true};
    }
    return {parse_custom_order_mode(s), false};
}

/* Inverse of parse_order_name, for printing/logging. */
inline std::string order_name(custom_order_mode mode, bool mod_minimizer) {
    if (!mod_minimizer) return custom_order_mode_name(mode);
    switch (mode) {
        case custom_order_mode::NATIVE:
            return "mod-minimizer-native-hash";
        case custom_order_mode::XOR_HASH:
            return "mod-minimizer-xor-hash";
        case custom_order_mode::XOR_PCG_HASH:
            return "mod-minimizer-pcg-xor-hash";
        case custom_order_mode::GREEDYMINI:
        default:
            return "mod-minimizer-greedymini";
    }
}

/*
    Replaces sshash's hasher_type (upstream: mixer_64). One binary, order
    chosen at runtime (build_configuration::order), same "no per-order build
    variant" posture as kraken2/kmc. NATIVE reproduces mixer_64's exact
    formula locally (rather than holding a sshash::mixer_64 member) purely to
    avoid this header depending on hash_util.hpp's ordering; behavior is
    byte-identical to upstream when no --order flag is passed.
*/
struct custom_order_hasher {
    custom_order_hasher() { seed(0); }
    explicit custom_order_hasher(const uint64_t s) { seed(s); }

    void seed(const uint64_t s) {
        m_seed = s;
        m_native_magic = pthash::xxhash_64::hash(s, 0).first();
    }

    /*
        Called once at build time (dictionary_builder::build), right after
        seed(). `l` is sshash's minimizer length m -- the width, in bases, of
        the value passed to hash(). No-op for NATIVE/XOR_HASH/XOR_PCG_HASH.
    */
    void configure(custom_order_mode mode, uint64_t l,
                   std::string const& greedy_mini_filename = "",
                   GmAlphabet gm_alphabet = GmAlphabet::BINARY) {
        m_mode = mode;
        m_l = (int)l;
        if (mode != custom_order_mode::GREEDYMINI) return;

        auto raw = load_gm_file_raw(greedy_mini_filename);
        if (raw.empty()) {
            throw std::runtime_error("failed to load GreedyMini table: " + greedy_mini_filename);
        }
        GreedyMiniTable table(raw, gm_alphabet);
        load_table(table);
    }

    inline uint64_t hash(uint64_t value) const {
        switch (m_mode) {
            case custom_order_mode::XOR_HASH:
                return value ^ m_seed;
            case custom_order_mode::XOR_PCG_HASH:
                return pcg_rxs_m_xs64(value ^ m_seed);
            case custom_order_mode::GREEDYMINI:
                return m_alphabet == GmAlphabet::DNA ? rank_dna_alphabet(value)
                                                     : rank_binary_dna(value);
            case custom_order_mode::NATIVE:
            default:
                return (value * 0x517cc1b727220a95ull) ^ m_native_magic;
        }
    }

    custom_order_mode mode() const { return m_mode; }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_mode);
        visitor.visit(t.m_seed);
        visitor.visit(t.m_native_magic);
        visitor.visit(t.m_l);
        visitor.visit(t.m_gm_k);
        visitor.visit(t.m_alphabet);
        visitor.visit(t.m_ext_bits);
        visitor.visit(t.m_gamma_bits);
        visitor.visit(t.m_key_bits);
        visitor.visit(t.m_gm_ranks);
    }

    void load_table(GreedyMiniTable const& table) {
        if (table.empty()) {
            throw std::invalid_argument("GreedyMini order requires a non-empty table");
        }
        if (table.gm_k() > m_l) {
            throw std::invalid_argument("GreedyMini gm_k exceeds minimizer length m");
        }
        m_alphabet = table.alphabet();
        m_gm_k = table.gm_k();
        m_key_bits = table.key_bits();
        if (m_alphabet == GmAlphabet::DNA) {
            m_ext_bits = 0;
            m_gamma_bits = 2 * (m_l - m_gm_k);
        } else {
            m_ext_bits = m_l - m_gm_k;
            m_gamma_bits = m_l;  // Gamma: 1 bit/base (DNA, bpc-1 = 1)
        }
        if (16 + m_key_bits + m_ext_bits + m_gamma_bits > 63) {
            throw std::invalid_argument(
                "GreedyMini tiebreaker width leaves no room to pack the primary rank");
        }
        m_gm_ranks.assign(table.data(), table.data() + table.size());
    }

    uint64_t tb_toggle() const {
        const uint64_t mask = (m_gamma_bits >= 64) ? ~uint64_t(0) : ((uint64_t(1) << m_gamma_bits) - 1);
        return m_seed & mask;
    }

    /*
        Dual-orientation lookup for the binary-alphabet backend (sigma-extend
        + Extend-k), ported from kraken2::GMRanker::rank_binary_dna, with one
        deviation kraken2 doesn't need: a canonical GreedyMini table
        deliberately assigns a *mirror pair* of Sigma keys the same rank (and
        every non-UHS key shares the same "infinity" sentinel), so `rank`
        alone is not injective in the key. kraken2's own scanner tolerates
        that, but sshash's minimizer_iterator assumes a hash tie implies a
        value tie (its incremental tie-break bookkeeping asserts exactly
        that). So the packed value here always folds the raw Sigma `key` in
        as a tie-break *below* rank: rank stays the primary, ascending sort
        key (kInfRank sorts last, so "any UHS beats any non-UHS" is
        preserved), and `key` deterministically disambiguates every rank tie,
        making hash() properly injective in the underlying m-mer.
    */
    uint64_t rank_binary_dna(uint64_t lmer) const {
        const uint16_t* ranks = m_gm_ranks.data();
        const uint64_t toggle = tb_toggle();

        auto probe = [&](uint64_t x) {
            const uint64_t atcg = gm_acgt_to_atcg(x);
            const uint32_t key = gm_extract_dna_key_atcg(atcg, m_l, m_gm_k);
            return std::make_tuple(ranks[key], key, atcg);
        };
        auto finish = [&](std::tuple<uint16_t, uint32_t, uint64_t> const& s) -> uint64_t {
            const uint16_t rank = std::get<0>(s);
            const uint32_t key = std::get<1>(s);
            const uint64_t atcg = std::get<2>(s);
            uint64_t mid = (((uint64_t)rank << m_key_bits) | key) << m_ext_bits;
            if (m_ext_bits > 0) mid |= gm_extract_suffix_sigma_atcg(atcg, m_l, m_gm_k);
            const uint64_t gamma = gm_extract_dna_tb_atcg(atcg, m_l) ^ toggle;
            return (mid << m_gamma_bits) | gamma;
        };

        const auto a = probe(lmer);
        const uint64_t rc = dna_reverse_complement(lmer, (uint8_t)m_l);
        if (rc == lmer) return finish(a);
        const auto b = probe(rc);
        return pick_better_primary(std::get<0>(a), std::get<0>(b), [&] { return finish(a); },
                                   [&] { return finish(b); });
    }

    /*
        Dual-orientation lookup for the DNA-alphabet backend (native 4^gm_k
        table, no bit remap), ported from kraken2::GMRanker::rank_dna_alphabet;
        see rank_binary_dna's comment for why `key` is always folded into the
        packed value below `rank`.
    */
    uint64_t rank_dna_alphabet(uint64_t lmer) const {
        const uint16_t* ranks = m_gm_ranks.data();
        const uint64_t toggle = tb_toggle();

        auto probe = [&](uint64_t x) {
            const uint32_t key = gm_extract_dna_alphabet_key(x, m_l, m_gm_k);
            return std::make_tuple(ranks[key], key, x);
        };
        auto finish = [&](std::tuple<uint16_t, uint32_t, uint64_t> const& s) -> uint64_t {
            const uint16_t rank = std::get<0>(s);
            const uint32_t key = std::get<1>(s);
            const uint64_t x = std::get<2>(s);
            const uint64_t primary = ((uint64_t)rank << m_key_bits) | key;
            const uint64_t tb = gm_extract_dna_alphabet_tb(x, m_l, m_gm_k) ^ toggle;
            return (primary << m_gamma_bits) | tb;
        };

        const auto a = probe(lmer);
        const uint64_t rc = dna_reverse_complement(lmer, (uint8_t)m_l);
        if (rc == lmer) return finish(a);
        const auto b = probe(rc);
        return pick_better_primary(std::get<0>(a), std::get<0>(b), [&] { return finish(a); },
                                   [&] { return finish(b); });
    }

    /*
        Smaller raw rank wins (kInfRank sorts last, so "any UHS beats any
        non-UHS" falls out of plain numeric comparison); equal raw rank ->
        smaller packed (injective) value breaks the tie deterministically.
    */
    template <typename FinishA, typename FinishB>
    static uint64_t pick_better_primary(uint16_t pa, uint16_t pb, FinishA finish_a,
                                        FinishB finish_b) {
        if (pa < pb) return finish_a();
        if (pb < pa) return finish_b();
        const uint64_t ra = finish_a();
        const uint64_t rb = finish_b();
        return ra < rb ? ra : rb;
    }

    custom_order_mode m_mode = custom_order_mode::NATIVE;
    uint64_t m_seed = 0;
    uint64_t m_native_magic = 0;
    int m_l = 0;
    int m_gm_k = 0;
    GmAlphabet m_alphabet = GmAlphabet::BINARY;
    int m_ext_bits = 0;
    int m_gamma_bits = 0;
    int m_key_bits = 0;
    std::vector<uint16_t> m_gm_ranks;  // embedded in the .sshash file; no sidecar needed
};

using hasher_type = custom_order_hasher;

}  // namespace sshash
