// Custom minimizer order integration: swaps sshash's own hasher_type
// (upstream: mixer_64, still reproduced unchanged below for the "native"
// order) for one of the orders implemented once in common/minimizer-orders/
// and shared by kraken2 and kmc -- see orders/README.md for the order
// reference and CLI-flag mapping.
//
// GreedyMini ranks an m-mer by its centre gm_k-mer y,
// m = |x| + gm_k + |z| with |x| = |z|, through the canonical table; the
// table is mirror symmetric, so hash(x) == hash(rc(x)) (strand_invariant())
// and the minimizer code hashes the RAW forward m-mer, canonicalizing only
// the value of a locus that becomes the minimum. hash() is computed from the
// Sigma/Gamma bit planes of the m-mer and of its reverse complement
// (hash_from_planes), which the strand-invariant iterator rolls along the
// k-mer; hash_reference() is the slow per-base specification used by tests.
//
// Works on sshash's own 2-bit packing, A=00, C=01, G=11, T=10 (upstream's default encoding; the
// traditional one, -D SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING, is not supported by this
// header): bit 0 of a base is its Sigma (A,T -> 0; C,G -> 1), bit 1 its Gamma (A,C -> 0; G,T -> 1),
// and complementing a base flips bit 1 only. The hash of a window is a function of the base
// sequence alone, so it is the same as with the packing greedymini_rank.h documents (ACGT = 0123).

#pragma once

#include "hash_util.hpp"

#include "../../common/minimizer-orders/greedymini_centre.h"
#include "../../common/minimizer-orders/greedymini_rank.h"
#include "../../common/minimizer-orders/pcg_rxs_m_xs.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
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
    DNA -> binary projection polarity for the GreedyMini
    backend. The Sigma plane (bit 0 of an sshash base: A=00, C=01, G=11, T=10)
    classifies {A,T} vs {C,G}; AT0_CG1 (default, the historical behaviour) maps AT->0, CG->1,
    AT1_CG0 flips every Sigma bit (AT->1, CG->0). Gamma (bit 1, the tie-break) is untouched by
    the flip.
*/
enum class gm_projection : uint8_t { AT0_CG1 = 0, AT1_CG0 = 1 };

inline const char* gm_projection_name(gm_projection p) {
    return p == gm_projection::AT1_CG0 ? "at1-cg0" : "at0-cg1";
}

inline gm_projection parse_gm_projection(std::string const& s) {
    if (s == "at0-cg1") return gm_projection::AT0_CG1;
    if (s == "at1-cg0") return gm_projection::AT1_CG0;
    throw std::invalid_argument("unknown gm-projection '" + s +
                                "' (choices: at0-cg1, at1-cg0)");
}

/*
    The order, as a compile-time constant: `dispatch_order(mode, f)` calls f with a
    std::integral_constant<custom_order_mode, mode>, so that the code of the hot paths (the
    minimizer iterators, util::compute_minimizer, streaming_query, the builder) is generated once
    per order, with the hash function inlined, instead of branching on the mode at every hash.
*/
template <custom_order_mode Mode>
using order_tag = std::integral_constant<custom_order_mode, Mode>;

template <typename F>
inline decltype(auto) dispatch_order(custom_order_mode mode, F&& f) {
    switch (mode) {
        case custom_order_mode::XOR_HASH:
            return f(order_tag<custom_order_mode::XOR_HASH>());
        case custom_order_mode::XOR_PCG_HASH:
            return f(order_tag<custom_order_mode::XOR_PCG_HASH>());
        case custom_order_mode::GREEDYMINI:
            return f(order_tag<custom_order_mode::GREEDYMINI>());
        case custom_order_mode::NATIVE:
        default:
            return f(order_tag<custom_order_mode::NATIVE>());
    }
}

/*
    The hash of a strand-dependent order (native, xor-hash, xor-pcg-hash; they expect the
    numerically canonical m-mer) as a one-word value type: native is upstream's mixer_64 with its
    magic precomputed, a multiply and a xor. The iterators keep one by value, next to their own
    state, instead of a pointer to (or copy of) the full hasher, which for GreedyMini carries a
    rank table and ~100 bytes of derived constants.
*/
template <custom_order_mode Mode>
struct plain_hasher {
    static_assert(Mode != custom_order_mode::GREEDYMINI);

    uint64_t key = 0;  // native: mixer_64's magic, xxhash_64(seed); xor-hash, xor-pcg-hash: the seed

    inline uint64_t hash(uint64_t value) const {
        if constexpr (Mode == custom_order_mode::NATIVE) {
            return (value * 0x517cc1b727220a95ull) ^ key;
        } else if constexpr (Mode == custom_order_mode::XOR_HASH) {
            return value ^ key;
        } else {
            return pcg_rxs_m_xs64(value ^ key);
        }
    }
};

/*
    Sigma/Gamma bit planes of an m-mer (m bits each, bit p describing the base at position p of
    the packed word, the bit order of `pext`), plus those of its reverse complement:
    Sb = Sigma(rc), Gb = Gamma(rc) (complementing a base flips Gamma only, Sigma = {A,T} vs
    {C,G} is strand symmetric, and reversing the word reverses the planes). Sigma carries the
    --gm-projection flip; Gamma does not carry the seed toggle (applied when the hash is packed).
*/
struct gm_planes {
    uint64_t S = 0, G = 0, Sb = 0, Gb = 0;
};

/*
    Replaces sshash's hasher_type (upstream: mixer_64). One binary, order
    chosen at runtime (build_configuration::order), same "no per-order build
    variant" posture as kraken2/kmc; the hot paths are instantiated per order (dispatch_order).
    NATIVE reproduces mixer_64's exact formula (plain_hasher<NATIVE>); behavior is
    byte-identical to upstream when no --order flag is passed.
*/
struct custom_order_hasher {
    custom_order_hasher() { seed(0); }
    explicit custom_order_hasher(const uint64_t s) { seed(s); }

    void seed(const uint64_t s) {
        m_seed = s;
        m_native_magic = pthash::xxhash_64::hash(s, 0).first();
        if (m_mode == custom_order_mode::GREEDYMINI) m_toggle = tb_toggle();  // seeded after configure()
    }

    /*
        Called once at build time (dictionary_builder::build), right after
        seed(). `l` is sshash's minimizer length m -- the width, in bases, of
        the value passed to hash(). No-op for NATIVE/XOR_HASH/XOR_PCG_HASH.
        For GreedyMini the table must be mirror symmetric, (l - gm_k) even and
        16 + 2l <= 63, else this throws.
    */
    void configure(custom_order_mode mode, uint64_t l,
                   std::string const& greedy_mini_filename = "",
                   gm_projection projection = gm_projection::AT0_CG1) {
        m_mode = mode;
        m_l = (int)l;
        m_projection = projection;
        if (mode != custom_order_mode::GREEDYMINI) {
            if (projection != gm_projection::AT0_CG1) {
                throw std::invalid_argument("--gm-projection only applies to -O greedymini");
            }
            return;
        }
#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
        throw std::invalid_argument(
            "GreedyMini needs sshash's default nucleotide encoding (A=00,C=01,G=11,T=10), not "
            "-DSSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING");
#endif

        auto raw = load_gm_file_raw(greedy_mini_filename);
        if (raw.empty()) {
            throw std::runtime_error("failed to load GreedyMini table: " + greedy_mini_filename);
        }
        GreedyMiniTable table(raw);
        load_table(table);
        derive();
    }

    /*
        True if hash(x) == hash(rc(x)) for every m-mer x, i.e. the order is
        strand-invariant, so callers may hash the raw forward m-mer instead of
        the numerically-canonical one (util::canonical_mmer_at) and only
        canonicalize the value of a locus that actually becomes the minimum.
        Holds for every GreedyMini configuration by construction; derived from
        the mode, so the on-disk format carries nothing for it.
    */
    inline bool strand_invariant() const { return m_mode == custom_order_mode::GREEDYMINI; }

    /* The hash of the order chosen at run time; the hot loops know the order at compile time
       and use plain<Mode>().hash() / hash_gm() directly (see dispatch_order). */
    inline uint64_t hash(uint64_t value) const {
        return dispatch_order(m_mode, [&](auto mode) {
            if constexpr (decltype(mode)::value == custom_order_mode::GREEDYMINI) {
                return hash_gm(value);
            } else {
                return plain<decltype(mode)::value>().hash(value);
            }
        });
    }

    /* The hash of a strand-dependent order (native, xor-hash, xor-pcg-hash), expecting the
       numerically canonical m-mer, as a value of its own. */
    template <custom_order_mode Mode>
    plain_hasher<Mode> plain() const {
        return plain_hasher<Mode>{Mode == custom_order_mode::NATIVE ? m_native_magic : m_seed};
    }

    /* The GreedyMini hash, which accepts either strand of the m-mer. */
    inline uint64_t hash_gm(uint64_t value) const { return hash_from_planes(planes_of_mmer(value)); }

    /*
        Bit-plane API. hash_gm(x) is *defined* as
        hash_from_planes(planes_of_mmer(x)); the strand-invariant minimizer iterator needs the
        full hash of a window only when ranks tie, and otherwise only its rank (centre_rank,
        from the Sigma plane of the k-mer, rolled one base per call).
    */
    /* Planes of the m-mer x (m = the hasher's l) from scratch: two pext and one 64-bit bit
       reversal (S | G << 32 reversed yields both reverse-complement planes). Bits above
       position m-1 are zero. */
    inline gm_planes planes_of_mmer(uint64_t x) const {
        gm_planes p;
        p.S = m_planes.sigma(x);
        p.G = m_planes.gamma(x);
        const uint64_t T = gm_centre::reverse_bits64(p.S | (p.G << 32));
        p.Sb = T >> (64 - m_l);
        p.Gb = ((T & 0xFFFFFFFFull) >> (32 - m_l)) ^ m_l_mask;
        return p;
    }

    /* Sigma plane only, which is all the rank of a window needs (see centre_rank): of a whole
       k-mer (k <= 63; low/high 64-bit halves of the packed word, bits above k-1 zero), and
       slid by one base: the base b (a 2-bit code of the packing) enters at the highest position
       k-1 of the word, the lowest position leaves. */
    inline uint64_t sigma_of_kmer(uint64_t lo, uint64_t hi, int k) const {
        return m_planes.sigma_wide(lo, hi, k);
    }

    inline uint64_t roll_sigma(uint64_t S, uint64_t b, int k) const {
        return gm_centre::roll_in_high(S, m_planes.base_sigma(b), k);
    }

    /* The primary (rank) field of the hash of the m-mer at positions i .. i+m-1 of a k-mer
       with Sigma plane S: one table load. hash >> 2m == centre_rank, so comparing ranks first
       and the full hashes only on equal ranks orders the windows exactly as the hashes do.
       GreedyMini only. */
    inline uint16_t centre_rank(uint64_t S, int i) const {
        return m_order.rank_at(S, i);
    }

    inline uint64_t hash_from_planes(gm_planes const& p) const {
        return hash_from_planes(p.S, p.G, p.Sb, p.Gb);
    }

    /*
        The GreedyMini hash of an m-mer from its planes (S,G) and those of its reverse
        complement (Sb,Gb), exactly m bits each. GreedyMini only.

        ONE table probe: the centre gm_k-mer of rc(x) is the mirror of the centre gm_k-mer of
        x, so a mirror-symmetric table gives the same rank (checked at configure time). The
        tie-break min over the two frames is branchless; it is the low 2l bits of the value,
        below the rank, so the packed value is a total order consistent with
        (rank, min Sigma, Gamma) and injective on canonical m-mers (the winning frame's
        (Sigma, Gamma) determine the DNA).
    */
    inline uint64_t hash_from_planes(uint64_t S, uint64_t G, uint64_t Sb, uint64_t Gb) const {
        return m_order.value_of_both_strands(S, G, Sb, Gb);
    }

    /*
        Slow, per-base specification of the hash that hash_gm() must equal bit for bit (used
        by the tests, independent of the bit tricks of hash_from_planes): builds the Sigma/Gamma strings of the m-mer and of its reverse
        complement base by base. hash = (rank[Sigma(y)] << 2l) | min(F(x), F(rc(x))) with
        F(v) = (Sigma(v) << l) | (Gamma(v) ^ toggle), y = the middle gm_k bases of x.
    */
    uint64_t hash_reference(uint64_t value) const {
        const int l = m_l;
        /* per-position 0/1 planes, position 0 = the first base (the highest bits of the word);
           the packing is A=0,C=1,G=3,T=2 */
        auto planes = [&](uint64_t v, std::vector<int>& sigma, std::vector<int>& gamma) {
            sigma.assign(l, 0);
            gamma.assign(l, 0);
            for (int p = 0; p != l; ++p) {
                const uint64_t b = (v >> (2 * (l - 1 - p))) & 3;
                const int sig = (b == 1 or b == 3) ? 1 : 0;  // {C,G} vs {A,T}
                sigma[p] = sig ^ (m_projection == gm_projection::AT1_CG0 ? 1 : 0);
                gamma[p] = (b == 2 or b == 3) ? 1 : 0;  // {A,C} vs {G,T}
            }
        };
        /* reverse complement, base by base: A <-> T (0 <-> 2), C <-> G (1 <-> 3) */
        uint64_t rc = 0;
        for (int p = 0; p != l; ++p) {
            const uint64_t b = (value >> (2 * p)) & 3;
            rc |= uint64_t(b == 0 ? 2 : b == 2 ? 0 : b == 1 ? 3 : 1) << (2 * (l - 1 - p));
        }
        std::vector<int> s1, g1, s2, g2;
        planes(value, s1, g1);
        planes(rc, s2, g2);
        return gm_centre::Reference(*m_table, l, m_seed & ((uint64_t(1) << l) - 1),
                                    /*injective=*/true)
            .value(s1, g1, s2, g2);
    }

    custom_order_mode mode() const { return m_mode; }
    gm_projection projection() const { return m_projection; }
    int gm_k() const { return m_gm_k; }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
        /* the derived, non-serialized members follow the loaded ones */
        if (m_mode == custom_order_mode::GREEDYMINI and !m_gm_ranks.empty()) derive();
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
        visitor.visit(t.m_gm_ranks);
        visitor.visit(t.m_projection);
    }

    /* Preconditions: (l - gm_k) even, 16 + 2l <= 63 and a mirror-symmetric table (the single
       table probe is strand invariant only then); the table checks are the shared
       gm_centre::validate_table. */
    void validate(GreedyMiniTable const& table) const {
        if (16 + 2 * m_l > 63) {
            throw std::invalid_argument("GreedyMini: minimizer length too large to pack");
        }
        gm_centre::validate_table(table, m_l, /*require_rank_pairs=*/false);
    }

    void load_table(GreedyMiniTable const& table) {
        if (table.empty()) {
            throw std::invalid_argument("GreedyMini order requires a non-empty table");
        }
        if (table.gm_k() > m_l) {
            throw std::invalid_argument("GreedyMini gm_k exceeds minimizer length m");
        }
        validate(table);
        m_gm_k = table.gm_k();
        m_gm_ranks.assign(table.data(), table.data() + table.size());
    }

    /* seed toggle of the Gamma plane (1 bit per base) */
    uint64_t tb_toggle() const { return m_seed & ((uint64_t(1) << m_l) - 1); }

    /* ---- derived (non-serialized) state ---- */

    /*
        Recompute everything hash() hoists out of the per-call path (toggle,
        masks, shifts) from the serialized members only; called after
        configure() and after a load (see visit()).
    */
    void derive() {
        const int l = m_l;
        m_toggle = tb_toggle();
        m_l_mask = (uint64_t(1) << l) - 1;
        m_planes = gm_centre::TwoBitPlanes<false>(l, m_projection == gm_projection::AT1_CG0);
        /* the shared order, over a table rebuilt from the serialized ranks (infinity = 0xFFFF) */
        std::vector<uint64_t> raw(m_gm_ranks.size());
        for (size_t i = 0; i != raw.size(); ++i) {
            raw[i] = m_gm_ranks[i] == GreedyMiniTable::kInfRank ? raw.size() : m_gm_ranks[i];
        }
        m_table = std::make_shared<const GreedyMiniTable>(raw);
        m_order = gm_centre::Order(m_table, l, m_seed, /*injective=*/true);
    }

    custom_order_mode m_mode = custom_order_mode::NATIVE;
    uint64_t m_seed = 0;
    uint64_t m_native_magic = 0;
    int m_l = 0;
    int m_gm_k = 0;
    gm_projection m_projection = gm_projection::AT0_CG1;
    std::vector<uint16_t> m_gm_ranks;  // embedded in the .sshash file; no sidecar needed

    /* derived by derive(); never serialized */
    uint64_t m_toggle = 0, m_l_mask = 0;
    gm_centre::TwoBitPlanes<false> m_planes;  // sshash's encoding: Sigma = bit 0, Gamma = bit 1
    std::shared_ptr<const GreedyMiniTable> m_table;
    gm_centre::Order m_order;
};

using hasher_type = custom_order_hasher;

/* whether a hasher (a plain_hasher<Mode> or the GreedyMini hasher itself) is strand-invariant */
template <typename Hasher>
inline constexpr bool strand_invariant_hasher_v = std::is_same_v<Hasher, custom_order_hasher>;

}  // namespace sshash
