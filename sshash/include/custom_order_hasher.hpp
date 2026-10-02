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
    backend. The Sigma plane (MSB of the ATCG-remapped base, A=00,T=01,C=10,
    G=11) classifies {A,T} vs {C,G}; AT0_CG1 (default, the historical
    behaviour) maps AT->0, CG->1, AT1_CG0 flips every Sigma bit (AT->1, CG->0).
    Gamma (the LSB plane, the tie-break) is untouched by the flip.
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
    The hash of the strand-dependent orders (native, xor-hash, xor-pcg-hash) -- three words,
    no table -- as a small value type. minimizer_iterator keeps one by value, next to its own
    state, instead of a pointer to (or copy of) the full hasher, which for GreedyMini carries a
    rank table and ~100 bytes of derived constants.
*/
struct plain_hasher {
    custom_order_mode mode = custom_order_mode::NATIVE;
    uint64_t seed = 0;
    uint64_t native_magic = 0;

    inline uint64_t hash(uint64_t value) const {
        switch (mode) {
            case custom_order_mode::XOR_HASH:
                return value ^ seed;
            case custom_order_mode::XOR_PCG_HASH:
                return pcg_rxs_m_xs64(value ^ seed);
            case custom_order_mode::NATIVE:
            default:
                return (value * 0x517cc1b727220a95ull) ^ native_magic;
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
        m_sigma_flip = 0;
        if (mode != custom_order_mode::GREEDYMINI) {
            if (projection != gm_projection::AT0_CG1) {
                throw std::invalid_argument("--gm-projection only applies to -O greedymini");
            }
            return;
        }
        if (projection == gm_projection::AT1_CG0) {
            /* Sigma = odd bits of the ATCG word; flip them for the l bases. */
            const uint64_t odd = 0xAAAAAAAAAAAAAAAAull;
            m_sigma_flip = (l >= 32) ? odd : (odd & ((uint64_t(1) << (2 * l)) - 1));
        }

        auto raw = load_gm_file_raw(greedy_mini_filename);
        if (raw.empty()) {
            throw std::runtime_error("failed to load GreedyMini table: " + greedy_mini_filename);
        }
        GreedyMiniTable table(raw, GmAlphabet::BINARY);
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

    inline uint64_t hash(uint64_t value) const {
        return m_mode == custom_order_mode::GREEDYMINI ? hash_gm(value) : hash_plain(value);
    }

    /* hash() split by strand-invariance, so that hot loops that know the order's kind
       (minimizer_iterator, util::compute_minimizer) inline only the half they need:
       hash_plain() for native / xor-hash / xor-pcg-hash (which expect the numerically
       canonical m-mer), hash_gm() for GreedyMini (which accepts either strand). */
    inline uint64_t hash_plain(uint64_t value) const { return plain().hash(value); }

    plain_hasher plain() const { return plain_hasher{m_mode, m_seed, m_native_magic}; }

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
        p.S = gm_pext_even(x ^ (x >> 1), m_even_l) ^ m_flip_l;
        p.G = gm_pext_odd(x, m_odd_l);
        const uint64_t T = reverse_bits64(p.S | (p.G << 32));
        p.Sb = T >> (64 - m_l);
        p.Gb = ((T & 0xFFFFFFFFull) >> (32 - m_l)) ^ m_l_mask;
        return p;
    }

    /* Sigma plane only, which is all the rank of a window needs (see centre_rank): of a whole
       k-mer (k <= 63; low/high 64-bit halves of the packed word, bits above k-1 zero), and
       slid by one base: the base b (A=0,C=1,G=2,T=3) enters at the highest position k-1 of
       the word, the lowest position leaves. */
    inline uint64_t sigma_of_kmer(uint64_t lo, uint64_t hi, int k) const {
        const uint64_t evens = 0x5555555555555555ull;
        const uint64_t kmask = (uint64_t(1) << k) - 1;
        uint64_t S = (gm_pext_even(lo ^ (lo >> 1), evens) | (gm_pext_even(hi ^ (hi >> 1), evens) << 32)) &
                     kmask;
        return m_sigma_flip ? S ^ kmask : S;
    }

    inline uint64_t roll_sigma(uint64_t S, uint64_t b, int k) const {
        return (S >> 1) | ((((b ^ (b >> 1)) & 1) ^ m_flip_bit) << (k - 1));
    }

    /* The primary (rank) field of the hash of the m-mer at positions i .. i+m-1 of a k-mer
       with Sigma plane S: one table load. hash >> 2m == centre_rank, so comparing ranks first
       and the full hashes only on equal ranks orders the windows exactly as the hashes do.
       GreedyMini only. */
    inline uint16_t centre_rank(uint64_t S, int i) const {
        return m_gm_ranks[(S >> (i + m_centre_shift)) & m_key_mask];
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
        const uint64_t pa = (S << m_l) | (G ^ m_toggle);
        const uint64_t pb = (Sb << m_l) | (Gb ^ m_toggle);
        const uint64_t tb = pa < pb ? pa : pb;
        return ((uint64_t)m_gm_ranks[(S >> m_centre_shift) & m_key_mask] << m_rank_shift) | tb;
    }

    /*
        Slow, per-base specification of the hash that hash_gm() must equal bit for bit (used
        by the tests, independent of the bit tricks of hash_from_planes): builds the Sigma/Gamma strings of the m-mer and of its reverse
        complement base by base. hash = (rank[Sigma(y)] << 2l) | min(F(x), F(rc(x))) with
        F(v) = (Sigma(v) << l) | (Gamma(v) ^ toggle), y = the middle gm_k bases of x.
    */
    uint64_t hash_reference(uint64_t value) const {
        const int l = m_l, gk = m_gm_k, i = (l - gk) / 2;
        const uint64_t toggle = m_seed & ((uint64_t(1) << l) - 1);
        auto planes = [&](uint64_t v, uint64_t& sigma, uint64_t& gamma) {
            sigma = gamma = 0;
            for (int p = 0; p != l; ++p) {  // base p = 0 is the one in the highest bits
                const uint64_t b = (v >> (2 * (l - 1 - p))) & 3;  // A=0,C=1,G=2,T=3
                const uint64_t sig = (b == 1 or b == 2) ? 1 : 0;  // {C,G} vs {A,T}
                const uint64_t gam = (b >> 1) & 1;                // {A,C} vs {G,T}
                sigma = (sigma << 1) | (sig ^ (m_projection == gm_projection::AT1_CG0 ? 1 : 0));
                gamma = (gamma << 1) | gam;
            }
        };
        uint64_t s1, g1, s2, g2;
        planes(value, s1, g1);
        planes(dna_reverse_complement(value, (uint8_t)l), s2, g2);
        uint64_t key = 0;
        for (int p = i; p != i + gk; ++p) key = (key << 1) | ((s1 >> (l - 1 - p)) & 1);
        const uint64_t f1 = (s1 << l) | (g1 ^ toggle);
        const uint64_t f2 = (s2 << l) | (g2 ^ toggle);
        return ((uint64_t)m_gm_ranks[key] << (2 * l)) | (f1 < f2 ? f1 : f2);
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
        visitor.visit(t.m_sigma_flip);
    }

    /* Preconditions: (l - gm_k) even, 16 + 2l <= 63, and a mirror-symmetric
       table (rank[key] == rank[reverse_bits(key)] for every key) -- the single table probe
       is strand invariant only then. */
    void validate() const {
        if (((m_l - m_gm_k) & 1) != 0) {
            throw std::invalid_argument(
                "GreedyMini requires (m - gm_k) to be even (m = " + std::to_string(m_l) +
                ", gm_k = " + std::to_string(m_gm_k) + "): the table is extended by the same "
                "number of bases on both sides of the gm_k-mer");
        }
        if (16 + 2 * m_l > 63) {
            throw std::invalid_argument("GreedyMini: minimizer length too large to pack");
        }
        const uint64_t n = uint64_t(1) << m_gm_k;
        if (m_gm_ranks.size() != n) {
            throw std::invalid_argument("GreedyMini: unexpected table size");
        }
        for (uint64_t key = 0; key != n; ++key) {
            const uint64_t mirror = reverse_bits64(key) >> (64 - m_gm_k);
            if (m_gm_ranks[key] != m_gm_ranks[mirror]) {
                throw std::invalid_argument(
                    "GreedyMini requires a mirror-symmetric (canonical) table: rank[" +
                    std::to_string(key) + "] != rank[" + std::to_string(mirror) + "]");
            }
        }
    }

    void load_table(GreedyMiniTable const& table) {
        if (table.empty()) {
            throw std::invalid_argument("GreedyMini order requires a non-empty table");
        }
        if (table.gm_k() > m_l) {
            throw std::invalid_argument("GreedyMini gm_k exceeds minimizer length m");
        }
        m_gm_k = table.gm_k();
        m_gm_ranks.assign(table.data(), table.data() + table.size());
        validate();
    }

    /* seed toggle of the Gamma plane (1 bit per base) */
    uint64_t tb_toggle() const { return m_seed & ((uint64_t(1) << m_l) - 1); }

    /* ---- derived (non-serialized) state ---- */

    static inline uint64_t reverse_bits64(uint64_t v) {
#if defined(__has_builtin)
#if __has_builtin(__builtin_bitreverse64)
        return __builtin_bitreverse64(v);
#endif
#endif
        v = __builtin_bswap64(v);
        v = ((v >> 4) & 0x0F0F0F0F0F0F0F0Full) | ((v & 0x0F0F0F0F0F0F0F0Full) << 4);
        v = ((v >> 2) & 0x3333333333333333ull) | ((v & 0x3333333333333333ull) << 2);
        v = ((v >> 1) & 0x5555555555555555ull) | ((v & 0x5555555555555555ull) << 1);
        return v;
    }

    /*
        Recompute everything hash() hoists out of the per-call path (toggle,
        masks, shifts) from the serialized members only; called after
        configure() and after a load (see visit()).
    */
    void derive() {
        const int l = m_l, gk = m_gm_k;
        m_toggle = tb_toggle();
        m_l_mask = (uint64_t(1) << l) - 1;
        m_flip_l = m_sigma_flip ? m_l_mask : 0;
        m_flip_bit = m_sigma_flip ? 1 : 0;
        const uint64_t evens = 0x5555555555555555ull;
        m_even_l = evens & ((uint64_t(1) << (2 * l)) - 1);
        m_odd_l = m_even_l << 1;
        m_key_mask = (uint64_t(1) << gk) - 1;
        m_centre_shift = (l - gk) / 2;
        m_rank_shift = 2 * l;
    }

    custom_order_mode m_mode = custom_order_mode::NATIVE;
    uint64_t m_seed = 0;
    uint64_t m_native_magic = 0;
    int m_l = 0;
    int m_gm_k = 0;
    gm_projection m_projection = gm_projection::AT0_CG1;
    uint64_t m_sigma_flip = 0;  // XOR mask on the ATCG word: flips every Sigma bit, or 0
    std::vector<uint16_t> m_gm_ranks;  // embedded in the .sshash file; no sidecar needed

    /* derived by derive(); never serialized */
    uint64_t m_flip_bit = 0, m_toggle = 0, m_l_mask = 0, m_flip_l = 0, m_key_mask = 0;
    uint64_t m_even_l = 0, m_odd_l = 0;
    int m_rank_shift = 0, m_centre_shift = 0;
};

using hasher_type = custom_order_hasher;

}  // namespace sshash
