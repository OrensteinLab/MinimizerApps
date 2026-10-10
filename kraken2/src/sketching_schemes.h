/*
 * Sketching + ordering for Kraken 2 minimizers.
 *
 * SketchingScheme — how the window representative is chosen
 * (CLI --sketching):
 *   MINIMIZER — "minimizer": classic leftmost min over w l-mers
 *   MOD       — "mod-minimizer": strand-symmetric mod-sampling (mod-sampling,
 *               Groot Koerkamp & Pibiri 2024, with a per-window preferred
 *               strand; see mod_minimizer.h)
 *
 * Order — how substrings are ranked (O_l for MINIMIZER, O_t for MOD)
 * (CLI --order):
 *   RANDOM_XOR — "xor": XOR toggle-mask
 *   GM         — "greedymini": GreedyMini UHS table
 *   RANDOM_PCG — "xor-pcg-hash": PCG RXS-M-XS scramble of (canonical ^ toggle)
 *
 * GmProjection — the DNA -> binary projection of the GreedyMini order
 * (CLI --gm-projection): AT0_CG1 "at0-cg1" (A,T -> 0, C,G -> 1; default) or
 * AT1_CG0 "at1-cg0" (flipped), as in the sshash and minimap2 forks.
 *
 * The numeric ids are persisted in opts.k2d; do not renumber. Sketching id 1
 * was the earlier forward-only (not strand-symmetric) mod-minimizer: databases
 * built with it are rejected (rebuild them).
 */

#ifndef KRAKEN2_SKETCHING_SCHEMES_H_
#define KRAKEN2_SKETCHING_SCHEMES_H_

#include <cstdint>
#include <stdexcept>
#include <string>

namespace kraken2 {

enum class SketchingScheme : uint8_t {
  MINIMIZER = 0,
  MOD       = 2  // 1: retired forward-only mod-minimizer
};

enum class Order : uint8_t {
  RANDOM_XOR = 0,
  GM         = 1,
  RANDOM_PCG = 2
};

enum class GmProjection : uint8_t {
  AT0_CG1 = 0,
  AT1_CG0 = 1
};

inline const char *gm_projection_name(GmProjection p) {
  return p == GmProjection::AT1_CG0 ? "at1-cg0" : "at0-cg1";
}

inline GmProjection parse_gm_projection(const std::string &s) {
  if (s == "at0-cg1")
    return GmProjection::AT0_CG1;
  if (s == "at1-cg0")
    return GmProjection::AT1_CG0;
  throw std::invalid_argument("unknown gm-projection '" + s +
                              "' (choices: at0-cg1, at1-cg0)");
}

inline const char *sketching_scheme_name(SketchingScheme s) {
  return s == SketchingScheme::MOD ? "mod-minimizer" : "minimizer";
}

inline const char *order_name(Order o) {
  switch (o) {
    case Order::GM:         return "greedymini";
    case Order::RANDOM_PCG: return "xor-pcg-hash";
    case Order::RANDOM_XOR:
    default:                return "xor";
  }
}

// Scheme label (sketching + order): xor, xor-pcg-hash, greedymini,
// mod-minimizer-xor, mod-minimizer-xor-pcg-hash, mod-minimizer-greedymini.
// (The benchmark does not use mod-minimizer-xor-pcg-hash.)
inline std::string scheme_label(SketchingScheme s, Order o) {
  if (s == SketchingScheme::MOD)
    return std::string("mod-minimizer-") + order_name(o);
  return order_name(o);
}

// Exact lowercase names only; no aliases.
inline SketchingScheme parse_sketching_scheme(const std::string &s) {
  if (s == "minimizer")
    return SketchingScheme::MINIMIZER;
  if (s == "mod-minimizer")
    return SketchingScheme::MOD;
  throw std::invalid_argument(
      "unknown sketching '" + s + "' (choices: minimizer, mod-minimizer)");
}

inline Order parse_order(const std::string &s) {
  if (s == "xor")
    return Order::RANDOM_XOR;
  if (s == "greedymini")
    return Order::GM;
  if (s == "xor-pcg-hash")
    return Order::RANDOM_PCG;
  throw std::invalid_argument(
      "unknown order '" + s + "' (choices: xor, xor-pcg-hash, greedymini)");
}

}  // namespace kraken2

#endif
