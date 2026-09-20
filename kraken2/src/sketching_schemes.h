/*
 * Sketching + ordering for Kraken 2 minimizers.
 *
 * SketchingScheme — how the window representative is chosen:
 *   MINIMIZER — classic leftmost min over w ℓ-mers
 *   MOD       — mod-sampling (Groot Koerkamp & Pibiri 2024)
 *
 * Order — how substrings are ranked (O_ℓ for MINIMIZER, O_t for MOD):
 *   RANDOM_XOR — XOR toggle-mask
 *   GM         — GreedyMini UHS table
 *   RANDOM_PCG — PCG RXS-M-XS scramble of (canonical ^ toggle)
 */

#ifndef KRAKEN2_SKETCHING_SCHEMES_H_
#define KRAKEN2_SKETCHING_SCHEMES_H_

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace kraken2 {

enum class SketchingScheme : uint8_t {
  MINIMIZER = 0,
  MOD       = 1
};

enum class Order : uint8_t {
  RANDOM_XOR = 0,
  GM         = 1,
  RANDOM_PCG = 2
};

inline const char *sketching_scheme_name(SketchingScheme s) {
  return s == SketchingScheme::MOD ? "Mod" : "Minimizer";
}

inline const char *order_name(Order o) {
  switch (o) {
    case Order::GM:         return "GreedyMini";
    case Order::RANDOM_PCG: return "RandomPcg";
    case Order::RANDOM_XOR:
    default:                return "RandomXor";
  }
}

inline std::string to_lower_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

inline SketchingScheme parse_sketching_scheme(const std::string &s) {
  const std::string lower = to_lower_copy(s);
  if (lower == "minimizer")
    return SketchingScheme::MINIMIZER;
  if (lower == "mod" || lower == "modminimizer")
    return SketchingScheme::MOD;
  throw std::invalid_argument(
      "unknown sketching scheme '" + s + "' (choices: Minimizer, Mod)");
}

inline Order parse_order(const std::string &s) {
  const std::string lower = to_lower_copy(s);
  if (lower == "random" || lower == "randomxor" || lower == "random_xor" ||
      lower == "xor")
    return Order::RANDOM_XOR;
  if (lower == "greedymini" || lower == "gm")
    return Order::GM;
  if (lower == "randompcg" || lower == "random_pcg" || lower == "pcg")
    return Order::RANDOM_PCG;
  throw std::invalid_argument(
      "unknown order '" + s +
      "' (choices: RandomXor, GreedyMini, RandomPcg)");
}

}  // namespace kraken2

#endif
