// Tests of mod_minimizer_rule.h: le_revcomp_planes (O(1)) against the
// per-base window_le_revcomp on random and low-complexity windows of every
// length 1..64, and CenterPick against a direct evaluation of the center rule.
// Build and run: make test

#include "mod_minimizer_rule.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

static int failures = 0;

static void check(bool ok, const char *what) {
  if (!ok && failures++ < 10)
    std::cerr << "FAILED: " << what << "\n";
}

int main() {
  std::mt19937_64 rng(11);
  for (int len = 1; len <= 64; ++len) {
    for (int it = 0; it < 4000; ++it) {
      std::vector<int> b((size_t)len);
      const int kind = it % 4;  // i.i.d., {A,C}, palindromes, one change
      for (int j = 0; j < len; ++j)
        b[(size_t)j] = kind == 1 ? (int)(rng() & 1) : (int)(rng() & 3);
      if (kind >= 2)
        for (int j = 0; j < len / 2; ++j)
          b[(size_t)(len - 1 - j)] = 3 - b[(size_t)j];
      if (kind == 3 && len > 1)
        b[rng() % (size_t)len] ^= 1;
      if (kind == 2 && (len & 1))
        b[(size_t)(len / 2)] = 0;  // odd length: no exact palindrome, keep it near one
      uint64_t hi = rng(), lo = rng();  // garbage above len must be ignored
      for (int j = 0; j < len; ++j) {
        hi = (hi << 1) | (uint64_t)(b[(size_t)j] >> 1);
        lo = (lo << 1) | (uint64_t)(b[(size_t)j] & 1);
      }
      const bool want = mod_rule::window_le_revcomp(len, [&](int j) { return b[(size_t)j]; });
      check(mod_rule::le_revcomp_planes(hi, lo, len) == want, "le_revcomp_planes == window_le_revcomp");
    }
  }
  for (uint64_t w = 1; w <= 40; ++w) {
    for (int it = 0; it < 2000; ++it) {
      std::vector<uint64_t> offs;
      const int n = 1 + (int)(rng() % 4);
      for (int i = 0; i < n; ++i)
        offs.push_back(rng() % w);
      mod_rule::CenterPick c;
      for (uint64_t p : offs)
        c.add(p, w);
      for (bool le : {false, true}) {
        uint64_t best = ~0ULL, lo = 0, hi = 0;
        for (uint64_t p : offs) {
          const uint64_t d = (uint64_t)std::llabs(2 * (long long)p - (long long)(w - 1));
          if (d < best) best = d, lo = hi = p;
          else if (d == best) lo = std::min(lo, p), hi = std::max(hi, p);
        }
        bool called = false;
        const uint64_t got = c.pick([&]() { called = true; return le; });
        check(got == (lo == hi || le ? lo : hi), "CenterPick picks the center-closest offset");
        check(called == (lo != hi), "v <= RC(v) is evaluated only on a two-way tie");
        check(lo == hi || lo + hi == w - 1, "two equally close offsets are mirror images");
      }
    }
  }
  check(mod_rule::t_is_canonical(15, 10, 5) && mod_rule::t_is_canonical(19, 10, 9) &&
            !mod_rule::t_is_canonical(15, 10, 6) && mod_rule::t_is_canonical(31, 5, 6),
        "t_is_canonical");
  std::cout << (failures ? "FAILED" : "OK") << "\n";
  return failures ? 1 : 0;
}
