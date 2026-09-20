/*
 * Correctness tests for GreedyMini UHS early-out and fresh pext extract.
 */

#include "greedymini_rank.h"
#include "minimizer_rankers.h"
#include "dna_kmer.h"
#include "utilities.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace kraken2;

static int failures = 0;

static void check(bool cond, const char *msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << "\n";
    failures++;
  }
}

int main() {
  // Closed-form masks match bit-set definition for a few (l, gm_k).
  for (int l = 8; l <= 31; l += 7) {
    for (int gk = 4; gk <= l && gk <= 15; gk += 5) {
      uint64_t expect_key = 0, expect_tb = 0;
      const int base0 = l - gk;
      for (int i = 0; i < gk; i++)
        expect_key |= 1ULL << (2 * (base0 + i) + 1);
      for (int i = 0; i < l; i++)
        expect_tb |= 1ULL << (2 * i);
      uint64_t expect_suf = 0;
      for (int i = 0; i < l - gk; i++)
        expect_suf |= 1ULL << (2 * i + 1);
      check(gm_dna_key_pext_mask(l, gk) == expect_key, "key mask closed form");
      check(gm_dna_tb_pext_mask(l) == expect_tb, "tb mask closed form");
      check(gm_dna_suffix_sigma_pext_mask(l, gk) == expect_suf,
            "suffix Σ mask closed form");
    }
  }

  // Synthetic table: gm_k=4 → 16 keys; UHS = {0,1,2,5} with ranks 0..3
  std::vector<uint64_t> raw(16, 16);
  raw[0] = 0;
  raw[1] = 1;
  raw[2] = 2;
  raw[5] = 3;

  GreedyMiniTable table(raw);
  check(table.gm_k() == 4, "gm_k");
  check(table.uhs_size() == 4, "uhs_size");
  check(table.in_uhs(0) && table.in_uhs(5), "in_uhs members");
  check(!table.in_uhs(3) && !table.in_uhs(15), "in_uhs non-members");

  const int l = 8, bpc = 2;
  const uint64_t lmer_mask = (1ULL << (2 * l)) - 1ULL;
  auto shared = std::make_shared<const GreedyMiniTable>(table);
  GMRanker ranker(l, bpc, /*toggle*/0, lmer_mask, true, shared);

  check(ranker.Rank(0) == 0, "all-zero finite rank");

  auto pack_orient = [&](uint64_t mer) -> uint64_t {
    uint32_t k = 0;
    uint64_t g = 0;
    gm_extract_dna(mer, l, 4, k, g);
    if (!table.in_uhs(k))
      return UINT64_MAX;
    const uint64_t atcg = gm_acgt_to_atcg(mer);
    const uint64_t suf = gm_extract_suffix_sigma_atcg(atcg, l, 4);
    return ((((uint64_t)table.primary(k) << (l - 4)) | suf) << l) | g;
  };
  auto expect_gm_orient = [&](uint64_t mer) -> uint64_t {
    const uint64_t rc =
        dna_reverse_complement(mer, (uint8_t)l, DNA_REVCOM_VERSION);
    const uint32_t kf = gm_extract_dna_key(mer, l, 4);
    const uint32_t kr = gm_extract_dna_key(rc, l, 4);
    const bool f_ok = table.in_uhs(kf);
    const bool r_ok = table.in_uhs(kr);
    if (!f_ok && !r_ok)
      return UINT64_MAX;
    if (!f_ok)
      return pack_orient(rc);
    if (!r_ok || rc == mer)
      return pack_orient(mer);
    const uint16_t pf = table.primary(kf);
    const uint16_t pr = table.primary(kr);
    if (pf < pr)
      return pack_orient(mer);
    if (pr < pf)
      return pack_orient(rc);
    const uint64_t rf = pack_orient(mer);
    const uint64_t rr = pack_orient(rc);
    return rf < rr ? rf : rr;
  };

  uint64_t lmer = 0;
  // ACGT bases A A G T … → ATCG MSBs 0,0,1,0 = key 2 (UHS rank 2)
  const int bases[8] = {0, 0, 2, 3, 0, 0, 0, 0};
  for (int i = 0; i < l; i++)
    lmer = (lmer << 2) | (uint64_t)bases[i];
  check(gm_extract_dna_key(lmer, l, 4) == 2, "extracted AT/CG key==2");
  check(ranker.Rank(lmer) == expect_gm_orient(lmer),
        "canonical gm_k then Extend-k pack");

  // Non-UHS on both strands → infinity
  uint64_t lmer_inf = 0;
  const int bases_inf[8] = {0, 1, 2, 3, 0, 0, 0, 0};  // ACGTAAAA
  for (int i = 0; i < l; i++)
    lmer_inf = (lmer_inf << 2) | (uint64_t)bases_inf[i];
  check(gm_extract_dna_key(lmer_inf, l, 4) == 6, "fwd key==6");
  check(ranker.Rank(lmer_inf) == expect_gm_orient(lmer_inf),
        "orientation matches probe rule");

  // gm_k primary beats DNA-lex: find mer where DNA-lex prefers RC but UHS
  // primary prefers fwd (or vice versa).
  {
    bool found = false;
    for (uint64_t mer = 0; mer < (1ULL << (2 * l)) && !found; mer++) {
      const uint64_t rc =
          dna_reverse_complement(mer, (uint8_t)l, DNA_REVCOM_VERSION);
      if (rc == mer)
        continue;
      const uint64_t dna_can = mer < rc ? mer : rc;
      const uint32_t k_can = gm_extract_dna_key(dna_can, l, 4);
      const uint32_t k_other =
          gm_extract_dna_key(dna_can == mer ? rc : mer, l, 4);
      if (!table.in_uhs(k_can) || !table.in_uhs(k_other))
        continue;
      if (table.primary(k_can) <= table.primary(k_other))
        continue;  // DNA-lex strand already has better/equal primary
      // DNA-lex strand has worse primary → Rank must follow other strand.
      const uint64_t got = ranker.Rank(dna_can);
      const uint64_t want = pack_orient(dna_can == mer ? rc : mer);
      check(got == want, "gm_k primary wins over DNA-lex strand");
      found = true;
    }
    check(found, "found DNA-lex vs gm_k conflict case");
  }

  // AND each strand independently; never Rank(RC(x & mask)).
  {
    check(MakeSpacedSeedMask(31, 7, 2) == 0x3ffffffff3333333ULL,
          "NT s=7 ℓ=31 expanded mask");
    const uint64_t mask = lmer_mask & ~0xFULL;
    GMRanker spaced(l, bpc, 0, lmer_mask, true, shared, mask);
    bool found_not_and_then_rc = false;
    for (uint64_t mer = 0; mer < (1ULL << (2 * l)); mer++) {
      if (spaced.Rank(mer) != spaced.Rank(mer & mask)) {
        found_not_and_then_rc = true;
        break;
      }
    }
    check(found_not_and_then_rc, "GM Rank is not AND-then-RC");
    auto pack_m = [&](uint64_t mer) -> uint64_t {
      return pack_orient(mer & mask);
    };
    auto expect_pair = [&](uint64_t mer) -> uint64_t {
      const uint64_t rc =
          dna_reverse_complement(mer, (uint8_t)l, DNA_REVCOM_VERSION);
      const uint64_t a = mer & mask;
      const uint64_t b = rc & mask;
      const uint32_t kf = gm_extract_dna_key(a, l, 4);
      const uint32_t kr = gm_extract_dna_key(b, l, 4);
      const bool f_ok = table.in_uhs(kf);
      const bool r_ok = table.in_uhs(kr);
      if (!f_ok && !r_ok)
        return UINT64_MAX;
      if (!f_ok)
        return pack_m(rc);
      if (!r_ok || rc == mer)
        return pack_m(mer);
      const uint16_t pf = table.primary(kf);
      const uint16_t pr = table.primary(kr);
      if (pf < pr)
        return pack_m(mer);
      if (pr < pf)
        return pack_m(rc);
      const uint64_t rf = pack_m(mer);
      const uint64_t rr = pack_m(rc);
      return rf < rr ? rf : rr;
    };
    check(spaced.Rank(lmer) == expect_pair(lmer),
          "spaced Rank matches per-strand AND");
  }

  // revcom_version=0 must drive GM strand RC, not hardcoded v1.
  {
    GMRanker r0(l, bpc, 0, lmer_mask, true, shared, 0, /*revcom*/0);
    GMRanker r1(l, bpc, 0, lmer_mask, true, shared, 0, /*revcom*/1);
    auto expect_at = [&](uint64_t mer, int ver) -> uint64_t {
      const uint64_t rc =
          dna_reverse_complement(mer, (uint8_t)l, ver);
      const uint32_t kf = gm_extract_dna_key(mer, l, 4);
      const uint32_t kr = gm_extract_dna_key(rc, l, 4);
      const bool f_ok = table.in_uhs(kf);
      const bool r_ok = table.in_uhs(kr);
      if (!f_ok && !r_ok)
        return UINT64_MAX;
      if (!f_ok)
        return pack_orient(rc);
      if (!r_ok || rc == mer)
        return pack_orient(mer);
      const uint16_t pf = table.primary(kf);
      const uint16_t pr = table.primary(kr);
      if (pf < pr)
        return pack_orient(mer);
      if (pr < pf)
        return pack_orient(rc);
      const uint64_t rf = pack_orient(mer);
      const uint64_t rr = pack_orient(rc);
      return rf < rr ? rf : rr;
    };
    bool found_rc_diff = false;
    bool found_rank_diff = false;
    for (uint64_t mer = 0; mer < (1ULL << (2 * l)); mer++) {
      const uint64_t rc0 = dna_reverse_complement(mer, (uint8_t)l, 0);
      const uint64_t rc1 = dna_reverse_complement(mer, (uint8_t)l, 1);
      if (rc0 != rc1)
        found_rc_diff = true;
      check(r0.Rank(mer) == expect_at(mer, 0), "v0 Rank matches v0 RC");
      check(r1.Rank(mer) == expect_at(mer, 1), "v1 Rank matches v1 RC");
      if (r0.Rank(mer) != r1.Rank(mer))
        found_rank_diff = true;
    }
    check(found_rc_diff, "v0 and v1 RC encodings differ");
    check(found_rank_diff, "GM Rank depends on revcom_version");
  }

  // ACGT→ATCG: A=00,C=01,G=10,T=11 → A=00,T=01,C=10,G=11
  check(gm_acgt_to_atcg(0x0) == 0x0, "remap A");
  check(gm_acgt_to_atcg(0x1) == 0x2, "remap C");
  check(gm_acgt_to_atcg(0x2) == 0x3, "remap G");
  check(gm_acgt_to_atcg(0x3) == 0x1, "remap T");
  // AT/CG key of RC equals reverse of forward key (canonical UHS property).
  {
    const int k = 4;
    uint64_t fwd = 0;
    const int fb[4] = {0, 1, 2, 3};  // ACGT
    for (int i = 0; i < k; i++)
      fwd = (fwd << 2) | (uint64_t)fb[i];
    // DNA RC of ACGT = ACGT (palindrome under RC for this encoding? A↔T,C↔G
    // ACGT → RC complement reverse: ACGT rev=TGCA, comp=ACGT — yes palindrome)
    // Use AAGT → RC = ACTT
    fwd = 0;
    const int fb2[4] = {0, 0, 2, 3};  // AAGT
    for (int i = 0; i < k; i++)
      fwd = (fwd << 2) | (uint64_t)fb2[i];
    // RC: rev T G A A → comp A C T T = 0,1,3,3
    uint64_t rc = 0;
    const int rb[4] = {0, 1, 3, 3};
    for (int i = 0; i < k; i++)
      rc = (rc << 2) | (uint64_t)rb[i];
    const uint32_t fk = gm_extract_dna_key(fwd, k, k);
    const uint32_t rk = gm_extract_dna_key(rc, k, k);
    uint32_t rev = 0;
    for (int i = 0; i < k; i++)
      rev = (rev << 1) | ((fk >> i) & 1u);
    check(rk == rev, "AT/CG key(RC) == reverse(key(fwd))");
  }

  GMRanker stream(l, bpc, 0, lmer_mask, true, shared);
  std::mt19937_64 rng(42);
  uint64_t lm = 0;
  for (int i = 0; i < l; i++)
    lm = ((lm << 2) | (rng() & 3ULL)) & lmer_mask;
  for (int step = 0; step < 20000; step++) {
    check(stream.Rank(lm) == expect_gm_orient(lm), "stream matches orient");
    lm = ((lm << 2) | (rng() & 3ULL)) & lmer_mask;
  }

  const uint64_t toggle = 0xe37e28c4271b5a2dULL & ((1ULL << l) - 1ULL);
  GMRanker tog(l, bpc, toggle, lmer_mask, true, shared);
  uint64_t x = 0;
  for (int i = 0; i < l; i++)
    x = ((x << 2) | (rng() & 3ULL)) & lmer_mask;
  uint64_t prev = tog.Rank(x);
  for (int step = 0; step < 5000; step++) {
    x = ((x << 2) | (rng() & 3ULL)) & lmer_mask;
    const uint64_t r = tog.Rank(x);
    prev = r;
  }
  check(prev == tog.Rank(x), "toggle rank deterministic");

  auto exported = gm_table_to_raw(table);
  check(exported == raw, "export matches input encoding");

  // ---- DNA alphabet (4^k table, contiguous key + suffix tiebreaker) ----
  {
    // gm_k=2 → 16 entries; UHS keys 0,1,5 with ranks 0,1,2; sentinel=16
    std::vector<uint64_t> dna_raw(16, 16);
    dna_raw[0] = 0;
    dna_raw[1] = 1;
    dna_raw[5] = 2;
    GreedyMiniTable dna_table(dna_raw, GmAlphabet::DNA);
    check(dna_table.alphabet() == GmAlphabet::DNA, "dna alphabet");
    check(dna_table.gm_k() == 2, "dna gm_k bases");
    check(dna_table.key_bits() == 4, "dna key_bits");
    check(dna_table.uhs_size() == 3, "dna uhs");

    // Loading DNA-sized table as binary would treat gm_k=log2(16)=4 — allowed
    // size-wise, but sentinel would be 16 either way. Require DNA flag for
    // 4^k semantics: loading with BINARY yields gm_k=4.
    GreedyMiniTable as_bin(dna_raw, GmAlphabet::BINARY);
    check(as_bin.gm_k() == 4, "same bytes as binary → gm_k=4");

    // Odd log2 size must reject DNA.
    bool threw = false;
    try {
      std::vector<uint64_t> odd(8, 8);
      odd[0] = 0;
      GreedyMiniTable bad(odd, GmAlphabet::DNA);
      (void)bad;
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    check(threw, "DNA rejects non-4^k size");

    const int dl = 6, dbpc = 2;
    const uint64_t dmask = (1ULL << (2 * dl)) - 1ULL;
    auto dna_shared = std::make_shared<const GreedyMiniTable>(dna_table);
    GMRanker dna_ranker(dl, dbpc, 0, dmask, true, dna_shared);

    // lmer ACGTAC: key=AC=1 → primary 1; tb=GTAC=0xb1; pack=(1<<8)|0xb1
    uint64_t dlmer = 0;
    const int dbases[6] = {0, 1, 2, 3, 0, 1};  // A C G T A C
    for (int i = 0; i < dl; i++)
      dlmer = (dlmer << 2) | (uint64_t)dbases[i];
    const uint32_t dkey = gm_extract_dna_alphabet_key(dlmer, dl, 2);
    const uint64_t dtb = gm_extract_dna_alphabet_tb(dlmer, dl, 2);
    check(dkey == 1, "dna alphabet key==1");
    check(dtb == 0xb1, "dna alphabet tb==GTAC");
    const uint64_t dr = dna_ranker.Rank(dlmer);
    check(dr == ((1ULL << 8) | 0xb1), "dna alphabet packs primary|tb");

    // Same primary, different suffix → different packed rank (tiebreak works).
    uint64_t dlmer_tb = 0;
    const int dbases_tb[6] = {0, 1, 0, 0, 0, 0};  // A C A A A A
    for (int i = 0; i < dl; i++)
      dlmer_tb = (dlmer_tb << 2) | (uint64_t)dbases_tb[i];
    const uint64_t dr_tb = dna_ranker.Rank(dlmer_tb);
    check(gm_extract_dna_alphabet_key(dlmer_tb, dl, 2) == 1, "same key");
    check(dr_tb == (1ULL << 8), "dna tb all-zero suffix");
    check(dr != dr_tb, "dna tiebreaker distinguishes same key");

    // first 2 bases GG (2,2) → key = 0b1010 = 10 → non-UHS
    uint64_t dlmer2 = 0;
    const int dbases2[6] = {2, 2, 0, 0, 0, 0};
    for (int i = 0; i < dl; i++)
      dlmer2 = (dlmer2 << 2) | (uint64_t)dbases2[i];
    check(gm_extract_dna_alphabet_key(dlmer2, dl, 2) == 10, "dna key==10");
    check(dna_ranker.Rank(dlmer2) == UINT64_MAX, "dna non-UHS → INF");

    check(dna_ranker.Rank(dlmer) == dna_ranker.Rank(
              dna_reverse_complement(dlmer, (uint8_t)dl, DNA_REVCOM_VERSION)),
          "dna alphabet Rank is strand-independent");

    // AA is DNA-lex smaller than TT, but give TT the better UHS primary.
    {
      std::vector<uint64_t> orient_raw(16, 16);
      orient_raw[0] = 2;   // AA
      orient_raw[15] = 0;  // TT best
      GreedyMiniTable ot(orient_raw, GmAlphabet::DNA);
      auto osh = std::make_shared<const GreedyMiniTable>(ot);
      GMRanker ork(dl, dbpc, 0, dmask, true, osh);
      uint64_t all_a = 0, all_t = 0;
      for (int i = 0; i < dl; i++) {
        all_a = (all_a << 2) | 0ULL;
        all_t = (all_t << 2) | 3ULL;
      }
      check(all_a < all_t, "AAAAAA DNA-lex-smaller than TTTTTT");
      check(ot.primary(0) > ot.primary(15), "AA worse primary than TT");
      const uint64_t want =
          ((uint64_t)ot.primary(15) << 8) |
          gm_extract_dna_alphabet_tb(all_t, dl, 2);
      check(ork.Rank(all_a) == want, "dna gm_k primary wins over DNA-lex strand");
      check(ork.Rank(all_a) == ork.Rank(all_t),
            "dna alphabet Rank(AAAAAA)==Rank(TTTTTT)");
    }
  }

  // Prefer env overrides; else ship-in-repo tables relative to cwd (src/ via make test).
  auto resolve_gm = [](const char *env_key, const char *rel_name) -> std::string {
    const char *env = std::getenv(env_key);
    if (env && *env)
      return env;
    const std::string candidates[] = {
        std::string("../../orders/GreedyMini/") + rel_name,  // cwd = kraken2/src (make test)
        std::string("../orders/GreedyMini/") + rel_name,     // cwd = kraken2/
        std::string("orders/GreedyMini/") + rel_name,        // cwd = repo root
    };
    for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); i++) {
      std::ifstream in(candidates[i], std::ios::binary);
      if (in)
        return candidates[i];
    }
    return {};
  };

  const std::string gm_path =
      resolve_gm("GM_FILE", "w5_k15_canonical.gm");
  if (!gm_path.empty()) {
    auto loaded = load_gm_table(gm_path, GmAlphabet::BINARY);
    if (loaded) {
      check(loaded->gm_k() == 15, "real gm_k");
      check(loaded->alphabet() == GmAlphabet::BINARY, "real binary alphabet");
      check(loaded->uhs_size() == 9076, "canonical UHS size");
      check(loaded->size() == 32768, "real table size");
      std::cout << "Validated w5_k15_canonical: UHS=" << loaded->uhs_size()
                << "/" << loaded->size() << "\n";

      // Monte Carlo density factor dens×(w+1): Extend-k keeps ℓ=31 ≈ ℓ=15.
      auto dens_factor = [&](int ell) -> double {
        const int w = 5;
        const uint64_t mask =
            (ell >= 32) ? ~0ULL : ((1ULL << (2 * ell)) - 1ULL);
        GMRanker r(ell, 2, 0, mask, true, loaded);
        const int N = 80000;
        std::vector<uint8_t> s((size_t)(N + ell));
        std::mt19937_64 g(7);
        for (auto &b : s)
          b = (uint8_t)(g() & 3ULL);
        auto pack = [&](int pos) -> uint64_t {
          uint64_t x = 0;
          for (int i = 0; i < ell; i++)
            x = (x << 2) | s[(size_t)(pos + i)];
          return x;
        };
        auto win_min = [&](int start) -> uint64_t {
          uint64_t best = UINT64_MAX;
          for (int i = 0; i < w; i++) {
            const uint64_t rk = r.Rank(pack(start + i));
            if (rk < best)
              best = rk;
          }
          return best;
        };
        uint64_t distinct = 0, total = 0;
        uint64_t last = UINT64_MAX;
        bool have = false;
        for (int i = 0; i + w <= N; i++) {
          const uint64_t m = win_min(i);
          total++;
          if (!have || m != last) {
            distinct++;
            last = m;
            have = true;
          }
        }
        const double dens = total ? (double)distinct / (double)total : 0.0;
        return dens * (double)(w + 1);
      };
      const double f15 = dens_factor(15);
      const double f31 = dens_factor(31);
      std::cout << "dens factor ℓ=15: " << f15 << "  ℓ=31: " << f31 << "\n";
      check(f15 > 1.35 && f15 < 1.65, "canonical dens factor ℓ=15 ~1.47");
      check(f31 > 1.35 && f31 < 1.70, "Extend-k dens factor ℓ=31 ~ℓ=15");
      check(f31 < f15 + 0.12, "ℓ=31 not much worse than ℓ=15");
    }
  }

  const std::string dna_gm_path =
      resolve_gm("GM_DNA_FILE", "w5_k7_improved_w5_dna.gm");
  if (!dna_gm_path.empty()) {
    auto dna_loaded = load_gm_table(dna_gm_path, GmAlphabet::DNA);
    if (dna_loaded) {
      check(dna_loaded->alphabet() == GmAlphabet::DNA, "real dna alphabet");
      check(dna_loaded->gm_k() == 7, "real dna gm_k");
      check(dna_loaded->size() == 16384, "real dna table size 4^7");
      check(dna_loaded->uhs_size() == 5248, "real dna UHS size");
      std::cout << "Validated w5_k7_dna: UHS=" << dna_loaded->uhs_size()
                << "/" << dna_loaded->size() << "\n";
    }
  }

  if (failures) {
    std::cerr << failures << " failure(s)\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
