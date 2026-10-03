/*
 * Correctness tests for GreedyMini UHS early-out and fresh pext extract.
 */

#include "greedymini_rank.h"
#include "../../common/minimizer-orders/greedymini_centre.h"
#include "minimizer_rankers.h"
#include "dna_kmer.h"
#include "kraken2_data.h"
#include "utilities.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
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

  // ---- Binary GreedyMini: l == gm_k (old behaviour) and centre extension ----
  {
    // Mirror-symmetric synthetic table, gm_k = 4: ranks are shared by a key
    // and its mirror: {0}:0 {1,8}:1 {2,4}:2 {5,10}:3.
    std::vector<uint64_t> sraw(16, 16);
    sraw[0] = 0;
    sraw[1] = sraw[8] = 1;
    sraw[2] = sraw[4] = 2;
    sraw[5] = sraw[10] = 3;
    GreedyMiniTable stable(sraw);
    check(stable.gm_k() == 4, "gm_k");
    check(stable.uhs_size() == 7, "uhs_size");
    check(stable.in_uhs(0) && stable.in_uhs(5), "in_uhs members");
    check(!stable.in_uhs(3) && !stable.in_uhs(15), "in_uhs non-members");
    check(gm_table_to_raw(stable) == sraw, "export matches input encoding");
    auto sshared = std::make_shared<const GreedyMiniTable>(stable);

    // l == gm_k: the pre-centre rule (orient by the better primary, then
    // (primary << l) | Gamma), with the spaced mask applied per strand.
    for (uint64_t toggle : {0ULL, 0x5ULL}) {
      for (uint64_t mask_bits : {~0ULL, ~0xCULL}) {
        const int lk = 4;
        const uint64_t lm = (1ULL << (2 * lk)) - 1ULL;
        const uint64_t mask = (mask_bits == ~0ULL) ? 0 : (lm & mask_bits);
        GMRanker r(lk, 2, toggle, lm, true, sshared, mask);
        check(!r.centre_extension(), "l == gm_k is not an extension");
        auto pack = [&](uint64_t w) -> uint64_t {
          uint32_t k;
          uint64_t g;
          gm_extract_dna(mask ? (w & mask) : w, lk, lk, k, g);
          if (!stable.in_uhs(k))
            return UINT64_MAX;
          return ((uint64_t)stable.primary(k) << lk) | (g ^ (toggle & lm));
        };
        bool all_ok = true;
        for (uint64_t mer = 0; mer < (1ULL << (2 * lk)); mer++) {
          const uint64_t rc =
              dna_reverse_complement(mer, (uint8_t)lk, DNA_REVCOM_VERSION);
          const uint64_t pa = pack(mer), pb = pack(rc);
          const uint16_t ka = stable.primary(
              gm_extract_dna_key(mask ? (mer & mask) : mer, lk, lk));
          const uint16_t kb = stable.primary(
              gm_extract_dna_key(mask ? (rc & mask) : rc, lk, lk));
          uint64_t want;
          if (ka != kb && ka != GreedyMiniTable::kInfRank &&
              kb != GreedyMiniTable::kInfRank)
            want = ka < kb ? pa : pb;
          else
            want = pa < pb ? pa : pb;  // inf == UINT64_MAX; tie -> smaller
          if (r.Rank(mer) != want)
            all_ok = false;
        }
        check(all_ok, "l == gm_k Rank equals orient-then-pack rule");
      }
    }
  }

  // Per-base test adapter: kraken2's ACGT=0123 words -> per-position Sigma /
  // Gamma bits (A,T -> Sigma 0; Gamma 1 for G,T), strand by strand. The spaced
  // mask (applied to each strand's own word) coarsens only the tie-break
  // fields; the rank lookup sees the unmasked Sigma. Fed to the shared slow
  // reference gm_centre::Reference.
  struct CentreRef {
    const GreedyMiniTable &table;
    int l, gk, i, ver;
    uint64_t toggle, mask;  // mask 0 = none
    struct P {
      std::vector<int> sig, sigm, gamm;
    };
    P planes(uint64_t w) const {
      P p;
      p.sig.assign((size_t)l, 0);
      p.sigm.assign((size_t)l, 0);
      p.gamm.assign((size_t)l, 0);
      const uint64_t wm = mask ? (w & mask) : w;
      for (int q = 0; q < l; q++) {
        const int b = (int)((w >> (2 * (l - 1 - q))) & 3);  // A C G T
        const int bm = (int)((wm >> (2 * (l - 1 - q))) & 3);
        p.sig[(size_t)q] = (b == 1 || b == 2);
        p.sigm[(size_t)q] = (bm == 1 || bm == 2);
        p.gamm[(size_t)q] = (bm >= 2);
      }
      return p;
    }
    gm_centre::Reference ref() const {
      return gm_centre::Reference(table, l, toggle);
    }
    uint64_t rank(uint64_t w) const {
      const uint64_t rc = dna_reverse_complement(w, (uint8_t)l, ver);
      const P a = planes(w), b = planes(rc);
      return ref().value(a.sig, a.gamm, b.sig, b.gamm, &a.sigm, &b.sigm);
    }
    bool tuple(uint64_t w, std::vector<int> &out) const {
      const uint64_t rc = dna_reverse_complement(w, (uint8_t)l, ver);
      const P a = planes(w), b = planes(rc);
      return ref().tuple(a.sig, a.gamm, b.sig, b.gamm, &out, &a.sigm, &b.sigm);
    }
  };

  auto rand_lmer = [](std::mt19937_64 &g, int ll, int bias) -> uint64_t {
    // bias 0: uniform; 1: AT-rich; 2: CG-rich; 3: low-entropy repeat;
    // 4: Sigma-palindromic (position p and l-1-p in the same AT/CG class,
    // so both strand frames tie on Sigma and only Gamma separates them)
    uint64_t x = 0;
    if (bias == 4) {
      std::vector<uint64_t> b((size_t)ll);
      for (int p = 0; p < ll; p++) {
        b[(size_t)p] = g() & 3ULL;
        if (p >= ll - 1 - p)
          b[(size_t)p] = (b[(size_t)(ll - 1 - p)] == 0 || b[(size_t)(ll - 1 - p)] == 3)
                             ? ((g() & 1ULL) ? 0 : 3)
                             : ((g() & 1ULL) ? 1 : 2);
      }
      for (int p = 0; p < ll; p++)
        x = (x << 2) | b[(size_t)p];
      return x;
    }
    if (bias == 3) {
      const int unit = 1 + (int)(g() % 4);
      std::vector<uint64_t> u(unit);
      for (auto &c : u)
        c = g() & 3ULL;
      for (int p = 0; p < ll; p++)
        x = (x << 2) | u[(size_t)(p % unit)];
      return x;
    }
    for (int p = 0; p < ll; p++) {
      uint64_t b = g() & 3ULL;  // 0 A,1 C,2 G,3 T
      if (bias == 1 && (g() % 4) != 0)
        b = (g() & 1ULL) ? 0 : 3;
      if (bias == 2 && (g() % 4) != 0)
        b = (g() & 1ULL) ? 1 : 2;
      x = (x << 2) | b;
    }
    return x;
  };
  auto lmask_of = [](int ll) { return (1ULL << (2 * ll)) - 1ULL; };

  // Centre extension on the synthetic table, exhaustively.
  {
    std::vector<uint64_t> sraw(16, 16);
    sraw[0] = 0;
    sraw[1] = sraw[8] = 1;
    sraw[2] = sraw[4] = 2;
    sraw[5] = sraw[10] = 3;
    auto sshared = std::make_shared<const GreedyMiniTable>(sraw);
    for (int ll : {6, 8, 10}) {
      for (uint64_t toggle : {0ULL, 0x2d3ULL}) {
        const uint64_t lm = lmask_of(ll);
        GMRanker r(ll, 2, toggle & ((1ULL << ll) - 1ULL), lm, true, sshared);
        check(r.centre_extension(), "l > gm_k uses the centre extension");
        CentreRef ref{*sshared, ll, 4, (ll - 4) / 2, DNA_REVCOM_VERSION,
                      toggle & ((1ULL << ll) - 1ULL), 0};
        bool eq = true, inv = true;
        std::vector<std::pair<uint64_t, uint64_t>> fin;  // (rank, canonical)
        for (uint64_t mer = 0; mer < (1ULL << (2 * ll)); mer++) {
          const uint64_t rc =
              dna_reverse_complement(mer, (uint8_t)ll, DNA_REVCOM_VERSION);
          const uint64_t got = r.Rank(mer);
          if (got != ref.rank(mer))
            eq = false;
          if (got != r.Rank(rc))
            inv = false;
          if (got != UINT64_MAX && mer <= rc)
            fin.push_back({got, mer});
        }
        check(eq, "centre Rank == per-base reference (exhaustive)");
        check(inv, "centre Rank(x) == Rank(rc(x)) (exhaustive)");
        std::sort(fin.begin(), fin.end());
        bool inj = true;
        for (size_t j = 1; j < fin.size(); j++)
          if (fin[j].first == fin[j - 1].first)
            inj = false;
        check(inj && !fin.empty(), "centre injective on canonical (exh.)");
      }
    }
    // Spaced seed: unmasked lookup, masked tie-break fields; strand invariant.
    {
      const int ll = 8;
      const uint64_t lm = lmask_of(ll);
      const uint64_t mask = lm & ~0xF0ULL;
      GMRanker r(ll, 2, 0, lm, true, sshared, mask);
      CentreRef ref{*sshared, ll, 4, 2, DNA_REVCOM_VERSION, 0, mask};
      bool eq = true, inv = true;
      for (uint64_t mer = 0; mer < (1ULL << (2 * ll)); mer++) {
        const uint64_t rc =
            dna_reverse_complement(mer, (uint8_t)ll, DNA_REVCOM_VERSION);
        if (r.Rank(mer) != ref.rank(mer))
          eq = false;
        if (r.Rank(mer) != r.Rank(rc))
          inv = false;
      }
      check(eq, "spaced centre Rank == reference (mask only in tie-break)");
      check(inv, "spaced centre Rank strand invariant");
    }
    // revcom_version 0 drives the strand RC, not a hardcoded v1.
    {
      const int ll = 8;
      const uint64_t lm = lmask_of(ll);
      GMRanker r0(ll, 2, 0, lm, true, sshared, 0, 0);
      GMRanker r1(ll, 2, 0, lm, true, sshared, 0, 1);
      CentreRef ref0{*sshared, ll, 4, 2, 0, 0, 0};
      CentreRef ref1{*sshared, ll, 4, 2, 1, 0, 0};
      bool m0 = true, m1 = true, diff = false;
      for (uint64_t mer = 0; mer < (1ULL << (2 * ll)); mer++) {
        m0 = m0 && r0.Rank(mer) == ref0.rank(mer);
        m1 = m1 && r1.Rank(mer) == ref1.rank(mer);
        diff = diff || r0.Rank(mer) != r1.Rank(mer);
      }
      check(m0, "v0 Rank matches v0 RC");
      check(m1, "v1 Rank matches v1 RC");
      check(diff, "GM Rank depends on revcom_version");
    }

    // Rejections.
    auto throws = [](const std::function<void()> &f) {
      try {
        f();
      } catch (const std::invalid_argument &) {
        return true;
      }
      return false;
    };
    check(throws([&] { GMRanker(7, 2, 0, lmask_of(7), true, sshared); }),
          "odd l - gm_k rejected");
    check(throws([&] { GMRanker(9, 2, 0, lmask_of(9), true, sshared); }),
          "odd l - gm_k rejected (9)");
    check(throws([&] { GMRanker(3, 2, 0, lmask_of(3), true, sshared); }),
          "gm_k > l rejected");
    {
      std::vector<uint64_t> asym(16, 16);  // not mirror symmetric
      asym[0] = 0;
      asym[1] = 1;
      asym[2] = 2;
      asym[5] = 3;
      auto ash = std::make_shared<const GreedyMiniTable>(asym);
      check(throws([&] { GMRanker(8, 2, 0, lmask_of(8), true, ash); }),
            "asymmetric table rejected for the centre extension");
      // ... but l == gm_k does not need symmetry (unchanged behaviour)
      check(!throws([&] { GMRanker(4, 2, 0, lmask_of(4), true, ash); }),
            "asymmetric table still accepted at l == gm_k");
    }
    {
      std::vector<uint64_t> shared_rank(16, 16);  // rank shared by non-mirrors
      shared_rank[1] = shared_rank[8] = 0;
      shared_rank[2] = shared_rank[4] = 0;
      auto sh2 = std::make_shared<const GreedyMiniTable>(shared_rank);
      check(throws([&] { GMRanker(8, 2, 0, lmask_of(8), true, sh2); }),
            "rank shared by non-mirror keys rejected");
    }
  }

  // ---- shared header greedymini_centre.h: helpers ----
  {
    check(gm_centre::reverse_low_bits(0b0001101ULL, 7) == 0b1011000ULL,
          "reverse_low_bits");
    check(gm_centre::mirror_key(0b0011, 4) == 0b1100, "mirror_key");
    check(gm_centre::centre_key(0b101101ULL, 1, 4) == 0b0110,
          "centre_key picks bits [i, i+gm_k) from the MSB");
    bool t1 = false, t2 = false, ok3 = true;
    try { gm_centre::check_length(8, 4); } catch (...) { ok3 = false; }
    try { gm_centre::check_length(7, 4); } catch (const std::invalid_argument &) { t1 = true; }
    try { gm_centre::check_length(3, 4); } catch (const std::invalid_argument &) { t2 = true; }
    check(t1 && t2 && ok3, "check_length");
    // value_of_planes (rc planes by bit reversal) == value_of_two_strands.
    std::vector<uint64_t> sraw(16, 16);
    sraw[0] = 0;
    sraw[1] = sraw[8] = 1;
    sraw[2] = sraw[4] = 2;
    sraw[5] = sraw[10] = 3;
    auto st = std::make_shared<const GreedyMiniTable>(sraw);
    const int ll = 10;
    gm_centre::Order ord(st, ll, 0x2a5);
    check(ord.value_bits() == 2 + 2 * 3 + 1 + ll,
          "value_bits");
    std::mt19937_64 g(3);
    bool same = true;
    for (int it = 0; it < 200000; it++) {
      const uint64_t x = g() & lmask_of(ll);
      const uint64_t rc = dna_reverse_complement(x, (uint8_t)ll, 1);
      const uint64_t a = gm_acgt_to_atcg(x), b = gm_acgt_to_atcg(rc);
      const uint64_t sa = gm_extract_dna_key_atcg(a, ll, ll);
      const uint64_t ga = gm_extract_dna_tb_atcg(a, ll);
      const uint64_t sb = gm_extract_dna_key_atcg(b, ll, ll);
      const uint64_t gb = gm_extract_dna_tb_atcg(b, ll);
      same = same && ord.value_of_planes(sa, ga) ==
                         ord.value_of_two_strands(sa, ga, sb, gb) &&
             ord.value_of_planes(sa, ga) == ord.value_of_planes(sb, gb);
    }
    check(same, "value_of_planes == two-strand value, strand invariant");
    // density_factor vs a brute-force leftmost-minimum count.
    std::vector<uint64_t> rk(5000);
    for (auto &v : rk)
      v = (g() % 7 == 0) ? UINT64_MAX : (g() % 40);
    for (int w : {1, 3, 5}) {
      uint64_t sel = 0;
      long long last = -1;
      for (size_t i = 0; i + (size_t)w <= rk.size(); i++) {
        size_t best = i;
        for (int j = 1; j < w; j++)
          if (rk[i + (size_t)j] < rk[best])
            best = i + (size_t)j;
        sel += (long long)best != last;
        last = (long long)best;
      }
      const double want = (double)sel / (double)(rk.size() - (size_t)w + 1) * (w + 1);
      check(std::abs(gm_centre::density_factor(rk, w) - want) < 1e-12,
            "density_factor == brute force");
    }
  }

  // DB guard: version-1 Minimizer+binary GM databases with l > gm_k are legacy.
  {
    IndexOptions o = {0};
    o.l = 31;
    o.k = 35;
    o.dna_db = true;
    o.gm_k = 15;
    o.gm_alphabet = 0;
    EncodeIndexSchemes(o, SketchingScheme::MINIMIZER, Order::GM);
    check(o.db_version == KRAKEN2_DB_VERSION &&
              KRAKEN2_DB_VERSION >= KRAKEN2_DB_VERSION_GM_CENTRE,
          "new databases are written with the centre-extension version");
    check(!IndexUsesLegacyGmPrefix(o), "current GM database is not legacy");
    o.db_version = KRAKEN2_DB_VERSION_SCHEMES;
    check(IndexUsesLegacyGmPrefix(o), "v1 GM minimizer db, l > gm_k: legacy");
    o.l = 15;
    check(!IndexUsesLegacyGmPrefix(o), "v1 GM db with l == gm_k: fine");
    o.l = 31;
    o.sketching_scheme = (uint8_t)SketchingScheme::MOD;
    o.gm_k = 6;
    check(!IndexUsesLegacyGmPrefix(o), "v1 Mod+GM db is unaffected");
    o.sketching_scheme = (uint8_t)SketchingScheme::MINIMIZER;
    o.reserved = (uint32_t)Order::RANDOM_XOR;
    check(!IndexUsesLegacyGmPrefix(o), "v1 random-order db is unaffected");
    o.db_version = KRAKEN2_DB_VERSION_ORIGINAL;
    check(!IndexUsesLegacyGmPrefix(o), "v0 db is unaffected");
  }

  // ACGT->ATCG: A=00,C=01,G=10,T=11 -> A=00,T=01,C=10,G=11
  check(gm_acgt_to_atcg(0x0) == 0x0, "remap A");
  check(gm_acgt_to_atcg(0x1) == 0x2, "remap C");
  check(gm_acgt_to_atcg(0x2) == 0x3, "remap G");
  check(gm_acgt_to_atcg(0x3) == 0x1, "remap T");
  // AT/CG key of RC equals reverse of forward key (canonical UHS property).
  {
    const int k = 4;
    const int fb2[4] = {0, 0, 2, 3};  // AAGT
    const int rb[4] = {0, 1, 3, 3};   // RC = ACTT
    uint64_t fwd = 0, rc = 0;
    for (int i = 0; i < k; i++) {
      fwd = (fwd << 2) | (uint64_t)fb2[i];
      rc = (rc << 2) | (uint64_t)rb[i];
    }
    const uint32_t fk = gm_extract_dna_key(fwd, k, k);
    const uint32_t rk = gm_extract_dna_key(rc, k, k);
    uint32_t rev = 0;
    for (int i = 0; i < k; i++)
      rev = (rev << 1) | ((fk >> i) & 1u);
    check(rk == rev, "AT/CG key(RC) == reverse(key(fwd))");
  }

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

      // --- centre extension on the real table, l = 17 .. 31 ---
      for (int ll : {17, 21, 25, 31}) {
        const uint64_t lm = lmask_of(ll);
        const uint64_t tog = 0xe37e28c4271b5a2dULL & ((1ULL << ll) - 1ULL);
        for (uint64_t mask : {(uint64_t)0, (uint64_t)MakeSpacedSeedMask(ll, ll / 4, 2)}) {
          GMRanker r(ll, 2, tog, lm, true, loaded, mask);
          CentreRef ref{*loaded, ll, 15, (ll - 15) / 2, DNA_REVCOM_VERSION,
                        tog, mask};
          std::mt19937_64 g(1000 + (uint64_t)ll);
          const int n_iter = (mask == 0 && ll == 31) ? 10000000 : 1500000;
          uint64_t bad_ref = 0, bad_inv = 0, finite = 0;
          for (int it = 0; it < n_iter; it++) {
            const uint64_t x = rand_lmer(g, ll, it % 5);
            const uint64_t rc =
                dna_reverse_complement(x, (uint8_t)ll, DNA_REVCOM_VERSION);
            const uint64_t got = r.Rank(x);
            if (got != ref.rank(x))
              bad_ref++;
            if (got != r.Rank(rc))
              bad_inv++;
            finite += got != UINT64_MAX;
          }
          check(bad_ref == 0, "real table: centre Rank == per-base reference");
          check(bad_inv == 0, "real table: Rank(x) == Rank(rc(x))");
          check(finite > 0, "real table: some finite ranks");
        }
        // Pairwise order == lexicographic (rank, Sigma(xyz), Gamma^toggle),
        // on l-mers sharing a centre as well (the interesting ties).
        {
          GMRanker r(ll, 2, tog, lm, true, loaded);
          CentreRef ref{*loaded, ll, 15, (ll - 15) / 2, DNA_REVCOM_VERSION,
                        tog, 0};
          std::mt19937_64 g(77 + (uint64_t)ll);
          const int i = (ll - 15) / 2;
          bool order_ok = true;
          for (int it = 0; it < 200000 && order_ok; it++) {
            uint64_t a = rand_lmer(g, ll, 0), b = rand_lmer(g, ll, 0);
            if (it & 1) {  // same centre y, different x / z
              const int sh = 2 * (ll - i - 15);
              const uint64_t cmask = ((1ULL << 30) - 1ULL) << sh;
              b = (b & ~cmask) | (a & cmask);
            }
            std::vector<int> ta, tb;
            const bool fa = ref.tuple(a, ta), fb = ref.tuple(b, tb);
            const uint64_t ra = r.Rank(a), rb = r.Rank(b);
            if (!fa || !fb) {
              order_ok = (fa == (ra != UINT64_MAX)) && (fb == (rb != UINT64_MAX));
              continue;
            }
            order_ok = (ra < rb) == (ta < tb) && (ra == rb) == (ta == tb);
          }
          check(order_ok, "real table: order == (rank, Sigma, Gamma) lex");
        }
        // Injectivity on canonical l-mers (random + shared-centre families).
        {
          GMRanker r(ll, 2, 0, lm, true, loaded);
          std::mt19937_64 g(5 + (uint64_t)ll);
          std::vector<std::pair<uint64_t, uint64_t>> fin;
          const int i = (ll - 15) / 2;
          const int sh = 2 * (ll - i - 15);
          const uint64_t cmask = ((1ULL << 30) - 1ULL) << sh;
          uint64_t centre = 0;
          for (int it = 0; it < 1500000; it++) {
            uint64_t x = rand_lmer(g, ll, it % 5);
            if ((it & 7) < 4) {
              if ((it & 7) == 0)
                centre = x;
              x = (x & ~cmask) | (centre & cmask);  // shared centre
            }
            const uint64_t rc =
                dna_reverse_complement(x, (uint8_t)ll, DNA_REVCOM_VERSION);
            const uint64_t rk = r.Rank(x);
            if (rk != UINT64_MAX)
              fin.push_back({rk, x < rc ? x : rc});
          }
          std::sort(fin.begin(), fin.end());
          fin.erase(std::unique(fin.begin(), fin.end()), fin.end());
          bool inj = true;
          for (size_t j = 1; j < fin.size(); j++)
            if (fin[j].first == fin[j - 1].first)
              inj = false;
          check(inj && fin.size() > 200000,
                "real table: injective on canonical l-mers");
        }
      }
      // l == gm_k on the real table: orient-then-pack, as before.
      {
        const int lk = 15;
        GMRanker r(lk, 2, 0, lmask_of(lk), true, loaded);
        std::mt19937_64 g(15);
        bool ok = true;
        for (int it = 0; it < 1000000 && ok; it++) {
          const uint64_t x = rand_lmer(g, lk, it % 5);
          const uint64_t rc =
              dna_reverse_complement(x, (uint8_t)lk, DNA_REVCOM_VERSION);
          uint32_t ka, kb;
          uint64_t ga, gb;
          gm_extract_dna(x, lk, lk, ka, ga);
          gm_extract_dna(rc, lk, lk, kb, gb);
          const uint16_t pa = loaded->primary(ka), pb = loaded->primary(kb);
          uint64_t want = UINT64_MAX;
          if (pa != GreedyMiniTable::kInfRank)
            want = ((uint64_t)pa << lk) | ga;
          if (pb != GreedyMiniTable::kInfRank) {
            const uint64_t wb = ((uint64_t)pb << lk) | gb;
            if (wb < want)
              want = wb;
          }
          ok = r.Rank(x) == want;
        }
        check(ok, "real table: l == gm_k Rank equals old orient+pack");
      }
      // Invalid lengths for gm_k = 15 are rejected with a clear error.
      for (int ll : {16, 20, 30}) {
        bool threw = false;
        try {
          GMRanker r(ll, 2, 0, lmask_of(ll), true, loaded);
        } catch (const std::invalid_argument &ex) {
          threw = std::string(ex.what()).find("even") != std::string::npos;
        }
        check(threw, "l - gm_k odd rejected for the default table");
      }

      // Monte Carlo density factor dens x (w+1) (one stream; the full study
      // is benchmark/kraken2/density_sweep.sh): the centre extension keeps the
      // design value (~1.47) for every valid l, not only l == gm_k.
      auto dens_factor = [&](int ell) -> double {
        const int w = 5;
        GMRanker r(ell, 2, 0, lmask_of(ell), true, loaded);
        const int N = 300000;
        std::mt19937_64 g(7);
        std::vector<uint64_t> rk((size_t)N + w - 1);
        uint64_t x = 0;
        for (int i = 0; i < ell - 1; i++)
          x = ((x << 2) | (g() & 3ULL)) & lmask_of(ell);
        for (auto &v : rk) {
          x = ((x << 2) | (g() & 3ULL)) & lmask_of(ell);
          v = r.Rank(x);
        }
        return gm_centre::density_factor(rk, w);
      };
      const double f15 = dens_factor(15);
      const double f31 = dens_factor(31);
      std::cout << "dens factor l=15: " << f15 << "  l=31: " << f31 << "\n";
      check(f15 > 1.40 && f15 < 1.54, "canonical dens factor l=15 ~1.47");
      check(f31 > 1.40 && f31 < 1.54, "centre dens factor l=31 ~ l=15");
      check(std::abs(f31 - f15) < 0.04, "l=31 matches l=15 (centre extension)");
    }
  }

  // ---- injective total-order mode (sshash's order) vs the reference formula
  // key(x) = (rank_or_inf, min(Sigma_x << l | Gamma_x ^ t,
  //                            Sigma_rc << l | Gamma_rc ^ t)),
  // l = 21, both Sigma polarities, random toggles, tables w11_k13 / w11_k15.
  for (const char *tname : {"w11_k13_canonical.gm", "w11_k15_canonical.gm"}) {
    const std::string tpath = resolve_gm("GM_FILE_INJ", tname);
    if (tpath.empty())
      continue;
    auto tbl = load_gm_table(tpath, GmAlphabet::BINARY);
    check(tbl != nullptr, "injective test table loads");
    if (!tbl)
      continue;
    const int ll = 21;
    const int gk = tbl->gm_k();
    const uint64_t smask = (1ULL << ll) - 1ULL;  // l-bit plane mask
    // Per-base planes of a word (A,T -> Sigma 0 unless flipped; Gamma 1 for G,T).
    auto per_base = [&](uint64_t w, bool flip, uint64_t &S, uint64_t &G) {
      S = G = 0;
      for (int q = 0; q < ll; q++) {
        const int b = (int)((w >> (2 * (ll - 1 - q))) & 3);  // A C G T
        S = (S << 1) | (uint64_t)((b == 1 || b == 2) ^ (int)flip);
        G = (G << 1) | (uint64_t)(b >= 2);
      }
    };
    std::mt19937_64 g(900 + (uint64_t)gk);
    uint64_t pairs = 0, bad_bits = 0, bad_order = 0, bad_eq = 0, bad_two = 0,
             inf_pairs = 0, tie_rank = 0, pal_keys = 0;
    for (int flip = 0; flip < 2; flip++) {
      for (int tg = 0; tg < 3; tg++) {
        const uint64_t toggle = g() & smask;
        gm_centre::Order ord(tbl, ll, toggle, /*injective=*/true);
        auto ref_key = [&](uint64_t w, uint64_t &fullval, bool &pal) {
          const uint64_t rc =
              dna_reverse_complement(w, (uint8_t)ll, DNA_REVCOM_VERSION);
          uint64_t S, G, Sr, Gr;
          per_base(w, flip, S, G);
          per_base(rc, flip, Sr, Gr);
          const uint32_t key = (uint32_t)(S >> ((ll - gk) / 2)) & ((1u << gk) - 1u);
          const uint64_t a = (S << ll) | (G ^ toggle);
          const uint64_t b = (Sr << ll) | (Gr ^ toggle);
          pal = gm_centre::mirror_key(key, gk) == key;
          fullval = ((uint64_t)tbl->primary(key) << (2 * ll)) | (a < b ? a : b);
          // the rank_or_inf, min tuple is exactly this packed integer
          return (uint64_t)tbl->primary(key);
        };
        auto order_val = [&](uint64_t w) {
          uint64_t S, G;
          per_base(w, flip, S, G);  // caller handles polarity
          return ord.value_of_planes(S, G);
        };
        for (int it = 0; it < 1700000; it++) {
          uint64_t a = rand_lmer(g, ll, it % 5), b = rand_lmer(g, ll, (it / 5) % 5);
          if (it & 1) {  // same centre, different tails -> equal ranks
            const int i = (ll - gk) / 2;
            const uint64_t cmask = ((1ULL << (2 * gk)) - 1ULL) << (2 * (ll - i - gk));
            b = (b & ~cmask) | (a & cmask);
          }
          uint64_t fa, fb;
          bool pa, pb;
          const uint64_t ra = ref_key(a, fa, pa), rb = ref_key(b, fb, pb);
          const uint64_t va = order_val(a), vb = order_val(b);
          pairs++;
          pal_keys += pa;
          inf_pairs += (ra == GreedyMiniTable::kInfRank);
          tie_rank += (ra == rb);
          if (va != fa || vb != fb)
            bad_bits++;
          if ((va < vb) != (fa < fb))
            bad_order++;
          if ((va == vb) != (fa == fb))
            bad_eq++;
          if ((va >> (2 * ll)) != ra)
            bad_bits++;
          // rank alone, and strand-by-strand with both frames given
          const uint64_t rc =
              dna_reverse_complement(a, (uint8_t)ll, DNA_REVCOM_VERSION);
          uint64_t S, G, Sr, Gr;
          per_base(a, flip, S, G);
          per_base(rc, flip, Sr, Gr);
          if (ord.rank_of_sigma(S) != ra ||
              ord.value_of_two_strands(S, G, Sr, Gr) != va)
            bad_two++;
        }
      }
    }
    std::cout << "injective " << tname << ": " << pairs << " pairs ("
              << inf_pairs << " infinite, " << tie_rank << " equal-rank, "
              << pal_keys << " palindromic-key)\n";
    check(pairs >= 10000000, "injective test pair count");
    check(inf_pairs > 100000 && tie_rank > 100000 && pal_keys > 1000,
          "injective test covers infinite / tied / palindromic windows");
    check(bad_bits == 0, "injective value == (rank << 2l) | min-frame field");
    check(bad_order == 0, "injective order == reference formula");
    check(bad_eq == 0, "injective equality == reference equality");
    check(bad_two == 0, "injective two-strand value / rank_of_sigma");

    // rank_at: window of a longer Sigma plane.
    {
      gm_centre::Order ord(tbl, ll, 0, true);
      bool ok = true;
      for (int it = 0; it < 100000 && ok; it++) {
        const uint64_t P = g() & ((1ULL << 40) - 1ULL);
        for (int sh = 0; sh <= 40 - ll; sh++)
          ok = ok && ord.rank_at(P, sh) == ord.rank_of_sigma((P >> sh) & smask);
      }
      check(ok, "rank_at == rank_of_sigma of the shifted window");
    }
    // Injective on canonical l-mers (including infinite windows).
    {
      gm_centre::Order ord(tbl, ll, 0x15a5a5, true);
      std::vector<std::pair<uint64_t, uint64_t>> v;
      for (int it = 0; it < 1500000; it++) {
        const uint64_t x = rand_lmer(g, ll, it % 5);
        const uint64_t rc =
            dna_reverse_complement(x, (uint8_t)ll, DNA_REVCOM_VERSION);
        uint64_t S, G;
        per_base(x, false, S, G);
        v.push_back({ord.value_of_planes(S, G), x < rc ? x : rc});
      }
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      bool inj = true;
      for (size_t j = 1; j < v.size(); j++)
        inj = inj && v[j].first != v[j - 1].first;
      check(inj, "injective mode: injective on canonical l-mers (inf too)");
    }
    // l too large for the injective packing; compressed packing still fits.
    {
      bool t = false;
      try {
        gm_centre::Order bad(tbl, 31, 0, true);
      } catch (const std::invalid_argument &) {
        t = true;
      }
      check(t, "injective mode rejects 16 + 2l > 64");
      gm_centre::Order ok(tbl, 31, 0, false);
      check(ok.value_bits() <= 64, "compressed packing fits l = 31");
    }
  }

  // ---- injective mode == sshash's hash: the 384 Python-reference vectors of
  // sshash/test/gm_centre_vectors.inc (sshash packs A=0,C=1,G=3,T=2; seed 1234;
  // cfg 0/1 = w11_k13 at0-cg1/at1-cg0, cfg 2/3 = w11_k15). Read-only include.
  {
    static const struct { int cfg; uint64_t x, h; } vectors[] = {
#include "../../sshash/test/gm_centre_vectors.inc"
    };
    const char *tn[2] = {"w11_k13_canonical.gm", "w11_k15_canonical.gm"};
    const int ll = 21;
    const uint64_t smask = (1ULL << ll) - 1ULL;
    std::shared_ptr<const GreedyMiniTable> tb[2];
    for (int t = 0; t < 2; t++) {
      const std::string tp = resolve_gm("GM_FILE_INJ", tn[t]);
      if (!tp.empty())
        tb[t] = load_gm_table(tp, GmAlphabet::BINARY);
    }
    if (tb[0] && tb[1]) {
      gm_centre::Order ords[4] = {
          gm_centre::Order(tb[0], ll, 1234 & smask, true),
          gm_centre::Order(tb[0], ll, 1234 & smask, true),
          gm_centre::Order(tb[1], ll, 1234 & smask, true),
          gm_centre::Order(tb[1], ll, 1234 & smask, true)};
      uint64_t n = 0, bad = 0;
      for (const auto &v : vectors) {
        uint64_t S = 0, G = 0;
        for (int q = 0; q < ll; q++) {
          const int sb = (int)((v.x >> (2 * (ll - 1 - q))) & 3);  // A C T G
          const int b = (sb == 3) ? 2 : (sb == 2) ? 3 : sb;       // -> A C G T
          S = (S << 1) | (uint64_t)(((b == 1 || b == 2) ? 1 : 0) ^ (v.cfg & 1));
          G = (G << 1) | (uint64_t)(b >= 2);
        }
        bad += ords[v.cfg].value_of_planes(S, G) != v.h;
        n++;
      }
      std::cout << n << " sshash Python-reference vectors checked\n";
      check(n == 384 && bad == 0, "injective Order == sshash Python vectors");
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
