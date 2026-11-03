// Expected density of a KMC3 signature scheme (see driver.h). Counts as KMC3's
// own n_super_kmers and total k-mers do (CKmerBinCollector::PutExtendedKmer:
// one superwindow per call, n - k + 1 k-mers): the loop below is
// CSplitter::ProcessReadsT (kmc/kmc_core/splitter.cpp) with each
// PutExtendedKmer replaced by the count, on KMC's own CMmer and signature
// orders (kmc/kmc_api/mmer.*, custom_order.*), so the tie rule (an equal
// signature is adopted, the last one wins on a rescan) and the cut of a
// superwindow at k + 255 bases are KMC's.
//
// Options: --k K (the k-mer length, kmc -k) --p P (the signature length, kmc -p)
// --scheme default|lex|xor|xor-pcg-hash|greedymini [--gm FILE] [--gm-lut binary|dna]
// [--projection at0-cg1|at1-cg0]; --seeds are kmc -os.

#include "driver.h"

#include "../../../kmc/kmc_api/custom_order.h"
#include "../../../kmc/kmc_api/mmer.h"

namespace {

typedef uint32_t uint32;

std::mutex g_order_mu;
bool g_installed = false;
uint64_t g_installed_seed = 0;

struct Kmc3Backend {
  uint32 kmer_len, signature_len;
  CCustomOrder::Config cfg;
  std::vector<char> codes;

  explicit Kmc3Backend(const expdens::Args &a)
      : kmer_len((uint32)a.need_i("k")), signature_len((uint32)a.need_i("p")) {
    const std::string s = a.need("scheme");
    cfg.signature_len = signature_len;
    if (s == "default") {
      cfg.mode = CCustomOrder::Mode::DEFAULT;
    } else if (s == "lex" || s == "xor" || s == "xor-pcg-hash") {
      cfg.mode = CCustomOrder::ParseMode(s == "lex" ? "xor" : s);
      cfg.label = s;
    } else if (s == "greedymini") {
      cfg.mode = CCustomOrder::Mode::TABLE;
      cfg.file = a.need("gm");
      cfg.lut_backend = a.str("gm-lut", "binary") == "dna" ? CCustomOrder::LutBackend::DnaExpanded
                                                           : CCustomOrder::LutBackend::BinaryRuntime;
      cfg.gm_at1 = a.str("projection", "at0-cg1") == "at1-cg0";
    } else {
      throw std::runtime_error("unknown kmc3 scheme " + s);
    }
    lex_ = s == "lex";
  }

  // KMC's signature order is process-wide (static CMmer tables): install it once per seed.
  void set_seed(uint64_t seed) {
    std::lock_guard<std::mutex> g(g_order_mu);
    if (g_installed && g_installed_seed == seed) return;
    cfg.seed = lex_ ? 0 : seed;
    CCustomOrder::Install(cfg);
    g_installed = true;
    g_installed_seed = seed;
  }

  void count(const char *text, size_t seq_size, uint64_t &n_kmers, uint64_t &n_selected) {
    codes.resize(seq_size);
    for (size_t j = 0; j < seq_size; j++) {
      switch (text[j]) {
        case 'A': codes[j] = 0; break;
        case 'C': codes[j] = 1; break;
        case 'G': codes[j] = 2; break;
        case 'T': codes[j] = 3; break;
        default: codes[j] = -1;
      }
    }
    switch (CMmer::LutMode()) {
      case CMmer::LUT_BINARY: return process<CMmer::LUT_BINARY>(codes.data(), (uint32)seq_size, n_kmers, n_selected);
      case CMmer::LUT_DNA: return process<CMmer::LUT_DNA>(codes.data(), (uint32)seq_size, n_kmers, n_selected);
      default: return process<CMmer::LUT_DEFAULT>(codes.data(), (uint32)seq_size, n_kmers, n_selected);
    }
  }

 private:
  bool lex_ = false;

  // CSplitter::ProcessReadsT for one read; put(len) stands for PutExtendedKmer(seq + ..., len).
  template <int Mode>
  void process(char *seq, uint32 seq_size, uint64_t &n_kmers, uint64_t &n_selected) {
    auto put = [&](uint32 n) {
      n_selected++;
      n_kmers += n - kmer_len + 1;
    };
    uint32 signature_start_pos;
    CMmer current_signature(signature_len), end_mmer(signature_len);
    uint32 i = 0;
    uint32 len = 0;  // length of extended kmer
    while (i + kmer_len - 1 < seq_size) {
      bool contains_N = false;
      for (uint32 j = 0; j < signature_len; ++j, ++i)
        if (seq[i] < 0) {
          contains_N = true;
          break;
        }
      if (contains_N) {
        ++i;
        continue;
      }
      len = signature_len;
      signature_start_pos = i - signature_len;
      current_signature.insert(seq + signature_start_pos);
      end_mmer.set(current_signature);
      for (; i < seq_size; ++i) {
        if (seq[i] < 0) {
          if (len >= kmer_len) put(len);
          len = 0;
          ++i;
          break;
        }
        end_mmer.insert_sym<Mode>(seq[i]);
        if (end_mmer < current_signature) {
          if (len >= kmer_len) {
            put(len);
            len = kmer_len - 1;
          }
          current_signature.set(end_mmer);
          signature_start_pos = i - signature_len + 1;
        } else if (end_mmer == current_signature) {
          current_signature.set(end_mmer);
          signature_start_pos = i - signature_len + 1;
        } else if (signature_start_pos + kmer_len - 1 < i) {
          put(len);
          len = kmer_len - 1;
          ++signature_start_pos;
          end_mmer.insert(seq + signature_start_pos);
          current_signature.set(end_mmer);
          for (uint32 j = signature_start_pos + signature_len; j <= i; ++j) {
            end_mmer.insert_sym<Mode>(seq[j]);
            if (end_mmer <= current_signature) {
              current_signature.set(end_mmer);
              signature_start_pos = j - signature_len + 1;
            }
          }
        }
        ++len;
        if (len == kmer_len + 255) {
          put(len);
          i -= kmer_len - 2;
          len = 0;
          break;
        }
      }
    }
    if (len >= kmer_len) put(len);
  }
};

}  // namespace

int main(int argc, char **argv) { return expdens::run<Kmc3Backend>(argc, argv); }
