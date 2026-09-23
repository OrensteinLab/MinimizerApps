// Monotone deque per residue class for streaming mod-sampling (Groot
// Koerkamp & Pibiri 2024): maintaining w = k-l+1 of these (one per t-mer
// locus residue i mod w) and combining their fronts gives the amortized-O(1)
// incremental analogue of mod-sampling's own O(k-t+1) brute-force scan.
// Shared by kraken2 (kraken2::ModSampler, src/mod_minimizer.h) and sshash
// (mod_minimizer_iterator, sshash/include/mod_minimizer_iterator.hpp); each
// app keeps its own tool-specific glue (bit-packing/ranking for kraken2,
// Kmer/hasher_type for sshash) around this shared primitive, since the
// call-site shapes differ.

#ifndef MINIMIZER_ORDERS_MOD_RESIDUE_DEQUE_H_
#define MINIMIZER_ORDERS_MOD_RESIDUE_DEQUE_H_

#include <cstdint>
#include <sys/types.h>  // ssize_t

struct ModTmerData {
  uint64_t rank;
  ssize_t pos;  // absolute t-mer start index since last reset
};

// Fixed-capacity monotone deque for one residue class (no heap allocs).
struct ModResidueDeque {
  static constexpr int kMaxCap = 16;

  ModTmerData data[kMaxCap];
  int head = 0;
  int size = 0;
  int cap = kMaxCap;

  void clear() { head = size = 0; }

  bool empty() const { return size == 0; }

  const ModTmerData &front() const { return data[head]; }

  const ModTmerData &back() const {
    return data[(head + size - 1) % cap];
  }

  void pop_front() {
    if (!size)
      return;
    head = (head + 1) % cap;
    --size;
  }

  void pop_back() {
    if (size)
      --size;
  }

  void push_back(const ModTmerData &v) {
    const int tail = (head + size) % cap;
    data[tail] = v;
    ++size;
  }

  void push_monotone(const ModTmerData &v) {
    while (!empty() && back().rank > v.rank)
      pop_back();
    push_back(v);
  }
};

#endif  // MINIMIZER_ORDERS_MOD_RESIDUE_DEQUE_H_
