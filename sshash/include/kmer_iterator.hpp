#pragma once

#include "util.hpp"

namespace sshash {

template <typename Kmer, typename BitVector>
struct kmer_iterator  //
{
    kmer_iterator() {}

    kmer_iterator(BitVector const& bv, const uint64_t k)
        : m_bv(&bv), m_uint_kmer_bits(Kmer::bits_per_char * k), m_pos(0), m_avail(0), m_buff(0) {}

    kmer_iterator(BitVector const& bv, const uint64_t k, const uint64_t pos)
        : kmer_iterator(bv, k)  //
    {
        at(pos);
    }

    void at(uint64_t pos) {
        m_pos = pos;
        m_avail = 0;
        m_buff = 0;
    }

    Kmer get() {
        if (m_avail < m_uint_kmer_bits) fill_buff();
        auto kmer = m_buff;
        kmer.take(m_uint_kmer_bits);
        return kmer;
    }

    Kmer get_reverse() {
        if (m_avail < m_uint_kmer_bits) fill_buff_reverse();
        auto kmer = m_buff;
        kmer.drop(Kmer::uint_kmer_bits - m_uint_kmer_bits);
        return kmer;
    }

    void next() {
        if (m_avail < Kmer::bits_per_char) fill_buff();
        m_buff.drop_char();
        m_avail -= Kmer::bits_per_char;
        m_pos += Kmer::bits_per_char;
    }

    void next_reverse() {
        if (m_avail < Kmer::bits_per_char) fill_buff_reverse();
        m_buff.pad_char();
        m_avail -= Kmer::bits_per_char;
        m_pos -= Kmer::bits_per_char;
    }

    inline uint64_t get_next_char() {
        if (m_avail < Kmer::bits_per_char) fill_buff_inlined();
        m_avail -= Kmer::bits_per_char;
        m_pos += Kmer::bits_per_char;
        return m_buff.pop_char();
    }

    inline uint64_t position() const { return m_pos; }

private:
    /* The refill itself, always inlined into `get_next_char`, the per-k-mer step of the string
       iterators (`sshash bench`'s iterator, check): there, in a translation unit as big as
       tools/sshash.cpp, the compiler's size heuristics otherwise leave the call out of line and
       the iteration gets 70% slower. The other callers (the streaming query's extension steps)
       use `fill_buff`, left to the compiler: forcing it there makes the streaming lookup 8%
       slower. */
    __attribute__((always_inline)) inline void fill_buff_inlined() {
        static_assert(Kmer::uint_kmer_bits % 64 == 0);
        for (int i = Kmer::uint_kmer_bits - 64; i >= 0; i -= 64) {
            m_buff.append64(m_bv->get_word64(m_pos + i));
        }
        m_avail = Kmer::uint_kmer_bits;
    }

    inline void fill_buff() { fill_buff_inlined(); }

    inline void fill_buff_reverse() {
        static_assert(Kmer::uint_kmer_bits % 64 == 0);
        /* fill the buffer with the bits at [m_pos - uint_kmer_bits, m_pos),
           in order: `append64` appends at the low end, so words must be
           appended from the highest address down to the lowest */
        for (int i = 64; i <= Kmer::uint_kmer_bits; i += 64) {
            m_buff.append64(m_bv->get_word64(std::max<uint64_t>(m_pos, Kmer::uint_kmer_bits) - i));
        }
        m_avail = std::min<uint64_t>(m_pos, Kmer::uint_kmer_bits);
        m_buff.pad(Kmer::uint_kmer_bits - m_avail);
    }

    BitVector const* m_bv;
    uint64_t m_uint_kmer_bits;
    uint64_t m_pos, m_avail;
    Kmer m_buff;
};

}  // namespace sshash
