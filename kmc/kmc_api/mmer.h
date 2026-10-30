/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/
// Modified by the Orenstein Lab, 2026: ranking signatures by a custom order. See kmc/MODIFICATIONS.md.

#ifndef _MMER_H
#define _MMER_H
#include <cinttypes>
#include "../../common/minimizer-schemes/greedymini_rank.h"
#include "../../common/minimizer-schemes/dna_revcomp.h"
#ifndef MIN
#define MIN(x,y)	((x) < (y) ? (x) : (y))
#endif
#ifndef uchar
using uchar = unsigned char;
#endif

// *************************************************************************
// *************************************************************************


class CMmer
{
	uint32_t str;
	uint32_t str_rc;
	uint32_t mask;
	uint64_t current_val;
	mutable uint32_t sig_val;
	mutable bool sig_valid;
	uint32_t bin_key;
	uint32_t bin_key_mask;
	uint32_t rc_shift;
	uint32_t* norm;
	uint32_t len;

	// GM DNA LUT: 4^p sort index in gm_dna_sort_lut; canonical DNA in norm[].
	static const uint32_t* gm_dna_sort_lut;
	// GM binary LUT: 2^p uint32 ranks over the projected DNA (-oat1 selects at1-cg0; default at0-cg1).
	static bool gm_binary_lut;
	static const uint32_t* gm_binary_rank_lut;
	static uint32_t gm_tie_seed;
	static uint32_t norm5[1 << 10];
	static uint32_t norm6[1 << 12];
	static uint32_t norm7[1 << 14];
	static uint32_t norm8[1 << 16];
	static uint32_t norm9[1 << 18];
	static uint32_t norm10[1 << 20];
	static uint32_t norm11[1 << 22];

	static bool is_allowed(uint32_t mmer, uint32_t len)
	{
		if ((mmer & 0x3f) == 0x3f)            // TTT suffix
			return false;
		if ((mmer & 0x3f) == 0x3b)            // TGT suffix
			return false;
		if ((mmer & 0x3c) == 0x3c)            // TG* suffix !!!! consider issue #152
			return false;

		for (uint32_t j = 0; j < len - 3; ++j)
			if ((mmer & 0xf) == 0)                // AA inside
				return false;
			else
				mmer >>= 2;

		if (mmer == 0)            // AAA prefix
			return false;
		if (mmer == 0x04)        // ACA prefix
			return false;
		if ((mmer & 0xf) == 0)    // *AA prefix
			return false;
	
		return true;
	}

	static inline uint32_t DnaBaseToBinaryBit(uint32_t base2)
	{
		return (base2 ^ (base2 >> 1)) & 1u;
	}

	friend class CSignatureMapper;
	struct _si
	{			
		static void init_norm(uint32_t* norm, uint32_t len)
		{
			uint32_t special = 1 << len * 2;
			for(uint32_t i = 0 ; i < special ; ++i)
			{				
				uint32_t rev = DnaReverseComplement(i, len);
				uint32_t str_val = is_allowed(i, len) ? i : special;
				uint32_t rev_val = is_allowed(rev, len) ? rev : special;
				norm[i] = MIN(str_val, rev_val);				
			}
		}

		_si()
		{
			init_norm(norm5, 5);
			init_norm(norm6, 6);
			init_norm(norm7, 7);
			init_norm(norm8, 8);
			init_norm(norm9, 9);
			init_norm(norm10, 10);
			init_norm(norm11, 11);
		}

	}static _init;

	inline void pack_binary()
	{
		const uint32_t tie_dna = str ^ gm_tie_seed;
		const uint32_t tie_rc = str_rc ^ gm_tie_seed;
		current_val = ((uint64_t)gm_binary_rank_lut[bin_key] << 32)
			| (uint64_t)(tie_dna < tie_rc ? tie_dna : tie_rc);
	}

public:
	enum { LUT_DEFAULT = 0, LUT_DNA = 1, LUT_BINARY = 2 };
	static int LutMode()
	{
		if (gm_binary_lut) return LUT_BINARY;
		if (gm_dna_sort_lut) return LUT_DNA;
		return LUT_DEFAULT;
	}

	CMmer(uint32_t _len);
	static void SetGmDnaLut(const uint32_t* sort_lut, const uint32_t* canonical_lut, uint32_t p);
	static void SetGmBinaryLut(const uint32_t* rank_lut, uint32_t tie_seed);
	static void ClearCustomOrder();
	// Shared with kraken2's own binary-alphabet GreedyMini key derivation
	// (common/minimizer-schemes/greedymini_rank.h); verified bit-identical to
	// this function's previous from-scratch loop over DnaBaseToBinaryBit.
	static inline uint32_t DnaToBinaryKey(uint32_t dna, uint32_t p)
	{
		return gm_extract_dna_key((uint64_t)dna, (int)p, (int)p);
	}
	// Shared with kraken2's own reverse complement (dna_revcomp.h); verified
	// bit-identical to this function's previous from-scratch implementation.
	static inline uint32_t DnaReverseComplement(uint32_t v, uint32_t p)
	{
		return (uint32_t)dna_reverse_complement((uint64_t)v, (uint8_t)p, DNA_REVCOM_VERSION);
	}
	inline void insert(uchar symb);
	template<int Mode>
	inline void insert_sym(uchar symb);
	inline uint32_t get() const;
	inline bool operator==(const CMmer& x);
	inline bool operator<(const CMmer& x);
	inline void clear();
	inline bool operator<=(const CMmer& x);
	inline void set(const CMmer& x);
	inline void insert(const char* seq);
	
};

//--------------------------------------------------------------------------
template<int Mode>
inline void CMmer::insert_sym(uchar symb)
{
	str <<= 2;
	str += symb;
	str &= mask;

	if (Mode == LUT_BINARY)
	{
		str_rc = ((str_rc >> 2) | ((symb ^ 3u) << rc_shift)) & mask;
		bin_key = ((bin_key << 1) | DnaBaseToBinaryBit(symb)) & bin_key_mask;
		pack_binary();
	}
	else if (Mode == LUT_DNA)
		current_val = gm_dna_sort_lut[str];
	else
		current_val = norm[str];
	sig_valid = false;
}

inline void CMmer::insert(uchar symb)
{
	if (gm_binary_lut)
		insert_sym<LUT_BINARY>(symb);
	else if (gm_dna_sort_lut)
		insert_sym<LUT_DNA>(symb);
	else
		insert_sym<LUT_DEFAULT>(symb);
}

//--------------------------------------------------------------------------
inline uint32_t CMmer::get() const
{
	if (!sig_valid)
	{
		if (gm_binary_lut)
			sig_val = (uint32_t)current_val ^ gm_tie_seed;
		else if (gm_dna_sort_lut)
			sig_val = norm[str];
		else
			sig_val = (uint32_t)current_val;
		sig_valid = true;
	}
	return sig_val;
}

//--------------------------------------------------------------------------
inline bool CMmer::operator==(const CMmer& x)
{
	return current_val == x.current_val;
}

//--------------------------------------------------------------------------
inline bool CMmer::operator<(const CMmer& x)
{
	return current_val < x.current_val;
}

//--------------------------------------------------------------------------
inline void CMmer::clear()
{
	str = 0;
	str_rc = 0;
	bin_key = 0;
	current_val = 0;
	sig_val = 0;
	sig_valid = false;
}

//--------------------------------------------------------------------------
inline bool CMmer::operator<=(const CMmer& x)
{
	return current_val <= x.current_val;
}

//--------------------------------------------------------------------------
inline void CMmer::set(const CMmer& x)
{
	str = x.str;
	str_rc = x.str_rc;
	current_val = x.current_val;
	bin_key = x.bin_key;
	sig_val = x.sig_val;
	sig_valid = x.sig_valid;
}

//--------------------------------------------------------------------------
inline void CMmer::insert(const char* seq)
{
	switch (len)
	{
	case 5: 
		str = (seq[0] << 8) + (seq[1] << 6) + (seq[2] << 4) + (seq[3] << 2) + (seq[4]);
		break;
	case 6:
		str = (seq[0] << 10) + (seq[1] << 8) + (seq[2] << 6) + (seq[3] << 4) + (seq[4] << 2) + (seq[5]);
		break;
	case 7:
		str = (seq[0] << 12) + (seq[1] << 10) + (seq[2] << 8) + (seq[3] << 6) + (seq[4] << 4 ) + (seq[5] << 2) + (seq[6]);
		break;
	case 8:
		str = (seq[0] << 14) + (seq[1] << 12) + (seq[2] << 10) + (seq[3] << 8) + (seq[4] << 6) + (seq[5] << 4) + (seq[6] << 2) + (seq[7]);
		break;
	case 9:
		str = (seq[0] << 16) + (seq[1] << 14) + (seq[2] << 12) + (seq[3] << 10) + (seq[4] << 8) + (seq[5] << 6) + (seq[6] << 4) + (seq[7] << 2) + (seq[8]);
		break;
	case 10:
		str = (seq[0] << 18) + (seq[1] << 16) + (seq[2] << 14) + (seq[3] << 12) + (seq[4] << 10) + (seq[5] << 8) + (seq[6] << 6) + (seq[7] << 4) + (seq[8] << 2) + (seq[9]);
		break;
	case 11:
		str = (seq[0] << 20) + (seq[1] << 18) + (seq[2] << 16) + (seq[3] << 14) + (seq[4] << 12) + (seq[5] << 10) + (seq[6] << 8) + (seq[7] << 6) + (seq[8] << 4) + (seq[9] << 2) + (seq[10]);
		break;
	default:
		break;
	}

	if (gm_binary_lut)
	{
		bin_key = DnaToBinaryKey(str, len);
		str_rc = DnaReverseComplement(str, len);
		pack_binary();
	}
	else if (gm_dna_sort_lut)
		current_val = gm_dna_sort_lut[str];
	else
		current_val = norm[str];
	sig_valid = false;
}


#endif
