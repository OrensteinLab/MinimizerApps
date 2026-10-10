/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 3.2.4
  Date   : 2024-02-09
*/

#include "../kmc_api/mmer.h"


uint32_t CMmer::norm5[];
uint32_t CMmer::norm6[];
uint32_t CMmer::norm7[];
uint32_t CMmer::norm8[];
uint32_t CMmer::norm9[];
uint32_t CMmer::norm10[];
uint32_t CMmer::norm11[];

const uint32_t* CMmer::gm_dna_sort_lut = nullptr;
bool CMmer::gm_binary_lut = false;
const uint32_t* CMmer::gm_binary_rank_lut = nullptr;
uint32_t CMmer::gm_tie_seed = 0;

CMmer::_si CMmer::_init;

void CMmer::SetGmDnaLut(const uint32_t* sort_lut, const uint32_t* canonical_lut, uint32_t p)
{
	gm_dna_sort_lut = sort_lut;
	gm_binary_lut = false;
	gm_binary_rank_lut = nullptr;
	gm_tie_seed = 0;
	uint32_t* tables[] = { nullptr, nullptr, nullptr, nullptr, nullptr, norm5, norm6, norm7, norm8, norm9, norm10, norm11 };
	if (p < 5 || p > 11)
		return;
	uint32_t* norm_tbl = tables[p];
	const uint32_t n = 1u << (p * 2);
	for (uint32_t i = 0; i < n; ++i)
		norm_tbl[i] = canonical_lut[i];
}

void CMmer::SetGmBinaryLut(const uint32_t* rank_lut, uint32_t tie_seed)
{
	gm_dna_sort_lut = nullptr;
	gm_binary_lut = true;
	gm_binary_rank_lut = rank_lut;
	gm_tie_seed = tie_seed & 0xFFFFFFFFu;
}

void CMmer::ClearCustomOrder()
{
	gm_dna_sort_lut = nullptr;
	gm_binary_lut = false;
	gm_binary_rank_lut = nullptr;
	gm_tie_seed = 0;
	_si::init_norm(norm5, 5);
	_si::init_norm(norm6, 6);
	_si::init_norm(norm7, 7);
	_si::init_norm(norm8, 8);
	_si::init_norm(norm9, 9);
	_si::init_norm(norm10, 10);
	_si::init_norm(norm11, 11);
}

//--------------------------------------------------------------------------
CMmer::CMmer(uint32_t _len)
{
	switch (_len)
	{
	case 5:
		norm = norm5;
		break;
	case 6:
		norm = norm6;
		break;
	case 7:
		norm = norm7;
		break;
	case 8:
		norm = norm8;
		break;
	case 9:
		norm = norm9;
		break;
	case 10:
		norm = norm10;
		break;
	case 11:
		norm = norm11;
		break;
	default:
		break;
	}
	len = _len;
	mask = (1 << _len * 2) - 1;
	bin_key_mask = (1u << _len) - 1u;
	rc_shift = 2u * (_len - 1u);
	str = 0;
	str_rc = 0;
	bin_key = 0;
	current_val = 0;
	sig_val = 0;
	sig_valid = false;
}

//--------------------------------------------------------------------------
