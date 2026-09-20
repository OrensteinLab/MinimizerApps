/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Custom signature order: GreedyMini rank table, random hash (xor seed),
  or random PCG (xor seed then RXS-M-XS scramble). Default KMC lex order
  is unchanged until Install().
*/

#ifndef _CUSTOM_ORDER_H
#define _CUSTOM_ORDER_H

#include <string>
#include <cstdint>

class CCustomOrder
{
public:
	enum : uint8_t { kFlagDefault = 0, kFlagCustom = 1 };

	enum class Mode : uint8_t
	{
		DEFAULT = 0,
		TABLE = 1,
		HASH = 2,
		PCG = 3,
	};

	enum class LutBackend : uint8_t
	{
		DnaExpanded = 0,
		BinaryRuntime = 1,
	};

	struct Config
	{
		Mode mode = Mode::DEFAULT;
		LutBackend lut_backend = LutBackend::BinaryRuntime;
		std::string file;
		std::string label;
		uint32_t signature_len = 0;
		uint64_t seed = 0;
	};

	struct Info
	{
		Mode mode = Mode::DEFAULT;
		std::string mode_name;
		std::string label;
		std::string file;
		uint64_t seed = 0;
		uint32_t k = 0;
	};

	static bool Active();
	static Info GetInfo();
	static Mode ParseMode(const std::string& name);

	static void Install(const Config& cfg);
	static void Clear();

	static void WriteSidecar(const std::string& db_prefix);
	static bool InstallFromSidecar(const std::string& db_prefix, uint32_t signature_len);

private:
	static void check_k(uint32_t k, uint32_t signature_len);
	static void activate_gm_dna_lut(uint32_t p);
	static void activate_gm_binary_lut(uint32_t p, uint32_t tie_seed);
	static const char* mode_name(Mode mode);
	static std::string default_label(const std::string& path);
	static void install_table(const std::string& path, const std::string& label, uint32_t signature_len, uint64_t seed, LutBackend lut_backend);
	static void install_random(Mode mode, const std::string& label, uint32_t signature_len, uint64_t seed);
};

#endif
