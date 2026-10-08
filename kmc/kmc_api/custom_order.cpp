#include "custom_order.h"
#include "mmer.h"
#include "../../common/minimizer-orders/pcg_rxs_m_xs.h"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace {

constexpr uint64_t KMC_DNA_GM_MAGIC = 0x314D444E414D4B43ull; // "KMCANDM1"

std::vector<uint64_t> g_gm_binary_ranks;
std::vector<uint32_t> g_gm_binary_ranks32;
std::vector<uint32_t> g_gm_dna_sort;
std::vector<uint32_t> g_gm_dna_canonical;
std::vector<uint64_t> g_gm_rank_by_dna;
std::string g_file;
std::string g_label;
uint64_t g_seed = 0;
CCustomOrder::Mode g_mode = CCustomOrder::Mode::DEFAULT;
uint32_t g_k = 0;
bool g_active = false;
bool g_gm_dna_prebuilt = false;
bool g_gm_binary_lut = false;

uint32_t log2_pow2(uint64_t n)
{
	uint32_t k = 0;
	while (n > 1)
	{
		if (n & 1)
			return 0;
		n >>= 1;
		++k;
	}
	return k;
}

uint64_t read_u64_le(std::ifstream& in)
{
	unsigned char b[8];
	in.read(reinterpret_cast<char*>(b), 8);
	uint64_t x = 0;
	for (int i = 7; i >= 0; --i)
		x = (x << 8) | b[i];
	return x;
}

void write_u64_le(std::ofstream& out, uint64_t x)
{
	unsigned char b[8];
	for (int i = 0; i < 8; ++i)
	{
		b[i] = (unsigned char)(x & 0xFF);
		x >>= 8;
	}
	out.write(reinterpret_cast<char*>(b), 8);
}

uint32_t read_u32_le(std::ifstream& in)
{
	unsigned char b[4];
	in.read(reinterpret_cast<char*>(b), 4);
	return (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

void write_u32_le(std::ofstream& out, uint32_t x)
{
	unsigned char b[4];
	b[0] = (unsigned char)(x & 0xFF);
	b[1] = (unsigned char)((x >> 8) & 0xFF);
	b[2] = (unsigned char)((x >> 16) & 0xFF);
	b[3] = (unsigned char)((x >> 24) & 0xFF);
	out.write(reinterpret_cast<char*>(b), 4);
}

void gm_pick_canonical_dna(uint32_t i, uint32_t p, uint32_t seed, const uint64_t* rank_by_dna, uint32_t& canonical)
{
	const uint32_t rc = CMmer::DnaReverseComplement(i, p);
	const uint64_t ri = rank_by_dna[i];
	const uint64_t rr = rank_by_dna[rc];
	if (ri != rr)
	{
		canonical = ri < rr ? i : rc;
		return;
	}
	const uint32_t tie_dna = i ^ seed;
	const uint32_t tie_rc = rc ^ seed;
	if (tie_dna != tie_rc)
		canonical = tie_dna < tie_rc ? i : rc;
	else
		canonical = i < rc ? i : rc;
}

void build_gm_dna_luts(uint32_t p, uint32_t seed, const uint64_t* rank_by_dna)
{
	const uint32_t n = 1u << (2 * p);
	seed &= 0xFFFFFFFFu;
	g_gm_dna_sort.resize(n);
	g_gm_dna_canonical.resize(n);

	for (uint32_t i = 0; i < n; ++i)
		gm_pick_canonical_dna(i, p, seed, rank_by_dna, g_gm_dna_canonical[i]);

	std::vector<uint32_t> order(n);
	std::iota(order.begin(), order.end(), 0u);
	std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
		const uint64_t ra = rank_by_dna[a];
		const uint64_t rb = rank_by_dna[b];
		if (ra != rb)
			return ra < rb;
		return (a ^ seed) < (b ^ seed);
	});

	std::vector<uint32_t> packed(n);
	for (uint32_t pos = 0; pos < n; ++pos)
		packed[order[pos]] = pos;
	for (uint32_t i = 0; i < n; ++i)
		g_gm_dna_sort[i] = packed[g_gm_dna_canonical[i]];
}

void build_gm_dna_luts_from_binary(const std::vector<uint64_t>& ranks, uint32_t p, uint32_t seed)
{
	const uint32_t n = 1u << (2 * p);
	g_gm_rank_by_dna.resize(n);
	for (uint32_t i = 0; i < n; ++i)
		g_gm_rank_by_dna[i] = ranks[CMmer::DnaToBinaryKey(i, p)];
	build_gm_dna_luts(p, seed, g_gm_rank_by_dna.data());
}

void build_gm_rank_by_dna(CCustomOrder::Mode mode, uint32_t p)
{
	const uint32_t n = 1u << (2 * p);
	g_gm_rank_by_dna.resize(n);
	for (uint32_t i = 0; i < n; ++i)
		g_gm_rank_by_dna[i] = mode == CCustomOrder::Mode::HASH ? (uint64_t)i ^ g_seed : pcg_rxs_m_xs64((uint64_t)i ^ g_seed);
}

bool is_dna_gm_file(const std::string& path)
{
	std::ifstream in(path, std::ios::binary);
	if (!in)
		return false;
	return in && read_u64_le(in) == KMC_DNA_GM_MAGIC;
}

void load_dna_gm_file(const std::string& path, uint64_t& seed, uint32_t& k)
{
	std::ifstream in(path, std::ios::binary);
	if (!in)
		throw std::runtime_error("Cannot open custom order file: " + path);
	if (read_u64_le(in) != KMC_DNA_GM_MAGIC)
		throw std::runtime_error("Not a KMC DNA .gm file: " + path);
	const uint64_t k64 = read_u64_le(in);
	seed = read_u64_le(in);
	const uint64_t n = read_u64_le(in);
	if (k64 < 5 || k64 > 11 || n != (1ull << (2 * k64)))
		throw std::runtime_error("Invalid DNA custom order header in " + path);
	k = (uint32_t)k64;
	g_gm_dna_sort.resize((size_t)n);
	g_gm_dna_canonical.resize((size_t)n);
	for (uint64_t i = 0; i < n; ++i)
	{
		g_gm_dna_sort[(size_t)i] = read_u32_le(in);
		if (!in)
			throw std::runtime_error("Truncated DNA sort table in " + path);
	}
	for (uint64_t i = 0; i < n; ++i)
	{
		g_gm_dna_canonical[(size_t)i] = read_u32_le(in);
		if (!in)
			throw std::runtime_error("Truncated DNA canonical table in " + path);
	}
}

// Classic .gm format ([u64 len][len x u64]); shared with kraken2's own
// GreedyMini table loader (common/minimizer-orders/greedymini_rank.h).
std::vector<uint64_t> load_gm_file(const std::string& path)
{
	auto ranks = load_gm_file_raw(path);
	if (ranks.empty())
		throw std::runtime_error("Cannot open or read custom order file: " + path);
	if (ranks.size() > (1ull << 16))
		throw std::runtime_error("Invalid custom order table size in " + path);
	return ranks;
}

}  // namespace

const char* CCustomOrder::mode_name(Mode mode)
{
	switch (mode)
	{
	case Mode::TABLE:
		if (g_gm_binary_lut)
			return "table_binary";
		return g_gm_dna_prebuilt ? "table_dna" : "table";
	case Mode::HASH: return "hash";
	case Mode::PCG: return "pcg";
	default: return "default";
	}
}

std::string CCustomOrder::default_label(const std::string& path)
{
	std::string base = path;
	auto slash = base.find_last_of("/\\");
	if (slash != std::string::npos)
		base = base.substr(slash + 1);
	auto dot = base.rfind('.');
	if (dot != std::string::npos)
		base = base.substr(0, dot);
	return base;
}

bool CCustomOrder::Active() { return g_active; }

CCustomOrder::Info CCustomOrder::GetInfo()
{
	Info info;
	info.mode = g_mode;
	info.mode_name = mode_name(g_mode);
	info.label = g_label;
	info.file = g_file;
	info.seed = g_seed;
	info.k = g_k;
	return info;
}

CCustomOrder::Mode CCustomOrder::ParseMode(const std::string& name)
{
	if (name.empty())
		return Mode::TABLE;
	if (name == "hash" || name == "random_hash" || name == "xor")
		return Mode::HASH;
	if (name == "pcg" || name == "random_pcg")
		return Mode::PCG;
	if (name == "table" || name == "table_dna" || name == "table_binary")
		return Mode::TABLE;
	if (name == "default")
		return Mode::DEFAULT;
	throw std::runtime_error("Unknown custom order mode: " + name);
}

void CCustomOrder::check_k(uint32_t k, uint32_t signature_len)
{
	if (k != signature_len)
	{
		std::ostringstream ostr;
		ostr << "Custom order k=" << k << " does not match signature length -p" << signature_len;
		throw std::runtime_error(ostr.str());
	}
	if (k < 5 || k > 11)
		throw std::runtime_error("Custom order k must be in 5..11 (KMC signature length)");
}

void CCustomOrder::activate_gm_dna_lut(uint32_t p)
{
	CMmer::SetGmDnaLut(g_gm_dna_sort.data(), g_gm_dna_canonical.data(), p);
	g_k = p;
	g_active = true;
}

void CCustomOrder::activate_gm_binary_lut(uint32_t p, uint32_t tie_seed)
{
	g_gm_binary_ranks32.resize(g_gm_binary_ranks.size());
	for (size_t i = 0; i < g_gm_binary_ranks.size(); ++i)
	{
		if (g_gm_binary_ranks[i] > 0xFFFFFFFFull)
			throw std::runtime_error("GM binary rank does not fit in 32 bits");
		g_gm_binary_ranks32[i] = (uint32_t)g_gm_binary_ranks[i];
	}
	CMmer::SetGmBinaryLut(g_gm_binary_ranks32.data(), tie_seed);
	g_k = p;
	g_active = true;
	g_gm_binary_lut = true;
}

void CCustomOrder::install_table(const std::string& path, const std::string& label, uint32_t signature_len, uint64_t seed, LutBackend lut_backend, bool gm_at1)
{
	g_file = path;
	g_label = label;
	g_seed = seed;
	g_mode = Mode::TABLE;
	g_gm_binary_lut = false;
	g_gm_dna_prebuilt = is_dna_gm_file(path);

	if (g_gm_dna_prebuilt)
	{
		// lut_backend only chooses a runtime representation for a *binary*-alphabet
		// source table; a _dna.gm file is already DNA-native, so it always installs
		// via the DNA-expanded path regardless of the requested backend (relevant
		// now that BinaryRuntime is the ambient default rather than an explicit ask).
		if (gm_at1)
			throw std::runtime_error("-oat1 needs a binary-alphabet .gm table, not a DNA table (" + path + ")");
		uint32_t p = 0;
		load_dna_gm_file(path, g_seed, p);
		g_gm_binary_ranks.clear();
		g_gm_binary_ranks32.clear();
		g_gm_rank_by_dna.clear();
		check_k(p, signature_len);
		activate_gm_dna_lut(p);
		return;
	}

	g_gm_dna_sort.clear();
	g_gm_dna_canonical.clear();
	g_gm_binary_ranks = load_gm_file(path);
	const uint32_t p = log2_pow2(g_gm_binary_ranks.size());
	if (p == 0 || ((uint64_t)1 << p) != g_gm_binary_ranks.size())
		throw std::runtime_error("Custom order table size is not a power of 2");
	check_k(p, signature_len);
	if (gm_at1)
	{
		// The at1-cg0 projection complements every key (A,T -> 1, C,G -> 0), so
		// rank the projection x by the table's rank of ~x: permute the table once.
		// Both backends (binary runtime, DNA-expanded) and the sidecar use it.
		const uint64_t key_mask = ((uint64_t)1 << p) - 1;
		std::vector<uint64_t> flipped(g_gm_binary_ranks.size());
		for (uint64_t key = 0; key < flipped.size(); ++key)
			flipped[key] = g_gm_binary_ranks[key ^ key_mask];
		g_gm_binary_ranks.swap(flipped);
	}

	if (lut_backend == LutBackend::BinaryRuntime)
	{
		activate_gm_binary_lut(p, (uint32_t)g_seed);
		return;
	}

	build_gm_dna_luts_from_binary(g_gm_binary_ranks, p, (uint32_t)g_seed);
	activate_gm_dna_lut(p);
}

void CCustomOrder::install_random(Mode mode, const std::string& label, uint32_t signature_len, uint64_t seed)
{
	g_gm_binary_ranks.clear();
	g_gm_binary_ranks32.clear();
	g_file.clear();
	g_label = label;
	g_seed = seed;
	g_mode = mode;
	g_gm_dna_prebuilt = false;
	check_k(signature_len, signature_len);
	build_gm_rank_by_dna(mode, signature_len);
	build_gm_dna_luts(signature_len, (uint32_t)g_seed, g_gm_rank_by_dna.data());
	activate_gm_dna_lut(signature_len);
}

void CCustomOrder::Install(const Config& cfg)
{
	if (cfg.mode == Mode::DEFAULT)
	{
		Clear();
		return;
	}

	switch (cfg.mode)
	{
	case Mode::TABLE:
	{
		if (cfg.file.empty())
			throw std::runtime_error("Custom order TABLE mode requires a file");
		const std::string label = cfg.label.empty() ? default_label(cfg.file) : cfg.label;
		install_table(cfg.file, label, cfg.signature_len, cfg.seed, cfg.lut_backend, cfg.gm_at1);
		break;
	}
	case Mode::HASH:
	case Mode::PCG:
		if (cfg.label.empty())
			throw std::runtime_error("Custom order label required for random modes");
		install_random(cfg.mode, cfg.label, cfg.signature_len, cfg.seed);
		break;
	default:
		Clear();
		break;
	}
}

void CCustomOrder::Clear()
{
	g_gm_binary_ranks.clear();
	g_gm_binary_ranks32.clear();
	g_gm_dna_sort.clear();
	g_gm_dna_canonical.clear();
	g_gm_rank_by_dna.clear();
	g_file.clear();
	g_label.clear();
	g_seed = 0;
	g_mode = Mode::DEFAULT;
	g_k = 0;
	g_active = false;
	g_gm_dna_prebuilt = false;
	g_gm_binary_lut = false;
	CMmer::ClearCustomOrder();
}

void CCustomOrder::WriteSidecar(const std::string& db_prefix)
{
	if (!g_active)
		return;
	if (g_mode == Mode::TABLE && !g_gm_binary_ranks.empty())
	{
		std::ofstream out(db_prefix + ".kmc_co", std::ios::binary);
		if (!out)
			throw std::runtime_error("Cannot write " + db_prefix + ".kmc_co");
		write_u64_le(out, (uint64_t)g_gm_binary_ranks.size());
		for (uint64_t r : g_gm_binary_ranks)
			write_u64_le(out, r);
	}
	else if (g_mode == Mode::TABLE)
	{
		std::ofstream out(db_prefix + ".kmc_co", std::ios::binary);
		if (!out)
			throw std::runtime_error("Cannot write " + db_prefix + ".kmc_co");
		const uint64_t n = g_gm_dna_sort.size();
		write_u64_le(out, KMC_DNA_GM_MAGIC);
		write_u64_le(out, g_k);
		write_u64_le(out, g_seed);
		write_u64_le(out, n);
		for (uint32_t v : g_gm_dna_sort)
			write_u32_le(out, v);
		for (uint32_t v : g_gm_dna_canonical)
			write_u32_le(out, v);
	}
	else
	{
		std::ofstream out(db_prefix + ".kmc_co", std::ios::binary);
		if (!out)
			throw std::runtime_error("Cannot write " + db_prefix + ".kmc_co");
		write_u64_le(out, g_seed);
	}
	{
		std::ofstream out(db_prefix + ".kmc_col");
		if (!out)
			throw std::runtime_error("Cannot write " + db_prefix + ".kmc_col");
		out << g_label << "\n"
			<< mode_name(g_mode) << "\n"
			<< g_seed << "\n";
	}
}

bool CCustomOrder::InstallFromSidecar(const std::string& db_prefix, uint32_t signature_len)
{
	const std::string col_path = db_prefix + ".kmc_col";
	std::ifstream lab(col_path);
	if (!lab)
		return false;

	Config cfg;
	cfg.signature_len = signature_len;
	std::string mode_line;
	std::string seed_line;
	std::getline(lab, cfg.label);
	std::getline(lab, mode_line);
	std::getline(lab, seed_line);
	if (cfg.label.empty())
		cfg.label = "custom";
	if (!seed_line.empty())
		cfg.seed = std::stoull(seed_line);

	try
	{
		cfg.mode = ParseMode(mode_line);
		if (cfg.mode == Mode::TABLE)
			cfg.file = db_prefix + ".kmc_co";
		Install(cfg);
	}
	catch (...)
	{
		Clear();
		return false;
	}
	return Active();
}
