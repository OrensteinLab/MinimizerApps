CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11
COMMON := ../../../common/minimizer-schemes
KMC_API := ../../../kmc/kmc_api

.PHONY: all clean

all: gm_sigma_extend

gm_sigma_extend: gm_sigma_extend.cpp $(COMMON)/greedymini_rank.h $(COMMON)/dna_revcomp.h
	$(CXX) $(CXXFLAGS) -o $@ gm_sigma_extend.cpp

# Signature microbenchmark (experiments/profile_run.sh), built like kmc itself (-O3 -fsigned-char, no -march)
bench_signature_impls: bench_signature_impls.cpp $(KMC_API)/mmer.cpp $(KMC_API)/custom_order.cpp
	$(CXX) -O3 -std=c++17 -fsigned-char -o $@ bench_signature_impls.cpp $(KMC_API)/mmer.cpp $(KMC_API)/custom_order.cpp

clean:
	rm -f gm_sigma_extend bench_signature_impls
