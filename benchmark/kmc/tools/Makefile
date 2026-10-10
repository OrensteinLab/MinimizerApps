CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11
COMMON := ../../../common/minimizer-orders

.PHONY: all clean

all: gm_sigma_extend

gm_sigma_extend: gm_sigma_extend.cpp $(COMMON)/greedymini_rank.h $(COMMON)/dna_revcomp.h
	$(CXX) $(CXXFLAGS) -o $@ gm_sigma_extend.cpp

clean:
	rm -f gm_sigma_extend
