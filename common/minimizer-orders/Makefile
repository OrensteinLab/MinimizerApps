CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11

.PHONY: all clean test

all: gm_serializer

test: test_mod_minimizer_rule
	./test_mod_minimizer_rule

test_mod_minimizer_rule: test_mod_minimizer_rule.cpp mod_minimizer_rule.h
	$(CXX) $(CXXFLAGS) -o $@ test_mod_minimizer_rule.cpp

gm_serializer: gm_serializer.cpp greedymini_rank.h
	$(CXX) $(CXXFLAGS) -o $@ gm_serializer.cpp

clean:
	rm -f gm_serializer test_mod_minimizer_rule
