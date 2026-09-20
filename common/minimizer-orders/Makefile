CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11

.PHONY: all clean

all: gm_serializer

gm_serializer: gm_serializer.cpp greedymini_rank.h
	$(CXX) $(CXXFLAGS) -o $@ gm_serializer.cpp

clean:
	rm -f gm_serializer
