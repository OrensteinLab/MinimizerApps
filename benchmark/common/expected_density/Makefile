# One expected-density tool per app (driver.h + <app>_backend), each built
# against the app's own selection code with the app's own compiler flags.
REPO := $(abspath ../../..)
K2 := $(REPO)/kraken2/src
MM2 := $(REPO)/minimap2
KMCAPI := $(REPO)/kmc/kmc_api
SSH := $(REPO)/sshash
# the pthash submodule (git submodule update --init sshash/external/pthash)
SSH_EXT ?= $(SSH)/external
BUILD := build

.PHONY: all clean
all: $(BUILD)/kraken2_expected_density $(BUILD)/minimap2_expected_density $(BUILD)/kmc3_expected_density $(BUILD)/sshash_expected_density

$(BUILD):
	mkdir -p $@

# kraken2: mmscanner.cc as in kraken2/src/Makefile (mmscanner.o, without profiling)
K2_CXXFLAGS = -fopenmp -Wall -std=c++11 -O3 -g -march=native -DLINEAR_PROBING
$(BUILD)/k2_mmscanner.o: $(K2)/mmscanner.cc $(wildcard $(K2)/*.h) | $(BUILD)
	$(CXX) $(K2_CXXFLAGS) -I$(K2) -c -o $@ $<
$(BUILD)/kraken2_expected_density: kraken2_backend.cc driver.h $(BUILD)/k2_mmscanner.o
	$(CXX) $(K2_CXXFLAGS) -pthread -I$(K2) -o $@ kraken2_backend.cc $(BUILD)/k2_mmscanner.o

# minimap2: sketch.c, kalloc.c and scheme.cpp with minimap2/Makefile's flags
MM2_CFLAGS = -g -Wall -O2 -DHAVE_KALLOC
MM2_CXXFLAGS = -g -Wall -O2 -std=c++11 -march=native -DHAVE_KALLOC
$(BUILD)/mm2_%.o: $(MM2)/%.c $(wildcard $(MM2)/*.h) | $(BUILD)
	$(CC) $(MM2_CFLAGS) -I$(MM2) -c -o $@ $<
$(BUILD)/mm2_scheme.o: $(MM2)/scheme.cpp $(wildcard $(MM2)/*.h) | $(BUILD)
	$(CXX) $(MM2_CXXFLAGS) -I$(MM2) -I$(REPO) -c -o $@ $<
MM2_OBJS = $(BUILD)/mm2_sketch.o $(BUILD)/mm2_kalloc.o $(BUILD)/mm2_misc.o $(BUILD)/mm2_scheme.o
$(BUILD)/minimap2_expected_density: minimap2_backend.cpp driver.h $(MM2_OBJS)
	$(CXX) $(MM2_CXXFLAGS) -pthread -I$(MM2) -o $@ minimap2_backend.cpp $(MM2_OBJS) -lz -lm

# kmc3: the signature orders of kmc_api with kmc/Makefile's flags
KMC_CXXFLAGS = -Wall -O3 -fsigned-char -m64 -pthread -std=c++14
$(BUILD)/kmc_%.o: $(KMCAPI)/%.cpp $(wildcard $(KMCAPI)/*.h) | $(BUILD)
	$(CXX) $(KMC_CXXFLAGS) -c -o $@ $<
KMC_OBJS = $(BUILD)/kmc_custom_order.o $(BUILD)/kmc_mmer.o
$(BUILD)/kmc3_expected_density: kmc3_backend.cpp driver.h $(KMC_OBJS)
	$(CXX) $(KMC_CXXFLAGS) -o $@ kmc3_backend.cpp $(KMC_OBJS)

# sshash: header-only iterators and hashers, with sshash/CMakeLists.txt's flags
SSH_INCLUDES = $(addprefix -I$(SSH_EXT)/pthash/,include external/bits/include external/fastmod \
	external/bits/external/essentials/include external/xxHash external/mm_file/include) -I$(SSH) -I$(SSH)/include -I$(SSH_EXT)/..
SSH_CXXFLAGS = -std=c++17 -O3 -mbmi2 -mavx2 -pthread -Wall
$(BUILD)/sshash_expected_density: sshash_backend.cpp driver.h $(wildcard $(SSH)/include/*.hpp) | $(BUILD)
	$(CXX) $(SSH_CXXFLAGS) $(SSH_INCLUDES) -o $@ sshash_backend.cpp

clean:
	rm -rf $(BUILD)
