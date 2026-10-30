# Modified by the Orenstein Lab, 2026: build changes (dynamic link by default, Python autodetect, raduls flags) and the custom signature orders. See MODIFICATIONS.md.
dummy := $(shell git submodule update --init --recursive)

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
UNAME_P := $(shell uname -p)

KMC_MAIN_DIR = kmc_core
KMC_CLI_DIR = kmc_CLI
KMC_API_DIR = kmc_api
KMC_DUMP_DIR = kmc_dump
KMC_TOOLS_DIR = kmc_tools
PY_KMC_API_DIR = py_kmc_api

OUT_BIN_DIR = bin
OUT_INCLUDE_DIR = include

D_OS =
D_ARCH = 

ifeq ($(UNAME_S),Darwin)
	D_OS=MACOS
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	else
		D_ARCH=X64
	endif
else
	D_OS=LINUX
	D_ARCH=X64
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	endif
	ifeq ($(UNAME_M),aarch64)
		D_ARCH=ARM64
	endif
endif

CPU_FLAGS =
STATIC_CFLAGS = 
STATIC_LFLAGS = 
PY_FLAGS =
LINK_STATIC ?= 0

# Python for py_kmc_api (needs python*-devel). Override: make PYTHON=python3.9 py_kmc_api
PYTHON ?= $(shell \
	for p in python3.11 python3.10 python3.9 python3; do \
		$$p -c 'import sysconfig,os; exit(0 if os.path.isfile(os.path.join(sysconfig.get_path("include"),"Python.h")) else 1)' 2>/dev/null && echo $$p && break; \
	done)
PYTHON_INCLUDES := $(shell $(PYTHON)-config --includes 2>/dev/null)
PYTHON_EXT_SUFFIX := $(shell $(PYTHON)-config --extension-suffix 2>/dev/null)

ifneq ($(strip $(PYTHON)),)
ALL_TARGETS := kmc kmc_dump kmc_tools py_kmc_api
else
ALL_TARGETS := kmc kmc_dump kmc_tools
$(info note: py_kmc_api skipped — no Python.h found; install python3-devel or set PYTHON=)
endif

ifeq ($(D_OS),MACOS)
CC = g++-11

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8.4-a
	else
		CPU_FLAGS = -m64
	endif
	STATIC_CFLAGS = -static-libgcc -static-libstdc++ -pthread
	STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread	
	PY_FLAGS = -Wl,-undefined,dynamic_lookup -fPIC 
else
CC = g++

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8-a
		STATIC_CFLAGS = -pthread
		STATIC_LFLAGS = -pthread
	else
		CPU_FLAGS = -m64
		STATIC_CFLAGS = -pthread
		ifeq ($(LINK_STATIC),1)
			STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
		else
			STATIC_LFLAGS = -lpthread
		endif
	endif
	PY_FLAGS = -fPIC
endif


CFLAGS	= -Wall -O3 -fsigned-char $(CPU_FLAGS) $(STATIC_CFLAGS) -std=c++14
CLINK	= -lm $(STATIC_LFLAGS) -O3 -std=c++14
PY_KMC_API_CFLAGS = $(PY_FLAGS) -Wall -shared -std=c++14 -O3

KMC_CLI_OBJS = \
$(KMC_CLI_DIR)/kmc.o

KFF_OBJS = \
$(KMC_MAIN_DIR)/kff_writer.o

KMC_CORE_OBJS = \
$(KMC_MAIN_DIR)/mem_disk_file.o \
$(KMC_MAIN_DIR)/rev_byte.o \
$(KMC_MAIN_DIR)/bkb_writer.o \
$(KMC_MAIN_DIR)/cpu_info.o \
$(KMC_MAIN_DIR)/bkb_reader.o \
$(KMC_MAIN_DIR)/fastq_reader.o \
$(KMC_MAIN_DIR)/timer.o \
$(KMC_MAIN_DIR)/develop.o \
$(KMC_MAIN_DIR)/kb_completer.o \
$(KMC_MAIN_DIR)/kb_storer.o \
$(KMC_MAIN_DIR)/kmer.o \
$(KMC_MAIN_DIR)/splitter.o \
$(KMC_MAIN_DIR)/kb_collector.o \
$(KMC_MAIN_DIR)/kmc_runner.o

ifeq ($(UNAME_S),Darwin)
ifeq ($(D_ARCH),ARM64)
	RADULS_OBJS = \
	$(KMC_MAIN_DIR)/raduls_neon.o
else
	RADULS_OBJS =
endif
else
ifeq ($(D_ARCH),ARM64)
	RADULS_OBJS = \
	$(KMC_MAIN_DIR)/raduls_neon.o
else
	RADULS_OBJS = \
	$(KMC_MAIN_DIR)/raduls_sse2.o \
	$(KMC_MAIN_DIR)/raduls_sse41.o \
	$(KMC_MAIN_DIR)/raduls_avx2.o \
	$(KMC_MAIN_DIR)/raduls_avx.o
endif
endif

LIB_ZLIB=3rd_party/cloudflare/libz.a
LIB_KMC_CORE = $(OUT_BIN_DIR)/libkmc_core.a


KMC_DUMP_OBJS = \
$(KMC_DUMP_DIR)/nc_utils.o \
$(KMC_DUMP_DIR)/kmc_dump.o

KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/custom_order.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

KMC_API_SRC_FILES = $(wildcard $(KMC_API_DIR)/*.cpp)
PY_KMC_API_OBJS = $(patsubst $(KMC_API_DIR)/%.cpp,$(PY_KMC_API_DIR)/%.o,$(KMC_API_SRC_FILES))

KMC_TOOLS_OBJS = \
$(KMC_TOOLS_DIR)/kmer_file_header.o \
$(KMC_TOOLS_DIR)/kmc_tools.o \
$(KMC_TOOLS_DIR)/nc_utils.o \
$(KMC_TOOLS_DIR)/parameters_parser.o \
$(KMC_TOOLS_DIR)/parser.o \
$(KMC_TOOLS_DIR)/tokenizer.o \
$(KMC_TOOLS_DIR)/fastq_filter.o \
$(KMC_TOOLS_DIR)/fastq_reader.o \
$(KMC_TOOLS_DIR)/fastq_writer.o \
$(KMC_TOOLS_DIR)/percent_progress.o \
$(KMC_TOOLS_DIR)/kff_info_reader.o

$(LIB_ZLIB):
	cd 3rd_party/cloudflare; ./configure; make libz.a

$(KMC_CLI_OBJS) $(KMC_CORE_OBJS) $(KMC_DUMP_OBJS) $(KMC_API_OBJS) $(KFF_OBJS) $(KMC_TOOLS_OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -I 3rd_party/cloudflare -c $< -o $@

$(KMC_MAIN_DIR)/raduls_sse2.o: $(KMC_MAIN_DIR)/raduls_sse2.cpp
	$(CC) $(CFLAGS) -msse2 -mno-sse3 -mno-sse4.1 -mno-avx -c $< -o $@
$(KMC_MAIN_DIR)/raduls_sse41.o: $(KMC_MAIN_DIR)/raduls_sse41.cpp
	$(CC) $(CFLAGS) -msse4.1 -mno-avx -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx.o: $(KMC_MAIN_DIR)/raduls_avx.cpp
	$(CC) $(CFLAGS) -mavx -mno-avx2 -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx2.o: $(KMC_MAIN_DIR)/raduls_avx2.cpp
	$(CC) $(CFLAGS) -mavx2 -c $< -o $@

$(KMC_MAIN_DIR)/raduls_neon.o: $(KMC_MAIN_DIR)/raduls_neon.cpp
	$(CC) $(CFLAGS) -c $< -o $@


$(LIB_KMC_CORE): $(KMC_CORE_OBJS) $(RADULS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS)
	-mkdir -p $(OUT_INCLUDE_DIR)
	cp $(KMC_MAIN_DIR)/kmc_runner.h $(OUT_INCLUDE_DIR)/kmc_runner.h
	-mkdir -p $(OUT_BIN_DIR)
	ar rcs $@ $^

kmc: $(KMC_CLI_OBJS) $(LIB_KMC_CORE) $(LIB_ZLIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) $(CLINK) -o $(OUT_BIN_DIR)/$@ $^

kmc_dump: $(KMC_DUMP_OBJS) $(KMC_API_OBJS)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) $(CLINK) -o $(OUT_BIN_DIR)/$@ $^

kmc_tools: $(KMC_TOOLS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS) $(LIB_ZLIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) $(CLINK) -I 3rd_party/cloudflare -o $(OUT_BIN_DIR)/$@ $^

$(PY_KMC_API_DIR)/%.o: $(KMC_API_DIR)/%.cpp
	$(CC) -c -fPIC -Wall -O3 $(CPU_FLAGS) -std=c++14 $^ -o $@

py_kmc_api: $(PY_KMC_API_OBJS)
	@test -n "$(PYTHON)" || (echo "py_kmc_api: set PYTHON= to a Python with development headers (python3-devel)" && exit 1)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) $(PY_KMC_API_CFLAGS) $(PY_KMC_API_DIR)/py_kmc_api.cpp $(PY_KMC_API_OBJS) \
	-I $(KMC_API_DIR) \
	-I $(PY_KMC_API_DIR)/libs/pybind11/include \
	$(PYTHON_INCLUDES) \
	-o $(OUT_BIN_DIR)/$@$(PYTHON_EXT_SUFFIX)

all: $(ALL_TARGETS)

clean:
	-rm -f $(KMC_MAIN_DIR)/*.o
	-rm -f $(KMC_API_DIR)/*.o
	-rm -f $(KMC_DUMP_DIR)/*.o
	-rm -f $(KMC_TOOLS_DIR)/*.o
	-rm -f $(PY_KMC_API_DIR)/*.o
	-rm -f $(PY_KMC_API_DIR)/*.so
	-rm -rf $(OUT_BIN_DIR)
	-rm -rf $(OUT_INCLUDE_DIR)
	cd 3rd_party/cloudflare; make clean;
