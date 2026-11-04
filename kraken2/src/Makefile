CXX ?= g++
KRAKEN2_SKIP_FOPENMP ?= -fopenmp
# -march=native enables BMI2 pext used by greedymini_rank.h when available.
CXXFLAGS = $(KRAKEN2_SKIP_FOPENMP) -Wall -std=c++11 -O3 -fPIC -g -march=native -flto
CXXFLAGS += -DLINEAR_PROBING
LDFLAGS ?=
LDFLAGS += -flto
CFLAGS = -Wall -std=c99 -O2 -g -march=native

.PHONY: all clean install test check_prof_mode

PROGS = estimate_capacity build_db classify dump_table lookup_accession_numbers k2mask blast_to_fasta merge libtax.so benchmark_minimizers
TESTS = test_scheme_rankers test_mod_minimizer

# UHS hit counters are always updated in GMRanker (cheap); BENCH only gates
# NextMinimizer RDTSC phase timing in mmscanner.cc.
BENCH_CPPFLAGS = -DKRAKEN2_MINIMIZER_BENCHMARK

# Phase timing and lookup counts in classify (env KRAKEN2_MINIMIZER_PROF path;
# benchmark/kraken2/experiments/profile_run.sh).
# Stamp tracks ENABLE_MINIMIZER_PROF so toggling the flag rebuilds classify.
PROF_STAMP = .minimizer_prof_mode
ifeq ($(ENABLE_MINIMIZER_PROF),1)
  PROF_CPPFLAGS = -DKRAKEN2_MINIMIZER_PROF
  PROF_MODE = 1
else
  PROF_CPPFLAGS =
  PROF_MODE = 0
endif

# Headers whose layout/inlines are merged under -flto; keep deps tight to avoid
# -Wlto-type-mismatch from stale TUs (e.g. classify.o vs mmscanner.o).
MM_HEADERS = mmscanner.h minimizer_rankers.h sketching_schemes.h \
	mod_minimizer.h dna_kmer.h greedymini_rank.h \
	../../common/minimizer-schemes/greedymini_rank.h \
	../../common/minimizer-schemes/greedymini_center.h \
	../../common/minimizer-schemes/mod_minimizer_rule.h

all: check_prof_mode $(PROGS)

# Always ensure stamp exists before objects that list it as a prerequisite
# (direct `make classify` after clean must not fail).
$(PROF_STAMP):
	@echo "$(PROF_MODE)" > "$@"

check_prof_mode: $(PROF_STAMP)
	@desired="$(PROF_MODE)"; \
	if [ "$$(cat "$(PROF_STAMP)" 2>/dev/null)" != "$$desired" ]; then \
	  echo "$$desired" > "$(PROF_STAMP)"; \
	  rm -f classify.o mmscanner.o classify; \
	fi

install: $(PROGS)
	cp $(PROGS) "$(KRAKEN2_DIR)/"

test: $(TESTS)
	./test_scheme_rankers
	./test_mod_minimizer

clean:
	rm -f *.o $(PROGS) $(TESTS) scan_dump $(PROF_STAMP)

taxonomy.o: taxonomy.cc taxonomy.h mmap_file.h
hyperloglogplus.o: hyperloglogplus.cc hyperloglogplus.h
mmap_file.o: mmap_file.cc mmap_file.h
compact_hash.o: compact_hash.cc compact_hash.h kv_store.h mmap_file.h kraken2_data.h
mmscanner.o: mmscanner.cc $(MM_HEADERS) $(PROF_STAMP)
	$(CXX) $(CXXFLAGS) $(PROF_CPPFLAGS) -c -o $@ mmscanner.cc
seqreader.o: seqreader.cc seqreader.h kseq.h
omp_hack.o: omp_hack.cc omp_hack.h
reports.o: reports.cc reports.h kraken2_data.h
aa_translate.o: aa_translate.cc aa_translate.h
utilities.o: utilities.cc utilities.h

classify.o: classify.cc kraken2_data.h kv_store.h taxonomy.h seqreader.h \
	$(MM_HEADERS) compact_hash.h aa_translate.h reports.h utilities.h \
	readcounts.h $(PROF_STAMP)
	$(CXX) $(CXXFLAGS) $(PROF_CPPFLAGS) -c -o $@ classify.cc
dump_table.o: dump_table.cc compact_hash.h taxonomy.h $(MM_HEADERS) kraken2_data.h reports.h
estimate_capacity.o: estimate_capacity.cc kv_store.h $(MM_HEADERS) seqreader.h utilities.h
build_db.o: build_db.cc build_db.h taxonomy.h $(MM_HEADERS) seqreader.h compact_hash.h kv_store.h kraken2_data.h utilities.h
lookup_accession_numbers.o: lookup_accession_numbers.cc mmap_file.h utilities.h
k2mask.o: k2mask.cc gzstream.h threadpool.h
merge.o: merge.cc taxonomy.cc reports.cc mmap_file.cc utilities.cc hyperloglogplus.cc

blast_defline.o: blast_defline.c
blast_to_fasta.o: blast_to_fasta.c
blast_utils.o: blast_utils.c

build_db: build_db.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ build_db.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o utilities.o $(LDFLAGS)

classify: classify.o reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ classify.o reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o $(LDFLAGS)

estimate_capacity: estimate_capacity.o seqreader.o mmscanner.o omp_hack.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ estimate_capacity.o seqreader.o mmscanner.o omp_hack.o utilities.o $(LDFLAGS)

dump_table: dump_table.o mmap_file.o compact_hash.o omp_hack.o taxonomy.o reports.o hyperloglogplus.o
	$(CXX) $(CXXFLAGS) -o $@ dump_table.o mmap_file.o compact_hash.o omp_hack.o taxonomy.o reports.o hyperloglogplus.o $(LDFLAGS)

lookup_accession_numbers: lookup_accession_numbers.o mmap_file.o omp_hack.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ lookup_accession_numbers.o mmap_file.o omp_hack.o utilities.o $(LDFLAGS)

k2mask: k2mask.o seqreader.o
	$(CXX) $(CXXFLAGS) -o $@ k2mask.o seqreader.o $(LDFLAGS) -lz -pthread

merge: merge.o taxonomy.o reports.o mmap_file.o utilities.o hyperloglogplus.o
	$(CXX) $(CXXFLAGS) -o $@ merge.o taxonomy.o reports.o mmap_file.o utilities.o hyperloglogplus.o -lc

blast_to_fasta: blast_to_fasta.o blast_defline.o blast_utils.o
	$(CC) $(CFLAGS) -o  $@  blast_to_fasta.o blast_utils.o blast_defline.o $(LDFLAGS) -lm -lc

libtax.o: libtax.cc

libtax.so: taxonomy.o mmap_file.o libtax.o
	$(CXX) $(CXXFLAGS) -shared -o libtax.so taxonomy.o mmap_file.o libtax.o $(LDFLAGS)

benchmark_minimizers.o: benchmark_minimizers.cc $(MM_HEADERS)
	$(CXX) $(CXXFLAGS) $(BENCH_CPPFLAGS) -c -o $@ benchmark_minimizers.cc

# Separate mmscanner TU: RDTSC phase timing only in the bench binary.
mmscanner_bench.o: mmscanner.cc $(MM_HEADERS)
	$(CXX) $(CXXFLAGS) $(BENCH_CPPFLAGS) -c -o $@ mmscanner.cc

benchmark_minimizers: benchmark_minimizers.o mmscanner_bench.o
	$(CXX) $(CXXFLAGS) -o $@ benchmark_minimizers.o mmscanner_bench.o $(LDFLAGS)

test_scheme_rankers.o: test_scheme_rankers.cc greedymini_rank.h minimizer_rankers.h \
	kraken2_data.h $(MM_HEADERS)
test_scheme_rankers: test_scheme_rankers.o
	$(CXX) $(CXXFLAGS) -o $@ test_scheme_rankers.o $(LDFLAGS)

test_mod_minimizer.o: test_mod_minimizer.cc $(MM_HEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ test_mod_minimizer.cc

test_mod_minimizer: test_mod_minimizer.o mmscanner.o
	$(CXX) $(CXXFLAGS) -o $@ test_mod_minimizer.o mmscanner.o $(LDFLAGS)

# Test/benchmark tool (not part of `all`): key-stream dump and scan timing via
# the public MinimizerScanner API only. Build: make -C kraken2/src scan_dump
scan_dump.o: scan_dump.cc $(MM_HEADERS) utilities.h
	$(CXX) $(CXXFLAGS) -c -o $@ scan_dump.cc

scan_dump: scan_dump.o mmscanner.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ scan_dump.o mmscanner.o utilities.o $(LDFLAGS)
