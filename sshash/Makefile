.PHONY: all clean test

BUILD_DIR = build

all:
	mkdir -p $(BUILD_DIR) && cd $(BUILD_DIR) && cmake .. -DCMAKE_BUILD_TYPE=Release && $(MAKE) -j

clean:
	rm -rf $(BUILD_DIR)

test: all
	cd $(BUILD_DIR) && ./test_alphabet 31 && ./test_minimizer
	# GreedyMini / mod-minimizer tests (need rank tables):
	#   ./test_gm_center DIR  ./test_gm_rolling DIR  ./test_streaming_scheme INPUT.fa TABLE.gm  ./test_mod_minimizer
