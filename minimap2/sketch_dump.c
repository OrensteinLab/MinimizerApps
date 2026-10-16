// Test driver (not built by `make all`): runs one sketch over a sequence read
// from stdin and prints each emitted minimizer's raw end position, one per
// line. validate_mod_minimizer.py compares this against a brute-force
// reference, so the check covers the compiled order.cpp, not a port of it.
//
// Usage: sketch_dump <order> <w> <k> <seed> <t> [gm_file [at0-cg1|at1-cg0]] < seq.txt
// (order: native, xor, xor-pcg-hash, greedymini, mod-minimizer-xor; HPC off;
// t = 0 for non-mod orders)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minimap.h"
#include "mmpriv.h"
#include "order.h"

int main(int argc, char *argv[])
{
	mm_order_conf_t conf;
	mm_order_t *ord = 0;
	mm128_v p = {0, 0, 0};
	char errbuf[256], *seq = 0;
	size_t cap = 0, len = 0, i;
	int c, w, k;
	if (argc < 6) {
		fprintf(stderr, "Usage: sketch_dump <order> <w> <k> <seed> <t> < seq\n");
		return 1;
	}
	w = atoi(argv[2]), k = atoi(argv[3]);
	while ((c = getchar()) != EOF) {
		if (c == '\n' || c == '\r') continue;
		if (len + 1 >= cap) cap = cap ? cap << 1 : 1 << 16, seq = (char*)realloc(seq, cap);
		seq[len++] = (char)c;
	}
	if (len == 0) return 0;
	memset(&conf, 0, sizeof(conf));
	if (strcmp(argv[1], "native") != 0) {
		if (mm_order_parse_kind(argv[1], &conf.kind) != 0) {
			fprintf(stderr, "unknown order '%s'\n", argv[1]);
			return 1;
		}
		conf.seed = strtoull(argv[4], 0, 0);
		conf.t = atoi(argv[5]);
		conf.k = k;
		conf.w = w;
		if (argc > 6) conf.greedymini_file = argv[6];
		if (argc > 7 && strcmp(argv[7], "at1-cg0") == 0) conf.gm_projection = MM_GM_PROJ_AT1_CG0;
		if ((ord = mm_order_init(&conf, errbuf, sizeof errbuf)) == 0) {
			fprintf(stderr, "%s\n", errbuf);
			return 1;
		}
		mm_sketch_ordered(ord, 0, seq, (int)len, w, k, 0, 0, &p);
	} else mm_sketch(0, seq, (int)len, w, k, 0, 0, &p);
	for (i = 0; i < p.n; ++i)
		printf("%u\n", (uint32_t)p.a[i].y >> 1);
	free(p.a); free(seq);
	mm_order_destroy(ord);
	return 0;
}
