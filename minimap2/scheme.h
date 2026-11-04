// Custom minimizer order integration (lab fork, not part of upstream
// minimap2). See common/minimizer-schemes/README.md (orders/README.md at the
// monorepo root) for the shared order implementation this wraps, and this
// file's companion scheme.cpp for the design rationale (rank-vs-index-key
// separation, allocator/threading rules, HPC/palindrome handling).
//
// scheme.cpp is the only C++ translation unit in minimap2/; everything else
// stays plain C. This header is extern-"C" safe (like minimap.h already is)
// so it can be included from both.

#ifndef MM_SCHEME_H
#define MM_SCHEME_H

#include <stdint.h>
#include <stdio.h>
#include "minimap.h" // mm128_v, mm_scheme_kind_t, mm_scheme_conf_t, mm_scheme_t

#ifdef __cplusplus
extern "C" {
#endif

// Materializes conf into a ready-to-use order context: loads the GreedyMini
// table from disk for the greedymini schemes (MM_SCHEME_GREEDYMINI, and
// MM_SCHEME_MOD_GREEDYMINI, whose table ranks the t-mers), validates conf->t
// for MOD_* kinds. Returns NULL and writes a human-readable message into
// errbuf (always NUL-terminated) on failure; conf->greedymini_file only
// needs to stay valid for the duration of this call.
mm_scheme_t *mm_scheme_init(const mm_scheme_conf_t *conf, char *errbuf, size_t errbuf_len);

// minimap2's own native invertible hash (sketch.c), exported (was file-local
// `static inline hash64`) so scheme.cpp can derive the *stored* index key the
// same way for every order -- see scheme.cpp's top-of-file design comment.
uint64_t mm_hash64(uint64_t key, uint64_t mask);

void mm_scheme_destroy(mm_scheme_t *ord);

// A fresh, independently-destroyable copy sharing the same (read-only)
// GreedyMini table when present. NULL-safe (clone of NULL is NULL). Used to
// give each part of a multi-part index its own mm_idx_t::order while a
// single mm_scheme_t is only ever materialized/deserialized once (see
// mm_idx_reader_read in index.c).
mm_scheme_t *mm_scheme_clone(const mm_scheme_t *ord);

mm_scheme_kind_t mm_scheme_kind(const mm_scheme_t *ord);

// One call per sequence -- same contract as mm_sketch(): km/allocator rules
// (kalloc vs libc realloc, matching whichever the call site already uses),
// output appended to *p (caller may set p->n = 0 first). index.c and map.c
// call this instead of mm_sketch() whenever mi->scheme != NULL.
void mm_sketch_scheme(const mm_scheme_t *ord, void *km, const char *str, int len,
                        int w, int k, uint32_t rid, int is_hpc, mm128_v *p);

// .mmi persistence: write the order section ([u64 len][blob], including the
// GreedyMini table itself, so a reload never needs --greedymini-file) at the
// file position, or read one back (no filesystem access besides fp). Both
// stream: no full-size temporary buffer. 0 / non-NULL on success; the loader
// writes a message to errbuf on failure.
int mm_scheme_dump(const mm_scheme_t *ord, FILE *fp);
mm_scheme_t *mm_scheme_load(FILE *fp, char *errbuf, size_t errbuf_len);

// Checks that ord can sketch the index's k-mers with window w (k <= 28; for
// mod-minimizer 0 < t < k and t == k (mod w); for greedymini a table built for
// length k), as mm_idx_load does for an order read from an .mmi. Returns 0,
// or -1 with a message in errbuf.
int mm_scheme_validate(const mm_scheme_t *ord, int k, int w, char *errbuf, size_t errbuf_len);

// 1 if conf (the command-line order) describes ord (an index's stored order):
// same kind, seed, mod t, projection and, when conf names one, the same
// GreedyMini table. Used to warn when a prebuilt index overrides the options.
int mm_scheme_matches_conf(const mm_scheme_t *ord, const mm_scheme_conf_t *conf);

// The greedymini sketch (sketch_gm.c; mm_sketch_scheme calls it for
// MM_SCHEME_GREEDYMINI) and what it needs of the scheme, set when the scheme
// is created or loaded and read-only afterwards: ranks, the GreedyMini table
// (0xFFFF for a key outside the UHS); l, the ranked length (the k-mer length
// k = gm_k + 2i); i, the center offset; toggle, the Gamma tie-break toggle
// (the low l bits of the seed); flip, the Sigma polarity of the projection
// (0 = at0-cg1, 1 = at1-cg0). Same contract as mm_sketch().
typedef struct {
	const uint16_t *ranks;
	int l, i, gm_k;
	uint64_t toggle, flip;
} mm_gm_params_t;
void mm_sketch_gm(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc,
                  const mm_gm_params_t *gp, mm128_v *p);

// CLI helpers (main.c): name <-> kind, both directions.
const char *mm_scheme_kind_name(mm_scheme_kind_t kind);
int mm_scheme_parse_kind(const char *s, mm_scheme_kind_t *out); // 0 on success, -1 on unknown name

#ifdef __cplusplus
}
#endif

#endif // MM_SCHEME_H
