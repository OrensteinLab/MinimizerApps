#!/usr/bin/env bash
# Scheme flag builder for the minimap2 benchmark harness, mirroring
# benchmark/sshash/lib/schemes.sh's shape (see orders/README.md for the
# CLI-flag-mapping table across all four apps).
#
# Schemes are the canonical names (benchmark/common/scheme_names.py) plus the
# baseline `native`. The GreedyMini parameters are NOT part of the scheme name:
# the table (GM_TABLE, its k is gm_k) and the AT/CG projection (PROJECTION,
# --gm-projection at0-cg1|at1-cg0) are separate parameters. `sweep` runs
# greedymini once per projection in PROJECTIONS, each its own run (own
# directory, own rows). mod-minimizer-xor is the strand-symmetric
# mod-minimizer with xor as O_t (minimap2/scheme.cpp).
ALL_SCHEMES="native,xor,xor-pcg-hash,greedymini,mod-minimizer-xor"
DEFAULT_SCHEMES="$ALL_SCHEMES"

# GreedyMini table stem under $GM_DIR (its k, gm_k, is parsed from the name);
# the projections `sweep` covers, and the one a single `run` uses.
GM_TABLE=${GM_TABLE:-w10_k15_canonical}
PROJECTIONS=${PROJECTIONS:-at0-cg1,at1-cg0}
PROJECTION=${PROJECTION:-at0-cg1}

scheme_is_known() {
  case ",$ALL_SCHEMES," in
    *",$1,"*) return 0 ;;
    *) return 1 ;;
  esac
}

scheme_needs_gm() {
  [ "$1" = greedymini ]
}

scheme_is_mod() {
  case "$1" in
    mod-minimizer-*) return 0 ;;
    *) return 1 ;;
  esac
}

scheme_gm_file() {
  # Only meaningful when scheme_needs_gm.
  echo "$GM_DIR/$GM_TABLE.gm"
}

# Stdlib-only helpers from benchmark/common/scheme_names.py (no parsing here).
_scheme_names() {
  PYTHONPATH="$REPO_ROOT/benchmark/common${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -c "$@"
}

# GreedyMini table k, from the table name (w10_k15_canonical -> 15).
scheme_gm_k() {
  _scheme_names 'import sys; from scheme_names import parse_table_name; print(parse_table_name(sys.argv[1])["gm_k"])' "$GM_TABLE"
}

# Normalised projection ('at0-cg1' | 'at1-cg0'); fails on anything else.
scheme_projection() {
  _scheme_names 'import sys; from scheme_names import projection_from_flag; print(projection_from_flag(sys.argv[1]))' "${1:-$PROJECTION}"
}

# Run-directory segment(s): <scheme>[/k<gm_k>_<projection>_t<mod_t>] (scheme_names.scheme_dir_name).
# $1 scheme, $2 mod_t. Uses the current GM_TABLE / PROJECTION.
scheme_dir_name() {
  local scheme="$1" mod_t="${2:-}" gm_k="" proj=""
  if scheme_needs_gm "$scheme"; then
    gm_k=$(scheme_gm_k) || return 1
    proj=$(scheme_projection) || return 1
  fi
  _scheme_names '
import sys
from scheme_names import scheme_dir_name
s, k, p, t = sys.argv[1:5]
print(scheme_dir_name(s, int(k) if k else None, p or None, int(t) if t else None))
' "$scheme" "$gm_k" "$proj" "$mod_t"
}

# Fills the global array SCHEME_FLAGS with the minimap2 CLI args for one
# (scheme, seed, mod_t) combination. `native` needs no flags at all; greedymini
# needs --greedymini-file and --gm-projection ($PROJECTION), its --order-seed
# toggles the Gamma tie-break (as kraken2's and sshash's greedymini seeds); a
# mod-minimizer-* scheme additionally needs --mod-t.
build_scheme_flags() {
  local scheme="$1" seed="$2" mod_t="$3"
  SCHEME_FLAGS=()
  [ "$scheme" = "native" ] && return 0
  if scheme_needs_gm "$scheme"; then
    local proj
    proj=$(scheme_projection) || return 1
    SCHEME_FLAGS+=(--scheme greedymini --greedymini-file "$(scheme_gm_file "$scheme")" --order-seed "$seed" --gm-projection "$proj")
    return 0
  fi
  SCHEME_FLAGS+=(--scheme "$scheme" --order-seed "$seed")
  if scheme_is_mod "$scheme"; then
    SCHEME_FLAGS+=(--mod-t "$mod_t")
  fi
}

# Every scheme, native included, is swept across the full seed list (minimap2_bench.sh sweep). native
# has no order-dependent behavior, so its seeds build identical indexes and its repeats only measure
# run-to-run timing noise; they give the baseline the same number of runs as the other schemes, which
# the significance tests need. The other schemes take --order-seed (greedymini: the Gamma tie-break toggle).
