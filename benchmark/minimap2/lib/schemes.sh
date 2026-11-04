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
# directory, own rows). mod-minimizer-xor / mod-minimizer-greedymini are the
# strand-symmetric mod-minimizer with xor, respectively a GreedyMini table of
# t-mers, as O_t (minimap2/scheme.cpp); mod-minimizer-greedymini runs once per
# projection too, and only where the table's gm_k is the dataset's t
# (run_one skips it otherwise: w20_k9 fits k=19, w=10, t=9, not map-ont's k=15).
# The scheme lists come from benchmark/common/scheme_names.py (one Python call at
# source time): every scheme (the baseline native + CANONICAL), the greedymini
# schemes (a GreedyMini table) and the mod-minimizer schemes (--mod-t).
IFS=';' read -r ALL_SCHEMES _GM_SCHEMES _MOD_SCHEMES < <(
  PYTHONPATH="$REPO_ROOT/benchmark/common${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -c '
from scheme_names import CANONICAL, GREEDYMINI_SCHEMES, MOD_SCHEMES
print(";".join(",".join(x) for x in (("native",) + CANONICAL, GREEDYMINI_SCHEMES, MOD_SCHEMES)))')
# An empty list (no usable $PYTHON, scheme_names.py not importable) would make every sweep a silent no-op.
if [ -z "$ALL_SCHEMES" ] || [ -z "$_GM_SCHEMES" ] || [ -z "$_MOD_SCHEMES" ]; then
  echo "error: could not read the scheme lists from benchmark/common/scheme_names.py (PYTHON=$PYTHON)" >&2
  exit 1
fi
DEFAULT_SCHEMES="$ALL_SCHEMES"

# GreedyMini table stems under $GM_DIR (their k, gm_k, is parsed from the name):
# greedymini's table of k-mers (center-extended to k) and mod-minimizer-greedymini's
# table of t-mers (t = gm_k; the w20_k9 table is KMC3's); the projections `sweep`
# covers, and the one a single `run` uses.
GM_TABLE=${GM_TABLE:-w10_k15_canonical}
MOD_GM_TABLE=${MOD_GM_TABLE:-w20_k9_canonical}
PROJECTIONS=${PROJECTIONS:-at0-cg1,at1-cg0}
PROJECTION=${PROJECTION:-at0-cg1}

_in_list() {  # <item> <comma-separated list>
  case ",$2," in
    *",$1,"*) return 0 ;;
    *) return 1 ;;
  esac
}

scheme_is_known() { _in_list "$1" "$ALL_SCHEMES"; }

scheme_needs_gm() { _in_list "$1" "$_GM_SCHEMES"; }

# The GreedyMini table stem of a scheme (only meaningful when scheme_needs_gm): a mod-minimizer
# ranks t-mers (MOD_GM_TABLE), any other scheme k-mers (GM_TABLE).
scheme_gm_table() {
  if scheme_is_mod "$1"; then echo "$MOD_GM_TABLE"; else echo "$GM_TABLE"; fi
}

scheme_is_mod() { _in_list "$1" "$_MOD_SCHEMES"; }

scheme_gm_file() {
  # Only meaningful when scheme_needs_gm.
  echo "$GM_DIR/$(scheme_gm_table "$1").gm"
}

# Stdlib-only helpers from benchmark/common/scheme_names.py (no parsing here).
_scheme_names() {
  PYTHONPATH="$REPO_ROOT/benchmark/common${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -c "$@"
}

# GreedyMini table k of a scheme, from the table name (w10_k15_canonical -> 15).
scheme_gm_k() {
  _scheme_names 'import sys; from scheme_names import parse_table_name; print(parse_table_name(sys.argv[1])["gm_k"])' "$(scheme_gm_table "$1")"
}

# Whether the scheme's GreedyMini table fits what it ranks (scheme_names.gm_table_fits):
# $1 scheme, $2 k, $3 mod_t. Returns 0 if it fits (always, for a scheme without a
# table), 1 if not, 2 on an error (e.g. an unparsable table name), so a caller can
# skip a scheme that does not apply without also skipping a misconfigured one.
scheme_table_fits() {
  scheme_needs_gm "$1" || return 0
  _scheme_names '
import os, sys
from scheme_names import gm_table_fits, parse_table_name
scheme, gm_file, k, t = sys.argv[1:5]
try:
    if not os.path.isfile(gm_file):
        raise FileNotFoundError("no GreedyMini table %s" % gm_file)
    fits = gm_table_fits(scheme, parse_table_name(gm_file)["gm_k"], k, t or None)
except Exception as e:
    sys.exit("scheme_table_fits: %s" % e)  # exit status 1: mapped to 2 below
sys.exit(0 if fits else 3)
' "$1" "$(scheme_gm_file "$1")" "$2" "$3"
  case $? in
    0) return 0 ;;
    3) return 1 ;;
    *) return 2 ;;
  esac
}

# Normalised projection ('at0-cg1' | 'at1-cg0'); fails on anything else.
scheme_projection() {
  _scheme_names 'import sys; from scheme_names import projection_from_flag; print(projection_from_flag(sys.argv[1]))' "${1:-$PROJECTION}"
}

# Run-directory segment(s): <scheme>[/k<gm_k>_<projection>_t<mod_t>] (scheme_names.scheme_dir_name).
# $1 scheme, $2 mod_t. Uses the scheme's table and the current PROJECTION.
scheme_dir_name() {
  local scheme="$1" mod_t="${2:-}" gm_k="" proj=""
  if scheme_needs_gm "$scheme"; then
    gm_k=$(scheme_gm_k "$scheme") || return 1
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
# (scheme, seed, mod_t) combination. `native` needs no flags at all; the
# greedymini schemes need --greedymini-file and --gm-projection ($PROJECTION),
# their --order-seed toggles the Gamma tie-break (as kraken2's and sshash's
# greedymini seeds); a mod-minimizer-* scheme additionally needs --mod-t.
build_scheme_flags() {
  local scheme="$1" seed="$2" mod_t="$3"
  SCHEME_FLAGS=()
  [ "$scheme" = "native" ] && return 0
  SCHEME_FLAGS+=(--scheme "$scheme" --order-seed "$seed")
  if scheme_needs_gm "$scheme"; then
    local proj
    proj=$(scheme_projection) || return 1
    SCHEME_FLAGS+=(--greedymini-file "$(scheme_gm_file "$scheme")" --gm-projection "$proj")
  fi
  if scheme_is_mod "$scheme"; then
    SCHEME_FLAGS+=(--mod-t "$mod_t")
  fi
}

# Every scheme, native included, is swept across the full seed list (minimap2_bench.sh sweep). native
# has no order-dependent behavior, so its seeds build identical indexes and its repeats only measure
# run-to-run timing noise; they give the baseline the same number of runs as the other schemes, which
# the significance tests need. The other schemes take --order-seed (greedymini: the Gamma tie-break toggle).
