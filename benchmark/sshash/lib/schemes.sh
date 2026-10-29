# Scheme name -> sshash build flags (see ../../../orders/README.md and
# scripts/schemes.py, the Python-side twin of this table). Every scheme
# (including native) takes a seed with sshash -- unlike kmc's lex, there is
# no seed-less baseline to special-case here.
#
# Schemes carry only their canonical name (benchmark/common/scheme_names.py). A
# GreedyMini scheme runs once per variant, a (table, projection) pair from
# config.yaml's `greedymini` / `mod_minimizer_greedymini` sections (override with
# --gm-tables / --gm-projections); a variant's gm_k, projection, gm_table and
# mod_t travel as separate fields (see expand_scheme_variants), never in the name.

ALL_SCHEMES="native,xor,xor-pcg-hash,greedymini,mod-minimizer-xor,mod-minimizer-greedymini"
# (GM_TABLES / GM_PROJECTIONS, set in lib/constants.sh: comma-separated overrides, empty = config.yaml's)

# Run-directory segment of a variant: "<scheme>[/k<gm_k>_<projection>_t<mod_t>]"
# (benchmark/common/scheme_names.py:scheme_dir_name, called through python so the
# two never drift). Args: scheme [gm_k projection mod_t] ("" or "-" = none).
scheme_dir_name() {
  local scheme="$1" gm_k="${2:--}" projection="${3:--}" mod_t="${4:--}"
  local args=("$scheme")
  [[ "$gm_k" == "-" || -z "$gm_k" ]] || args+=(--gm-k "$gm_k")
  [[ "$projection" == "-" || -z "$projection" ]] || args+=(--projection "$projection")
  [[ "$mod_t" == "-" || -z "$mod_t" ]] || args+=(--mod-t "$mod_t")
  bench_py schemes dir-name "${args[@]}"
}

# Fills the SCHEME_VARIANTS array with one tab-separated line per run variant of the
# comma-separated scheme list `$1` at k=`$2`, m=`$3`: scheme, gm_file, projection,
# gm_k, gm_table, mod_t ("-" = empty). Split a line with unpack_variant.
expand_scheme_variants() {
  local schemes="$1" k="$2" m="$3" line
  local args=("$schemes" --k "$k" --m "$m")
  [[ -z "$GM_TABLES" ]] || args+=(--gm-tables "$GM_TABLES")
  [[ -z "$GM_PROJECTIONS" ]] || args+=(--gm-projections "$GM_PROJECTIONS")
  SCHEME_VARIANTS=()
  while IFS= read -r line; do
    SCHEME_VARIANTS+=("$line")
  done < <(bench_py schemes variants "${args[@]}") || true
  [[ "${#SCHEME_VARIANTS[@]}" -gt 0 ]] || die "no runnable schemes in '$schemes' (see orders/README.md)"
}

# Sets V_SCHEME, V_GM_FILE, V_PROJECTION, V_GM_K, V_GM_TABLE, V_MOD_T from one
# expand_scheme_variants line ("-" -> empty).
unpack_variant() {
  local f
  IFS=$'\t' read -r V_SCHEME V_GM_FILE V_PROJECTION V_GM_K V_GM_TABLE V_MOD_T <<<"$1"
  for f in V_GM_FILE V_PROJECTION V_GM_K V_GM_TABLE V_MOD_T; do
    [[ "${!f}" != "-" ]] || printf -v "$f" '%s' ""
  done
}

# Fills SCHEME_FLAGS with sshash build's -s/--scheme/--greedy-mini-file/--mod-t/--gm-projection/--regular
# flags for scheme `$1`, seed `$2` (GreedyMini table `$3`, modality `$4`:
# "canonical" default or "regular", sshash build --regular; projection `$5`,
# default at0-cg1; mod-minimizer t `$6`, default: derived from k, m).
#
# `$2` comes in as config.yaml's hex literal (e.g. "0xDEADBEEF") -- sshash's
# CLI parser (external/pthash/external/cmd_line_parser) reads -s via
# strtoull(value, nullptr, 10), i.e. BASE 10 ONLY, so a "0x..." string
# silently parses as 0 (strtoull stops at the non-decimal 'x'). $((seed))
# converts it to a decimal string sshash actually understands correctly.
build_scheme_flags() {
  local scheme="$1" seed="$2" gm_file="$3" modality="${4:-canonical}" projection="${5:-}" t="${6:-}"
  SCHEME_FLAGS=(-s "$((seed))")
  case "$scheme" in
    native) ;;
    xor|xor-pcg-hash) SCHEME_FLAGS+=(--scheme "$scheme") ;;
    greedymini|mod-minimizer-greedymini)
      [[ -f "$gm_file" ]] || die "GreedyMini table not found: $gm_file"
      SCHEME_FLAGS+=(--scheme "$scheme" --greedy-mini-file "$gm_file")
      if [[ "$scheme" == mod-minimizer-greedymini ]]; then
        [[ -n "$t" ]] || t="$(bench_py schemes mod-t "$scheme" "${k:-$K}" "${m:-$M}" --gm-file "$gm_file")"
        SCHEME_FLAGS+=(--mod-t "$t")
      fi
      # Only pass the flag for the non-default polarity, so default runs use
      # exactly the flags they always did.
      [[ -z "$projection" || "$projection" == "at0-cg1" ]] || SCHEME_FLAGS+=(--gm-projection "$projection")
      ;;
    mod-minimizer-xor)
      # k, m: run_one's (the caller's) values, else the defaults
      [[ -n "$t" ]] || t="$(bench_py schemes mod-t "$scheme" "${k:-$K}" "${m:-$M}")"
      SCHEME_FLAGS+=(--scheme "$scheme" --mod-t "$t")
      ;;
    *) die "unknown scheme: $scheme (see orders/README.md)" ;;
  esac
  if [[ "$modality" == "regular" ]]; then
    [[ "$scheme" != mod-minimizer-* ]] || die "the mod-minimizer schemes have no regular modality"
    SCHEME_FLAGS+=(--regular)
  fi
}

resolve_schemes() {
  if [[ "$SCHEMES" == "all" ]]; then
    echo "$ALL_SCHEMES"
  else
    echo "$SCHEMES"
  fi
}
