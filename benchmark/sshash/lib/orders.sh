# Order name -> sshash build flags (see ../../../orders/README.md and
# scripts/orders.py, the Python-side twin of this table). Every order
# (including native) takes a seed with sshash -- unlike kmc's lex, there is
# no seed-less baseline to special-case here.
#
# Orders carry only their canonical name (benchmark/common/order_names.py). A
# GreedyMini order runs once per variant, a (table, projection) pair from
# config.yaml's `greedymini` / `mod_minimizer_greedymini` sections (override with
# --gm-tables / --gm-projections); a variant's gm_k, projection, gm_table and
# mod_t travel as separate fields (see expand_order_variants), never in the name.

ALL_ORDERS="native,xor,xor-pcg-hash,greedymini,mod-minimizer-xor,mod-minimizer-greedymini"
# (GM_TABLES / GM_PROJECTIONS, set in lib/constants.sh: comma-separated overrides, empty = config.yaml's)

# Run-directory segment of a variant: "<order>[/k<gm_k>_<projection>_t<mod_t>]"
# (benchmark/common/order_names.py:order_dir_name, called through python so the
# two never drift). Args: order [gm_k projection mod_t] ("" or "-" = none).
order_dir_name() {
  local order="$1" gm_k="${2:--}" projection="${3:--}" mod_t="${4:--}"
  local args=("$order")
  [[ "$gm_k" == "-" || -z "$gm_k" ]] || args+=(--gm-k "$gm_k")
  [[ "$projection" == "-" || -z "$projection" ]] || args+=(--projection "$projection")
  [[ "$mod_t" == "-" || -z "$mod_t" ]] || args+=(--mod-t "$mod_t")
  bench_py orders dir-name "${args[@]}"
}

# Fills the ORDER_VARIANTS array with one tab-separated line per run variant of the
# comma-separated order list `$1` at k=`$2`, m=`$3`: order, gm_file, projection,
# gm_k, gm_table, mod_t ("-" = empty). Split a line with unpack_variant.
expand_order_variants() {
  local orders="$1" k="$2" m="$3" line
  local args=("$orders" --k "$k" --m "$m")
  [[ -z "$GM_TABLES" ]] || args+=(--gm-tables "$GM_TABLES")
  [[ -z "$GM_PROJECTIONS" ]] || args+=(--gm-projections "$GM_PROJECTIONS")
  ORDER_VARIANTS=()
  while IFS= read -r line; do
    ORDER_VARIANTS+=("$line")
  done < <(bench_py orders variants "${args[@]}") || true
  [[ "${#ORDER_VARIANTS[@]}" -gt 0 ]] || die "no runnable orders in '$orders' (see orders/README.md)"
}

# Sets V_ORDER, V_GM_FILE, V_PROJECTION, V_GM_K, V_GM_TABLE, V_MOD_T from one
# expand_order_variants line ("-" -> empty).
unpack_variant() {
  local f
  IFS=$'\t' read -r V_ORDER V_GM_FILE V_PROJECTION V_GM_K V_GM_TABLE V_MOD_T <<<"$1"
  for f in V_GM_FILE V_PROJECTION V_GM_K V_GM_TABLE V_MOD_T; do
    [[ "${!f}" != "-" ]] || printf -v "$f" '%s' ""
  done
}

# Fills ORDER_FLAGS with sshash build's -s/-O/--greedy-mini-file/--mod-t/--gm-projection/--regular
# flags for order `$1`, seed `$2` (GreedyMini table `$3`, modality `$4`:
# "canonical" default or "regular", sshash build --regular; projection `$5`,
# default at0-cg1; mod-minimizer t `$6`, default: derived from k, m).
#
# `$2` comes in as config.yaml's hex literal (e.g. "0xDEADBEEF") -- sshash's
# CLI parser (external/pthash/external/cmd_line_parser) reads -s via
# strtoull(value, nullptr, 10), i.e. BASE 10 ONLY, so a "0x..." string
# silently parses as 0 (strtoull stops at the non-decimal 'x'). $((seed))
# converts it to a decimal string sshash actually understands correctly.
build_order_flags() {
  local order="$1" seed="$2" gm_file="$3" modality="${4:-canonical}" projection="${5:-}" t="${6:-}"
  ORDER_FLAGS=(-s "$((seed))")
  case "$order" in
    native) ;;
    xor|xor-pcg-hash) ORDER_FLAGS+=(-O "$order") ;;
    greedymini|mod-minimizer-greedymini)
      [[ -f "$gm_file" ]] || die "GreedyMini table not found: $gm_file"
      ORDER_FLAGS+=(-O "$order" --greedy-mini-file "$gm_file")
      if [[ "$order" == mod-minimizer-greedymini ]]; then
        [[ -n "$t" ]] || t="$(bench_py orders mod-t "$order" "${k:-$K}" "${m:-$M}" --gm-file "$gm_file")"
        ORDER_FLAGS+=(--mod-t "$t")
      fi
      # Only pass the flag for the non-default polarity, so default runs use
      # exactly the flags they always did.
      [[ -z "$projection" || "$projection" == "at0-cg1" ]] || ORDER_FLAGS+=(--gm-projection "$projection")
      ;;
    mod-minimizer-xor)
      # k, m: run_one's (the caller's) values, else the defaults
      [[ -n "$t" ]] || t="$(bench_py orders mod-t "$order" "${k:-$K}" "${m:-$M}")"
      ORDER_FLAGS+=(-O "$order" --mod-t "$t")
      ;;
    *) die "unknown order: $order (see orders/README.md)" ;;
  esac
  if [[ "$modality" == "regular" ]]; then
    [[ "$order" != mod-minimizer-* ]] || die "the mod-minimizer orders have no regular modality"
    ORDER_FLAGS+=(--regular)
  fi
}

resolve_orders() {
  if [[ "$ORDERS" == "all" ]]; then
    echo "$ALL_ORDERS"
  else
    echo "$ORDERS"
  fi
}
