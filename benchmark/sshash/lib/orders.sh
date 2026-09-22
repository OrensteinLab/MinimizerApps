# Order name -> sshash build flags (see ../../../orders/README.md and
# scripts/orders.py, the Python-side twin of this table). Every order
# (including native) takes a seed with sshash -- unlike kmc's lex, there is
# no seed-less baseline to special-case here.
#
# The 4 mod-minimizer-* orders (Groot Koerkamp & Pibiri, WABI 2024) reuse the
# 4 base orders as O_t (the order on t-mers), plus a `t` length: sshash's own
# `-O` flag accepts all 8 names directly and decomposes them internally (see
# sshash/include/custom_order_hasher.hpp::parse_order_name), so this file
# just needs to also thread `--t-length` through for them.

ALL_ORDERS="native,xor-hash,xor-pcg-hash,greedymini,mod-minimizer-native-hash,mod-minimizer-xor-hash,mod-minimizer-pcg-xor-hash,mod-minimizer-greedymini"

order_dir_name() {
  case "$1" in
    native) echo "native" ;;
    xor-hash) echo "xor_hash" ;;
    xor-pcg-hash) echo "xor_pcg_hash" ;;
    greedymini) echo "greedymini" ;;
    mod-minimizer-native-hash) echo "mod_minimizer_native_hash" ;;
    mod-minimizer-xor-hash) echo "mod_minimizer_xor_hash" ;;
    mod-minimizer-pcg-xor-hash) echo "mod_minimizer_pcg_xor_hash" ;;
    mod-minimizer-greedymini) echo "mod_minimizer_greedymini" ;;
    *) die "unknown order: $1 (see orders/README.md)" ;;
  esac
}

# t values to sweep for order `$1`, space-separated, or the sentinel "-" for
# orders with no t at all (the 4 base orders):
#   - mod-minimizer-greedymini is fixed at t=11 (matches w21_k11_canonical.gm's
#     own gm_k exactly, no Extend-k needed there).
#   - the 3 hash-based mod-minimizer orders sweep both t=10 (reproduces the
#     mod-minimizer paper's own published Human/GRCh38 numbers) and t=11
#     (matches the GreedyMini t-order table's length, for a fair comparison
#     against mod-minimizer-greedymini).
order_t_values() {
  case "$1" in
    mod-minimizer-greedymini) echo "11" ;;
    mod-minimizer-*) echo "10 11" ;;
    *) echo "-" ;;
  esac
}

# GreedyMini table to use for order `$1`. `$2` is the CLI/config override for
# plain greedymini only (--greedy-mini-file) -- mod-minimizer-greedymini's
# O_t table always comes from config.yaml's mod_t_gm_file (no CLI override,
# to avoid ambiguity when both orders run in the same sweep).
order_gm_file() {
  case "$1" in
    greedymini) echo "${2:-$DEFAULT_GM_FILE}" ;;
    mod-minimizer-greedymini) echo "$MOD_T_GM_FILE" ;;
    *) echo "" ;;
  esac
}

# Fills the global ORDER_FLAGS array with sshash build's -O/-s/
# --greedy-mini-file/--t-length flags for order `$1` with seed `$2` (GreedyMini
# table `$3`, t-length `$4` for mod-minimizer-* orders).
#
# `$2` comes in as config.yaml's hex literal (e.g. "0xDEADBEEF") -- sshash's
# CLI parser (external/pthash/external/cmd_line_parser) reads -s via
# strtoull(value, nullptr, 10), i.e. BASE 10 ONLY, so a "0x..." string
# silently parses as 0 (strtoull stops at the non-decimal 'x'). $((seed))
# converts it to a decimal string sshash actually understands correctly.
build_order_flags() {
  local order="$1" seed="$2" gm_file="$3" t_len="${4:-}"
  ORDER_FLAGS=(-s "$((seed))")
  case "$order" in
    native) ;;
    xor-hash|xor-pcg-hash) ORDER_FLAGS+=(-O "$order") ;;
    greedymini)
      [[ -n "$gm_file" ]] || gm_file="$DEFAULT_GM_FILE"
      [[ -f "$gm_file" ]] || die "GreedyMini table not found: $gm_file"
      ORDER_FLAGS+=(-O greedymini --greedy-mini-file "$gm_file")
      ;;
    mod-minimizer-native-hash|mod-minimizer-xor-hash|mod-minimizer-pcg-xor-hash)
      [[ -n "$t_len" ]] || die "order $order requires a t value (see order_t_values)"
      ORDER_FLAGS+=(-O "$order" --t-length "$t_len")
      ;;
    mod-minimizer-greedymini)
      [[ -n "$t_len" ]] || die "order $order requires a t value (see order_t_values)"
      [[ -n "$gm_file" ]] || gm_file="$MOD_T_GM_FILE"
      [[ -f "$gm_file" ]] || die "GreedyMini table not found: $gm_file"
      ORDER_FLAGS+=(-O mod-minimizer-greedymini --greedy-mini-file "$gm_file" --t-length "$t_len")
      ;;
    *) die "unknown order: $order (see orders/README.md)" ;;
  esac
}

resolve_orders() {
  if [[ "$ORDERS" == "all" ]]; then
    echo "$ALL_ORDERS"
  else
    echo "$ORDERS"
  fi
}
