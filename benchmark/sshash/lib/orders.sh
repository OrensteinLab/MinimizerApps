# Order name -> sshash build flags (see ../../../orders/README.md and
# scripts/orders.py, the Python-side twin of this table). Every order
# (including native) takes a seed with sshash -- unlike kmc's lex, there is
# no seed-less baseline to special-case here.

ALL_ORDERS="native,xor-hash,xor-pcg-hash,greedymini"

order_dir_name() {
  case "$1" in
    native) echo "native" ;;
    xor-hash) echo "xor_hash" ;;
    xor-pcg-hash) echo "xor_pcg_hash" ;;
    greedymini) echo "greedymini" ;;
    *) die "unknown order: $1 (use native, xor-hash, xor-pcg-hash, or greedymini)" ;;
  esac
}

# Fills the global ORDER_FLAGS array with sshash build's -O/-s/--greedy-mini-file
# flags for order `$1` with seed `$2` (and, for greedymini, GreedyMini table `$3`).
#
# `$2` comes in as config.yaml's hex literal (e.g. "0xDEADBEEF") -- sshash's
# CLI parser (external/pthash/external/cmd_line_parser) reads -s via
# strtoull(value, nullptr, 10), i.e. BASE 10 ONLY, so a "0x..." string
# silently parses as 0 (strtoull stops at the non-decimal 'x'). $((seed))
# converts it to a decimal string sshash actually understands correctly.
build_order_flags() {
  local order="$1" seed="$2" gm_file="$3"
  ORDER_FLAGS=(-s "$((seed))")
  case "$order" in
    native) ;;
    xor-hash|xor-pcg-hash) ORDER_FLAGS+=(-O "$order") ;;
    greedymini)
      [[ -n "$gm_file" ]] || gm_file="$DEFAULT_GM_FILE"
      [[ -f "$gm_file" ]] || die "GreedyMini table not found: $gm_file"
      ORDER_FLAGS+=(-O greedymini --greedy-mini-file "$gm_file")
      ;;
    *) die "unknown order: $order (use native, xor-hash, xor-pcg-hash, or greedymini)" ;;
  esac
}

resolve_orders() {
  if [[ "$ORDERS" == "all" ]]; then
    echo "$ALL_ORDERS"
  else
    echo "$ORDERS"
  fi
}
