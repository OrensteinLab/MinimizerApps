#!/usr/bin/env bash
# Order flag builder for the minimap2 benchmark harness, mirroring
# benchmark/sshash/lib/orders.sh's shape (see orders/README.md for the
# CLI-flag-mapping table across all four apps).

# No mod-minimizer-greedymini for minimap2 (unlike kraken2/sshash) -- would
# need a second, differently-sized GreedyMini table (gm_k == the
# mod-minimizer t, not k) just for that one variant; dropped instead.
ALL_ORDERS="native,xor-hash,xor-pcg-hash,greedymini,mod-minimizer-xor-hash,mod-minimizer-pcg-xor-hash"
DEFAULT_ORDERS="$ALL_ORDERS"

order_needs_gm() {
  case "$1" in
    greedymini) return 0 ;;
    *) return 1 ;;
  esac
}

order_is_mod() {
  case "$1" in
    mod-minimizer-*) return 0 ;;
    *) return 1 ;;
  esac
}

order_gm_file() {
  # $1: order name. Only meaningful when order_needs_gm.
  case "$1" in
    greedymini) echo "$GM_DIR/w10_k15_canonical.gm" ;;
    *) echo "" ;;
  esac
}

# Fills the global array ORDER_FLAGS with the minimap2 CLI args for one
# (order, seed, mod_t) combination. `native` needs no flags at all; a
# greedymini-backed order needs --greedymini-file instead of --order-seed
# (the table's own ranks are already deterministic); a mod-minimizer-*
# order additionally needs --mod-t.
build_order_flags() {
  local order="$1" seed="$2" mod_t="$3"
  ORDER_FLAGS=()
  [ "$order" = "native" ] && return 0
  ORDER_FLAGS+=(--minimizer-order "$order")
  if order_needs_gm "$order"; then
    ORDER_FLAGS+=(--greedymini-file "$(order_gm_file "$order")")
  else
    ORDER_FLAGS+=(--order-seed "$seed")
  fi
  if order_is_mod "$order"; then
    ORDER_FLAGS+=(--mod-t "$mod_t")
  fi
}

order_dir_name() { echo "$1"; }

# mod-minimizer (Groot Koerkamp & Pibiri, WABI 2024) is only motivated when
# the sampled substring (k) is large *relative to* the window (w) -- that's
# the regime its density bound improves on a classic minimizer in. At w>=k
# (e.g. w=50 with our own k=15) there's no such rationale, so the sweep
# skips mod-minimizer-* there rather than report numbers for a use outside
# the scheme's own design regime (independent of, and in addition to, the
# cross-strand-symmetry limitation noted in minimap2/order.cpp).
order_relevant_for_w() {
  local order="$1"
  if order_is_mod "$order" && [ "$W" -ge "$K" ]; then
    return 1
  fi
  return 0
}
