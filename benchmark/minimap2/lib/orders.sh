#!/usr/bin/env bash
# Order flag builder for the minimap2 benchmark harness, mirroring
# benchmark/sshash/lib/orders.sh's shape (see orders/README.md for the
# CLI-flag-mapping table across all four apps).

# greedymini-at1 is greedymini with the flipped AT/CG projection
# (--gm-projection at1-cg0, as in sshash); mod-minimizer-xor is the
# strand-symmetric mod-minimizer with xor as O_t (minimap2/order.cpp).
ALL_ORDERS="native,xor,xor-pcg-hash,greedymini,greedymini-at1,mod-minimizer-xor"
DEFAULT_ORDERS="$ALL_ORDERS"

order_needs_gm() {
  case "$1" in
    greedymini|greedymini-at1) return 0 ;;
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
    greedymini|greedymini-at1) echo "$GM_DIR/w10_k15_canonical.gm" ;;
    *) echo "" ;;
  esac
}

# Fills the global array ORDER_FLAGS with the minimap2 CLI args for one
# (order, seed, mod_t) combination. `native` needs no flags at all; a
# greedymini-backed order needs --greedymini-file instead of --order-seed
# (the table's own ranks are already deterministic), and greedymini-at1
# also --gm-projection at1-cg0; a mod-minimizer-* order additionally needs
# --mod-t.
build_order_flags() {
  local order="$1" seed="$2" mod_t="$3"
  ORDER_FLAGS=()
  [ "$order" = "native" ] && return 0
  if order_needs_gm "$order"; then
    ORDER_FLAGS+=(--minimizer-order greedymini --greedymini-file "$(order_gm_file "$order")")
    [ "$order" = "greedymini-at1" ] && ORDER_FLAGS+=(--gm-projection at1-cg0)
    return 0
  fi
  ORDER_FLAGS+=(--minimizer-order "$order" --order-seed "$seed")
  if order_is_mod "$order"; then
    ORDER_FLAGS+=(--mod-t "$mod_t")
  fi
}

order_dir_name() { echo "$1"; }

# native has no order-dependent behavior at all; greedymini's ranks come
# from its table file, not a seed. Every other order takes --order-seed and
# should be swept across the full seed list, not just one.
order_is_seed_dependent() {
  case "$1" in
    native) return 1 ;;
    *) ! order_needs_gm "$1" ;;
  esac
}
