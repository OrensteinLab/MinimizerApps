#!/usr/bin/env bash
# Order flag builder for the minimap2 benchmark harness, mirroring
# benchmark/sshash/lib/orders.sh's shape (see orders/README.md for the
# CLI-flag-mapping table across all four apps).
#
# Orders are the canonical names (benchmark/common/order_names.py) plus the
# baseline `native`. The GreedyMini parameters are NOT part of the order name:
# the table (GM_TABLE, its k is gm_k) and the AT/CG projection (PROJECTION,
# --gm-projection at0-cg1|at1-cg0) are separate parameters. `sweep` runs
# greedymini once per projection in PROJECTIONS, each its own run (own
# directory, own rows). mod-minimizer-xor is the strand-symmetric
# mod-minimizer with xor as O_t (minimap2/order.cpp).
ALL_ORDERS="native,xor,xor-pcg-hash,greedymini,mod-minimizer-xor"
DEFAULT_ORDERS="$ALL_ORDERS"

# GreedyMini table stem under $GM_DIR (its k, gm_k, is parsed from the name);
# the projections `sweep` covers, and the one a single `run` uses.
GM_TABLE=${GM_TABLE:-w10_k15_canonical}
PROJECTIONS=${PROJECTIONS:-at0-cg1,at1-cg0}
PROJECTION=${PROJECTION:-at0-cg1}

order_is_known() {
  case ",$ALL_ORDERS," in
    *",$1,"*) return 0 ;;
    *) return 1 ;;
  esac
}

order_needs_gm() {
  [ "$1" = greedymini ]
}

order_is_mod() {
  case "$1" in
    mod-minimizer-*) return 0 ;;
    *) return 1 ;;
  esac
}

order_gm_file() {
  # Only meaningful when order_needs_gm.
  echo "$GM_DIR/$GM_TABLE.gm"
}

# Stdlib-only helpers from benchmark/common/order_names.py (no parsing here).
_order_names() {
  PYTHONPATH="$REPO_ROOT/benchmark/common${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON" -c "$@"
}

# GreedyMini table k, from the table name (w10_k15_canonical -> 15).
order_gm_k() {
  _order_names 'import sys; from order_names import parse_table_name; print(parse_table_name(sys.argv[1])["gm_k"])' "$GM_TABLE"
}

# Normalised projection ('at0-cg1' | 'at1-cg0'); fails on anything else.
order_projection() {
  _order_names 'import sys; from order_names import projection_from_flag; print(projection_from_flag(sys.argv[1]))' "${1:-$PROJECTION}"
}

# Run-directory segment(s): <order>[/k<gm_k>_<projection>_t<mod_t>] (order_names.order_dir_name).
# $1 order, $2 mod_t. Uses the current GM_TABLE / PROJECTION.
order_dir_name() {
  local order="$1" mod_t="${2:-}" gm_k="" proj=""
  if order_needs_gm "$order"; then
    gm_k=$(order_gm_k) || return 1
    proj=$(order_projection) || return 1
  fi
  _order_names '
import sys
from order_names import order_dir_name
o, k, p, t = sys.argv[1:5]
print(order_dir_name(o, int(k) if k else None, p or None, int(t) if t else None))
' "$order" "$gm_k" "$proj" "$mod_t"
}

# Fills the global array ORDER_FLAGS with the minimap2 CLI args for one
# (order, seed, mod_t) combination. `native` needs no flags at all; greedymini
# needs --greedymini-file and --gm-projection ($PROJECTION), its --order-seed
# toggles the Gamma tie-break (as kraken2's and sshash's greedymini seeds); a
# mod-minimizer-* order additionally needs --mod-t.
build_order_flags() {
  local order="$1" seed="$2" mod_t="$3"
  ORDER_FLAGS=()
  [ "$order" = "native" ] && return 0
  if order_needs_gm "$order"; then
    local proj
    proj=$(order_projection) || return 1
    ORDER_FLAGS+=(--minimizer-order greedymini --greedymini-file "$(order_gm_file "$order")" --order-seed "$seed" --gm-projection "$proj")
    return 0
  fi
  ORDER_FLAGS+=(--minimizer-order "$order" --order-seed "$seed")
  if order_is_mod "$order"; then
    ORDER_FLAGS+=(--mod-t "$mod_t")
  fi
}

# native has no order-dependent behavior at all. Every other order takes
# --order-seed (greedymini: the Gamma tie-break toggle) and is swept across
# the full seed list, not just one.
order_is_seed_dependent() {
  [ "$1" != native ]
}
