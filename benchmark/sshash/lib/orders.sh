# Order name -> sshash build flags (see ../../../orders/README.md and
# scripts/orders.py, the Python-side twin of this table). Every order
# (including native) takes a seed with sshash -- unlike kmc's lex, there is
# no seed-less baseline to special-case here.

ALL_ORDERS="native,xor-hash,xor-pcg-hash,greedymini-k13,greedymini-k15"

# Splits a plain-greedymini order name "greedymini-<label>[-at1]" into
# GM_LABEL (table from config.yaml's greedymini_labels, e.g. k13) and
# GM_PROJECTION (sshash --gm-projection: "at0-cg1" default, "at1-cg0" for the
# "-at1" suffix: the flipped DNA->binary polarity).
order_gm_split() {
  local o="$1"
  GM_PROJECTION="at0-cg1"
  if [[ "$o" == *-at1 ]]; then
    GM_PROJECTION="at1-cg0"
    o="${o%-at1}"
  fi
  GM_LABEL="${o#greedymini-}"
  [[ "$GM_LABEL" != "$o" && -n "$GM_LABEL" ]] || die "unknown order: $1 (greedymini-<label>[-at1], see config.yaml's greedymini_labels)"
  return 0
}

order_dir_name() {
  case "$1" in
    native) echo "native" ;;
    xor-hash) echo "xor_hash" ;;
    xor-pcg-hash) echo "xor_pcg_hash" ;;
    # plain greedymini table/projection variants (see order_gm_split)
    greedymini-*) echo "${1//-/_}" ;;
    *) die "unknown order: $1 (see orders/README.md)" ;;
  esac
}

# GreedyMini table to use for order `$1`: greedymini-<label> reads config.yaml's
# greedymini_labels; empty for the other orders.
order_gm_file() {
  case "$1" in
    greedymini-*)
      order_gm_split "$1"
      greedymini_label_file "$GM_LABEL"
      ;;
    *) echo "" ;;
  esac
}

# Fills ORDER_FLAGS with sshash build's -s/-O/--greedy-mini-file/--gm-projection/--regular
# flags for order `$1`, seed `$2` (GreedyMini table `$3`, modality `$4`:
# "canonical" default or "regular", sshash build --regular).
#
# `$2` comes in as config.yaml's hex literal (e.g. "0xDEADBEEF") -- sshash's
# CLI parser (external/pthash/external/cmd_line_parser) reads -s via
# strtoull(value, nullptr, 10), i.e. BASE 10 ONLY, so a "0x..." string
# silently parses as 0 (strtoull stops at the non-decimal 'x'). $((seed))
# converts it to a decimal string sshash actually understands correctly.
build_order_flags() {
  local order="$1" seed="$2" gm_file="$3" modality="${4:-canonical}"
  ORDER_FLAGS=(-s "$((seed))")
  case "$order" in
    native) ;;
    xor-hash|xor-pcg-hash) ORDER_FLAGS+=(-O "$order") ;;
    greedymini-*)
      [[ -f "$gm_file" ]] || die "GreedyMini table not found: $gm_file"
      order_gm_split "$order"
      ORDER_FLAGS+=(-O greedymini --greedy-mini-file "$gm_file")
      # Only pass the flag for the non-default polarity, so default runs use
      # exactly the flags they always did.
      [[ "$GM_PROJECTION" == "at0-cg1" ]] || ORDER_FLAGS+=(--gm-projection "$GM_PROJECTION")
      ;;
    *) die "unknown order: $order (see orders/README.md)" ;;
  esac
  if [[ "$modality" == "regular" ]]; then
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
