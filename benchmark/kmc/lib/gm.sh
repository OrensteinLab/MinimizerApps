# GreedyMini table paths: w = k - p + 1 (e.g. k=28, p=9 -> w20_k9_canonical*.gm)

signature_len_effective() {
  echo "${SIGNATURE_LEN:-9}"
}

gm_window_for_k() {
  local k="$1" p
  p="$(signature_len_effective)"
  echo $((k - p + 1))
}

gm_table_stem() {
  local k="$1" p w
  p="$(signature_len_effective)"
  w="$(gm_window_for_k "$k")"
  echo "w${w}_k${p}_canonical"
}

# Label passed to kmc -ol (logs / sidecar only); the projection is its own
# meta.json column, not part of the label.
gm_label_for_k() {
  echo "gm_$(gm_table_stem "$1")"
}

resolve_gm_path() {
  local k="$1" seed="${2:-0}" mode="${GM_LUT:-binary}"
  bench_py gm resolve --k "$k" -p "$(signature_len_effective)" --seed "$((seed))" --lut-mode "$mode"
}
