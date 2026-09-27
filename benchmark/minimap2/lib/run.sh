#!/usr/bin/env bash
# Single-config runner: build index -> map reads -> score accuracy -> meta.json.
set -uo pipefail

run_one() {
  local order="$1" seed="${2:-0}" seed_idx="${3:-0}"
  local run_dir="$MINIMAP2_RESULTS_DIR/runs/$DATASET/k${K}_w${W}/$order/$seed_idx"
  mkdir -p "$run_dir"

  local data_dir="$MINIMAP2_DATA_DIR/$DATASET"
  local genome="$data_dir/genome.fa"
  local reads="$data_dir/reads.fa"
  local truth="$data_dir/truth.tsv"
  if [ ! -f "$genome" ] || [ ! -f "$reads" ]; then
    echo "error: dataset '$DATASET' not found under $data_dir (run: $0 make-dataset)" >&2
    return 1
  fi

  build_order_flags "$order" "$seed" "$MOD_T"

  local mmi="$run_dir/index.mmi"
  local build_time="$run_dir/build_time.txt"
  local build_log="$run_dir/build.log"
  local paf="$run_dir/map.paf"
  local map_time="$run_dir/map_time.txt"
  local map_log="$run_dir/map.log"
  local acc_json="$run_dir/accuracy.json"
  local meta_json="$run_dir/meta.json"

  # -x map-pb: the preset both Winnowmap's own paper (Jain et al. 2020) and
  # minimap2's own docs recommend for PacBio reads -- every dataset here is
  # PacBio-simulated (see lib/pbsim_reads.py). Beyond scoring/chaining
  # defaults tuned for PacBio's error profile, this enables HPC
  # (homopolymer-compressed) k-mers (`-H`), which we'd otherwise silently
  # never turn on despite benchmarking exactly the read type it's meant
  # for. -k/-w passed *after* still win (confirmed via `-d` build log's own
  # "kmer size: $K; skip: $W" report) -- only k/w/H need to be at build
  # time (baked into the index); the other -x map-pb effects (scoring,
  # chaining bandwidth, etc.) are map-time parameters and must be repeated
  # on the map invocation too, where the loaded index's own k/w/H silently
  # win over the preset's (a harmless, expected stderr warning).
  echo "== [$order] building index (k=$K w=$W ${ORDER_FLAGS[*]:-}) ==" >&2
  if ! /usr/bin/time -v -o "$build_time" \
      "$MINIMAP2_BIN" -d "$mmi" -x map-pb -k "$K" -w "$W" "${ORDER_FLAGS[@]}" -t "$THREADS" "$genome" 2> "$build_log"; then
    echo "  build FAILED, see $build_log" >&2
    tail -5 "$build_log" >&2
    return 1
  fi

  echo "== [$order] mapping ${reads##*/} ==" >&2
  if ! /usr/bin/time -v -o "$map_time" \
      "$MINIMAP2_BIN" -x map-pb -c -t "$THREADS" "$mmi" "$reads" > "$paf" 2> "$map_log"; then
    echo "  map FAILED, see $map_log" >&2
    tail -5 "$map_log" >&2
    return 1
  fi

  # "distinct minimizers" is NOT the density_factor numerator -- it collapses
  # every repeated occurrence of the same minimizer key down to one count,
  # so it always undercounts total selected positions on a genome with any
  # repeat content (see mm_idx_stat in index.c: it separately tracks `sum`,
  # the true total-occurrences count, and prints it only via "average
  # occurrences"/"average spacing", not directly). Use "average spacing"
  # (genome_len / total_occurrences) for density_factor instead --
  # write_meta.py's own comment has the full reasoning.
  local density
  density=$(grep -o 'distinct minimizers: [0-9]*' "$build_log" | grep -o '[0-9]*' | head -1)
  local avg_occ
  avg_occ=$(grep -o 'average occurrences: [0-9.]*' "$build_log" | grep -o '[0-9.]*$' | head -1)
  local avg_spacing
  avg_spacing=$(grep -o 'average spacing: [0-9.]*' "$build_log" | grep -o '[0-9.]*$' | head -1)

  "$PYTHON" "$HERE_DIR/lib/score_accuracy.py" --paf "$paf" --truth "$truth" --json-out "$acc_json" >&2

  "$PYTHON" "$HERE_DIR/lib/write_meta.py" \
    --order "$order" --seed "$seed" --seed-idx "$seed_idx" --k "$K" --w "$W" --mod-t "$MOD_T" --dataset "$DATASET" \
    --build-time "$build_time" --map-time "$map_time" \
    --density "${density:-0}" --avg-occurrences "${avg_occ:-0}" --avg-spacing "${avg_spacing:-0}" \
    --accuracy "$acc_json" --out "$meta_json" >/dev/null

  echo "  -> $meta_json"
}
