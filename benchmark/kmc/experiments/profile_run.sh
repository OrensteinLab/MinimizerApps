#!/usr/bin/env bash
# Signature-cost profile: how much computing the signature of every k-mer costs per base with each
# order KMC3 runs, as in its splitter. tools/bench_signature_impls (built like kmc: -O3 -fsigned-char)
# runs kmc's own signature code (CMmer) with the default order and with the GreedyMini order through
# the DNA-alphabet and the binary table, and the minimum over each window of signatures, on the first
# N_BASES bases of a dataset's reads, one core, REPS times (median reported). One run per dataset and
# seed of the full run (the seed selects the GreedyMini DNA table and the binary table's tie-break),
# written to
#   <results>/profile/signature_<dataset>_<seed>.txt
# Then it re-exports (python -m scripts export), which adds the signature_cost sheet to
# kmc3_results.xlsx (mean / std over seeds, ratio to the default order). The runs follow the
# measurement protocol (node gate, pinning, memory binding); no run of the full run is touched.
#
# It works on a mirror of the full run, RESULTS_NAME (default results-prof), made from SOURCE_RESULTS
# (default results) by common/mirror_results.py: large files hard-linked, small ones copied, so the
# profiles and the re-export are written there and the source result set is never modified.
#
#   experiments/profile_run.sh
#   DATASETS="f_vesca h_sapiens_2" N_BASES=2e8 REPS=7 FORCE=1 NO_EXPORT=1 ...
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$HERE"
if [ -z "${KMC3_BENCH_ROOT:-}" ] && [ -f .env ]; then
  KMC3_BENCH_ROOT="$(sed -n 's/^KMC3_BENCH_ROOT=//p' .env | tr -d "\"' ")"
fi
: "${KMC3_BENCH_ROOT:?KMC3_BENCH_ROOT is not set (benchmark/kmc/.env)}"
export KMC3_BENCH_ROOT
RES="$KMC3_BENCH_ROOT/${RESULTS_NAME:-results-prof}"
export KMC3_RESULTS_DIR="$RES"
DATA="${KMC3_DATA_DIR:-$KMC3_BENCH_ROOT/data}"
DATASETS="${DATASETS:-f_vesca g_gallus m_balbisiana h_sapiens_2 h_sapiens_3}"
SEEDS="${SEEDS:-0xDEADBEEF 0xC0FFEE01 0x31415926 0x27182818 0x16180339}"   # the full run's 5 (config.yaml)
N_BASES="${N_BASES:-2e8}"
REPS="${REPS:-7}"
P=9; L=28
GM="$REPO/orders/GreedyMini/w20_k9_canonical"
PYTHON="${PYTHON:-$REPO/kmc/.venv/bin/python}"
"$PYTHON" "$HERE/../common/mirror_results.py" "$KMC3_BENCH_ROOT/${SOURCE_RESULTS:-results}" "$RES" || exit 1

# shellcheck source=../common/wait_low_load.sh
source "$HERE/../common/wait_low_load.sh"
# shellcheck source=../common/pin_cpus.sh
source "$HERE/../common/pin_cpus.sh"

make -C tools bench_signature_impls >&2 || { echo "error: build of bench_signature_impls failed" >&2; exit 1; }
BIN="$HERE/tools/bench_signature_impls"
mkdir -p "$RES/profile"
tmp="$(mktemp -d "${TMPDIR:-/tmp}/kmc_profile.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

fail=0
for ds in $DATASETS; do
  # the first N_BASES bases of the dataset's first input file (plain FASTQ for the tool)
  first="$(head -1 "$DATA/$ds/files_gz.lst")"
  awk -v n="$N_BASES" 'NR%4==2{b+=length($0)} {print} NR%4==0 && b>=n{exit}' < <(zcat "$first") > "$tmp/$ds.fq"
  for seed in $SEEDS; do
    out="$RES/profile/signature_${ds}_${seed}.txt"
    [ -s "$out" ] && grep -q "^TIMING" "$out" && [ "${FORCE:-0}" != 1 ] && continue
    dna="${GM}_dna_s$(printf '0x%x' "$seed").gm"
    [ "$seed" = 0xDEADBEEF ] && dna="${GM}_dna.gm"      # the seed baked into the unsuffixed table
    [ -f "$dna" ] || { echo "error: no DNA table $dna" >&2; fail=1; continue; }
    wait_low_load 4 0 2
    BENCH_PIN_SPARE=1 bench_pin_cpus 1
    ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "$BIN" "$P" "$GM.gm" "$dna" "$seed" "$L" "$tmp/$ds.fq" "$N_BASES" "$REPS" \
      > "$out" 2> "$out.log" || fail=1
    echo "$(date +%T) $ds $seed" >&2
  done
  rm -f "$tmp/$ds.fq"
done

if [ "${NO_EXPORT:-0}" != 1 ]; then
  (cd "$HERE" && "$PYTHON" -m scripts export --results-dir "$RES") >&2 || fail=1
fi
exit "$fail"
