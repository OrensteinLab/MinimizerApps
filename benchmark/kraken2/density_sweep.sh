#!/usr/bin/env bash
# density_sweep.sh - minimizer density factor (selected positions / windows x
# (w+1)) of xor and greedymini over a range of minimizer lengths l, using
# kraken2/src/benchmark_minimizers (so it is exactly the density that tool
# reports: expected = i.i.d. random DNA, particular = library FASTAs under -d).
#
#   density_sweep.sh -b BIN -g TABLE.gm -v LABEL -o OUT.csv [options]
#     -b BIN     benchmark_minimizers binary (build: make -C kraken2/src benchmark_minimizers;
#                to compare against an older revision build that revision's binary)
#     -g FILE    binary greedymini table (.gm)
#     -v LABEL   variant label written to the CSV (e.g. prefix, centre)
#     -o FILE    CSV to append to (header written when the file is new)
#     -L "l.."   minimizer lengths (default "15 17 19 21 23 25 27 29 31"); k = l + W - 1
#     -W INT     window w = k - l + 1 (default 5, the table's design window)
#     -n INT     random bases per stream (default 2000000)
#     -S INT     number of independent streams; stream i uses --seed i and a
#                different tie-break/hash toggle (default 5)
#     -s INT     minimizer spaces (default 0; kraken2-build NT default is 7)
#     -d DIR     Kraken DB dir (library/ FASTAs) for the particular density
#                (default: a tiny dummy FASTA, i.e. expected density only; use
#                -S 1 -n 100000 for particular-only runs)
#     -P PROJ    greedymini projection(s), comma-separated: at0-cg1 | at1-cg0 (--gm-projection);
#                default both, one invocation per projection (rows carry the projection
#                in the variant label when more than one is given)
#     -R         also run the xor baselines (xor = kraken2 default, xor-pcg-hash)
# Runs under nice -n 19 with one thread per invocation unless -j is given.
# Aggregate with density_summary.py.
set -euo pipefail
BIN=""; GM=""; LABEL=""; OUT=""; LS="15 17 19 21 23 25 27 29 31"; W=5; PROJ=at0-cg1,at1-cg0
N=2000000; SEEDS=5; SPACES=0; DB=""; WITH_RANDOM=0; JOBS=1
while getopts "b:g:v:o:L:W:n:S:s:d:j:P:R" o; do
  case "$o" in
    b) BIN=$OPTARG;; g) GM=$OPTARG;; v) LABEL=$OPTARG;; o) OUT=$OPTARG;;
    L) LS=$OPTARG;; W) W=$OPTARG;; n) N=$OPTARG;; S) SEEDS=$OPTARG;;
    s) SPACES=$OPTARG;; d) DB=$OPTARG;; j) JOBS=$OPTARG;; P) PROJ=$OPTARG;; R) WITH_RANDOM=1;;
    *) sed -n '2,28p' "$0"; exit 1;;
  esac
done
[[ -x $BIN && -f $GM && -n $LABEL && -n $OUT ]] || { sed -n '2,28p' "$0"; exit 1; }
[[ $JOBS -le 4 ]] || { echo "refusing more than 4 jobs" >&2; exit 1; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
if [[ -z $DB ]]; then
  DB=$TMP/db; mkdir -p "$DB"
  awk 'BEGIN { srand(1); print ">dummy"; s = "";
               for (i = 0; i < 5000; i++) s = s substr("ACGT", int(rand() * 4) + 1, 1);
               print s }' >"$DB/dummy.fna"
fi
SCHEMES="greedymini"; [[ $WITH_RANDOM == 1 ]] && SCHEMES="xor xor-pcg-hash greedymini"
[[ -s $OUT ]] || echo "variant,spaces,scheme,seed_idx,toggle_hex,l,k,w,dens_expected,dens_particular,fac_expected,fac_particular,ns_per_minimizer,sm_mean_exp,sm_max_exp,sm_mean_part,sm_max_part" >"$OUT"
IFS=, read -r -a PROJS <<<"$PROJ"
for l in $LS; do
 for P in "${PROJS[@]}"; do
  VL=$LABEL; [[ ${#PROJS[@]} -gt 1 ]] && VL=$LABEL-$P
  for ((i = 1; i <= SEEDS; i++)); do
    tog=$(printf '0x%x' $(( 0xe37e28c4271b5a2d ^ (i * 0x1d8e4e27c47d124f) )))
    nice -n 19 "$BIN" -S $SCHEMES -g "$GM" --gm-projection "$P" -l "$l" -k $((l + W - 1)) -n "$N" \
      -j "$JOBS" -s "$SPACES" -T "$tog" --seed "$i" --seed-idx "$i" \
      -d "$DB" -o "$TMP/o.csv" >/dev/null 2>"$TMP/err" \
      || { echo "l=$l seed=$i failed:"; tail -3 "$TMP/err"; continue; }
    tail -n +2 "$TMP/o.csv" | sed "s/^/$VL,$SPACES,/" >>"$OUT"
  done
 done
done
