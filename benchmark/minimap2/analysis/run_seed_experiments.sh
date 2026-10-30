#!/usr/bin/env bash
# Mechanism + timing experiments for the scheme analysis, on a fixed random
# read sample per dataset, single thread, seed-0 index of every scheme:
#   seeds/<ds>/<cap>/<variant>.csv      one pass with --dbg-seed-occ (seeds, anchors)
#   timing/<ds>/<cap>/<variant>.r<i>.csv  REPS passes without it (per-read time),
#                                      variants shuffled within each replicate
# A variant is a scheme's run-directory name (scheme_names.scheme_dir_name), e.g.
# native, xor, greedymini/k15_at1-cg0, mod-minimizer-xor/t9,
# mod-minimizer-greedymini/k9_at0-cg1_t9; in file names "/" is written "__". Every
# variant with a seed-0 index under the dataset's run directory is used (the
# greedymini schemes once per projection). Existing outputs are kept, so a re-run
# only adds the variants that lack them.
# Usage: run_seed_experiments.sh <dataset> <n_reads> <reps> <cap|default>...
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
# Same .env / PYTHON / MINIMAP2_BIN resolution as minimap2_bench.sh.
source "$HERE/lib/env.sh"
source "$HERE/lib/constants.sh"
source "$HERE/../common/wait_low_load.sh"
source "$HERE/../common/pin_cpus.sh"   # every timed step is pinned (seed_stats.py runs minimap2 -t 1)
load_bench_env
DATA="$MINIMAP2_DATA_DIR"
RESULTS="$MINIMAP2_RESULTS_DIR"
BIN="$MINIMAP2_BIN"
PY="$PYTHON"

ds=$1 n_reads=$2 reps=$3; shift 3
# per-read csvs live under results/csvs/analysis/<ds>; the sampled-reads fasta
# (non-csv) under results/analysis/<ds> -- the one extra dir under results/.
out="$RESULTS/csvs/analysis/$ds"
sdir="$RESULTS/analysis/$ds"
mkdir -p "$out" "$sdir"
sample="$sdir/sample_${n_reads}.fa"
preset=map-pb
[ -f "$DATA/$ds/preset" ] && preset=$(tr -d '[:space:]' < "$DATA/$ds/preset")
if [ ! -s "$sample" ]; then
  "$PY" - "$DATA/$ds/reads.fa" "$n_reads" "$sample" <<'EOF'
import random, sys
src, n, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
names, seqs, cur = [], {}, None
for line in open(src):
    if line.startswith(">"):
        cur = line[1:].split()[0]; names.append(cur); seqs[cur] = []
    else:
        seqs[cur].append(line.strip())
pick = sorted(random.Random(7).sample(names, min(n, len(names))))  # all reads when the set is small (D3)
with open(dst, "w") as f:
    for r in pick:
        f.write(f">{r}\n{''.join(seqs[r])}\n")
EOF
fi

# The dataset's single k<k>_w<w> run dir (k is its preset's own).
cond=$(cd "$RESULTS/runs/$ds" && ls -d k*_w* | grep -v _f | head -1)
[ -n "$cond" ] || { echo "error: no runs under $RESULTS/runs/$ds" >&2; exit 1; }

# schemes in canonical sequence first (stats.variant_key: scheme_names.sort_index), then the parameter directory
VARIANTS=$(cd "$RESULTS/runs/$ds/$cond" && find . -path '*/0/index.mmi' | sed 's#^\./##; s#/0/index.mmi$##' \
  | PYTHONPATH="$HERE/analysis" "$PY" -c '
import sys
from stats import variant_key
print("\n".join(sorted(sys.stdin.read().split(), key=variant_key)))')
[ -n "$VARIANTS" ] || { echo "error: no seed-0 indexes under $RESULTS/runs/$ds/$cond" >&2; exit 1; }

for cap in "$@"; do
  flag=""
  [ "$cap" != default ] && flag="--cap $cap"
  for v in $VARIANTS; do
    idx="$RESULTS/runs/$ds/$cond/$v/0/index.mmi"   # same index for both caps
    f="$out/seeds/$cap/${v//\//__}.csv"
    [ -s "$f" ] || { wait_low_load; bench_pin_cpus 1; ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "$PY" "$HERE/analysis/seed_stats.py" --minimap2 "$BIN" --index "$idx" --reads "$sample" --out "$f" --seeds --preset "$preset" $flag; }
  done
  for r in $(seq 1 "$reps"); do
    for v in $(printf '%s\n' $VARIANTS | shuf --random-source=<(yes $r)); do
      idx="$RESULTS/runs/$ds/$cond/$v/0/index.mmi"
      f="$out/timing/$cap/${v//\//__}.r$r.csv"
      [ -s "$f" ] || { wait_low_load; bench_pin_cpus 1; ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "$PY" "$HERE/analysis/seed_stats.py" --minimap2 "$BIN" --index "$idx" --reads "$sample" --out "$f" --preset "$preset" $flag; }
    done
  done
done
echo "SEED_EXPERIMENTS_DONE $ds"
