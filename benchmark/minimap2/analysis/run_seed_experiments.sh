#!/usr/bin/env bash
# Mechanism + timing experiments for the scheme analysis, on a fixed random
# read sample per dataset, single thread, seed-0 index of every order:
#   seeds/<ds>/<cap>/<order>.csv      one pass with --dbg-seed-occ (seeds, anchors)
#   timing/<ds>/<cap>/<order>.r<i>.csv  REPS passes without it (per-read time),
#                                      orders shuffled within each replicate
# Usage: run_seed_experiments.sh <dataset> <n_reads> <reps> <cap|default>...
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
ROOT=${MINIMAP2_BENCH_ROOT:-/localdata/shay_sirek/minimap2}
BIN="$REPO/minimap2/minimap2"
PY=${PYTHON:-/usr/bin/python3.11}
ORDERS="native xor xor-pcg-hash greedymini greedymini-at1 mod-minimizer-xor"

ds=$1 n_reads=$2 reps=$3; shift 3
out="$ROOT/results/analysis/$ds"
mkdir -p "$out"
sample="$out/sample_${n_reads}.fa"
preset=map-pb
[ -f "$ROOT/data/$ds/preset" ] && preset=$(tr -d '[:space:]' < "$ROOT/data/$ds/preset")
if [ ! -s "$sample" ]; then
  "$PY" - "$ROOT/data/$ds/reads.fa" "$n_reads" "$sample" <<'EOF'
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
cond=$(cd "$ROOT/results/runs/$ds" && ls -d k*_w* | grep -v _f | head -1)
[ -n "$cond" ] || { echo "error: no runs under $ROOT/results/runs/$ds" >&2; exit 1; }

for cap in "$@"; do
  flag=""
  [ "$cap" != default ] && flag="--cap $cap"
  for o in $ORDERS; do
    idx="$ROOT/results/runs/$ds/$cond/$o/0/index.mmi"   # same index for both caps
    f="$out/seeds/$cap/$o.csv"
    [ -s "$f" ] || "$PY" "$HERE/analysis/seed_stats.py" --minimap2 "$BIN" --index "$idx" --reads "$sample" --out "$f" --seeds --preset "$preset" $flag
  done
  for r in $(seq 1 "$reps"); do
    for o in $(printf '%s\n' $ORDERS | shuf --random-source=<(yes $r)); do
      idx="$ROOT/results/runs/$ds/$cond/$o/0/index.mmi"
      f="$out/timing/$cap/$o.r$r.csv"
      [ -s "$f" ] || "$PY" "$HERE/analysis/seed_stats.py" --minimap2 "$BIN" --index "$idx" --reads "$sample" --out "$f" --preset "$preset" $flag
    done
  done
done
echo "SEED_EXPERIMENTS_DONE $ds"
