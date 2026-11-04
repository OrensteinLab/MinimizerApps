#!/usr/bin/env bash
# Shell defaults, mirroring benchmark/sshash/lib/constants.sh's shape.
# config.yaml is the single source of truth for the seed list (the same
# five seeds as kmc/sshash); load_bench_config() below pulls just the fields this harness needs
# out of it with awk, same approach as sshash's own constants.sh.

HERE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$HERE_DIR/../.." && pwd)"

# Python >= 3.9: prefer the `minimap2_bench` conda env's once set up (setup.sh/environment.yml
# -- `conda env list` shows where it actually lands, typically
# ~/.conda/envs/minimap2_bench), else fall back to the newest system python3.
_CONDA_ENV_DIR="$(command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | awk '$1=="minimap2_bench"{print $NF}')"
if [ -z "$_CONDA_ENV_DIR" ] && [ -x "$HOME/.conda/envs/minimap2_bench/bin/python3" ]; then
  _CONDA_ENV_DIR="$HOME/.conda/envs/minimap2_bench"
fi
if [ -n "$_CONDA_ENV_DIR" ] && [ -x "$_CONDA_ENV_DIR/bin/python3" ]; then
  PYTHON="${PYTHON:-$_CONDA_ENV_DIR/bin/python3}"
  PBSIM="${PBSIM:-$_CONDA_ENV_DIR/bin/pbsim}"
  PBSIM1="${PBSIM1:-$_CONDA_ENV_DIR/opt/pbsim1/bin/pbsim}"   # PBSIM 1.0.3 (setup.sh)
  SAMTOOLS="${SAMTOOLS:-$_CONDA_ENV_DIR/bin/samtools}"
elif command -v /usr/bin/python3.11 >/dev/null 2>&1; then
  PYTHON="${PYTHON:-/usr/bin/python3.11}"
else
  PYTHON="${PYTHON:-python3}"
fi
# Without the conda env, fall back to whatever is on PATH (empty if absent:
# cmd_pbsim_reads / export report a missing tool themselves) so `set -u`
# callers never hit an unbound variable.
SAMTOOLS="${SAMTOOLS:-$(command -v samtools 2>/dev/null || true)}"
PBSIM="${PBSIM:-$(command -v pbsim 2>/dev/null || true)}"
PBSIM1="${PBSIM1:-}"
MINIMAP2_BIN="${MINIMAP2_BIN:-$REPO_ROOT/minimap2/minimap2}"

# k defaults to the dataset preset's own (preset_default_k), as the Winnowmap
# paper ran minimap2 with `-cx map-pb` / `-cx map-ont`; K=... overrides it.
# mod-minimizer's t defaults to the paper's r=4 formula t = 4 + (k-4) % w
# (k=19,w=10 -> 9; k=15 -> 5); MOD_T=... overrides it.
K=${K:-}
W=${W:-10}
MOD_T=${MOD_T:-}
# Map-time repeat-filter cap (minimap2 -f INT, absolute). Empty = minimap2's
# own per-index default (a quantile of that index's own occurrence counts, so
# it drifts with the scheme); set e.g. MID_OCC=1244 (native's own value on
# chm13_chrX) to hold it fixed across schemes as a control.
MID_OCC=${MID_OCC:-}

# minimap2 preset for the current dataset: $PRESET if set, else the one-word
# data/<dataset>/preset file (e.g. map-ont for ONT reads), else map-pb (PacBio
# CLR, which also turns on HPC k-mers). Used for both build and map.
dataset_preset() {
  local f="$MINIMAP2_DATA_DIR/$DATASET/preset"
  if [ -n "${PRESET:-}" ]; then echo "$PRESET"
  elif [ -f "$f" ]; then tr -d '[:space:]' < "$f"
  else echo map-pb
  fi
}

# minimap2's own k for a preset (minimap2/options.c mm_set_opt: map-pb and
# the other HPC PacBio presets k=19, map-ont/default k=15).
preset_default_k() {
  case "$1" in
    map-pb|map10k|ava-pb) echo 19 ;;
    *) echo 15 ;;
  esac
}
DATASET=${DATASET:-smoke}
GM_DIR=${GM_DIR:-"$REPO_ROOT/orders/GreedyMini"}
# minimap2 -t for build and map; the pinning reserves THREADS + 2 cores (16 = the paper's runs).
THREADS=${THREADS:-16}
FORCE=${FORCE:-0}

load_bench_config() {
  local cfg="$HERE_DIR/config.yaml"
  [ -f "$cfg" ] || return 0
  # First 5 "  - 0x..." lines under `seeds:` -> a bash array of seeds.
  mapfile -t ALL_SEEDS < <(awk '/^seeds:/{f=1;next} f&&/^[a-z]/{f=0} f&&/- /{gsub(/^ *- */,"");print}' "$cfg")
}
