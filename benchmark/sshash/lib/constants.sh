# Shared runtime defaults for sshash_bench.sh (seeds/paths come from config.yaml via load_bench_config).

CURL_RESUME_CODES=(18 26 28 52 55 56 7)
ENA_FASTQ="https://ftp.sra.ebi.ac.uk/vol1/fastq"
TIME_BIN="/usr/bin/time"
REPO_ROOT="$(cd "$HERE/../../sshash" && pwd)"
# The harness's own venv lives under benchmark/sshash/ (not sshash/, which is
# upstream-tracked source we don't touch -- see Makefile's `venv` target).
if [[ -x "$HERE/.venv/bin/python" ]]; then
  PYTHON="${PYTHON:-$HERE/.venv/bin/python}"
else
  PYTHON="${PYTHON:-python3}"
fi

CMD="run"
DATASETS=""
K=31
M=21
ORDERS="all"
MODALITY="canonical"
GM_FILE=""
THREADS=8
N_SEEDS=5
FORCE=0
KEEP_INDEX=0
RUN_TAG=""
SSHASH=""
DATA_DIR=""
RESULTS_DIR=""
DOWNLOAD_TIMEOUT=600
SKIP_MISSING=0
EXPERIMENT=""
DATASET_FILE=""

# Python helper: run `python -m scripts` from the benchmark directory.
bench_py() {
  (cd "$HERE" && "$PYTHON" -m scripts "$@")
}

# Load shared constants from config.yaml (seeds, paths, dataset list) via a
# small bash-native parse -- same approach as benchmark/kmc/lib/constants.sh
# (config.yaml's shape is fixed and simple enough that spawning python for a
# handful of scalars is unnecessary indirection).
_config_scalar() {
  local section="$1:" key="$2:"
  awk -v section="$section" -v key="$key" '
    $0 == section { insec = 1; next }
    insec && /^[^[:space:]]/ { insec = 0 }
    insec {
      line = $0
      sub(/^[[:space:]]+/, "", line)
      if (index(line, key) == 1) {
        val = substr(line, length(key) + 1)
        sub(/^[[:space:]]+/, "", val)
        sub(/[[:space:]]*(#.*)?$/, "", val)
        gsub(/^"|"$/, "", val)
        print val
        exit
      }
    }
  ' "$HERE/config.yaml"
}

load_bench_config() {
  BENCH_SEEDS=()
  local seed
  while IFS= read -r seed; do
    BENCH_SEEDS+=("$seed")
  done < <(awk '
    /^seeds:/ { insec = 1; next }
    insec && /^[^[:space:]]/ { insec = 0 }
    insec {
      line = $0
      sub(/^[[:space:]]*-[[:space:]]*/, "", line)
      sub(/[[:space:]]*(#.*)?$/, "", line)
      if (line != "") print line
    }
  ' "$HERE/config.yaml")

  local data_rel
  data_rel="$(_config_scalar paths data)"
  data_rel="${data_rel:-data}"
  DATASETS_CSV="$BENCHMARK_DIR/$data_rel/datasets.csv"

  ALL_DATASETS="$(_config_scalar defaults datasets)"
  DATASETS="${DATASETS:-$ALL_DATASETS}"

  local gm_dir
  gm_dir="$(_config_scalar paths gm_dir)"
  gm_dir="${gm_dir:-../../orders/GreedyMini}"
  GM_DIR_ABS="$(cd "$BENCHMARK_DIR" && cd "$gm_dir" && pwd)"

  local default_gm
  default_gm="$(_config_scalar defaults gm_file)"
  default_gm="${default_gm:-w11_k15_canonical.gm}"
  DEFAULT_GM_FILE="$GM_DIR_ABS/$default_gm"

  # mod-minimizer-greedymini's O_t table for t=11, separate from plain
  # greedymini's own (Extend-k'd) table above -- see config.yaml.
  local mod_t_gm
  mod_t_gm="$(_config_scalar defaults mod_t_gm_file)"
  mod_t_gm="${mod_t_gm:-w21_k11_canonical.gm}"
  MOD_T_GM_FILE="$GM_DIR_ABS/$mod_t_gm"

  # ...and for t=10.
  local mod_t_gm_t10
  mod_t_gm_t10="$(_config_scalar defaults mod_t_gm_file_t10)"
  mod_t_gm_t10="${mod_t_gm_t10:-w22_k10_canonical.gm}"
  MOD_T_GM_FILE_T10="$GM_DIR_ABS/$mod_t_gm_t10"
}

# gm file for a "mod-minimizer-greedymini-<label>" variant, from config.yaml's
# top-level mod_t_gm_labels map -- register a new variant by adding a line
# there, not by touching this function.
mod_t_gm_label_file() {
  local label="$1" file
  file="$(_config_scalar mod_t_gm_labels "$label")"
  [[ -n "$file" ]] || die "unknown GreedyMini label '$label' (add it to config.yaml's mod_t_gm_labels)"
  echo "$GM_DIR_ABS/$file"
}
