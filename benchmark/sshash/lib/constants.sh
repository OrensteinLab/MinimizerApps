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

# Single source of the run defaults (k, m, seeds; the order list is ALL_ORDERS in lib/orders.sh);
# config.yaml only holds the dataset list, seeds, paths and the GreedyMini table labels.
CMD="run"
DATASETS=""
K=31
M=21
ORDERS="all"
MODALITY="canonical"
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
}

# gm file for a plain "greedymini-<label>[-at1]" table variant (the label is
# the part between "greedymini-" and an optional trailing "-at1"), from
# config.yaml's top-level greedymini_labels map -- register a new table
# variant by adding a line there, not by touching this function.
greedymini_label_file() {
  local label="$1" file
  file="$(_config_scalar greedymini_labels "$label")"
  [[ -n "$file" ]] || die "unknown GreedyMini label '$label' (add it to config.yaml's greedymini_labels)"
  echo "$GM_DIR_ABS/$file"
}
