# Shared runtime defaults for kmc3_bench.sh (seeds/paths come from config.yaml via load_bench_config).

CURL_RESUME_CODES=(18 26 28 52 55 56 7)
ENA_FASTQ="https://ftp.sra.ebi.ac.uk/vol1/fastq"
TIME_BIN="/usr/bin/time"
REPO_ROOT="$(cd "$HERE/../../kmc" && pwd)"
if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
  PYTHON="${PYTHON:-$REPO_ROOT/.venv/bin/python}"
else
  PYTHON="${PYTHON:-python3}"
fi

CMD="run"
DATASETS=""
KS="28,55"
INPUTS="gzip"
THREADS=12
MEM_OVERRIDE=""
BINS=""
STRICT_MEM=0
RUN_TAG=""
KMC=""
SIGNATURE_ORDER="all"
ORDER_LABEL=""
CUSTOM_ORDER=""
N_SEEDS=1
SIGNATURE_LEN=""
FORCE=0
KEEP_TMP=0
KEEP_DB=0
DATA_DIR=""
RESULTS_DIR=""
DOWNLOAD_TIMEOUT=600
SKIP_MISSING=0
EXPERIMENT=""
RL_SOURCE=""
RL_OUT_ID=""
RL_JOIN=""
RL_TARGET_KMERS=""
RL_AGAINST="f_vesca"
RL_PROBE_READS=200000

# Python helper: run `python -m scripts` from the benchmark directory.
bench_py() {
  (cd "$HERE" && "$PYTHON" -m scripts "$@")
}

# Load shared constants from config.yaml (seeds, paths, dataset list) via a
# small bash-native parse -- config.yaml's shape is fixed and simple enough
# that spawning python (and round-tripping through its own YAML parser) for
# ~4 scalars is unnecessary indirection.
_config_scalar() {
  # top-level "section:\n  key: value" pair
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
  PAPER_CSV="$BENCHMARK_DIR/$data_rel/paper_kmc3.csv"

  ALL_DATASETS="$(_config_scalar defaults datasets)"
  DATASETS="${DATASETS:-$ALL_DATASETS}"
}
