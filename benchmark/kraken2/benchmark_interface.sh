# Shared scheme / DB-size CLI for run_benchmark.sh and benchmark_minimizers.sh.
# Source after benchmark_setup.sh with CODEBASE set:
#   source "${SCRIPT_DIR}/benchmark_interface.sh"
#   bench_iface_init

# shellcheck shell=bash

_bench_iface_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BENCH_SCHEME_CHOICES="Random | GreedyMini | ModRandom | ModGreedyMini | RandomPcg"

bench_iface_init() {
    local gm_dir="${CODEBASE}/../orders/GreedyMini"
    BASE_DIR="${KRAKEN2_BENCHMARK_DIR:-${_bench_iface_dir}/benchmark_data}"
    GM_FILE="${GM_FILE:-${gm_dir}/w5_k15_canonical.gm}"
    MOD_GM_FILE="${MOD_GM_FILE:-${gm_dir}/w5_k6_canonical.gm}"
    DEFAULT_SCHEMES=("Random" "GreedyMini" "ModRandom" "ModGreedyMini" "RandomPcg")
    SCHEME_FILTER_SET=0
    USE_SAMPLE=0
    USE_SMOKE_SAMPLE=0
    SIZE_EXPLICIT=0
    SKIP_ARGOS=0
    USE_FULL_NCBI=0
    USE_STREX_LIB=0
}

bench_scheme_help_line() {
    echo "  ${BENCH_SCHEME_CHOICES}"
}

# Common -o/-d/-S/-g/-G/-M/-m/-z option text (callers add their own extras).
bench_common_options_help() {
    cat <<EOF
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCHMARK_DIR or ./benchmark_data)
  -d DATA_SRC   (ignored; kept for compat)
  -S SCHEME...  Run only these schemes (space-separated; repeatable):
                  ${BENCH_SCHEME_CHOICES}
                Examples: -S Random ModRandom
                          -S "Random ModRandom RandomPcg"
                          -S Random -S RandomPcg
  -g GM_FILE    GreedyMini binary ℓ-mer rank file (.gm) for GreedyMini
  -M MOD_GM     binary w5_k6_canonical.gm for ModGreedyMini (or env MOD_GM_FILE)
  -m            NCBI ~10% subsample (db/standard_p10) → runs/standard_p10/
  -z            NCBI ~0.5% subsample (db/standard_p0.5) → runs/standard_p0.5/
  -N            NCBI 100% (db/standard) → runs/standard/ (same as omitting -m/-z)
  -E            Wood 2019 library (db/strain_exclusion) → runs/strain_exclusion/
                Mason + FDA-ARGOS queries (Tables S1 and S5); not with -m/-z/-N
  -a            With -E, skip FDA-ARGOS (Mason / Table S1 only)
                Default library: db/standard → runs/standard/ (unlabeled speed
                FASTQs under data/speed/). Paper accuracy is -E only.
EOF
}

# Append scheme name(s) to DEFAULT_SCHEMES; first -S replaces the default list.
bench_append_schemes() {
    local -a names=()
    local name
    if [[ "${SCHEME_FILTER_SET}" -eq 0 ]]; then
        DEFAULT_SCHEMES=()
        SCHEME_FILTER_SET=1
    fi
    read -r -a names <<< "$*"
    for name in "${names[@]}"; do
        [[ -z "${name}" ]] && continue
        case "${name}" in
            RandomXor) name="Random" ;;
            ModRandomXor) name="ModRandom" ;;
            RandomPCG|randompcg|random_pcg) name="RandomPcg" ;;
        esac
        case "${name}" in
            Random|GreedyMini|ModRandom|ModGreedyMini|RandomPcg) ;;
            *)
                echo "Error: Unknown scheme '${name}'." >&2
                bench_scheme_help_line >&2
                exit 1
                ;;
        esac
        DEFAULT_SCHEMES+=("${name}")
    done
}

# Handle shared getopts letters except -S. Returns 0 if handled.
bench_handle_common_opt() {
    case "$1" in
        o) BASE_DIR="${OPTARG}" ;;
        d) : ;; # compat no-op
        g) GM_FILE="${OPTARG}" ;;
        M) MOD_GM_FILE="${OPTARG}" ;;
        m) USE_SAMPLE=1; SIZE_EXPLICIT=1; USE_FULL_NCBI=0 ;;
        z) USE_SMOKE_SAMPLE=1; SIZE_EXPLICIT=1; USE_FULL_NCBI=0 ;;
        N) USE_FULL_NCBI=1; SIZE_EXPLICIT=1 ;;
        E) USE_STREX_LIB=1 ;;
        a) SKIP_ARGOS=1 ;;
        *) return 1 ;;
    esac
    return 0
}

# Handle -S: call as bench_handle_S_opt "$@" from inside the getopts loop
# (OPTARG / OPTIND from getopts; "$@" must be the script's original args).
bench_handle_S_opt() {
    bench_append_schemes "${OPTARG}"
    while [[ ${OPTIND} -le $# ]]; do
        local next="${!OPTIND}"
        [[ "${next}" == -* ]] && break
        bench_append_schemes "${next}"
        OPTIND=$((OPTIND + 1))
    done
}

bench_validate_scheme_filter() {
    if [[ "${SCHEME_FILTER_SET}" -eq 1 && "${#DEFAULT_SCHEMES[@]}" -eq 0 ]]; then
        echo "Error: -S requires at least one scheme." >&2
        exit 1
    fi
}

# Library ids: standard | standard_p10 | standard_p0.5 | strain_exclusion
# Shared FASTA: db/<id>/   Scheme indexes: runs/<id>/{scheme}/{seed}/

BENCH_LIB_IDS=(standard standard_p10 standard_p0.5 strain_exclusion)

# Mason + FDA-ARGOS accuracy mode (-E).
bench_is_paper_accuracy_mode() {
    [[ "${USE_STREX_LIB}" -eq 1 ]]
}

bench_is_kraken_lib_dir() {
    [[ -d "$1/library" || -d "$1/taxonomy" ]]
}

# Prefer canonical if it is already a Kraken library; else first legacy hit; else canonical.
bench_pick_kraken_lib() {
    local canonical="$1"
    shift
    if bench_is_kraken_lib_dir "${canonical}"; then
        echo "${canonical}"
        return
    fi
    local d
    for d in "$@"; do
        if bench_is_kraken_lib_dir "${d}"; then
            echo "${d}"
            return
        fi
    done
    echo "${canonical}"
}

bench_pick_dir() {
    local canonical="$1"
    shift
    if [[ -d "${canonical}" ]]; then
        echo "${canonical}"
        return
    fi
    local d
    for d in "$@"; do
        if [[ -d "${d}" ]]; then
            echo "${d}"
            return
        fi
    done
    echo "${canonical}"
}

bench_canonical_shared_db() {
    echo "${BASE_DIR}/db/$1"
}

bench_canonical_run_root() {
    echo "${BASE_DIR}/runs/$1"
}

bench_shared_db_for() {
    case "$1" in
        standard)
            bench_pick_kraken_lib "${BASE_DIR}/db/standard" "${BASE_DIR}/db"
            ;;
        standard_p10)
            bench_pick_kraken_lib "${BASE_DIR}/db/standard_p10" "${BASE_DIR}/sample_db"
            ;;
        standard_p0.5)
            bench_pick_kraken_lib "${BASE_DIR}/db/standard_p0.5" "${BASE_DIR}/sample_db_smoke"
            ;;
        strain_exclusion)
            bench_pick_kraken_lib "${BASE_DIR}/db/strain_exclusion" \
                "${BASE_DIR}/db_strain_exclusion"
            ;;
        *)
            echo "Error: unknown library id '$1'" >&2
            return 1
            ;;
    esac
}

bench_run_root_for() {
    case "$1" in
        standard)
            bench_pick_dir "${BASE_DIR}/runs/standard" \
                "${BASE_DIR}/runs/standard_p100" \
                "${BASE_DIR}/standard/p100" "${BASE_DIR}/p100"
            ;;
        standard_p10)
            bench_pick_dir "${BASE_DIR}/runs/standard_p10" \
                "${BASE_DIR}/standard/p10" "${BASE_DIR}/p10"
            ;;
        standard_p0.5)
            bench_pick_dir "${BASE_DIR}/runs/standard_p0.5" \
                "${BASE_DIR}/standard/p0.5" "${BASE_DIR}/p0.5"
            ;;
        strain_exclusion)
            bench_pick_dir "${BASE_DIR}/runs/strain_exclusion" \
                "${BASE_DIR}/strain-exclusion"
            ;;
        *)
            echo "Error: unknown library id '$1'" >&2
            return 1
            ;;
    esac
}

bench_lib_id_from_flags() {
    if [[ "${USE_STREX_LIB}" -eq 1 ]]; then
        echo "strain_exclusion"
    elif [[ "${USE_SMOKE_SAMPLE}" -eq 1 ]]; then
        echo "standard_p0.5"
    elif [[ "${USE_SAMPLE}" -eq 1 ]]; then
        echo "standard_p10"
    else
        echo "standard"
    fi
}

# Sets LIB_ID, SHARED_DB, RUN_LAYOUT (same as LIB_ID).
bench_resolve_shared_db() {
    local db="${SHARED_DB:-}"
    if [[ -n "${db}" ]]; then
        case "${db%/}" in
            */db/strain_exclusion|*/db_strain_exclusion) USE_STREX_LIB=1 ;;
        esac
    fi
    if bench_is_paper_accuracy_mode && [[ "${SIZE_EXPLICIT}" -eq 1 ]]; then
        echo "Error: -E (strain-exclusion library) cannot be combined with -m/-z/-N." >&2
        exit 1
    fi
    LIB_ID="$(bench_lib_id_from_flags)"
    RUN_LAYOUT="${LIB_ID}"
    SHARED_DB="${SHARED_DB:-$(bench_shared_db_for "${LIB_ID}")}"
}

bench_run_root() {
    bench_run_root_for "${LIB_ID}"
}

# SCHEME, SCHEME_BASE, SCHEME_ARGS, SCHEME_NEEDS_GM for one -S name.
bench_configure_scheme() {
    local scheme="$1"
    case "${scheme}" in
        Random|RandomXor)
            SCHEME="Random"
            SCHEME_BASE="Random"
            SCHEME_ARGS="--sketching-scheme Minimizer --order RandomXor"
            SCHEME_NEEDS_GM=0
            ;;
        RandomPcg|RandomPCG)
            SCHEME="RandomPcg"
            SCHEME_BASE="RandomPcg"
            SCHEME_ARGS="--sketching-scheme Minimizer --order RandomPcg"
            SCHEME_NEEDS_GM=0
            ;;
        GreedyMini)
            SCHEME="GreedyMini"
            SCHEME_BASE="GreedyMini_$(basename "${GM_FILE}" .gm)"
            SCHEME_ARGS="--sketching-scheme Minimizer --order GreedyMini --greedy-mini-file ${GM_FILE}"
            SCHEME_NEEDS_GM=1
            ;;
        ModRandom|ModRandomXor)
            SCHEME="ModRandom"
            SCHEME_BASE="ModRandom"
            SCHEME_ARGS="--sketching-scheme Mod --order RandomXor"
            SCHEME_NEEDS_GM=0
            ;;
        ModGreedyMini)
            SCHEME="ModGreedyMini"
            SCHEME_BASE="ModGreedyMini_$(basename "${MOD_GM_FILE}" .gm)"
            SCHEME_ARGS="--sketching-scheme Mod --order GreedyMini --greedy-mini-file ${MOD_GM_FILE}"
            SCHEME_NEEDS_GM=1
            ;;
        *)
            echo "Error: Unknown scheme '${scheme}'." >&2
            bench_scheme_help_line >&2
            return 1
            ;;
    esac
}

# Sets NEED_GM / NEED_MOD_GM from DEFAULT_SCHEMES.
bench_scheme_gm_needs() {
    NEED_GM=0
    NEED_MOD_GM=0
    local scheme
    for scheme in "${DEFAULT_SCHEMES[@]}"; do
        case "${scheme}" in
            GreedyMini) NEED_GM=1 ;;
            ModGreedyMini) NEED_MOD_GM=1 ;;
        esac
    done
}

# Build a .gm from its matching .txt (via the shared serializer) if missing.
bench_ensure_gm_file() {
    local gm_path="$1" txt_path tool_dir="${CODEBASE}/../common/minimizer-orders"
    [[ -n "${gm_path}" && -f "${gm_path}" ]] && return 0
    [[ -n "${gm_path}" ]] || return 1
    txt_path="${gm_path%.gm}.txt"
    [[ -f "${txt_path}" ]] || return 1
    echo "Building ${gm_path} from ${txt_path}" >&2
    make -C "${tool_dir}" >&2
    "${tool_dir}/gm_serializer" --canonical "${txt_path}" -o "${gm_path}" >&2
}

# Exit if selected schemes need missing .gm tables.
bench_require_gm_files() {
    bench_scheme_gm_needs
    bench_ensure_gm_file "${GM_FILE}"
    bench_ensure_gm_file "${MOD_GM_FILE}"
    if [[ "${NEED_GM}" -eq 1 && ! -f "${GM_FILE}" ]]; then
        echo "Error: GreedyMini requires -g / GM_FILE (.gm): ${GM_FILE}" >&2
        exit 1
    fi
    if [[ "${NEED_MOD_GM}" -eq 1 && ! -f "${MOD_GM_FILE}" ]]; then
        echo "Error: ModGreedyMini requires -M / MOD_GM_FILE (.gm): ${MOD_GM_FILE}" >&2
        exit 1
    fi
}
