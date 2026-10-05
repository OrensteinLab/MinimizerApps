#!/usr/bin/env bash
# run_benchmark — minimizer-scheme benchmark (speed + accuracy).
#
# Prep:
#   ./benchmark/kraken2/download_data.sh      # data/strain-exclusion + data/fda-argos (Table S5)
#   ./benchmark/kraken2/prepare_cami_data.sh  # data/cami (CAMI II marine sample; optional)
#   ./benchmark/kraken2/prepare_mock_data.sh  # data/mock (ZymoBIOMICS D6300 run; optional)
#   ./benchmark/kraken2/download_library.sh   # db/standard + db/standard_p{10,0.5}
#   ./benchmark/kraken2/run_benchmark.sh -m   # → runs/standard_p10
#   ./benchmark/kraken2/run_benchmark.sh -E   # → runs/strain_exclusion
#
# Layout:
#   db/{standard,strain_exclusion,standard_p10,standard_p0.5}/
#   data/{strain-exclusion,fda-argos,cami,mock}/
#   runs/{standard,strain_exclusion,standard_p10,standard_p0.5}/
#     {scheme}/{seed}/{db,results}

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${SCRIPT_DIR}/benchmark_setup.sh"
# shellcheck source=benchmark_interface.sh
source "${SCRIPT_DIR}/benchmark_interface.sh"

CODEBASE="$(cd "${SCRIPT_DIR}/../../kraken2" && pwd)"

# GNU time is required for -v (wall + max RSS). Bash builtin `time` cannot.
TIME_BIN="${TIME_BIN:-}"
if [[ -z "${TIME_BIN}" ]]; then
    if [[ -x /usr/bin/time ]]; then
        TIME_BIN=/usr/bin/time
    elif command -v time >/dev/null 2>&1 && time -v true >/dev/null 2>&1; then
        TIME_BIN="$(command -v time)"
    fi
fi
if [[ -z "${TIME_BIN}" || ! -x "${TIME_BIN}" ]]; then
    echo "Error: GNU time not found (need /usr/bin/time with -v)." >&2
    echo "  Install: sudo apt-get install -y time   (or re-run ./benchmark/kraken2/setup.sh)" >&2
    exit 1
fi

for script in k2 kraken2 kraken2-build; do
    path="${CODEBASE}/scripts/${script}"
    if [[ ! -f "${path}" ]]; then
        echo "Error: missing ${path}" >&2
        exit 1
    fi
    echo "Found ${script} in ${path}"
done

bench_iface_init
SCHEME=""
SKIP_CAMI="${SKIP_CAMI:-0}"
SKIP_MOCK="${SKIP_MOCK:-0}"
EXTRA_SETS="${EXTRA_SETS:-0}"
SKIP_DENSITY="${SKIP_DENSITY:-0}"
DENSITY_JOBS="${DENSITY_JOBS:-0}"   # benchmark_minimizers -j (0 = all cores)
COMPARE_ONLY="${COMPARE_ONLY:-0}"

N_SEEDS="${N_SEEDS:-5}"
N_CLASSIFY="${N_CLASSIFY:-50}"
PARALLEL_JOBS="${PARALLEL_JOBS:-1}"
SCHEME_THREADS="${SCHEME_THREADS:-16}"
# shellcheck source=bench_seeds.sh
source "${SCRIPT_DIR}/bench_seeds.sh"

show_help() {
    cat <<EOF
Usage: $0 [options]

Build + classify over N seeds, then score.

Paper pairing (Wood et al. 2019): one strain-exclusion library (-E) and two
query sets — Mason held-out genomes (Table S1 sens/PPV/F1) and FDA-ARGOS
(Table S5 genus concordance). NCBI -m/-z/-N never classify those sets.
Two extra query sets are used when prepared (prepare_cami_data.sh,
prepare_mock_data.sh): CAMI II marine reads (data/cami; library-aware genus /
species sens/PPV/F1) and the ZymoBIOMICS D6300 mock (data/mock; profile
precision/recall, L1, Bray-Curtis). With -E they are classified automatically;
with -m/-z/-N only when -Q is given.

Classify order: seed → dataset → cold + timed repeats → shuffled schemes.
Seed = toggle mask (Random* order; GreedyMini Γ / Mod tb).

Options:
$(bench_common_options_help)
  -n N          Number of seeds (default: ${N_SEEDS})
  -r N          Timed classification repeats per seed after 1 cold run
                (default: ${N_CLASSIFY}; cold run is discarded from stats)
  -j N          Max concurrent scheme builds per seed (default: ${PARALLEL_JOBS});
                classify always runs schemes sequentially
  -t N          Threads for k2 build (default: ${SCHEME_THREADS});
                classify always uses 1 thread
  -C            Skip the CAMI query set (data/cami)
  -K            Skip the mock-community query set (data/mock)
  -Q            With -m/-z/-N also classify data/cami and data/mock
                (default for them: speed set only; -E always includes them)
  -D            Skip the density microbench (benchmark_minimizers.sh). By default
                it runs after the builds for every -S scheme x seed 1..N that has
                no row yet in runs/<lib>/tables/minimizer_density.csv (density
                factor on random DNA and on the library FASTAs, same toggle as
                the seed's databases); DENSITY_JOBS sets its -j (default 0 = all
                cores)
  -c            Compare only (existing results + manifest; also plots).
                Reuses -m/-z/-N/-E when given; otherwise auto-detects a single
                tree under runs/ that has a manifest
  -h, --help    Show help
EOF
    exit 0
}

while getopts ":o:d:S:g:M:n:r:j:t:mzaNEcATCKQDh-" opt; do
    case "${opt}" in
        S)
            bench_handle_S_opt "$@"
            SCHEME="${DEFAULT_SCHEMES[*]}"
            ;;
        n) N_SEEDS="${OPTARG}" ;;
        r) N_CLASSIFY="${OPTARG}" ;;
        j) PARALLEL_JOBS="${OPTARG}" ;;
        t) SCHEME_THREADS="${OPTARG}" ;;
        c) COMPARE_ONLY=1 ;;
        C) SKIP_CAMI=1 ;;
        K) SKIP_MOCK=1 ;;
        Q) EXTRA_SETS=1 ;;
        D) SKIP_DENSITY=1 ;;
        A) : ;; # ARGOS on by default in paper mode (compat with old run_eval.sh -A)
        T) : ;; # seed list already includes the paper toggle (old run_eval.sh -T)
        h) show_help ;;
        -)
            case "${OPTARG}" in
                help) show_help ;;
                *) echo "Unknown option --${OPTARG}" >&2; exit 1 ;;
            esac
            ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
        *)
            if ! bench_handle_common_opt "${opt}"; then
                echo "Invalid option: -${opt}" >&2
                exit 1
            fi
            ;;
    esac
done
shift $((OPTIND - 1)) || true

bench_validate_scheme_filter
bench_require_gm_files

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && show_help

if [[ -n "${SEEDS_STR:-}" ]]; then
    read -r -a SEEDS <<< "${SEEDS_STR}"
fi

DATA_DIR="${DATA_DIR:-${BASE_DIR}/data}"
TABLES="${SCRIPT_DIR}/tables"
bench_resolve_shared_db

# Bind summary/manifest/figures paths to this library's run root.
set_run_root() {
    RUN_ROOT="$(bench_run_root)"
    SUMMARY_FILE="${RUN_ROOT}/benchmark_comparison.txt"
    RUN_MANIFEST="${RUN_ROOT}/.benchmark_run_manifest"
    FIGURES_DIR="${RUN_ROOT}/figures"
    ROUND_TABLES="${RUN_ROOT}/classify_round_tables.txt"
}

set_run_root

# For -c without -m/-z/-N/-E: pick the only run root that has a manifest.
_manifest_at() {
    [[ -f "$1/.benchmark_run_manifest" ]]
}

resolve_compare_run_root() {
    if [[ -f "${RUN_MANIFEST}" ]]; then
        return 0
    fi
    if bench_is_paper_accuracy_mode; then
        echo "Error: -c requires existing ${RUN_MANIFEST}" >&2
        echo "  (-E selected $(bench_run_root_for strain_exclusion))" >&2
        exit 1
    fi
    if [[ "${SIZE_EXPLICIT}" -eq 1 ]]; then
        echo "Error: -c requires existing ${RUN_MANIFEST}" >&2
        echo "  (flag selected ${LIB_ID}; check -m/-z/-N)" >&2
        exit 1
    fi
    local -a found=()
    local id root
    for id in "${BENCH_LIB_IDS[@]}"; do
        root="$(bench_run_root_for "${id}")"
        if _manifest_at "${root}"; then
            found+=("${id}")
        fi
    done
    if [[ ${#found[@]} -eq 1 ]]; then
        echo "[*] -c auto-detected run root: ${found[0]}"
        case "${found[0]}" in
            strain_exclusion) USE_STREX_LIB=1 ;;
            standard_p10) USE_SAMPLE=1 ;;
            standard_p0.5) USE_SMOKE_SAMPLE=1 ;;
            standard) ;;
        esac
        SHARED_DB=""
        bench_resolve_shared_db
        set_run_root
        return 0
    fi
    if [[ ${#found[@]} -eq 0 ]]; then
        echo "Error: -c found no .benchmark_run_manifest under ${BASE_DIR}/runs/{$(IFS=,; echo "${BENCH_LIB_IDS[*]}")}" >&2
        exit 1
    fi
    echo "Error: -c found multiple run roots with manifests: ${found[*]}" >&2
    echo "  Pass -m (standard_p10), -z (standard_p0.5), -N (standard), or -E (strain_exclusion)." >&2
    exit 1
}
if ! [[ "${N_SEEDS}" =~ ^[0-9]+$ ]] || [[ "${N_SEEDS}" -lt 1 ]]; then
    echo "Error: -n / N_SEEDS must be a positive integer." >&2
    exit 1
fi
if ! [[ "${N_CLASSIFY}" =~ ^[0-9]+$ ]] || [[ "${N_CLASSIFY}" -lt 1 ]]; then
    echo "Error: -r / N_CLASSIFY must be a positive integer." >&2
    exit 1
fi
if ! [[ "${PARALLEL_JOBS}" =~ ^[0-9]+$ ]] || [[ "${PARALLEL_JOBS}" -lt 1 ]]; then
    echo "Error: -j / PARALLEL_JOBS must be a positive integer." >&2
    exit 1
fi
N_SCHEMES="${#DEFAULT_SCHEMES[@]}"
if [[ "${PARALLEL_JOBS}" -gt "${N_SCHEMES}" ]]; then
    PARALLEL_JOBS="${N_SCHEMES}"
fi
if ! [[ "${SCHEME_THREADS}" =~ ^[0-9]+$ ]] || [[ "${SCHEME_THREADS}" -lt 1 ]]; then
    echo "Error: -t / SCHEME_THREADS must be a positive integer." >&2
    exit 1
fi
if [[ "${#SEEDS[@]}" -lt "${N_SEEDS}" ]]; then
    echo "Error: need at least ${N_SEEDS} seeds; have ${#SEEDS[@]}." >&2
    exit 1
fi
SEEDS=("${SEEDS[@]:0:${N_SEEDS}}")

configure_scheme() {
    bench_configure_scheme "$1"
}

# Seed = toggle mask for all schemes (Random order; GM Γ / Mod tb XOR).
configure_seed_paths() {
    local seed_idx="$1"
    local seed="${SEEDS[$((seed_idx - 1))]}"
    SEED="${seed}"
    SEED_TAG="${seed_idx}"
    SCHEME_CONTEXT="${RUN_LAYOUT}/${SCHEME_BASE}/${SEED_TAG}"
    DB_DIR="${RUN_ROOT}/${SCHEME_BASE}/${SEED_TAG}/db"
    OUTPUT_DIR="${RUN_ROOT}/${SCHEME_BASE}/${SEED_TAG}/results"
    TOGGLE_ARGS=(--toggle-mask "${seed}")
}

# Paper (Wood et al. Methods): Mason + FDA-ARGOS are both classified against
# the strain-exclusion library only. NCBI -m/-z/-N use an unlabeled speed set.
# Prefer data/{strain-exclusion,fda-argos,speed}/R1.fastq; accept older names.
pe_pair_in_dir() {
    local dir="$1"
    local r1 r2
    [[ -d "${dir}" ]] || return 1
    for r1 in \
        "${dir}/R1.fastq" \
        "${dir}/accuracy_1.fq" \
        "${dir}/argos_R1.fastq" \
        "${dir}/benchmark_R1.fastq"
    do
        [[ -s "${r1}" ]] || continue
        r2="${r1/R1.fastq/R2.fastq}"
        r2="${r2/_1.fq/_2.fq}"
        r2="${r2/_R1.fastq/_R2.fastq}"
        if [[ -s "${r2}" ]]; then
            printf '%s\t%s\n' "${r1}" "${r2}"
            return 0
        fi
    done
    return 1
}

add_dataset() {
    local id="$1" r1="$2" r2="$3" truth="${4:-}"
    DATASET_IDS+=("${id}")
    DATASET_R1["${id}"]="${r1}"
    DATASET_R2["${id}"]="${r2}"
    DATASET_TRUTH["${id}"]="${truth}"
}

# CAMI (data/cami, truth.tsv) and mock community (data/mock, composition.csv):
# added when prepared and not skipped with -C / -K.
add_extra_datasets() {
    local pair r1 r2
    if [[ "${SKIP_CAMI}" -eq 0 ]] && pair="$(pe_pair_in_dir "${DATA_DIR}/cami")"; then
        IFS=$'\t' read -r r1 r2 <<< "${pair}"
        add_dataset "cami" "${r1}" "${r2}" "${DATA_DIR}/cami/truth.tsv"
    fi
    if [[ "${SKIP_MOCK}" -eq 0 ]] && pair="$(pe_pair_in_dir "${DATA_DIR}/mock")"; then
        IFS=$'\t' read -r r1 r2 <<< "${pair}"
        add_dataset "mock" "${r1}" "${r2}" "${DATA_DIR}/mock/composition.csv"
    fi
}

setup_datasets() {
    DATASET_IDS=()
    declare -gA DATASET_R1=()
    declare -gA DATASET_R2=()
    declare -gA DATASET_TRUTH=()
    local pair r1 r2 se_dir arg_dir
    se_dir="${DATA_DIR}/strain-exclusion"
    arg_dir="${DATA_DIR}/fda-argos"

    if bench_is_paper_accuracy_mode; then
        if pair="$(pe_pair_in_dir "${se_dir}")"; then
            IFS=$'\t' read -r r1 r2 <<< "${pair}"
            add_dataset "strain-exclusion" "${r1}" "${r2}" "${se_dir}/truth.tsv"
        elif pair="$(pe_pair_in_dir "${DATA_DIR}/simulated")"; then
            IFS=$'\t' read -r r1 r2 <<< "${pair}"
            add_dataset "strain-exclusion" "${r1}" "${r2}" "${DATA_DIR}/simulated/truth.tsv"
        elif pair="$(pe_pair_in_dir "${BASE_DIR}/evaluation/reads")"; then
            IFS=$'\t' read -r r1 r2 <<< "${pair}"
            add_dataset "strain-exclusion" "${r1}" "${r2}" "${BASE_DIR}/evaluation/reads/truth.tsv"
        fi
        if [[ "${SKIP_ARGOS}" -eq 0 ]]; then
            if pair="$(pe_pair_in_dir "${arg_dir}")"; then
                IFS=$'\t' read -r r1 r2 <<< "${pair}"
                add_dataset "fda-argos" "${r1}" "${r2}" "${arg_dir}/truth.tsv"
            elif pair="$(pe_pair_in_dir "${BASE_DIR}/evaluation/reads/argos")"; then
                IFS=$'\t' read -r r1 r2 <<< "${pair}"
                add_dataset "fda-argos" "${r1}" "${r2}" \
                    "${BASE_DIR}/evaluation/reads/argos/truth.tsv"
            fi
        fi
        add_extra_datasets
        if [[ "${#DATASET_IDS[@]}" -eq 0 ]]; then
            echo "Error: -E needs Mason, FDA-ARGOS, CAMI and/or mock FASTQs." >&2
            echo "  ./benchmark/kraken2/download_data.sh -o ${BASE_DIR}" >&2
            echo "  ./benchmark/kraken2/prepare_strex_data.sh -o ${BASE_DIR}" >&2
            echo "  ./benchmark/kraken2/prepare_cami_data.sh -o ${BASE_DIR}" >&2
            echo "  ./benchmark/kraken2/prepare_mock_data.sh -o ${BASE_DIR}" >&2
            exit 1
        fi
        return 0
    fi

    # NCBI -m/-z/-N: unlabeled speed corpus only (not Table S1 / S5).
    if pair="$(pe_pair_in_dir "${DATA_DIR}/speed")"; then
        IFS=$'\t' read -r r1 r2 <<< "${pair}"
        add_dataset "speed" "${r1}" "${r2}" ""
    elif pair="$(pe_pair_in_dir "${DATA_DIR}")"; then
        IFS=$'\t' read -r r1 r2 <<< "${pair}"
        add_dataset "speed" "${r1}" "${r2}" ""
    fi
    if [[ "${EXTRA_SETS}" -eq 1 ]]; then
        add_extra_datasets
    fi
    if [[ "${#DATASET_IDS[@]}" -eq 0 ]]; then
        echo "Error: NCBI -m/-z/-N do not classify strain-exclusion or FDA-ARGOS." >&2
        echo "  Paper accuracy (Tables S1 and S5): ./benchmark/kraken2/run_benchmark.sh -E" >&2
        echo "  Speed on this library: put unlabeled PE FASTQs in ${DATA_DIR}/speed/" >&2
        exit 1
    fi
}

perf_stem() {
    echo "classification_performance_$1"
}

nextmin_stem() {
    echo "nextminimizer_timing_$1"
}

class_out_name() {
    echo "${1}_classification.txt"
}

report_out_name() {
    echo "${1}_report.txt"
}

initialize_shared_environment() {
    if [[ ! -d "${SHARED_DB}/taxonomy" || ! -d "${SHARED_DB}/library" ]]; then
        echo "Error: missing ${SHARED_DB}/{taxonomy,library}." >&2
        echo "Run ./benchmark/kraken2/download_library.sh first (use -m or -z for subsampled libraries)." >&2
        echo "Strain-exclusion library: ./benchmark/kraken2/prepare_strex_data.sh then -E" >&2
        exit 1
    fi
    setup_datasets
    mkdir -p "${RUN_ROOT}"
    echo "========================================================="
    echo " Kraken 2 scheme benchmark"
    echo "========================================================="
    echo "  Library:          ${RUN_LAYOUT}  (${SHARED_DB})"
    echo "  Base directory:   ${BASE_DIR}"
    echo "  Run root:         ${RUN_ROOT}"
    echo "  Datasets:         ${DATASET_IDS[*]}"
    echo "  Layout:           ${RUN_LAYOUT}/{scheme}/{seed}/{db,results}"
    echo "  Build order:      seed × scheme (${PARALLEL_JOBS} concurrent builds/seed)"
    echo "  Classify order:   seed → dataset → rep i → shuffled schemes (sequential)"
    echo "  Schemes:          ${DEFAULT_SCHEMES[*]}"
    echo "  Seeds:            ${N_SEEDS} (build ${SCHEME_THREADS} threads)"
    echo "  Classify:         1 thread; ${N_CLASSIFY} timed + 1 cold (discarded) per seed"
    echo "  Toggle seeds:     ${SEEDS[*]}"
    echo "  GM (Minimizer):   ${GM_FILE}"
    echo "  GM (Mod):         ${MOD_GM_FILE}"
    echo "========================================================="
}

build_database() {
    if [[ ! -f "${DB_DIR}/hash.k2d" ]]; then
        ln -sfn "${SHARED_DB}/taxonomy" "${DB_DIR}/taxonomy"
        ln -sfn "${SHARED_DB}/library" "${DB_DIR}/library"
        echo "[*] Building ${SCHEME_CONTEXT} (toggle ${SEED})..."
        local build_output="${OUTPUT_DIR}/build_performance.txt"
        if ! "${TIME_BIN}" -v "${K2:-k2}" build \
            --db "${DB_DIR}" \
            --threads "${SCHEME_THREADS}" \
            ${SCHEME_ARGS} \
            "${TOGGLE_ARGS[@]}" \
            2> "${build_output}"; then
            echo "Error: k2 build failed for ${SCHEME_CONTEXT}" >&2
            return 1
        fi
    else
        echo "[*] Existing hash.k2d for ${SCHEME_CONTEXT}; skipping build."
    fi

    if [[ ! -f "${DB_DIR}/hash.k2d" || ! -f "${DB_DIR}/opts.k2d" ]]; then
        echo "Error: database build incomplete for ${SCHEME_CONTEXT}" >&2
        return 1
    fi

    if [[ "${SCHEME_NEEDS_GM}" -eq 1 ]] && [[ ! -f "${DB_DIR}/gm.k2d" ]]; then
        echo "Error: gm.k2d missing after ${SCHEME} build in ${DB_DIR}" >&2
        return 1
    fi
}

# Parse GNU time -v file → wall seconds and max RSS in GB (stdout: "wall rss").
# Label contains colons ("h:mm:ss or m:ss"), so strip by fixed prefix — not ".*:".
parse_time_v() {
    local f="$1"
    local elapsed rss_kb wall="?" rss_gb="?"
    elapsed=$(grep -E 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):' "${f}" 2>/dev/null \
        | sed -E 's/^[[:space:]]*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):[[:space:]]*//' \
        | tr -d '[:space:]') || true
    rss_kb=$(grep -E 'Maximum resident set size \(kbytes\):' "${f}" 2>/dev/null \
        | sed -E 's/.*Maximum resident set size \(kbytes\):[[:space:]]*//' \
        | tr -d '[:space:]') || true
    if [[ -n "${elapsed}" ]]; then
        # Formats: h:mm:ss(.xx) | m:ss(.xx) | ss(.xx)
        wall=$(awk -v t="${elapsed}" 'BEGIN {
            n = split(t, a, ":")
            if (n == 3)      printf "%.3f", a[1]*3600 + a[2]*60 + a[3]
            else if (n == 2) printf "%.3f", a[1]*60 + a[2]
            else if (n == 1) printf "%.3f", a[1]+0
            else             print "?"
        }')
    fi
    if [[ -n "${rss_kb}" && "${rss_kb}" =~ ^[0-9]+$ ]]; then
        rss_gb=$(awk -v k="${rss_kb}" 'BEGIN { printf "%.3f", k / 1048576 }')
    fi
    printf '%s %s' "${wall}" "${rss_gb}"
}

# Parse "N sequences classified (XX.XX%)" from a classify stderr/time file.
parse_classified_pct() {
    local f="$1"
    local pct
    pct=$(grep -E 'sequences classified \([0-9.]+%\)' "${f}" 2>/dev/null \
        | head -1 \
        | sed -E 's/.*\(([0-9.]+)%\).*/\1/') || true
    if [[ -n "${pct}" ]]; then
        printf '%s' "${pct}"
    else
        printf '?'
    fi
}

classify_once() {
    local time_file="$1"
    local nextmin_file="$2"
    local class_out="$3"
    local report_out="$4"
    local r1="$5"
    local r2="$6"
    local -a classify_extra=()
    # Paper / 2.0.8-beta parity: native kraken2 had no hit-group threshold; k2
    # defaults to --minimum-hit-groups 2 (post-2.0.8). Disable for -E scoring.
    if bench_is_paper_accuracy_mode; then
        classify_extra+=(--minimum-hit-groups 0)
    fi
    env KRAKEN2_MINIMIZER_PROF="${nextmin_file}" \
        "${TIME_BIN}" -v "${K2:-k2}" classify \
        --db "${DB_DIR}" \
        --threads 1 \
        --paired \
        "${classify_extra[@]}" \
        --output "${class_out}" \
        --report "${report_out}" \
        "${r1}" "${r2}" \
        2> "${time_file}"
}

# Print + append a wall/RSS table for one (seed, rep) across schemes.
# Args: seed_idx rep_label  then triples: scheme wall_s rss_gb ...
log_classify_round_table() {
    local seed_idx="$1"
    local ds="$2"
    local rep_label="$3"
    shift 3
    local table_file="${ROUND_TABLES}"
    {
        echo ""
        echo "---- seed#${seed_idx}  dataset=${ds}  rep=${rep_label} ----"
        printf '%-24s %12s %10s %10s\n' "Scheme" "Wall(s)" "RSS(GB)" "Class%"
        printf '%-24s %12s %10s %10s\n' "------------------------" "------------" "----------" "----------"
        while (( $# >= 4 )); do
            printf '%-24s %12s %10s %10s\n' "$1" "$2" "$3" "$4"
            shift 4
        done
    } | tee -a "${table_file}"
}

clear_scheme_classify_outputs() {
    rm -f "${OUTPUT_DIR}/classification_performance.txt"
    rm -f "${OUTPUT_DIR}/classification_performance_cold.txt"
    rm -f "${OUTPUT_DIR}"/classification_performance_*.txt
    rm -f "${OUTPUT_DIR}/nextminimizer_timing.txt"
    rm -f "${OUTPUT_DIR}/nextminimizer_timing_cold.txt"
    rm -f "${OUTPUT_DIR}"/nextminimizer_timing_*.txt
}

build_seed_scheme() {
    local scheme="$1"
    local seed_idx="$2"
    configure_scheme "${scheme}" || return 1
    configure_seed_paths "${seed_idx}"

    mkdir -p "${DB_DIR}" "${OUTPUT_DIR}"

    echo ""
    echo "---------------------------------------------------------"
    echo " Build: ${SCHEME}  seed#${seed_idx}  toggle=${SEED}"
    echo " DB:    ${DB_DIR}"
    echo " Out:   ${OUTPUT_DIR}"
    echo " Threads: ${SCHEME_THREADS}"
    echo "---------------------------------------------------------"

    build_database || return 1
    clear_scheme_classify_outputs
    echo "${N_CLASSIFY}" > "${OUTPUT_DIR}/n_classify_repeats.txt"

    {
        flock -x 200
        printf '%s\t%s\t%s\t%s\t%d\t%d\n' \
            "${SCHEME}" "${SCHEME_CONTEXT}" "${OUTPUT_DIR}" "${SEED}" \
            "${seed_idx}" "${N_CLASSIFY}" \
            >> "${RUN_MANIFEST}"
    } 200>"${RUN_MANIFEST}.lock"
}

build_seed_schemes() {
    local seed_idx="$1"
    local -a pids=()
    local scheme fail=0 pid
    for scheme in "${DEFAULT_SCHEMES[@]}"; do
        while (( ${#pids[@]} >= PARALLEL_JOBS )); do
            wait -n || fail=1
            local -a still=()
            for pid in "${pids[@]}"; do
                if kill -0 "${pid}" 2>/dev/null; then
                    still+=("${pid}")
                fi
            done
            pids=("${still[@]}")
        done
        build_seed_scheme "${scheme}" "${seed_idx}" &
        pids+=("$!")
    done
    for pid in "${pids[@]}"; do
        wait "${pid}" || fail=1
    done
    return "${fail}"
}

# Fisher–Yates shuffle of named array (in place).
shuffle_array() {
    local -n _arr="$1"
    local i j tmp
    for ((i = ${#_arr[@]} - 1; i > 0; i--)); do
        j=$((RANDOM % (i + 1)))
        tmp="${_arr[i]}"
        _arr[i]="${_arr[j]}"
        _arr[j]="${tmp}"
    done
}

# One classify pass over all schemes for (seed_idx, rep_label).
# Scheme order is shuffled each round to dilute first-/last-position bias.
# rep_label: "cold" or "01".."NN". Writes per-scheme timing files.
classify_seed_round() {
    local seed_idx="$1"
    local rep_label="$2"
    local ds scheme time_file nextmin_file class_out report_out pstem nstem
    local wall rss class_pct r1 r2
    local -a row=()
    local -a schemes

    for ds in "${DATASET_IDS[@]}"; do
        schemes=("${DEFAULT_SCHEMES[@]}")
        shuffle_array schemes
        row=()
        echo ""
        echo "[*] Classify seed#${seed_idx} dataset=${ds} rep=${rep_label}  order: ${schemes[*]}"
        r1="${DATASET_R1[${ds}]}"
        r2="${DATASET_R2[${ds}]}"
        pstem="$(perf_stem "${ds}")"
        nstem="$(nextmin_stem "${ds}")"
        for scheme in "${schemes[@]}"; do
            configure_scheme "${scheme}" || return 1
            configure_seed_paths "${seed_idx}"
            class_out="${OUTPUT_DIR}/$(class_out_name "${ds}")"
            report_out="${OUTPUT_DIR}/$(report_out_name "${ds}")"
            if [[ "${rep_label}" == "cold" ]]; then
                time_file="${OUTPUT_DIR}/${pstem}_cold.txt"
                nextmin_file="${OUTPUT_DIR}/${nstem}_cold.txt"
            else
                time_file="${OUTPUT_DIR}/${pstem}_${rep_label}.txt"
                nextmin_file="${OUTPUT_DIR}/${nstem}_${rep_label}.txt"
            fi
            echo "    ${SCHEME} (${ds}) ..."
            if ! classify_once "${time_file}" "${nextmin_file}" \
                    "${class_out}" "${report_out}" "${r1}" "${r2}"; then
                echo "Error: k2 classify failed for ${SCHEME_CONTEXT} ${ds} (rep ${rep_label})" >&2
                return 1
            fi
            read -r wall rss <<< "$(parse_time_v "${time_file}")"
            class_pct="$(parse_classified_pct "${time_file}")"
            row+=("${SCHEME}" "${wall}" "${rss}" "${class_pct}")
        done
        log_classify_round_table "${seed_idx}" "${ds}" "${rep_label}" "${row[@]}"
    done
}

score_dataset_output() {
    local ds="$1"
    local nodes="${SHARED_DB}/taxonomy/nodes.dmp"
    local kraken="${OUTPUT_DIR}/$(class_out_name "${ds}")"
    local truth="${DATASET_TRUTH[${ds}]:-}"
    if [[ -z "${truth}" || ! -s "${truth}" ]]; then
        return 0
    fi
    if [[ "${ds}" == "strain-exclusion" ]]; then
        python3 "${SCRIPT_DIR}/score_kraken.py" \
            --nodes "${nodes}" \
            --truth "${truth}" \
            --kraken "${kraken}" \
            --scheme "${SCHEME}" \
            --dataset strain-exclusion \
            --out-tsv "${OUTPUT_DIR}/strain_exclusion_metrics.tsv" \
            --expected-s1 "${TABLES}/table_s1.csv"
    elif [[ "${ds}" == "fda-argos" ]]; then
        python3 "${SCRIPT_DIR}/score_argos.py" \
            --nodes "${nodes}" \
            --truth "${truth}" \
            --kraken "${kraken}" \
            --scheme "${SCHEME}" \
            --out-tsv "${OUTPUT_DIR}/fda_argos_metrics.tsv" \
            --expected "${TABLES}/table_s5_kraken2_expected.csv"
    elif [[ "${ds}" == "cami" ]]; then
        # Library-aware: reads whose true genus/species is absent from this
        # library are counted as excluded, not as misses.
        python3 "${SCRIPT_DIR}/score_cami.py" \
            --nodes "${nodes}" \
            --truth "${truth}" \
            --library-dir "${SHARED_DB}" \
            --kraken "${kraken}" \
            --scheme "${SCHEME}" \
            --dataset cami \
            --out-tsv "${OUTPUT_DIR}/cami_metrics.tsv"
    elif [[ "${ds}" == "mock" ]]; then
        # Profile scoring from the kraken2 report; "truth" is composition.csv.
        python3 "${SCRIPT_DIR}/score_mock.py" \
            --nodes "${nodes}" \
            --report "${OUTPUT_DIR}/$(report_out_name "${ds}")" \
            --composition "${truth}" \
            --scheme "${SCHEME}" \
            --dataset mock \
            --out-tsv "${OUTPUT_DIR}/mock_metrics.tsv"
    fi
}

finalize_seed_scheme_outputs() {
    local seed_idx="$1"
    local scheme ds pad pstem nstem class_out report_out
    printf -v pad "%02d" "${N_CLASSIFY}"
    for scheme in "${DEFAULT_SCHEMES[@]}"; do
        configure_scheme "${scheme}" || return 1
        configure_seed_paths "${seed_idx}"
        for ds in "${DATASET_IDS[@]}"; do
            class_out="${OUTPUT_DIR}/$(class_out_name "${ds}")"
            report_out="${OUTPUT_DIR}/$(report_out_name "${ds}")"
            if [[ ! -f "${class_out}" || ! -f "${report_out}" ]]; then
                echo "Error: classification outputs missing for ${SCHEME_CONTEXT} (${ds})" >&2
                return 1
            fi
            pstem="$(perf_stem "${ds}")"
            nstem="$(nextmin_stem "${ds}")"
            ln -sfn "${pstem}_${pad}.txt" "${OUTPUT_DIR}/${pstem}.txt"
            if [[ -f "${OUTPUT_DIR}/${nstem}_${pad}.txt" ]]; then
                ln -sfn "${nstem}_${pad}.txt" "${OUTPUT_DIR}/${nstem}.txt"
            fi
            score_dataset_output "${ds}" || return 1
        done
    done
}

run_seed_with_schemes() {
    local seed_idx="$1"
    local seed="${SEEDS[$((seed_idx - 1))]}"
    local i pad

    echo ""
    echo "#########################################################"
    echo " Seed #${seed_idx} (${seed}) — schemes: ${DEFAULT_SCHEMES[*]}"
    echo " Build: ${PARALLEL_JOBS} concurrent (${SCHEME_THREADS} threads)"
    echo " Classify: seed → dataset → i → shuffled schemes (1 thread; ${N_CLASSIFY} timed + 1 cold)"
    echo "#########################################################"

    build_seed_schemes "${seed_idx}" || return 1

    # Cold round (excluded from compare stats), then timed interleaved rounds.
    classify_seed_round "${seed_idx}" "cold" || return 1
    for ((i = 1; i <= N_CLASSIFY; i++)); do
        printf -v pad "%02d" "${i}"
        classify_seed_round "${seed_idx}" "${pad}" || return 1
    done
    finalize_seed_scheme_outputs "${seed_idx}" || return 1
}

run_all_scheme_benchmarks() {
    mkdir -p "${RUN_ROOT}"
    # Several invocations may share one run root (e.g. -S Random GreedyMini first,
    # later -S RandomPcg ModRandom ModGreedyMini): only the manifest lines of the
    # schemes of THIS invocation are replaced; other schemes' lines (and the
    # round tables they wrote) stay, so compare reads all orders of the library.
    local keep_other=0
    if [[ -s "${RUN_MANIFEST}" ]]; then
        local tmp_manifest="${RUN_MANIFEST}.tmp.$$" sch
        local -a names=()
        for sch in "${DEFAULT_SCHEMES[@]}"; do
            configure_scheme "${sch}" || exit 1
            names+=("${SCHEME}")
        done
        awk -F'\t' -v names="$(IFS=,; echo "${names[*]}")" \
            'BEGIN { n = split(names, a, ","); for (i = 1; i <= n; i++) drop[a[i]] = 1 }
             !($1 in drop)' "${RUN_MANIFEST}" > "${tmp_manifest}"
        mv -f "${tmp_manifest}" "${RUN_MANIFEST}"
        if [[ -s "${RUN_MANIFEST}" ]]; then keep_other=1; fi
    else
        : > "${RUN_MANIFEST}"
    fi
    if [[ "${keep_other}" -eq 0 ]]; then
        : > "${ROUND_TABLES}"
    fi
    local idx fail=0
    for ((idx = 1; idx <= N_SEEDS; idx++)); do
        run_seed_with_schemes "${idx}" || fail=1
    done
    rm -f "${RUN_MANIFEST}.lock"
    if [[ "${fail}" -ne 0 ]]; then
        echo "Error: one or more seed/scheme benchmarks failed." >&2
        exit 1
    fi
}

compare_results() {
    local compare_py="${SCRIPT_DIR}/compare_benchmark.py"
    if [[ ! -f "${compare_py}" ]]; then
        echo "Error: missing ${compare_py}" >&2
        exit 1
    fi
    if [[ ! -f "${RUN_MANIFEST}" ]]; then
        echo "Error: missing run manifest ${RUN_MANIFEST}" >&2
        exit 1
    fi
    echo ""
    echo "[*] Statistical comparison across seeds (reference: Random)..."
    python3 "${compare_py}" \
        --base-dir "${RUN_ROOT}" \
        --manifest "${RUN_MANIFEST}" \
        --out "${SUMMARY_FILE}"
}

# Density microbench (expected = random DNA, particular = the library FASTAs
# SHARED_DB) per seed with the seed's toggle; fills only the (scheme, seed)
# rows missing from ${RUN_ROOT}/tables/minimizer_density.csv, which
# compare_benchmark.py reads. Failures warn and never abort the benchmark.
run_density() {
    if [[ "${SKIP_DENSITY}" -eq 1 ]]; then
        echo "[*] Density microbench skipped (-D)."
        return 0
    fi
    local script="${SCRIPT_DIR}/benchmark_minimizers.sh"
    local bin="${CODEBASE}/src/benchmark_minimizers"
    if [[ ! -x "${bin}" ]]; then
        echo "Warning: missing ${bin}; density microbench skipped (make -C kraken2/src benchmark_minimizers)." >&2
        return 0
    fi
    echo ""
    echo "[*] Density microbench (${LIB_ID}; library ${SHARED_DB})..."
    if ! SHARED_DB="${SHARED_DB}" OUT_FILE="${RUN_ROOT}/tables/minimizer_density.csv" \
        bash "${script}" -o "${BASE_DIR}" -S "${DEFAULT_SCHEMES[@]}" \
        -n "${N_SEEDS}" -j "${DENSITY_JOBS}" \
        -g "${GM_FILE}" -M "${MOD_GM_FILE}" \
        > "${RUN_ROOT}/density_microbench.log" 2>&1; then
        echo "Warning: density microbench failed; see ${RUN_ROOT}/density_microbench.log" >&2
        return 0
    fi
    grep '^\[\*\] seed' "${RUN_ROOT}/density_microbench.log" || true
}

plot_histograms() {
    local plot_py="${SCRIPT_DIR}/plot_benchmark_histograms.py"
    if [[ ! -f "${plot_py}" ]]; then
        echo "Warning: missing ${plot_py}; skipping histograms." >&2
        return 0
    fi
    echo ""
    echo "[*] Plotting classify / NextMinimizer histograms..."
    if ! python3 "${plot_py}" \
        --base-dir "${RUN_ROOT}" \
        --manifest "${RUN_MANIFEST}" \
        --out-dir "${FIGURES_DIR}"; then
        echo "Warning: histogram plotting failed (matplotlib missing?)." >&2
        return 0
    fi
}

if [[ "${COMPARE_ONLY}" -eq 1 ]]; then
    resolve_compare_run_root
    run_density
    compare_results
    plot_histograms
    echo ""
    echo "Comparison complete."
    echo "  Summary:  ${SUMMARY_FILE}"
    echo "  Tables:   ${SUMMARY_FILE%.txt}.xlsx  and  $(dirname "${SUMMARY_FILE}")/tables/"
    echo "  Figures:  ${FIGURES_DIR}/"
    exit 0
fi

initialize_shared_environment
run_all_scheme_benchmarks
run_density
compare_results
plot_histograms

echo ""
echo "Benchmark complete."
echo "  Summary:  ${SUMMARY_FILE}"
echo "  Tables:   ${SUMMARY_FILE%.txt}.xlsx  and  $(dirname "${SUMMARY_FILE}")/tables/"
echo "  Rounds:   ${ROUND_TABLES}"
echo "  Figures:  ${FIGURES_DIR}/"
