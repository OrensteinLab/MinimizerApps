# Shared scheme / DB-size CLI for run_benchmark.sh and benchmark_minimizers.sh.
# Source after benchmark_setup.sh with CODEBASE set:
#   source "${SCRIPT_DIR}/benchmark_interface.sh"
#   bench_iface_init

# shellcheck shell=bash

_bench_iface_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BENCH_SCHEME_CHOICES="xor | xor-pcg-hash | greedymini | mod-minimizer-xor | mod-minimizer-greedymini"

bench_iface_init() {
    local gm_dir="${CODEBASE}/../orders/GreedyMini"
    BASE_DIR="${KRAKEN2_BENCH_ROOT:-${_bench_iface_dir}/benchmark_data}"
    GM_FILE="${GM_FILE:-${gm_dir}/w5_k15_canonical.gm}"
    MOD_GM_FILE="${MOD_GM_FILE:-${gm_dir}/w5_k6_canonical.gm}"
    GM_PROJECTION="${GM_PROJECTION:-at0-cg1,at1-cg0}"   # comma-separated list: one run variant per projection
    DEFAULT_SCHEMES=("xor" "xor-pcg-hash" "greedymini" "mod-minimizer-xor" "mod-minimizer-greedymini")
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

# Common -o/-S/-g/-M/-P/-m/-z/-N/-E/-a option text (callers add their own extras).
bench_common_options_help() {
    cat <<EOF
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCH_ROOT or ./benchmark_data)
  -S SCHEME...  Run only these schemes (space-separated; repeatable):
                  ${BENCH_SCHEME_CHOICES}
                Examples: -S xor mod-minimizer-xor
                          -S "xor xor-pcg-hash greedymini"
                          -S xor -S xor-pcg-hash
  -g GM_FILE    binary w5_k15 table (.gm) for greedymini (or env GM_FILE)
  -M MOD_GM     binary w5_k6 table (.gm) for mod-minimizer-greedymini
                (or env MOD_GM_FILE)
  -P PROJ[,PROJ] greedymini DNA->binary projection(s) for both greedymini schemes:
                at0-cg1 (A,T -> 0) | at1-cg0 (flipped); comma-separated list =
                one run variant per projection (default: both, or env
                GM_PROJECTION). xor, xor-pcg-hash, mod-minimizer-xor run once.
  -m            NCBI ~10% subsample (db/standard_p10) → results/runs/standard_p10/
  -z            NCBI ~0.5% subsample (db/standard_p0.5) → results/runs/standard_p0.5/
  -N            NCBI 100% (db/standard) → results/runs/standard/ (same as omitting -m/-z)
  -E            Wood 2019 library (db/strain_exclusion) → results/runs/strain_exclusion/
                Mason + FDA-ARGOS queries (Tables S1 and S5); not with -m/-z/-N
  -a            With -E, skip FDA-ARGOS (Mason / Table S1 only)
                Default library: db/standard → results/runs/standard/ (unlabeled speed
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
            xor|xor-pcg-hash|greedymini|mod-minimizer-xor|mod-minimizer-greedymini) ;;
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
        g) GM_FILE="${OPTARG}" ;;
        M) MOD_GM_FILE="${OPTARG}" ;;
        P)
            GM_PROJECTION="${OPTARG}"
            ;;
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
# Shared FASTA: db/<id>/   Scheme indexes: results/runs/<id>/{scheme}/{seed}/

BENCH_LIB_IDS=(standard standard_p10 standard_p0.5 strain_exclusion)

# Mason + FDA-ARGOS accuracy mode (-E).
bench_is_paper_accuracy_mode() {
    [[ "${USE_STREX_LIB}" -eq 1 ]]
}

bench_canonical_shared_db() {
    echo "${BASE_DIR}/db/$1"
}

# <ROOT>/results/runs (KRAKEN2_RESULTS_DIR overrides <ROOT>/results)
bench_runs_dir() {
    echo "${KRAKEN2_RESULTS_DIR:-${BASE_DIR}/results}/runs"
}

bench_canonical_run_root() {
    echo "$(bench_runs_dir)/$1"
}

bench_shared_db_for() {
    case "$1" in
        standard|standard_p10|standard_p0.5|strain_exclusion) bench_canonical_shared_db "$1" ;;
        *)
            echo "Error: unknown library id '$1'" >&2
            return 1
            ;;
    esac
}

bench_run_root_for() {
    case "$1" in
        standard|standard_p10|standard_p0.5|strain_exclusion) bench_canonical_run_root "$1" ;;
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
            */db/strain_exclusion) USE_STREX_LIB=1 ;;
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

# For one variant ("<scheme>" or "<scheme>@<projection>", see bench_expand_variants) sets
#   SCHEME        the scheme (label in tables / manifest)
#   SCHEME_BASE   run directory (benchmark/common/scheme_names.scheme_dir_name:
#                 <scheme>[/k<gm_k>_<projection>_t<mod_t>])
#   SCHEME_GM_K / SCHEME_PROJECTION / SCHEME_GM_TABLE / SCHEME_MOD_T  id fields
#   SCHEME_ARGS, SCHEME_NEEDS_GM
bench_configure_scheme() {
    local scheme="${1%%@*}" gm_file="" proj=""
    [[ "$1" == *@* ]] && proj="${1#*@}"
    [[ -n "${proj}" ]] || proj="${GM_PROJECTIONS[0]:-${GM_PROJECTION%%,*}}"
    case "${scheme}" in
        xor)
            SCHEME_ARGS="--sketching minimizer --order xor"
            SCHEME_NEEDS_GM=0
            ;;
        xor-pcg-hash)
            SCHEME_ARGS="--sketching minimizer --order xor-pcg-hash"
            SCHEME_NEEDS_GM=0
            ;;
        greedymini)
            gm_file="${GM_FILE}"
            SCHEME_ARGS="--sketching minimizer --order greedymini --greedymini-file ${GM_FILE} --gm-projection ${proj}"
            SCHEME_NEEDS_GM=1
            ;;
        mod-minimizer-xor)
            SCHEME_ARGS="--sketching mod-minimizer --order xor"
            SCHEME_NEEDS_GM=0
            ;;
        mod-minimizer-greedymini)
            gm_file="${MOD_GM_FILE}"
            SCHEME_ARGS="--sketching mod-minimizer --order greedymini --greedymini-file ${MOD_GM_FILE} --gm-projection ${proj}"
            SCHEME_NEEDS_GM=1
            ;;
        *)
            echo "Error: Unknown scheme '${scheme}'." >&2
            bench_scheme_help_line >&2
            return 1
            ;;
    esac
    local ident
    ident="$("${PYTHON:-python3}" "${_bench_iface_dir}/bench_const.py" scheme-id \
        "${scheme}" "${gm_file}" "${gm_file:+${proj}}" "${BENCH_L:-31}" "${BENCH_K:-35}")" || return 1
    IFS='|' read -r SCHEME_BASE SCHEME_GM_K SCHEME_PROJECTION SCHEME_GM_TABLE SCHEME_MOD_T <<< "${ident}"
    SCHEME="${scheme}"
}

# GM_PROJECTION (comma list) -> GM_PROJECTIONS (validated, de-duplicated) and
# VARIANTS = DEFAULT_SCHEMES with every GreedyMini scheme expanded to
# "<scheme>@<projection>" per projection (other schemes appear once, bare).
bench_expand_variants() {
    local -a raw=() 
    local p s
    GM_PROJECTIONS=()
    IFS=',' read -r -a raw <<< "${GM_PROJECTION}"
    for p in "${raw[@]}"; do
        case "${p}" in
            at0-cg1|at1-cg0) ;;
            "") continue ;;
            *) echo "Error: -P / GM_PROJECTION must be a comma list of at0-cg1, at1-cg0 (got '${p}')." >&2; exit 1 ;;
        esac
        [[ " ${GM_PROJECTIONS[*]:-} " == *" ${p} "* ]] || GM_PROJECTIONS+=("${p}")
    done
    if [[ ${#GM_PROJECTIONS[@]} -eq 0 ]]; then
        echo "Error: empty -P / GM_PROJECTION." >&2
        exit 1
    fi
    VARIANTS=()
    for s in "${DEFAULT_SCHEMES[@]}"; do
        case "${s}" in
            greedymini|mod-minimizer-greedymini)
                for p in "${GM_PROJECTIONS[@]}"; do VARIANTS+=("${s}@${p}"); done ;;
            *) VARIANTS+=("${s}") ;;
        esac
    done
}

# Sets NEED_GM / NEED_MOD_GM from DEFAULT_SCHEMES.
bench_scheme_gm_needs() {
    NEED_GM=0
    NEED_MOD_GM=0
    local scheme
    for scheme in "${DEFAULT_SCHEMES[@]}"; do
        case "${scheme}" in
            greedymini) NEED_GM=1 ;;
            mod-minimizer-greedymini) NEED_MOD_GM=1 ;;
        esac
    done
}

# Build a .gm from its matching .txt (via the shared serializer) if missing.
bench_ensure_gm_file() {
    local gm_path="$1" txt_path tool_dir="${CODEBASE}/../common/minimizer-schemes"
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
        echo "Error: greedymini requires -g / GM_FILE (.gm): ${GM_FILE}" >&2
        exit 1
    fi
    if [[ "${NEED_MOD_GM}" -eq 1 && ! -f "${MOD_GM_FILE}" ]]; then
        echo "Error: mod-minimizer-greedymini requires -M / MOD_GM_FILE (.gm): ${MOD_GM_FILE}" >&2
        exit 1
    fi
}
