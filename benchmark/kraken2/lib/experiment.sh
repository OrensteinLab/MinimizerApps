# Run a multi-config experiment from benchmark/kraken2/experiments/*.yaml.
# Mirrors benchmark/kmc/lib/experiment.sh's pattern, backed by the same
# shared benchmark/common engine.

# shellcheck shell=bash

bench_py_experiments() {
    "${PYTHON:-python3}" "${SCRIPT_DIR}/experiments.py" "$@"
}

reset_run_options() {
    SCHEME=""
    N_SEEDS=1
    SCHEME_THREADS=16
    PARALLEL_JOBS=1
    LIBRARY="standard"
}

library_flag() {
    case "${LIBRARY}" in
        standard_p10) echo "-m" ;;
        standard_p0.5) echo "-z" ;;
        strain_exclusion) echo "-E" ;;
        standard|*) echo "" ;;
    esac
}

run_experiment() {
    local spec="${EXPERIMENT:?experiment name or path required}"
    local -a env_lines
    mapfile -t env_lines < <(bench_py_experiments shell-env-all "${spec}")
    local n="${#env_lines[@]}" i flag

    echo "=== experiment ${spec} (${n} configs) ==="
    for ((i = 0; i < n; i++)); do
        reset_run_options
        eval "${env_lines[$i]}"
        flag="$(library_flag)"
        echo "=== [$((i + 1))/${n}] scheme=${SCHEME} n_seeds=${N_SEEDS} threads=${SCHEME_THREADS} jobs=${PARALLEL_JOBS} library=${LIBRARY} $(date) ==="
        "${SCRIPT_DIR}/run_benchmark.sh" ${flag:+"${flag}"} \
            -S "${SCHEME}" -n "${N_SEEDS}" \
            -j "${PARALLEL_JOBS}" -t "${SCHEME_THREADS}"
    done
}

list_experiments() {
    bench_py_experiments list
}
