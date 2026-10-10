# Shared toggle-mask seeds (source after SCRIPT_DIR is set).
# Canonical list: benchmark/kraken2/bench_const.py DEFAULT_SEEDS
# shellcheck shell=bash

if [[ -z "${SCRIPT_DIR:-}" ]]; then
    echo "Error: bench_seeds.sh requires SCRIPT_DIR" >&2
    return 1 2>/dev/null || exit 1
fi

mapfile -t DEFAULT_SEEDS < <(python3 "${SCRIPT_DIR}/bench_const.py" seeds)
if [[ "${#DEFAULT_SEEDS[@]}" -lt 1 ]]; then
    echo "Error: bench_const.py seeds produced no toggles." >&2
    return 1 2>/dev/null || exit 1
fi
SEEDS=("${DEFAULT_SEEDS[@]}")
