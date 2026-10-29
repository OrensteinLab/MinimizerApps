#!/usr/bin/env bash
# Main-configuration kmc benchmark (k=28 p=9 w=20, gzip, bins 1024, threads 24, no -sm; -m12, and
# h_sapiens_3 -m32 as in datasets.csv) into a fresh results dir, each timed kmc invocation waiting
# for a quiet machine (benchmark/common/wait_low_load.sh), then export.
#
#   RESULTS_NAME=results experiments/full_run.sh
#
# Per dataset 30 runs, same set as the existing results (config dirs gzip_n1024_t24_m{12,32} and gm_proj_n1024_t24):
#   default x5; xor x5; xor-pcg-hash x5; greedymini lut_dna at0-cg1 x5   (one `-S all -n 5` call, interleaved)
#   greedymini lut_binary at0-cg1 x5; greedymini lut_binary at1-cg0 x5   (--run-tag gm_proj_n1024_t24)
# Re-running skips runs that already have a kmc.json. KMC temp dirs live in <results>/runs/<run>/tmp
# (data filesystem) and are removed after every run. Density caches (deterministic, not timings) are
# copied from the old results' csvs because export only reads them; set OLD_RESULTS to change the source.
set -uo pipefail
fail=0   # a failed step is reported and the run goes on; the exit status says so (full_run_all.sh)
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
if [ -z "${KMC3_BENCH_ROOT:-}" ] && [ -f .env ]; then
  KMC3_BENCH_ROOT="$(sed -n 's/^KMC3_BENCH_ROOT=//p' .env | tr -d "\"' ")"
fi
: "${KMC3_BENCH_ROOT:?KMC3_BENCH_ROOT is not set (benchmark/kmc/.env)}"
export KMC3_BENCH_ROOT
RES="$KMC3_BENCH_ROOT/${RESULTS_NAME:-results}"
OLD="${OLD_RESULTS:-$KMC3_BENCH_ROOT/results}"
export KMC3_RESULTS_DIR="$RES"
DATASETS="${DATASETS:-f_vesca g_gallus m_balbisiana h_sapiens_2 h_sapiens_3}"
mkdir -p "$RES/csvs"
echo "results: $RES  $(date -Is)"

for f in kmc3_density_expected.csv kmc3_density_expected_schemes.csv kmc3_density_particular.csv; do
  if [ -f "$OLD/csvs/$f" ] && [ ! -f "$RES/csvs/$f" ]; then cp "$OLD/csvs/$f" "$RES/csvs/$f"; fi
done

common=(-k 28 --input gzip -t 24 --bins 1024 --results-dir "$RES")
for ds in $DATASETS; do
  echo "=== $ds $(date -Is)"
  NO_EXPORT=1 ./kmc3_bench.sh run -d "$ds" "${common[@]}" -S all -n 5 \
    || { echo "run failed: $ds main (re-run resumes)" >&2; fail=1; }
  for proj in at0-cg1 at1-cg0; do
    NO_EXPORT=1 ./kmc3_bench.sh run -d "$ds" "${common[@]}" -S greedymini -n 5 \
      --gm-lut binary --projection "$proj" --run-tag gm_proj_n1024_t24 \
      || { echo "run failed: $ds $proj (re-run resumes)" >&2; fail=1; }
  done
done
./kmc3_bench.sh export --results-dir "$RES" || { echo "export failed" >&2; fail=1; }
echo "=== done $(date -Is)$([ "$fail" = 0 ] || echo "  (some step failed, see above)")"
exit "$fail"
