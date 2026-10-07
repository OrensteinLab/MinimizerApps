#!/usr/bin/env bash
# One-time conda env setup, matching kraken2's own CONDA_ENV pattern.
# `pbsim2` (HMM error models) and PBSIM 1.0.3 (built below; the Winnowmap
# paper's simulator) are the long-read simulators lib/pbsim_reads.py drives;
# `samtools` indexes genomes (faidx).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
conda env create -f "$HERE/environment.yml" || conda env update -f "$HERE/environment.yml"

# PBSIM 1.0.3 (Ono et al. 2013), the simulator the Winnowmap paper used for its
# D1/D2, built from the official release into the env's opt/pbsim1 (bioconda's
# `pbsim` package is the same tarball + the same one-line compile fix).
ENV_DIR="$(conda env list | awk '$1=="minimap2_bench"{print $NF}')"
PB1="$ENV_DIR/opt/pbsim1"
if [ ! -x "$PB1/bin/pbsim" ]; then
  tmp="$(mktemp -d)"
  curl -sSL -o "$tmp/pbsim-1.0.3.tar.gz" \
    https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/pbsim/pbsim-1.0.3.tar.gz
  echo "c070394d42b6e4cdae0bb2aaf4856782  $tmp/pbsim-1.0.3.tar.gz" | md5sum -c -
  tar xzf "$tmp/pbsim-1.0.3.tar.gz" -C "$tmp"
  sed -i "s/argv\[optind\] == '\\\\0'/argv[optind] == nullptr/" "$tmp/pbsim-1.0.3/src/pbsim.cpp"
  (cd "$tmp/pbsim-1.0.3" && ./configure --prefix="$PB1" && make && make install)
  mkdir -p "$PB1/data" && cp "$tmp/pbsim-1.0.3/data/model_qc_clr" "$tmp/pbsim-1.0.3/data/model_qc_ccs" "$PB1/data/"
  rm -rf "$tmp"
fi
echo "Activate with: conda activate minimap2_bench"
echo "Then point PYTHON at its python3 in lib/constants.sh (or export PYTHON=...)."
