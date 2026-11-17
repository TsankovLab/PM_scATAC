#!/bin/bash
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
set -euo pipefail
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD=${SCRATCH_DIR}/meso_R2Q14; mkdir -p "$WD"
bsub -J q14_hist -P ${LSF_ACCOUNT} -q premium -n 2 -W 1:00 \
  -R "rusage[mem=32000]" -R "span[hosts=1]" -o "$WD/hist.out" -e "$WD/hist.err" \
  /bin/bash -c "export OPENBLAS_NUM_THREADS=2; source ${CONDA_INIT}; conda activate meso_scatac; Rscript '$SD/R2_Q14_histology_disentangle.R'"
echo "WD=$WD"
