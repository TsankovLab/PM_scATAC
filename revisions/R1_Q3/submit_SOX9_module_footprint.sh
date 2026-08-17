#!/bin/bash
# submit_SOX9_module_footprint.sh
# Observed-Tn5 footprint of SOX9 motifs (chromVAR positions) in SOX9-high vs
# SOX9-low module-gated cells, via ArchR getFootprints/plotFootprints.
# Usage: bash submit_SOX9_module_footprint.sh
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR=/sc/arion/scratch/giottb01/meso_SOX9_module_footprint
mkdir -p "${WORKDIR}"

bsub \
  -J sox9_module_footprint \
  -P acc_Tsankov_Normal_Lung \
  -q premium \
  -n 8 -W 6:00 \
  -R "rusage[mem=16000]" -R "span[hosts=1]" \
  -o "${WORKDIR}/footprint.out" \
  -e "${WORKDIR}/footprint.err" \
  /bin/bash -c "
source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
conda activate meso_scatac
Rscript '${SCRIPTDIR}/R1_Q3_SOX9_module_footprint_archr.R'
"
echo "Submitted. Monitor: bjobs -u \$USER ; logs in ${WORKDIR}/footprint.{out,err}"
