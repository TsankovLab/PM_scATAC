#!/bin/bash
# Submits ChromBPNet training for the purity-gated SOX9+ / SOX9- cell sets
# (module-score gated: SOX9pos module>=0.20 n=335, SOX9neg module<=-0.10 n=631).
#
# Step 1: prep job — subsets scatac_ArchR_SOX9_P23 to the gated cells, adds the
#         SOX9_module_gate column, saves to scatac_ArchR_SOX9_module_P23/
# Step 2: launch job (dependent on prep) — calls run_chrombpnet.sh, which submits
#         peak-calling then per-celltype ChromBPNet master jobs.
#
# Mirrors the previous iteration (submit_SOX9_P23_chrombpnet.sh) but uses the
# module-gated project/column from the adapted config.sh.
#
# Usage: bash submit_SOX9_module_chrombpnet.sh

set -euo pipefail

PIPELINEDIR="/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chrombpnet_hdma"
source "${PIPELINEDIR}/config.sh"

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_SCRIPT="${SCRIPTDIR}/R1_Q3_SOX9_module_prep_archr.R"

echo "=== SOX9 module-gated ChromBPNet submission ==="
echo "  WORKDIR        : ${WORKDIR}"
echo "  ARCHR_PROJECT  : ${ARCHR_PROJECT}"
echo "  CELLTYPE_COLUMN: ${CELLTYPE_COLUMN}"
echo "  CELLTYPES      : ${CELLTYPES[*]}"

mkdir -p "${WORKDIR}"

# Step 1: prep — subset ArchR project and assign SOX9_module_gate labels
PREP_JOB=$(bsub \
    -J cbp_SOX9_module_prep \
    -P "${HPC_PROJECT}" \
    -q premium \
    -n 4 -W 4:00 \
    -R "rusage[mem=32000]" -R "span[hosts=1]" \
    -o "${WORKDIR}/prep_archr.out" \
    -e "${WORKDIR}/prep_archr.err" \
    /bin/bash -c "
source '${CONDA_INIT}'
conda activate '${CONDA_ENV_R}'
Rscript '${PREP_SCRIPT}'
" | awk '{print $2}' | tr -d '<>')

echo "  Prep job submitted: ${PREP_JOB}"

# Step 2: launch — runs run_chrombpnet.sh (peak-calling + master jobs) after prep
bsub \
    -J cbp_SOX9_module_launch \
    -P "${HPC_PROJECT}" \
    -q premium \
    -n 1 -W 1:00 \
    -R "rusage[mem=8000]" -R "span[hosts=1]" \
    -w "done(${PREP_JOB})" \
    -o "${WORKDIR}/launch_chrombpnet.out" \
    -e "${WORKDIR}/launch_chrombpnet.err" \
    /bin/bash -c "bash '${PIPELINEDIR}/run_chrombpnet.sh'"

echo "  Launch job submitted (depends on ${PREP_JOB})"
echo "=== Done. Monitor with: bjobs -u \$USER ==="
