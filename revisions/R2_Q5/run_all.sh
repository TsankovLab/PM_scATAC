#!/bin/bash
###############################################################################
# Reproduce the R2_Q5 figures (joint scATAC + scRNA embedding).
#
# Step 1 reads the ArchR project and the Seurat object and needs a few GB.
# Step 8 is the expensive one (CCA co-embedding, ~1-2 h, 40 GB): run it on a compute
# node, e.g. with ./submit_08_joint_embedding.sh, then run step 9 (seconds).
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
./run_step.sh 01_annotations.R 2>&1 | tee logs/01_annotations.log
./run_step.sh 08_joint_embedding.R 2>&1 | tee logs/08_joint_embedding.log   # ~1-2 h; or submit_08_joint_embedding.sh
./run_step.sh 09_figure_joint.R 2>&1 | tee logs/09_figure_joint.log
echo; echo "ALL DONE  $(date '+%H:%M:%S')"
