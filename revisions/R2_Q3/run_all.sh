#!/bin/bash
###############################################################################
# Reproduce every figure downstream of the epiAneufinder runs (steps 3, 4, 9, 10, 11).
#
# Steps 1 and 2 are the expensive ones and are NOT run here -- they only need to be
# run once, and out_5Mb/ already holds their output. To redo them from scratch:
#     ./run_step.sh 01_extract_fragments.R                       # ~10 min, ~11 GB
#     for S in P1 P4 P5 P8 P10 P11 P12 P14 P23; do
#       bsub -P acc_Tsankov_Normal_Lung -q premium -n 8 -R "rusage[mem=8000] span[hosts=1]" \
#            -W 6:00 -J epi_$S -o logs/epi_$S.out -e logs/epi_$S.err \
#            "$PWD/run_step.sh 02_run_epianeufinder.R $S"
#     done
#
# Step 3 is the memory hungry one (an 11331 x 11331 distance matrix for P23) -- run on a
# node with >=8 GB, not on a login node. Step 11 pulls the chromVAR z matrix out of the
# ArchR project on its first run (cached to chromvar_z_cache.rds).
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
run() {  # run <log name> <script> ; environment variables set by the caller are passed on
  echo "=================================================================="
  echo ">>> $1   $(date '+%H:%M:%S')"
  echo "=================================================================="
  ./run_step.sh "$2" 2>&1 | tee "logs/$1.log"
}
run 03_clone_calls                 03_clone_calls.R
run 04_figure_subclones_overview   04_figure_subclones_overview.R
run 09_P4_visium_chr8q             09_P4_visium_chr8q.R
MAL_CALL=label run 09_P4_visium_chr8q_label 09_P4_visium_chr8q.R
MAL_CALL=label SCORE=relative run 09_P4_visium_chr8q_label_relative 09_P4_visium_chr8q.R
run 10_P4_visium_chr8q_clustering  10_P4_visium_chr8q_clustering.R
MAL_CALL=label run 10_P4_visium_chr8q_clustering_label 10_P4_visium_chr8q_clustering.R
run 11_tf_subclone_wilcoxon        11_tf_subclone_wilcoxon.R
echo; echo "ALL DONE  $(date '+%H:%M:%S')"
