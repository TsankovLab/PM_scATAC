#!/bin/bash
###############################################################################
# Reproduce everything downstream of the epiAneufinder runs (steps 3-9).
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
# Steps 3-9 take ~7 min in total on a single interactive core. Step 3 is the memory
# hungry one (an 11331 x 11331 distance matrix for P23) -- run on a node with >=8 GB,
# not on a login node.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
for s in 03_clone_calls.R \
         04_figure_subclones_overview.R \
         05_figure_circular_tree.R \
         06_P4_chr8q_validation.R \
         07_chromvar_tf_variability.R \
         08_figure_tf_variability.R \
         09_P4_visium_chr8q.R ; do
  echo "=================================================================="
  echo ">>> $s   $(date '+%H:%M:%S')"
  echo "=================================================================="
  ./run_step.sh "$s" 2>&1 | tee "logs/${s%.R}.log"
done
echo; echo "ALL DONE  $(date '+%H:%M:%S')"
