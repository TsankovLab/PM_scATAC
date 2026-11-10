#!/bin/bash
#BSUB -J max_dist_SE
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=240000] span[hosts=1]"
#BSUB -W 24:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/max_dist_SE_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/max_dist_SE_%J.err
mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs
module purge
/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7/max_dist_sweep_SE.R
