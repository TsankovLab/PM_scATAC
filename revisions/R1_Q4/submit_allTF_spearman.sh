#!/bin/bash
#BSUB -J allTF_spearman_R1Q4
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 4
#BSUB -R "rusage[mem=32000]"
#BSUB -W 4:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/allTF_spearman_R1Q4_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/allTF_spearman_R1Q4_%J.err

mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

export OPENBLAS_NUM_THREADS=1

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/plot_allTF_spearman_per_sample.R
