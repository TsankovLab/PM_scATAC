#!/bin/bash
#BSUB -J SOX9_allHDMA
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 2
#BSUB -R "rusage[mem=72000]"
#BSUB -W 2:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/logs/SOX9_allHDMA_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/logs/SOX9_allHDMA_%J.err

module purge

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/R1_Q3_SOX9_vs_all_HDMA_patterns.R
