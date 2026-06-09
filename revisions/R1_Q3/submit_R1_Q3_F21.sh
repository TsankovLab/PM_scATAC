#!/bin/bash
#BSUB -J R1Q3_F21
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 2
#BSUB -R "rusage[mem=16000]"
#BSUB -W 0:30
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/logs/R1Q3_F21_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/logs/R1Q3_F21_%J.err

module purge

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/R1_Q3_F21_SOX9_TGFB_correlation.R
