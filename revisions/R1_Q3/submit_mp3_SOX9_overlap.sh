#!/bin/bash
#BSUB -J mp3_SOX9_ov
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 2
#BSUB -R "rusage[mem=48000]"
#BSUB -W 1:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/logs/mp3_SOX9_ov_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/logs/mp3_SOX9_ov_%J.err

module purge

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/R1_Q3_mp3_SOX9_overlap.R
