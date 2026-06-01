#!/bin/bash
#BSUB -J plot_harmonized_motifs
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 4
#BSUB -W 4:00
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/plot_harmonized_motifs.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/plot_harmonized_motifs.err

ml anaconda3/latest
source activate meso_scatac

Rscript /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/plot_harmonized_motifs.R
