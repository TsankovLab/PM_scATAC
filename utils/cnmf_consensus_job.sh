#!/bin/bash
#BSUB -J cnmf_consensus
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 1
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -W 3:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml anaconda3/2022.10

source activate cnmf

#source /broad/software/scripts/useuse
#use .anaconda3-5.3.1
#use Anaconda

#source activate /ahg/regevdata/projects/ICA_Lung/Bruno/conda/cnmf

projdir=${1}
cnmf_out=${2}
k_selection=${3}

cd $projdir$cnmf_out

cnmf consensus --output-dir ./ --name cnmf --components $k_selection --local-density-threshold 0.3 --show-clustering
