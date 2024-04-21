#!/bin/bash
#$ -M bgiotti@broadinstitute.org
#$ -N CNMF_consensus
#$ -cwd
#$ -q broad
#$ -l h_vmem=8g
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=24:00:00
#$ -e cnmf_consensus.err

source /broad/software/scripts/useuse
#use .anaconda3-5.3.1
use Anaconda

source activate /ahg/regevdata/projects/ICA_Lung/Bruno/conda/cnmf
projdir=${1}
cnmf_out=${2}
k_selection=${3}

cd $projdir$cnmf_out

cnmf consensus --output-dir ./ --name cnmf --components $k_selection --local-density-threshold 0.3 --show-clustering
