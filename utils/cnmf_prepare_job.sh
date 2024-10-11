#!/bin/bash
#BSUB -J cnmf_prepare
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
k_list=${2}
nfeat=${3}
cnmf_out=${4}
cores=${5}

# Print variables
echo $projdir
echo $k_list
echo $nfeat
echo $cnmf_out
echo $cores

#genes_file=${3}
#cnmf_dir="cNMF_${k_list}_${nfeat}"

#k_list=' ' read -ra vector_array <<< "$k_list"
cd ${projdir}${cnmf_out}

# Run cNMF
#cnmf prepare --output-dir ./ --name cnmf -c ../counts_nmf_${nfeat}.txt -k $k_list --n-iter 100 --seed 14 --tpm ../norm_nmf_${nfeat}.txt --numgenes $nfeat --total-workers $cores #--genes $genes_file #
cnmf prepare --output-dir ./ --name cnmf -c ../counts_nmf_${nfeat}.txt -k $k_list --n-iter 100 --seed 14 --numgenes $nfeat --total-workers $cores #--genes $genes_file #




## Using UGER
# logdir = '.'
# uger_factorize_cmd = "qsub -cwd -b y -l h_vmem=2g,h_rt=3:00:00 -o . -e . -N cnmf -t 1-4 'python ./cnmf.py factorize --output-dir example_PBMC/cNMF --name pbmc_cNMF --worker-index $SGE_TASK_ID'"
# print('Factorize command to simultaneously run factorization over %d nodes of a compute cluster using UGER:\n%s' % (numworkers, uger_factorize_cmd))
# #!{cmd}
