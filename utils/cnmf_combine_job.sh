#!/bin/bash
#BSUB -J cnmf_combine
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 1
#BSUB -R rusage[mem=2000]
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

#env >> ~/job_env
#set >> ~/job_env

# Run cNMF
# cnmf prepare --output-dir ./ --name $cnmf_out -c counts_nmf_${nfeat}.txt -k $k_list --n-iter 100 --seed 14 --tpm norm_nmf_${nfeat}.txt --numgenes $nfeat --total-workers 4 #--genes $genes_file #
#echo $SGE_TASK_ID
cd ${projdir}${cnmf_out}

# cnmf factorize --output-dir ./ --name $cnmf_out --worker-index $SGE_TASK_ID
# echo 'done!'
cnmf combine --output-dir ./ --name cnmf
cnmf k_selection_plot --output-dir ./ --name cnmf

# cnmf consensus --output-dir ./ --name $cnmf_out --components $k_selection --local-density-threshold 2 --show-clustering




## Using UGER
# logdir = '.'
# uger_factorize_cmd = "qsub -cwd -b y -l h_vmem=2g,h_rt=3:00:00 -o . -e . -N cnmf -t 1-4 'python ./cnmf.py factorize --output-dir example_PBMC/cNMF --name pbmc_cNMF --worker-index $SGE_TASK_ID'"
# print('Factorize command to simultaneously run factorization over %d nodes of a compute cluster using UGER:\n%s' % (numworkers, uger_factorize_cmd))
#!{cmd}



