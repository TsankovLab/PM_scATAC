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

projdir=${1}
cnmf_out=${2}
repodir=${3}
nfeat=${4}
k_list=${5}
cores=${6}
cnmf_name=${7}

# Run cNMF
cd ${projdir}/${cnmf_out}

echo "prepare inputs..."
cnmf prepare --output-dir ./ --name cnmf -c ../counts_nmf_${nfeat}_${cnmf_name}.txt -k $k_list --n-iter 100 --seed 14 --numgenes $nfeat --total-workers $cores #--genes $genes_file #

echo "submit array job for cNMF factorization..."
chmod +x $repodir/cnmf_factorization_parallel.sh
job_id=$(bsub -P acc_Tsankov_Normal_Lung -J "cnmf_factorization[1-$cores]" -R rusage[mem=64000] -W 5:00 $repodir/cnmf_factorization_parallel.sh $projdir $cnmf_out $cores | awk '{print $2}' | sed 's/<//g' | sed 's/>//g')

### While loop to wait until factorization job array is completed
while true; do
  job_status=$(bjobs $job_id 2>&1)
  if [[ $job_status == *"DONE"* ]] || [[ $job_status == *"EXIT"* ]]; then
    echo "Job $job_id finished."
    break
  else
    echo "Waiting for job $job_id to finish..."
    sleep 10  # Wait for 10 seconds before checking again
  fi
done

echo "combine K iterations..."
cnmf combine --output-dir ./ --name cnmf
cnmf k_selection_plot --output-dir ./ --name cnmf

cnmf consensus --output-dir ./ --name cnmf --components $k_list --local-density-threshold 0.3 --show-clustering	

echo "cNMF jobs completed!"

## Using UGER
# logdir = '.'
# uger_factorize_cmd = "qsub -cwd -b y -l h_vmem=2g,h_rt=3:00:00 -o . -e . -N cnmf -t 1-4 'python ./cnmf.py factorize --output-dir example_PBMC/cNMF --name pbmc_cNMF --worker-index $SGE_TASK_ID'"
# print('Factorize command to simultaneously run factorization over %d nodes of a compute cluster using UGER:\n%s' % (numworkers, uger_factorize_cmd))
#!{cmd}



