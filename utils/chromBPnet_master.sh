#!/bin/bash
#BSUB -J chrBP_model
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=2
#BSUB -R v100
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

# ml anaconda3/2022.10
# ml cuda/11.7.0
# ml cudnn/8.9.5-11
# ml proxies
# ml java/11.0.2
# ml tensorrt/8.5.3.1

ml anaconda3/2020.11
source activate chrombpnet

chromBPdir=${1}
echo $chromBPdir
grefdir=${2}
echo $grefdir
repodir=${3}
echo $repodir
celltype=${4}
echo $celltype


#mkdir $chromBPdir
cd $chromBPdir

chmod +x ${repodir}/utils/chrBPnet_training_new.sh
chmod +x ${repodir}/utils/average_CNT_scores.py
chmod +x ${repodir}/utils/TFmodisco_counts.sh
chmod +x ${repodir}/utils/TFmodisco_profile.sh

job_ids=""
echo "run training model and contribution scores"
for fold_number in 0 1 2 3 4; do
    job_id=$(bsub -J ${celltype}_cBP \
         -P acc_Tsankov_Normal_Lung \
         -q gpu \
         -n 8 \
         -W 96:00 \
         -gpu num=2 \
         -R h100nvl \
         -R rusage[mem=32000] \
         -R span[hosts=1] \
         -o ${chromBPdir}/chormBPtraining_${celltype}_f${fold_number}.out \
         -e ${chromBPdir}/chormBPtraining_${celltype}_f${fold_number}.err \
         ${repodir}/utils/chrBPnet_training_new.sh "$chromBPdir" "$grefdir" "$celltype" "$fold_number" \
         | awk '{print $2}' | sed 's/<//;s/>//')
    
    # Append the job ID to the job_ids string
    if [ -z "$job_ids" ]; then
        job_ids="done(${job_id})"
    else
        job_ids="$job_ids && done(${job_id})"
    fi

done

chromBPct_dir=${chromBPdir}/${celltype}_model
echo $chromBPct_dir

# Wait for all jobs to complete
bsub -J wait_jobs -P acc_Tsankov_Normal_Lung -w "$job_ids"  -o wait.log -e wait.err /bin/bash -c "echo 'All jobs completed.'"

### Combine contribution scores 
echo "combine contribution scores"
ml anaconda3/2020.11
source deactivate
source activate h5py # activate another environment with hdf5plugin installed to read h5 files
# Explicitly set PATH
# Fix environment variables
export PATH=/sc/arion/work/giottb01/conda/envs/h5py/bin:$PATH
unset PYTHONPATH
export LD_LIBRARY_PATH=/sc/arion/work/giottb01/conda/envs/h5py/lib:$LD_LIBRARY_PATH

conda list | grep hdf5plugin
/sc/arion/work/giottb01/conda/envs/h5py/bin/python -c "import hdf5plugin; print('hdf5plugin is installed')"
/sc/arion/work/giottb01/conda/envs/h5py/bin/python $repodir/utils/average_CNT_scores.py $chromBPct_dir $celltype



bsub -J ${celltype}_TFmd_c \
    -P acc_Tsankov_Normal_Lung \
    -q premium \
    -n 8 \
    -W 96:00 \
    -R rusage[mem=64000] \
    -R span[hosts=1] \
    -o ${chromBPdir}/${celltype}_TFmodisco_counts.out \
    -e ${chromBPdir}/${celltype}_TFmodisco_counts.err \
    ${repodir}/utils/TFmodisco_counts.sh $chromBPdir $celltype

bsub -J ${celltype}_TFmd_p \
    -P acc_Tsankov_Normal_Lung \
    -q premium \
    -n 8 \
    -W 96:00 \
    -R rusage[mem=64000] \
    -R span[hosts=1] \
    -o ${chromBPdir}/${celltype}_TFmodisco_profiles.out \
    -e ${chromBPdir}/${celltype}_TFmodisco_profiles.err \
    ${repodir}/utils/TFmodisco_profile.sh $chromBPdir $celltype

# bsub -J ${celltype}_combS \
#      -P acc_Tsankov_Normal_Lung \
#      -o ${chromBPdir}/${celltype}_combine_scores.out \
#      -e ${chromBPdir}/${celltype}_combine_scores.err \
#      -R rusage[mem=64000] \
#      -w "$job_ids" \
#      /bin/bash -c "source activate chrombpnet && python $repodir/utils/average_CNT_scores.py $chromBPct_dir"
     

## Troubleshoot numpy (version installed should be 1.23.4)
# python
# import deepdish as dd
# import numpy as np
# d = {"test": np.array([0,0])}
# dd.io.save("test.h5", d, compression='blosc')


# !!!!!!ATTENTION!!!!!HOW TO USE THE NEW GPUS ADDED!!!!!ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!
# To submit jobs to NVLinked H100 GPU nodes, flag "-R h100nvl" is required.
# PLEASE ADD  #BSUB -R h100nvl to your LSF script or -R h100nvl to your LSF command line

# To submit jobs to L40S GPU nodes, flag "-R l40s" is required. 
# PLEASE ADD #BSUB -R l40s to your LSF script or -R l40s to your LSF command line. 