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

for fold_number in 1 2 3 4; do
    bsub -J ${celltype}_cBP \
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
         ${repodir}/utils/chrBPnet_training_new.sh "$chromBPdir" "$grefdir" "$celltype" "$fold_number"
done


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