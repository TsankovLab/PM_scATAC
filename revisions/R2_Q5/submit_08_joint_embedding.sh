#!/bin/bash
# Submit the joint CCA co-embedding (step 8) to LSF.
cd "$(dirname "$0")"
mkdir -p logs
bsub -P acc_Tsankov_Normal_Lung -q premium -n 4 \
     -R "rusage[mem=40000] span[hosts=1]" -W 12:00 \
     -J R2Q5_joint -o logs/08_joint.out -e logs/08_joint.err \
     "$PWD/run_step.sh 08_joint_embedding.R"
