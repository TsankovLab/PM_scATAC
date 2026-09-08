#!/bin/bash
# Submit the cross-modality label transfer (step 2) to LSF.
cd "$(dirname "$0")"
mkdir -p logs
bsub -P acc_Tsankov_Normal_Lung -q premium -n 4 \
     -R "rusage[mem=32000] span[hosts=1]" -W 12:00 \
     -J R2Q5_transfer -o logs/02_transfer.out -e logs/02_transfer.err \
     "$PWD/run_step.sh 02_label_transfer.R"
