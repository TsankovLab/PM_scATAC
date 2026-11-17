#!/bin/bash
# Submit the joint CCA co-embedding (step 8) to LSF.
source "${MESO_CONFIG:-$(for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$PWD"; do while [ ! -f "$d/config.sh" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; [ -f "$d/config.sh" ] && { echo "$d/config.sh"; break; }; done)}"   # data / results / reference roots: see config.sh at the repo root
cd "$(dirname "$0")"
mkdir -p logs
bsub -P ${LSF_ACCOUNT} -q premium -n 4 \
     -R "rusage[mem=40000] span[hosts=1]" -W 12:00 \
     -J R2Q5_joint -o logs/08_joint.out -e logs/08_joint.err \
     "$PWD/run_step.sh 08_joint_embedding.R"
