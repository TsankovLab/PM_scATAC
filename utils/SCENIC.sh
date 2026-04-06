#!/bin/bash
#BSUB -J scenic_multi_cpu       # job name
#BSUB -n 128                    # total CPU cores across nodes
#BSUB -R "rusage[mem=32GB]"     # memory per core
#BSUB -W 12:00                  # walltime
#BSUB -q normal
#BSUB -o scenic_multi_cpu.%J.out
#BSUB -e scenic_multi_cpu.%J.err
#BSUB -R "span[ptile=16]"       # 16 cores per node

# -------------------------------
# Load modules / activate conda
# -------------------------------

ml anaconda3/2022.10
source activate /sc/arion/work/giottb01/conda/envs/pyscenic

# -------------------------------
# Import variables
# -------------------------------
projdir_SC=${1}
motifs_tss=${2}
motifs_weights=${3}
TFs=${4}
expr_mat=${5}
scrna_pipeline_dir=${6}
vg=${7}
motif_window=${8}

# -------------------------------
# Start Dask scheduler (only on first node)
# -------------------------------
if [ $LSB_JOBINDEX -eq 0 ]; then
    echo "Starting Dask scheduler..."
    dask-scheduler --scheduler-file dask_scheduler.json &
    sleep 10
fi

# -------------------------------
# Start Dask workers on all nodes
# -------------------------------
echo "Starting Dask workers..."
# Each worker uses 1 thread to avoid Numba/TBB issues, memory limit per worker 16GB
dask-worker --scheduler-file dask_scheduler.json \
            --nthreads 1 \
            --memory-limit 16GB \
            --nworkers 8 &  # 8 workers per node (adjust based on cores)

sleep 10  # give workers time to connect

# -------------------------------
# Run SCENIC / GRNBoost2 Python script
# -------------------------------
/sc/arion/work/giottb01/conda/envs/pyscenic/bin/python "${scrna_pipeline_dir}/SCENIC.py" \
    $projdir_SC \
    $motifs_tss \
    $motifs_weights \
    $TFs \
    $expr_mat \
    $vg \
    $motif_window \
    dask_scheduler.json

# -------------------------------
# Optional cleanup
# -------------------------------
pkill -f dask-worker
pkill -f dask-scheduler


#$ -M bgiotti@broadinstitute.org
#$ -N SCENIC
#$ -cwd
#$ -q broad
#$ -l h_vmem=8g
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=24:00:00
#$ -e pySCENIC.err\n
#$ -m bea


# projdir_SC=${1}
# motifs_tss=${2}
# motifs_weights=${3}
# TFs=${4}
# expr_mat=${5}
# scrna_pipeline_dir=${6}
# vg=${7}
# motif_window=${8}

# need to add the full path of the conda environment to call the correct python version
#/sc/arion/work/giottb01/conda/envs/pyscenic/bin/python "${scrna_pipeline_dir}/SCENIC.py" $projdir_SC $motifs_tss $motifs_weights $TFs $expr_mat $vg $motif_window