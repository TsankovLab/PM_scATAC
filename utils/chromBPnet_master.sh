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

ml anaconda3/2022.10
source activate chrombpnet

chromBPdir=${1}
grefdir=${2}
repodir=${3}
celltype=${4}
biasdir=${5}

chromBPct_dir=${chromBPdir}/${celltype}
mkdir -p $chromBPct_dir
cd $chromBPct_dir

chmod +x ${repodir}/utils/chromBPnet_training.sh
chmod +x ${repodir}/utils/chromBPnet_contribution_scores.sh
chmod +x ${repodir}/utils/chromBPnet_average_CNT_scores.py
chmod +x ${repodir}/utils/TFmodisco_counts.sh
chmod +x ${repodir}/utils/TFmodisco_profile.sh
chmod +x ${repodir}/utils/finemo_motif_calls.sh

echo "=== Run training model and contribution scores ==="

all_contrib_jobs=""

for fold_number in 0 1 2 3 4; do
    MODEL_H5=no_bias_model/fold_${fold_number}/models/chrombpnet_nobias.h5
    count_scores_file=no_bias_model/fold_${fold_number}/contribution_scores.count_scores.h5
    profile_scores_file=no_bias_model/fold_${fold_number}/contribution_scores.profile_scores.h5

    train_job_id=""
    contrib_job_id=""

    # --- Training job ---
    if [ ! -f "${MODEL_H5}" ]; then
        echo "[Fold ${fold_number}] Training job needed."
        rm -rf no_bias_model/fold_${fold_number}
        train_job_id=$(bsub -J ${celltype}_CBPtrain_f${fold_number} \
             -P acc_Tsankov_Normal_Lung \
             -q gpu \
             -n 1 \
             -W 24:00 \
             -gpu num=1 \
             -R a100 \
             -R rusage[mem=64000] \
             -R span[hosts=1] \
             -o ${chromBPdir}/chromBPtraining_${celltype}_f${fold_number}.out \
             -e ${chromBPdir}/chromBPtraining_${celltype}_f${fold_number}.err \
             ${repodir}/utils/chromBPnet_training.sh "$chromBPdir" "$grefdir" "$celltype" "$fold_number" "$biasdir" \
             | awk '{print $2}' | sed 's/<//;s/>//')
    else
        echo "[Fold ${fold_number}] Model already exists, skipping training."
        train_job_id=$(bsub -J dummy_train_${celltype}_f${fold_number} \
             -o logs/dummy_train_${celltype}_f${fold_number}.out \
             -e logs/dummy_train_${celltype}_f${fold_number}.err \
             /bin/bash -c "echo 'Model exists for ${celltype}, fold ${fold_number}'" \
             | awk '{print $2}' | sed 's/<//;s/>//')
    fi

    # --- Contribution job ---
    if [ ! -f "${count_scores_file}" ] || [ ! -f "${profile_scores_file}" ]; then
        echo "[Fold ${fold_number}] Contribution job needed."
        contrib_job_id=$(bsub -J ${celltype}_CBPcontrib_f${fold_number} \
             -P acc_Tsankov_Normal_Lung \
             -w "done(${train_job_id})" \
             -q gpu \
             -n 1 \
             -W 48:00 \
             -gpu num=1 \
             -R a100 \
             -R rusage[mem=64000] \
             -R span[hosts=1] \
             -o ${chromBPdir}/chromBPnet_contribution_scores_${celltype}_f${fold_number}.out \
             -e ${chromBPdir}/chromBPnet_contribution_scores_${celltype}_f${fold_number}.err \
             ${repodir}/utils/chromBPnet_contribution_scores.sh "$chromBPdir" "$grefdir" "$celltype" "$fold_number" "$biasdir" "$MODEL_H5" \
             | awk '{print $2}' | sed 's/<//;s/>//')
    else
        echo "[Fold ${fold_number}] Contribution already exists, skipping."
        contrib_job_id=$(bsub -J dummy_contrib_${celltype}_f${fold_number} \
             -o logs/dummy_contrib_${celltype}_f${fold_number}.out \
             -e logs/dummy_contrib_${celltype}_f${fold_number}.err \
             /bin/bash -c "echo 'Contribution exists for ${celltype}, fold ${fold_number}'" \
             | awk '{print $2}' | sed 's/<//;s/>//')
    fi

    # Collect for global dependency
    if [ -z "$all_contrib_jobs" ]; then
        all_contrib_jobs="done(${contrib_job_id})"
    else
        all_contrib_jobs="$all_contrib_jobs && done(${contrib_job_id})"
    fi
done


echo "=== Combine contribution scores ==="

source activate h5py
export PATH=/sc/arion/work/giottb01/conda/envs/h5py/bin:$PATH
unset PYTHONPATH
export LD_LIBRARY_PATH=/sc/arion/work/giottb01/conda/envs/h5py/lib:$LD_LIBRARY_PATH

# Combine HDF5 scores
/sc/arion/work/giottb01/conda/envs/h5py/bin/python $repodir/utils/chromBPnet_average_CNT_scores.py $chromBPct_dir $celltype

# Average bigwig files
wiggletools mean no_bias_model/fold_{0..4}/contribution_scores.counts_scores.bw \
    > no_bias_model/temp_contribution_counts_scores.wig
wigToBigWig no_bias_model/temp_contribution_counts_scores.wig ${grefdir}/hg38.chrom.sizes no_bias_model/${celltype}_averaged_contribution_scores_counts.bw
rm no_bias_model/temp_contribution_counts_scores.wig

wiggletools mean no_bias_model/fold_{0..4}/contribution_scores.profile_scores.bw \
    > no_bias_model/temp_contribution_profile_scores.wig
wigToBigWig no_bias_model/temp_contribution_profile_scores.wig ${grefdir}/hg38.chrom.sizes no_bias_model/${celltype}_averaged_contribution_scores_profile.bw
rm no_bias_model/temp_contribution_profile_scores.wig

# Cleanup fold-level contribution files
avg_contribution_file=no_bias_model/${celltype}_averaged_contributions_counts.h5
if [ -f "${avg_contribution_file}" ]; then
    echo "Average contribution file exists. Cleaning up fold-level contribution files..."
    rm -rf no_bias_model/fold_{0..4}/contribution*
fi

echo "=== Run TF-MoDISco on averaged contributions (parallel) ==="
TFmd_c_id=$(bsub -J ${celltype}_TFmd_c \
    -P acc_Tsankov_Normal_Lung \
    -q premium \
    -n 8 \
    -W 72:00 \
    -R rusage[mem=64000] \
    -R span[hosts=1] \
    -o ${chromBPdir}/${celltype}_TFmodisco_counts.out \
    -e ${chromBPdir}/${celltype}_TFmodisco_counts.err \
    ${repodir}/utils/TFmodisco_counts.sh $chromBPdir $celltype | awk '{print $2}' | sed 's/<//;s/>//')

TFmd_p_id=$(bsub -J ${celltype}_TFmd_p \
    -P acc_Tsankov_Normal_Lung \
    -q premium \
    -n 8 \
    -W 72:00 \
    -R rusage[mem=64000] \
    -R span[hosts=1] \
    -o ${chromBPdir}/${celltype}_TFmodisco_profiles.out \
    -e ${chromBPdir}/${celltype}_TFmodisco_profiles.err \
    ${repodir}/utils/TFmodisco_profile.sh $chromBPdir $celltype | awk '{print $2}' | sed 's/<//;s/>//')

echo "=== Submit finemo job (depends on TF-MoDISco) ==="
finemo_job_id=$(bsub -J ${celltype}_finemo \
    -P acc_Tsankov_Normal_Lung \
    -q gpu \
    -n 8 \
    -W 72:00 \
    -gpu num=2 \
    -R a100 \
    -R rusage[mem=32000] \
    -R span[hosts=1] \
    -o ${chromBPdir}/finemo_${celltype}.out \
    -e ${chromBPdir}/finemo_${celltype}.err \
    -w "done(${TFmd_c_id}) && done(${TFmd_p_id})" \
    ${repodir}/utils/finemo_motif_calls.sh "$chromBPdir" "$celltype" \
    | awk '{print $2}' | sed 's/<//;s/>//')
echo "Submitted finemo job with ID: $finemo_job_id"

# Wait for finemo to finish
echo "Waiting for finemo job..."
bwait -w "done(${finemo_job_id})"

echo "=== Run R script for finemo motif labels ==="
source activate meso_scatac
Rscript $repodir/utils/chromBPnet_finemo_motif_labels.R $chromBPdir $celltype
echo "R script execution completed."

# bsub -J ${celltype}_combS \
#      -P acc_Tsankov_Normal_Lung \
#      -o ${chromBPdir}/${celltype}_combine_scores.out \
#      -e ${chromBPdir}/${celltype}_combine_scores.err \
#      -R rusage[mem=64000] \
#      -w "$job_ids" \
#      /bin/bash -c "source activate chrombpnet && python $repodir/utils/average_CNT_scores.py $chromBPct_dir"

### Remove contribution scores h5 files for each fold in each cell types provided the average file is found ###



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