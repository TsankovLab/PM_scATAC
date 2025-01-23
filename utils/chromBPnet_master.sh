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
biasdir=${5}
echo $biasdir

#mkdir $chromBPdir
chromBPct_dir=${chromBPdir}/${celltype}
echo $chromBPct_dir

mkdir $chromBPct_dir
cd $chromBPct_dir

chmod +x ${repodir}/utils/chrBPnet_training.sh
chmod +x ${repodir}/utils/chromBPnet_average_CNT_scores.py
chmod +x ${repodir}/utils/TFmodisco_counts.sh
chmod +x ${repodir}/utils/TFmodisco_profile.sh
chmod +x ${repodir}/utils/finemo_motif_calls.sh

### Create background regions file - Go in the output chromBPnet folder ####
# Remove regions overlapping black listed regions
bedtools slop -i ${grefdir}/blacklist.bed.gz -g ${grefdir}/hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a ../MACS2_${celltype}/${celltype}_peaks_capped.narrowPeak -b temp.bed  > ${celltype}_peakset_all_no_blacklist.bed
wc -l ${celltype}_peakset_all_no_blacklist.bed # # Make sure number of peaks is not more than 250K

# Generate training validation and test chromosome sets
head -n 23  ${grefdir}/hg38.chrom.sizes >  hg38.chrom.subset.sizes

# Generate background regions
for fold_number in 0 1 2 3 4; do
    
    negatives_file=no_bias_model/output_negatives_f${fold_number}_negatives.bed
    if [ ! -f "${negatives_file}" ]; then
    echo "negatives file not found. Identifying background peaks..."
    #rm -r output_auxiliary
    chrombpnet prep nonpeaks \
        -g ${grefdir}/genome_references/hg38.genome.fa \
        -p ${celltype}_peakset_all_no_blacklist.bed \
        -c ${grefdir}/hg38.chrom.sizes \
        -fl ${grefdir}/folds/fold_${fold_number}.json \
        -br ${grefdir}/blacklist.bed.gz \
        -o no_bias_model/output_negatives_f${fold_number}
    fi
done


job_ids=""
echo "Run training model and contribution scores"
for fold_number in 0 1 2 3 4; do
    job_id=$(bsub -J ${celltype}_cBP \
         -P acc_Tsankov_Normal_Lung \
         -q gpu \
         -n 8 \
         -W 72:00 \
         -gpu num=2 \
         -R h100nvl \
         -R rusage[mem=32000] \
         -R span[hosts=1] \
         -o ${chromBPdir}/chormBPtraining_${celltype}_f${fold_number}.out \
         -e ${chromBPdir}/chormBPtraining_${celltype}_f${fold_number}.err \
         ${repodir}/utils/chrBPnet_training.sh "$chromBPdir" "$grefdir" "$celltype" "$fold_number" "$biasdir" \
         | awk '{print $2}' | sed 's/<//;s/>//')
    
    # Append the job ID to the job_ids string
    if [ -z "$job_ids" ]; then
        job_ids="done(${job_id})"
    else
        job_ids="$job_ids && done(${job_id})"
    fi
done


# Wait for all jobs to complete
wait_job_id=$(bsub -J wait_jobs \
    -P acc_Tsankov_Normal_Lung \
    -w "$job_ids" \
    -o wait.log \
    -e wait.err \
    /bin/bash -c "echo 'All jobs completed.'" | awk '{print $2}' | sed 's/<//;s/>//')

# Actively wait for the wait_jobs to complete
bwait -w "done(${wait_job_id})"




### Combine contribution scores 
echo "combine contribution scores"
#ml anaconda3/2020.11
#source deactivate
source activate h5py # activate another environment with hdf5plugin installed to read h5 files

# Explicitly set PATH to detect hdf5plugin
# Fix environment variables
export PATH=/sc/arion/work/giottb01/conda/envs/h5py/bin:$PATH
unset PYTHONPATH
export LD_LIBRARY_PATH=/sc/arion/work/giottb01/conda/envs/h5py/lib:$LD_LIBRARY_PATH

conda list | grep hdf5plugin
/sc/arion/work/giottb01/conda/envs/h5py/bin/python -c "import hdf5plugin; print('hdf5plugin is installed')"
/sc/arion/work/giottb01/conda/envs/h5py/bin/python $repodir/utils/chromBPnet_average_CNT_scores.py $chromBPct_dir $celltype

echo "Take average of bigwig files counts"
wiggletools mean no_bias_model/fold_0/contribution_scores.counts_scores.bw \
no_bias_model/fold_1/contribution_scores.counts_scores.bw \
no_bias_model/fold_2/contribution_scores.counts_scores.bw \
no_bias_model/fold_3/contribution_scores.counts_scores.bw \
no_bias_model/fold_4/contribution_scores.counts_scores.bw \
> no_bias_model/temp_contribution_counts_score.wig
wigToBigWig no_bias_model/temp_contribution_counts_score.wig ${grefdir}/hg38.chrom.sizes no_bias_model/${celltype}_averaged_contribution_scores_counts.bw

echo "Take average of bigwig files profile"
wiggletools mean no_bias_model/fold_0/contribution_scores.profile_scores.bw \
no_bias_model/fold_1/contribution_scores.profile_scores.bw \
no_bias_model/fold_2/contribution_scores.profile_scores.bw \
no_bias_model/fold_3/contribution_scores.profile_scores.bw \
no_bias_model/fold_4/contribution_scores.profile_scores.bw \
> no_bias_model/temp_contribution_profile_score.wig
wigToBigWig no_bias_model/temp_contribution_profile_score.wig ${grefdir}/hg38.chrom.sizes no_bias_model/${celltype}_averaged_contribution_scores_profile.bw

echo "Run TFmodisco on averaged h5 contribution counts and profile files"
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


# Submit the finemo job
echo "Run finemo for motif calls"
finemo_job_id=$(bsub -J ${celltype}_finemo \
    -P acc_Tsankov_Normal_Lung \
    -q gpu \
    -n 8 \
    -W 72:00 \
    -gpu num=2 \
    -R h100nvl \
    -R rusage[mem=32000] \
    -R span[hosts=1] \
    -o ${chromBPdir}/finemo_${celltype}.out \
    -e ${chromBPdir}/finemo_${celltype}.err \
    -w "done(${TFmd_c_id}) && done(${TFmd_p_id})" \
    ${repodir}/utils/finemo_motif_calls.sh "$chromBPdir" "$celltype" \
    | awk '{print $2}' | sed 's/<//;s/>//')

echo "Submitted finemo job with ID: $finemo_job_id"

# Wait for the finemo job to complete
echo "Waiting for the finemo job to finish..."
bwait -w "done(${finemo_job_id})"

# Run the R script
source activate meso_scatac
echo "Running R script for finemo motif labels..."
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
cd $chromBPct_dir

avg_contribution_file=no_bias_model/${celltype}_averaged_contributions_counts.h5
if [ -f "${avg_contribution_file}" ]; then
echo "average contribution file found. Delete fold contribution files"
rm -r fold_0/contribution*
rm -r fold_1/contribution*
rm -r fold_2/contribution*
rm -r fold_3/contribution*
rm -r fold_4/contribution*
fi



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