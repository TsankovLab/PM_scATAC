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

### Create background regions file - Go in the output chromBPnet folder ####
mkdir $chromBPdir/$celltype
cd $chromBPdir/$celltype

chmod +x ${repodir}/utils/bias_training.sh


# Remove regions overlapping black listed regions
bedtools slop -i ${grefdir}/blacklist.bed.gz -g ${grefdir}/hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a peakset_all.bed -b temp.bed  > peakset_all_no_blacklist.bed
wc -l peakset_all_no_blacklist.bed # # Make sure number of peaks is not more than 250K

# Generate training validation and test chromosome sets
head -n 23  ${grefdir}/hg38.chrom.sizes >  hg38.chrom.subset.sizes

# mkdir splits
# chrombpnet prep splits -c hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op splits/fold_0

# Train chrombpnet bias model
for fold_number in 0 1 2 3 4; do
    
    negatives_file=${celltype}_bias_model/fold_${fold_number}/output_negatives_${fold_number}
    if [ ! -f "${negatives_file}" ]; then
    echo "negatives file not found. Identifying background peaks..."
    #rm -r output_auxiliary
    chrombpnet prep nonpeaks \
    	-g ${grefdir}/genome_references/hg38.genome.fa \
        -p peakset_all_no_blacklist.bed \
    	-c ${grefdir}/hg38.chrom.sizes \
    	-fl ${grefdir}/folds/fold_${fold_number}.json \
    	-br ${grefdir}/blacklist.bed.gz \
    	-o bias_model/output_negatives_f${fold_number}
    fi
done


# Train chrombpnet bias model
job_ids=""
for fold_number in 0 1 2 3 4; do
    BIAS_MODEL_H5=bias_model/fold_${fold_number}/models/model_bias.h5
    if [ ! -f "${BIAS_MODEL_H5}" ]; then
    rm -r bias_model/fold_$fold_number
    echo "run training bias model"
     job_id=$(
        bsub -J ${celltype}_bias \
         -P acc_Tsankov_Normal_Lung \
         -q gpu \
         -n 8 \
         -W 48:00 \
         -gpu num=2 \
         -R h100nvl \
         -R rusage[mem=32000] \
         -R span[hosts=1] \
         -o ${chromBPdir}/bias_training_${celltype}_f${fold_number}.out \
         -e ${chromBPdir}/bias_training_${celltype}_f${fold_number}.err \
         ${repodir}/utils/bias_training.sh "$chromBPdir" "$grefdir" "$celltype" "$fold_number" \
         | awk '{print $2}' | sed 's/<//;s/>//')
    
    # Append the job ID to the job_ids string
        if [ -z "$job_ids" ]; then
            job_ids="done(${job_id})"
        else
            job_ids="$job_ids && done(${job_id})"
        fi
 
    else
    echo "model_bias.h5 file found!"
    fi
done

# Wait for all jobs to complete
bsub -J wait_jobs -P acc_Tsankov_Normal_Lung -w "$job_ids"  -o wait.log -e wait.err /bin/bash -c "echo 'All jobs completed.'"

### Combine contribution scores 
echo "combine bias models"
#ml anaconda3/2020.11
#source deactivate
source activate h5py # activate another environment with hdf5plugin installed to read h5 files

chromBPct_dir=${chromBPdir}/${celltype}
echo $chromBPct_dir

# Explicitly set PATH to detect hdf5plugin
# Fix environment variables
export PATH=/sc/arion/work/giottb01/conda/envs/h5py/bin:$PATH
unset PYTHONPATH
export LD_LIBRARY_PATH=/sc/arion/work/giottb01/conda/envs/h5py/lib:$LD_LIBRARY_PATH

conda list | grep hdf5plugin
/sc/arion/work/giottb01/conda/envs/h5py/bin/python -c "import hdf5plugin; print('hdf5plugin is installed')"
/sc/arion/work/giottb01/conda/envs/h5py/bin/python $repodir/utils/average_bias_models.py $chromBPct_dir
