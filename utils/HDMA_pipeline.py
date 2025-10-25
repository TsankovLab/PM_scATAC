conda activate chrombpnet

cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA
#output_dir=/sc/arion/scratch/giottb01/chromBPnet
chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet
chromBPdir=/sc/arion/scratch/giottb01/chromBPnet

model_head=profile
mkdir ${chromBPdir}/modisco_merged_${model_head}

python code/03-chrombpnet/02-compendium/01-modisco_to_pfm.py \
  -c ${chromBPdir}/chrombpnet_models_${model_head}_paths.tsv \
  -o ${chromBPdir}/modisco_merged_${model_head}/pfms_pos.txt \
  -p pos_patterns \
  -t 0.3 \
  -ml 3


python code/03-chrombpnet/02-compendium/01-modisco_to_pfm.py \
  -c ${chromBPdir}/chrombpnet_models_${model_head}_paths.tsv \
  -o ${chromBPdir}/modisco_merged_${model_head}/pfms_neg.txt \
  -p neg_patterns \
  -t 0.3 \
  -ml 3


# Positive motifs
gimme cluster \
  ${chromBPdir}/modisco_merged_${model_head}/pfms_pos.txt \
  ${chromBPdir}/modisco_merged_${model_head}/cluster_pos \
  -t 0.8 \
  -N 8

# Negative motifs
gimme cluster \
  ${chromBPdir}/modisco_merged_${model_head}/pfms_neg.txt \
  ${chromBPdir}/modisco_merged_${model_head}/cluster_neg \
  -t 0.8 \
  -N 8

#output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/modisco_merged_counts



### Merge modisco motifs per model head ###
#BATCH_SIZE=5
#bias_params="Heart_c0_thresh0.4"
#out_dir=${modisco_merged_dir}
#model_head="counts"


awk -F'\t' 'BEGIN{OFS="\t"} {n=int((NR-1)/3)+1; print $1, "pos_patterns", $2, n}' ${chromBPdir}/modisco_merged_${model_head}/cluster_pos/cluster_key.txt > ${chromBPdir}/modisco_merged_${model_head}/cluster_pos/cluster_key_modified_pos.txt
awk -F'\t' 'BEGIN{OFS="\t"} {n=int((NR-1)/3)+1; print $1, "neg_patterns", $2, n}' ${chromBPdir}/modisco_merged_${model_head}/cluster_neg/cluster_key.txt > ${chromBPdir}/modisco_merged_${model_head}/cluster_neg/cluster_key_modified_neg.txt
cat ${chromBPdir}/modisco_merged_${model_head}/cluster_pos/cluster_key_modified_pos.txt ${chromBPdir}/modisco_merged_${model_head}/cluster_neg/cluster_key_modified_neg.txt > ${chromBPdir}/modisco_merged_${model_head}/cluster_key_combined.txt


cluster_key=${chromBPdir}/modisco_merged_${model_head}/cluster_key_combined.txt


modisco_dir=${chromBPdir}
contribs_dir=${chromBPdir}

# set inputs
#cluster_key="${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv"

# SUBMIT JOBS ------------------------------------------------------------------

# get batches
batches=$(cut -f 4 $cluster_key | sort | uniq)

chmod +x ../git_repo/utils/modisco_merge_job.sh

for batch in ${batches[@]}; do

bsub -J modisco_merge \
-P acc_Tsankov_Normal_Lung \
-q premium \
-n 8 \
-W 12:00 \
-R rusage[mem=80000] \
-R span[hosts=1] \
-o ${chromBPdir}/modisco_merged_${model_head}/modisco_merge_${batch}.out \
-e ${chromBPdir}/modisco_merged_${model_head}/modisco_merge_${batch}.err \
../git_repo/utils/modisco_merge_job.sh "$chromBPdir" "$model_head" "$cluster_key" "$modisco_dir" "$contribs_dir" "$batch"

done

### Compile merged modisco motifs ####
mkdir $chromBPdir/modisco_merged_${model_head}/compiled

python ../git_repo/utils/04a-compile_modisco_obj.py \
  --output_dir $chromBPdir/modisco_merged_${model_head} \
  --cluster_key_path $cluster_key \
  --compiled_modisco_h5_path ${chromBPdir}/modisco_merged_${model_head}/compiled/modisco_compiled.h5 \
  --compiled_modisco_tsv_path ${chromBPdir}/modisco_merged_${model_head}/compiled/modisco_compiled.tsv \
  --modisco_merged_dir ${chromBPdir}/modisco_merged_${model_head}/merged_motifs



### Find known motifs matching modisco motifs using TOMTOM tool ####
conda activate chrombpnet

compiled_h5=${chromBPdir}/modisco_merged_${model_head}/compiled/modisco_compiled.h5
motif_db=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
tomtom_dir=${chromBPdir}/modisco_merged_${model_head}/compiled_tomtom
mkdir $tomtom_dir
python -u code/03-chrombpnet/02-compendium/04b-get_tomtom_matches.py --modisco-h5 $compiled_h5 \
    --out-dir  $tomtom_dir \
    --meme-db $motif_db \
    --verbose True
                    
#echo "done."

### Run Fi-Nemo on averaged motifs 
chmod +x ../git_repo/utils/finemo_motif_calls_on_merged.sh
chromBPdir=/sc/arion/scratch/giottb01/chromBPnet
#output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/finemo_on_merged_${model_head}
#finemo_counts_file=no_bias_model/finemo_out_counts/hits.bed
#if [ -f "${modisco_counts_file}" ] && [ -f "${modisco_profile_file}" ] && [ ! -f "$finemo_counts_file" ]; then
celltypes=("Myeloid" "Malignant" "Fibroblasts" "Endothelial" "B_cells" \
           "Mesothelium" "SmoothMuscle" "T_cells" "NK" "Plasma" \
           "pDCs" "Alveolar")

# celltypes=("NK" "Plasma" \
#            "pDCs" "Alveolar")

for celltype in ${celltypes[@]}; do

    echo "=== Submit finemo job ==="
    bsub -J ${celltype}_finemo \
        -P acc_Tsankov_Normal_Lung \
        -q gpu \
        -n 1 \
        -W 24:00 \
        -gpu num=1 \
        -R a100 \
        -R rusage[mem=64000] \
        -R span[hosts=1] \
        -o ${chromBPdir}/finemo_${celltype}.out \
        -e ${chromBPdir}/finemo_${celltype}.err \
        ../git_repo/utils/finemo_motif_calls_on_merged.sh "$chromBPdir" "$celltype" "$model_head" \
        | awk '{print $2}' | sed 's/<//;s/>//'
    echo "Submitted finemo job with ID: $finemo_job_id"
    
    # Wait for finemo to finish
    #echo "Waiting for finemo job..."
    #bwait -w "done(${finemo_job_id})"
done

echo "=== Run R script for finemo motif labels ==="
conda activate meso_scatac

for celltype in ${celltypes[@]}; do

Rscript ../git_repo/utils/chromBPnet_finemo_motif_labels.R $chromBPdir $celltype

done

echo "R script execution completed."

### Copy chromBPnet results files (Excluding large files !!) 
chromBPdir_dest=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet
cd $chromBPdir
for celltype in ${celltypes[@]}; do
    mkdir -p "${chromBPdir_dest}/${celltype}/no_bias_model"
    rsync -av \
        --include="*/" \
        --include="*.bw" \
        --include="*.tsv" \
        --include="fold_0/***" \
        --include="fold_1/***" \
        --include="fold_2/***" \
        --include="fold_3/***" \
        --include="fold_4/***" \
        --exclude="*" \
        "${celltype}/no_bias_model/" \
        "${chromBPdir_dest}/${celltype}/no_bias_model/"
done

cp -r ${chromBPdir}/modisco_merged_counts ${chromBPdir_dest}/
cp -r ${chromBPdir}/modisco_merged_profile ${chromBPdir_dest}/

### Remove log files ####
cd $chromBPdir_dest
find . -type f \( -name "*.out" -o -name "*.err" \) -delete

