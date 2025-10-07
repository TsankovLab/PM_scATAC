cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA

python code/03-chrombpnet/02-compendium/01-modisco_to_pfm.py \
  -c ../main/scatac_ArchR/chrombpnet_models_paths.tsv \
  -o ../main/scatac_ArchR/chromBPnet/pfms.txt \
  -p pos_patterns \
  -t 0.3 \
  -ml 3


python code/03-chrombpnet/02-compendium/01-modisco_to_pfm.py \
  -c ../main/scatac_ArchR/chrombpnet_models_paths.tsv \
  -o ../main/scatac_ArchR/chromBPnet/pfms_neg.txt \
  -p neg_patterns \
  -t 0.3 \
  -ml 3


# Positive motifs
gimme cluster \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/pfms.txt \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/cluster_pos \
  -t 0.8 \
  -N 8

# Negative motifs
gimme cluster \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/pfms_neg.txt \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/cluster_neg \
  -t 0.8 \
  -N 8

output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/modisco_merged_counts
model_head=counts



BATCH_SIZE=5
bias_params="Heart_c0_thresh0.4"
out_dir=${modisco_merged_dir}
model_head="counts"

output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet
mkdir $output_dir

model_head=counts

awk -F'\t' 'BEGIN{OFS="\t"} {print $1, "pos_patterns", $2, 1}' /sc/arion/scratch/giottb01/chromBPnet/cluster_pos/cluster_key.txt > /sc/arion/scratch/giottb01/chromBPnet/cluster_pos/cluster_key_modified_pos.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, "pos_patterns", $2, 1}' /sc/arion/scratch/giottb01/chromBPnet/cluster_neg/cluster_key.txt > /sc/arion/scratch/giottb01/chromBPnet/cluster_neg/cluster_key_modified_neg.txt
cat /sc/arion/scratch/giottb01/chromBPnet/cluster_pos/cluster_key_modified_pos.txt /sc/arion/scratch/giottb01/chromBPnet/cluster_neg/cluster_key_modified_neg.txt > ${output_dir}/cluster_key_combined.txt


cluster_key=${output_dir}/cluster_key_combined.txt


modisco_dir=/sc/arion/scratch/giottb01/chromBPnet
contribs_dir=/sc/arion/scratch/giottb01/chromBPnet
batch=1

# set inputs
#cluster_key="${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv"

# SUBMIT JOBS ------------------------------------------------------------------

# get batches
batches=( $(cut -f 4 $cluster_key | sort | uniq) )

for batch in ${batches[@]}; do
# Merge modisco clusters (?)
python -u ../git_repo/utils/03-merge_modisco.py --out-dir ${output_dir} \
                    --model-head ${model_head} \
                    --cluster-key ${cluster_key} \
                    --modisco-dir ${modisco_dir} \
                    --contribs-dir ${contribs_dir} \
                    --batch ${batch}

done                      
