### Set up chromBPnet #### from https://github.com/kundajelab/chrombpnet?tab=readme-ov-file

# Request GPU cores
bsub -q interactive -P acc_Tsankov_Normal_Lung -n 1 -W 12:00 -gpu num=1 -R v100 -R rusage[mem=64000] -R span[hosts=1] -XF -Is /bin/bash

# Install conda 
conda create -n chrombpnet python=3.8
conda activate chrombpnet

# Install supplementary packages 
conda install -y -c conda-forge -c bioconda samtools bedtools ucsc-bedgraphtobigwig pybigwig meme

# install chromBPnet
pip install chrombpnet

conda activate chrombpnet
ml cuda
ml cudnn/8.2.0
ml proxies


# Create background regions file - Go in the output chromBPnet folder
cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet

# Remove regions overlapping black listed regions
bedtools slop -i ../../../../../chromBPnet/blacklist.bed.gz -g ../../../../../chromBPnet/hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a peakset_all.bed -b temp.bed  > peakset_all_no_blacklist.bed
wc -l peakset_all_no_blacklist.bed # # Make sure number of peaks is not more than 250K

# Generate training validation and test chromosome sets
head -n 23  ../../../../../chromBPnet/hg38.chrom.sizes >  hg38.chrom.subset.sizes

mkdir splits
chrombpnet prep splits -c hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op splits/fold_0

chrombpnet prep nonpeaks -g ../../../../../genome_references/hg38.genome.fa -p peakset_all_no_blacklist.bed -c  ../../../../../chromBPnet/hg38.chrom.sizes -fl splits/fold_0.json -br ../../../../../chromBPnet/blacklist.bed.gz -o output


# shuf -n 200000 peakset_all_no_blacklist.bed > peakset_all_no_blacklist_shuffled.bed
# wc -l peakset_all_no_blacklist_shuffled.bed

# chrombpnet prep nonpeaks -g ../../../../../genome_references/hg38.genome.fa -p peakset_all_no_blacklist_shuffled.bed -c  ../../../../../chromBPnet/hg38.chrom.sizes -fl splits/fold_0.json -br ../../../../../chromBPnet/blacklist.bed.gz -o output2

# Train bias model
chrombpnet bias pipeline \
        -ifrag fragments_NKT_cells.tsv \
        -d "ATAC" \
        -g ../../../../../genome_references/hg38.genome.fa \
        -c ../../../../../chromBPnet/hg38.chrom.sizes \
        -p peakset_all_no_blacklist.bed \
        -n output_negatives.bed \
        -fl splits/fold_0.json \
        -b 0.5 \
        -o bias_model/ \
        -fp tnk

# Run bias model as a job
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombpnet_bias_NKT.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombpnet_bias_NKT.sh



# Train CD8 exhausted chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext.sh

# Train CD8 exhausted chrombpnet model F1
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f1.sh

# Train CD8 exhausted chrombpnet model F2
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f2.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f2.sh

# Train CD8 exhausted chrombpnet model F3
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f3.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f3.sh

# Train CD8 exhausted chrombpnet model F4
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_ext_f4.sh




# Train NK KRLC1 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1.sh

# Train NK KRLC1 chrombpnet model F1
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f1.sh

# Train NK KRLC1 chrombpnet model F2
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f2.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f2.sh

# Train NK KRLC1 chrombpnet model F3
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh

# Train NK KRLC1 chrombpnet model F4
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f4.sh



# Train NK FGFBP2 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2.sh

# Train NK FGFBP2 chrombpnet model F1
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f1.sh

# Train NK FGFBP2 chrombpnet model F2
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f2.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f2.sh

# Train NK FGFBP2 chrombpnet model F3
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f3.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f3.sh

# Train NK FGFBP2 chrombpnet model F4
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_FGFBP2_f4.sh

# Train CD8 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8.sh

# Train CD8 chrombpnet model F1
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f1.sh

# Train CD8 chrombpnet model F2
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f2.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f2.sh

# Train CD8 chrombpnet model F3
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f3.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f3.sh

# Train CD8 chrombpnet model F4
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD8_f4.sh

# Train CD4 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4.sh

# Train CD4 chrombpnet model F1
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f1.sh

# Train CD4 chrombpnet model F2
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f2.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f2.sh

# Train CD4 chrombpnet model F3
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f3.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f3.sh

# Train CD4 chrombpnet model F4
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_CD4_f4.sh

# Train Tregs chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs.sh

# Train Tregs chrombpnet model F1
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f1.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f1.sh

# Train Tregs chrombpnet model F2
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f2sh.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f2sh.sh

# Train Tregs chrombpnet model F3
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f3.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f3.sh

# Train Tregs chrombpnet model F4
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f4.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_Tregs_f4.sh



# Get contribution scores from NK KRLC1 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_KLRC1_contribution.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_KLRC1_contribution.sh

# Get contribution scores from CD8 ext chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_CD8ext_contribution.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_CD8ext_contribution.sh

# Get contribution scores from NK Tregs chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_Tregs_contribution.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_Tregs_contribution.sh

# Get contribution scores from CD8 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_CD8_contribution.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_CD8_contribution.sh

# Get contribution scores from CD4 chrombpnet model
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_CD4_contribution.sh
bsub </sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chromBPnet_NKT_CD4_contribution.sh


# Get contribution score bigwigs
cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet
MODEL_H5=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/KRLC1_model/models/chrombpnet.h5
REGIONS=peakset_NK_KLRC1.bed
GENOME=../../../../../genome_references/hg38.genome.fa
CHROM_SIZES=../../../../../chromBPnet/hg38.chrom.sizes
OUTPUT_PREFIX=KLRC1_contribution_scores
chrombpnet contribs_bw -m $MODEL_H5 -r $REGIONS -g $GENOME -c $CHROM_SIZES -op $OUTPUT_PREFIX 





