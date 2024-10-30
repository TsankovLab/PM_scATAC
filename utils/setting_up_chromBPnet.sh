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

# install TF-modisco lite
pip install modisco-lite

# Activate conda env with GPU modules 
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

# Train bias model #### I did this only on fold 0 !!
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






