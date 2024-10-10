
#### Reproduce tutorial ####
conda activate chrombpnet
mkdir /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/tutorial/data/download
cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/tutorial/

# Download bam files
wget https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam -O data/download/rep1.bam
wget https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam -O data/download/rep2.bam
wget https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam -O data/download/rep3.bam

# download overlap peaks (default peaks on ENCODE)
wget https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz -O overlap.bed.gz

# Ensure that the peak regions do not intersect with the blacklist regions (extended by 1057 bp on both sides) by running the command below.
bedtools slop -i blacklist.bed.gz -g ../hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a overlap.bed.gz -b temp.bed  > peaks_no_blacklist.bed

head -n 24  ../hg38.chrom.sizes >  hg38.chrom.subset.sizes

mkdir data/splits
chrombpnet prep splits -c hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op data/splits/fold_0

# create non-peaks regions
chrombpnet prep nonpeaks -g ../../genome_references/hg38.genome.fa -p peaks_no_blacklist.bed -c  ../hg38.chrom.sizes -fl data/splits/fold_0.json -br ../blacklist.bed.gz -o data/output

## Try replacing values in 4:9 columns with . as is in my peakset
cp peaks_no_blacklist.bed peaks_no_blacklist2.bed
awk -v FS='[[:space:]]+' 'BEGIN{OFS="\t"} {$4=".";$1=$1} 1' peaks_no_blacklist.bed > peaks_no_blacklist2.bed
awk -v FS='[[:space:]]+' 'BEGIN{OFS="\t"} {$5=".";$1=$1} 1' peaks_no_blacklist2.bed > peaks_no_blacklist3.bed
awk -v FS='[[:space:]]+' 'BEGIN{OFS="\t"} {$6=".";$1=$1} 1' peaks_no_blacklist3.bed > peaks_no_blacklist4.bed
awk -v FS='[[:space:]]+' 'BEGIN{OFS="\t"} {$7=".";$1=$1} 1' peaks_no_blacklist4.bed > peaks_no_blacklist5.bed
awk -v FS='[[:space:]]+' 'BEGIN{OFS="\t"} {$8=".";$1=$1} 1' peaks_no_blacklist5.bed > peaks_no_blacklist6.bed
awk -v FS='[[:space:]]+' 'BEGIN{OFS="\t"} {$9=".";$1=$1} 1' peaks_no_blacklist6.bed > peaks_no_blacklist7.bed
less peaks_no_blacklist7.bed | head

# create non-peaks regions
chrombpnet prep nonpeaks -g ../../genome_references/hg38.genome.fa -p peaks_no_blacklist7.bed -c  ../hg38.chrom.sizes -fl data/splits/fold_0.json -br ../blacklist.bed.gz -o data/output2


