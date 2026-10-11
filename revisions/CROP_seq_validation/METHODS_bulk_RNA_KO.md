# Methods — bulk RNA-seq of TF knockouts in mesothelioma cell lines

Items in [brackets] are wet-lab details that are not recorded in the analysis files and need to be filled in.

## Samples and sequencing

Four mesothelioma cell lines (NCI-H2052, MSTO-211H, NCI-H2452, NCI-H28) were transduced with guides against eight transcription factors from the CROP-seq screen (TCF3, TEAD4, TWIST1, PITX1, SOX9, TEAD2, BPTF, HMGA1), three additional factors (RUNX1, RUNX2, SNAI2) or non-targeting controls (NTC), with two to three replicates per condition and line [guide sequences, Cas9 system, selection, time of harvest, RNA extraction and library kit]. Libraries were dual-indexed and sequenced paired-end, 2 × 151 bp, on one lane of an Illumina [NovaSeq X] flow cell.

## Pre-processing

Base calls were demultiplexed with BCL Convert v4.4.6, allowing one mismatch per index. Read quality was assessed with FastQC v0.12.1 and summarised with MultiQC v1.35. Reads were aligned to the human genome (hg38) with STAR v2.7.11b, using an index built with the GENCODE v46 annotation (`--sjdbOverhang 150`) and default alignment parameters; no adapter trimming was applied. Gene-level counts were obtained in the same run (`--quantMode GeneCounts`). Library strandedness was determined per sample from the totals of the three STAR count columns; all libraries were reverse-stranded and the corresponding column was used. Counts of Ensembl genes sharing a gene symbol were summed (61,471 symbols).

## Quality control and sample sets

Every sample carrying an NTC or one of the eight screen guides in the four lines (76 samples) was audited on: input read pairs, uniquely mapped fraction, reads assigned to genes, fraction of reads assigned, genes detected, mitochondrial fraction and the share of counts in the 100 most abundant genes. Samples were flagged as low depth (< 5 million reads assigned to genes) or as QC outliers (robust z-score, median/MAD within the cell line, beyond 3 on the unfavourable side for any metric), and as expression outliers (robust z < −3 for the median Spearman correlation with the other samples of the line, or > 3 for the distance to the line centroid in the first five principal components). Principal component analysis was run per line on TMM-normalised log-CPM, without and with regression of the cell-cycle score and sequencing depth, and the variance of each component explained by KO, replicate, depth, assigned fraction and cell cycle was recorded.

Cell-line identity was verified from the RNA-seq reads independently of expression. SNPs were called jointly over the exons of the 3,000 protein-coding genes with the highest minimum expression across lines (bcftools v1.22 `mpileup`/`call`; uniquely mapped reads, base quality ≥ 20, biallelic SNPs with QUAL ≥ 30). A genotype was assigned at sites covered by at least 10 reads (alternative-allele fraction < 0.1, 0.25–0.75 or > 0.9), and each sample was assigned to the line whose reference genotype (majority genotype of the other samples named after that line) it matched best. Samples whose genotype disagreed with their name were relabelled to the genotyped line, and one sample showing a mixture of two lines was excluded.

The main analysis set comprised 70 samples (NTC and the eight screen KOs; 16–19 per line; NTC: 1 in NCI-H2052, 3 in MSTO-211H, 3 in NCI-H2452, 2 in NCI-H28). In this set the median depth was 37.9 million read pairs (range 7.5–159.2), with a median of 80.9% uniquely mapped (57.0–90.4%) and 27.3 million read pairs assigned to genes (5.2–128.9; 73.0% of reads). Results were also recomputed on a stricter set of 60 samples excluding those flagged in the QC and batch review.

The i5 library index group, the only batch-related variable in the sample sheet, was recorded for every sample. Within a line both replicates of a KO share the same group, so it could not be included as a covariate.

## Normalisation and KO effects

Analyses were done in R v4.3.3 with edgeR v4.0.16, limma v3.58.1 and Seurat v5.2.1. Genes passing `filterByExpr` in all four lines were retained (15,630 genes). Counts were normalised within each line by TMM and expressed as log2 CPM (prior count 2). Where stated, expression was additionally quantile-normalised across all samples (`limma::normalizeQuantiles`).

For each line, a linear model with one coefficient per group (NTC and each KO) was fitted with limma, optionally including the cell-cycle score and the fraction of reads assigned to genes as covariates. Two contrasts were tested per KO with empirical-Bayes moderated t-statistics: NTC versus the KO, and the mean of the other KOs versus the KO; the second does not use the NTC samples and removes the response common to all KOs. Knockdown was assessed from the log2 ratio of the target transcript between NTC and the KO in each line.

## Reproducibility across cell lines

For every pair of lines, the KO profiles (moderated t-statistics, centred across the eight KOs within each line) were compared by Spearman correlation over the 500 genes with the largest KO signal (maximum absolute t over KOs) in the two lines not included in the pair, so that gene selection was independent of the pair being scored. Reproducibility of a KO was the mean correlation between its profiles in different lines (six line pairs), compared with its correlations to the other KOs. Significance was obtained by permuting the KO labels independently within each line (20,000 permutations), per KO and for the mean over KOs, with Benjamini–Hochberg correction over the eight KOs. Robustness to analysis choices was evaluated over a grid of sample set, normalisation, covariates, statistic, feature space and correlation method (12,672 configurations).

## Gene programme (Cm) scores and comparison with CROP-seq

The 20 gene programmes (Cm1–Cm20) defined from single-cell RNA-seq of pleural mesothelioma were scored in each sample with Seurat `AddModuleScore` on the top 50 genes of each programme; as an alternative, the mean of within-line gene z-scores over the top 30, 50 or 100 genes was used. Scores were z-scored within each line, replicates were averaged per group, and the effect of a KO on a programme was computed as the KO minus NTC or as the KO minus the mean of the other eight groups. Effects were combined across lines with equal weight per line. For a given KO–programme pair, significance was assessed by permuting the group label within each line (10,000 permutations), with Benjamini–Hochberg correction over the 160 pairs.

In the CROP-seq screen (NCI-H2052), the same programmes were scored per cell and the same contrasts computed from group means. Agreement between screen and bulk was measured by Spearman correlation over the 160 KO–programme effects, for each line and for the four-line mean, and separately for the component shared by all KOs (the mean effect per programme) and the KO-specific remainder. The reliability of each data set was estimated from the correlation between random halves of the cells of each guide in the screen (300 splits) and between replicates within each bulk line, with Spearman–Brown correction; the upper bound of the attainable agreement was taken as the square root of the product of the two reliabilities. The effect of each processing step on these measures was evaluated over 1,280 combinations of sample set, normalisation (CP10K, TMM, TMM with quantile normalisation, RLE, upper quartile), programme scoring, covariates, within-line scaling, reference and statistic.

## Proliferation score

Proliferation was scored per sample as the sum of the S-phase and G2/M scores from Seurat `CellCycleScoring` (default gene lists) on log-normalised counts, with or without prior quantile normalisation across samples. The difference between the mean of the NTC samples and the mean of the KO replicates was computed within each line. For display, scores were z-scored within each line across its NTC and KO samples and pooled, and each KO was compared with NTC by a two-sided Wilcoxon rank-sum test. The same score was computed per cell in the CROP-seq screen, z-scored across cells, and compared between the cells of each KO and the NTC cells in the same way.
