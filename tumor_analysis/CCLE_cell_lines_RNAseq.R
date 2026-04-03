# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

suppressPackageStartupMessages(library(data.table))

release_name <- "DepMap Public 26Q1"
depmap_dir <- "depmap_26Q1"
dir.create(depmap_dir, showWarnings = FALSE, recursive = TRUE)

get_depmap_manifest <- function() {
  manifest_url <- "https://depmap.org/portal/api/download/files"
  fread(manifest_url, data.table = FALSE)
}

find_depmap_file <- function(manifest, release_name, filename_candidates) {
  nms <- names(manifest)

  release_col <- nms[grepl("release", nms, ignore.case = TRUE)][1]
  file_col    <- nms[grepl("filename|file.*name|name", nms, ignore.case = TRUE)][1]
  url_col     <- nms[grepl("^url$|download.*url|signed.*url", nms, ignore.case = TRUE)][1]

  if (any(is.na(c(release_col, file_col, url_col)))) {
    stop("Could not identify needed columns in manifest:\n",
         paste(nms, collapse = ", "))
  }

  sub <- manifest[tolower(trimws(manifest[[release_col]])) == tolower(trimws(release_name)), , drop = FALSE]

  for (fname in filename_candidates) {
    hit <- sub[tolower(trimws(sub[[file_col]])) == tolower(trimws(fname)), , drop = FALSE]
    if (nrow(hit) > 0) {
      return(list(
        filename = hit[[file_col]][1],
        url = hit[[url_col]][1]
      ))
    }
  }

  stop("No matching filename found for release ", release_name,
       "\nCandidates: ", paste(filename_candidates, collapse = ", "))
}

download_depmap_manifest_file <- function(manifest, release_name, filename_candidates,
                                          out_dir = ".", overwrite = TRUE) {
  hit <- find_depmap_file(manifest, release_name, filename_candidates)
  dest <- file.path(out_dir, hit$filename)

  if (file.exists(dest) && overwrite) unlink(dest)

  message("Downloading: ", hit$filename)
  message("URL: ", hit$url)

  utils::download.file(
    url = hit$url,
    destfile = dest,
    mode = "wb",
    method = "libcurl"
  )

  # sanity check: reject html pages
  first_lines <- tryCatch(readLines(dest, n = 5, warn = FALSE), error = function(e) character())
  first_text <- paste(first_lines, collapse = "\n")

  if (grepl("<html|<!DOCTYPE html|gtag\\(", first_text, ignore.case = TRUE)) {
    stop("Downloaded file is HTML, not CSV: ", dest)
  }

  dest
}

strip_depmap_gene_suffix <- function(x) {
  vapply(strsplit(x, "\\.\\.", perl = TRUE), `[`, character(1), 1)
}

manifest <- get_depmap_manifest()

# inspect manifest columns if needed
print(names(manifest))

ccle_meta_file <- download_depmap_manifest_file(
  manifest,
  release_name,
  c("Models.csv", "Model.csv"),
  out_dir = depmap_dir,
  overwrite = TRUE
)

ccle_rna_file <- download_depmap_manifest_file(
  manifest,
  release_name,
  c("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv"),
  out_dir = depmap_dir,
  overwrite = TRUE
)

cnv_file <- download_depmap_manifest_file(
  manifest,
  release_name,
  c("OmicsCNGeneWGS.csv"),
  out_dir = depmap_dir,
  overwrite = TRUE
)

# verify first few lines before reading
readLines(ccle_meta_file, n = 3)
readLines(ccle_rna_file, n = 3)
readLines(cnv_file, n = 3)

# now read
ccle_meta <- fread(ccle_meta_file, data.table = FALSE)
ccle_rna  <- fread(ccle_rna_file, data.table = FALSE)
cnv       <- fread(cnv_file, data.table = FALSE)


# RNA: set rownames from first column
rna_id_col <- names(ccle_rna)[1]
rownames(ccle_rna) <- ccle_rna[[rna_id_col]]
ccle_rna[[rna_id_col]] <- NULL

# RNA subset
cellLine_ids = c('NCI-H2052','MSTO-211H','NCI-H2452','NCI-H28')
ccle_meta = ccle_meta[ccle_meta$CellLineName %in% cellLine_ids,]
ccle_rna = ccle_rna[match(ccle_meta$ModelID, ccle_rna$ModelID),]
rownames(ccle_rna) = ccle_rna$ModelID
ccle_rna = ccle_rna[,7:ncol(ccle_rna)]
ccle_rna = t(ccle_rna)
colnames (ccle_rna) = ccle_meta$CellLineName[match(colnames (ccle_rna),ccle_meta$ModelID)]
rownames(ccle_rna) <- sub("\\s*\\(\\d+\\)$", "", rownames(ccle_rna))

# CNV subset
cellLine_ids = c('NCI-H2052','MSTO-211H','NCI-H2452','NCI-H28')
cnv = cnv[match(ccle_meta$ModelID, cnv$ModelID),]
rownames(cnv) = cnv$ModelID
cnv = cnv[,7:ncol(cnv)]
cnv = t(cnv)
colnames (cnv) = ccle_meta$CellLineName[match(colnames (cnv),ccle_meta$ModelID)]
rownames(cnv) <- sub("\\s*\\(\\d+\\)$", "", rownames(cnv))



check_genes = readRDS ('selected_TF.rds')
target_TFs = c('PITX1','MEF2A','MEF2D','BPTF','TCF3','HMGA1','SOX9','TFAP2A','RUNX1','TEAD1','TEAD2','TEAD3','TEAD4','TWIST1')
target_TFs_recovered = c('TCF3','TEAD2','PITX1','SOX9','SNAI2','MEF2A','MEF2D','BPTF','TEAD4','HMGA1')
target_TFs_not_in_list = target_TFs[!target_TFs %in% check_genes]


ccle_rna2 = as.data.frame (ccle_rna[,'NCI-H2052',drop=F])
ccle_rna2$gene = rownames(ccle_rna2)
colnames(ccle_rna2)[1] = 'expression'
ccle_rna2$gene = factor (ccle_rna2$gene, levels = ccle_rna2$gene[order(-ccle_rna2$expression)])
ccle_rna2$target = ccle_rna2$gene %in% target_TFs
ccle_rna2$target[ccle_rna2$gene %in% target_TFs_recovered] = 'recovered'
ccle_rna2$target[ccle_rna2$gene %in% target_TFs_not_in_list] = 'not_in_list' 
ccle_rna2 = ccle_rna2[check_genes,]
# ccle_p1 = ggplot(ccles[ccles$viestra_cls != 'clNA',], aes(fill=cell_line, y=expression, x=gene)) + 
#     geom_bar(position="dodge", stat="identity") + scale_fill_viridis (discrete = T) + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#     facet_wrap (~viestra_cls, drop=TRUE, scales = 'free_x', ncol=12) + ggtitle ("TF in Viestra clusters") + 

ccle_p = ggplot(ccle_rna2, 
    aes(fill=target, y=log2(expression+1), x=gene)) +
    geom_bar(position="dodge", stat="identity") + 
    scale_fill_viridis (discrete = T)  + 
    theme_classic() +
    # geom_bar_pattern(position="dodge", stat="identity", pattern_fill = "black",
    #                fill = "white", colour = "black", pattern_spacing = 0.01,
    #                pattern_frequency = 5, pattern_angle = 45) +
    ggtitle ("genes to check") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf (file.path('Plots','CCLE_tf_expression_barplots.pdf'), width=12, height=3)
wrap_plots (ccle_p, ncol=1)
dev.off()


#### Import inferCNV results to generate inferCNV plot but averaging by sample ####
library (biomaRt)
library (GenomicRanges)
library (BSgenome.Hsapiens.UCSC.hg38)
library (dplyr)
# cluster infercnv heatmap not by cluster
# library (infercnv)
# infercnv_result = readRDS ('infercnv/all_samples_subsampled_Inf_ref_broad_meso/infercnv.results.obj.Rds')
# icnf_exp = infercnv_result@expr.data

#source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scATAC_functions.R')

makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6,
  gene_level_annotation = TxDb.Hsapiens.UCSC.hg38.knownGene){
  chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  mcols(windows)$wSeq <- as.character(seqnames(windows))
    mcols(windows)$wStart <- start(windows)
    mcols(windows)$wEnd <- end(windows)
  message("Subtracting Blacklist...")
  windowsBL <- lapply(seq_along(windows), function(x){
      if(x %% 100 == 0){
        message(sprintf("%s of %s", x, length(windows)))
      }
      gr <- GenomicRanges::setdiff(windows[x,], blacklist)
      mcols(gr) <- mcols(windows[x,])
      return(gr)
    })
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  windowsBL <- unlist(GRangesList(windowsBL), use.names = TRUE)
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split (windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x){
    message(sprintf("%s of %s", x, length(windowSplit)))
      chrSeq <- Biostrings::getSeq(genome,chromSizes[which(as.character(seqnames(chromSizes))==names(windowSplit)[x])])
      grx <- windowSplit[[x]]
      aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
      mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
      mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
      return(grx)
    }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  # get gene density
  gene_density = genes (gene_level_annotation)
  mcols(windowNuc)$gene_density = countOverlaps (windowNuc, gene_density)
  windowNuc
}

# Get GRanges of bins excluding black list regions
ws = 10000000
ss = ws
windowSize = ws
slidingSize = ss
genome = BSgenome.Hsapiens.UCSC.hg38
chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  
#region_df = do.call (rbind, region)
# specify the database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# loop through rows, get genes, then paste with collapse,
# and finally bind back with data d.

# best current pattern
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"   # try "uswest" or "asia" if needed
)
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
  filters    = "hgnc_symbol",
  values     = rownames(cnv),
  mart       = ensembl
)

gene_info_filtered = gene_info[gene_info$chromosome_name %in% c(1:22),]
colnames (gene_info_filtered) = c('chr','start','end')
gene_info_filtered$chr = paste0('chr',gene_info_filtered$chr)
gene_info_gr = makeGRangesFromDataFrame (gene_info_filtered, keep.extra.columns=T)

ov = findOverlaps (windows, gene_info_gr)
qhits = queryHits (ov)
shits = subjectHits (ov)
icnv_regions = lapply (unique(qhits), function(x) 
    colMeans(cnv[gene_info_gr@elementMetadata[,1][shits[which (qhits == x)]],,drop=F]))

region_names = as.character(windows)[unique(qhits)]
chr_names = sapply (region_names, function(x) unlist(strsplit (x,'\\:'))[1])
icnv_regions_df = do.call (cbind, icnv_regions)
sample_row = sapply (rownames(icnv_regions_df), function(x) unlist(strsplit (x,'\\_'))[1])
icnv_regions_df_sample = lapply (unique(sample_row), function(x) colMeans (icnv_regions_df[sample_row == x,, drop=F]))
icnv_regions_df_sample = do.call (cbind, icnv_regions_df_sample)
colnames (icnv_regions_df_sample) = unique(sample_row)
#sample_to_keep = c('p786neg','p811','p826','p846','p848','p4','p8neg','p7','p9','p11','p12','p13')
#icnv_regions_df_sample = icnv_regions_df_sample[,sample_to_keep]
rownames (icnv_regions_df_sample) = region_names
icnv_regions_df_sample = t(icnv_regions_df_sample)
colnames (icnv_regions_df_sample) = chr_names
chr_names= factor (chr_names, levels = unique(chr_names))
#rownames(icnv_regions_df_sample) = c('p786','p811','p826','p846','p848','p4','p8','p7','p9','p11','p12','p13')

icnv_regions_df_sample_c = icnv_regions_df_sample
#icnv_regions_df_sample_c[icnv_regions_df_sample_c > 2] = 2
#icnv_regions_df_sample_c[icnv_regions_df_sample_c < 0] = 0
#palette_cnv = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging",3))

library (ComplexHeatmap)
library (circlize)
ht = Heatmap (
  icnv_regions_df_sample_c,
  #row_order = match(levels(srt$sampleID), rownames(icnv_regions_df_sample)),
  column_gap = unit(.2,'mm'),
  column_names_gp = gpar(fontsize = 0),
  cluster_column_slices = FALSE,
  column_split = chr_names, 
  border=T,
  cluster_columns=F)
  #colorRamp2(c(0.8, 1, 1.2), c(palette_cnv[1], palette_cnv[2], palette_cnv[3])))
pdf (file.path('Plots','cnv_sample_mean_heatmap.pdf'),height=3, width=13)
ht
dev.off()


### Redo with expression data
gene_info <- getBM (attributes = c("chromosome_name", "start_position", "end_position","external_gene_name"),
                   filters = "external_gene_name",
                   values = rownames(ccle_rna),
                   mart = ensembl)

gene_info_filtered = gene_info[gene_info$chromosome_name %in% c(1:22),]
colnames (gene_info_filtered) = c('chr','start','end')
gene_info_filtered$chr = paste0('chr',gene_info_filtered$chr)
gene_info_gr = makeGRangesFromDataFrame (gene_info_filtered, keep.extra.columns=T)

ov = findOverlaps (windows, gene_info_gr)
qhits = queryHits (ov)
shits = subjectHits (ov)
icnv_regions = lapply (unique(qhits), function(x) 
    colMeans(log2(ccle_rna+1)[gene_info_gr@elementMetadata[,1][shits[which (qhits == x)]],,drop=F]))

region_names = as.character(windows)[unique(qhits)]
chr_names = sapply (region_names, function(x) unlist(strsplit (x,'\\:'))[1])
icnv_regions_df = do.call (cbind, icnv_regions)
sample_row = sapply (rownames(icnv_regions_df), function(x) unlist(strsplit (x,'\\_'))[1])
icnv_regions_df_sample = lapply (unique(sample_row), function(x) colMeans (icnv_regions_df[sample_row == x,, drop=F]))
icnv_regions_df_sample = do.call (cbind, icnv_regions_df_sample)
colnames (icnv_regions_df_sample) = unique(sample_row)
#sample_to_keep = c('p786neg','p811','p826','p846','p848','p4','p8neg','p7','p9','p11','p12','p13')
#icnv_regions_df_sample = icnv_regions_df_sample[,sample_to_keep]
rownames (icnv_regions_df_sample) = region_names
icnv_regions_df_sample = t(icnv_regions_df_sample)
colnames (icnv_regions_df_sample) = chr_names
chr_names= factor (chr_names, levels = unique(chr_names))
#rownames(icnv_regions_df_sample) = c('p786','p811','p826','p846','p848','p4','p8','p7','p9','p11','p12','p13')

icnv_regions_df_sample_c = icnv_regions_df_sample
#icnv_regions_df_sample_c[icnv_regions_df_sample_c > 2] = 2
#icnv_regions_df_sample_c[icnv_regions_df_sample_c < 0] = 0
#palette_cnv = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging",3))

library (ComplexHeatmap)
library (circlize)
icnv_regions_df_sample_c2 = icnv_regions_df_sample_c
ht = Heatmap (
  scale(icnv_regions_df_sample_c2),
  #row_order = match(levels(srt$sampleID), rownames(icnv_regions_df_sample)),
  column_gap = unit(.2,'mm'),
  column_names_gp = gpar(fontsize = 0),
  cluster_column_slices = FALSE,
  column_split = chr_names, 
  border=T,
  cluster_columns=F)
  #colorRamp2(c(0.8, 1, 1.2), c(palette_cnv[1], palette_cnv[2], palette_cnv[3])))
pdf (file.path('Plots','cnv_sample_expression_mean_heatmap.pdf'),height=3, width=13)
ht
dev.off()
