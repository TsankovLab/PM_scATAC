###############################################################
# R2_Q3 : PER-SAMPLE CLONALITY from scRNA inferCNV (malignant cells)
#
# Modality 2 = scRNA. inferCNV already produced a per-cell, per-gene residual
# CNV signal (expr.data, genes ordered by genome) using a normal reference.
# Here, for each tumor sample we (1) recompute clonal structure with the SAME
# method used for scATAC (hierarchical clustering on 1 - Pearson across cells,
# #subclones by silhouette) so the two modalities are directly comparable, and
# (2) collapse the gene-level signal into 1-Mb bins to produce a per-sample
# pseudobulk CNV profile that can be correlated with the scATAC profile.
#
# NOTE on its role in the current pipeline: downstream, only the per-sample CELL LISTS
# it writes (scRNA_cnv/<S>_subclones.csv, column `cell`) are consumed -- R2_Q3_arm_level_cnv.R
# uses them to know which scRNA cells belong to which tumour. The CNV values and clone
# labels here are superseded by the arm-level calls. Run it first regardless: step 2
# fails without those files.
###############################################################
suppressMessages({ library(infercnv); library(ComplexHeatmap); library(circlize); library(cluster) })
set.seed(1)

SC   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
IC   <- file.path(SC, "tumor_compartment", "scrna", "infercnv")
OUT  <- file.path(SC, "git_repo_claude", "R2_Q3")
setwd(OUT); dir.create("Plots", showWarnings=FALSE); dir.create("scRNA_cnv", showWarnings=FALSE)

BIN      <- 1e6          # pseudobulk bin width
MINCELLS <- 50           # a sample below this is not worth clustering
AUTOSOMES<- paste0("chr", 1:22)

## ---- STEP 1: load inferCNV and put genes in genome order -------------------
## expr.data is the per-gene CNV RESIDUAL (centred near 1), already denoised against
## inferCNV's own normal reference. Autosomes only, then sorted chr1..chr22 by start so
## that column position on every heatmap is genomic position.
ic <- readRDS(file.path(IC, "run.final.infercnv_obj"))
E  <- ic@expr.data                                  # genes x cells (CNV residual, ~1-centered)
go <- ic@gene_order
gchr   <- as.character(go$chr); gstart <- as.integer(go$start)
keepg  <- gchr %in% AUTOSOMES
E <- E[keepg, , drop=FALSE]; gchr <- gchr[keepg]; gstart <- gstart[keepg]
gord   <- order(match(gchr, AUTOSOMES), gstart)
E <- E[gord, , drop=FALSE]; gchr <- gchr[gord]; gstart <- gstart[gord]
cat("inferCNV genes (autosomal):", nrow(E), " cells:", ncol(E), "\n")

## ---- STEP 2: sample membership ---------------------------------------------
## inferCNV stores which columns belong to which tumour group; "reference" is the
## normal set it was run against and is not a tumour.
obs <- ic@observation_grouped_cell_indices          # sample -> column indices
samples <- setdiff(names(obs), "reference")
cat("tumor samples:", paste(samples, collapse=", "), "\n")

## gene -> 1 Mb bin key, zero-padded so alphabetical order is genomic order
gbin <- paste0(gchr, ":", sprintf("%05d", (gstart %/% BIN)))

summ <- list(); pbulk <- list()
for (s in samples){
  idx <- obs[[s]]; cells <- colnames(E)[idx]
  if (length(cells) < MINCELLS){ cat(sprintf("[%s] %d cells < %d -> skip\n", s, length(cells), MINCELLS)); next }
  M <- E[, idx, drop=FALSE]                          # genes x cells for this sample
  ## ---- STEP 3: per-cell centring -------------------------------------------
  ## express as deviation from 1 (inferCNV's neutral value), then subtract each cell's
  ## own median so a cell with globally shifted signal does not dominate the distance.
  ## Same convention as the scATAC side, which is what makes the two comparable.
  Mc <- sweep(M - 1, 2, apply(M - 1, 2, median), "-")
  cat(sprintf("\n[%s] %d cells\n", s, ncol(Mc)))

  ## ---- STEP 4: clone calling, identical recipe to the scATAC side ----------
  ## distance = 1 - Pearson between cells; ward.D2; k chosen by silhouette over k=2..6,
  ## a candidate k is only admissible if >=2 clusters hold max(20, 5%) cells, and a best
  ## silhouette < 0.10 means CLONAL. Two heterogeneity summaries are recorded alongside:
  ##   het_pair = mean pairwise (1 - correlation)   [cell-cell dissimilarity]
  ##   het_bin  = mean per-gene SD across cells     [per-feature spread]
  ## Neither is comparable ACROSS modalities (scATAC is inflated by sparsity) -- use the
  ## discrete call and the pseudobulk profile for that.
  cormat <- cor(Mc)
  d  <- as.dist(1 - cormat)
  hc <- hclust(d, method="ward.D2")
  het_pair <- mean(1 - cormat[upper.tri(cormat)])
  het_bin  <- mean(apply(Mc, 1, sd))
  best_k <- 1; best_sil <- -1; assign <- rep(1, ncol(Mc))
  for (k in 2:min(6, ncol(Mc)-1)){
    cl  <- cutree(hc, k=k)
    sil <- mean(silhouette(cl, d)[, "sil_width"])
    tab <- table(cl); frac_ok <- sum(tab >= max(20, 0.05*ncol(Mc)))
    if (sil > best_sil && frac_ok >= 2){ best_sil <- sil; best_k <- k; assign <- cl }
  }
  ## k_eff = number of clusters that are actually robust (>= max(20,5%) cells). The raw
  ## cutree label count can exceed this, so downstream code must use k_eff, never the
  ## number of distinct values in `assign`.
  tb_as <- table(assign); robust <- names(tb_as)[tb_as >= max(20, 0.05*ncol(Mc))]
  k_eff <- length(robust); if (best_sil < 0.10) k_eff <- 1
  largest <- max(tb_as)/ncol(Mc)
  cls <- if (k_eff <= 1) "clonal" else sprintf("subclonal (k=%d)", k_eff)
  cat(sprintf("[%s] het_pair=%.3f het_bin=%.3f best_k=%d sil=%.2f k_eff=%d largest=%.2f -> %s\n",
              s, het_pair, het_bin, best_k, best_sil, k_eff, largest, cls))

  ## ---- STEP 5: the file the rest of the pipeline needs ---------------------
  ## `cell` is what R2_Q3_arm_level_cnv.R reads; `subclone` is the raw cutree label.
  write.csv(data.frame(cell=colnames(Mc), subclone=assign),
            file.path("scRNA_cnv", paste0(s,"_subclones.csv")), row.names=FALSE)
  summ[[s]] <- data.frame(sample=s, modality="scRNA", n_cells=ncol(Mc),
                          het_pair=het_pair, het_bin=het_bin, k_subclones=k_eff,
                          best_sil=best_sil, largest_frac=largest, class=cls, stringsAsFactors=FALSE)

  ## ---- STEP 6: pseudobulk 1 Mb profile --------------------------------------
  ## average over cells first, then over genes within a bin, giving one CNV landscape per
  ## tumour that can be correlated against the scATAC landscape.
  gmean <- rowMeans(Mc)
  agg <- tapply(gmean, factor(gbin, levels=unique(gbin)), mean)
  bk  <- names(agg)
  pbulk[[s]] <- data.frame(chr=sub(":.*","",bk),
                           start=as.integer(sub(".*:","",bk))*BIN,
                           mean_cnv=as.numeric(agg), row.names=NULL)
  names(pbulk[[s]])[3] <- s

  ## ---- STEP 7: per-sample heatmap, cells grouped by subclone ---------------
  ## rows = cells (subsampled to 1500 for rendering), columns = genes in genome order,
  ## split by chromosome; neither axis is clustered so position is interpretable.
  pc  <- if (ncol(Mc) > 1500) sample(ncol(Mc), 1500) else seq_len(ncol(Mc))
  ordc<- pc[order(assign[pc])]
  ca  <- factor(assign[ordc])
  col_fun <- colorRamp2(c(-0.15,0,0.15), c("#2166ac","white","#b2182b"))
  ha  <- HeatmapAnnotation(subclone=ca,
             col=list(subclone=setNames(scales::hue_pal()(nlevels(ca)), levels(ca))),
             which="row", show_annotation_name=FALSE)
  chsplit <- factor(gchr, levels=AUTOSOMES)
  ht <- Heatmap(t(Mc[, ordc]), name="CNV", col=col_fun,
                cluster_rows=FALSE, cluster_columns=FALSE,
                show_row_names=FALSE, show_column_names=FALSE,
                column_split=chsplit, column_gap=unit(0,"mm"),
                column_title_gp=gpar(fontsize=6), border=TRUE, left_annotation=ha,
                row_title=sprintf("%s | %d cells | %s", s, ncol(Mc), cls),
                use_raster=TRUE, raster_quality=2)
  pdf(file.path("Plots", paste0("scRNA_cnv_heatmap_", s, ".pdf")), width=11, height=5.5)
  draw(ht, column_title=sprintf("scRNA inferCNV — %s", s))
  dev.off()
}

## ---- STEP 8: cohort summary + merged pseudobulk table ----------------------
summ_df <- do.call(rbind, summ)
write.csv(summ_df, "scRNA_clonality_summary.csv", row.names=FALSE)
cat("\n================ scRNA inferCNV per-sample clonality summary ================\n")
print(summ_df, row.names=FALSE, digits=3)

pb <- Reduce(function(a,b) merge(a,b, by=c("chr","start"), all=TRUE), pbulk)
write.csv(pb, "scRNA_pseudobulk_bin_profiles.csv", row.names=FALSE)
cat("\nDONE\n")
