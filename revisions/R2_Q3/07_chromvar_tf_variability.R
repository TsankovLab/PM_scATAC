###############################################################################
# STEP 7 -- TF activity between the epiAneufinder subclones.
#
# Question: do the CNV clones correspond to distinct regulatory states, and is there a
# TF that moves between clones in tumour after tumour?
#
# Only the tumours that actually split (step 3) can be tested; the ones the ARM_MIN
# criterion left as a single clone are skipped, and the "across samples" summary is a
# median over the tumours that remain.
#
#   chromVAR z  per-cell deviation z-score for each of 869 cisBP motifs, taken from the
#               MotifMatrix of the tumour-compartment ArchR project. Pulling it out is
#               the slow step, so it is cached to chromvar_z_cache.rds (~128 MB).
#   per sample  |difference| of mean z between the two clones -- the variability of that
#               motif across the clones, deliberately UNSIGNED: which clone is "c1" is an
#               arbitrary label, so a direction would not mean anything,
#               eta^2 (fraction of that motif's variance the clone split explains),
#               Wilcoxon + BH, and the ARI of an INDEPENDENT chromVAR-only k = 2
#               clustering against the CNV clones -- does TF activity recover the same
#               partition when it is not told about it?
#   across      MEDIAN |difference| with the IQR. Median, not mean: samples differ
#               several-fold in how far apart their clones sit, so one tumour must not
#               set the ranking.
#   gene score  the SAME contrast on a second ATAC readout: mean ArchR GeneScoreMatrix
#               value of each motif's TF gene, log2(score + 1), computed PER CLONE, so
#               that |gene score c2 - c1| is the gene-level twin of |dz|. Cached per
#               (tumour, clone, gene) to archr_tf_genescore.csv, with the gene's
#               coordinates so a gene sitting on the arm that defines the split can be
#               flagged. CAVEATS: gene score comes from the same fragments as the
#               chromVAR z, AND a clone defined by a copy-number gain has more fragments
#               across that whole arm, so genes there gain score for reasons that have
#               nothing to do with regulation. Both are consistency checks, not
#               independent evidence.
#
#   flags       AP-1/bZIP + SMARCC1 are marked "technical": that axis tracks per-cell
#               fragment count in this dataset and floors in every comparison, so it is
#               not clonal biology. Immune/interferon motifs are flagged as a class and
#               tested for enrichment among the most variable motifs.
#
# NOTE the FDR column ranks motifs WITHIN a sample; it scales with cell count and must
#      not be compared between samples. Use `diff` and `eta2` for that.
#
# Input : epi_clone_labels.csv (step 3); ArchR project (first run only)
# Output: chromvar_z_cache.rds, archr_tf_genescore.csv,
#         archr_tf_genescore_clonediff.csv, epi_chromvar_diff_<S>.csv,
#         epi_chromvar_mirror.csv, epi_chromvar_recurrence.csv
###############################################################################
suppressMessages({ library(data.table) })
source("00_common.R")
ZCACHE <- "chromvar_z_cache.rds"

## AP-1/bZIP accessibility axis -- technical, see note above
TECH <- c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","BACH1","BACH2","NFE2",
          "SMARCC1","JDP2","BATF","BATF3","NFE2L1","NFE2L2","MAF","MAFF","MAFK","MAFG")
IMMUNE <- c(paste0("IRF", 1:9), paste0("STAT", c(1,2,3,4,"5A","5B",6)),
            "NFKB1","NFKB2","REL","RELA","RELB","SPI1","SPIB","SPIC",
            paste0("IKZF", 1:4), paste0("NFATC", 1:4), "TBX21","EOMES","PRDM1",
            "RUNX3","ETS1","ELF4","BCL11B","TCF7","LEF1","FOXP3","POU2F2","POU2AF1",
            "CIITA","NLRC5")

## ---- chromVAR z (cached) ----------------------------------------------------
## Row names carry a numeric suffix (TFAP2B_1); strip it to the TF symbol.
if (file.exists(ZCACHE)){
  Z <- readRDS(ZCACHE); cat("chromVAR z from cache\n")
} else {
  suppressMessages(library(ArchR)); addArchRThreads(1)
  proj <- loadArchRProject(ARCHR, showLogo = FALSE)
  se <- getMatrixFromProject(proj, useMatrix = "MotifMatrix", verbose = FALSE)
  Z <- as.matrix(SummarizedExperiment::assays(se)$z)
  rownames(Z) <- make.unique(sub("_[0-9]+$", "", SummarizedExperiment::rowData(se)$name))
  saveRDS(Z, ZCACHE)
}
cat("chromVAR z:", nrow(Z), "motifs x", ncol(Z), "cells\n")

LAB <- fread("epi_clone_labels.csv")
cat("clone labels:", nrow(LAB), "cells\n")

## a tumour that no longer splits must not leave its old per-sample table behind
invisible(file.remove(Sys.glob("epi_chromvar_diff_*.csv")))

## ---- per sample -------------------------------------------------------------
diffs <- list(); metr <- list()
for (S in SAMPLES){
  l <- LAB[sample == S]; v <- setNames(l$clone, l$cell)
  if (length(unique(v)) < 2){ cat("skip", S, "- one clone (primary split below ARM_MIN)\n"); next }
  cells <- intersect(names(v), colnames(Z))
  if (length(cells) < 40){ cat("skip", S, "- <40 cells in the motif matrix\n"); next }
  v <- v[cells]; zt <- Z[, cells, drop = FALSE]
  i1 <- which(v == "c1"); i2 <- which(v == "c2"); n1 <- length(i1); n2 <- length(i2)
  if (min(n1, n2) < 20){ cat("skip", S, "- a clone has <20 cells with motif data\n"); next }

  m1 <- rowMeans(zt[, i1, drop = FALSE]); m2 <- rowMeans(zt[, i2, drop = FALSE])
  gm <- rowMeans(zt)
  eta2 <- (n1*(m1-gm)^2 + n2*(m2-gm)^2) / rowSums((zt - gm)^2)
  pw  <- apply(zt, 1, function(x) tryCatch(wilcox.test(x[i1], x[i2])$p.value,
                                           error = function(e) NA))
  fdr <- p.adjust(pw, "BH")
  ## absdiff, not (m2 - m1): the clone labels carry no order, so the variability of a
  ## motif across the clones is a magnitude and nothing else
  dd  <- data.frame(sample = S, TF = rownames(zt), mean_c1 = m1, mean_c2 = m2,
                    absdiff = abs(m2 - m1), eta2 = eta2, p = pw, fdr = fdr, row.names = NULL)
  dd$abs_pct <- rank(dd$absdiff) / nrow(dd)            # within-sample percentile
  dd$tech    <- dd$TF %in% TECH
  diffs[[S]] <- dd
  write.csv(dd[order(dd$fdr), ], sprintf("epi_chromvar_diff_%s.csv", S), row.names = FALSE)

  ## independent chromVAR-only clustering, on the 100 most variable motifs
  tv <- order(apply(zt, 1, var), decreasing = TRUE)[seq_len(min(100, nrow(zt)))]
  km <- kmeans(scale(t(zt[tv, , drop = FALSE])), centers = 2, nstart = 10)$cluster
  bio <- dd[!dd$tech, ]
  metr[[S]] <- data.frame(sample = S, n_c1 = n1, n_c2 = n2,
    n_sig_fdr05 = sum(fdr < 0.05, na.rm = TRUE),
    mean_eta2 = mean(eta2, na.rm = TRUE), max_eta2 = max(eta2, na.rm = TRUE),
    ARI_chromvar_vs_CNV = ARI(km, v),
    top_TF_any = dd$TF[which.max(dd$absdiff)],
    top_TF_nontech = bio$TF[which.max(bio$absdiff)], stringsAsFactors = FALSE)
}
metr <- do.call(rbind, metr); rownames(metr) <- NULL
write.csv(metr, "epi_chromvar_mirror.csv", row.names = FALSE)
cat("\n=== does TF activity mirror the epiAneufinder clone split? ===\n")
print(transform(metr, mean_eta2 = round(mean_eta2, 4), max_eta2 = round(max_eta2, 3),
                ARI_chromvar_vs_CNV = round(ARI_chromvar_vs_CNV, 3)), row.names = FALSE)

## ---- across samples ---------------------------------------------------------
D <- do.call(rbind, diffs)
rec <- do.call(rbind, lapply(split(D, D$TF), function(x)
  data.frame(TF = x$TF[1], n_samples = nrow(x), n_sig = sum(x$fdr < 0.05, na.rm = TRUE),
             n_strong = sum(x$abs_pct >= 0.95, na.rm = TRUE),
             median_absdiff = median(x$absdiff),
             q25 = quantile(x$absdiff, .25, names = FALSE),
             q75 = quantile(x$absdiff, .75, names = FALSE),
             mean_absdiff = mean(x$absdiff),
             med_eta2 = median(x$eta2, na.rm = TRUE), tech = x$tech[1],
             immune = x$TF[1] %in% IMMUNE, stringsAsFactors = FALSE)))
rec <- rec[order(-rec$median_absdiff), ]
write.csv(rec, "epi_chromvar_recurrence.csv", row.names = FALSE)

cat("\n=== most variable motifs between epiAneufinder clones (median |diff|) ===\n")
print(head(rec[, c("TF","n_samples","n_sig","n_strong","median_absdiff","med_eta2","tech","immune")], 20),
      row.names = FALSE)
cat("\n--- excluding the AP-1/bZIP technical axis ---\n")
print(head(rec[!rec$tech, c("TF","n_sig","n_strong","median_absdiff","q25","q75","immune")], 20),
      row.names = FALSE)
for (q in c(0.90, 0.95)){
  hit <- rec$median_absdiff >= quantile(rec$median_absdiff, q)
  tb <- table(immune = rec$immune, top = hit); ft <- fisher.test(tb)
  cat(sprintf("\nimmune/interferon in the top %.0f%%: %d of %d (%.0f%%) vs %.0f%% background | OR %.2f, p = %.2g\n",
      100*(1-q), tb["TRUE","TRUE"], sum(hit), 100*tb["TRUE","TRUE"]/sum(hit),
      100*mean(rec$immune), ft$estimate, ft$p.value))
}
## ---- gene score of the same TFs, same cells, PER CLONE (cached) -------------
## ArchR gene score for each motif's TF, averaged within each clone of each tumour, so
## the between-clone difference can be compared with the chromVAR one.  Gene score is
## derived from the SAME fragments as the chromVAR deviations and inherits the clones'
## copy-number difference directly, so this is an internal consistency check, not
## independent evidence -- see the driver-arm flag in step 8.
GSCACHE <- "archr_tf_genescore.csv"
if (file.exists(GSCACHE)){
  GS <- read.csv(GSCACHE, stringsAsFactors = FALSE); cat("\ngene scores from cache\n")
} else {
  suppressMessages({ library(ArchR); library(Matrix) }); addArchRThreads(1)
  proj <- loadArchRProject(ARCHR, showLogo = FALSE)
  se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = FALSE)
  G  <- SummarizedExperiment::assays(se)[[1]]
  ## NB: do NOT assign rownames() on the dgCMatrix -- in this environment that leaves an
  ## object whose `[` no longer dispatches ("object of type 'S4' is not subsettable").
  ## Keep the names in plain vectors and index the sparse matrix by POSITION only.
  rn <- SummarizedExperiment::rowData(se)$name
  cn <- colnames(se)
  gi <- match(intersect(unique(rownames(Z)), rn), rn)   # first row per TF symbol
  ci <- which(cn %in% LAB$cell)                         # exactly the clone-labelled cells
  gene <- rn[gi]; cells <- cn[ci]
  rd   <- SummarizedExperiment::rowData(se)
  gchr <- as.character(rd$seqnames)[gi]; gmb <- rd$start[gi] / 1e6
  G   <- log2(as.matrix(G[gi, ci, drop = FALSE]) + 1)   # ArchR scales gene scores to 10k
  key <- paste(LAB$sample, LAB$clone)[match(cells, LAB$cell)]
  cat("\nlabelled cells with a gene score, per clone:\n"); print(table(key))
  GS <- do.call(rbind, lapply(sort(unique(key)), function(k){
    i <- which(key == k)
    data.frame(sample = sub(" .*", "", k), clone = sub(".* ", "", k),
               gene = gene, chr = gchr, mb = gmb, n_cells = length(i),
               mean_score = rowMeans(G[, i, drop = FALSE]),
               pct_score  = rowMeans(G[, i, drop = FALSE] > 0),
               row.names = NULL) }))
  write.csv(GS, GSCACHE, row.names = FALSE)
  rm(se, G, proj); gc(FALSE)
}
cat(sprintf("gene scores: %d TF genes x %d clones in %d tumours\n",
            length(unique(GS$gene)), nrow(unique(GS[, c("sample","clone")])),
            length(unique(GS$sample))))
cat("all tumours that split are covered:", all(metr$sample %in% GS$sample), "\n")

## the gene-level twin of the chromVAR contrast, for the tumours with two clones
DEC <- read.csv("epi_clone_decision.csv", stringsAsFactors = FALSE)
GD <- do.call(rbind, lapply(metr$sample, function(S){
  x  <- GS[GS$sample == S, ]
  c1 <- x[x$clone == "c1", ]; c2 <- x[x$clone == "c2", ][match(c1$gene, x$gene[x$clone == "c2"]), ]
  data.frame(sample = S, gene = c1$gene,
             score_c1 = c1$mean_score, score_c2 = c2$mean_score,
             absdiff_gs = abs(c2$mean_score - c1$mean_score),
             mean_score = (c1$mean_score + c2$mean_score) / 2,
             arm = armof(c1$chr, c1$mb),
             on_split_arm = armof(c1$chr, c1$mb) == DEC$arm[DEC$sample == S],
             row.names = NULL) }))
write.csv(GD, "archr_tf_genescore_clonediff.csv", row.names = FALSE)
cat("\n=== |gene score c2 - c1| per tumour ===\n")
print(do.call(rbind, lapply(split(GD, GD$sample), function(x) data.frame(
  sample = x$sample[1], n_genes = nrow(x), median = round(median(x$absdiff_gs), 3),
  max = round(max(x$absdiff_gs), 3), top = x$gene[which.max(x$absdiff_gs)],
  n_on_split_arm = sum(x$on_split_arm),
  median_on_arm = round(median(x$absdiff_gs[x$on_split_arm]), 3),
  median_off_arm = round(median(x$absdiff_gs[!x$on_split_arm]), 3))))[
    match(metr$sample, sort(metr$sample)), ], row.names = FALSE)

cat("\nDONE\n")
