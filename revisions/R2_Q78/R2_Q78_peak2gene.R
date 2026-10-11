###############################################################################
# R2_Q78 -- are cHub2711 and cHub3727 / cHub3728 linked to HAND2 and MEF2C by peak-to-gene correlation?
#
# ArchR addPeak2GeneLinks on the main project, set up as in the manuscript (git_repo/main_analysis/
# scatac_main_hub_analysis.R): useMatrix = GeneScoreMatrix, reducedDims = IterativeLSI. The project has
# no integrated RNA, so "gene" here is the gene activity score. One change: maxDist is raised from
# 250 kb to 500 kb, because the MEF2C hubs lie 363-387 kb from the MEF2C TSS and would not be tested
# at all with the manuscript setting. cHub2711 is ~115 kb from the HAND2 TSS.
# The gene score uses accessibility in the gene body and up to 100 kb around it, so these hubs do not
# contribute to the score they are correlated with (no built-in circularity at these distances).
#
# The main project is not modified: the output directory is redirected in memory to scratch, where
# addPeak2GeneLinks writes its (multi-GB) aggregate matrices.
#
# Reported for every peak of the four hubs (2711, 3722, 3727, 3728): correlation and FDR with the
# expected gene, its rank among all genes tested for that peak, and an empirical p against all
# peak-gene pairs genome-wide at a comparable distance. HUB3722 is kept as a negative example: it is
# 14.5 Mb from MEF2C and cannot be tested against it.
#
# Output: R2_Q78_peak2gene_hub_peaks.csv, R2_Q78_peak2gene_all_genes_for_hub_peaks.csv,
#         R2_Q78_peak2gene_summary.csv, Plots/R2_Q78_cHubs_peak2gene.pdf
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges); library(ggplot2); library(ggrepel); library(grid); library(gtable)
                                 library(RColorBrewer); library(paletteer); library(circlize) })
addArchRThreads(8); addArchRGenome("hg38"); set.seed(1)
ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT    <- file.path(ROOT, "git_repo_claude", "R2_Q78")
SCR    <- "/sc/arion/scratch/giottb01/R2_Q78_peak2gene"
HUBDIR <- file.path(ROOT, "main", "scatac_ArchR", "hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5")
MAXDIST <- 500000
setwd(OUT); dir.create("Plots", showWarnings = FALSE); dir.create(SCR, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "git_repo", "utils", "scATAC_functions.R")); source(file.path(ROOT, "git_repo", "utils", "palettes.R"))
EXPECT <- c(HUB2711 = "HAND2", HUB3722 = "MEF2C", HUB3727 = "MEF2C", HUB3728 = "MEF2C")
hubs <- readRDS(file.path(HUBDIR, "global_hubs_obj.rds"))
proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE)
proj@projectMetadata$outputDirectory <- SCR                      # in memory only: nothing is written to main/scatac_ArchR

P2G_RDS <- file.path(SCR, "peak2gene_links_500kb.rds")
if (file.exists(P2G_RDS)) { p2g <- readRDS(P2G_RDS) } else {
  proj <- addPeak2GeneLinks(proj, useMatrix = "GeneScoreMatrix", reducedDims = "IterativeLSI", maxDist = MAXDIST)
  p2g <- metadata(proj@peakSet)$Peak2GeneLinks; saveRDS(p2g, P2G_RDS) }
pk <- metadata(p2g)$peakSet; gn <- metadata(p2g)$geneSet
L <- as.data.frame(p2g); L$gene <- gn$name[L$idxRNA]; L$peak <- sprintf("%s:%d-%d", as.character(seqnames(pk))[L$idxATAC], start(pk)[L$idxATAC], end(pk)[L$idxATAC])
L$dist <- abs((start(pk)[L$idxATAC] + end(pk)[L$idxATAC]) / 2 - start(gn)[L$idxRNA])
cat(sprintf("peak-gene pairs tested within %d kb: %d ; with correlation >= 0.2 and FDR < 1e-4: %d (%.1f%%)\n", MAXDIST / 1e3, nrow(L), sum(L$Correlation >= 0.2 & L$FDR < 1e-4, na.rm = TRUE), 100 * mean(L$Correlation >= 0.2 & L$FDR < 1e-4, na.rm = TRUE)))

## peaks of the four hubs
hp <- do.call(rbind, lapply(names(EXPECT), function(h) { d <- hubs$hubsMerged[[match(h, hubs$hubs_id)]]; data.frame(hub = h, chr = as.character(d$seqnames), start = d$start, end = d$end) }))
ov <- findOverlaps(GRanges(hp$chr, IRanges(hp$start, hp$end)), pk); hp$idxATAC <- NA; hp$idxATAC[queryHits(ov)] <- subjectHits(ov)
A <- merge(hp, L, by = "idxATAC"); A <- A[order(A$hub, A$start, -A$Correlation), ]
write.csv(A[, c("hub", "peak", "gene", "dist", "Correlation", "FDR", "VarQATAC", "VarQRNA")], "R2_Q78_peak2gene_all_genes_for_hub_peaks.csv", row.names = FALSE)

res <- do.call(rbind, lapply(seq_len(nrow(hp)), function(i) { h <- hp$hub[i]; g <- EXPECT[[h]]; a <- A[A$idxATAC == hp$idxATAC[i] & A$hub == h, ]; x <- a[a$gene == g, ]
  if (!nrow(x)) return(data.frame(hub = h, peak = sprintf("%s:%d-%d", hp$chr[i], hp$start[i], hp$end[i]), expected_gene = g, tested = FALSE, dist_kb = NA, Correlation = NA, FDR = NA, rank_among_genes = NA, n_genes_tested = nrow(a), top_gene = if (nrow(a)) a$gene[1] else NA, top_cor = if (nrow(a)) a$Correlation[1] else NA, emp_p = NA))
  bg <- L$Correlation[abs(L$dist - x$dist) < 50000]; bg <- bg[is.finite(bg)]
  data.frame(hub = h, peak = x$peak, expected_gene = g, tested = TRUE, dist_kb = round(x$dist / 1e3), Correlation = x$Correlation, FDR = x$FDR, rank_among_genes = sum(a$Correlation > x$Correlation, na.rm = TRUE) + 1, n_genes_tested = nrow(a),
             top_gene = a$gene[1], top_cor = a$Correlation[1], emp_p = mean(bg >= x$Correlation)) }))
write.csv(res, "R2_Q78_peak2gene_hub_peaks.csv", row.names = FALSE)
options(width = 230); cat("\n=== every hub peak vs the expected gene (GeneScore-based peak2gene) ===\n"); print(format(res, digits = 3), row.names = FALSE)
sm <- do.call(rbind, lapply(names(EXPECT), function(h) { r <- res[res$hub == h, ]; data.frame(hub = h, expected_gene = EXPECT[[h]], n_peaks = nrow(r), n_tested = sum(r$tested), n_linked_cor0.2_FDR1e4 = sum(r$Correlation >= 0.2 & r$FDR < 1e-4, na.rm = TRUE),
  n_linked_cor0.45 = sum(r$Correlation >= 0.45 & r$FDR < 1e-4, na.rm = TRUE), max_cor = suppressWarnings(max(r$Correlation, na.rm = TRUE)), median_cor = median(r$Correlation, na.rm = TRUE), n_expected_is_top_gene = sum(r$rank_among_genes == 1, na.rm = TRUE)) }))
write.csv(sm, "R2_Q78_peak2gene_summary.csv", row.names = FALSE); cat("\n=== per hub ===\n"); print(format(sm, digits = 3), row.names = FALSE)
cat("\n=== strongest gene links of each hub (any gene; correlation >= 0.2, FDR < 1e-4), counted over its peaks ===\n")
for (h in names(EXPECT)) { a <- A[A$hub == h & A$Correlation >= 0.2 & A$FDR < 1e-4 & !is.na(A$Correlation), ]; tb <- sort(table(a$gene), decreasing = TRUE)
  cat(sprintf("%s: %s\n", h, if (length(tb)) paste(sprintf("%s (%d peaks, max r %.2f)", names(tb), as.integer(tb), sapply(names(tb), function(g) max(a$Correlation[a$gene == g]))), collapse = "; ") else "none")) }

## figure: gene-context windows with the peak2gene links of the hub peaks drawn as loops
genes <- getGenes(proj); celltype_order <- c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs')
CT_PAL <- palette_celltype_lv1; names(CT_PAL)[names(CT_PAL) == "Smooth Muscle"] <- "SmoothMuscle"; CT_PAL <- CT_PAL[celltype_order]
show <- c("HUB2711", "HUB3727", "HUB3728"); a <- A[A$hub %in% show & A$Correlation >= 0.2 & A$FDR < 1e-4 & !is.na(A$Correlation), ]
loopGR <- GRanges(as.character(seqnames(pk))[a$idxATAC], IRanges(pmin((start(pk)[a$idxATAC] + end(pk)[a$idxATAC]) %/% 2, start(gn)[a$idxRNA]), pmax((start(pk)[a$idxATAC] + end(pk)[a$idxATAC]) %/% 2, start(gn)[a$idxRNA])), value = a$Correlation)
ctx <- do.call(c, lapply(c("HAND2", "MEF2C"), function(g0) { h <- intersect(names(EXPECT)[EXPECT == g0], show); b <- hp[hp$hub %in% h, ]; g <- genes[which(genes$symbol == g0)][1]; lr <- loopGR[as.character(seqnames(loopGR)) == b$chr[1]]
  r <- GRanges(b$chr[1], IRanges(min(c(b$start, start(g), start(lr))), max(c(b$end, end(g), end(lr))))); r$label <- sprintf("%s: peak-to-gene links", paste(sub("HUB", "cHub", h), collapse = " / ")); r }))
hub_pk <- GRanges(hp$chr[hp$hub %in% show], IRanges(hp$start[hp$hub %in% show], hp$end[hp$hub %in% show]))
plotTracksSideBySide(proj, ctx, groupBy = "celltype_lv1", groups = celltype_order, pal = CT_PAL, features = GRangesList(`cHub peaks` = hub_pk), loops = GRangesList(`peak2gene` = loopGR),
                     file = "Plots/R2_Q78_cHubs_peak2gene.pdf", padding = 0.05, panelWidth = 4.4, tileSize = 500, outlineProfiles = FALSE)
cat("\nDONE\n")
