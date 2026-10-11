###############################################################################
# R2_Q78 -- cHub links and peak-to-gene links together, peak-to-gene at the manuscript distance (250 kb)
#
# Question: within 250 kb, does ANY peak-to-gene link reach HAND2 or MEF2C, and do the links of
# cHub2711 / cHub3727 / cHub3728 connect to those genes?
# Peak-to-gene: ArchR addPeak2GeneLinks, GeneScoreMatrix + IterativeLSI (as scatac_main_hub_analysis.R).
# The pairs were computed once out to 500 kb (R2_Q78_peak2gene.R); the correlation of a pair does not
# depend on the distance limit, so the 250 kb set is the subset with |peak centre - TSS| <= 250 kb and
# the FDR is recomputed (BH) over that subset exactly as ArchR does it (t statistic from the correlation
# over the 500 cell aggregates). A link = correlation >= 0.2 and FDR < 1e-4 (manuscript cut-off 0.2).
#
# Figure: three panels side by side, one per hub, each window = hub +/- 250 kb (widened to 30 kb past the
# TSS of HAND2 / MEF2C where that lies further away), with
#   cHub peaks | cHub links (pairwise links between the hub's peaks) | peak-to-gene links | genes.
# Peak-to-gene links drawn: every link from a peak of the hub, plus every link to HAND2 / MEF2C that
# lies inside the window. MEF2C's TSS is 364-399 kb from the two MEF2C hubs, i.e. outside both the
# 250 kb limit and the window, so no link from those hubs to MEF2C can exist at this setting.
#
# Output: R2_Q78_p2g250_links_to_HAND2_MEF2C.csv, R2_Q78_p2g250_hub_peak_links.csv,
#         Plots/R2_Q78_cHubs_chub_and_p2g_links_250kb.pdf
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges); library(ggplot2); library(ggrepel); library(grid); library(gtable)
                                 library(RColorBrewer); library(paletteer); library(circlize) })
addArchRThreads(8); addArchRGenome("hg38")
ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT    <- file.path(ROOT, "git_repo_claude", "R2_Q78"); SCR <- "/sc/arion/scratch/giottb01/R2_Q78_peak2gene"
HUBDIR <- file.path(ROOT, "main", "scatac_ArchR", "hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5")
MAXDIST <- 250000; COR <- 0.2; FDRC <- 1e-4; NAGG <- 500
setwd(OUT); source(file.path(ROOT, "git_repo", "utils", "scATAC_functions.R")); source(file.path(ROOT, "git_repo", "utils", "palettes.R"))
HUBS <- c("HUB2711", "HUB3727", "HUB3728"); GOI <- c("HAND2", "MEF2C")
hubs <- readRDS(file.path(HUBDIR, "global_hubs_obj.rds")); p2g <- readRDS(file.path(SCR, "peak2gene_links_500kb.rds"))
pk <- metadata(p2g)$peakSet; gn <- metadata(p2g)$geneSet; L <- as.data.frame(p2g)
L$gene <- gn$name[L$idxRNA]; L$chr <- as.character(seqnames(pk))[L$idxATAC]; L$pmid <- (start(pk)[L$idxATAC] + end(pk)[L$idxATAC]) %/% 2; L$tss <- start(gn)[L$idxRNA]; L$dist <- abs(L$pmid - L$tss)
L <- L[L$dist <= MAXDIST & is.finite(L$Correlation), ]
tt <- L$Correlation / sqrt(pmax(1 - L$Correlation^2, 1e-12) / (NAGG - 2)); L$FDR250 <- p.adjust(2 * pt(-abs(tt), NAGG - 2), "fdr"); L$link <- L$Correlation >= COR & L$FDR250 < FDRC
cat(sprintf("pairs within %d kb: %d ; links (r >= %.1f, FDR < %g): %d (%.1f%%)\n", MAXDIST / 1e3, nrow(L), COR, FDRC, sum(L$link), 100 * mean(L$link)))

hp <- do.call(rbind, lapply(HUBS, function(h) { d <- hubs$hubsMerged[[match(h, hubs$hubs_id)]]; data.frame(hub = h, chr = as.character(d$seqnames), start = d$start, end = d$end) }))
ov <- findOverlaps(pk, GRanges(hp$chr, IRanges(hp$start, hp$end))); hubOfPeak <- rep(NA_character_, length(pk)); hubOfPeak[queryHits(ov)] <- hp$hub[subjectHits(ov)]; L$hub <- hubOfPeak[L$idxATAC]

## ---- 1. any peak-to-gene link to HAND2 / MEF2C within 250 kb? ---------------------------
G <- L[L$gene %in% GOI, ]; G <- G[order(G$gene, -G$Correlation), ]
write.csv(G[, c("gene", "chr", "pmid", "tss", "dist", "Correlation", "FDR250", "link", "hub")], "R2_Q78_p2g250_links_to_HAND2_MEF2C.csv", row.names = FALSE)
options(width = 220); cat("\n=== peak-to-gene links to HAND2 and MEF2C within 250 kb ===\n")
for (g in GOI) { x <- G[G$gene == g, ]; lk <- x[x$link, ]
  cat(sprintf("%s (TSS %s:%s): %d peaks tested, %d linked; correlation of the links %.2f-%.2f (median %.2f); farthest link %.0f kb; links from the three hubs' peaks: %d\n", g, x$chr[1], format(x$tss[1], big.mark = ","), nrow(x), nrow(lk),
              min(lk$Correlation), max(lk$Correlation), median(lk$Correlation), max(lk$dist) / 1e3, sum(!is.na(lk$hub))))
  cat("  top 5 links:\n"); print(format(head(lk[, c("chr", "pmid", "dist", "Correlation", "FDR250", "hub")], 5), digits = 3), row.names = FALSE) }
cat("\nhub peaks tested against the genes of interest at 250 kb:\n"); print(table(hub = L$hub[!is.na(L$hub)], gene_of_interest = ifelse(L$gene[!is.na(L$hub)] %in% GOI, L$gene[!is.na(L$hub)], "other")))
H <- L[!is.na(L$hub), ]; H <- H[order(H$hub, H$gene, -H$Correlation), ]; write.csv(H[, c("hub", "chr", "pmid", "gene", "tss", "dist", "Correlation", "FDR250", "link")], "R2_Q78_p2g250_hub_peak_links.csv", row.names = FALSE)
cat("\n=== genes linked to each hub's peaks at 250 kb (n peaks linked / n tested, max r) ===\n")
for (h in HUBS) { x <- H[H$hub == h, ]; if (!nrow(x)) { cat(h, ": no gene TSS within 250 kb of its peaks -> nothing tested\n"); next }
  sp <- split(x, x$gene); s <- sapply(sp, function(d) sprintf("%s %d/%d (%.2f)", d$gene[1], sum(d$link), nrow(d), max(d$Correlation))); cat(h, ":", paste(s[order(-sapply(sp, function(d) max(d$Correlation)))], collapse = "; "), "\n") }

## ---- 2. tracks: hub +/- 250 kb, cHub links and peak-to-gene links together -----------------
proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE); proj@projectMetadata$outputDirectory <- SCR
celltype_order <- c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs')
CT_PAL <- palette_celltype_lv1; names(CT_PAL)[names(CT_PAL) == "Smooth Muscle"] <- "SmoothMuscle"; CT_PAL <- CT_PAL[celltype_order]
## window = hub +/- 250 kb, widened where needed to reach 30 kb past the TSS of the gene of interest, so that links to that gene are in view
EXPECT <- c(HUB2711 = "HAND2", HUB3727 = "MEF2C", HUB3728 = "MEF2C")
reg <- do.call(c, lapply(HUBS, function(h) { b <- hp[hp$hub == h, ]; t0 <- start(gn)[match(EXPECT[[h]], gn$name)]
  r <- GRanges(b$chr[1], IRanges(min(min(b$start) - MAXDIST, t0 - 30000), max(max(b$end) + MAXDIST, t0 + 30000))); r$label <- sprintf("%s and %s", sub("HUB", "cHub", h), EXPECT[[h]]); r }))
hub_pk <- GRanges(hp$chr, IRanges(hp$start, hp$end)); hubWin <- reduce(resize(hub_pk, width(hub_pk) + 2000, fix = "center"), min.gapwidth = 15000)
chub <- subsetByOverlaps(hubs$peakLinks[[1]], hubWin, type = "within")
S <- L[L$link & (!is.na(L$hub) | L$gene %in% GOI), ]
p2gr <- GRanges(S$chr, IRanges(pmin(S$pmid, S$tss), pmax(S$pmid, S$tss)), value = S$Correlation); p2gr <- subsetByOverlaps(p2gr, reg, type = "within")
cat(sprintf("\nloops drawn: %d cHub links, %d peak-to-gene links (from hub peaks: %d; to HAND2 / MEF2C: %d)\n", length(chub), length(p2gr), sum(S$link & !is.na(S$hub)), sum(S$gene %in% GOI)))
plotTracksSideBySide(proj, reg, groupBy = "celltype_lv1", groups = celltype_order, pal = CT_PAL, features = GRangesList(`cHub peaks` = hub_pk),
                     loops = GRangesList(`cHub links` = chub, `peak2gene (250 kb)` = p2gr), file = "Plots/R2_Q78_cHubs_chub_and_p2g_links_250kb.pdf",
                     padding = 0.02, panelWidth = 4.0, tileSize = 500, outlineProfiles = FALSE, sizes = c(10, 0.45, 2.2, 1.9))
cat("\nDONE\n")
