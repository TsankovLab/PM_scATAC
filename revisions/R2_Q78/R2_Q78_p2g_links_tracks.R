###############################################################################
# R2_Q78 -- peak-to-gene links of the three cHubs, with the cHub regions as an annotation layer
#
# Two loci side by side: cHub2711 / HAND2, and cHub3727 + cHub3728 / MEF2C (the two MEF2C hubs are 11 kb apart, so one panel
# shows both). Tracks, top to bottom:
#   coverage per cell type
#   cHub annotation   every cHub of the hub set that falls in the window, as blocks on ONE row, each labelled with its id number
#                     (ArchR's feature track has no labels, so .featureTracks is wrapped for the call and restored afterwards)
#   coverage profiles are outlined in black
#   peak-to-gene      links only (the pairwise cHub links are not drawn: the hub is shown by its region)
#   genes
# Peak-to-gene distance: widened from the manuscript's 250 kb to 400 kb, the smallest round limit at which the
# two MEF2C hubs (364-399 kb from the MEF2C TSS) can be tested against MEF2C at all. Same limit in all panels.
# Pairs come from the 500 kb run of R2_Q78_peak2gene.R (GeneScoreMatrix, IterativeLSI); the 400 kb set is the
# subset by distance with the FDR recomputed (BH) over it. Link = correlation >= 0.2, FDR < 1e-4.
# Links drawn: all links from a peak of a focal hub, and all links to HAND2 / MEF2C in the window; colour = correlation.
#
# Output: Plots/R2_Q78_cHubs_p2g_links.pdf, R2_Q78_p2g400_hub_peak_links.csv, R2_Q78_p2g400_links_to_HAND2_MEF2C.csv
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges); library(ggplot2); library(ggrepel); library(grid); library(gtable)
                                 library(RColorBrewer); library(paletteer); library(circlize) })
addArchRThreads(8); addArchRGenome("hg38")
ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT    <- file.path(ROOT, "git_repo_claude", "R2_Q78"); SCR <- "/sc/arion/scratch/giottb01/R2_Q78_peak2gene"
HUBDIR <- file.path(ROOT, "main", "scatac_ArchR", "hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5")
MAXDIST <- 400000; FLANK <- 250000; COR <- 0.2; FDRC <- 1e-4; NAGG <- 500
setwd(OUT); source(file.path(ROOT, "git_repo", "utils", "scATAC_functions.R")); source(file.path(ROOT, "git_repo", "utils", "palettes.R"))
EXPECT <- c(HUB2711 = "HAND2", HUB3727 = "MEF2C", HUB3728 = "MEF2C"); HUBS <- names(EXPECT); GOI <- unique(EXPECT)
hubs <- readRDS(file.path(HUBDIR, "global_hubs_obj.rds")); bed <- read.table(file.path(HUBDIR, "hub_regions.bed"), sep = "\t", stringsAsFactors = FALSE); colnames(bed)[1:6] <- c("chr", "start", "end", "id", "strand", "genes")
p2g <- readRDS(file.path(SCR, "peak2gene_links_500kb.rds")); pk <- metadata(p2g)$peakSet; gn <- metadata(p2g)$geneSet; L <- as.data.frame(p2g)
L$gene <- gn$name[L$idxRNA]; L$chr <- as.character(seqnames(pk))[L$idxATAC]; L$pmid <- (start(pk)[L$idxATAC] + end(pk)[L$idxATAC]) %/% 2; L$tss <- start(gn)[L$idxRNA]; L$dist <- abs(L$pmid - L$tss)
L <- L[L$dist <= MAXDIST & is.finite(L$Correlation), ]; tt <- L$Correlation / sqrt(pmax(1 - L$Correlation^2, 1e-12) / (NAGG - 2)); L$FDR400 <- p.adjust(2 * pt(-abs(tt), NAGG - 2), "fdr"); L$link <- L$Correlation >= COR & L$FDR400 < FDRC
hp <- do.call(rbind, lapply(HUBS, function(h) { d <- hubs$hubsMerged[[match(h, hubs$hubs_id)]]; data.frame(hub = h, chr = as.character(d$seqnames), start = d$start, end = d$end) }))
ov <- findOverlaps(pk, GRanges(hp$chr, IRanges(hp$start, hp$end))); hubOfPeak <- rep(NA_character_, length(pk)); hubOfPeak[queryHits(ov)] <- hp$hub[subjectHits(ov)]; L$hub <- hubOfPeak[L$idxATAC]
options(width = 220); cat(sprintf("pairs within %d kb: %d ; links: %d (%.1f%%)\n", MAXDIST / 1e3, nrow(L), sum(L$link), 100 * mean(L$link)))
H <- L[!is.na(L$hub), ]; H <- H[order(H$hub, H$gene, -H$Correlation), ]; write.csv(H[, c("hub", "chr", "pmid", "gene", "tss", "dist", "Correlation", "FDR400", "link")], "R2_Q78_p2g400_hub_peak_links.csv", row.names = FALSE)
G <- L[L$gene %in% GOI, ]; write.csv(G[order(G$gene, -G$Correlation), c("gene", "chr", "pmid", "tss", "dist", "Correlation", "FDR400", "link", "hub")], "R2_Q78_p2g400_links_to_HAND2_MEF2C.csv", row.names = FALSE)
cat("\n=== focal hub peaks vs the expected gene at 400 kb ===\n")
for (h in HUBS) { x <- H[H$hub == h & H$gene == EXPECT[[h]], ]; cat(sprintf("%s -> %s: %d peaks tested, %d linked; correlation %.2f to %.2f; distance %.0f-%.0f kb\n", h, EXPECT[[h]], nrow(x), sum(x$link), min(x$Correlation), max(x$Correlation), min(x$dist) / 1e3, max(x$dist) / 1e3)) }
cat("\ngenes linked to each hub (n peaks linked / tested, max r):\n")
for (h in HUBS) { x <- H[H$hub == h, ]; sp <- split(x, x$gene); s <- sapply(sp, function(d) sprintf("%s %d/%d (%.2f)", d$gene[1], sum(d$link), nrow(d), max(d$Correlation))); cat(" ", h, ":", paste(s[order(-sapply(sp, function(d) max(d$Correlation)))], collapse = "; "), "\n") }
for (g in GOI) { x <- G[G$gene == g, ]; cat(sprintf("%s: %d peaks tested within 400 kb, %d linked, of which %d are peaks of the focal hubs\n", g, nrow(x), sum(x$link), sum(x$link & !is.na(x$hub)))) }

## ---- windows, annotation layer, combined link layer -----------------------------------------
## one panel per locus: cHub3727 and cHub3728 are 11 kb apart, so they share the MEF2C panel
reg <- do.call(c, lapply(GOI, function(g) { h <- HUBS[EXPECT == g]; b <- hp[hp$hub %in% h, ]; t0 <- start(gn)[match(g, gn$name)]
  r <- GRanges(b$chr[1], IRanges(min(min(b$start) - FLANK, t0 - 30000), max(max(b$end) + FLANK, t0 + 30000))); r$label <- sprintf("%s and %s", paste(sub("HUB", "cHub", h), collapse = " / "), g); r }))
allhub <- GRanges(bed$chr, IRanges(bed$start, bed$end), id = bed$id); inwin <- subsetByOverlaps(allhub, reg)
cat(sprintf("\ncHubs in the two windows: %d (focal: %s)\n", length(unique(inwin$id)), paste(HUBS, collapse = ", "))); print(as.data.frame(inwin)[, c("seqnames", "start", "end", "id")], row.names = FALSE)
inwin$num <- sub("HUB", "", inwin$id); feat <- GRangesList(cHubs = inwin)
S <- L[L$link & (!is.na(L$hub) | L$gene %in% GOI), ]
pg <- subsetByOverlaps(GRanges(S$chr, IRanges(pmin(S$pmid, S$tss), pmax(S$pmid, S$tss)), value = S$Correlation), reg, type = "within")
cat(sprintf("peak-to-gene links drawn: %d (correlation %.2f-%.2f)\n", length(pg), min(pg$value), max(pg$value)))

proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE); proj@projectMetadata$outputDirectory <- SCR
celltype_order <- c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs')
CT_PAL <- palette_celltype_lv1; names(CT_PAL)[names(CT_PAL) == "Smooth Muscle"] <- "SmoothMuscle"; CT_PAL <- CT_PAL[celltype_order]
## label every cHub block with its id number: wrap ArchR's feature-track builder for this call only
ft_orig <- get(".featureTracks", envir = asNamespace("ArchR"))
ft_lab  <- function(features = NULL, region = NULL, ...) { p <- ft_orig(features = features, region = region, ...)
  f <- subsetByOverlaps(features[[1]], region, ignore.strand = TRUE); if (!length(f)) return(p)
  d <- data.frame(x = (pmax(start(f), start(region)) + pmin(end(f), end(region))) / 2, name = names(features)[1], lab = f$num)
  p + ggrepel::geom_text_repel(data = d, ggplot2::aes(x = x, y = name, label = lab), inherit.aes = FALSE, direction = "x", nudge_y = 0.42, size = 2.7, colour = "black",
                               segment.size = 0.2, min.segment.length = 0, box.padding = 0.08, point.padding = 0, max.overlaps = Inf, seed = 1) +
    ## blocks sit at the bottom of the row (labels above them), so nothing but the blocks separates this track from the links below
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.12, 0.9))) }
environment(ft_lab) <- environment(); assignInNamespace(".featureTracks", ft_lab, ns = "ArchR")
tryCatch(
  plotTracksSideBySide(proj, reg, groupBy = "celltype_lv1", groups = celltype_order, pal = CT_PAL, features = feat, loops = GRangesList(`peak2gene` = pg),
                       file = "Plots/R2_Q78_cHubs_p2g_links.pdf", padding = 0.02, panelWidth = 4.0, tileSize = 1000, sizes = c(10, 0.62, 1.6, 1.9),
                       outlineProfiles = TRUE, outlineWidth = 0.03),          # thin outline: at 0.5-0.7 Mb per panel a heavier one hides the fill
  finally = assignInNamespace(".featureTracks", ft_orig, ns = "ArchR"))
cat("\nDONE\n")
