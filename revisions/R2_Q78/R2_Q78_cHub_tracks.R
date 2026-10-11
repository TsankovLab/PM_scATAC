###############################################################################
# R2_Q78 -- browser tracks of cHub2711, cHub3722 and cHub3728, side by side
#
# The reviewer asks to see three cHubs: 2711 (expected near HAND2) and 3722 / 3728 (expected
# near MEF2C). Hubs are those of the main analysis:
#   main/scatac_ArchR/hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5  (6,719 hubs)
# Tracks: main ArchR project, grouped by celltype_lv1, drawn with plotTracksSideBySide()
# (git_repo/utils/scATAC_functions.R): coverage per cell type, the hub's peaks, the pairwise
# peak links of the hub, and the gene track.
#
# Before plotting, the script checks WHERE each hub is relative to the gene it is supposed to
# be near, because the hub ids are easy to mistype. In this hub set:
#   HUB2711  chr4:173.41 Mb  nearest gene SCRG1, HAND2 ~0.1 Mb away      -> consistent
#   HUB3728  chr5:89.29 Mb   annotated to MEF2C                           -> consistent
#   HUB3722  chr5:74.25 Mb   LINC01335, ~14.5 Mb from MEF2C               -> NOT near MEF2C
#   HUB3727  chr5:89.27 Mb   annotated to MEF2C, next to HUB3728          -> the likely intended hub
# So two versions are written: the three ids exactly as requested, and the same figure with
# HUB3727 in place of HUB3722. A third figure widens each window to include the gene, with
# the hub peaks marked, so the position of the hub relative to HAND2 / MEF2C is visible.
#
# Output (Plots/):
#   R2_Q78_cHubs_2711_3722_3728.pdf           as requested
#   R2_Q78_cHubs_2711_3727_3728.pdf           with the two MEF2C hubs
#   R2_Q78_cHubs_gene_context.pdf             wide windows: hub + HAND2 / MEF2C
#   R2_Q78_hub_positions.csv                  coordinates, peaks, nearest genes, distance to HAND2 / MEF2C
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges); library(ggplot2); library(ggrepel); library(grid); library(gtable)
                                 library(RColorBrewer); library(paletteer); library(circlize) })
addArchRThreads(8); addArchRGenome("hg38")
ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT    <- file.path(ROOT, "git_repo_claude", "R2_Q78")
HUBDIR <- file.path(ROOT, "main", "scatac_ArchR", "hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5")
setwd(OUT); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "scATAC_functions.R"))
source(file.path(ROOT, "git_repo", "utils", "palettes.R"))

REQUESTED <- c("HUB2711", "HUB3722", "HUB3728")
CORRECTED <- c("HUB2711", "HUB3727", "HUB3728")
EXPECT    <- c(HUB2711 = "HAND2", HUB3722 = "MEF2C", HUB3727 = "MEF2C", HUB3728 = "MEF2C")
celltype_order <- c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial',
                    'Myeloid','T_cells','NK','B_cells','Plasma','pDCs')
## palettes.R names this cell type "Smooth Muscle"; the project metadata uses "SmoothMuscle"
CT_PAL <- palette_celltype_lv1; names(CT_PAL)[names(CT_PAL) == "Smooth Muscle"] <- "SmoothMuscle"
stopifnot(all(celltype_order %in% names(CT_PAL))); CT_PAL <- CT_PAL[celltype_order]

hubs <- readRDS(file.path(HUBDIR, "global_hubs_obj.rds"))
bed  <- read.table(file.path(HUBDIR, "hub_regions.bed"), sep = "\t", stringsAsFactors = FALSE)
colnames(bed)[1:6] <- c("chr", "start", "end", "id", "strand", "genes")
proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE)

## ---- where is each hub relative to the gene it should be near? ---------------------------
genes <- getGenes(proj)
ids   <- union(REQUESTED, CORRECTED)
pos <- do.call(rbind, lapply(ids, function(h) { b <- bed[bed$id == h, ]; g <- genes[which(genes$symbol == EXPECT[[h]])][1]
  hg <- GRanges(b$chr, IRanges(b$start, b$end)); same <- as.character(seqnames(g)) == b$chr
  near <- genes[nearest(hg, genes)]
  data.frame(hub = h, chr = b$chr, start = b$start, end = b$end, width_kb = round((b$end - b$start) / 1e3, 1),
             n_peaks = nrow(hubs$hubsMerged[[match(h, hubs$hubs_id)]]), annotated_genes = b$genes, nearest_gene = near$symbol,
             expected_gene = EXPECT[[h]], gene_chr = as.character(seqnames(g)), gene_start = start(g), gene_end = end(g),
             distance_to_expected_kb = if (same) round(distance(hg, g) / 1e3, 1) else NA) }))
write.csv(pos, "R2_Q78_hub_positions.csv", row.names = FALSE)
options(width = 220); cat("=== hub positions relative to the expected gene ===\n"); print(pos, row.names = FALSE)

## ---- 1 + 2: the hubs themselves -----------------------------------------------------------
hdr <- function(i, regions) sprintf("c%s  %s\n%s:%s-%s", sub("HUB", "Hub", regions$label[i]), ifelse(is.null(regions$genes), "", regions$genes[i]),
                                    as.character(seqnames(regions))[i], format(start(regions)[i], big.mark = ","), format(end(regions)[i], big.mark = ","))
for (set in list(list(ids = REQUESTED, file = "Plots/R2_Q78_cHubs_2711_3722_3728.pdf"),
                 list(ids = CORRECTED, file = "Plots/R2_Q78_cHubs_2711_3727_3728.pdf"))) {
  message("plotting ", paste(set$ids, collapse = ", "))
  plotTracksSideBySide(proj, set$ids, hubs = hubs, hub_bed = bed, groupBy = "celltype_lv1", groups = celltype_order,
                       pal = CT_PAL, file = set$file, padding = 0.25, panelWidth = 2.6)
}

## ---- 3: gene context -- window from the hub to the expected gene ---------------------------
ctx_ids <- c("HUB2711", "HUB3727", "HUB3728")
ctx <- do.call(c, lapply(c("HAND2", "MEF2C"), function(gn) { h <- names(EXPECT)[EXPECT == gn & names(EXPECT) %in% ctx_ids]; b <- bed[bed$id %in% h, ]; g <- genes[which(genes$symbol == gn)][1]
  r <- GRanges(b$chr[1], IRanges(min(c(b$start, start(g))), max(c(b$end, end(g)))))
  r$label <- sprintf("%s and %s", paste(sub("HUB", "cHub", h), collapse = " / "), gn); r }))
hub_pk <- do.call(c, lapply(ctx_ids, function(h) { d <- hubs$hubsMerged[[match(h, hubs$hubs_id)]]; GRanges(as.character(d$seqnames), IRanges(d$start, d$end)) }))
plotTracksSideBySide(proj, ctx, groupBy = "celltype_lv1", groups = celltype_order, pal = CT_PAL,
                     features = GRangesList(`cHub peaks` = hub_pk),
                     ## only the links of these three hubs: the window also contains other hubs (e.g. over the MEF2C gene body)
                     loops = GRangesList(`cHub links` = subsetByOverlaps(hubs$peakLinks[[1]], reduce(resize(hub_pk, width(hub_pk) + 2000, fix = "center"), min.gapwidth = 15000), type = "within")),
                     ## no profile outline here: at 0.1-0.7 Mb per panel the outline would cover the cell-type colour
                     file = "Plots/R2_Q78_cHubs_gene_context.pdf", padding = 0.08, panelWidth = 4.2, tileSize = 500, outlineProfiles = FALSE)
cat("\nDONE\n")
