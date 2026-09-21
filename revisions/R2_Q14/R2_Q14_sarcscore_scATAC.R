###############################################################################
# R2_Q14 -- sarcomatoid score computed FROM THE scATAC DATA, so P23 is included.
#
# Problem this fixes.  scATAC_BAP1_vs_sarcscore.pdf plotted 9 of the 11 scATAC
# tumours.  The y axis came from tumor_compartment/scrna/cnmf20_sarcomatoid_sample_
# order.csv, a scRNA-derived score, and that file has no row for P10 or P23:
#   * P23 was profiled by scATAC ONLY -- it has no scRNA at all, so it can never
#     have an scRNA-derived score.  It also carries the HIGHEST BAP1 gene score in
#     the cohort (1.010), so dropping it removes the most influential point.
#   * P10 has 4,625 scRNA cells but ZERO of them are malignant, and the sarcomatoid
#     programme is a malignant-cell programme -- so an scRNA-derived score could not
#     exist for it either.  It DOES have 976 malignant cells in the scATAC compartment
#     (more than P8, P1, P14, P3 or P13, all of which are plotted), so once the score
#     is computed from scATAC that obstacle is gone and P10 is included.  Its lack of
#     malignant scRNA is instead carried into the figure as an open point: P10 and P23
#     are the two tumours whose malignant compartment has no RNA corroboration.
#
# Fix.  Score the same cNMF programme directly on the scATAC malignant cells with
# ArchR::addModuleScore over the GeneScoreMatrix -- the recipe already used in
# git_repo/tumor_analysis/scatac_tumor_SOX9.R (cNMF20 of the 25-programme malignant
# spectra, nBin 25, nBgd 100).  This gives every scATAC tumour a score, P23
# included, from the modality that actually measured it.
#
# Module scores are z-scored ACROSS the 25 programmes WITHIN each cell before the
# per-sample summary.  Without that, a cell with globally high gene scores (deeper
# library, or a copy-number gain over the programme's genes) scores high on every
# programme, and the per-sample means would track coverage rather than phenotype.
#
# Validation: the new scATAC score is correlated against the old scRNA score over
# the 9 tumours both have.  A method that disagrees with the published score on
# those 9 has no business assigning one to P23.
#
# NOTHING is written back into tumor_compartment/scatac_ArchR: addModuleScore works
# on the in-memory project, the output directory is redirected here, and
# saveArchRProject() is never called.
#
# Input : tumor_compartment/scatac_ArchR,
#         git_repo/files/malignant_cnmf_genelist_25_nfeat_5000.rds,
#         tumor_compartment/scrna/cnmf20_sarcomatoid_sample_order.csv (validation),
#         BAP1_genescore_per_sample.csv
# Output: scATAC_sarcscore_per_sample.csv, scATAC_sarcscore_vs_scRNA.csv,
#         Plots/scATAC_BAP1_vs_sarcscore_withP23.pdf
###############################################################################
suppressPackageStartupMessages({
  library(ArchR); library(ggplot2); library(ggpubr); library(ggrepel)
})
addArchRGenome("hg38"); addArchRThreads(4)

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
TCOMP  <- file.path(ROOT, "tumor_compartment")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

SARC_MODULE <- "cNMF20"    # the sarcomatoid programme (VIM, AXL, CAV1, FBN1, THY1, ...)
TOP_GENES   <- 30          # top genes per programme, as in scatac_tumor_SOX9.R
NBIN        <- 25
NBGD        <- 100
MIN_GENES   <- 5           # a programme needs this many genes in the gene score matrix
DROP        <- c("normal1", "P11_HOX")   # P10 is kept -- see header
## tumours whose malignant compartment has no matched malignant scRNA
NO_RNA_MALIG <- c("P10", "P23")

## ---- project ------------------------------------------------------------------
proj <- loadArchRProject(file.path(TCOMP, "scatac_ArchR"), showLogo = FALSE)
proj@projectMetadata$outputDirectory <- file.path(OUTDIR, "archr_scratch")
dir.create(file.path(OUTDIR, "archr_scratch"), showWarnings = FALSE)
cd <- getCellColData(proj)
grp <- if ("Sample3" %in% colnames(cd)) "Sample3" else "Sample"
cat("malignant project:", nCells(proj), "cells\n"); print(table(cd[[grp]]))

keep <- rownames(cd)[!cd[[grp]] %in% DROP & !grepl("^RPL", cd[[grp]])]
proj <- proj[keep, ]
cat("\nafter dropping", paste(DROP, collapse = ", "), "and the RPL normals:",
    nCells(proj), "cells\n")
print(table(getCellColData(proj)[[grp]]))

## ---- the cNMF programmes, restricted to genes the gene score matrix has -------
spectra <- readRDS(file.path(ROOT, "git_repo", "files",
                             "malignant_cnmf_genelist_25_nfeat_5000.rds"))
gnames  <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
spectra <- lapply(spectra, function(x) { x <- head(x, TOP_GENES); x[x %in% gnames] })
## cNMF16 keeps 1 gene of 30 -- addModuleScore drops to a vector and fails
## ("'x' must be an array of at least two dimensions").  A one-gene programme is not
## a programme, so it is removed rather than worked around.  Every other programme
## keeps 28-30 of 30, cNMF20 included.
drop_prog <- names(which(lengths(spectra) < MIN_GENES))
spectra   <- spectra[lengths(spectra) >= MIN_GENES]
stopifnot(SARC_MODULE %in% names(spectra))
cat("\ngenes kept per programme (of", TOP_GENES, "):",
    paste(range(lengths(spectra)), collapse = "-"),
    "| programmes dropped for <", MIN_GENES, "genes:",
    if (length(drop_prog)) paste(drop_prog, collapse = ", ") else "none",
    "| programmes used:", length(spectra), "\n")
cat(SARC_MODULE, "genes:", paste(spectra[[SARC_MODULE]], collapse = ", "), "\n")

## ---- module scores -------------------------------------------------------------
## All 25 programmes, because the within-cell z-score below needs the full set.
proj@cellColData <- proj@cellColData[, !colnames(proj@cellColData) %in%
                                     c(names(spectra), drop_prog)]
proj <- addModuleScore(ArchRProj = proj, useMatrix = "GeneScoreMatrix", name = "",
                       features = spectra, nBin = NBIN, nBgd = NBGD, seed = 1,
                       threads = getArchRThreads(),
                       logFile = createLogFile("addModuleScore"))
colnames(proj@cellColData) <- gsub("^\\.", "", colnames(proj@cellColData))
cd <- as.data.frame(getCellColData(proj))
stopifnot(SARC_MODULE %in% colnames(cd))

M  <- as.matrix(cd[, names(spectra)])
Mz <- t(scale(t(M)))                       # z across the 25 programmes, within a cell
cd$sarc_cell <- Mz[, SARC_MODULE]

## ---- per-sample score ------------------------------------------------------------
## Median over the sample's malignant cells: the distributions are skewed and a few
## very high cells would otherwise set a whole tumour's score.
S <- do.call(rbind, lapply(split(cd, cd[[grp]]), function(d) data.frame(
  sample = d[[grp]][1], n_cells = nrow(d),
  sarc_score_atac = median(d$sarc_cell),
  sarc_mean = mean(d$sarc_cell), sarc_sd = sd(d$sarc_cell),
  stringsAsFactors = FALSE)))
S <- S[order(-S$sarc_score_atac), ]
write.csv(S, "scATAC_sarcscore_per_sample.csv", row.names = FALSE)
cat("\n=== sarcomatoid score from scATAC (median z of", SARC_MODULE, ") ===\n")
print(transform(S, sarc_score_atac = round(sarc_score_atac, 3),
                sarc_mean = round(sarc_mean, 3), sarc_sd = round(sarc_sd, 3)),
      row.names = FALSE)

## ---- validation against the published scRNA score ---------------------------------
sc <- read.csv(file.path(TCOMP, "scrna", "cnmf20_sarcomatoid_sample_order.csv"),
               stringsAsFactors = FALSE)
V  <- merge(S[, c("sample", "n_cells", "sarc_score_atac")],
            data.frame(sample = sc$sampleID, sarc_score_rna = sc$x), by = "sample")
ct <- suppressWarnings(cor.test(V$sarc_score_atac, V$sarc_score_rna, method = "spearman"))
V  <- V[order(-V$sarc_score_rna), ]
write.csv(V, "scATAC_sarcscore_vs_scRNA.csv", row.names = FALSE)
cat(sprintf("\n=== validation: scATAC vs published scRNA score, n = %d tumours ===\n", nrow(V)))
print(transform(V, sarc_score_atac = round(sarc_score_atac, 3),
                sarc_score_rna = round(sarc_score_rna, 3)), row.names = FALSE)
cat(sprintf("Spearman rho = %.3f, p = %.3g\n", ct$estimate, ct$p.value))
cat("samples scored here but absent from the scRNA file:",
    paste(setdiff(S$sample, V$sample), collapse = ", "), "\n")

## ---- the figure, now with P23 -------------------------------------------------------
gs <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
D  <- merge(S[, c("sample", "n_cells", "sarc_score_atac")],
            data.frame(sample = gs$Sample, BAP1_status = gs$BAP1_status,
                       BAP1_gs = gs$BAP1_genescore), by = "sample")
D$BAP1_status <- factor(D$BAP1_status, levels = c("retained", "lost"))
D$rna_malig   <- !D$sample %in% NO_RNA_MALIG

## Reported both ways.  P10 and P23 are the tumours the earlier figure could not
## show, and they sit at opposite ends of the BAP1 axis, so the effect of adding
## them is stated rather than buried.
for (nm in list(c("all tumours", ""), c("excluding P10", "P10"),
                c("excluding P10 and P23", "P10|P23"))) {
  d <- if (nzchar(nm[2])) D[!grepl(paste0("^(", nm[2], ")$"), D$sample), ] else D
  ct2 <- suppressWarnings(cor.test(d$BAP1_gs, d$sarc_score_atac, method = "spearman"))
  cat(sprintf("BAP1 gene score vs scATAC sarc score, %-22s n = %2d, rho = %6.3f, p = %.3g\n",
              nm[1], nrow(d), ct2$estimate, ct2$p.value))
}
cs <- suppressWarnings(cor.test(D$BAP1_gs, D$sarc_score_atac, method = "spearman"))
cat("\n=== per-tumour values ===\n")
print(D[order(-D$sarc_score_atac), ], row.names = FALSE, digits = 3)

sp <- ggplot(D, aes(BAP1_gs, sarc_score_atac)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey55", linewidth = .4,
              linetype = "dashed", formula = y ~ x) +
  geom_point(aes(color = BAP1_status, shape = rna_malig), size = 2.4, stroke = .9) +
  geom_text_repel(aes(label = sample), size = 2.4, max.overlaps = 100, seed = 1) +
  stat_cor(method = "spearman", size = 2.6) +
  scale_color_manual(values = c(retained = "#2471a3", lost = "#c0392b"),
                     labels = c(retained = "BAP1-retained", lost = "BAP1-lost"),
                     name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "malignant scRNA available",
                                `FALSE` = "no malignant scRNA"), name = NULL) +
  gtheme_no_rot + xlab("BAP1 gene score (scATAC)") +
  ylab(paste0("sarcomatoid score (", SARC_MODULE, ", scATAC)")) +
  ggtitle(sprintf("scATAC: BAP1 vs histology axis (n = %d)", nrow(D)))
pdf(file.path("Plots", "scATAC_BAP1_vs_sarcscore_withP23.pdf"), width = 4.8, height = 3.8)
print(sp); dev.off()

cat("\nDONE\n")
