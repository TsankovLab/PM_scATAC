#!/usr/bin/env Rscript
# R1_Q3_SOX9_cluster_purity.R
#
# Is the SOX9-high pseudobulk (P23 cluster C9) actually pure for SOX9, or is it
# contaminated with SOX9-negative cells (diluting the ChromBPNet SOX9 signal)?
#
# Per-cell SOX9 gene activity (GeneScoreMatrix) and SOX9 motif activity
# (chromVAR MotifMatrix z) in C9 (high) vs C4 (low), with a contamination
# estimate: fraction of C9 cells whose SOX9 level looks "low".

suppressMessages(library(ArchR))
addArchRThreads(1)
OUT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
proj <- loadArchRProject(
  "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23",
  showLogo = FALSE)

grp <- proj$SOX9_cluster                      # SOX9_high_P23 / SOX9_low_P23
message("group sizes: "); print(table(grp))

## ---- SOX9 gene activity ---------------------------------------------------
gsm <- getMatrixFromProject(proj, "GeneScoreMatrix")
gs_names <- rowData(gsm)$name
gi <- which(gs_names == "SOX9")
gsv <- as.numeric(assay(gsm)[gi, ])           # per cell, aligned to colnames(gsm)
gs <- data.frame(cell = colnames(gsm), gene_SOX9 = gsv)
gs$grp <- grp[match(gs$cell, proj$cellNames)]

## ---- SOX9 motif (chromVAR) if available -----------------------------------
mz <- rep(NA_real_, nrow(gs))
ok <- tryCatch({
  mm <- getMatrixFromProject(proj, "MotifMatrix")
  zi <- grep("SOX9", rownames(mm), ignore.case = TRUE, value = TRUE)
  zi <- zi[grepl("z:|deviations", zi) | TRUE][1]
  za <- assays(mm)[["z"]]
  ri <- grep("SOX9", rownames(za), value = TRUE)[1]
  message("chromVAR SOX9 row: ", ri)
  mz <- as.numeric(za[ri, match(gs$cell, colnames(za))]); TRUE
}, error = function(e) { message("MotifMatrix unavailable: ", conditionMessage(e)); FALSE })
gs$motif_SOX9_z <- mz

write.csv(gs, file.path(OUT, "SOX9_C9_purity_percell.csv"), row.names = FALSE)

## ---- summaries ------------------------------------------------------------
hi <- gs$gene_SOX9[gs$grp == "SOX9_high_P23"]
lo <- gs$gene_SOX9[gs$grp == "SOX9_low_P23"]
cat("\n#### SOX9 gene activity per group ####\n")
print(rbind(
  SOX9_high_C9 = summary(hi),
  SOX9_low_C4  = summary(lo)))

# contamination: C9 cells at/below the C4 (low) median or below a global gate
lo_med <- median(lo)
gate   <- quantile(gs$gene_SOX9, 0.5)         # global median as a simple SOX9+/- gate
cat(sprintf("\nC4(low) median SOX9 gene score = %.3f\nglobal median (gate) = %.3f\n", lo_med, gate))
cat(sprintf("C9(high) cells at/below C4 median        : %d / %d (%.1f%%)\n",
            sum(hi <= lo_med), length(hi), 100*mean(hi <= lo_med)))
cat(sprintf("C9(high) cells below global-median gate   : %d / %d (%.1f%%)  <- putative SOX9- contaminants\n",
            sum(hi < gate), length(hi), 100*mean(hi < gate)))
cat(sprintf("C9(high) cells with ZERO SOX9 gene score  : %d / %d (%.1f%%)\n",
            sum(hi == 0), length(hi), 100*mean(hi == 0)))

if (ok) {
  hz <- gs$motif_SOX9_z[gs$grp=="SOX9_high_P23"]
  lz <- gs$motif_SOX9_z[gs$grp=="SOX9_low_P23"]
  cat("\n#### SOX9 chromVAR motif z per group ####\n")
  print(rbind(SOX9_high_C9 = summary(hz), SOX9_low_C4 = summary(lz)))
  cat(sprintf("C9 cells with SOX9 motif z <= 0 : %.1f%%\n", 100*mean(hz <= 0, na.rm=TRUE)))
}
message("\nWrote ", file.path(OUT, "SOX9_C9_purity_percell.csv"))
