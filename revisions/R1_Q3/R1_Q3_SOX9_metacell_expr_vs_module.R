#!/usr/bin/env Rscript
# R1_Q3_SOX9_metacell_expr_vs_module.R
# Metacell-level coherence of SOX9 expression with its program. Two module
# definitions, both correlated with SOX9 expression across tumor metacells:
#   (1) "ours"   : the top-20 P4 cluster-32 SOX9 module we used for the ATAC gate
#                  (computed incl. AND excl. SOX9 itself, to avoid self-correlation)
#   (2) "SCENIC" : the actual SOX9(+) regulon AUC, aggregated from per-cell SCENIC
#                  AUC to metacells via the cells_merged membership
# RUNX2 included as a comparator (driver). NOTE: this is a TRANSCRIPTIONAL
# coherence test (does SOX9 mRNA track its co-expressed program), NOT a DNA-binding
# test -- the footprint/ChromBPNet analyses address occupancy.

suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT  <- file.path(ROOT, "git_repo_claude/R1_Q3")
SCR  <- file.path(ROOT, "tumor_compartment/scrna")
MC   <- file.path(SCR, "metacells.rds")
MODF <- file.path(OUT, "cluster32_top20_genes.csv")
AUCF <- file.path(SCR, "SCENIC/vg_5000_mw_tss500bp/tumor_programs/auc_mtx.csv")
TUMOR <- c("P1","P11","P11_HOX","P12","P13","P4","P5","P8")

## ---- metacells ----
message("Loading metacells ...")
mc <- readRDS(MC)
mc <- mc[, mc$sampleID3 %in% TUMOR]
message("Tumor metacells: ", ncol(mc), " | samples: ", paste(sort(unique(mc$sampleID3)),collapse=", "))
mc <- NormalizeData(mc, verbose = FALSE)
expr <- GetAssayData(mc, assay = "RNA", layer = "data")

## ---- module (ours) ----
modg <- read.csv(MODF, stringsAsFactors=FALSE)$gene
modg <- modg[modg %in% rownames(mc)]
modg_no <- setdiff(modg, c("SOX9"))
mc <- AddModuleScore(mc, features=list(modg),    name="ourMod",      seed=1)
mc <- AddModuleScore(mc, features=list(modg_no), name="ourModNoSOX9", seed=1)

## ---- SCENIC AUC aggregated to metacells ----
message("Loading SCENIC AUC ...")
auc <- read.csv(AUCF, row.names=1, check.names=FALSE)
colnames(auc) <- gsub("\\(\\+\\)|\\(\\-\\)", "", colnames(auc))
agg_auc <- function(reg) {
  v <- setNames(auc[[reg]], rownames(auc))
  vapply(mc$cells_merged, function(cm){ bc<-strsplit(cm,",")[[1]]; mean(v[bc], na.rm=TRUE) }, numeric(1))
}
has <- function(reg) reg %in% colnames(auc)
sox9_auc <- if (has("SOX9")) agg_auc("SOX9") else rep(NA, ncol(mc))
runx2_auc<- if (has("RUNX2")) agg_auc("RUNX2") else rep(NA, ncol(mc))
cat("metacells with non-NA SOX9 AUC: ", sum(!is.na(sox9_auc)), "/", ncol(mc), "\n")

## ---- assemble + correlate ----
df <- data.frame(
  sample      = mc$sampleID3,
  SOX9_expr   = as.numeric(expr["SOX9", ]),
  RUNX2_expr  = if ("RUNX2" %in% rownames(expr)) as.numeric(expr["RUNX2", ]) else NA,
  ourMod      = mc$ourMod1,
  ourModNoSOX9= mc$ourModNoSOX91,
  SOX9_auc    = sox9_auc,
  RUNX2_auc   = runx2_auc)
write.csv(df, file.path(OUT, "SOX9_metacell_expr_module_table.csv"), row.names=FALSE)

corr <- function(x,y,lab){
  ok<-complete.cases(x,y)
  data.frame(comparison=lab, n=sum(ok),
             pearson=round(cor(x[ok],y[ok]),3),
             spearman=round(cor(x[ok],y[ok],method="spearman"),3))
}
summ <- rbind(
  corr(df$SOX9_expr, df$ourMod,       "SOX9 expr ~ ourModule(incl SOX9)"),
  corr(df$SOX9_expr, df$ourModNoSOX9, "SOX9 expr ~ ourModule(excl SOX9)"),
  corr(df$SOX9_expr, df$SOX9_auc,     "SOX9 expr ~ SCENIC SOX9 regulon AUC"),
  corr(df$RUNX2_expr,df$ourModNoSOX9, "RUNX2 expr ~ ourModule(excl SOX9)  [comparator]"),
  corr(df$RUNX2_expr,df$RUNX2_auc,    "RUNX2 expr ~ SCENIC RUNX2 regulon AUC [comparator]"))
cat("\n=== metacell correlations (n tumor metacells) ===\n"); print(summ)
write.csv(summ, file.path(OUT, "SOX9_metacell_expr_module_correlations.csv"), row.names=FALSE)

## ---- plots ----
mk <- function(yvar, ylab, rr){
  d <- df[complete.cases(df[,c("SOX9_expr",yvar)]),]
  ggplot(d, aes(.data[["SOX9_expr"]], .data[[yvar]], color=sample)) +
    geom_point(size=1, alpha=.7) + geom_smooth(method="lm", se=TRUE, color="black", linewidth=.5) +
    labs(title=paste0("SOX9 expr vs ", ylab, " (metacells)"),
         subtitle=sprintf("Pearson r=%.2f  Spearman=%.2f  (n=%d metacells)",
                          rr$pearson, rr$spearman, rr$n),
         x="SOX9 expression (metacell, log-norm)", y=ylab) +
    theme_bw(base_size=11)
}
p1 <- mk("ourModNoSOX9","our SOX9 module (excl SOX9)", summ[2,])
p2 <- mk("SOX9_auc","SCENIC SOX9 regulon AUC", summ[3,])
pdf(file.path(OUT,"SOX9_metacell_expr_vs_module.pdf"), width=6, height=5); print(p1); print(p2); dev.off()
cat("\nWrote table + correlations + plot\n")
