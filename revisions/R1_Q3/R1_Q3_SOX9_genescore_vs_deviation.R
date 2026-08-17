#!/usr/bin/env Rscript
# R1_Q3_SOX9_genescore_vs_deviation.R
# Does SOX9 gene-score correlate with SOX9 chromVAR deviation, and is that
# correlation independent of the active-enhancer program (AP-1/RUNX) the cells
# were gated on? Tests "real regulator" vs "colinear passenger".
#
#   - per-cell SOX9 gene-score vs SOX9 deviation z (Pearson + Spearman)
#   - same self-correlation for bona-fide drivers (FOS/JUN/RUNX) as a yardstick
#   - PARTIAL correlation of SOX9 gs~dev controlling for AP-1 deviation and for
#     the SOX9 module score (the gating axis)

suppressPackageStartupMessages({library(ArchR); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
proj <- loadArchRProject(PROJ, showLogo = FALSE)

## ---- matrices ----
gsm <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gs_genes <- rowData(gsm)$name
GS <- assays(gsm)[[1]]; rownames(GS) <- gs_genes
mm  <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
Z   <- assays(mm)[["z"]]
# align cells
cells <- intersect(colnames(GS), colnames(Z))
GS <- GS[, cells]; Z <- Z[, cells]
gate <- setNames(proj$SOX9_module_gate, proj$cellNames)[cells]

motif_of <- function(sym) { h <- grep(paste0("^", sym, "_[0-9]"), rownames(Z)); if(length(h)) rownames(Z)[h[1]] else NA }
getgs <- function(sym) if (sym %in% rownames(GS)) as.numeric(GS[sym, ]) else rep(NA, length(cells))
getz  <- function(sym) { m <- motif_of(sym); if (is.na(m)) rep(NA, length(cells)) else as.numeric(Z[m, ]) }

## ---- self gene-score ~ deviation correlations ----
tfs <- c("SOX9","FOS","FOSL1","FOSL2","JUN","JUNB","JUND","RUNX1","RUNX2","CTCF","NFIC","TEAD1")
rows <- list()
for (tf in tfs) {
  g <- getgs(tf); z <- getz(tf)
  if (all(is.na(g)) || all(is.na(z))) next
  rp <- suppressWarnings(cor(g, z, method="pearson",  use="complete.obs"))
  rs <- suppressWarnings(cor(g, z, method="spearman", use="complete.obs"))
  rows[[tf]] <- data.frame(TF=tf, motif=motif_of(tf), pearson=round(rp,3), spearman=round(rs,3))
}
selfcor <- do.call(rbind, rows)
cat("=== self gene-score ~ own-motif deviation (across",length(cells),"cells) ===\n")
print(selfcor)
write.csv(selfcor, file.path(OUT, "SOX9_genescore_vs_deviation_selfcor.csv"), row.names=FALSE)

## ---- partial correlation: SOX9 gs~dev controlling for AP-1 program / module score ----
pcor <- function(x, y, z) {  # partial cor of x,y given z
  ok <- complete.cases(x,y,z); x<-x[ok]; y<-y[ok]; z<-z[ok]
  rxy<-cor(x,y); rxz<-cor(x,z); ryz<-cor(y,z)
  (rxy - rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
}
sox_g <- getgs("SOX9"); sox_z <- getz("SOX9")
ap1_z <- rowMeans(cbind(getz("FOSL2"), getz("JUN"), getz("JUNB")), na.rm=TRUE)  # AP-1 program activity
runx_z<- rowMeans(cbind(getz("RUNX1"), getz("RUNX2")), na.rm=TRUE)
# module score (gating axis) from the per-cell program file if available
modf <- file.path(OUT, "SOX9_program_module_P23_percell.csv")
mod <- NULL
if (file.exists(modf)) {
  md <- read.csv(modf, stringsAsFactors=FALSE)
  bc <- md[[1]]; sc_col <- grep("module|score", names(md), ignore.case=TRUE, value=TRUE)[1]
  mod <- setNames(md[[sc_col]], bc)[cells]
}
cat("\n=== SOX9 gene-score ~ SOX9 deviation: raw vs partial ===\n")
cat(sprintf("  raw Pearson                         : %.3f\n", cor(sox_g, sox_z, use="complete.obs")))
cat(sprintf("  partial | AP-1 deviation            : %.3f\n", pcor(sox_g, sox_z, ap1_z)))
cat(sprintf("  partial | RUNX deviation            : %.3f\n", pcor(sox_g, sox_z, runx_z)))
if (!is.null(mod) && sum(!is.na(mod))>10)
  cat(sprintf("  partial | SOX9 module score (gate)  : %.3f\n", pcor(sox_g, sox_z, as.numeric(mod))))

## ---- scatter plot SOX9 gs vs deviation ----
df <- data.frame(genescore=sox_g, deviation=sox_z, gate=gate)
df <- df[complete.cases(df) & df$gate %in% c("SOX9neg","SOX9pos"), ]
p <- ggplot(df, aes(genescore, deviation, color=gate)) +
  geom_point(size=.5, alpha=.5) + geom_smooth(method="lm", se=TRUE, color="black", linewidth=.6) +
  scale_color_manual(values=c(SOX9neg="#377EB8", SOX9pos="#E41A1C")) +
  labs(title="SOX9 gene-score vs chromVAR deviation",
       subtitle=sprintf("Pearson r = %.3f (n=%d malignant cells)", cor(sox_g,sox_z,use="complete.obs"), nrow(df)),
       x="SOX9 gene-score", y="SOX9 deviation z") +
  theme_bw(base_size=12)
ggsave(file.path(OUT, "SOX9_genescore_vs_deviation_scatter.pdf"), p, width=5, height=4.5)
cat("\nWrote selfcor CSV + scatter PDF\n")
