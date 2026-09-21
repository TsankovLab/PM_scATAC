###############################################################################
# R2_Q14 -- histology-adjusted BAP1 comparison in all three bulk cohorts, using a
#           MOLECULAR definition of BAP1 status.
#
#     expression ~ BAP1_group + sarc_score        (limma, per cohort)
#
# Why a molecular definition here.  Only MESOMICS has clinical BAP1 IHC; Bueno and
# TCGA have no BAP1 status at all.  Their genomic proxies are unusable:
#   Bueno FISH.chrom3  161 of 211 blank, and "del 3p" tumours have the same BAP1
#                      mRNA as the rest (3.72 vs 3.73) -- no signal.
#   TCGA  STATUS_3P    5 Lost, 4 Gained, the rest not called; BAP1 mRNA is HIGHER
#                      in "Lost" than in the unlabelled group -- wrong direction.
# So the only molecular definition with any discriminating power is BAP1 mRNA
# itself, split at tertiles (the split R2_Q14_BAP1expr_proxy_epithelioid.R already
# uses): BAP1-low = bottom third, BAP1-high = top third, middle third dropped.
#
# THE CIRCULARITY THIS INTRODUCES IS NOT MINOR, and it is the reason the clinical
# analysis is kept as the reference:
#   * Groups are defined by expression and then expression is tested, so any gene
#     correlated with BAP1 mRNA for reasons other than BAP1 status -- shared
#     regulation, tumour purity, RNA quality, chr3p copy number dragging neighbours
#     -- separates the groups by construction.
#   * BAP1 itself is guaranteed significant; it is excluded from the results.
#   * Discarding the middle tertile inflates the apparent effect relative to a
#     status label that would have put those tumours in one group or the other.
# The MESOMICS IHC calibration below quantifies how far the proxy is from the real
# thing, and it is the number to read before believing anything downstream.
#
# Input : bulkRNA_meso/bulk_RNA_studies.rds, bulk_RNA_studies_metadata.rds
#         BAP1_bulk_limma_sarcadjusted.csv (MESOMICS clinical, for comparison)
# Output: BAP1_bulk_molecular_limma.csv, BAP1_bulk_molecular_crosscohort.csv,
#         BAP1_proxy_vs_IHC_calibration.csv,
#         Plots/BAP1_bulk_molecular_volcano.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(ggrepel) })

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
BULK   <- file.path(ROOT, "bulkRNA_meso")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

RUNX    <- c("RUNX1", "RUNX2", "RUNX3")
IRF     <- paste0("IRF", 1:9)
STUDIES <- c("bueno", "tcga", "mesomics")
NTOP    <- 20

bl   <- readRDS(file.path(BULK, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))

get <- function(s){
  E  <- as.matrix(bl[[s]])
  md <- meta[[s]]
  md <- if (all(colnames(E) %in% rownames(md))) md[colnames(E), , drop = FALSE]
        else md[match(colnames(E), md$Sample), , drop = FALSE]
  list(E = E, md = md)
}

## ---- how well does the proxy recover the clinical label? -------------------------
d <- get("mesomics"); b <- as.numeric(d$E["BAP1", ]); ih <- as.character(d$md$IHC.BAP1)
q <- quantile(b, c(1/3, 2/3), na.rm = TRUE)
prox <- ifelse(b <= q[1], "BAP1low", ifelse(b >= q[2], "BAP1high", NA))
tb <- table(proxy = prox, IHC = ih)[, c("NO", "YES"), drop = FALSE]
cat("=== calibration: BAP1 mRNA tertile proxy vs clinical IHC (MESOMICS) ===\n")
cat("   IHC NO = BAP1 lost, YES = retained\n"); print(tb)
acc <- (tb["BAP1low", "NO"] + tb["BAP1high", "YES"]) / sum(tb)
cat(sprintf("concordance where the proxy commits: %.1f%% (%d of %d tumours)\n",
            100 * acc, tb["BAP1low","NO"] + tb["BAP1high","YES"], sum(tb)))
cat(sprintf("of %d IHC-lost tumours the proxy calls %d low, %d high, %d middle (dropped)\n",
            sum(ih == "NO", na.rm = TRUE), tb["BAP1low","NO"], tb["BAP1high","NO"],
            sum(ih == "NO" & is.na(prox), na.rm = TRUE)))
write.csv(as.data.frame.matrix(tb), "BAP1_proxy_vs_IHC_calibration.csv")

## ---- per cohort ---------------------------------------------------------------------
run <- function(s){
  d <- get(s); E <- d$E; md <- d$md
  bap1 <- as.numeric(E["BAP1", ]); sarc <- as.numeric(md$sarc_score)
  q <- quantile(bap1, c(1/3, 2/3), na.rm = TRUE)
  grp <- ifelse(bap1 <= q[1], "lost", ifelse(bap1 >= q[2], "retained", NA))
  ok  <- !is.na(grp) & is.finite(sarc)
  E <- E[, ok, drop = FALSE]
  G <- factor(grp[ok], levels = c("retained", "lost")); S <- sarc[ok]
  keep <- apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0 && mean(x > 0) > 0.5)
  E <- E[keep & rownames(E) != "BAP1", , drop = FALSE]     # BAP1 is circular by design
  cat(sprintf("\n--- %s: %d tumours (%d lost / %d retained), %d genes ---\n",
              s, ncol(E), sum(G == "lost"), sum(G == "retained"), nrow(E)))
  cat(sprintf("    sarc score %.2f lost vs %.2f retained\n",
              mean(S[G == "lost"]), mean(S[G == "retained"])))
  fa <- eBayes(lmFit(E, model.matrix(~ G + S)), trend = TRUE, robust = TRUE)
  fu <- eBayes(lmFit(E, model.matrix(~ G)),     trend = TRUE, robust = TRUE)
  a <- topTable(fa, coef = "Glost", number = Inf, sort.by = "none")
  u <- topTable(fu, coef = "Glost", number = Inf, sort.by = "none")
  h <- topTable(fa, coef = "S",     number = Inf, sort.by = "none")
  i <- match(rownames(a), rownames(u)); j <- match(rownames(a), rownames(h))
  R <- data.frame(cohort = s, gene = rownames(a),
                  logFC_adj = a$logFC, t_adj = a$t, P_adj = a$P.Value, FDR_adj = a$adj.P.Val,
                  t_unadj = u$t[i], FDR_unadj = u$adj.P.Val[i],
                  t_sarc = h$t[j], FDR_sarc = h$adj.P.Val[j], stringsAsFactors = FALSE)
  for (aa in c(0.01, 0.05, 0.10)) cat(sprintf(
    "    FDR < %.2f : adjusted %5d | unadjusted %5d | histology %5d\n", aa,
    sum(R$FDR_adj < aa), sum(R$FDR_unadj < aa), sum(R$FDR_sarc < aa)))
  attr(R, "camera") <- cameraPR(fa$t[, "Glost"],
                                list(IRF = which(rownames(E) %in% IRF),
                                     RUNX = which(rownames(E) %in% RUNX)))
  R
}
res <- lapply(setNames(STUDIES, STUDIES), run)
ALL <- do.call(rbind, res); rownames(ALL) <- NULL
write.csv(ALL, "BAP1_bulk_molecular_limma.csv", row.names = FALSE)

for (s in STUDIES){
  cat(sprintf("\n=== %s: cameraPR on the adjusted BAP1 t ===\n", s))
  print(attr(res[[s]], "camera"), digits = 3)
}
cat("\n=== RUNX, all cohorts (molecular proxy, histology-adjusted) ===\n")
print(ALL[ALL$gene %in% RUNX, c("cohort","gene","logFC_adj","t_adj","P_adj","FDR_adj","t_sarc")],
      row.names = FALSE, digits = 3)
cat("\n=== IRF, all cohorts ===\n")
print(ALL[ALL$gene %in% IRF, c("cohort","gene","logFC_adj","t_adj","P_adj","FDR_adj")],
      row.names = FALSE, digits = 3)

## ---- reproducible across cohorts? -----------------------------------------------------
w <- reshape(ALL[, c("cohort","gene","t_adj")], idvar = "gene",
             timevar = "cohort", direction = "wide")
names(w) <- sub("^t_adj\\.", "", names(w))
w <- w[complete.cases(w), ]
cat(sprintf("\n=== %d genes measured in all three cohorts ===\n", nrow(w)))
for (p in combn(STUDIES, 2, simplify = FALSE))
  cat(sprintf("  %-9s vs %-9s Spearman rho = %6.3f\n", p[1], p[2],
              cor(w[[p[1]]], w[[p[2]]], method = "spearman")))
w$n_consistent <- rowSums(sign(w[, STUDIES]) == sign(w[[STUDIES[1]]]))
FDR <- reshape(ALL[, c("cohort","gene","FDR_adj")], idvar = "gene",
               timevar = "cohort", direction = "wide")
names(FDR) <- sub("^FDR_adj\\.", "", names(FDR))
m <- merge(w, FDR, by = "gene", suffixes = c("_t", "_fdr"))
hit <- m[rowSums(m[, paste0(STUDIES, "_fdr")] < 0.05, na.rm = TRUE) == 3 &
         m$n_consistent == 3, ]
cat(sprintf("genes at FDR < 0.05 in ALL THREE cohorts with the same sign: %d\n", nrow(hit)))
if (nrow(hit)) print(head(hit[order(-abs(hit$mesomics_t)),
                             c("gene", paste0(STUDIES, "_t"))], 25),
                     row.names = FALSE, digits = 3)
write.csv(m, "BAP1_bulk_molecular_crosscohort.csv", row.names = FALSE)

## ---- volcano, one panel per cohort ------------------------------------------------------
ALL$key <- paste(ALL$cohort, ALL$gene)
ALL$rank <- ave(ALL$P_adj, ALL$cohort, FUN = function(x) rank(x, ties.method = "min"))
ALL$lab  <- ifelse(ALL$rank <= NTOP, ALL$gene, NA)
pv <- ggplot(ALL, aes(logFC_adj, -log10(P_adj))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey78", size = .7, alpha = .55) +
  geom_point(data = subset(ALL, FDR_adj < 0.05), colour = "grey45", size = .8) +
  geom_point(data = subset(ALL, rank <= NTOP), colour = "grey15", size = 1.3) +
  geom_text_repel(aes(label = lab), size = 1.9, colour = "grey15", max.overlaps = 30,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  facet_wrap(~ cohort, scales = "free", nrow = 1) +
  gtheme_no_rot +
  xlab("BAP1 effect, low - high BAP1 mRNA (log2, histology-adjusted)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle("Bulk RNA, MOLECULAR BAP1 proxy (mRNA tertiles) adjusted for sarcomatoid score")
pdf(file.path("Plots", "BAP1_bulk_molecular_volcano.pdf"), width = 11, height = 4)
print(pv); dev.off()

cat("\nDONE\n")
