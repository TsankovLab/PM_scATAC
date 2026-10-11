# Spec 50 (comparison direction: NTC - KO by default => NES>0 = higher in NTC / lower after KO; DIFF_DIR=KO_NTC for the reverse): is there an antigen-presentation (AP) signal specific to the TEAD4 KO? Preranked GSEA (fgsea) on the FULL KO-vs-NTC statistic
# (bulk DESeq2 Wald stat per line; CROP-seq screen Welch t) for AP gene sets (KEGG, Reactome class I / class II, GO BP/CC, curated MHC-I / MHC-II /
# APM core list) plus IFN-gamma/alpha Hallmark as context. "Specific" = TEAD4's NES compared with the NES of every OTHER KO in the same source
# (null = other KOs vs NTC), per source. Also a gene-level view of the curated AP genes in TEAD4 KO vs NTC.
# needs: spec1_prepare_stats.R (spec_bulk_stats.rds, spec_screen_stats.rds)
# outputs: spec50_AP_fgsea_all_KOs.csv, spec50_AP_TEAD4_specificity.csv, spec50_AP_core_genes_TEAD4.csv,
#          fig4D_TEAD4_AP_NES_heatmap, fig4E_TEAD4_AP_vs_other_KOs, fig4F_TEAD4_AP_core_genes
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })
KO <- "TEAD4"
DIFF_DIR <- Sys.getenv("DIFF_DIR", "NTC_KO"); SGN <- if (DIFF_DIR == "NTC_KO") -1 else 1; DLAB <- if (DIFF_DIR == "NTC_KO") "NTC-KO" else "KO-NTC"; DSUF <- if (DIFF_DIR == "NTC_KO") "" else "_KOminusNTC"   # default NTC - KO (positive = higher in NTC, i.e. lower after KO)
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))
SRC <- c("CROP-seq (H2052)", paste0(LINES, " bulk")); short <- c("Screen", "H2052", "MSTO", "H2452", "H28"); names(short) <- SRC
stats <- c(list(SGN * ss$t0), lapply(LINES, function(cl) SGN * bs[[cl]]$stat0)); names(stats) <- SRC
KOS <- intersect(ORDER, Reduce(intersect, lapply(stats, colnames)))   # KOs present in every source (TEAD4 + the other shared KOs)
cat("KOs compared:", paste(KOS, collapse = ","), "\n")

## ---------------------------------------------------------------- gene sets
m <- msigdbr(species = "Homo sapiens"); sets_all <- split(m$gene_symbol, m$gs_name); sets_all <- lapply(sets_all, unique)
WANT <- c(KEGG_AP = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", REACTOME_MHC_II = "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
          REACTOME_MHC_I = "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
          REACTOME_MHC_I_loading = "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
          GOBP_AP = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", GOBP_MHC_I = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
          GOBP_MHC_II = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
          GOCC_MHC_I_complex = "GOCC_MHC_CLASS_I_PROTEIN_COMPLEX", GOCC_MHC_II_complex = "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX",
          HALLMARK_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE", HALLMARK_IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE")
stopifnot(all(WANT %in% names(sets_all))); sets <- sets_all[WANT]; names(sets) <- names(WANT)
MHC1 <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M"); APM <- c("TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "NLRC5")
MHC2 <- c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CD74", "CIITA", "CTSS", "CTSL", "IFI30")
sets[["CURATED_MHC_I_APM"]] <- c(MHC1, APM); sets[["CURATED_MHC_II"]] <- MHC2
SETLAB <- c(KEGG_AP = "KEGG antigen processing & presentation", REACTOME_MHC_II = "Reactome MHC class II antigen presentation", REACTOME_MHC_I = "Reactome class I MHC antigen processing",
            REACTOME_MHC_I_loading = "Reactome class I peptide loading", GOBP_AP = "GO BP antigen processing & presentation", GOBP_MHC_I = "GO BP peptide antigen via MHC I",
            GOBP_MHC_II = "GO BP exogenous peptide via MHC II", GOCC_MHC_I_complex = "GO CC MHC class I complex", GOCC_MHC_II_complex = "GO CC MHC class II complex",
            HALLMARK_IFNG = "Hallmark IFN-gamma response (context)", HALLMARK_IFNA = "Hallmark IFN-alpha response (context)",
            CURATED_MHC_I_APM = "Curated MHC-I + APM genes", CURATED_MHC_II = "Curated MHC-II genes")

## ---------------------------------------------------------------- fgsea, every KO x source
set.seed(1); res <- list()
for (s in SRC) for (k in KOS) {
  x <- stats[[s]][, k]; x <- x[is.finite(x)]; x <- x[!duplicated(names(x))]
  r <- fgseaMultilevel(sets, x, minSize = 3, maxSize = 1500, eps = 0)
  res[[length(res) + 1]] <- data.frame(source = s, KO = k, set = r$pathway, size = r$size, NES = r$NES, pval = r$pval, padj = r$padj, stringsAsFactors = FALSE)
}
res <- do.call(rbind, res); write.csv(res, file.path(OUT, paste0("spec50_AP_fgsea_all_KOs", DSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- specificity: TEAD4 NES vs the other KOs' NES (same source, same set)
spec <- do.call(rbind, lapply(split(res, list(res$source, res$set), drop = TRUE), function(d) {
  t4 <- d[d$KO == KO, ]; ot <- d[d$KO != KO, ]; if (!nrow(t4) || !nrow(ot)) return(NULL)
  data.frame(source = t4$source, set = t4$set, TEAD4_NES = t4$NES, TEAD4_p = t4$pval, TEAD4_padj = t4$padj, other_median_NES = median(ot$NES, na.rm = TRUE),
             n_other = sum(!is.na(ot$NES)), rank_of_TEAD4_by_NES = sum(ot$NES < t4$NES, na.rm = TRUE) + 1,   # 1 = most negative among KOs
             rank_by_abs = sum(abs(ot$NES) > abs(t4$NES), na.rm = TRUE) + 1, empirical_z = (t4$NES - mean(ot$NES, na.rm = TRUE)) / sd(ot$NES, na.rm = TRUE))
})); rownames(spec) <- NULL
write.csv(spec, file.path(OUT, paste0("spec50_AP_TEAD4_specificity", DSUF, ".csv")), row.names = FALSE)
options(width = 220)
cat(sprintf("\n=== TEAD4 %s: AP gene sets, NES (positive = higher in %s) | rank among 8 KOs by NES (1 = most negative) | empirical z vs other KOs ===\n", DLAB, if (SGN < 0) "NTC" else "KO"))
for (st in names(SETLAB)) { d <- spec[spec$set == st, ]; d$source <- short[d$source]
  cat(sprintf("\n-- %s (n genes in set=%d)\n", SETLAB[st], length(sets[[st]]))); print(format(d[match(short, d$source), c("source", "TEAD4_NES", "TEAD4_padj", "other_median_NES", "rank_of_TEAD4_by_NES", "empirical_z")], digits = 2), row.names = FALSE) }

## ---------------------------------------------------------------- fig4D: TEAD4 NES heatmap (sets x sources) with padj stars
save_fig(paste0("fig4D_TEAD4_AP_NES_heatmap", DSUF), 8.5, 6.4, function() {
  M <- sapply(SRC, function(s) sapply(names(SETLAB), function(st) spec$TEAD4_NES[spec$source == s & spec$set == st][1]))
  P <- sapply(SRC, function(s) sapply(names(SETLAB), function(st) spec$TEAD4_padj[spec$source == s & spec$set == st][1]))
  vmax <- max(abs(M), na.rm = TRUE); par(mar = c(4.5, 17, 3.6, 1), family = "sans")
  image(seq_along(SRC), seq_len(nrow(M)), t(M)[, nrow(M):1], col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(nrow(M))) for (j in seq_along(SRC)) { p <- P[i, j]; v <- M[i, j]; if (is.na(v)) next
    text(j, nrow(M) - i + 1, sprintf("%.2f%s", v, if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""), cex = 0.7, col = INK) }
  axis(1, at = seq_along(SRC), labels = FALSE); text(seq_along(SRC), 0.3, short[SRC], srt = 45, adj = 1, xpd = NA, cex = 0.8)
  axis(2, at = nrow(M):1, labels = SETLAB[rownames(M)], las = 1, tick = FALSE, cex.axis = 0.7)
  abline(v = 1.5, col = "white", lwd = 3); abline(h = 2.5, col = "white", lwd = 3); box(col = MUTED, lwd = 0.6)
  mtext(sprintf("AP gene sets, TEAD4 %s (preranked GSEA)", DLAB), side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext(sprintf("colour/number = NES (red = higher in %s, blue = lower); * padj<0.05 ** <0.01 *** <0.001", if (SGN < 0) "NTC than TEAD4 KO" else "TEAD4 KO than NTC"), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
})

## ---------------------------------------------------------------- fig4E: TEAD4 NES vs other KOs' NES (null), per source, for selected AP sets
SEL <- c("KEGG_AP", "REACTOME_MHC_II", "GOBP_MHC_I", "CURATED_MHC_I_APM", "CURATED_MHC_II", "HALLMARK_IFNG")
save_fig(paste0("fig4E_TEAD4_AP_vs_other_KOs", DSUF), 12.5, 7, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (st in SEL) {
    d <- res[res$set == st, ]; d$s <- factor(short[d$source], levels = short); par(mar = c(4, 4, 2.2, 0.8))
    ylim <- range(c(d$NES, 0), na.rm = TRUE); plot(NA, xlim = c(0.5, 5.5), ylim = ylim + c(-0.1, 0.1) * diff(ylim), xaxt = "n", xlab = "", ylab = paste0("NES (", DLAB, ")"))
    axis(1, at = 1:5, labels = FALSE); text(1:5, par("usr")[3] - 0.04 * diff(par("usr")[3:4]), short, srt = 45, adj = 1, xpd = NA, cex = 0.75)
    abline(h = 0, lty = 2, col = INK); abline(v = 1.5, col = "grey70", lty = 3)
    for (i in 1:5) { o <- d[d$s == short[i] & d$KO != KO, ]; set.seed(2); points(jitter(rep(i, nrow(o)), amount = 0.15), o$NES, pch = 16, col = "grey55", cex = 0.9) }
    for (i in 1:5) { t4 <- d[d$s == short[i] & d$KO == KO, ]; points(i, t4$NES, pch = 23, bg = KO_COL[[KO]], col = INK, cex = 2, lwd = 1.2)
      if (is.finite(t4$padj) && t4$padj < 0.05) text(i, t4$NES, "*", pos = if (t4$NES > 0) 3 else 1, cex = 1.4, font = 2) }
    mtext(SETLAB[st], side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85)
  }
  mtext(sprintf("Antigen presentation, TEAD4 (orange diamond) vs the other 7 KOs (grey), NES per source, %s", DLAB), outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1.0)
  mtext(sprintf("each KO, %s, same source; * = TEAD4 padj<0.05", DLAB), outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.55, col = INK2)
})

## ---------------------------------------------------------------- fig4F: curated AP genes, TEAD4 KO vs NTC statistic per source (+ per-line DESeq2 log2FC)
core <- c(MHC1, APM, MHC2)
GM <- sapply(SRC, function(s) { x <- stats[[s]][, KO]; x[core] }); rownames(GM) <- core
getv <- function(M) { v <- setNames(rep(NA_real_, length(core)), core); ok <- intersect(core, rownames(M)); v[ok] <- SGN * M[ok, KO]; v }
LM <- cbind(Screen = getv(ss$diff0), sapply(LINES, function(cl) getv(bs[[cl]]$lfc0))); colnames(LM) <- short
grp <- ifelse(core %in% MHC1, "MHC-I", ifelse(core %in% APM, "APM", "MHC-II"))
write.csv(data.frame(gene = core, group = grp, setNames(as.data.frame(GM), paste0("stat_", short)), setNames(as.data.frame(LM), paste0("lfc_", short))), file.path(OUT, paste0("spec50_AP_core_genes_TEAD4", DSUF, ".csv")), row.names = FALSE)
cat(sprintf("\n=== curated AP gene sets: mean TEAD4 %s stat per source (t or Wald), genes detected ===\n", DLAB))
for (g in c("MHC-I", "APM", "MHC-II")) cat(sprintf("%-7s %s\n", g, paste(sprintf("%s %+.2f (n=%d)", short, colMeans(GM[grp == g, , drop = FALSE], na.rm = TRUE), colSums(is.finite(GM[grp == g, , drop = FALSE]))), collapse = " | ")))
save_fig(paste0("fig4F_TEAD4_AP_core_genes", DSUF), 7.5, 9.5, function() {
  n <- length(core); vmax <- 4; Z <- pmax(pmin(GM, vmax), -vmax); par(mar = c(4.5, 6, 3.6, 1), family = "sans")
  image(seq_along(SRC), seq_len(n), t(Z)[, n:1], col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  na_i <- which(is.na(GM), arr.ind = TRUE); if (nrow(na_i)) rect(na_i[, 2] - 0.5, n - na_i[, 1] + 0.5, na_i[, 2] + 0.5, n - na_i[, 1] + 1.5, col = "#efefef", border = NA)
  axis(1, at = seq_along(SRC), labels = FALSE); text(seq_along(SRC), 0.2, short[SRC], srt = 45, adj = 1, xpd = NA, cex = 0.8)
  text(0.45, n:1, core, adj = 1, xpd = NA, cex = 0.6, col = c(`MHC-I` = "#eb6834", APM = "#4a3aa7", `MHC-II` = "#1baf7a")[grp])
  gb <- which(grp[-1] != grp[-n]); for (g in gb) abline(h = n - g + 0.5, col = "white", lwd = 2.5); abline(v = 1.5, col = "white", lwd = 3); box(col = MUTED, lwd = 0.6)
  mtext(sprintf("Curated antigen-presentation genes, TEAD4 %s", DLAB), side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext(sprintf("colour = %s statistic (screen: Welch t; bulk: Wald), capped at +-4; red = higher in %s; grey = gene not tested", DLAB, if (SGN < 0) "NTC" else "KO"), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
})
