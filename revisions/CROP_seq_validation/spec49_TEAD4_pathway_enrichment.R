# Spec 49 (comparison direction: NTC - KO by default, DIFF_DIR=KO_NTC for the reverse; UP = higher in NTC than TEAD4 KO): pathway enrichment of the top TEAD4-KO-vs-NTC genes in EACH cell line (bulk DESeq2 Wald stat, 4 lines) and in the CROP-seq screen
# (Welch t, NCI-H2052), split into UP (KO>NTC) and DOWN (KO<NTC) genes; then agreement (gene-level and pathway-level overlap) and identity
# (which terms are shared) across sources. GO Biological Process (clusterProfiler::enrichGO) + MSigDB Hallmark (enricher), universe = genes tested
# in that source. Env: TOPN (default 150 genes per direction), PADJ (default 0.05 significance for term overlap).
# needs: spec1_prepare_stats.R (spec_bulk_stats.rds, spec_screen_stats.rds)
# outputs: spec49_TEAD4_top_genes.csv, spec49_TEAD4_enrichment_{GO,Hallmark}.csv, spec49_TEAD4_pathway_agreement.csv,
#          fig4A_TEAD4_GO_dotplot, fig4B_TEAD4_Hallmark_dotplot, fig4C_TEAD4_gene_and_pathway_overlap
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr) })
TOPN <- as.integer(Sys.getenv("TOPN", "150")); PADJ <- as.numeric(Sys.getenv("PADJ", "0.05")); KO <- "TEAD4"
DIFF_DIR <- Sys.getenv("DIFF_DIR", "NTC_KO"); SGN <- if (DIFF_DIR == "NTC_KO") -1 else 1; DLAB <- if (DIFF_DIR == "NTC_KO") "NTC-KO" else "KO-NTC"; DSUF <- if (DIFF_DIR == "NTC_KO") "" else "_KOminusNTC"   # default NTC - KO (positive = higher in NTC, i.e. lower after KO)

ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))
SRC <- c("CROP-seq (H2052)", paste0(LINES, " bulk")); short <- c("Screen", "H2052", "MSTO", "H2452", "H28"); names(short) <- SRC
stat <- c(list(SGN * ss$t0[, KO]), lapply(LINES, function(cl) SGN * bs[[cl]]$stat0[, KO])); names(stat) <- SRC
univ <- c(list(ss$genes), lapply(LINES, function(cl) bs[[cl]]$genes)); names(univ) <- SRC

## ---------------------------------------------------------------- top genes per source and direction
top <- do.call(rbind, lapply(SRC, function(s) { x <- stat[[s]]; x <- x[is.finite(x)]
  rbind(data.frame(source = s, dir = "UP", gene = names(sort(x, decreasing = TRUE))[1:TOPN]),
        data.frame(source = s, dir = "DOWN", gene = names(sort(x))[1:TOPN])) }))
top$stat <- mapply(function(s, g) stat[[s]][g], top$source, top$gene)
write.csv(top, file.path(OUT, paste0("spec49_TEAD4_top_genes", DSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- enrichment
hm <- msigdbr(species = "Homo sapiens", category = "H"); hm2g <- hm[, c("gs_name", "gene_symbol")]; hm2g$gs_name <- sub("HALLMARK_", "", hm2g$gs_name)
run_enr <- function(kind) do.call(rbind, lapply(SRC, function(s) do.call(rbind, lapply(c("UP", "DOWN"), function(d) {
  g <- top$gene[top$source == s & top$dir == d]
  e <- tryCatch(if (kind == "GO") enrichGO(g, universe = univ[[s]], OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
                else enricher(g, universe = univ[[s]], TERM2GENE = hm2g, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500), error = function(e) NULL)
  if (is.null(e) || !nrow(e@result)) return(NULL)
  r <- e@result; r$source <- s; r$dir <- d; r$k <- as.numeric(sub("/.*", "", r$GeneRatio)); r$gene_ratio <- r$k / as.numeric(sub(".*/", "", r$GeneRatio)); r }))))
GO <- run_enr("GO"); HM <- run_enr("HM")
write.csv(GO, file.path(OUT, paste0("spec49_TEAD4_enrichment_GO", DSUF, ".csv")), row.names = FALSE); write.csv(HM, file.path(OUT, paste0("spec49_TEAD4_enrichment_Hallmark", DSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- agreement: gene-level and pathway-level Jaccard (sig terms) between every source pair, per direction
jac <- function(a, b) if (!length(a) && !length(b)) NA else length(intersect(a, b)) / length(union(a, b))
pairs <- combn(SRC, 2, simplify = FALSE); rows <- list()
for (d in c("UP", "DOWN")) for (p in pairs) {
  gj <- jac(top$gene[top$source == p[1] & top$dir == d], top$gene[top$source == p[2] & top$dir == d])
  gn <- length(intersect(top$gene[top$source == p[1] & top$dir == d], top$gene[top$source == p[2] & top$dir == d]))
  sg <- function(R, s) unique(R$ID[R$source == s & R$dir == d & R$p.adjust < PADJ])
  rows[[length(rows) + 1]] <- data.frame(dir = d, a = short[p[1]], b = short[p[2]], genes_shared = gn, gene_jaccard = gj,
    GO_jaccard = jac(sg(GO, p[1]), sg(GO, p[2])), GO_shared = length(intersect(sg(GO, p[1]), sg(GO, p[2]))),
    HM_jaccard = jac(sg(HM, p[1]), sg(HM, p[2])), HM_shared = length(intersect(sg(HM, p[1]), sg(HM, p[2]))))
}
agree <- do.call(rbind, rows); write.csv(agree, file.path(OUT, paste0("spec49_TEAD4_pathway_agreement", DSUF, ".csv")), row.names = FALSE)
options(width = 200)
cat(sprintf("\n=== TEAD4 %s: top%d UP / DOWN genes per source; significant = padj<%g ===\n", DLAB, TOPN, PADJ))
cat("# significant terms per source:\n")
print(do.call(rbind, lapply(SRC, function(s) data.frame(source = s, GO_up = sum(GO$p.adjust < PADJ & GO$source == s & GO$dir == "UP"), GO_down = sum(GO$p.adjust < PADJ & GO$source == s & GO$dir == "DOWN"),
  HM_up = sum(HM$p.adjust < PADJ & HM$source == s & HM$dir == "UP"), HM_down = sum(HM$p.adjust < PADJ & HM$source == s & HM$dir == "DOWN")))), row.names = FALSE)
cat("\n# pairwise overlap:\n"); print(format(agree, digits = 2), row.names = FALSE)
cat("\n# mean pairwise Jaccard among the 4 bulk lines vs screen-to-bulk:\n")
agree$type <- ifelse(agree$a == "Screen", "screen-bulk", "bulk-bulk")
print(aggregate(cbind(gene_jaccard, GO_jaccard, HM_jaccard) ~ dir + type, agree, mean), row.names = FALSE)
for (d in c("UP", "DOWN")) { cat(sprintf("\n# Hallmark terms significant (padj<%g) in >=2 sources, %s:\n", PADJ, d))
  s <- HM[HM$dir == d & HM$p.adjust < PADJ, ]; tb <- sort(table(s$ID), decreasing = TRUE); tb <- tb[tb >= 2]
  for (t in names(tb)) cat(sprintf("  %-32s %d/5: %s\n", t, tb[t], paste(short[s$source[s$ID == t]], collapse = ",")))
  cat(sprintf("# GO terms significant in >=3 sources, %s:\n", d))
  s <- GO[GO$dir == d & GO$p.adjust < PADJ, ]; tb <- sort(table(s$Description), decreasing = TRUE); tb <- tb[tb >= 3]
  for (t in head(names(tb), 25)) cat(sprintf("  %-60s %d/5: %s\n", t, tb[t], paste(short[s$source[s$Description == t]], collapse = ",")))
}

## ---------------------------------------------------------------- dotplots: terms selected by the best rank across sources, columns = source, panels = UP / DOWN
dotplot_fig <- function(R, name, title, nterm, wrap, lab_col) {
  R$lab <- R[[lab_col]]
  save_fig(name, 12.5, 4.2 + nterm * 0.3, function() {
    layout(matrix(1:2, nrow = 1)); par(oma = c(0, 0, 3.2, 0), family = "sans")
    vmax <- max(-log10(pmax(R$p.adjust, 1e-30)), na.rm = TRUE)
    for (d in c("UP", "DOWN")) {
      Rd <- R[R$dir == d, ]
      # terms: ranked by number of sources significant, then by best padj; keep top nterm
      agg <- do.call(rbind, lapply(split(Rd, Rd$lab), function(x) data.frame(lab = x$lab[1], nsig = sum(x$p.adjust < PADJ), best = min(x$p.adjust))))
      agg <- agg[order(-agg$nsig, agg$best), ]; terms <- head(agg$lab[agg$nsig > 0], nterm)
      par(mar = c(4, 20, 2.2, 1.5))
      plot(NA, xlim = c(0.5, length(SRC) + 0.5), ylim = c(0.5, max(length(terms), 1) + 0.5), axes = FALSE, xlab = "", ylab = "")
      if (!length(terms)) { title(paste(d, "- no significant terms")); next }
      for (ti in seq_along(terms)) for (si in seq_along(SRC)) {
        x <- Rd[Rd$lab == terms[ti] & Rd$source == SRC[si], ]; if (!nrow(x)) next
        nl <- -log10(x$p.adjust[1]); sigf <- x$p.adjust[1] < PADJ
        points(si, length(terms) - ti + 1, pch = 21, cex = 0.6 + 4 * min(x$gene_ratio[1], 0.3),
               bg = if (sigf) DIV[round(51 + 50 * min(nl / vmax, 1))] else "grey92", col = if (sigf) INK else "grey70", lwd = if (sigf) 1.2 else 0.5)
      }
      axis(1, at = seq_along(SRC), labels = FALSE, tick = FALSE); text(seq_along(SRC), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), short[SRC], srt = 45, adj = 1, xpd = NA, cex = 0.8)
      axis(2, at = length(terms):1, labels = strtrim(terms, 60), las = 1, tick = FALSE, cex.axis = 0.62)
      abline(v = 1.5, col = "grey60", lty = 3); box(col = MUTED, lwd = 0.6)
      mtext(sprintf("%s: %s", d, if ((d == "UP") == (SGN < 0)) "lower in TEAD4 KO" else "higher in TEAD4 KO"), side = 3, adj = 0, line = 0.3, font = 2, cex = 0.9, col = if (d == "UP") "#eb6834" else "#2a78d6")
    }
    mtext(title, outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1.0)
    mtext(sprintf("dot colour = -log10(padj), size = gene ratio, filled+black outline = padj<%g, grey = not significant; rows sorted by number of sources where significant", PADJ), outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.55, col = INK2)
  })
}
dotplot_fig(GO, paste0("fig4A_TEAD4_GO_dotplot", DSUF), sprintf("GO Biological Process: TEAD4 %s, screen and each bulk line", DLAB), 25, 60, "Description")
dotplot_fig(HM, paste0("fig4B_TEAD4_Hallmark_dotplot", DSUF), sprintf("MSigDB Hallmark: TEAD4 %s, screen and each bulk line", DLAB), 20, 40, "ID")

## ---------------------------------------------------------------- fig4C: source x source heatmaps (gene Jaccard, GO Jaccard), UP and DOWN
save_fig(paste0("fig4C_TEAD4_gene_and_pathway_overlap", DSUF), 13, 7.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (d in c("UP", "DOWN")) for (m in c("gene_jaccard", "GO_jaccard", "HM_jaccard")) {
    M <- matrix(NA, 5, 5, dimnames = list(short, short)); diag(M) <- 1
    for (i in which(agree$dir == d)) { M[agree$a[i], agree$b[i]] <- agree[[m]][i]; M[agree$b[i], agree$a[i]] <- agree[[m]][i] }
    vmax <- max(M[lower.tri(M)], na.rm = TRUE); par(mar = c(3.5, 4, 2.2, 1))
    image(1:5, 1:5, t(M)[, 5:1], col = colorRampPalette(c("#f7f7f7", "#1baf7a"))(100), zlim = c(0, max(vmax, 0.01)), axes = FALSE, xlab = "", ylab = "")
    for (i in 1:5) for (j in 1:5) if (i != j) text(j, 6 - i, if (is.na(M[i, j])) "NA" else sprintf("%.2f", M[i, j]), cex = 0.75, col = INK)
    axis(1, 1:5, short, tick = FALSE, cex.axis = 0.8); axis(2, 5:1, short, tick = FALSE, las = 1, cex.axis = 0.8); box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s: %s", d, c(gene_jaccard = "top-gene Jaccard", GO_jaccard = "GO BP term Jaccard", HM_jaccard = "Hallmark Jaccard")[m]), side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85)
  }
  mtext(sprintf("Agreement of TEAD4 (%s) responses between the screen and each bulk line (top%d genes per direction; terms padj<%g)", DLAB, TOPN, PADJ), outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.0)
})
