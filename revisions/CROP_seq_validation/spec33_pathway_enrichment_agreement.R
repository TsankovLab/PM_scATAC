# Spec 33: for each KO, take the top 100 DE genes (by |DESeq2 Wald stat| vs NTC) IN EACH CELL LINE SEPARATELY, run GO Biological
# Process enrichment (clusterProfiler::enrichGO, org.Hs.eg.db + GO.db, hypergeometric test, universe = all genes tested in that line),
# and check whether the enriched GO terms AGREE across the 4 cell lines for the same KO (Jaccard overlap of significant-term sets).
# needs: spec1_prepare_stats.R (spec_bulk_stats.rds)
# outputs: spec33_top100_genes_per_KO_per_line.csv, spec33_GO_enrichment_all.csv, spec33_KO_pathway_agreement_summary.csv,
#          fig2N_pathway_enrichment_agreement.{png,pdf} (top GO term per KO per line), fig2O_pathway_overlap_all_combinations.{png,pdf}
#          (one combined heatmap, KO x line-pair, all 6 pairwise combinations, diagonal removed)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(GO.db) })

bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))

## ---------------------------------------------------------------- top 100 genes per KO per line (ranked by |Wald stat| vs NTC)
top100 <- list()
for (cl in LINES) {
  st <- bs[[cl]]$stat0
  for (k in ORDER) {
    if (!(k %in% colnames(st))) next
    x <- st[, k]; x <- x[is.finite(x)]
    g <- names(sort(abs(x), decreasing = TRUE))[1:min(100, length(x))]
    top100[[length(top100) + 1]] <- data.frame(line = cl, KO = k, gene = g, stat = x[g])
  }
}
top100 <- do.call(rbind, top100); write.csv(top100, file.path(OUT, "spec33_top100_genes_per_KO_per_line.csv"), row.names = FALSE)

## ---------------------------------------------------------------- GO BP enrichment, per KO per line
enr <- list()
for (cl in LINES) {
  universe <- bs[[cl]]$genes
  for (k in ORDER) {
    g <- top100$gene[top100$line == cl & top100$KO == k]
    if (!length(g)) next
    e <- tryCatch(enrichGO(gene = g, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                            pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500), error = function(e) NULL)
    if (is.null(e) || nrow(e@result) == 0) next
    d <- e@result; d$line <- cl; d$KO <- k
    enr[[length(enr) + 1]] <- d
  }
}
enr <- do.call(rbind, enr); rownames(enr) <- NULL
write.csv(enr, file.path(OUT, "spec33_GO_enrichment_all.csv"), row.names = FALSE)
cat(sprintf("GO BP enrichment tested: %d KO x line combos with >=1 result\n", length(unique(paste(enr$KO, enr$line)))))

## ---------------------------------------------------------------- per-KO agreement: pairwise Jaccard of SIGNIFICANT (padj<0.25) GO terms across lines
sig <- enr[enr$p.adjust < 0.25, ]
pairs <- combn(LINES, 2, simplify = FALSE)
agree <- list(); jac_by_ko <- list()
for (k in ORDER) {
  sets <- lapply(LINES, function(cl) unique(sig$ID[sig$KO == k & sig$line == cl])); names(sets) <- LINES
  n_lines_any <- sum(sapply(sets, length) > 0)
  JM <- matrix(NA, length(LINES), length(LINES), dimnames = list(LINES, LINES)); diag(JM) <- 1
  for (p in pairs) { a <- sets[[p[1]]]; b <- sets[[p[2]]]; j <- if (!length(a) && !length(b)) NA else length(intersect(a, b)) / length(union(a, b)); JM[p[1], p[2]] <- j; JM[p[2], p[1]] <- j }
  jac_by_ko[[k]] <- JM
  jac <- JM[lower.tri(JM)]
  shared_all_terms <- table(unlist(sets)); shared2plus <- names(shared_all_terms)[shared_all_terms >= 2]
  desc <- setNames(sig$Description, sig$ID)[shared2plus]
  agree[[k]] <- data.frame(KO = k, n_lines_with_hit = n_lines_any, mean_jaccard = mean(jac, na.rm = TRUE),
                            n_terms_shared_by_2plus_lines = length(shared2plus),
                            shared_terms = paste(unique(desc), collapse = "; "))
}
agree <- do.call(rbind, agree); agree <- agree[order(-agree$mean_jaccard), ]
write.csv(agree, file.path(OUT, "spec33_KO_pathway_agreement_summary.csv"), row.names = FALSE)
options(width = 200); cat("\n=== per-KO cross-line GO BP agreement (padj<0.25, top100 genes by |stat|) ===\n")
print(agree[, c("KO", "n_lines_with_hit", "mean_jaccard", "n_terms_shared_by_2plus_lines")], row.names = FALSE)
cat("\nshared GO terms per KO (first few):\n")
for (i in seq_len(nrow(agree))) if (nchar(agree$shared_terms[i])) cat(sprintf("  %-8s %s\n", agree$KO[i], substr(agree$shared_terms[i], 1, 160)))

## ---------------------------------------------------------------- fig2N: top GO term per KO per line
save_fig("fig2N_pathway_enrichment_agreement", 12, 7.6, function() {
  ko_ord <- agree$KO
  best_per_cell <- do.call(rbind, lapply(ko_ord, function(k) do.call(rbind, lapply(LINES, function(cl) {
    d <- enr[enr$KO == k & enr$line == cl, ]; if (!nrow(d)) return(data.frame(KO = k, line = cl, best = NA, nlp = NA))
    d <- d[order(d$p.adjust), ][1, ]; data.frame(KO = k, line = cl, best = d$Description, nlp = -log10(d$p.adjust))
  }))))
  M <- matrix(NA, length(ko_ord), length(LINES), dimnames = list(ko_ord, LINES)); LAB <- M
  for (i in seq_len(nrow(best_per_cell))) { M[best_per_cell$KO[i], best_per_cell$line[i]] <- best_per_cell$nlp[i]; LAB[best_per_cell$KO[i], best_per_cell$line[i]] <- best_per_cell$best[i] }
  vmax <- max(M, na.rm = TRUE); Mz <- M; Mz[] <- pmin(M, vmax)
  par(mar = c(9, 8, 3.6, 1), family = "sans")
  image(seq_len(ncol(M)), seq_len(nrow(M)), t(Mz)[, nrow(M):1, drop = FALSE], col = DIV[51:101], zlim = c(0, vmax), axes = FALSE, xlab = "", ylab = "")
  na_i <- which(is.na(M), arr.ind = TRUE); if (nrow(na_i)) rect(na_i[, 1] - 0.5, nrow(M) - na_i[, 2] + 0.5, na_i[, 1] + 0.5, nrow(M) - na_i[, 2] + 1.5, col = "#efefef", border = NA)
  axis(1, at = seq_len(ncol(M)), labels = colnames(M), las = 2, tick = FALSE, cex.axis = 0.8)
  for (ri in seq_len(nrow(M))) text(0.35, nrow(M) - ri + 1, rownames(M)[ri], adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[rownames(M)[ri]])
  for (ri in seq_len(nrow(M))) for (ci in seq_len(ncol(M))) {
    lab <- LAB[ri, ci]; if (is.na(lab)) next
    words <- strsplit(lab, " ")[[1]]
    wrapped <- paste(sapply(split(words, ceiling(seq_along(words) / 3)), paste, collapse = " "), collapse = "\n")
    text(ci, nrow(M) - ri + 1, wrapped, cex = 0.38, col = if (M[ri, ci] > vmax * 0.55) "white" else INK)
  }
  box(col = MUTED, lwd = 0.6)
  mtext("Top GO Biological Process term per KO per line (top100 genes by |Wald stat| vs NTC)", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("colour = -log10(padj) of that line's single best-matching GO BP term; grey = no term tested/found; text in each cell = that term's name", side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  hin <- par("din")[2]
  op <- par(fig = c(0.03, 0.30, 0, 0.9 / hin), new = TRUE, oma = c(0, 0, 0, 0), mar = c(2.0, 0.5, 1.0, 0.5), cex = 1); on.exit(par(op))
  image(seq(0, vmax, length.out = 51), 1, matrix(seq(0, vmax, length.out = 51), ncol = 1), col = DIV[51:101], zlim = c(0, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(1, cex.axis = 0.55, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext("-log10(padj), best GO BP term", side = 1, line = 1.2, cex = 0.55, col = INK2); box(lwd = 0.5)
})

## ---------------------------------------------------------------- fig2O: ONE combined heatmap, KO x line-pair (6 pairwise combinations, diagonal removed)
short <- function(l) sub("NCI-", "", sub("-211H", "", l))
ko_ord <- agree$KO
pair_lab <- sapply(pairs, function(p) paste(short(p), collapse = "-"))
MC <- matrix(NA, length(ko_ord), length(pairs), dimnames = list(ko_ord, pair_lab))
for (k in ko_ord) for (pi in seq_along(pairs)) MC[k, pi] <- jac_by_ko[[k]][pairs[[pi]][1], pairs[[pi]][2]]

save_fig("fig2O_pathway_overlap_all_combinations", 8.5, 6.6, function() {
  par(mar = c(9, 7, 3.6, 1), family = "sans")
  vmax <- max(MC, na.rm = TRUE)
  image(seq_len(ncol(MC)), seq_len(nrow(MC)), t(MC)[, nrow(MC):1, drop = FALSE], col = colorRampPalette(c("#f7f7f7", "#1baf7a"))(100), zlim = c(0, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_len(ncol(MC)), labels = colnames(MC), las = 2, tick = FALSE, cex.axis = 0.8)
  for (ri in seq_len(nrow(MC))) text(0.3, nrow(MC) - ri + 1, rownames(MC)[ri], adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[rownames(MC)[ri]])
  for (ri in seq_len(nrow(MC))) for (ci in seq_len(ncol(MC))) { v <- MC[ri, ci]; text(ci, nrow(MC) - ri + 1, sprintf("%.2f", v), cex = 0.68, col = if (v > vmax * 0.55) "white" else INK) }
  box(col = MUTED, lwd = 0.6)
  mtext("Cross-line overlap of significant GO Biological Process terms, per KO x line-pair", side = 3, adj = 0, line = 2.2, font = 2, cex = 0.95)
  mtext("Jaccard index of significant (padj<0.25) GO BP term sets between each pair of lines; diagonal (self-overlap = 1) removed; KOs ordered by mean pairwise Jaccard",
        side = 3, adj = 0, line = 1.0, cex = 0.55, col = INK2)
  hin <- par("din")[2]
  op <- par(fig = c(0.03, 0.30, 0, 0.9 / hin), new = TRUE, oma = c(0, 0, 0, 0), mar = c(2.0, 0.5, 1.0, 0.5), cex = 1); on.exit(par(op))
  image(seq(0, vmax, length.out = 51), 1, matrix(seq(0, vmax, length.out = 51), ncol = 1), col = colorRampPalette(c("#f7f7f7", "#1baf7a"))(100), zlim = c(0, vmax), axes = FALSE, xlab = "", ylab = "")
  # (colorbar drawn in its own small fig region below the main plot's bottom margin -- kept separate from the rotated column labels)
  axis(1, cex.axis = 0.55, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext("Jaccard index", side = 1, line = 1.2, cex = 0.55, col = INK2); box(lwd = 0.5)
})

## ---------------------------------------------------------------- fig2P: full (KO x line) x (KO x line) Jaccard overlap matrix, every unit vs every other unit (not just same-KO pairs)
units <- expand.grid(line = LINES, KO = ko_ord, stringsAsFactors = FALSE)
units <- units[order(match(units$KO, ko_ord), match(units$line, LINES)), ]
units$id <- paste0(units$KO, "_", short(units$line)); n <- nrow(units)
sets_all <- lapply(seq_len(n), function(i) unique(sig$ID[sig$KO == units$KO[i] & sig$line == units$line[i]]))
J <- matrix(NA, n, n, dimnames = list(units$id, units$id))
for (i in seq_len(n)) for (j in seq_len(n)) { a <- sets_all[[i]]; b <- sets_all[[j]]
  J[i, j] <- if (i == j) NA else if (!length(a) && !length(b)) NA else length(intersect(a, b)) / length(union(a, b)) }
write.csv(as.data.frame(J), file.path(OUT, "spec33_full_KOline_jaccard_matrix.csv"))

save_fig("fig2P_pathway_overlap_full_KO_line_matrix", 12.5, 12, function() {
  par(mar = c(9, 9, 1, 1), oma = c(0, 0, 3.6, 0), family = "sans")
  vmax <- max(J, na.rm = TRUE)
  image(seq_len(n), seq_len(n), t(J)[, n:1, drop = FALSE], col = colorRampPalette(c("#f7f7f7", "#1baf7a"))(100), zlim = c(0, vmax), axes = FALSE, xlab = "", ylab = "")
  na_i <- which(is.na(J), arr.ind = TRUE); if (nrow(na_i)) rect(na_i[, 1] - 0.5, n - na_i[, 2] + 0.5, na_i[, 1] + 0.5, n - na_i[, 2] + 1.5, col = "#efefef", border = NA)
  edges <- which(units$KO[-1] != units$KO[-n]); for (e in edges) { abline(v = e + 0.5, col = "white", lwd = 2.2); abline(h = n - e + 0.5, col = "white", lwd = 2.2) }
  axis(1, at = seq_len(n), labels = units$id, las = 2, tick = FALSE, cex.axis = 0.5, col.axis = INK2)
  axis(2, at = n:1, labels = units$id, las = 1, tick = FALSE, cex.axis = 0.5, col.axis = INK2)
  for (ki in seq_along(ko_ord)) { at <- mean(which(units$KO == ko_ord[ki])); mtext(ko_ord[ki], side = 3, at = at, line = 0.3, font = 2, cex = 0.62, col = KO_COL[ko_ord[ki]]) }
  box(col = MUTED, lwd = 0.6)
  mtext("Pathway overlap between every KO x line unit", outer = TRUE, side = 3, adj = 0, line = 2.1, font = 2, cex = 0.95)
  mtext("32 x 32 Jaccard matrix of significant GO BP terms (padj<0.25); diagonal removed; white lines = KO block boundaries (each block = one KO's 4 lines)", outer = TRUE, side = 3, adj = 0, line = 0.9, cex = 0.55, col = INK2)
  hin <- par("din")[2]
  op <- par(fig = c(0.03, 0.20, 0, 0.6 / hin), new = TRUE, oma = c(0, 0, 0, 0), mar = c(2.0, 0.5, 1.0, 0.5), cex = 1); on.exit(par(op))
  image(seq(0, vmax, length.out = 51), 1, matrix(seq(0, vmax, length.out = 51), ncol = 1), col = colorRampPalette(c("#f7f7f7", "#1baf7a"))(100), zlim = c(0, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(1, cex.axis = 0.55, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext("Jaccard index (GO BP terms, padj<0.25)", side = 1, line = 1.2, cex = 0.55, col = INK2); box(lwd = 0.5)
})
