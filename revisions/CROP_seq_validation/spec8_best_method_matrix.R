# Spec 8: correlation matrix bulk KO (rows, 4 lines) x CROP-seq KO (columns) for the best variant F4:
#   cell-cycle-adjusted test statistics (screen: Welch t after cc regression; bulk: DESeq2 Wald, design ~ cc + condition), all expressed genes,
#   each dataset centred across its KOs (removes the shared NTC / global response), Pearson correlation.
# Columns: the 6 screen KOs with >= 50 cells (used for the test), then the low-cell-count KOs (TEAD2 n=20, HMGA1 n=13), then the screen-only KOs.
# needs: spec1_prepare_stats.R
# outputs: spec_bestmethod_matrix.csv, fig10_best_method_matrix.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))

ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))
POWERED <- intersect(ORDER, names(ss$n)[ss$n >= 50]); LOWN <- setdiff(ORDER, POWERED); COLS <- c(POWERED, LOWN, EXTRA)
rows_lines <- c("NCI-H2052", "MSTO-211H", "NCI-H2452", "NCI-H28")          # screen line first
genes <- Reduce(intersect, c(list(ss$genes), lapply(bs[rows_lines], `[[`, "genes"))); cat("genes:", length(genes), "\n")
S <- ss$tcc[genes, COLS]; S <- prep_dataset(S, TRUE, FALSE)
rows <- list(); lab <- list()
for (cl in rows_lines) { m <- bs[[cl]]$statcc[genes, ]; m[!is.finite(m)] <- 0; m <- prep_dataset(m, TRUE, FALSE)
  for (ko in c(POWERED, LOWN)) if (ko %in% colnames(m)) { rows[[length(rows) + 1]] <- sapply(COLS, function(cc) cor(m[, ko], S[, cc])); lab[[length(lab) + 1]] <- c(cl, ko) } }
C <- do.call(rbind, rows); lab <- setNames(as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE), c("line", "ko")); rownames(C) <- paste(lab$line, lab$ko)
write.csv(C, file.path(OUT, "spec_bestmethod_matrix.csv"))

# statistics on the 6 well-powered screen KOs (exact permutation over all 6! KO assignments), for 3 lines and all 4
st <- lapply(list(`3 lines` = lab$line %in% THREE, `4 lines` = rep(TRUE, nrow(lab)), `NCI-H2052 only` = lab$line == "NCI-H2052"), function(sel) {
  s <- sel & lab$ko %in% POWERED; evaluate_specificity(C[s, POWERED, drop = FALSE], lab$ko[s], cols = POWERED) })
stt <- do.call(rbind, st); print(format(stt, digits = 3))
rank_pow <- sapply(seq_len(nrow(lab)), function(i) if (lab$ko[i] %in% POWERED) sum(C[i, POWERED] > C[i, lab$ko[i]]) + 1 else NA)

## ---------------------------------------------------------------- figure
n <- nrow(C); k <- ncol(C); vmax <- ceiling(quantile(abs(C), 0.99, names = FALSE) * 20) / 20
save_fig("fig10_best_method_matrix", 12.5, 12, function() {
  par(mar = c(5, 11, 14, 8.5), family = "sans")
  image(seq_len(k), seq_len(n), heat_z(pmax(pmin(C, vmax), -vmax)), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  ncell <- ss$n[COLS]
  pw <- COLS %in% POWERED
  axis(3, at = which(pw), labels = sprintf("%s (%d)", COLS, ss$n[COLS])[pw], las = 2, tick = FALSE, cex.axis = 0.95, line = -0.4, col.axis = INK, font.axis = 2)
  axis(3, at = which(!pw), labels = sprintf("%s (%d)", COLS, ss$n[COLS])[!pw], las = 2, tick = FALSE, cex.axis = 0.95, line = -0.4, col.axis = MUTED, font.axis = 1)
  axis(2, at = n:1, labels = lab$ko, las = 1, tick = FALSE, cex.axis = 0.8, line = -0.4)
  for (i in seq_len(n)) for (j in seq_len(k)) text(j, n - i + 1, sub("^([+-])0\\.", "\\1.", sprintf("%+.2f", C[i, j])), cex = 0.62, col = if (abs(C[i, j]) < 0.6 * vmax) INK2 else "white")
  for (i in seq_len(n)) { j <- match(lab$ko[i], COLS); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = INK, lwd = 1.8) }
  abline(v = c(length(POWERED), length(POWERED) + length(LOWN)) + 0.5, col = "white", lwd = 5)
  edges <- which(lab$line[-1] != lab$line[-n]); abline(h = n - edges + 0.5, col = "white", lwd = 5); bounds <- c(0, edges, n); par(xpd = NA)
  for (b in seq_along(rows_lines)) { y1 <- n - bounds[b + 1] + 0.5; y2 <- n - bounds[b] + 0.5
    rect(-1.15, y1, -0.85, y2, col = LINE_COL[[rows_lines[b]]], border = NA); text(-1.6, (y1 + y2) / 2, rows_lines[b], srt = 90, font = 2, cex = 0.95, col = LINE_COL[[rows_lines[b]]]) }
  text(k + 0.75, n:1, ifelse(is.na(rank_pow), "", sprintf("rank %d/6", rank_pow)), adj = 0, cex = 0.7, col = ifelse(!is.na(rank_pow) & rank_pow <= 2, "#b3400f", MUTED), font = ifelse(!is.na(rank_pow) & rank_pow <= 2, 2, 1))
  text(k + 0.75, n + 0.9, "matched KO rank\n(of 6 powered)", adj = 0, cex = 0.62, col = INK2)
  mtext("Bulk KO vs CROP-seq KO: best variant (cell-cycle adjusted test statistics, centred across KOs)", side = 3, adj = 0, line = 11.8, font = 2, cex = 1.05, at = -1.9)
  mtext("all expressed genes; rows = bulk KO (RUN1) per cell line, columns = CROP-seq KO (cells in brackets); matched pair outlined; bold columns = screen KOs with >= 50 cells", side = 3, adj = 0, line = 10.6, cex = 0.62, col = INK2, at = -1.9)
  mtext(sprintf("exact permutation test on those 6 screen KOs: 3 lines matched z = %.2f (p = %.3f) | 4 lines z = %.2f (p = %.3f) | NCI-H2052 only z = %.2f (p = %.3f)",
                stt$z_matched[1], stt$perm_p[1], stt$z_matched[2], stt$perm_p[2], stt$z_matched[3], stt$perm_p[3]), side = 3, adj = 0, line = 9.6, cex = 0.62, col = INK2, at = -1.9)
  mtext(sprintf("mean r matched %+.3f vs other KOs %+.3f (3 lines); r is small because only the KO-specific part of each profile is compared", stt$matched_r[1], stt$other_r[1]), side = 3, adj = 0, line = 8.6, cex = 0.62, col = INK2, at = -1.9)
  par(xpd = FALSE); draw_colorbar(c(0.35, 0.65), c(-vmax, vmax), "Pearson r between bulk KO profile and CROP-seq KO profile")
})
