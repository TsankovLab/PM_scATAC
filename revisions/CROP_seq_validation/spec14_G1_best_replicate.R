# Spec 14: CROP-seq screen restricted to G1 cells with the remaining cells still cell-cycle regressed (version D of spec13), versus cc-adjusted bulk, using
#          (a) replicates averaged, (b) ONLY the bulk replicate best matching the screen KO, (c) control: best replicate picked separately for every cell of the matrix.
# Baseline: the same procedure on the all-cell, cc-regressed screen (version A). Bulk is cc-adjusted (replicate-level pooled-variance t, or cc-regressed log difference).
# Variants: features (all expressed genes | screen variable genes) x statistic (t | log difference, RMS-scaled) x centred across KOs (no | yes) x Pearson.
# needs: spec1_prepare_stats.R
# outputs: spec14_results.csv, fig18_G1_best_replicate.{png,pdf}, fig18b_G1_best_replicate_matrix.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
ctx <- load_spec_context(LINES); ss <- ctx$ss
srt <- load_screen(); md <- srt@meta.data; keep <- md$Phase == "G1"; Y <- GetAssayData(srt, layer = "data")[ss$genes, ]
G1r <- screen_ko_stats(Y[, keep], md$merged_call[keep], md$cc[keep], TRUE); n_g1 <- G1r$n; rm(Y, srt)
SV <- list(`A all cells, cc regressed` = list(t = ss$tcc, diff = ss$diffcc, n = ss$n), `D G1 only, cc regressed` = list(t = G1r$t, diff = G1r$diff, n = n_g1))
genes <- ctx$F_all; ALLC <- c(ORDER, EXTRA)

## ---------------------------------------------------------------- bulk replicate-level cc-adjusted profiles (genes x KO replicates) per line
REP <- list()
for (cl in LINES) { m <- ctx$meta[ctx$meta$Cell_Line == cl, ]; s <- rownames(m); y <- ctx$lg[genes, s, drop = FALSE]; cc <- ctx$ccb$cc[match(s, ctx$ccb$sample)]
  g <- factor(m$Gene_Control, levels = c("NTC", sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix"))); X <- model.matrix(~ cc + g); q <- qr(X)
  s2 <- pmax(colSums(qr.resid(q, t(y))^2) / (length(s) - q$rank), 1e-4); yadj <- y - outer(qr.coef(q, t(y))["cc", ], cc - mean(cc))
  ntc_s <- s[m$Gene_Control == "NTC"]; ko_s <- s[m$Gene_Control != "NTC"]; D <- yadj[, ko_s, drop = FALSE] - rowMeans(yadj[, ntc_s, drop = FALSE])
  Tm <- D / sqrt(s2 * (1 + 1 / length(ntc_s))); Tm[!is.finite(Tm)] <- 0; D[!is.finite(D)] <- 0
  REP[[cl]] <- list(t = Tm, diff = D, ko = setNames(as.character(m[ko_s, "Gene_Control"]), ko_s)) }

## ---------------------------------------------------------------- rows (line x TF in ORDER) with replicate lists
rows <- do.call(rbind, lapply(LINES, function(cl) { k <- REP[[cl]]$ko; k <- k[k %in% ORDER]; data.frame(line = cl, ko = unique(k)) })); rows$key <- paste(rows$line, rows$ko)
run_variant <- function(sv, feats, ty, cen) {
  S <- SV[[sv]][[ty]][feats, ALLC]; S <- prep_dataset(S, cen, ty == "diff")[, ALLC, drop = FALSE]
  avg <- matrix(NA, nrow(rows), length(ALLC)); best <- ctl <- avg; pick <- rep(NA_character_, nrow(rows)); r_all <- list()
  for (cl in LINES) { M <- prep_dataset(REP[[cl]][[ty]][feats, , drop = FALSE], cen, ty == "diff"); ko <- REP[[cl]]$ko
    for (i in which(rows$line == cl)) { ss_k <- names(ko)[ko == rows$ko[i]]; cr <- cor(M[, ss_k, drop = FALSE], S); mcol <- match(rows$ko[i], ORDER)     # replicates x screen KOs
      avg[i, ] <- cor(rowMeans(M[, ss_k, drop = FALSE]), S)[1, ]; b <- which.max(cr[, mcol]); best[i, ] <- cr[b, ]; ctl[i, ] <- apply(cr, 2, max); pick[i] <- ss_k[b] } }
  colnames(avg) <- colnames(best) <- colnames(ctl) <- ALLC; list(avg = avg, best = best, ctl = ctl, pick = pick) }
row_sets <- list(ALL3 = rows$line %in% THREE, ALL4 = rep(TRUE, nrow(rows)), H2052 = rows$line == "NCI-H2052")
res <- list(); store <- list()
for (sv in names(SV)) for (fs in c("all_genes", "var_genes")) for (ty in c("t", "diff")) for (cen in c(FALSE, TRUE)) {
  feats <- if (fs == "all_genes") genes else ctx$F_var; R <- run_variant(sv, feats, ty, cen); store[[paste(sv, fs, ty, cen)]] <- R
  for (cs in c("all8", "powered")) { cols <- if (cs == "all8") ORDER else intersect(ORDER, names(SV[[sv]]$n)[SV[[sv]]$n >= 50])
    for (rs in names(row_sets)) { sel <- row_sets[[rs]] & rows$ko %in% cols; if (sum(sel) < 3) next
      for (md_ in c("avg", "best", "ctl")) res[[length(res) + 1]] <- data.frame(screen = sv, features = fs, stat = ty, centred = cen, cols = cs, rows = rs, mode = md_, evaluate_specificity(R[[md_]][sel, , drop = FALSE], rows$ko[sel], cols = cols)) } } }
res <- do.call(rbind, res); rownames(res) <- NULL; write.csv(res, file.path(OUT, "spec14_results.csv"), row.names = FALSE); options(width = 200)
cat("=== mean z over the 8 variants by screen version, replicate mode, column and row set ===\n"); a <- aggregate(z_matched ~ screen + mode + cols + rows, res, mean); a$mode <- factor(a$mode, c("avg", "best", "ctl")); a <- a[order(a$cols, a$rows, a$screen, a$mode), ]; print(format(a, digits = 3), row.names = FALSE)
cat("\n=== best variant per screen version (mode = best replicate), with its selection-bias control ===\n")
bb <- res[res$mode == "best", ]; top <- do.call(rbind, lapply(split(bb, list(bb$screen, bb$cols, bb$rows), drop = TRUE), function(d) d[which.max(d$z_matched), ]))
top$z_control <- mapply(function(sv, fs, ty, ce, cs, rs) res$z_matched[res$screen == sv & res$features == fs & res$stat == ty & res$centred == ce & res$cols == cs & res$rows == rs & res$mode == "ctl"], top$screen, top$features, top$stat, top$centred, top$cols, top$rows)
top$z_avg <- mapply(function(sv, fs, ty, ce, cs, rs) res$z_matched[res$screen == sv & res$features == fs & res$stat == ty & res$centred == ce & res$cols == cs & res$rows == rs & res$mode == "avg"], top$screen, top$features, top$stat, top$centred, top$cols, top$rows)
top <- top[order(top$cols, top$rows, top$screen), ]; print(format(top[, c("screen", "cols", "rows", "features", "stat", "centred", "n_rows", "z_matched", "perm_p", "top1", "z_avg", "z_control")], digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- figure 18: z by mode for both screen versions
save_fig("fig18_G1_best_replicate", 13, 8.4, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(0, 0, 4, 0.5), family = "sans"); svs <- names(SV); mds <- c("avg", "best", "ctl"); mcol <- c("#2a78d6", "#eb6834", GREY); mlab <- c("replicates averaged", "best replicate", "control: best per cell")
  for (cs in c("all8", "powered")) { d <- res[res$cols == cs & res$rows == "ALL4", ]; par(mar = c(6.5, 4.6, 3.6, 0.8)); plot(NA, xlim = c(0.4, 6.6), ylim = range(c(d$z_matched, 0)) + c(-0.05, 0.1), xaxt = "n", xlab = "", ylab = "row-standardised matched score z (4 lines)", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8); abline(h = 0, col = INK2)
    for (i in seq_along(svs)) for (j in seq_along(mds)) { v <- d$z_matched[d$screen == svs[i] & d$mode == mds[j]]; at <- (i - 1) * 3 + j; boxplot(v, at = at, add = TRUE, axes = FALSE, boxwex = 0.75, outline = FALSE, col = adjustcolor(mcol[j], 0.85), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1); points(at + runif(length(v), -0.22, 0.22), v, pch = 16, cex = 0.45, col = adjustcolor(INK, 0.5)) }
    axis(1, at = c(2, 5), labels = sub(", ", "\n", svs), tick = FALSE, cex.axis = 0.75, line = 1); legend("topleft", mlab, fill = adjustcolor(mcol, 0.85), bty = "n", cex = 0.65)
    mtext(sprintf("%s", ifelse(cs == "all8", "all 8 shared screen KOs", sprintf("screen KOs with >= 50 cells (A: %d, D: %d)", length(intersect(ORDER, names(ss$n)[ss$n >= 50])), length(intersect(ORDER, names(n_g1)[n_g1 >= 50]))))), side = 3, adj = 0, line = 0.7, font = 2, cex = 0.8) }
  mtext("G1-only screen (cc regressed) vs cc-adjusted bulk: averaged replicates, best replicate and selection control", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1)
  mtext("each box = 8 variants (features x statistic x centring); z = matched screen KO's row-standardised score (0 = no better than other KOs); all 4 cell lines", outer = TRUE, side = 3, adj = 0, line = 1.1, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- figure 18b: matrix of the best G1-only best-replicate variant (8 shared TFs ordered + screen-only KOs)
g <- top[top$screen == "D G1 only, cc regressed" & top$cols == "all8" & top$rows == "ALL4", ]; R <- store[[paste(g$screen, g$features, g$stat, g$centred)]]
Call <- R$best; Z <- t(scale(t(Call[, ORDER]))); Z[!is.finite(Z)] <- 0; jm <- match(rows$ko, ORDER); perf <- sapply(ORDER, function(k) mean(Z[cbind(seq_len(nrow(Z)), jm)][rows$ko == k])); cols <- names(sort(perf, decreasing = TRUE))
lab <- rows; lab$pick <- R$pick; lab$ko <- factor(lab$ko, cols); o <- order(match(lab$line, LINES), as.integer(lab$ko)); lab <- lab[o, ]; lab$ko <- as.character(lab$ko); C <- Call[o, cols]; Cx <- Call[o, EXTRA, drop = FALSE]
draw_ordered_matrix("fig18b_G1_best_replicate_matrix", C, Cx, lab, n_g1, sprintf("%s (%s)", lab$ko, sub(".*_", "", lab$pick)), "CROP-seq G1 cells (cc regressed) vs cc-adjusted bulk, best replicate per KO",
  sprintf("%s genes, %s statistic, %s; 8 shared TFs ordered best to worst, then screen-only KOs MEF2A / MEF2D; cells per KO (G1) in brackets; replicate in brackets after each row", g$features, ifelse(g$stat == "t", "test", "log-difference"), ifelse(g$centred, "centred across KOs", "uncentred")),
  sprintf("all 4 lines, 8 TFs: best replicate z = %.2f (exact permutation p = %.3f) | replicates averaged z = %.2f | selection-bias control z = %.2f", g$z_matched, g$perm_p, g$z_avg, g$z_control))
