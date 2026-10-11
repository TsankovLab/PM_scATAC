# Spec 16: CROP-seq KO profiles built ONLY from cells in which the target gene is undetected (lower expression), compared with cc-adjusted bulk.
# Per KO k: cells of KO k with target k undetected (log-normalised expression == 0) vs NTC cells; NTC reference = all NTC cells ("NTC all") or only NTC cells that also have gene k undetected
# ("NTC undetected", removes the detection/dropout confound). cc regressed pooled over the cells used. The 8 shared TFs + MEF2A/MEF2D genes are excluded from the features so the
# comparison is not driven by the target genes themselves. KOs need >= 15 undetected cells (HMGA1 is dropped: expressed in ~88% of cells). Baseline: the all-cell, cc-regressed screen on the SAME genes and KOs.
# Bulk: cc-adjusted replicate-level statistics (t | cc-regressed log difference); modes: replicates averaged, best replicate, control (best replicate per cell).
# needs: spec1_prepare_stats.R
# outputs: spec16_results.csv, fig20_undetected_cells.{png,pdf}, fig20b_undetected_matrix.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
ctx <- load_spec_context(LINES); ss <- ctx$ss; bs <- ctx$bs; ALLC <- c(ORDER, EXTRA)
srt <- load_screen(); md <- srt@meta.data; Yall <- GetAssayData(srt, layer = "data"); ccv <- md$cc
gx <- setdiff(Reduce(intersect, c(list(ss$genes), lapply(bs[LINES], `[[`, "genes"))), ALLC); Yx <- Yall[gx, ]; var_x <- intersect(ctx$F_var, gx)
ntc_i <- which(md$merged_call == "NTC")

## ---------------------------------------------------------------- screen profiles from undetected cells
MINC <- 15; prof <- list(`undetected, NTC all` = list(t = list(), diff = list()), `undetected, NTC undetected` = list(t = list(), diff = list())); cnt <- data.frame()
for (tf in ALLC) { ki <- which(md$merged_call == tf & Yall[tf, ] == 0); nu <- ntc_i[Yall[tf, ntc_i] == 0]
  cnt <- rbind(cnt, data.frame(KO = tf, cells_total = sum(md$merged_call == tf), cells_undetected = length(ki), NTC_undetected = length(nu))); if (length(ki) < MINC) next
  for (nv in names(prof)) { ref <- if (nv == "undetected, NTC all") ntc_i else nu; ci <- c(ki, ref); St <- screen_ko_stats(Yx[, ci], c(rep("KO", length(ki)), rep("NTC", length(ref))), ccv[ci], TRUE)
    prof[[nv]]$t[[tf]] <- St$t[, "KO"]; prof[[nv]]$diff[[tf]] <- St$diff[, "KO"] } }
write.csv(cnt, file.path(OUT, "spec16_cell_counts.csv"), row.names = FALSE); options(width = 200); print(format(cnt, digits = 3), row.names = FALSE)
avail <- names(prof[[1]]$t); avail_order <- intersect(ORDER, avail); cat("KOs used:", paste(avail, collapse = ", "), "\n")
SA <- screen_ko_stats(Yx, md$merged_call, ccv, TRUE)                                                                # baseline: all cells, cc regressed, same genes
SV <- list(`all cells, cc regressed` = list(t = SA$t[, avail], diff = SA$diff[, avail]))
for (nv in names(prof)) SV[[nv]] <- list(t = do.call(cbind, prof[[nv]]$t), diff = do.call(cbind, prof[[nv]]$diff))
n_und <- setNames(cnt$cells_undetected, cnt$KO); powered <- intersect(avail_order, names(n_und)[n_und >= 50]); cat("columns with >= 50 undetected cells:", paste(powered, collapse = ", "), "\n")

## ---------------------------------------------------------------- bulk replicate profiles (cc adjusted, genes gx)
REP <- list()
for (cl in LINES) { m <- ctx$meta[ctx$meta$Cell_Line == cl, ]; s <- rownames(m); y <- ctx$lg[gx, s, drop = FALSE]; cc <- ctx$ccb$cc[match(s, ctx$ccb$sample)]
  g <- factor(m$Gene_Control, levels = c("NTC", sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix"))); X <- model.matrix(~ cc + g); q <- qr(X)
  s2 <- pmax(colSums(qr.resid(q, t(y))^2) / (length(s) - q$rank), 1e-4); yadj <- y - outer(qr.coef(q, t(y))["cc", ], cc - mean(cc))
  ntc_s <- s[m$Gene_Control == "NTC"]; ko_s <- s[m$Gene_Control != "NTC"]; D <- yadj[, ko_s, drop = FALSE] - rowMeans(yadj[, ntc_s, drop = FALSE])
  Tm <- D / sqrt(s2 * (1 + 1 / length(ntc_s))); Tm[!is.finite(Tm)] <- 0; D[!is.finite(D)] <- 0; REP[[cl]] <- list(t = Tm, diff = D, ko = setNames(as.character(m[ko_s, "Gene_Control"]), ko_s)) }
rows <- do.call(rbind, lapply(LINES, function(cl) { k <- REP[[cl]]$ko; k <- k[k %in% avail_order]; data.frame(line = cl, ko = unique(k)) })); rows$key <- paste(rows$line, rows$ko)
run_variant <- function(sv, feats, ty, cen) {
  S <- prep_dataset(SV[[sv]][[ty]][feats, avail, drop = FALSE], cen, ty == "diff")[, avail, drop = FALSE]
  avg <- matrix(NA, nrow(rows), length(avail)); best <- ctl <- avg; pick <- rep(NA_character_, nrow(rows))
  for (cl in LINES) { M <- prep_dataset(REP[[cl]][[ty]][feats, , drop = FALSE], cen, ty == "diff"); ko <- REP[[cl]]$ko
    for (i in which(rows$line == cl)) { ss_k <- names(ko)[ko == rows$ko[i]]; cr <- cor(M[, ss_k, drop = FALSE], S); mcol <- match(rows$ko[i], avail)
      avg[i, ] <- cor(rowMeans(M[, ss_k, drop = FALSE]), S)[1, ]; b <- which.max(cr[, mcol]); best[i, ] <- cr[b, ]; ctl[i, ] <- apply(cr, 2, max); pick[i] <- ss_k[b] } }
  colnames(avg) <- colnames(best) <- colnames(ctl) <- avail; list(avg = avg, best = best, ctl = ctl, pick = pick) }
row_sets <- list(ALL3 = rows$line %in% THREE, ALL4 = rep(TRUE, nrow(rows)), H2052 = rows$line == "NCI-H2052")
res <- list(); store <- list()
for (sv in names(SV)) for (fs in c("all_genes", "var_genes")) for (ty in c("t", "diff")) for (cen in c(FALSE, TRUE)) {
  feats <- if (fs == "all_genes") gx else var_x; R <- run_variant(sv, feats, ty, cen); store[[paste(sv, fs, ty, cen)]] <- R
  for (cs in c("all", "powered")) { cols <- if (cs == "all") avail_order else powered
    for (rs in names(row_sets)) { sel <- row_sets[[rs]] & rows$ko %in% cols; if (sum(sel) < 3) next
      for (md_ in c("avg", "best", "ctl")) res[[length(res) + 1]] <- data.frame(screen = sv, features = fs, stat = ty, centred = cen, cols = cs, n_cols = length(cols), rows = rs, mode = md_, evaluate_specificity(R[[md_]][sel, , drop = FALSE], rows$ko[sel], cols = cols)) } } }
res <- do.call(rbind, res); rownames(res) <- NULL; write.csv(res, file.path(OUT, "spec16_results.csv"), row.names = FALSE)
cat("\n=== mean z over the 8 variants (features x statistic x centring) ===\n"); a <- aggregate(z_matched ~ screen + mode + cols + rows, res, mean); a$mode <- factor(a$mode, c("avg", "best", "ctl")); a <- a[order(a$cols, a$rows, a$screen, a$mode), ]
print(format(a[a$rows %in% c("ALL3", "ALL4"), ], digits = 3), row.names = FALSE)
cat("\n=== best variant per screen version: averaged replicates ===\n"); bb <- res[res$mode == "avg", ]; top <- do.call(rbind, lapply(split(bb, list(bb$screen, bb$cols, bb$rows), drop = TRUE), function(d) d[which.max(d$z_matched), ])); top <- top[order(top$cols, top$rows, top$screen), ]
print(format(top[top$rows %in% c("ALL3", "ALL4"), c("screen", "cols", "rows", "features", "stat", "centred", "n_rows", "matched_r", "other_r", "z_matched", "perm_p", "top1")], digits = 3), row.names = FALSE)
cat("\n=== best variant per screen version: best replicate (with its control) ===\n"); bb <- res[res$mode == "best", ]; top2 <- do.call(rbind, lapply(split(bb, list(bb$screen, bb$cols, bb$rows), drop = TRUE), function(d) d[which.max(d$z_matched), ]))
top2$z_control <- mapply(function(sv, fs, ty, ce, cs, rs) res$z_matched[res$screen == sv & res$features == fs & res$stat == ty & res$centred == ce & res$cols == cs & res$rows == rs & res$mode == "ctl"], top2$screen, top2$features, top2$stat, top2$centred, top2$cols, top2$rows); top2 <- top2[order(top2$cols, top2$rows, top2$screen), ]
print(format(top2[top2$rows %in% c("ALL3", "ALL4"), c("screen", "cols", "rows", "features", "stat", "centred", "n_rows", "z_matched", "perm_p", "top1", "z_control")], digits = 3), row.names = FALSE)
cat("\ntests:", nrow(res), "| p<0.05:", sum(res$perm_p < 0.05), "(expected ~", round(0.05 * nrow(res), 1), ")\n")

## ---------------------------------------------------------------- figure 20
save_fig("fig20_undetected_cells", 15, 8.6, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(0, 0, 4, 0.5), family = "sans"); svs <- names(SV); mds <- c("avg", "best", "ctl"); mcol <- c("#2a78d6", "#eb6834", GREY); mlab <- c("replicates averaged", "best replicate", "control: best per cell")
  for (cs in c("all", "powered")) { d <- res[res$cols == cs & res$rows == "ALL4", ]; par(mar = c(7, 4.6, 3.6, 0.8)); plot(NA, xlim = c(0.4, 3 * length(svs) + 0.6), ylim = range(c(d$z_matched, 0)) + c(-0.05, 0.12), xaxt = "n", xlab = "", ylab = "row-standardised matched score z (4 lines)", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8); abline(h = 0, col = INK2)
    for (i in seq_along(svs)) for (j in seq_along(mds)) { v <- d$z_matched[d$screen == svs[i] & d$mode == mds[j]]; at <- (i - 1) * 3 + j; boxplot(v, at = at, add = TRUE, axes = FALSE, boxwex = 0.75, outline = FALSE, col = adjustcolor(mcol[j], 0.85), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1); points(at + runif(length(v), -0.22, 0.22), v, pch = 16, cex = 0.45, col = adjustcolor(INK, 0.5)) }
    axis(1, at = seq(2, 3 * length(svs), 3), labels = gsub(", ", "\n", svs), tick = FALSE, cex.axis = 0.7, line = 1); legend("topleft", mlab, fill = adjustcolor(mcol, 0.85), bty = "n", cex = 0.65)
    mtext(sprintf("%s (%d KOs)", ifelse(cs == "all", "all KOs with >= 15 undetected cells", "KOs with >= 50 undetected cells"), length(if (cs == "all") avail_order else powered)), side = 3, adj = 0, line = 0.7, font = 2, cex = 0.8) }
  mtext("CROP-seq profiles from cells with the target gene undetected vs cc-adjusted bulk", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1)
  mtext(sprintf("each box = 8 variants (features x statistic x centring); target genes excluded from the features; KOs: %s; z = matched KO's row-standardised score (0 = no better than other KOs)", paste(avail_order, collapse = ", ")), outer = TRUE, side = 3, adj = 0, line = 1.1, cex = 0.58, col = INK2)
})
## matrix of the best undetected-cell variant (averaged replicates; shared TFs ordered best -> worst, then the screen-only KOs MEF2A / MEF2D)
g <- top[top$screen == "undetected, NTC all" & top$cols == "all" & top$rows == "ALL4", ]; R <- store[[paste(g$screen, g$features, g$stat, g$centred)]]; Call <- R$avg
Z <- t(scale(t(Call[, avail_order]))); Z[!is.finite(Z)] <- 0; jm <- match(rows$ko, avail_order); perf <- sapply(avail_order, function(k) mean(Z[cbind(seq_len(nrow(Z)), jm)][rows$ko == k])); cols <- names(sort(perf, decreasing = TRUE))
lab <- rows; lab$ko <- factor(lab$ko, cols); o <- order(match(lab$line, LINES), as.integer(lab$ko)); lab <- lab[o, ]; lab$ko <- as.character(lab$ko); C <- Call[o, cols]; Cx <- Call[o, intersect(EXTRA, avail), drop = FALSE]
draw_ordered_matrix("fig20b_undetected_matrix", C, Cx, lab, n_und, lab$ko, "CROP-seq cells with the target undetected vs cc-adjusted bulk (replicates averaged)",
  sprintf("%s genes, %s statistic, %s; shared TFs ordered best to worst, then screen-only KOs MEF2A / MEF2D; undetected cells per KO in brackets; target genes excluded", g$features, ifelse(g$stat == "t", "test", "log-difference"), ifelse(g$centred, "centred", "uncentred")),
  sprintf("all 4 lines, %d shared TFs: z = %.2f (exact permutation p = %.3f); best variant of this screen version, so optimistic", length(cols), g$z_matched, g$perm_p))
