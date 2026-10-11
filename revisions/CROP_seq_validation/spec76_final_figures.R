# Spec 76: final compact figures of the bulk KO agreement analysis. Direction NTC - KO (positive = lower after KO); DIFF_DIR=KO_NTC for the reverse.
# PRIMARY configuration (chosen from the marginal effect of every axis in spec73, not from the single best cell of the grid): sample set S70, TMM logCPM, limma moderated t of
# NTC - KO, profiles centred over the 8 KOs of the line (removes the shared-NTC term, spec74 C), Spearman. Shown without and with the cell-cycle index as covariate.
# Feature spaces: signal500 (the 500 genes with the largest KO signal in the 2 lines NOT being compared) | screen genes (top-100 per KO in CROP-seq) | Cm20_genes (mean t over the
# top-50 genes of each of the 20 Cm) | Hallmark (mean t per Hallmark set).
# p-values: per KO label permutation (N_PERM, KO labels permuted within every line), BH over the 8 KOs; global = mean same-KO correlation; family-wise p over the feature spaces x
# cell-cycle options shown (max statistic over the same permutations). Robustness = fraction of ALL spec73 configurations (S70, centred/z, t/t_rest) in which the KO has perm p < 0.05.
# fig8A: KO-by-line correlation heatmap (32 x 32) + per-KO same-KO vs other-KO boxplots with p-values.   fig8B: KO x Cm dot heatmap, screen beside bulk, preserved pairs outlined.
# fig8C: the preserved KO / Cm pairs + TEAD4 / Cm20: module score NTC vs KO in the screen and in every bulk line.   fig8D: per-KO summary over feature spaces (permutation z).
# needs: spec72_build_profiles.R, spec73_agreement_benchmark.R, spec75_Cm_preserved.R, spec1_prepare_stats.R
# outputs: spec76_primary_per_KO.csv, spec76_primary_global.csv, spec76_KOline_correlation_<space>.csv, fig8A..fig8D {png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })
SET <- Sys.getenv("SET", unset = "S70"); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); perms <- make_perms(N_PERM)
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_%s%s.rds", SET, DSUF))); genes <- pr$genes; top50 <- cm_top50(); CM <- names(top50); ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
hall <- { suppressPackageStartupMessages(library(msigdbr)); h <- msigdbr(species = "Homo sapiens", category = "H"); split(h$gene_symbol, sub("HALLMARK_", "", h$gs_name)) }
scr <- intersect(unique(unlist(lapply(KOS, function(k) ss$genes[order(abs(ss$tcc[, k]), decreasing = TRUE)[1:100]]))), genes)
setmean <- function(A, sets) { sets <- lapply(sets, intersect, genes); sets <- sets[sapply(sets, length) >= 10]; out <- array(NA_real_, c(length(sets), dim(A)[2:3]), dimnames = c(list(names(sets)), dimnames(A)[2:3])); for (i in seq_along(sets)) out[i, , ] <- colMeans(A[sets[[i]], , , drop = FALSE]); out }
SPACES <- c(signal500 = "top-500 KO-signal genes", screen = "CROP-seq top genes", Cm20_genes = "20 Cm modules", Hallmark = "Hallmark sets"); COVS <- c(none = "no covariate", cc = "cell cycle regressed")
cl_for <- function(cv, sp) { A <- pr$P$tmm[[cv]]$t; switch(sp, signal500 = cl_from_array(A, "signal", pr$gene_signal), screen = cl_from_array(A, scr), Cm20_genes = cl_from_array(setmean(A, top50)), Hallmark = cl_from_array(setmean(A, hall))) }

## ---------------------------------------------------------------- primary results + family-wise p
null_of <- function(CL) { li <- lapply(PAIRS, function(p) match(p, LINES)); V <- lapply(seq_along(CL), function(i) matrix(CL[[i]][cbind(as.vector(perms[[li[[i]][1]]]), as.vector(perms[[li[[i]][2]]]))], N_PERM, length(KOS)))
  Reduce(`+`, lapply(V, function(v) { v[is.na(v)] <- 0; v })) / Reduce(`+`, lapply(V, function(v) !is.na(v))) }            # NA-aware (KO absent from a line, RUN3); identical when complete
RES <- list(); perko <- list(); glob <- list(); zobs <- c(); znull <- NULL
for (cv in names(COVS)) for (sp in names(SPACES)) { CL <- cl_for(cv, sp); a <- agreement(CL, perms); id <- paste(cv, sp, sep = "|"); RES[[id]] <- list(CL = CL, a = a)
  perko[[id]] <- cbind(set = SET, cov = cv, space = sp, a$per_KO); glob[[id]] <- data.frame(set = SET, cov = cv, space = sp, t(a$global))
  nu <- rowMeans(null_of(CL), na.rm = TRUE); zobs[id] <- (mean(a$intra, na.rm = TRUE) - mean(nu)) / sd(nu); znull <- cbind(znull, (nu - mean(nu)) / sd(nu)) }
perko <- do.call(rbind, perko); glob <- do.call(rbind, glob); glob$familywise_p <- sapply(zobs, function(z) (1 + sum(apply(znull, 1, max) >= z)) / (N_PERM + 1))
f73 <- file.path(OUT, sprintf("spec73_grid_per_KO_%s.csv", SET)); g73 <- read.csv(if (file.exists(f73)) f73 else file.path(OUT, "spec73_grid_per_KO.csv")); g73 <- g73[g73$set == SET & g73$prep != "none" & g73$stat %in% c("t", "t_rest"), ]
perko$robust_frac_same_space <- mapply(function(k, sp) mean(g73$perm_p[g73$KO == k & g73$space == sp] < 0.05), perko$KO, perko$space)
perko$robust_frac_all_spaces <- sapply(perko$KO, function(k) mean(g73$perm_p[g73$KO == k] < 0.05))
write.csv(perko, file.path(OUT, sprintf("spec76_primary_per_KO_%s%s.csv", SET, DSUF)), row.names = FALSE); write.csv(glob, file.path(OUT, sprintf("spec76_primary_global_%s%s.csv", SET, DSUF)), row.names = FALSE)
options(width = 250); cat(sprintf("\n=== PRIMARY configuration, %s, %s, %d permutations ===\n", SET, DLAB, N_PERM)); print(glob, digits = 3, row.names = FALSE)
print(perko[, c("cov", "space", "KO", "mean_intra", "median_intra", "median_cross", "n_pos_pairs", "wilcox_p", "perm_p", "perm_fdr", "perm_z", "robust_frac_same_space")], digits = 3, row.names = FALSE)

## ---------------------------------------------------------------- fig8A: 32 x 32 KO-line heatmap + boxplots, 2 feature spaces (genes, modules), no covariate and cc regressed
big_matrix <- function(CL) { nm <- as.vector(t(outer(KOS, sub("NCI-|-211H", "", LINES), paste))); M <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm)); sl <- sub("NCI-|-211H", "", LINES)
  for (i in seq_along(PAIRS)) { a <- sub("NCI-|-211H", "", PAIRS[[i]][1]); b <- sub("NCI-|-211H", "", PAIRS[[i]][2]); C <- CL[[i]]; for (k in KOS) for (k2 in KOS) { M[paste(k, a), paste(k2, b)] <- C[k, k2]; M[paste(k2, b), paste(k, a)] <- C[k, k2] } }; M }
fig8A <- function(cv) { sps <- c("signal500", "Cm20_genes"); hl <- list()
  for (sp in sps) { M <- big_matrix(RES[[paste(cv, sp, sep = "|")]]$CL); write.csv(round(M, 4), file.path(OUT, sprintf("spec76_KOline_correlation_%s_%s_%s%s.csv", sp, cv, SET, DSUF)))
    ko <- factor(sub(" .*", "", rownames(M)), KOS); ln <- sub(".* ", "", rownames(M)); lim <- max(0.3, quantile(abs(M), 0.98, na.rm = TRUE))
    hl[[sp]] <- Heatmap(M, name = paste0("rho_", sp), col = colorRamp2(c(-lim, 0, lim), c("#2a78d6", "#f0efec", "#eb6834")), na_col = "#d5d4cf", cluster_rows = FALSE, cluster_columns = FALSE, row_split = ko, column_split = ko,
      row_gap = unit(1.2, "mm"), column_gap = unit(1.2, "mm"), row_title_rot = 0, row_title_gp = gpar(fontsize = 8, fontface = "bold"), column_title_gp = gpar(fontsize = 8, fontface = "bold"), show_row_names = FALSE, show_column_names = FALSE, border = FALSE,
      left_annotation = rowAnnotation(line = ln, col = list(line = setNames(LINE_COL[LINES], sub("NCI-|-211H", "", LINES))), show_annotation_name = FALSE, simple_anno_size = unit(2.2, "mm"), show_legend = FALSE),
      heatmap_legend_param = list(title = "Spearman rho", direction = "horizontal", title_position = "lefttop", legend_width = unit(3, "cm"), labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 7)),
      layer_fun = function(j, i, x, y, w, h, fill) if (all(ko[i] == ko[j])) grid.rect(gp = gpar(col = INK, fill = NA, lwd = 1.2)), width = unit(7.4, "cm"), height = unit(7.4, "cm"),
      column_title = NULL, row_title_side = "left") }
  bx <- function(sp) { a <- RES[[paste(cv, sp, sep = "|")]]$a; k <- a$per_KO; o <- order(k$perm_p, -k$mean_intra); pushViewport(viewport(layout = grid.layout(1, 1))); yr <- range(c(unlist(a$cross), a$intra), na.rm = TRUE); yr <- yr + c(-0.05, 0.22) * diff(yr)
    pushViewport(plotViewport(c(3.2, 3.2, 1.6, 0.4), xscale = c(0.4, 8.6), yscale = yr)); grid.yaxis(gp = gpar(fontsize = 7)); grid.lines(unit(c(0.4, 8.6), "native"), unit(c(0, 0), "native"), gp = gpar(lty = 2, col = INK2))
    bxp1 <- function(v, x, col) { v <- v[!is.na(v)]; if (length(v) < 2) return(invisible()); qq <- quantile(v, c(0.25, 0.5, 0.75)); w <- range(v[v >= qq[1] - 1.5 * diff(qq[c(1, 3)]) & v <= qq[3] + 1.5 * diff(qq[c(1, 3)])]); grid.segments(unit(x, "native"), unit(w[1], "native"), unit(x, "native"), unit(w[2], "native"), gp = gpar(col = INK2))
      grid.rect(unit(x, "native"), unit(qq[1], "native"), unit(0.3, "native"), unit(qq[3] - qq[1], "native"), just = c("centre", "bottom"), gp = gpar(fill = col, col = INK2, lwd = 0.6)); grid.segments(unit(x - 0.15, "native"), unit(qq[2], "native"), unit(x + 0.15, "native"), unit(qq[2], "native"), gp = gpar(col = INK, lwd = 1.2)) }
    for (i in seq_along(o)) { kk <- k$KO[o[i]]; bxp1(a$cross[[o[i]]], i - 0.19, "#d5d4cf"); bxp1(a$intra[kk, ], i + 0.19, adjustcolor(KO_COL[kk], 0.55)); grid.points(unit(rep(i + 0.19, 6) + seq(-0.1, 0.1, length.out = 6), "native"), unit(a$intra[kk, ], "native"), pch = 16, size = unit(1.6, "mm"), gp = gpar(col = INK2))
      grid.text(sprintf("p=%s\n%s", formatC(k$perm_p[o[i]], format = "g", digits = 2), stars(k$perm_fdr[o[i]])), unit(i, "native"), unit(yr[2], "native"), just = "top", gp = gpar(fontsize = 6.2, fontface = ifelse(k$perm_fdr[o[i]] < 0.05, "bold", "plain")))
      grid.text(kk, unit(i, "native"), unit(-0.6, "lines"), rot = 45, just = "right", gp = gpar(fontsize = 7.5)) }
    grid.text("Spearman rho between two lines", x = unit(-2.5, "lines"), rot = 90, gp = gpar(fontsize = 7.5))
    grid.text(sprintf("%s: global p = %s (family-wise %s)", SPACES[[sp]], formatC(a$global[["perm_p"]], format = "g", digits = 2), formatC(glob$familywise_p[glob$cov == cv & glob$space == sp], format = "g", digits = 2)), x = unit(0, "npc"), y = unit(1, "npc") + unit(0.9, "lines"), just = "left", gp = gpar(fontsize = 8, fontface = "bold"))
    popViewport(2) }
  function() { grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(2.3, 10, 6.4), "cm"))))
    grid.text(sprintf("Agreement of bulk KO profiles across cell lines (%s, %d samples, %s)", DLAB, length(pr$samples), COVS[[cv]]), x = unit(4, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    grid.text("TMM logCPM, limma moderated t, centred over the KOs of each line. Heatmaps: every KO of every line vs every KO of the other lines\n(grey = same line, not compared); boxed = same KO. Boxplots: grey = the KO vs other KOs in another line, colour = same KO (dots = the 6 line pairs);\np = label permutation, stars = BH FDR over the 8 KOs (* <0.05, ** <0.01, *** <0.001).", x = unit(4, "mm"), y = unit(1, "npc") - unit(10, "mm"), just = c("left", "top"), gp = gpar(fontsize = 6.4, col = INK2, lineheight = 1.05), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    draw(Legend(labels = sub("NCI-|-211H", "", LINES), legend_gp = gpar(fill = LINE_COL[LINES]), nrow = 1, title = "line", title_position = "leftcenter", labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 7, fontface = "bold"), grid_height = unit(3, "mm"), grid_width = unit(3, "mm")), x = unit(1, "npc") - unit(4, "mm"), y = unit(1, "npc") - unit(13, "mm"), just = c("right", "top"))
    for (j in 1:2) { pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j)); draw(hl[[sps[j]]], newpage = FALSE, heatmap_legend_side = "bottom", column_title = SPACES[[sps[j]]], column_title_gp = gpar(fontsize = 9, fontface = "bold")); popViewport()
      pushViewport(viewport(layout.pos.row = 3, layout.pos.col = j)); bx(sps[j]); popViewport() } } }
for (cv in names(COVS)) save_fig(sprintf("fig8A_KO_agreement_across_lines_%s_%s%s", cv, SET, DSUF), 8.6, 7.6, fig8A(cv))

## ---------------------------------------------------------------- fig8D: per-KO summary over feature spaces and covariate options
save_fig(sprintf("fig8D_KO_agreement_summary_%s%s", SET, DSUF), 8.2, 3.9, function() {
  Z <- sapply(names(RES), function(id) RES[[id]]$a$per_KO$perm_z); P <- sapply(names(RES), function(id) RES[[id]]$a$per_KO$perm_p); Fd <- sapply(names(RES), function(id) RES[[id]]$a$per_KO$perm_fdr); rownames(Z) <- KOS
  colnames(Z) <- SPACES[sub(".*\\|", "", colnames(Z))]; rb <- sapply(KOS, function(k) mean(g73$perm_p[g73$KO == k] < 0.05)); o <- order(-rowMeans(Z))
  ht <- Heatmap(Z[o, ], name = "z", col = colorRamp2(c(-4, 0, 4), c("#2a78d6", "#f0efec", "#eb6834")), cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(COVS[sub("\\|.*", "", names(RES))], COVS), rect_gp = gpar(col = "white", lwd = 2),
    cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%s\n%s", formatC(P[o, ][i, j], format = "g", digits = 2), sub("ns", "", stars(Fd[o, ][i, j]))), x, y, gp = gpar(fontsize = 6.5, fontface = ifelse(Fd[o, ][i, j] < 0.05, "bold", "plain"))),
    right_annotation = rowAnnotation(`configs with p<0.05` = anno_barplot(rb[o], gp = gpar(fill = INK2, col = NA), border = FALSE, width = unit(2.2, "cm"), axis_param = list(gp = gpar(fontsize = 6))), annotation_name_gp = gpar(fontsize = 7), annotation_name_rot = 0),
    row_names_side = "left", row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 8), column_names_rot = 30, column_title_gp = gpar(fontsize = 9, fontface = "bold"), width = unit(9.6, "cm"), height = unit(5.2, "cm"),
    heatmap_legend_param = list(title = "permutation z", labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 7)))
  draw(ht, column_title = sprintf("Same KO across cell lines: permutation z and p per TF (%s; number = p, stars = FDR)", DLAB), column_title_gp = gpar(fontsize = 9.5, fontface = "bold"), padding = unit(c(2, 2, 2, 2), "mm")) })

## ---------------------------------------------------------------- fig8B / fig8C: Cm modules, screen vs bulk
tab <- read.csv(file.path(OUT, sprintf("spec75_Cm_KO_table_%s%s.csv", SET, DSUF))); pl <- read.csv(file.path(OUT, sprintf("spec75_Cm_KO_per_line_%s%s.csv", SET, DSUF)))
mat <- function(v) { M <- matrix(NA_real_, length(KOS), length(CM), dimnames = list(KOS, CM)); M[cbind(tab$KO, tab$Cm)] <- tab[[v]]; M }
Zs <- sign(mat("screen_effect")) * pmin(-log10(mat("screen_fdr")), 10); Zb <- mat("bulk_z"); Pr <- mat("preserved_specific") == 1; Fs <- mat("screen_fdr"); Fb <- mat("bulk_fdr")
save_fig(sprintf("fig8B_Cm_KO_screen_vs_bulk_%s%s", SET, DSUF), 11.6, 4.3, function() {
  dot <- function(M, F, lim, nm, ttl, leg) Heatmap(M, name = nm, col = colorRamp2(c(-lim, 0, lim), c("#2a78d6", "#f0efec", "#eb6834")), rect_gp = gpar(type = "none"), cluster_rows = FALSE, cluster_columns = FALSE, column_title = ttl, column_title_gp = gpar(fontsize = 9, fontface = "bold"),
    cell_fun = function(j, i, x, y, w, h, fill) { grid.rect(x, y, w, h, gp = gpar(col = GRID, fill = NA, lwd = 0.4)); r <- 0.2 + 0.8 * min(abs(M[i, j]) / lim, 1); grid.circle(x, y, r = unit(1.9 * r, "mm"), gp = gpar(fill = fill, col = if (F[i, j] < 0.1) INK2 else NA, lwd = 0.5))
      if (Pr[i, j]) grid.rect(x, y, w, h, gp = gpar(col = INK, fill = NA, lwd = 1.6)) }, row_names_side = "left", row_names_gp = gpar(fontsize = 8.5), column_names_gp = gpar(fontsize = 8), column_names_rot = 45, width = unit(10.6, "cm"), height = unit(4.6, "cm"),
    heatmap_legend_param = list(title = leg, labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 7), direction = "horizontal", title_position = "lefttop", legend_width = unit(2.6, "cm")))
  draw(dot(Zs, Fs, 10, "screen", "CROP-seq screen (NCI-H2052)", "signed -log10 FDR") + dot(Zb, Fb, 5, "bulk", sprintf("bulk RNA-seq, %d cell lines combined", length(LINES)), "combined z"), heatmap_legend_side = "bottom", ht_gap = unit(8, "mm"),
       column_title = sprintf("Effect of each KO on the 20 Cm modules (%s; orange = module lower after KO): CROP-seq vs bulk", DLAB), column_title_gp = gpar(fontsize = 10, fontface = "bold"), padding = unit(c(2, 2, 2, 2), "mm"))
  grid.text("dot size = strength; thin outline = FDR < 0.1 (bulk) / < 0.05 shown by colour saturation (screen); thick box = PRESERVED and KO-specific (screen FDR<0.05, bulk FDR<0.1, same sign in >= 3 of 4 lines, top-3 KO for the module in both, and significant vs the other KOs in both)", x = unit(3, "mm"), y = unit(2.5, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 6.2, col = INK2)) })

sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); q <- load_audit(); meta <- load_full_audited_bulk()$meta; S <- apply(sc$score, 2, function(y) resid(lm(y ~ sc$cc)) + mean(y)); B <- pr$cm_scores$tmm
sel <- tab[tab$preserved_specific == 1, c("KO", "Cm", "bulk_p")]; sel <- sel[order(sel$bulk_p), ]; sel <- rbind(sel[!duplicated(sel$KO), ], head(sel[duplicated(sel$KO), ], 1)); sel <- sel[order(sel$bulk_p), ]; if (!any(sel$KO == "TEAD4" & sel$Cm == "Cm20")) sel <- rbind(sel, data.frame(KO = "TEAD4", Cm = "Cm20", bulk_p = NA))
save_fig(sprintf("fig8C_preserved_Cm_KO_pairs_%s%s", SET, DSUF), 2.35 * nrow(sel) + 0.6, 6.4, function() {
  layout(matrix(seq_len(2 * nrow(sel)), nrow = 2)); par(oma = c(0, 1.5, 3.6, 0.3), family = "sans")
  for (i in seq_len(nrow(sel))) { k <- sel$KO[i]; m <- sel$Cm[i]; r <- tab[tab$KO == k & tab$Cm == m, ]; a <- S[sc$group == "NTC", m]; b <- S[sc$group == k, m]
    par(mar = c(2.6, 3.6, 3.8, 0.5)); plot(NA, xlim = c(0.4, 2.6), ylim = quantile(c(a, b), c(0.005, 0.995)), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = if (i == 1) "module score, CROP-seq cells" else "", cex.axis = 0.7, cex.lab = 0.8)
    violin_at(a, 1, KO_COL[["NTC"]]); violin_at(b, 2, KO_COL[[k]]); segments(c(0.8, 1.8), c(median(a), median(b)), c(1.2, 2.2), lwd = 2, col = INK); axis(1, at = 1:2, labels = c("NTC", k), tick = FALSE, cex.axis = 0.8, line = -0.6)
    mtext(sprintf("%s KO, %s", k, m), side = 3, adj = 0, line = 2.2, font = 2, cex = 0.8); mtext(if (r$preserved_specific == 1) "preserved, KO-specific" else "NOT preserved in bulk", side = 3, adj = 0, line = 1.3, cex = 0.58, font = 2, col = if (r$preserved_specific == 1) INK2 else "#c0392b"); mtext(sprintf("screen: %s = %+.3f, FDR %s", DLAB, r$screen_effect, formatC(r$screen_fdr, format = "g", digits = 2)), side = 3, adj = 0, line = 0.5, cex = 0.55, col = INK2)
    d <- pl[pl$KO == k & pl$Cm == m, ]; par(mar = c(4.4, 3.6, 3, 0.5)); ys <- unlist(lapply(LINES, function(cl) { s <- pr$samples[meta[pr$samples, "Cell_Line"] == cl & meta[pr$samples, "Gene_Control"] %in% c("NTC", k)]; B[m, s] - mean(B[m, s]) }))
    plot(NA, xlim = c(0.5, 4.5), ylim = range(ys) + c(-0.05, 0.25) * diff(range(ys)), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = if (i == 1) "module score, bulk (centred per line)" else "", cex.axis = 0.7, cex.lab = 0.8); abline(h = 0, col = GRID)
    for (j in seq_along(LINES)) { cl <- LINES[j]; s <- pr$samples[meta[pr$samples, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; mu <- mean(B[m, s[g %in% c("NTC", k)]]); yn <- B[m, s[g == "NTC"]] - mu; yk <- B[m, s[g == k]] - mu
      segments(j - 0.18, mean(yn), j + 0.18, mean(yk), col = INK2, lwd = 1); points(rep(j - 0.18, length(yn)), yn, pch = 21, bg = KO_COL[["NTC"]], col = "white", cex = 1.5); points(rep(j + 0.18, length(yk)), yk, pch = 21, bg = KO_COL[[k]], col = "white", cex = 1.5)
      pp <- d$p[d$line == cl]; if (!length(pp)) { text(j, par("usr")[4], "KO not\nin RUN3", adj = c(0.5, 1), cex = 0.55, col = INK2, xpd = NA) } else text(j, par("usr")[4], sprintf("%s\np=%s", ifelse(d$effect[d$line == cl] > 0, "down", "up"), formatC(pp, format = "g", digits = 2)), adj = c(0.5, 1), cex = 0.55, font = ifelse(pp < 0.05, 2, 1), xpd = NA)
      text(j, par("usr")[3] - 0.05 * diff(par("usr")[3:4]), sub("NCI-|-211H", "", cl), srt = 45, adj = 1, xpd = NA, cex = 0.75) }
    mtext(sprintf("bulk combined: z = %+.2f, p = %s, FDR %s", r$bulk_z, formatC(r$bulk_p, format = "g", digits = 2), formatC(r$bulk_fdr, format = "g", digits = 2)), side = 3, adj = 0, line = 1.4, cex = 0.55, col = INK2)
    mtext(sprintf("vs other KOs: z = %+.2f, p = %s | camera z = %+.2f", r$bulk_rest_z, formatC(r$bulk_rest_p, format = "g", digits = 2), r$bulk_camera_z), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2) }
  mtext(sprintf("Cm modules whose KO effect is preserved between CROP-seq and bulk RNA-seq (%s)", SET), outer = TRUE, side = 3, adj = 0, line = 2.2, font = 2, cex = 0.95)
  mtext("top: CROP-seq cells (cell cycle regressed; bar = median). bottom: bulk samples per line, grey = NTC, colour = KO; 'down'/'up' = module after KO; p = limma per line (score ~ group + cell cycle + % assigned)", outer = TRUE, side = 3, adj = 0, line = 1, cex = 0.58, col = INK2) })
