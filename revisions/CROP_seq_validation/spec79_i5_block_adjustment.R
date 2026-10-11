# Spec 79: does accounting for the i5 library-index group (the only batch proxy in the sample sheet, spec77) improve agreement across lines and with CROP-seq? S70 set.
# Design fact (checked here): within a line the i5 group is nested in KO (both replicates of a KO always share it), so i5 cannot enter the limma model as a covariate
# (NCI-H2052: aliased with group; other lines: identifiable only from the 1-2 samples split from their KO-mates). What CAN remove an additive i5 shift is a WITHIN-BLOCK
# contrast: each KO profile centred over the KOs of the same i5 block of that line (= KO vs its block-mates; NTC cancels too). A KO alone in its block cannot be centred -> NA.
# Compared: "centred" (primary: over all 8 KOs of the line) vs "i5-block-centred", on the spec72 S70 profiles (TMM logCPM, limma moderated t NTC - KO), no covariate and cell cycle.
#   A across lines : agreement() (Spearman, KO-label permutation) on 4 feature spaces, as spec76
#   B vs CROP-seq  : per line, Spearman over the 8 KO x 20 Cm effects between bulk (module t, centred the same way) and the screen (cell-level module score, NTC - KO median,
#                    centred over the SAME set of KOs as the bulk profile of that line, so both sides estimate the same contrast); plus the 4-line mean
# needs: spec72_build_profiles.R (spec72_profiles_S70.rds), spec75_Cm_preserved.R (spec75_screen_scores_by_cell.rds), spec1_prepare_stats.R
# outputs: spec79_i5_blocks.csv, spec79_across_lines.csv, spec79_across_lines_per_KO.csv, spec79_vs_cropseq.csv, fig8F_i5_block_adjustment.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); perms <- make_perms(N_PERM)
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); genes <- pr$genes; top50 <- cm_top50(); CM <- names(top50); ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
meta <- load_full_audited_bulk()$meta; sheet <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv")); i5 <- setNames(sheet$i5_index_num, sheet$fastq_sample_id)

## ---------------------------------------------------------------- i5 block of every KO in every line (majority i5 of its S70 samples; ties -> the larger block)
blk <- do.call(rbind, lapply(LINES, function(cl) { s <- pr$samples[meta[pr$samples, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; b <- i5[s]
  do.call(rbind, lapply(c("NTC", KOS), function(k) { t <- table(b[g == k]); data.frame(line = cl, group = k, i5 = paste(names(t), collapse = "/"), n = sum(t), block = names(t)[which.max(t)]) })) }))
for (cl in LINES) { d <- blk$line == cl & blk$group != "NTC"; for (k in blk$group[d & grepl("/", blk$i5)]) { cand <- strsplit(blk$i5[d & blk$group == k], "/")[[1]]; sz <- sapply(cand, function(x) sum(blk$block[d] == x)); blk$block[d & blk$group == k] <- cand[which.max(sz)] } }
blk$block_size <- ave(seq_len(nrow(blk)), blk$line, blk$block, blk$group != "NTC", FUN = length); blk$block_size[blk$group == "NTC"] <- NA
write.csv(blk, file.path(OUT, "spec79_i5_blocks.csv"), row.names = FALSE); options(width = 200)
cat("i5 blocks (KOs of the same line and block are compared with each other):\n"); for (cl in LINES) { d <- blk[blk$line == cl, ]; cat(sprintf("  %-10s NTC in i5 %s | %s\n", cl, d$i5[d$group == "NTC"], paste(tapply(d$group[d$group != "NTC"], d$block[d$group != "NTC"], paste, collapse = ","), paste0("[i5 ", names(tapply(d$group[d$group != "NTC"], d$block[d$group != "NTC"], length)), "]"), collapse = "  "))) }
BLOCK <- lapply(setNames(LINES, LINES), function(cl) { d <- blk[blk$line == cl & blk$group != "NTC", ]; setNames(d$block, d$group)[KOS] })

centre_arr <- function(A, mode) { if (mode == "uncentred") return(A); for (cl in LINES) { X <- A[, , cl, drop = FALSE][, , 1]
    if (mode == "centred") X <- X - rowMeans(X, na.rm = TRUE) else for (b in unique(BLOCK[[cl]])) { k <- KOS[BLOCK[[cl]] == b]; X[, k] <- if (length(k) < 2) NA else X[, k, drop = FALSE] - rowMeans(X[, k, drop = FALSE]) }
    A[, , cl] <- X }; A }
setmean <- function(A, sets) { sets <- lapply(sets, intersect, dimnames(A)[[1]]); sets <- sets[sapply(sets, length) >= 10]; out <- array(NA_real_, c(length(sets), dim(A)[2:3]), dimnames = c(list(names(sets)), dimnames(A)[2:3])); for (i in seq_along(sets)) out[i, , ] <- colMeans(A[sets[[i]], , , drop = FALSE]); out }
hall <- { suppressPackageStartupMessages(library(msigdbr)); h <- msigdbr(species = "Homo sapiens", category = "H"); split(h$gene_symbol, sub("HALLMARK_", "", h$gs_name)) }
scr <- intersect(unique(unlist(lapply(KOS, function(k) ss$genes[order(abs(ss$tcc[, k]), decreasing = TRUE)[1:100]]))), genes)
SPACES <- c(signal500 = "KO-signal genes", screen = "CROP-seq genes", Cm20_genes = "Cm modules", Hallmark = "Hallmark"); MODES <- c(centred = "centred over all KOs (primary)", i5_block = "centred within i5 block")

## ---------------------------------------------------------------- A across lines
A_res <- list(); K_res <- list()
for (cv in c("none", "cc")) for (md in names(MODES)) { A <- centre_arr(pr$P$tmm[[cv]]$t, md)
  for (sp in names(SPACES)) { CL <- switch(sp, signal500 = cl_from_array(A, "signal", pr$gene_signal, centre = FALSE), screen = cl_from_array(A, scr, centre = FALSE),
                                          Cm20_genes = cl_from_array(setmean(A, top50), centre = FALSE), Hallmark = cl_from_array(setmean(A, hall), centre = FALSE))
    a <- agreement(CL, perms); A_res[[length(A_res) + 1]] <- data.frame(cov = cv, mode = md, space = sp, t(a$global))
    K_res[[length(K_res) + 1]] <- cbind(cov = cv, mode = md, space = sp, a$per_KO) } }
A_res <- do.call(rbind, A_res); K_res <- do.call(rbind, K_res)
write.csv(A_res, file.path(OUT, "spec79_across_lines.csv"), row.names = FALSE); write.csv(K_res, file.path(OUT, "spec79_across_lines_per_KO.csv"), row.names = FALSE)
cat("\n=== A across lines: global agreement (perm p, z, mean same-KO rho) ===\n"); w <- reshape(A_res[, c("cov", "space", "mode", "perm_p", "perm_z", "diag_mean")], idvar = c("cov", "space"), timevar = "mode", direction = "wide"); print(w, digits = 3, row.names = FALSE)
cat("\nper KO, KO-signal genes (perm p): \n"); for (cv in c("none", "cc")) { x <- K_res[K_res$cov == cv & K_res$space == "signal500", ]; print(cbind(cov = cv, reshape(x[, c("KO", "mode", "perm_p", "n_pairs")], idvar = "KO", timevar = "mode", direction = "wide")), digits = 3, row.names = FALSE) }

## ---------------------------------------------------------------- B vs CROP-seq (module level, per line; screen centred over the same KOs as the bulk profile)
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds"))
screen_eff <- function(cc) { S <- sc$score; if (cc) S <- apply(S, 2, function(y) resid(lm(y ~ sc$cc)) + mean(y)); n <- apply(S[sc$group == "NTC", ], 2, median)
  sapply(KOS, function(k) SGN * (n - apply(S[sc$group == k, , drop = FALSE], 2, median))) }                     # Cm x KO, NTC - KO
B_res <- do.call(rbind, lapply(c("none", "cc"), function(cv) { E <- screen_eff(cv == "cc"); E <- E[CM, ]
  do.call(rbind, lapply(c("uncentred", names(MODES)), function(md) { Mb <- centre_arr(pr$M$tmm[[cv]]$t[CM, , , drop = FALSE], md)
    do.call(rbind, lapply(LINES, function(cl) { Es <- E; if (md == "uncentred") { } else if (md == "centred") Es <- Es - rowMeans(Es) else for (b in unique(BLOCK[[cl]])) { k <- KOS[BLOCK[[cl]] == b]; Es[, k] <- if (length(k) < 2) NA else Es[, k, drop = FALSE] - rowMeans(Es[, k, drop = FALSE]) }
      x <- as.vector(Es); y <- as.vector(Mb[, , cl]); ok <- is.finite(x) & is.finite(y); ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      data.frame(cov = cv, mode = md, line = cl, n_pairs = sum(ok), rho = unname(ct$estimate), p = ct$p.value, sign_agree = mean(sign(x[ok]) == sign(y[ok]))) })) })) }))
write.csv(B_res, file.path(OUT, "spec79_vs_cropseq.csv"), row.names = FALSE)
cat("\n=== B bulk vs CROP-seq, 8 KO x 20 Cm per line (Spearman) ===\n"); print(reshape(B_res[, c("cov", "line", "mode", "rho", "p", "n_pairs")], idvar = c("cov", "line"), timevar = "mode", direction = "wide"), digits = 3, row.names = FALSE)
cat("\nmean rho over the 4 lines:\n"); print(round(tapply(B_res$rho, list(B_res$cov, B_res$mode), mean), 3))

## ---------------------------------------------------------------- figure
MCOL <- c(uncentred = "grey65", centred = "#2a78d6", i5_block = "#1baf7a")
save_fig("fig8F_i5_block_adjustment", 14, 5.4, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.3, 1.1, 1)); par(oma = c(0, 0, 3.4, 0.4), mar = c(6.6, 4.2, 2.6, 0.6), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
  ids <- expand.grid(space = names(SPACES), cov = c("none", "cc"), stringsAsFactors = FALSE)
  Z <- sapply(names(MODES), function(md) mapply(function(sp, cv) A_res$perm_z[A_res$space == sp & A_res$cov == cv & A_res$mode == md], ids$space, ids$cov))
  P <- sapply(names(MODES), function(md) mapply(function(sp, cv) A_res$perm_p[A_res$space == sp & A_res$cov == cv & A_res$mode == md], ids$space, ids$cov))
  bp <- barplot(t(Z), beside = TRUE, col = MCOL[names(MODES)], border = NA, las = 1, ylab = "global agreement across lines (permutation z)", names.arg = rep("", nrow(ids)), ylim = range(0, Z) * c(1.1, 1.2), cex.axis = 0.75); abline(h = qnorm(0.95), lty = 3, col = MUTED)
  text(bp, t(Z), ifelse(t(P) < 0.001, "***", ifelse(t(P) < 0.01, "**", ifelse(t(P) < 0.05, "*", ""))), pos = 3, cex = 0.8, xpd = NA)
  text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), paste0(SPACES[ids$space], "\n", ifelse(ids$cov == "none", "no covariate", "cell cycle regr.")), srt = 45, adj = 1, xpd = NA, cex = 0.62)
  legend("topright", MODES, fill = MCOL[names(MODES)], border = NA, bty = "n", cex = 0.7); ttl("A  agreement across cell lines")
  pk <- lapply(names(MODES), function(md) { d <- K_res[K_res$cov == "none" & K_res$space == "signal500" & K_res$mode == md, ]; setNames(d$perm_z, d$KO)[KOS] }); names(pk) <- names(MODES); o <- order(-pk$centred)
  par(mar = c(4.6, 4.2, 2.6, 0.6)); plot(NA, xlim = c(0.5, 8.5), ylim = range(c(unlist(pk), qnorm(0.95)), na.rm = TRUE) + c(-0.3, 0.3), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "same-KO agreement (permutation z)", cex.axis = 0.75)
  abline(h = c(0, qnorm(0.95)), lty = c(2, 3), col = c(INK, MUTED)); segments(1:8, pk$centred[o], 1:8, pk$i5_block[o], col = "grey60", lwd = 1.5)
  for (md in names(MODES)) points(1:8, pk[[md]][o], pch = 16, cex = 1.6, col = MCOL[[md]])
  text(1:8, par("usr")[3] - 0.04 * diff(par("usr")[3:4]), KOS[o], srt = 45, adj = 1, xpd = NA, cex = 0.75, font = 2, col = KO_COL[KOS[o]]); ttl("B  per TF (KO-signal genes, no covariate)")
  x <- B_res; par(mar = c(4.6, 4.2, 2.6, 0.6)); lv <- expand.grid(line = LINES, cov = c("none", "cc"), stringsAsFactors = FALSE)
  R <- sapply(names(MCOL), function(md) mapply(function(cl, cv) x$rho[x$line == cl & x$cov == cv & x$mode == md], lv$line, lv$cov))
  bp <- barplot(t(R), beside = TRUE, col = MCOL, border = NA, las = 1, ylab = "Spearman rho, bulk vs CROP-seq (8 KO x 20 Cm)", names.arg = rep("", nrow(lv)), cex.axis = 0.75); abline(h = 0, col = INK)
  text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), paste0(sub("NCI-|-211H", "", lv$line), "\n", ifelse(lv$cov == "none", "no cov.", "cc regr.")), srt = 45, adj = 1, xpd = NA, cex = 0.62)
  Pp <- sapply(names(MCOL), function(md) mapply(function(cl, cv) x$p[x$line == cl & x$cov == cv & x$mode == md], lv$line, lv$cov)); text(bp, pmax(t(R), 0), ifelse(t(Pp) < 0.05, "*", ""), pos = 3, cex = 0.9, xpd = NA)
  legend("topright", c("NTC - KO as is (not centred)", MODES), fill = MCOL, border = NA, bty = "n", cex = 0.62); ttl("C  matching with CROP-seq (module effects, per line)")
  mtext("Accounting for the i5 library-index batch (S70): KO profiles centred within their i5 block vs over all KOs", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0)
  mtext("i5 is nested in KO within every line, so it cannot be a model covariate; the within-block contrast cancels an additive i5 shift. TMM logCPM, limma t (NTC - KO); Spearman; permutation p; * p<0.05", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})
