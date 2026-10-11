# Spec 115: every processing step between counts and the KO x Cm effect, varied one factor at a time AND as a full grid, scored on three questions at once.
# Purpose: find which step (if any) is throwing the reproducibility away, rather than tuning one pipeline until something is significant.
#   bulk factors
#     set     S70 (70 samples) | RUN3 (S70 minus the samples flagged in the batch / QC review, agree_common.R)
#     norm    cp10k | tmm (per line) | tmmq (tmm + quantile normalisation across samples) | rle (DESeq-style median-of-ratios size factors) | uq (upper quartile)
#             (DESeq2 vst is not installed in the meso_scatac environment; rle is the same size-factor estimator on the log scale)
#     score   ams (Seurat AddModuleScore, top-50)  | mz30 | mz50 | mz100 (mean of the within-line gene z-scores over the top-N module genes: no control-gene subtraction,
#             every gene weighted equally, immune to a few highly expressed genes dominating)
#     cov     none | cc | pa (% reads assigned to genes, the main technical axis of spec71) | cc_pa      -- regressed out of the module scores with the KO groups protected
#     zline   module scores z-scored within the line (yes | no)
#     ref     ntc (KO - NTC) | rest (KO - mean of the other 8 groups)
#     stat    diff (difference of replicate-averaged group means) | t (limma moderated t of the same contrast)
#   screen  : AddModuleScore per cell (spec75), z over cells; cc regressed out when the bulk cov contains cc; same reference as the bulk side.
#   metrics (sign convention KO - reference: NEGATIVE = module lower in the KO)
#     M1 screen_rho   Spearman over the 160 KO x Cm effects, screen vs each line; mean of the 4 lines, H2052 alone, and screen vs the 4-line mean
#        screen_rho_rel   the same on the 6 KOs whose screen profile is reliable (spec114: HMGA1 13 cells, SOX9 one guide / 59 cells dropped)
#     M2 lines_dmo / lines_z   agreement() across the 4 bulk lines on the 20-module profiles: same-KO minus other-KO Spearman, KO-label permutation z; n KOs with p < 0.05
#     M3 down_sign    for the pairs the screen puts clearly DOWN in the KO (lowest quartile of the screen KO - rest effect, fixed set): fraction of line-level bulk effects
#                     that are negative (chance 0.5), and how many of those pairs are negative in all 4 lines (chance = n / 16)
# Selection is kept honest by construction: M2 never sees the screen, so a configuration picked on M2 can be validated on M1 / M3 (and the reverse).
# needs: spec72 (S70 sample list), spec71 audit, spec75 (screen scores)
# outputs: spec115_grid.csv, spec115_marginal_effects.csv, spec115_one_factor_at_a_time.csv, spec115_profiles.rds, fig27_stepwise_pipeline.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 230)
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "2000")); permsL <- make_perms(N_PERM); GR <- c("NTC", KOS); CM <- paste0("Cm", 1:20); REL_KO <- setdiff(KOS, c("HMGA1", "SOX9"))
fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit(); pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds"))
SETS <- list(S70 = pr$samples, RUN3 = setdiff(pr$samples, RUN3_EXCLUDE)); mods <- readRDS(CM_RDS); mods <- mods[CM]
cat("module gene list lengths:", range(lengths(mods)), "\n")
sp <- function(a, b) suppressWarnings(cor(as.vector(a), as.vector(b), method = "spearman", use = "pairwise.complete.obs"))

## ---------------------------------------------------------------- screen effects (fixed scores; cc on / off; ntc / rest)
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); keep <- sc$group %in% GR; gS <- sc$group[keep]; ccS <- sc$cc[keep]
scr_eff <- function(cc, ref) { z <- sc$score[keep, ]; if (cc) z <- apply(z, 2, function(v) resid(lm(v ~ ccS))); z <- scale(z); colnames(z) <- CM
  M <- sapply(GR, function(g) colMeans(z[gS == g, , drop = FALSE]))
  if (ref == "ntc") M[, KOS] - M[, "NTC"] else sapply(KOS, function(k) M[, k] - rowMeans(M[, setdiff(GR, k)])) }
SCR <- list(); for (cc in c(FALSE, TRUE)) for (ref in c("ntc", "rest")) SCR[[paste(cc, ref)]] <- scr_eff(cc, ref)
base_s <- SCR[["FALSE rest"]]; DOWN <- which(base_s <= quantile(base_s, 0.25), arr.ind = TRUE)          # fixed set of clearly-down pairs
cat("screen-down pairs (lowest quartile of KO - rest):", nrow(DOWN), "; per KO:", paste(sprintf("%s %d", KOS, tabulate(DOWN[, 2], 8)), collapse = " "), "\n")

## ---------------------------------------------------------------- bulk: normalised matrices and module scores per set x norm x scoring
norm_mat <- function(ids, nm) { ln <- meta[ids, "Cell_Line"]; det <- rownames(cs)[rowSums(cs[, ids] > 0) > 0]
  per_line <- function(f) do.call(cbind, lapply(LINES, function(cl) f(cs[det, ids[ln == cl]])))[, ids]
  switch(nm, cp10k = log2(sweep(cs[det, ids], 2, colSums(cs[det, ids]), "/") * 1e4 + 1), tmm = per_line(function(x) normalise(x, "tmm")$y),
         tmmq = { y <- per_line(function(x) normalise(x, "tmm")$y); z <- limma::normalizeQuantiles(y); dimnames(z) <- dimnames(y); z },
         rle = per_line(function(x) edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(x), method = "RLE"), log = TRUE, prior.count = 2)),
         uq = per_line(function(x) edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(x), method = "upperquartile"), log = TRUE, prior.count = 2)) ) }
score_mat <- function(y, ids, how) { ln <- meta[ids, "Cell_Line"]; y <- y[rowSums(is.finite(y)) == ncol(y), ]
  if (how == "ams") { y0 <- y - min(y); sb <- CreateSeuratObject(counts = cs[rownames(y0), ids]); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(y0))
    sb <- suppressWarnings(AddModuleScore(sb, features = lapply(mods, head, 50), name = "Cm_", seed = 1)); x <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", seq_along(CM))])); rownames(x) <- CM; return(x) }
  n <- as.integer(sub("mz", "", how)); out <- matrix(NA_real_, length(CM), length(ids), dimnames = list(CM, ids))
  for (cl in LINES) { j <- which(ln == cl); e <- y[rowMeans(cs[rownames(y), ids[j]]) >= 5, j]; zz <- t(scale(t(e))); zz <- zz[rowSums(is.finite(zz)) == ncol(zz), ]
    for (m in CM) { g <- head(intersect(mods[[m]], rownames(zz)), n); out[m, j] <- colMeans(zz[g, , drop = FALSE]) } }
  out }
bulk_eff <- function(Sx, ids, cov, zline, ref, stat) { ln <- meta[ids, "Cell_Line"]; gp <- meta[ids, "Gene_Control"]
  A <- array(NA_real_, c(length(CM), length(KOS), length(LINES)), dimnames = list(CM, KOS, LINES))
  for (cl in LINES) { j <- which(ln == cl); y <- Sx[, j]; g <- factor(gp[j], levels = GR[GR %in% gp[j]]); G <- model.matrix(~ 0 + g); colnames(G) <- levels(g)
    X <- switch(cov, none = NULL, cc = cbind(q[ids[j], "cc"]), pa = cbind(q[ids[j], "pct_assigned"]), cc_pa = cbind(q[ids[j], "cc"], q[ids[j], "pct_assigned"])); if (!is.null(X)) X <- scale(X, scale = FALSE)
    if (!is.null(X)) y <- limma::removeBatchEffect(y, covariates = X, design = G)
    if (zline) y <- t(scale(t(y)))
    kk <- intersect(KOS, levels(g)); L <- sapply(kk, function(k) { v <- setNames(numeric(nlevels(g)), levels(g)); v[k] <- 1
      if (ref == "ntc") v["NTC"] <- -1 else v[setdiff(levels(g), k)] <- -1 / (nlevels(g) - 1); v })
    if (stat == "diff") { M <- sapply(levels(g), function(l) rowMeans(y[, g == l, drop = FALSE])); A[, kk, cl] <- M %*% L }
    else A[, kk, cl] <- limma::eBayes(limma::contrasts.fit(limma::lmFit(y, G), L))$t }
  A }
metrics <- function(A, S) { per <- sapply(LINES, function(cl) sp(S, A[, , cl])); agg <- apply(A, c(1, 2), mean, na.rm = TRUE)
  ag <- agreement(cl_from_array(A, centre = TRUE, method = "spearman"), permsL); d <- A[cbind(rep(DOWN[, 1], 4), rep(DOWN[, 2], 4), rep(1:4, each = nrow(DOWN)))]; dm <- matrix(d, nrow(DOWN), 4)
  c(screen_rho_mean = mean(per), screen_rho_H2052 = per[["NCI-H2052"]], setNames(per[LINES[-1]], paste0("screen_rho_", c("MSTO", "H2452", "H28"))), screen_rho_agg = sp(S, agg),
    screen_rho_rel = mean(sapply(LINES, function(cl) sp(S[, REL_KO], A[, REL_KO, cl]))), screen_rho_rel_agg = sp(S[, REL_KO], agg[, REL_KO]),
    lines_dmo = ag$global[["dmo"]], lines_z = ag$global[["perm_z"]], lines_p = ag$global[["perm_p"]], lines_nKO05 = sum(ag$per_KO$perm_p < 0.05),
    down_sign = mean(d < 0, na.rm = TRUE), down_all4 = sum(rowSums(dm < 0, na.rm = TRUE) == rowSums(is.finite(dm)) & rowSums(is.finite(dm)) >= 3), down_n = nrow(DOWN), down_mean = mean(agg[DOWN]) ) }

F <- list(set = names(SETS), norm = c("cp10k", "tmm", "tmmq", "rle", "uq"), score = c("ams", "mz30", "mz50", "mz100"), cov = c("none", "cc", "pa", "cc_pa"), zline = c(TRUE, FALSE), ref = c("ntc", "rest"), stat = c("diff", "t"))
rows <- list(); PROF <- list(); t0 <- Sys.time()
for (st in F$set) for (nm in F$norm) { ids <- SETS[[st]]; y <- norm_mat(ids, nm)
  for (sco in F$score) { Sx <- score_mat(y, ids, sco)
    for (cv in F$cov) for (zl in F$zline) for (rf in F$ref) for (sa in F$stat) { A <- bulk_eff(Sx, ids, cv, zl, rf, sa); S <- SCR[[paste(grepl("cc", cv), rf)]]
      key <- paste(st, nm, sco, cv, zl, rf, sa); PROF[[key]] <- A
      rows[[length(rows) + 1]] <- data.frame(set = st, norm = nm, score = sco, cov = cv, zline = zl, ref = rf, stat = sa, t(metrics(A, S))) } }
  cat(sprintf("  %s %s done (%.1f min)\n", st, nm, as.numeric(difftime(Sys.time(), t0, units = "mins")))) }
G <- do.call(rbind, rows); write.csv(G, file.path(OUT, "spec115_grid.csv"), row.names = FALSE); saveRDS(list(profiles = PROF, screen = SCR, down = DOWN), file.path(OUT, "spec115_profiles.rds"))

MET <- c("screen_rho_mean", "screen_rho_H2052", "screen_rho_agg", "screen_rho_rel_agg", "lines_dmo", "lines_z", "lines_nKO05", "down_sign", "down_all4")
marg <- do.call(rbind, lapply(names(F), function(f) do.call(rbind, lapply(F[[f]], function(l) { s <- G[G[[f]] == l, ]; data.frame(factor = f, level = as.character(l), n = nrow(s), t(sapply(MET, function(m) median(s[[m]])))) }))))
write.csv(marg, file.path(OUT, "spec115_marginal_effects.csv"), row.names = FALSE)
cat(sprintf("\n=== %d configurations. Marginal effect of every step (median over all other choices) ===\n", nrow(G))); print(format(marg, digits = 3), row.names = FALSE)
base <- list(set = "S70", norm = "tmm", score = "ams", cov = "none", zline = TRUE, ref = "rest", stat = "diff")
ofat <- do.call(rbind, lapply(names(F), function(f) do.call(rbind, lapply(F[[f]], function(l) { b <- base; b[[f]] <- l; s <- G; for (n in names(b)) s <- s[s[[n]] == b[[n]], ]
  data.frame(changed = f, to = as.character(l), is_baseline = identical(l, base[[f]]), s[, MET]) }))))
write.csv(ofat, file.path(OUT, "spec115_one_factor_at_a_time.csv"), row.names = FALSE)
cat("\n=== one factor at a time from the baseline used in spec109-113 (S70, tmm, ams, no covariate, z within line, KO - rest, diff) ===\n"); print(format(ofat, digits = 3), row.names = FALSE)
cat("\n=== range of each metric over the whole grid (how much can processing move it at all?) ===\n"); print(round(t(sapply(MET, function(m) c(quantile(G[[m]], c(0, 0.05, 0.5, 0.95, 1)), baseline = ofat[ofat$is_baseline, m][1]))), 3))
cat("\n=== top 10 by cross-line agreement (lines_z; the screen plays no part in this ranking), with their screen metrics ===\n")
print(format(head(G[order(-G$lines_z), c(names(F), MET)], 10), digits = 3), row.names = FALSE)
cat("\n=== top 10 by agreement with the screen (screen_rho_agg), with their cross-line metrics ===\n")
print(format(head(G[order(-G$screen_rho_agg), c(names(F), MET)], 10), digits = 3), row.names = FALSE)
cat("\ndo the two criteria agree? Spearman over configurations between lines_z and screen_rho_agg:", round(cor(G$lines_z, G$screen_rho_agg, method = "s"), 3), "; lines_z vs down_sign:", round(cor(G$lines_z, G$down_sign, method = "s"), 3), "\n")

save_fig("fig27_stepwise_pipeline", 15, 8.4, function() {
  par(mfrow = c(2, 4), family = "sans", oma = c(0, 0, 3.2, 0.4), mar = c(6.2, 4.4, 3, 0.8)); show <- c(lines_z = "cross-line agreement (perm z)", screen_rho_agg = "screen vs bulk (4-line mean), rho", screen_rho_H2052 = "screen vs bulk H2052, rho", down_sign = "screen-down pairs: fraction negative in bulk")
  ref0 <- c(lines_z = 0, screen_rho_agg = 0, screen_rho_H2052 = 0, down_sign = 0.5); FC <- setNames(c("#4a3aa7", "#2a78d6", "#1baf7a", "#eda100", "#e87ba4", "#eb6834", "#898781"), names(F))
  for (m in names(show)) { lv <- paste(marg$factor, marg$level, sep = ":"); v <- marg[[m]]
    bp <- barplot(v - ref0[m], col = FC[marg$factor], border = NA, las = 2, names.arg = marg$level, cex.names = 0.6, ylab = show[m], offset = ref0[m], ylim = range(c(v, ref0[m], quantile(G[[m]], c(0.02, 0.98)))))
    abline(h = ref0[m], col = INK, lwd = 0.7); abline(h = median(G[[m]]), lty = 3, col = MUTED); mtext(show[m], side = 3, adj = 0, line = 1.3, font = 2, cex = 0.72); mtext("median over all other choices; dotted = grid median", side = 3, adj = 0, line = 0.3, cex = 0.52, col = INK2)
    if (m == "lines_z") legend("topright", names(F), fill = FC, border = NA, bty = "n", cex = 0.6, ncol = 2) }
  for (m in names(show)) { hist(G[[m]], breaks = 40, col = "#d5d4cf", border = "white", main = "", xlab = show[m], las = 1, cex.axis = 0.8); abline(v = ref0[m], col = INK, lwd = 0.8)
    abline(v = ofat[ofat$is_baseline, m][1], col = "#eb6834", lwd = 1.6); mtext(sprintf("all %d configurations", nrow(G)), side = 3, adj = 0, line = 1.3, font = 2, cex = 0.72); mtext("orange = the baseline of spec109-113; black = no signal", side = 3, adj = 0, line = 0.3, cex = 0.52, col = INK2) }
  mtext("Which processing step loses the reproducibility? Every choice from counts to the KO x Cm effect, on three questions", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.98)
  mtext("top: marginal effect of each level; bottom: how far the whole grid can move each metric. If no level stands out and the grid hugs the no-signal line, processing is not the cause", outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.62, col = INK2)
})
