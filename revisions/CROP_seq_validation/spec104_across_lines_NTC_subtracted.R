# Spec 104: fig13 panel C (bulk across the 4 lines, top-500 KO-signal genes, Spearman, no covariate) with plain NTC-subtracted profiles, i.e. NTC - KO WITHOUT the centring over the 8 KOs.
#   centred (fig13C) : limma t of NTC - KO, minus the mean over the 8 KOs of the line (the NTC cancels exactly)
#   NTC - KO, t      : limma moderated t of NTC - KO, as is
#   NTC - KO, logFC  : limma coefficient NTC - KO (log2 fold change), as is
#   mean of 9 groups - KO (t / logFC): no NTC subtraction; each KO against the mean of the 9 group means (NTC + 8 KOs) of its line
# uncentred profiles of the different KOs of a line share the NTC, so they correlate with each other (grey boxes shift up); the KO-label permutation p accounts for it
# needs: spec72 (S70), spec82
# outputs: spec104_across_lines_NTC_subtracted.csv, spec104_across_lines_NTC_subtracted_per_KO.csv, fig18_across_lines_NTC_subtracted.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); P <- pr$P$tmm$none
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
# no NTC subtraction: each KO against the mean of all 9 groups (NTC + 8 KOs) of the line -- limma contrast (mean of the 9 group means) - KO, same sign convention as NTC - KO
fab <- load_full_audited_bulk(); meta <- fab$meta; ids <- pr$samples; GR <- c("NTC", KOS); T9 <- F9 <- array(NA_real_, dim(P$t), dimnames = dimnames(P$t))
for (cl in LINES) { s <- ids[meta[ids, "Cell_Line"] == cl]; y <- normalise(fab$cs[pr$genes, s], "tmm")$y; g <- factor(meta[s, "Gene_Control"], levels = GR); des <- model.matrix(~ 0 + g); colnames(des) <- GR
  L <- sapply(KOS, function(k) { v <- setNames(rep(1 / 9, 9), GR); v[k] <- v[k] - 1; v }); ft <- limma::eBayes(limma::contrasts.fit(limma::lmFit(y, des), L)); T9[, , cl] <- ft$t; F9[, , cl] <- ft$coefficients }
V <- list(g9_t = list(A = T9, centre = FALSE, lab = "mean of 9 groups - KO, t"), g9_fc = list(A = F9, centre = FALSE, lab = "mean of 9 groups - KO, log fold change"), centred = list(A = P$t, centre = TRUE, lab = "KO-centred t (fig13C)"), ntc_t = list(A = P$t, centre = FALSE, lab = "NTC - KO, t"), ntc_fc = list(A = P$coef, centre = FALSE, lab = "NTC - KO, log fold change"))
AG <- lapply(V, function(v) agreement(cl_from_array(v$A, "signal", pr$gene_signal, centre = v$centre), permsL))
G <- do.call(rbind, lapply(names(V), function(n) data.frame(view = n, t(AG[[n]]$global), n_KO_p05 = sum(AG[[n]]$per_KO$perm_p < 0.05), n_KO_fdr05 = sum(AG[[n]]$per_KO$perm_fdr < 0.05)))); K <- do.call(rbind, lapply(names(V), function(n) cbind(view = n, AG[[n]]$per_KO)))
write.csv(G, file.path(OUT, "spec104_across_lines_NTC_subtracted.csv"), row.names = FALSE); write.csv(K, file.path(OUT, "spec104_across_lines_NTC_subtracted_per_KO.csv"), row.names = FALSE)
options(width = 200); print(G[, c("view", "diag_mean", "off_mean", "dmo", "perm_p", "perm_z", "n_KO_p05", "n_KO_fdr05")], digits = 3, row.names = FALSE)
print(reshape(transform(K[, c("KO", "view")], v = sprintf("%.4f (r %.2f)%s", K$perm_p, K$mean_intra, ifelse(K$perm_fdr < 0.05, " *", ""))), idvar = "KO", timevar = "view", direction = "wide"), row.names = FALSE)
o2 <- read.csv(file.path(OUT, "spec82_TF_agreement_ranked.csv"))$KO
save_fig("fig18_across_lines_NTC_subtracted", 16, 10.6, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3, 0.4), family = "sans")
  for (n in names(V)) { a <- AG[[n]]; k <- a$per_KO; yr <- range(c(unlist(a$cross), a$intra)); yr <- yr + c(-0.02, 0.17) * diff(yr); par(mar = c(4.6, 4.6, 3.4, 1))
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(o2)) { j <- match(o2[i], k$KO); v <- a$intra[o2[i], ]; boxplot(a$cross[[j]], at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
      boxplot(v, at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[o2[i]]], 0.65), border = INK2, outline = FALSE, axes = FALSE); points(i + 0.2 + seq(-0.1, 0.1, length.out = length(v)), v, pch = 16, cex = 0.75, col = INK2)
      text(i, yr[2], sprintf("p=%s\n%s", fp(k$perm_p[j]), ifelse(k$perm_fdr[j] < 0.01, "**", ifelse(k$perm_fdr[j] < 0.05, "*", ""))), cex = 0.68, adj = c(0.5, 1), font = ifelse(k$perm_p[j] < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.82, font = 2, col = KO_COL[o2])
    mtext(V[[n]]$lab, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.82); mtext(sprintf("same KO %.3f vs other KOs %.3f; global p = %s; * / ** FDR < 0.05 / 0.01", a$global[["diag_mean"]], a$global[["off_mean"]], fp(a$global[["perm_p"]])), side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2) }
  mtext("Bulk across the 4 lines: each KO vs the mean of all 9 groups (no NTC subtraction), vs the KO-centred profiles of fig13C, vs plain NTC subtraction", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.95)
  mtext("top-500 KO-signal genes from the 2 lines not compared; colour = same KO (6 line pairs), grey = other KOs; 20,000 KO-label permutations within lines; S70, no covariate", outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2)
})
