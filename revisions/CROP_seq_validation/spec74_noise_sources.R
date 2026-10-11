# Spec 74: where does the noise in the cross-line KO agreement come from? S70 sample set, TMM logCPM, NTC - KO moderated t (limma), centred, Spearman unless stated.
#   A variance fractions : per gene (top 5000 variable genes of the line) sequential R2 of technical (% assigned + log depth) -> cell cycle -> replicate number -> KO group;
#                          KO term compared with its chance level (KO labels permuted within the line)
#   B reliability ceiling: per-sample profiles (sample - mean NTC, centred over the KO samples of the line); Spearman between the 2 replicates of the same KO in the same line
#                          vs different KOs in the same line vs the same KO in another line. Cross-line agreement cannot exceed the within-line replicate agreement.
#   C shared NTC term    : mean correlation between DIFFERENT KOs of the same line (NTC - KO t, not centred); with 1 vs 4 NTC in NCI-H2052 (S70 vs S70N4); after centring
#   D read depth         : all samples thinned (binomial) to the depth of the shallowest S70 sample, agreement recomputed (3 seeds)
#   E knockdown          : own-target NTC - KO t of every KO x line vs how well that KO x line profile agrees with the same KO in the other lines
#   F lines              : global agreement with each line left out, and per line pair
# needs: spec71_sample_audit_76.R, spec72_build_profiles.R
# outputs: spec74_variance_fractions.csv, spec74_replicate_reliability.csv, spec74_shared_NTC.csv, spec74_depth_thinning.csv, spec74_knockdown_vs_agreement.csv,
#          spec74_leave_one_line_out.csv, fig7C_noise_sources.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "5000")); perms <- make_perms(N_PERM); set.seed(1)
fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit()
SET <- Sys.getenv("SET", unset = "S70"); OSUF <- if (SET == "S70") "" else paste0("_", SET)   # e.g. SET=RUN3 -> spec74_*_RUN3.csv
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_%s%s.rds", SET, DSUF))); pr4 <- readRDS(file.path(OUT, sprintf("spec72_profiles_%s%s.rds", if (SET == "S70") "S70N4" else "S70", DSUF))); ids <- pr$samples; genes <- pr$genes
A <- pr$P$tmm$none$t; sig <- pr$gene_signal; SP <- list(signal500 = "signal", all_genes = NULL)
glob <- function(A, sig, feats, pairs = PAIRS, lines_avail = LINES) agreement(cl_from_array(A, feats, sig, pairs = pairs, lines_avail = lines_avail), perms, pairs)

## ---------------------------------------------------------------- A variance fractions
r2seq <- function(y, terms) { n <- ncol(y); yc <- y - rowMeans(y); tss <- rowSums(yc^2); X <- matrix(1, n, 1); prev <- 0; out <- list()
  for (nm in names(terms)) { X <- cbind(X, terms[[nm]]); qrx <- qr(X); rss <- colSums(qr.resid(qrx, t(yc))^2); r2 <- 1 - rss / tss; out[[nm]] <- r2 - prev; prev <- r2 }
  out$residual <- 1 - prev; sapply(out, median) }
vf <- do.call(rbind, lapply(LINES, function(cl) { s <- ids[meta[ids, "Cell_Line"] == cl]; y <- normalise(cs[genes, s], "tmm")$y; y <- y[order(matrixStats::rowVars(y), decreasing = TRUE)[1:5000], ]; qq <- q[s, ]
  mm <- function(f) model.matrix(~ f)[, -1, drop = FALSE]; tech <- list(technical = cbind(qq$pct_assigned, log10(qq$assigned_M)), cell_cycle = qq$cc, replicate = mm(factor(qq$rep)))
  obs <- r2seq(y, c(tech, list(KO = mm(factor(qq$KO))))); nul <- rowMeans(replicate(20, r2seq(y, c(tech, list(KO = mm(factor(sample(qq$KO))))))))
  data.frame(line = cl, n = length(s), term = names(obs), median_R2 = obs, chance_R2 = nul, row.names = NULL) }))
write.csv(vf, file.path(OUT, paste0("spec74_variance_fractions", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- B replicate reliability (per-sample profiles)
gs <- rownames(sig)[order(rowMeans(sig), decreasing = TRUE)[1:500]]           # descriptive only (uses all lines)
D <- do.call(cbind, lapply(LINES, function(cl) { s <- ids[meta[ids, "Cell_Line"] == cl]; y <- normalise(cs[genes, s], "tmm")$y; ko <- s[meta[s, "Gene_Control"] != "NTC"]
  d <- rowMeans(y[, s[meta[s, "Gene_Control"] == "NTC"], drop = FALSE]) - y[, ko, drop = FALSE]; SGN * (d - rowMeans(d)) }))
rel <- do.call(rbind, lapply(list(signal500 = gs, all_genes = genes), function(f) { C <- cor(D[f, ], method = "spearman"); ut <- which(upper.tri(C), arr.ind = TRUE); a <- colnames(C)[ut[, 1]]; b <- colnames(C)[ut[, 2]]
  data.frame(a = a, b = b, rho = C[ut], same_KO = meta[a, "Gene_Control"] == meta[b, "Gene_Control"], same_line = meta[a, "Cell_Line"] == meta[b, "Cell_Line"], KO = meta[a, "Gene_Control"], line = meta[a, "Cell_Line"]) }))
rel$space <- sub("\\..*", "", rownames(rel)); rel$class <- with(rel, ifelse(same_line & same_KO, "replicates (same KO, same line)", ifelse(same_line, "different KO, same line", ifelse(same_KO, "same KO, other line", "different KO, other line"))))
write.csv(rel, file.path(OUT, paste0("spec74_replicate_reliability", OSUF, ".csv")), row.names = FALSE)
rel_sum <- aggregate(rho ~ space + class, rel, function(v) c(median = median(v), n = length(v)))

## ---------------------------------------------------------------- C shared NTC term
offd <- function(M) { C <- cor(M, method = "spearman"); mean(C[upper.tri(C)]) }
shared <- do.call(rbind, lapply(list(pr, pr4), function(p) do.call(rbind, lapply(LINES, function(cl) { M <- p$P$tmm$none$t[, , cl]; M <- M[, colSums(!is.na(M)) > 0]
  data.frame(set = p$set, line = cl, n_NTC = p$n[cl, "NTC"], within_line_cross_KO_rho = offd(M), after_centring = offd(M - rowMeans(M)), t_rest = offd(p$P$tmm$none$t_rest[, colnames(M), cl])) }))))
write.csv(shared, file.path(OUT, paste0("spec74_shared_NTC", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- D depth thinning
target <- min(colSums(cs[, ids])); thin <- do.call(rbind, lapply(0:3, function(seed) { x <- cs[, ids]; if (seed > 0) { set.seed(seed); x <- thinCounts(x, target.size = target); dimnames(x) <- dimnames(cs[, ids]) }
  tp <- t_profiles(x, meta, q, ids, genes); do.call(rbind, lapply(names(SP), function(sp) { a <- glob(tp$t, tp$sig, SP[[sp]])
    data.frame(depth = if (seed == 0) "full" else sprintf("thinned to %.1f M", target / 1e6), seed = seed, space = sp, t(a$global), n_KO_perm05 = sum(a$per_KO$perm_p < 0.05)) })) }))
write.csv(thin, file.path(OUT, paste0("spec74_depth_thinning", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- E knockdown vs agreement
CL <- cl_from_array(A, "signal", sig); kd <- do.call(rbind, lapply(KOS, function(k) do.call(rbind, lapply(LINES, function(cl) { i <- which(sapply(PAIRS, function(p) cl %in% p))
  data.frame(KO = k, line = cl, own_target_t = if (k %in% genes) A[k, k, cl] else NA, own_target_log2FC = if (k %in% genes) pr$P$tmm$none$coef[k, k, cl] else NA, n_KO = pr$n[cl, k], n_NTC = pr$n[cl, "NTC"],
             mean_rho_same_KO_other_lines = mean(sapply(CL[i], function(C) C[k, k]), na.rm = TRUE), n_strong_genes = sum(abs(A[, k, cl]) > 4)) }))))
write.csv(kd, file.path(OUT, paste0("spec74_knockdown_vs_agreement", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- F leave one line out / per pair
lolo <- do.call(rbind, lapply(c("none", LINES), function(h) { keep <- setdiff(LINES, h); pp <- PAIRS[sapply(PAIRS, function(p) all(p %in% keep))]
  do.call(rbind, lapply(names(SP), function(sp) { a <- glob(A, sig, SP[[sp]], pp, keep); data.frame(left_out = h, space = sp, n_pairs = length(pp), t(a$global), n_KO_perm05 = sum(a$per_KO$perm_p < 0.05), KOs_perm05 = paste(a$per_KO$KO[a$per_KO$perm_p < 0.05], collapse = ",")) })) }))
pairwise <- do.call(rbind, lapply(names(SP), function(sp) { C <- cl_from_array(A, SP[[sp]], sig); data.frame(space = sp, pair = PLAB, diag_mean = sapply(C, function(x) mean(diag(x), na.rm = TRUE)), off_mean = sapply(C, function(x) mean(x[row(x) != col(x)], na.rm = TRUE))) }))
write.csv(lolo, file.path(OUT, paste0("spec74_leave_one_line_out", OSUF, ".csv")), row.names = FALSE); write.csv(pairwise, file.path(OUT, paste0("spec74_pairwise_lines", OSUF, ".csv")), row.names = FALSE)

options(width = 250); cat("\n=== A variance fractions (median over the top 5000 variable genes) ===\n"); print(vf, digits = 2, row.names = FALSE)
cat("\n=== B replicate reliability (median Spearman of per-sample centred profiles) ===\n"); print(rel_sum, digits = 3)
cat("\n=== C shared NTC term ===\n"); print(shared, digits = 3, row.names = FALSE)
cat("\n=== D depth thinning ===\n"); print(thin, digits = 3, row.names = FALSE)
cat("\n=== E knockdown vs agreement ===\n"); print(kd, digits = 3, row.names = FALSE); ct <- cor.test(kd$own_target_t, kd$mean_rho_same_KO_other_lines, method = "spearman", exact = FALSE)   # KO absent from a line (RUN3): NA row, dropped
cat(sprintf("Spearman(own-target t, agreement) = %.3f, p = %.3g; by line mean agreement: %s\n", ct$estimate, ct$p.value, paste(LINES, round(tapply(kd$mean_rho_same_KO_other_lines, kd$line, mean, na.rm = TRUE)[LINES], 3), collapse = " | ")))
cat("\n=== F leave one line out ===\n"); print(lolo, digits = 3, row.names = FALSE); print(pairwise, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------- figure
save_fig(paste0("fig7C_noise_sources", OSUF), 15, 8.6, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3, 0.5), mar = c(5.2, 4.4, 2.6, 0.8), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.5, font = 2, cex = 0.8)
  tc <- c(technical = "#eda100", cell_cycle = "#1baf7a", replicate = "#4a3aa7", KO = "#eb6834", residual = "#d5d4cf"); M <- sapply(LINES, function(cl) vf$median_R2[vf$line == cl][match(names(tc), vf$term[vf$line == cl])])
  bp <- barplot(M, col = tc, border = "white", las = 1, ylab = "median fraction of variance", names.arg = sub("NCI-|-211H", "", LINES), cex.axis = 0.75, cex.names = 0.8, ylim = c(0, 1.25))
  ch <- sapply(LINES, function(cl) vf$chance_R2[vf$line == cl & vf$term == "KO"]); text(bp, colSums(M) + 0.04, sprintf("KO chance\n%.2f", ch), cex = 0.6, col = INK2); legend("top", names(tc), fill = tc, border = NA, bty = "n", ncol = 3, cex = 0.65); ttl("A  variance fractions per line")
  d <- rel[rel$space == "signal500", ]; cls <- c("replicates (same KO, same line)", "different KO, same line", "same KO, other line", "different KO, other line"); cc4 <- c("#eb6834", "#d5d4cf", "#2a78d6", "#d5d4cf")
  plot(NA, xlim = c(0.4, 4.6), ylim = range(d$rho), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman rho (per-sample profiles)", cex.axis = 0.75); abline(h = 0, col = INK, lty = 2)
  for (i in 1:4) { v <- d$rho[d$class == cls[i]]; boxplot(v, at = i, add = TRUE, boxwex = 0.6, col = adjustcolor(cc4[i], 0.6), border = INK2, outline = FALSE, axes = FALSE); if (length(v) < 200) points(jitter(rep(i, length(v)), amount = 0.12), v, pch = 16, cex = 0.5, col = INK2) }
  text(1:4, par("usr")[3] - 0.04 * diff(par("usr")[3:4]), c("replicates\n(same line)", "other KO\n(same line)", "same KO\n(other line)", "other KO\n(other line)"), xpd = NA, cex = 0.7, adj = c(0.5, 1)); ttl("B  reliability ceiling (top-500 signal genes)")
  x <- shared; xx <- barplot(rbind(x$within_line_cross_KO_rho, x$after_centring), beside = TRUE, col = c("#eb6834", "#2a78d6"), border = NA, las = 1, ylab = "mean rho between DIFFERENT KOs of a line", cex.axis = 0.75, names.arg = rep("", nrow(x)))
  text(colMeans(xx), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), sprintf("%s\n%s, %d NTC", sub("NCI-|-211H", "", x$line), x$set, x$n_NTC), srt = 45, adj = 1, xpd = NA, cex = 0.6); legend("topright", c("NTC - KO as is", "centred over KOs"), fill = c("#eb6834", "#2a78d6"), border = NA, bty = "n", cex = 0.7); ttl("C  shared NTC term")
  d <- thin[thin$space == "signal500", ]; plot(NA, xlim = c(0.5, 2.5), ylim = range(0, thin$perm_z), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "global agreement (permutation z)", cex.axis = 0.75); abline(h = qnorm(0.95), col = MUTED, lty = 3)
  for (sp in names(SP)) { d <- thin[thin$space == sp, ]; points(ifelse(d$depth == "full", 1, 2) + ifelse(sp == "signal500", -0.1, 0.1), d$perm_z, pch = 21, bg = ifelse(sp == "signal500", "#2a78d6", "#eda100"), col = "white", cex = 1.6) }
  axis(1, at = 1:2, labels = c("full depth", sprintf("thinned to %.1f M", target / 1e6)), tick = FALSE, cex.axis = 0.8); legend("bottomleft", names(SP), pt.bg = c("#2a78d6", "#eda100"), pch = 21, col = "white", bty = "n", cex = 0.7); ttl("D  read depth")
  plot(kd$own_target_t, kd$mean_rho_same_KO_other_lines, pch = 21, bg = LINE_COL[kd$line], col = "white", cex = 1.7, bty = "n", las = 1, xlab = sprintf("own target gene, %s moderated t", DLAB), ylab = "mean rho with the same KO in the other lines", cex.axis = 0.75, cex.lab = 0.8)
  abline(h = 0, v = 0, col = GRID); text(kd$own_target_t, kd$mean_rho_same_KO_other_lines, kd$KO, pos = 3, cex = 0.5, col = INK2); legend("topleft", LINES, pt.bg = LINE_COL[LINES], pch = 21, col = "white", bty = "n", cex = 0.65); ttl(sprintf("E  knockdown vs agreement (rho = %.2f, p = %.2g)", ct$estimate, ct$p.value))
  d <- lolo[lolo$space == "signal500", ]; bp <- barplot(d$perm_z, col = c(INK2, LINE_COL[LINES]), border = NA, las = 1, ylab = "global agreement (permutation z)", cex.axis = 0.75, names.arg = rep("", nrow(d))); abline(h = qnorm(0.95), col = MUTED, lty = 3)
  text(bp, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), c("all 4 lines", paste("without", sub("NCI-|-211H", "", LINES))), srt = 45, adj = 1, xpd = NA, cex = 0.7); text(bp, d$perm_z, sprintf("%d KO", d$n_KO_perm05), pos = 3, cex = 0.65, xpd = NA); ttl("F  leaving one line out")
  mtext(sprintf("Sources of noise in the cross-line agreement of bulk KO profiles (%s, %s)", SET, DLAB), outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1)
  mtext("TMM logCPM, limma moderated t; agreement = Spearman between centred KO profiles of two lines, KO labels permuted within line for the null; dotted = one-sided 5%", outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.62, col = INK2)
})
