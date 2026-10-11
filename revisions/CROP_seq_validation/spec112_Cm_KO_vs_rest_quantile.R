# Spec 112: the spec109 figure (KO vs rest, replicates averaged, 4 cell lines aggregated) recomputed on QUANTILE-NORMALISED expression,
# side by side with the current TMM scores.
#
# Normalisation arms, following the recipe already used in spec103:
#   tmm   TMM logCPM computed per cell line (as spec72 does)                                   [current]
#   tmmq  the same matrix, then limma::normalizeQuantiles across all 70 bulk samples, so every sample is forced to an identical
#         distribution. This removes any residual difference in the shape of the expression distribution between samples/lines, which
#         is the thing most likely to be adding spread to the module scores and so blurring the KO-vs-rest separation.
# Cm scores are rebuilt from scratch on each matrix with Seurat::AddModuleScore over the same top-50 genes per module (cm_top50), so the
# two arms differ ONLY in the normalisation.
#
# Everything downstream is identical to spec109: z per module within line -> replicate means per group -> the KO contributes one value
# per line (4), the rest contributes the 8 non-target groups x 4 lines (32) -> group-label permutation within line (10,000), two-sided.
# needs: repro_common.R, profiles.R, agree_common.R, spec_common.R, spec72_profiles_S70.rds
# outputs: spec112_quantile_vs_tmm_stats.csv, spec112_quantile_group_means.csv, fig24_Cm_KO_vs_rest_quantile.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3")
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000")); GR <- c("NTC", KOS)
PAIRS <- list(c("TEAD4", "Cm20"), c("TWIST1", "Cm8"), c("TWIST1", "Cm9"), c("TCF3", "Cm10"), c("TEAD2", "Cm10"), c("TEAD2", "Cm1"))

pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs
ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]
top50 <- cm_top50(); CM <- names(top50)
det <- rownames(cs)[rowSums(cs[, ids]) > 0]
score <- function(y) { y <- y - min(y); sb <- CreateSeuratObject(counts = cs[rownames(y), ids])
  sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(y))
  sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
  x <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", seq_along(CM))])); rownames(x) <- CM; x }
qn <- function(y) { z <- limma::normalizeQuantiles(y); dimnames(z) <- dimnames(y); z }
yT <- do.call(cbind, lapply(LINES, function(cl) normalise(cs[det, ids[line == cl]], "tmm")$y))[, ids]   # TMM per line, as spec72
cat("sample medians of TMM logCPM before quantile normalisation (range over 70 samples):", round(range(apply(yT, 2, median)), 3), "\n")
## NOT identical: logCPM has a large block of tied values at the detection floor, and normalizeQuantiles averages ties,
## so columns with different numbers of tied zeros keep slightly different medians. The spread narrows but does not vanish.
cat("after quantile normalisation:", round(range(apply(qn(yT), 2, median)), 3), "  (ties at the detection floor keep these from collapsing to one value)\n")
SC <- list(tmm = score(yT), tmmq = score(qn(yT)))
cat("\ncorrelation of the 20 module scores between the two normalisations (per module, over 70 samples):\n")
print(round(sapply(CM, function(m) cor(SC$tmm[m, ], SC$tmmq[m, ])), 3))

## within-line z -> replicate means per group, for each arm
Gof <- function(S) { Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j, drop = FALSE]))) }; Z[!is.finite(Z)] <- NA
  G <- array(NA_real_, c(nrow(Z), length(GR), length(LINES)), dimnames = list(rownames(Z), GR, LINES))
  for (cl in LINES) for (g in GR) { j <- which(line == cl & grp == g); if (length(j)) G[, g, cl] <- rowMeans(Z[, j, drop = FALSE], na.rm = TRUE) }; G }
GL <- lapply(SC, Gof)
write.csv(do.call(rbind, lapply(names(GL), function(nm) do.call(rbind, lapply(LINES, function(cl)
  data.frame(norm = nm, line = cl, group = GR, t(round(GL[[nm]][CM, , cl], 4))))))),
  file.path(OUT, "spec112_quantile_group_means.csv"), row.names = FALSE)

set.seed(1)
one <- function(G, cm, ko) { a <- G[cm, ko, ]; b <- as.vector(G[cm, setdiff(GR, ko), ])
  a <- a[is.finite(a)]; b <- b[is.finite(b)]; obs <- mean(a) - mean(b)
  null <- replicate(N_PERM, { d <- sapply(LINES, function(cl) { v <- G[cm, , cl]; v <- v[is.finite(v)]
      j <- sample(length(v), 1); c(v[j], mean(v[-j])) }); mean(d[1, ]) - mean(d[2, ]) })
  c(delta = obs, p = (1 + sum(abs(null) >= abs(obs))) / (N_PERM + 1),
    wilcox = wilcox.test(a, b)$p.value, n_same_sign = max(sum(G[cm, ko, ] > 0, na.rm = TRUE), sum(G[cm, ko, ] < 0, na.rm = TRUE))) }
res <- do.call(rbind, lapply(PAIRS, function(pp) { ko <- pp[1]; cm <- pp[2]
  t_ <- one(GL$tmm, cm, ko); q_ <- one(GL$tmmq, cm, ko)
  data.frame(KO = ko, Cm = cm, delta_tmm = t_[1], p_tmm = t_[2], wilcox_tmm = t_[3],
             delta_q = q_[1], p_q = q_[2], wilcox_q = q_[3], sign_lines_q = q_[4]) }))
res$p_tmm_fdr <- p.adjust(res$p_tmm, "BH"); res$p_q_fdr <- p.adjust(res$p_q, "BH")
write.csv(res, file.path(OUT, "spec112_quantile_vs_tmm_stats.csv"), row.names = FALSE)
options(width = 220); cat("\n=== KO vs rest (replicates averaged, 4 lines aggregated): TMM vs TMM+quantile ===\n")
print(format(res, digits = 3), row.names = FALSE)
cat("\nmedian permutation p:  tmm", round(median(res$p_tmm), 3), " tmmq", round(median(res$p_q), 3), "\n")
cat("pairs with p < 0.05:   tmm", sum(res$p_tmm < 0.05), " tmmq", sum(res$p_q < 0.05), "\n")
cat("pairs improved by quantile normalisation:", sum(res$p_q < res$p_tmm), "of", nrow(res), "\n")

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
LCOL <- LINE_COL[LINES]; short <- sub("NCI-|-211H", "", LINES)
save_fig("fig24_Cm_KO_vs_rest_quantile", 11, 7.2, function() {
  par(mfrow = c(2, 3), oma = c(3.2, 1.2, 4.0, 0.6), mar = c(3.0, 4.3, 3.8, 0.8), family = "sans")
  G <- GL$tmmq
  for (i in seq_len(nrow(res))) { ko <- res$KO[i]; cm <- res$Cm[i]
    a <- G[cm, ko, ]; b <- as.vector(G[cm, setdiff(GR, ko), ]); bl <- rep(LINES, each = length(GR) - 1)
    yr <- range(c(a, b), na.rm = TRUE) + c(-0.12, 0.42) * diff(range(c(a, b), na.rm = TRUE))
    boxplot(list(rest = b, KO = a), ylim = yr, col = c("#d5d4cf", adjustcolor(KO_COL[[ko]], 0.75)), border = INK2,
            outline = FALSE, boxwex = 0.55, las = 1, bty = "n", cex.axis = 0.8,
            names = c("rest\n(8 groups x 4 lines)", sprintf("%s\n(1 per line)", ko)),
            ylab = sprintf("%s score (quantile-norm, z within line)", cm))
    abline(h = 0, lty = 3, col = GRID)
    points(jitter(rep(1, length(b)), amount = 0.16), b, pch = 21, bg = adjustcolor(LCOL[bl], 0.6), col = INK2, cex = 0.7)
    points(seq(1.88, 2.12, length.out = length(a)), a, pch = 21, bg = LCOL[LINES], col = INK2, cex = 1.5)
    st <- ifelse(res$p_q_fdr[i] < 0.01, " **", ifelse(res$p_q_fdr[i] < 0.05, " *", ""))
    mtext(sprintf("%s in %s KO", cm, ko), side = 3, adj = 0, line = 2.2, font = 2, cex = 0.82, col = KO_COL[[ko]])
    mtext(sprintf("quantile: delta = %+.2f, p = %s%s", res$delta_q[i], fp(res$p_q[i]), st), side = 3, adj = 0, line = 1.2,
          cex = 0.6, col = INK2, font = ifelse(res$p_q[i] < 0.05, 2, 1))
    mtext(sprintf("TMM (current): delta = %+.2f, p = %s", res$delta_tmm[i], fp(res$p_tmm[i])), side = 3, adj = 0, line = 0.4, cex = 0.55, col = MUTED) }
  mtext("Cm module score after QUANTILE NORMALISATION: the KO vs the rest, replicates averaged, 4 cell lines aggregated", outer = TRUE, side = 3, adj = 0, line = 2.0, font = 2, cex = 0.92)
  mtext(sprintf("TMM logCPM per line, then limma::normalizeQuantiles across all 70 samples; module scores rebuilt with AddModuleScore (top-50 genes); group-label permutation within line (%s perms); grey line = the current TMM result", format(N_PERM, big.mark = ",")),
        outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.56, col = INK2)
  legend("bottom", short, pt.bg = LCOL, pch = 21, col = INK2, horiz = TRUE, bty = "n", xpd = NA, inset = c(0, -0.30), cex = 0.78, pt.cex = 1.2)
})
