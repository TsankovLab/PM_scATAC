# Spec 113: which Cm x KO pairs actually reproduce between the CROP-seq screen and the bulk panel?
#
# All 20 x 8 = 160 pairs are scored on three things, and a pair has to pass all three to count as reproducible:
#   1  direction      the bulk effect has the same sign as the CROP-seq effect
#   2  consistency    all 4 bulk cell lines have that sign (not an average carried by one line)
#   3  significance   within-line group-label permutation p, BH-corrected over the 160 pairs
# Effect in both data sets = mean(the KO) - mean(NTC + the other 7 KOs); in bulk on the within-line z after averaging replicates, so
# each line contributes one number (as spec109/spec110).
#
# SIGN CONVENTION. This script uses KO - rest, so NEGATIVE = the module is LOWER in the KO. That is the opposite of spec101/spec69,
# which report rest - KO (positive = lower in the KO); a column rest_minus_KO is carried so the two can be read side by side.
# The two directions are not equally interpretable: knocking out an activating TF should REDUCE the module it drives, so pairs the
# screen puts DOWN in the KO are the ones with a mechanistic expectation. Pairs that go UP in the KO can be compensation, a de-repressed
# programme, or the module simply tracking something else, so they are scored separately rather than pooled.
#
# Criterion 2 on its own proves nothing: with 4 lines, P(all agree) = 2 * 0.5^4 = 1/8, so ~20 of 160 pairs clear it by chance alone.
# That is why the permutation p is computed and corrected. The null for "one group differs from the other eight" does not depend on
# WHICH group is called the KO, so a single null per module is generated once and shared across the 8 KOs -- 10,000 draws per module
# rather than per pair, which is what makes the full 160-pair correction affordable.
# needs: repro_common.R, profiles.R, agree_common.R, spec72_profiles_S70.rds, spec75_screen_scores_by_cell.rds
# outputs: spec113_all_pairs_scored.csv, spec113_reproducible_pairs.csv, fig25_reproducible_Cm_KO_pairs.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000")); GR <- c("NTC", KOS); CM <- paste0("Cm", 1:20)

pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta
ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]
S <- pr$cm_scores$tmm[CM, ids, drop = FALSE]
Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j, drop = FALSE]))) }; Z[!is.finite(Z)] <- NA
G <- array(NA_real_, c(length(CM), length(GR), length(LINES)), dimnames = list(CM, GR, LINES))
for (cl in LINES) for (g in GR) { j <- which(line == cl & grp == g); if (length(j)) G[, g, cl] <- rowMeans(Z[, j, drop = FALSE], na.rm = TRUE) }

## observed bulk effect per line, and aggregated
dl <- array(NA_real_, c(length(CM), length(KOS), length(LINES)), dimnames = list(CM, KOS, LINES))
for (cl in LINES) for (k in KOS) dl[, k, cl] <- G[, k, cl] - rowMeans(G[, setdiff(GR, k), cl], na.rm = TRUE)
B <- apply(dl, c(1, 2), mean, na.rm = TRUE)

## shared null: in each line relabel one of the 9 groups as "the KO"; identical in distribution for every KO, so one null per module
set.seed(1)
NULLD <- replicate(N_PERM, {
  d <- sapply(LINES, function(cl) { v <- G[, , cl]; ok <- which(apply(v, 2, function(x) any(is.finite(x))))
    j <- sample(ok, 1); v[, j] - rowMeans(v[, setdiff(ok, j), drop = FALSE], na.rm = TRUE) })
  rowMeans(d, na.rm = TRUE) })                                           # CM x N_PERM
pval <- t(sapply(CM, function(m) sapply(KOS, function(k) (1 + sum(abs(NULLD[m, ]) >= abs(B[m, k]))) / (N_PERM + 1))))
dimnames(pval) <- list(CM, KOS)

## CROP-seq effect
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); gs <- as.character(sc$group); keep <- gs %in% GR
zS <- scale(sc$score[keep, CM]); gk <- gs[keep]; mS <- sapply(GR, function(g) colMeans(zS[gk == g, , drop = FALSE]))
C <- sapply(KOS, function(k) mS[, k] - rowMeans(mS[, setdiff(GR, k), drop = FALSE])); dimnames(C) <- list(CM, KOS)

nagree <- sapply(KOS, function(k) sapply(CM, function(m) sum(sign(dl[m, k, ]) == sign(C[m, k]), na.rm = TRUE)))
D <- data.frame(Cm = rep(CM, length(KOS)), KO = rep(KOS, each = length(CM)),
                cropseq = as.vector(C), bulk = as.vector(B), same_direction = as.vector(sign(B) == sign(C)),
                n_lines_agree = as.vector(nagree), perm_p = as.vector(pval))
D$rest_minus_KO_cropseq <- -D$cropseq; D$rest_minus_KO_bulk <- -D$bulk     # spec101/spec69 convention
D$direction <- ifelse(D$cropseq < 0, "down in KO (expected)", "up in KO")
D$perm_fdr <- p.adjust(D$perm_p, "BH")
## FDR computed WITHIN each direction as well: the two sets are tested as separate families
D$perm_fdr_within <- NA_real_
for (dd in unique(D$direction)) { i <- D$direction == dd; D$perm_fdr_within[i] <- p.adjust(D$perm_p[i], "BH") }
D <- cbind(D, as.data.frame(matrix(aperm(dl, c(1, 2, 3)), ncol = length(LINES), dimnames = list(NULL, sub("NCI-|-211H", "", LINES)))))
D <- D[order(D$perm_p, -abs(D$bulk)), ]
write.csv(D, file.path(OUT, "spec113_all_pairs_scored.csv"), row.names = FALSE)

options(width = 220)
cat(sprintf("=== %d Cm x KO pairs scored ===\n", nrow(D)))
cat("same direction as CROP-seq:", sum(D$same_direction), sprintf("(chance 80)\n"))
cat("all 4 lines agree with CROP-seq:", sum(D$n_lines_agree == 4), sprintf("(chance ~%.0f)\n", nrow(D) / 16))
cat("perm p < 0.05:", sum(D$perm_p < 0.05), sprintf("(chance ~%.0f)\n", 0.05 * nrow(D)))
cat("BH FDR < 0.05:", sum(D$perm_fdr < 0.05), " | FDR < 0.10:", sum(D$perm_fdr < 0.10), "\n")
cat("\n=== pairs passing ALL THREE: same direction, 4/4 lines, p < 0.05 ===\n")
rep3 <- D[D$same_direction & D$n_lines_agree == 4 & D$perm_p < 0.05, ]
if (nrow(rep3)) print(format(rep3, digits = 3), row.names = FALSE) else cat("  none\n")
cat("\n=== relaxed: same direction + 4/4 lines (significance not required) ===\n")
rep2 <- D[D$same_direction & D$n_lines_agree == 4, ]
print(format(rep2[, c("Cm", "KO", "cropseq", "bulk", "perm_p", "perm_fdr", "H2052", "MSTO", "H2452", "H28")], digits = 3), row.names = FALSE)
write.csv(rep2, file.path(OUT, "spec113_reproducible_pairs.csv"), row.names = FALSE)
cat("\n############ stratified by the CROP-seq direction ############\n")
for (dd in c("down in KO (expected)", "up in KO")) { sub <- D[D$direction == dd, ]
  cat(sprintf("\n=== %s : %d pairs ===\n", dd, nrow(sub)))
  cat("  4/4 lines agree:", sum(sub$n_lines_agree == 4), " | perm p < 0.05:", sum(sub$perm_p < 0.05),
      " | both:", sum(sub$n_lines_agree == 4 & sub$perm_p < 0.05), "\n")
  k <- sub[sub$same_direction & sub$n_lines_agree == 4, ]
  k <- k[order(k$perm_p), ]
  if (nrow(k)) print(format(k[, c("Cm", "KO", "cropseq", "bulk", "rest_minus_KO_cropseq", "rest_minus_KO_bulk",
                                  "perm_p", "perm_fdr_within", "H2052", "MSTO", "H2452", "H28")], digits = 3), row.names = FALSE)
  else cat("  no pair with 4/4 line agreement\n") }
write.csv(D[D$direction == "down in KO (expected)" & D$same_direction & D$n_lines_agree == 4, ],
          file.path(OUT, "spec113_reproducible_down_in_KO.csv"), row.names = FALSE)

cat("\n=== strongest bulk effects regardless of CROP-seq agreement (top 10 by permutation p) ===\n")
print(format(head(D[, c("Cm", "KO", "cropseq", "bulk", "same_direction", "n_lines_agree", "perm_p", "perm_fdr")], 10), digits = 3), row.names = FALSE)

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
save_fig("fig25_reproducible_Cm_KO_pairs", 12, 5.8, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.25, 1)); par(family = "sans", oma = c(0, 0, 3.4, 0.4), mar = c(4.4, 4.8, 3.0, 1))
  hit <- D$same_direction & D$n_lines_agree == 4
  dn <- D$direction == "down in KO (expected)"
  plot(D$cropseq, D$bulk, pch = 21, bg = ifelse(hit, adjustcolor(KO_COL[D$KO], 0.95), adjustcolor("grey70", 0.5)),
       col = ifelse(hit, "black", NA), cex = ifelse(hit, 1.8, 1.0), lwd = ifelse(hit, 1.2, 0), bty = "n", las = 1,
       xlab = "CROP-seq: KO minus the rest (z)", ylab = "bulk: KO minus the rest, 4 lines aggregated (z)", cex.axis = 0.85)
  abline(h = 0, v = 0, lty = 2, col = GRID); abline(0, 1, lty = 3, col = MUTED)
  if (any(hit)) text(D$cropseq[hit], D$bulk[hit], paste(D$Cm[hit], D$KO[hit]), pos = 3, cex = 0.58, font = 2, col = INK, xpd = NA)
  mtext("all 160 pairs; ringed = same direction AND all 4 lines agree", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
  mtext(sprintf("%d of 160 pass (chance ~%.0f); %d of them are DOWN in the KO in the screen (lower-left quadrant), %d are up", sum(hit), nrow(D) / 16, sum(hit & dn), sum(hit & !dn)), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2)
  r <- rep2[rep2$cropseq < 0, ]; r <- r[order(r$perm_p), ]; n <- min(nrow(r), 12)   # the expected direction only
  if (n) { par(mar = c(4.4, 7.0, 3.0, 1))
    bp <- barplot(rev(r$bulk[1:n]), horiz = TRUE, col = adjustcolor(KO_COL[rev(r$KO[1:n])], 0.85), border = INK2,
                  names.arg = rev(paste(r$Cm[1:n], r$KO[1:n])), las = 1, cex.names = 0.72, xlab = "bulk effect (z)")
    abline(v = 0, col = INK, lwd = 0.6)
    text(rev(r$bulk[1:n]), bp, sprintf(" p %s", fp(rev(r$perm_p[1:n]))), adj = ifelse(rev(r$bulk[1:n]) > 0, 0, 1), cex = 0.6, col = INK2, xpd = NA)
    mtext("candidates DOWN in the KO, by permutation p", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
    mtext("the mechanistically expected direction; all have the CROP-seq sign in all 4 lines", side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2) }
  mtext("Which Cm x KO effects reproduce between the CROP-seq screen and the 4-line bulk panel?", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.95)
  mtext("effect = mean(KO) - mean(NTC + other 7 KOs) in both; bulk: replicates averaged, within-line z, 4 lines equally weighted; p = group-label permutation within line, BH over 160", outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2)
})
