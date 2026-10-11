# Spec 114: how reproducible can the KO x Cm effects be? Reliability of each data set ON ITS OWN, and the ceiling that puts on any agreement between them.
# Two measurements cannot correlate better than sqrt(reliability_1 * reliability_2). If the screen's own KO-specific module effects do not reproduce between two halves of
# its cells, or a bulk line's do not reproduce between its replicates, no normalisation of the other side can recover agreement.
#   effect, both sides: module score of the KO minus a reference, 20 Cm x 8 KOs = 160 values.
#       ntc      = KO - NTC                       (contains the response shared by every KO)
#       specific = KO - mean of the other groups  (what distinguishes this KO; the quantity used in spec109-113)
#       shared   = the per-module mean of ntc over the 8 KOs (20 values)
#   A screen : cells of every group split at random into two halves (N_SPLIT times); half-vs-half Spearman over the 160 values, Spearman-Brown corrected to the full
#              data set; the same per KO over its 20 modules. Group sizes are very unequal (13 - 1722 cells), so the per-KO numbers are the point.
#   B bulk   : within a line, replicate A vs replicate B of every group (groups with one replicate drop out); same statistics.
#              CAVEAT: both replicates of a KO share the i5 library block (spec79), so this is an UPPER bound -- batch that is nested in KO counts as "KO signal" here.
#   C ceiling: sqrt(rel_screen * rel_bulk) per line, next to the correlation actually observed.
#   D modules: are the Cm programmes (defined in patient scRNA) coherent in the cell lines? Mean pairwise Spearman among the top-50 genes of a module across the samples
#              of a line (group means removed, so it is co-variation, not KO effect) against 200 random gene sets matched for expression.
# needs: spec72 (S70), spec75 (screen scores by cell)
# outputs: spec114_screen_reliability.csv, spec114_bulk_reliability.csv, spec114_ceiling.csv, spec114_module_coherence.csv, fig26_reliability_ceiling.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_SPLIT <- as.integer(Sys.getenv("N_SPLIT", unset = "300")); GR <- c("NTC", KOS); CM <- paste0("Cm", 1:20); options(width = 220); set.seed(1)
sb <- function(r) 2 * r / (1 + r)                                   # Spearman-Brown: half -> full length
eff <- function(M) { ntc <- M[, KOS, drop = FALSE] - M[, "NTC"]     # M = Cm x group means
  spec <- sapply(KOS, function(k) M[, k] - rowMeans(M[, setdiff(GR, k), drop = FALSE], na.rm = TRUE))
  list(ntc = ntc, specific = spec, shared = rowMeans(ntc, na.rm = TRUE)) }
sp <- function(a, b) suppressWarnings(cor(as.vector(a), as.vector(b), method = "spearman", use = "pairwise.complete.obs"))

## ---------------------------------------------------------------- A screen
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); keep <- sc$group %in% GR; zS <- scale(sc$score[keep, ]); colnames(zS) <- CM; gS <- sc$group[keep]
nS <- table(gS)[GR]; gm <- function(idx) sapply(GR, function(g) colMeans(zS[idx & gS == g, , drop = FALSE]))
ES <- eff(gm(rep(TRUE, length(gS))))
half <- replicate(N_SPLIT, { h <- logical(length(gS)); for (g in GR) { i <- which(gS == g); h[sample(i, floor(length(i) / 2))] <- TRUE }
  a <- eff(gm(h)); b <- eff(gm(!h))
  c(ntc = sp(a$ntc, b$ntc), specific = sp(a$specific, b$specific), shared = sp(a$shared, b$shared),
    setNames(sapply(KOS, function(k) sp(a$specific[, k], b$specific[, k])), paste0("spec_", KOS)),
    setNames(sapply(KOS, function(k) sp(a$ntc[, k], b$ntc[, k])), paste0("ntc_", KOS))) })
scr <- data.frame(what = rownames(half), half_rho = rowMeans(half), half_sd = apply(half, 1, sd)); scr$reliability <- sb(pmax(scr$half_rho, 0))
scr$n_cells <- NA; for (k in KOS) scr$n_cells[grepl(paste0("_", k, "$"), scr$what)] <- nS[k]
write.csv(scr, file.path(OUT, "spec114_screen_reliability.csv"), row.names = FALSE)
cat("=== A  CROP-seq screen: split-half reliability of the KO x Cm effects (", N_SPLIT, "random splits; Spearman-Brown to full size) ===\n"); print(format(scr, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- B bulk
pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; ids <- pr$samples
line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]; S <- pr$cm_scores$tmm[CM, ids]
Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j]))) }
EB <- list(); bulk <- list()
for (cl in LINES) { j <- line == cl
  M <- sapply(GR, function(g) rowMeans(Z[, j & grp == g, drop = FALSE])); EB[[cl]] <- eff(M)
  A <- B <- matrix(NA_real_, length(CM), length(GR), dimnames = list(CM, GR))
  for (g in GR) { s <- which(j & grp == g); if (length(s) >= 2) { A[, g] <- Z[, s[1]]; B[, g] <- if (length(s) == 2) Z[, s[2]] else rowMeans(Z[, s[-1], drop = FALSE]) } }
  ok <- GR[colSums(is.finite(A)) > 0]; a <- eff(A); b <- eff(B)
  per <- sapply(KOS, function(k) if (k %in% ok) sp(a$specific[, k], b$specific[, k]) else NA)
  bulk[[cl]] <- data.frame(line = cl, groups_with_2_reps = length(ok), NTC_split = "NTC" %in% ok,
    half_ntc = if ("NTC" %in% ok) sp(a$ntc, b$ntc) else NA, half_specific = sp(a$specific, b$specific), t(setNames(per, paste0("spec_", KOS)))) }
bulk <- do.call(rbind, bulk); bulk$rel_specific <- sb(pmax(bulk$half_specific, 0)); bulk$rel_ntc <- sb(pmax(bulk$half_ntc, 0))
write.csv(bulk, file.path(OUT, "spec114_bulk_reliability.csv"), row.names = FALSE)
cat("\n=== B  bulk: replicate-vs-replicate reliability within each line (upper bound: batch is nested in KO) ===\n"); print(format(bulk, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- C ceiling vs observed
relS <- setNames(scr$reliability, scr$what)
cei <- do.call(rbind, lapply(LINES, function(cl) { b <- bulk[bulk$line == cl, ]
  data.frame(line = cl, kind = c("specific (KO - rest)", "NTC-referenced (KO - NTC)", "shared (module mean over KOs)"),
    rel_screen = c(relS["specific"], relS["ntc"], relS["shared"]), rel_bulk = c(b$rel_specific, b$rel_ntc, NA),
    observed = c(sp(ES$specific, EB[[cl]]$specific), sp(ES$ntc, EB[[cl]]$ntc), sp(ES$shared, EB[[cl]]$shared))) }))
cei$ceiling <- sqrt(cei$rel_screen * cei$rel_bulk); cei$fraction_of_ceiling <- cei$observed / cei$ceiling
write.csv(cei, file.path(OUT, "spec114_ceiling.csv"), row.names = FALSE)
cat("\n=== C  observed screen-vs-bulk Spearman against the reliability ceiling ===\n"); print(format(cei, digits = 3), row.names = FALSE)
xl <- combn(LINES, 2, simplify = FALSE)
cat("\nbulk line vs bulk line, specific effects (160 values):", paste(sprintf("%s %.2f", sapply(xl, function(p) paste(sub("NCI-|-211H", "", p), collapse = "-")), sapply(xl, function(p) sp(EB[[p[1]]]$specific, EB[[p[2]]]$specific))), collapse = "  "), "\n")
cat("bulk line vs bulk line, shared response (20 values):  ", paste(sprintf("%s %.2f", sapply(xl, function(p) paste(sub("NCI-|-211H", "", p), collapse = "-")), sapply(xl, function(p) sp(EB[[p[1]]]$shared, EB[[p[2]]]$shared))), collapse = "  "), "\n")
cat("\nvariance of the 160 NTC-referenced effects explained by the shared (module-mean) term:\n")
vs <- function(E) { tot <- sum((E$ntc - mean(E$ntc))^2); 1 - sum(sweep(E$ntc, 1, E$shared)^2) / tot }
print(round(c(screen = vs(ES), sapply(EB, vs)), 3))

## ---------------------------------------------------------------- D module coherence
top50 <- cm_top50(); coh <- list()
for (cl in LINES) { s <- ids[line == cl]; y <- normalise(cs[rowSums(cs[, s] > 5) >= 3, s], "tmm")$y
  g <- factor(grp[line == cl]); y <- y - t(sapply(seq_len(nrow(y)), function(i) ave(y[i, ], g)))       # remove group means: co-variation beyond the KO effect
  mu <- rowMeans(normalise(cs[rownames(y), s], "tmm")$y); bin <- cut(rank(mu), 20, labels = FALSE)
  mc <- function(gs) { if (length(gs) < 5) return(NA); C <- suppressWarnings(cor(t(y[gs, ]), method = "spearman")); mean(C[upper.tri(C)], na.rm = TRUE) }
  for (m in CM) { gs <- intersect(top50[[m]], rownames(y)); obs <- mc(gs)
    nul <- replicate(200, mc(unlist(lapply(split(gs, bin[match(gs, rownames(y))]), function(x) sample(rownames(y)[bin == bin[match(x[1], rownames(y))]], length(x))))))
    coh[[length(coh) + 1]] <- data.frame(line = cl, Cm = m, n_expressed = length(gs), coherence = obs, null_mean = mean(nul, na.rm = TRUE), z = (obs - mean(nul, na.rm = TRUE)) / sd(nul, na.rm = TRUE)) } }
coh <- do.call(rbind, coh); write.csv(coh, file.path(OUT, "spec114_module_coherence.csv"), row.names = FALSE)
cat("\n=== D  module coherence in the bulk lines (top-50 genes; z vs expression-matched random sets) ===\n")
cat("genes of the top 50 expressed, per module x line:\n"); print(tapply(coh$n_expressed, list(coh$Cm, coh$line), mean)[CM, LINES])
cat("\ncoherence z:\n"); print(round(tapply(coh$z, list(coh$Cm, coh$line), mean)[CM, LINES], 1))
cat("\nmodules coherent (z > 3) in all 4 lines:", paste(CM[rowSums(tapply(coh$z, list(coh$Cm, coh$line), mean)[CM, LINES] > 3) == 4], collapse = " "), "\n")

save_fig("fig26_reliability_ceiling", 14, 5.4, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1.1, 1.15)); par(family = "sans", oma = c(0, 0, 3.2, 0.4), mar = c(5.2, 4.6, 3.2, 1))
  k <- scr[grepl("^spec_", scr$what), ]; k$KO <- sub("spec_", "", k$what); k <- k[order(k$n_cells), ]
  plot(k$n_cells, k$reliability, log = "x", pch = 21, bg = KO_COL[k$KO], col = INK2, cex = 2, ylim = c(-0.05, 1), bty = "n", las = 1, xlab = "cells with the guide in the CROP-seq screen (log)", ylab = "reliability of the KO-specific Cm profile (20 modules)")
  text(k$n_cells, k$reliability, k$KO, pos = 3, cex = 0.7, font = 2, col = KO_COL[k$KO], xpd = NA); abline(h = c(0, 0.5), lty = 3, col = GRID)
  mtext("A  screen: reliability depends on cell number", side = 3, adj = 0, line = 1.4, font = 2, cex = 0.82); mtext(sprintf("split-half, %d splits, Spearman-Brown", N_SPLIT), side = 3, adj = 0, line = 0.3, cex = 0.6, col = INK2)
  d <- cei[cei$kind != "shared (module mean over KOs)", ]; m <- rbind(ceiling = d$ceiling, observed = d$observed)
  bp <- barplot(m, beside = TRUE, col = c("#d5d4cf", "#2a78d6"), border = INK2, las = 2, ylim = c(min(-0.2, min(m, na.rm = TRUE) * 1.2), 1), cex.names = 0.6,
                names.arg = paste(sub("NCI-|-211H", "", d$line), ifelse(grepl("specific", d$kind), "specific", "vs NTC"), sep = "\n"), ylab = "Spearman over the 160 KO x Cm effects")
  abline(h = 0, col = INK, lwd = 0.6); legend("topright", c("ceiling = sqrt(rel screen x rel bulk)", "observed screen vs bulk"), fill = c("#d5d4cf", "#2a78d6"), border = NA, bty = "n", cex = 0.68)
  mtext("B  observed agreement vs what reliability allows", side = 3, adj = 0, line = 1.4, font = 2, cex = 0.82); mtext("bulk reliability is an upper bound (replicates share the library block)", side = 3, adj = 0, line = 0.3, cex = 0.6, col = INK2)
  Zc <- tapply(coh$z, list(coh$Cm, coh$line), mean)[CM, LINES]; par(mar = c(5.2, 4.6, 3.2, 3))
  image(seq_along(LINES), seq_along(CM), t(pmin(pmax(Zc, -6), 12))[, length(CM):1], col = colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834", "#7a1f00"))(101), zlim = c(-12, 12), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_along(LINES), labels = sub("NCI-|-211H", "", LINES), tick = FALSE, cex.axis = 0.8); axis(2, at = seq_along(CM), labels = rev(CM), las = 1, tick = FALSE, cex.axis = 0.65)
  for (i in seq_along(CM)) for (j in seq_along(LINES)) text(j, length(CM) - i + 1, sprintf("%.0f", Zc[i, j]), cex = 0.55, col = INK)
  mtext("C  are the Cm programmes coherent in each line?", side = 3, adj = 0, line = 1.4, font = 2, cex = 0.82); mtext("z of mean gene-gene Spearman (KO means removed) vs matched random sets", side = 3, adj = 0, line = 0.3, cex = 0.6, col = INK2)
  mtext("Why agreement is low: each data set's own reliability bounds what can be reproduced", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.98)
  mtext("KO x Cm effects (AddModuleScore, top-50 genes); specific = KO minus the mean of the other 8 groups; vs NTC = KO minus NTC", outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.62, col = INK2)
})
