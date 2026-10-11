# Spec 91: remove the dominant KO axes of each bulk line (top k PCs of the KO-centred NTC - KO profiles, k = 0..3; PC1 = the proliferation / inflammation axis of spec90)
# and test whether matching with CROP-seq and agreement across bulk lines improve.
#   profiles: bulk limma coef NTC - KO (TMM logCPM, S70, no covariate), centred over the 8 KOs of each line, minus the projection on that line's top k PCs;
#             screen: pseudobulk logCPM NTC - KO (spec88), centred over KOs (left as is)
#   genes   : screen top-5000 HVGs expressed in both (spec89) | top 2000 bulk variable genes | all shared genes
#   metrics : screen vs each line -- own-KO Spearman (diag) vs other KOs (off), exact permutation p of the bulk labels (8!), own-KO rank; RSA of the KO-pair similarities
#             across lines -- agreement() (Spearman, KO-label permutation) on the same residual profiles
# needs: spec72 (S70), spec88 (pseudobulk logCPM), spec1 (screen genes)
# outputs: spec91_vs_screen.csv, spec91_vs_screen_per_KO.csv, spec91_across_lines.csv, fig9M_remove_bulk_PCs.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000")); permsL <- make_perms(N_PERM)
BOTH <- Sys.getenv("STRIP_SCREEN", unset = "0") == "1"; OSUF <- if (BOTH) "_both" else ""   # STRIP_SCREEN=1: the top k PCs are removed from the screen profiles too
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); L <- as.matrix(read.csv(file.path(OUT, "spec88_screen_pseudobulk_logCPM.csv"), row.names = 1, check.names = FALSE))
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); Gall <- intersect(pr$genes, ss$genes)
SPACES <- list(hvg = intersect(rownames(L), pr$genes), top2000 = Gall[order(rowMeans(pr$gene_var[Gall, ]), decreasing = TRUE)][1:2000], all = Gall)
SPL <- c(hvg = "screen top-5000 HVGs", top2000 = "top 2000 bulk variable genes", all = "all shared genes")
cen <- function(F) F - rowMeans(F)
scr_prof <- function(f) { if (all(f %in% rownames(L))) cen(SGN * (L[f, "NTC"] - L[f, KOS])) else { S <- -ss$t0[f, KOS]; S[!is.finite(S)] <- 0; cen(S) } }   # pseudobulk FC where available, else Welch t
strip <- function(P, k) { P <- cen(P); if (k == 0) return(P); sv <- svd(P); P - sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], k) %*% t(sv$v[, 1:k, drop = FALSE]) }
perms8 <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
ut <- which(upper.tri(diag(8)))
VS <- list(); VK <- list(); XL <- list()
for (sp in names(SPACES)) { f <- SPACES[[sp]]; S <- scr_prof(f); Ssim <- cor(S, method = "spearman")
  for (k in 0:3) { A <- pr$P$tmm$none$coef[f, , ]; for (cl in LINES) A[, , cl] <- strip(A[, , cl], k)
    if (BOTH) { S <- strip(scr_prof(f), k); Ssim <- cor(S, method = "spearman") }
    for (cl in LINES) { C <- cor(S, A[, , cl], method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms8, 1, function(p) mean(C[cbind(1:8, p)]))
      B <- cor(A[, , cl], method = "spearman"); rs <- cor(Ssim[ut], B[ut], method = "spearman"); rnul <- apply(perms8, 1, function(p) cor(Ssim[ut], B[p, p][ut], method = "spearman"))
      VS[[length(VS) + 1]] <- data.frame(space = sp, k = k, line = cl, diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8), RSA = rs, RSA_p = mean(rnul >= rs - 1e-12))
      VK[[length(VK) + 1]] <- data.frame(space = sp, k = k, line = cl, KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank = sapply(1:8, function(i) rank(-C[i, ])[i])) }
    a <- agreement(cl_from_array(A, NULL, centre = FALSE), permsL); XL[[length(XL) + 1]] <- cbind(space = sp, k = k, a$per_KO, global_p = a$global[["perm_p"]], global_diag = a$global[["diag_mean"]], global_off = a$global[["off_mean"]]) } }
VS <- do.call(rbind, VS); VK <- do.call(rbind, VK); XL <- do.call(rbind, XL)
write.csv(VS, file.path(OUT, sub("\\.csv$", paste0(OSUF, ".csv"), "spec91_vs_screen.csv")), row.names = FALSE); write.csv(VK, file.path(OUT, sub("\\.csv$", paste0(OSUF, ".csv"), "spec91_vs_screen_per_KO.csv")), row.names = FALSE); write.csv(XL, file.path(OUT, sub("\\.csv$", paste0(OSUF, ".csv"), "spec91_across_lines.csv")), row.names = FALSE)
options(width = 220); cat("=== screen vs bulk line: diag - off (exact perm p), RSA ===\n"); print(transform(VS, dmo = diag - off)[, c("space", "k", "line", "diag", "off", "dmo", "perm_p", "top1", "RSA", "RSA_p")], digits = 3, row.names = FALSE)
cat("\n=== H2052, own-KO rho (rank) by k, screen HVGs ===\n"); x <- VK[VK$line == "NCI-H2052" & VK$space == "hvg", ]; print(reshape(transform(x[, c("KO", "k")], v = sprintf("%.2f (%d)", x$own, x$rank)), idvar = "KO", timevar = "k", direction = "wide"), row.names = FALSE)
cat("\n=== across lines: global p and per-KO p by k ===\n"); gx <- unique(XL[, c("space", "k", "global_diag", "global_off", "global_p")]); print(gx, digits = 3, row.names = FALSE)
x <- XL[XL$space == "hvg", ]; print(reshape(transform(x[, c("KO", "k")], v = sprintf("%.3f", x$perm_p)), idvar = "KO", timevar = "k", direction = "wide"), row.names = FALSE)

save_fig(paste0("fig9M_remove_bulk_PCs", OSUF), 15, 9.6, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0.4), mar = c(4.4, 4.4, 2.8, 0.8), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.78)
  for (sp in names(SPACES)) { d <- VS[VS$space == sp, ]; plot(NA, xlim = c(-0.2, 3.2), ylim = range(c(d$diag - d$off, 0)) + c(-0.02, 0.04), xaxt = "n", bty = "n", las = 1, xlab = "top PCs removed from each bulk line", ylab = "screen vs bulk: diag - off (Spearman)", cex.axis = 0.75)
    axis(1, at = 0:3, cex.axis = 0.75); abline(h = 0, lty = 2, col = INK)
    for (cl in LINES) { e <- d[d$line == cl, ]; lines(e$k, e$diag - e$off, col = LINE_COL[[cl]], lwd = 2); points(e$k, e$diag - e$off, pch = ifelse(e$perm_p < 0.05, 16, 1), col = LINE_COL[[cl]], cex = 1.4) }
    if (sp == "hvg") legend("topleft", c(sub("NCI-|-211H", "", LINES), "filled = perm p < 0.05"), col = c(LINE_COL[LINES], INK2), lwd = c(2, 2, 2, 2, NA), pch = c(NA, NA, NA, NA, 16), bty = "n", cex = 0.65)
    ttl(sprintf("A  match with CROP-seq per line -- %s", SPL[[sp]])) }
  x <- VK[VK$line == "NCI-H2052" & VK$space == "hvg", ]; M <- sapply(0:3, function(k) x$own[x$k == k][match(KOS, x$KO[x$k == k])]); rownames(M) <- KOS
  bp <- barplot(t(M), beside = TRUE, col = c("grey60", "#2a78d6", "#1baf7a", "#4a3aa7"), border = NA, las = 2, ylab = "own-KO Spearman, CROP-seq vs bulk H2052", cex.names = 0.75, cex.axis = 0.75); abline(h = 0, col = INK)
  rk <- sapply(0:3, function(k) x$rank[x$k == k][match(KOS, x$KO[x$k == k])]); text(bp, pmax(t(M), 0) + 0.02, ifelse(t(rk) == 1, "1st", ""), cex = 0.5, col = INK2)
  legend("topright", paste("PCs removed:", 0:3), fill = c("grey60", "#2a78d6", "#1baf7a", "#4a3aa7"), border = NA, bty = "n", cex = 0.62); ttl("B  H2052 per TF (screen HVGs); 1st = own KO best of 8")
  d <- VS[VS$space == "hvg", ]; plot(NA, xlim = c(-0.2, 3.2), ylim = range(c(d$RSA, 0)) + c(-0.05, 0.08), xaxt = "n", bty = "n", las = 1, xlab = "top PCs removed", ylab = "RSA: KO-pair similarities, screen vs bulk", cex.axis = 0.75); axis(1, at = 0:3, cex.axis = 0.75); abline(h = 0, lty = 2, col = INK)
  for (cl in LINES) { e <- d[d$line == cl, ]; lines(e$k, e$RSA, col = LINE_COL[[cl]], lwd = 2); points(e$k, e$RSA, pch = ifelse(e$RSA_p < 0.05, 16, 1), col = LINE_COL[[cl]], cex = 1.4) }; ttl("C  KO-similarity structure vs CROP-seq (screen HVGs)")
  gx <- unique(XL[, c("space", "k", "global_diag", "global_off", "global_p")]); plot(NA, xlim = c(-0.2, 3.2), ylim = range(c(gx$global_diag - gx$global_off, 0)) + c(-0.01, 0.02), xaxt = "n", bty = "n", las = 1, xlab = "top PCs removed", ylab = "across bulk lines: same-KO - other-KO rho", cex.axis = 0.75)
  axis(1, at = 0:3, cex.axis = 0.75); abline(h = 0, lty = 2, col = INK); sc <- c(hvg = "#eb6834", top2000 = "#2a78d6", all = "grey40")
  for (sp in names(SPACES)) { e <- gx[gx$space == sp, ]; lines(e$k, e$global_diag - e$global_off, col = sc[[sp]], lwd = 2); points(e$k, e$global_diag - e$global_off, pch = ifelse(e$global_p < 0.05, 16, 1), col = sc[[sp]], cex = 1.4) }
  legend("topright", SPL, col = sc, lwd = 2, bty = "n", cex = 0.62); ttl("D  reproducibility across bulk lines (filled = p < 0.05)")
  mtext("Removing the dominant KO axes (top PCs) of each bulk line before comparing with CROP-seq and across lines (S70, NTC - KO, centred over KOs)", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1)
  mtext(sprintf("PC1 = proliferation (E2F / G2M) vs inflammation (NF-kB / IFN) axis (spec90); A: mean own-KO minus mean other-KO Spearman, exact permutation p of bulk KO labels; screen profiles %s", if (BOTH) "ALSO stripped of their own top k PCs" else "not modified"), outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
})
