# Spec 90: what drives the two-block KO structure of bulk NCI-H2052 (TCF3 / PITX1 / TEAD2 vs TEAD4 / SOX9 / BPTF / HMGA1, spec89)?
#   1 axis      : PC1 of the KO-centred NTC - KO profiles (limma coef, TMM logCPM, S70, the spec89 genes) -> one score per KO (+ % variance); also projected on every sample
#   2 candidates: per sample technical QC (spec71 audit), cell cycle, i5 index group, replicate, the 20 Cm module scores, own-target knockdown;
#                 association with the axis at KO level (group means, n = 8, Spearman) and at sample level (n = 17, Pearson; KO group not modelled -- descriptive)
#   3 adjustment: re-derive the profiles with each candidate as covariate (limma) and measure how much of the block structure is left
#                 (PC1 % variance, mean |rho| between KOs) and the agreement with the screen (RSA, spec89)
#   4 biology   : Hallmark gene sets on the PC1 gene loadings (limma::cameraPR)
# needs: spec72 (S70), spec71 audit, spec88/89 outputs
# outputs: spec90_axis_KO.csv, spec90_covariate_association.csv, spec90_adjustment.csv, spec90_hallmark.csv, fig9L_H2052_axis_covariates.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
CL <- "NCI-H2052"; fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit()
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); L <- as.matrix(read.csv(file.path(OUT, "spec88_screen_pseudobulk_logCPM.csv"), row.names = 1, check.names = FALSE))
G <- intersect(rownames(L), pr$genes); s <- pr$samples[meta[pr$samples, "Cell_Line"] == CL]; g <- meta[s, "Gene_Control"]; y <- normalise(cs[pr$genes, s], "tmm")$y
sheet <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv")); i5 <- sheet$i5_index_num[match(s, sheet$fastq_sample_id)]
cmS <- pr$cm_scores$tmm[, s]; cen <- function(F) F - rowMeans(F)
screenSIM <- cor(cen(SGN * (L[G, "NTC"] - L[G, KOS])), method = "spearman"); ut <- which(upper.tri(screenSIM))
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })

## ---------------------------------------------------------------- 1 axis
P0 <- cen(pr$P$tmm$none$coef[G, , CL]); pc <- prcomp(t(P0), center = FALSE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
ko_score <- pc$x[, 1]; if (ko_score["TCF3"] < 0) { ko_score <- -ko_score; pc$rotation[, 1] <- -pc$rotation[, 1] }
load1 <- pc$rotation[, 1]; samp_score <- as.numeric(crossprod(y[G, ] - rowMeans(y[G, ]), load1)); names(samp_score) <- s
cat(sprintf("PC1 of H2052 KO profiles: %.0f%% of variance (PC2 %.0f%%)\n", ve[1], ve[2])); print(round(sort(ko_score), 1))
write.csv(data.frame(KO = KOS, PC1 = ko_score[KOS], PC2 = pc$x[KOS, 2]), file.path(OUT, "spec90_axis_KO.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 2 candidates
qc_cols <- c("pct_assigned", "assigned_M", "uniq_pct", "pct_unmapped", "pct_nofeature", "genes_detected", "pct_mito", "top100_share", "cc")
cand <- cbind(q[s, qc_cols], log10_depth = log10(q[s, "assigned_M"]), i5_group3 = as.numeric(i5 == 3), rep_R1 = as.numeric(grepl("_R1$", s)), t(cmS)); cand$assigned_M <- NULL
own_t <- sapply(KOS, function(k) if (k %in% rownames(pr$P$tmm$none$t)) pr$P$tmm$none$t[k, k, CL] else NA)
assoc <- do.call(rbind, lapply(colnames(cand), function(v) { x <- cand[[v]]; km <- tapply(x, g, mean)[KOS]
  a <- suppressWarnings(cor.test(km, ko_score[KOS], method = "spearman", exact = FALSE)); b <- suppressWarnings(cor.test(x, samp_score))
  data.frame(covariate = v, type = ifelse(v %in% rownames(cmS), "Cm module", ifelse(v %in% c("i5_group3", "rep_R1"), "batch / design", ifelse(v == "cc", "cell cycle", "technical QC"))),
             KO_level_rho = unname(a$estimate), KO_level_p = a$p.value, sample_level_r = unname(b$estimate), sample_level_p = b$p.value) }))
kd <- cor.test(own_t[KOS], ko_score[KOS], method = "spearman", exact = FALSE)
assoc <- rbind(assoc, data.frame(covariate = "own-target knockdown (t)", type = "knockdown", KO_level_rho = unname(kd$estimate), KO_level_p = kd$p.value, sample_level_r = NA, sample_level_p = NA))
assoc <- assoc[order(-abs(assoc$KO_level_rho)), ]; write.csv(assoc, file.path(OUT, "spec90_covariate_association.csv"), row.names = FALSE)
options(width = 200); cat("\n=== association with the H2052 bulk axis (sorted by |KO-level rho|) ===\n"); print(head(assoc, 20), digits = 3, row.names = FALSE)

## ---------------------------------------------------------------- 3 adjustment
block <- function(P) { C <- cor(cen(P), method = "spearman"); pcx <- prcomp(t(cen(P)), center = FALSE); obs <- cor(screenSIM[ut], C[ut], method = "spearman")
  nul <- apply(perms, 1, function(p) cor(screenSIM[ut], C[p, p][ut], method = "spearman"))
  c(PC1_pct = 100 * pcx$sdev[1]^2 / sum(pcx$sdev^2), mean_abs_rho = mean(abs(C[ut])), corr_with_axis = abs(cor(pcx$x[, 1], ko_score[rownames(pcx$x)])), RSA_screen = obs, RSA_p = mean(nul >= obs - 1e-12)) }
top_tech <- head(setdiff(assoc$covariate[assoc$type != "knockdown" & !assoc$covariate %in% rownames(cmS)], "cc"), 3); top_cm <- head(assoc$covariate[assoc$type == "Cm module"], 3)
adj_sets <- c(list(none = character(0), cc = "cc"), setNames(as.list(top_tech), top_tech), setNames(as.list(top_cm), top_cm), list(cc_pa = c("cc", "pct_assigned")))
adj <- do.call(rbind, lapply(names(adj_sets), function(n) { v <- adj_sets[[n]]; X <- if (length(v)) scale(as.matrix(cand[, v, drop = FALSE]), scale = FALSE) else NULL
  if (!is.null(X) && ncol(y) - 9 - ncol(X) < 2) return(NULL); e <- ko_effects(y[G, ], g, X); data.frame(adjusted_for = n, t(block(e$coef))) }))
write.csv(adj, file.path(OUT, "spec90_adjustment.csv"), row.names = FALSE); cat("\n=== block structure after adjusting the H2052 KO profiles ===\n"); print(adj, digits = 3, row.names = FALSE)

## ---------------------------------------------------------------- 4 biology of the axis
suppressPackageStartupMessages(library(msigdbr)); h <- msigdbr(species = "Homo sapiens", category = "H"); hs <- split(h$gene_symbol, sub("HALLMARK_", "", h$gs_name))
idx <- lapply(hs, function(x) which(names(load1) %in% x)); idx <- idx[sapply(idx, length) >= 10]
hm <- limma::cameraPR(load1, idx); hm$set <- rownames(hm); write.csv(hm, file.path(OUT, "spec90_hallmark.csv"), row.names = FALSE)
cat("\n=== Hallmark on the PC1 gene loadings (Up = higher in the TCF3 / PITX1 / TEAD2 side of NTC - KO, i.e. LOWER after those KOs) ===\n"); print(head(hm[, c("NGenes", "Direction", "PValue", "FDR")], 10), digits = 3)

## ---------------------------------------------------------------- figure
save_fig("fig9L_H2052_axis_covariates", 15, 9.2, function() {
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0.4), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
  par(mar = c(4.6, 4.4, 2.6, 0.8)); o <- names(sort(ko_score)); b <- barplot(ko_score[o], col = KO_COL[o], border = NA, las = 2, ylab = sprintf("PC1 score (%.0f%% of KO-profile variance)", ve[1]), cex.names = 0.75, cex.axis = 0.75); abline(h = 0, col = INK)
  ttl("A  the dominant axis of H2052 bulk KO profiles")
  par(mar = c(4.6, 9, 2.6, 0.8)); a <- head(assoc[!is.na(assoc$KO_level_rho), ], 14); a <- a[nrow(a):1, ]; tc <- c("technical QC" = "#eda100", "cell cycle" = "#1baf7a", "batch / design" = "#4a3aa7", "Cm module" = "#e87ba4", knockdown = "grey50")
  bp <- barplot(a$KO_level_rho, horiz = TRUE, names.arg = a$covariate, las = 1, col = tc[a$type], border = NA, xlim = c(-1, 1), xlab = "Spearman with KO PC1 score (n = 8 KOs)", cex.names = 0.7, cex.axis = 0.75); abline(v = 0, col = INK)
  text(ifelse(a$KO_level_rho > 0, a$KO_level_rho + 0.02, a$KO_level_rho - 0.02), bp, sprintf("p %.2g", a$KO_level_p), adj = ifelse(a$KO_level_rho > 0, 0, 1), cex = 0.55, xpd = NA)
  legend("bottomright", names(tc), fill = tc, border = NA, bty = "n", cex = 0.6); ttl("B  which covariate tracks the axis (KO level)")
  par(mar = c(4.6, 4.4, 2.6, 0.8)); best <- assoc$covariate[which(!is.na(assoc$sample_level_r) & assoc$type != "Cm module")[1]]; x <- cand[[best]]
  plot(x, samp_score, pch = 21, bg = KO_COL[g], col = "white", cex = 1.8, bty = "n", las = 1, xlab = best, ylab = "sample projection on the axis", cex.axis = 0.75); text(x, samp_score, sub("NCI-H2052_|H2052_|H28_", "", s), pos = 3, cex = 0.5, col = INK2, xpd = NA)
  r <- assoc[assoc$covariate == best, ]; ttl(sprintf("C  samples: %s (r %.2f, p %.2g)", best, r$sample_level_r, r$sample_level_p))
  par(mar = c(6.4, 4.4, 2.6, 0.8)); bp <- barplot(adj$PC1_pct, names.arg = adj$adjusted_for, las = 2, col = ifelse(adj$adjusted_for == "none", "grey50", "#2a78d6"), border = NA, ylab = "PC1 % of KO-profile variance", cex.names = 0.7, cex.axis = 0.75)
  abline(h = adj$PC1_pct[adj$adjusted_for == "none"], lty = 2, col = INK2); ttl("D  block strength after adjusting")
  bp <- barplot(adj$RSA_screen, names.arg = adj$adjusted_for, las = 2, col = ifelse(adj$adjusted_for == "none", "grey50", "#1baf7a"), border = NA, ylab = "RSA with the screen (KO-pair similarities)", cex.names = 0.7, cex.axis = 0.75, ylim = range(c(0, adj$RSA_screen)) + c(-0.05, 0.1))
  abline(h = 0, col = INK); text(bp, pmax(adj$RSA_screen, 0) + 0.02, sprintf("p %.2g", adj$RSA_p), cex = 0.55, adj = c(0.5, 0)); ttl("E  similarity structure vs CROP-seq after adjusting")
  par(mar = c(4.6, 12, 2.6, 0.8)); hh <- head(hm, 10); hh <- hh[nrow(hh):1, ]; z <- -log10(hh$PValue) * ifelse(hh$Direction == "Up", 1, -1)
  barplot(z, horiz = TRUE, names.arg = hh$set, las = 1, col = ifelse(z > 0, "#eb6834", "#2a78d6"), border = NA, xlab = "signed -log10 p (cameraPR on PC1 loadings)", cex.names = 0.62, cex.axis = 0.75); abline(v = 0, col = INK)
  ttl("F  biology of the axis (Hallmark)")
  mtext("What drives the two-block KO structure in bulk NCI-H2052? (PC1 of the KO-centred NTC - KO profiles)", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext("A: KO scores; B: covariate group means vs KO score; C: best sample-level covariate; D-E: profiles re-derived with the covariate in the limma model; F: positive = genes on the TCF3 / PITX1 / TEAD2 side", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.6, col = INK2)
})
