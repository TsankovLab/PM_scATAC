# Spec 75: which KO effects on the 20 Cm modules are preserved between the CROP-seq screen and bulk RNA-seq, and across the bulk cell lines? Every KO x Cm pair (8 x 20) gets the same test.
# Direction NTC - KO (positive = module lower after KO); DIFF_DIR=KO_NTC for the reverse.
#   bulk, per line : (a) module score (AddModuleScore, top-50 genes, TMM logCPM) ~ group + cc + % assigned, limma moderated t of NTC - KO over the 20 modules;
#                    (b) gene level: limma-voom with sample quality weights, same design, then limma::camera on the top-50 genes of the module (competitive gene-set test that
#                        accounts for inter-gene correlation) for the contrast NTC - KO.
#   bulk, combined : weighted Stouffer z over the lines (weight = sqrt(1 / (1/n_NTC + 1/n_KO))), two-sided p, BH over the 160 pairs; sign consistency over the lines;
#                    Cochran Q heterogeneity p of the module-score effects. Done with all 4 lines and without NCI-H2052 (so the screen line is not needed for the bulk claim).
#   screen         : AddModuleScore on cells (top-50 genes, as spec45), NTC vs KO cells, Wilcoxon, median difference; raw and with the cell-cycle index regressed out of the score; BH over 160.
#   specificity    : rank of the KO among the 8 KOs for that module (1 = strongest effect in the direction of the screen effect), in the screen and in bulk combined.
#   KO-specific    : same tests against the OTHER KOs instead of NTC (bulk: mean(other KOs) - KO per line, Stouffer; screen: cells of the other KOs vs KO cells). NTC is not used, so a
#                    shared-NTC artefact (spec74 C) or an effect common to every KO cannot produce it. preserved_specific = PRESERVED and significant with the same sign in both.
# PRESERVED = screen FDR < 0.05 AND bulk combined FDR < 0.1 (module score) AND same sign AND same sign in >= 3 of the 4 lines AND rank <= 3 in both.
# env: SET (S70 default | S70N4 | S75), CC (1 default = cell cycle adjusted on both sides, 0 = not)
# needs: spec71_sample_audit_76.R, spec72_build_profiles.R (module scores), git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen())
# outputs: spec75_Cm_KO_table_<SET>.csv, spec75_Cm_KO_per_line_<SET>.csv, spec75_preserved_<SET>.csv, spec75_screen_scores_by_cell.rds
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "agree_common.R"))
SET <- Sys.getenv("SET", unset = "S70"); CC <- Sys.getenv("CC", unset = "1") == "1"; SUF <- paste0(SET, if (CC) "" else "_noCC", DSUF)
top50 <- cm_top50(); CM <- names(top50); q <- load_audit(); fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_%s%s.rds", SET, DSUF))); ids <- pr$samples; sc <- pr$cm_scores$tmm

## ---------------------------------------------------------------- bulk per line
per_line <- do.call(rbind, lapply(LINES, function(cl) { s <- ids[meta[ids, "Cell_Line"] == cl]; kos_l <- KOS[KOS %in% meta[s, "Gene_Control"]]   # KOs present in this line (RUN3 drops some)
  g <- factor(meta[s, "Gene_Control"], levels = c("NTC", kos_l)); G <- model.matrix(~ 0 + g); colnames(G) <- levels(g)
  X <- covar_matrix(q[s, ], if (CC) "cc_pa" else "none"); if (!CC) X <- scale(cbind(pa = q[s, "pct_assigned"]), scale = FALSE); des <- cbind(G, X)
  d <- calcNormFactors(DGEList(cs[pr$genes, s])); v <- voomWithQualityWeights(d, des, plot = FALSE); fg <- lmFit(v, des)
  fm <- lmFit(sc[, s], des, weights = matrix(v$targets$sample.weights, nrow(sc), length(s), byrow = TRUE)); idx <- lapply(top50, function(x) which(rownames(v) %in% x))
  do.call(rbind, lapply(kos_l, function(k) { L <- setNames(numeric(ncol(des)), colnames(des)); L["NTC"] <- SGN; L[k] <- -SGN
    f <- eBayes(contrasts.fit(fm, L)); cam <- camera(v, idx, des, contrast = L)[CM, ]
    L2 <- setNames(numeric(ncol(des)), colnames(des)); L2[setdiff(kos_l, k)] <- SGN / (length(kos_l) - 1); L2[k] <- -SGN; f2 <- eBayes(contrasts.fit(fm, L2))      # mean(other KOs) - KO: NTC not used
    data.frame(line = cl, KO = k, Cm = CM, n_NTC = sum(g == "NTC"), n_KO = sum(g == k), effect = f$coefficients[, 1], se = f$stdev.unscaled[, 1] * sqrt(f$s2.post), t = f$t[, 1], p = f$p.value[, 1],
               camera_dir = ifelse(cam$Direction == "Up", 1, -1), camera_p = cam$PValue, n_genes = cam$NGenes, rest_effect = f2$coefficients[, 1], rest_t = f2$t[, 1], rest_p = f2$p.value[, 1], row.names = NULL) })) }))
per_line$rest_z <- sign(per_line$rest_t) * qnorm(per_line$rest_p / 2, lower.tail = FALSE)
per_line$z <- sign(per_line$t) * qnorm(per_line$p / 2, lower.tail = FALSE); per_line$camera_z <- per_line$camera_dir * qnorm(per_line$camera_p / 2, lower.tail = FALSE)
per_line$w <- sqrt(1 / (1 / per_line$n_NTC + 1 / per_line$n_KO)); write.csv(per_line, file.path(OUT, sprintf("spec75_Cm_KO_per_line_%s.csv", SUF)), row.names = FALSE)

combine <- function(d, tag) { r <- do.call(rbind, lapply(split(d, list(d$KO, d$Cm), drop = TRUE), function(x) { st <- function(z) sum(x$w * z) / sqrt(sum(x$w^2)); zs <- st(x$z); zc <- st(x$camera_z); zr <- sum(x$rest_z) / sqrt(nrow(x))
    wi <- 1 / x$se^2; eb <- sum(wi * x$effect) / sum(wi); Q <- sum(wi * (x$effect - eb)^2)
    data.frame(KO = x$KO[1], Cm = x$Cm[1], effect = eb, z = zs, p = 2 * pnorm(-abs(zs)), camera_z = zc, camera_p = 2 * pnorm(-abs(zc)), rest_z = zr, rest_p = 2 * pnorm(-abs(zr)), rest_n_pos = sum(x$rest_effect > 0), n_lines = nrow(x), n_pos = sum(x$effect > 0), n_neg = sum(x$effect < 0),
               het_p = if (nrow(x) > 1) pchisq(Q, nrow(x) - 1, lower.tail = FALSE) else NA) }))
  r$fdr <- p.adjust(r$p, "BH"); r$camera_fdr <- p.adjust(r$camera_p, "BH"); r$rest_fdr <- p.adjust(r$rest_p, "BH"); names(r)[-(1:2)] <- paste0(tag, "_", names(r)[-(1:2)]); r }
b4 <- combine(per_line, "bulk"); b3 <- combine(per_line[per_line$line != "NCI-H2052", ], "bulk3")
wide <- reshape(per_line[, c("KO", "Cm", "line", "effect", "p")], idvar = c("KO", "Cm"), timevar = "line", direction = "wide")

## ---------------------------------------------------------------- screen
suppressPackageStartupMessages(library(Seurat)); f_sc <- file.path(OUT, "spec75_screen_scores_by_cell.rds")
if (!file.exists(f_sc)) { srt <- load_screen(); srt <- AddModuleScore(srt, features = top50, name = "Cm_", seed = 1)
  S <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(S) <- CM; saveRDS(list(score = S, group = srt$merged_call, cc = srt$cc), f_sc) }
scr <- readRDS(f_sc); S <- scr$score; if (CC) S <- apply(S, 2, function(y) resid(lm(y ~ scr$cc)) + mean(y))
screen <- do.call(rbind, lapply(KOS, function(k) do.call(rbind, lapply(CM, function(m) { a <- S[scr$group == "NTC", m]; b <- S[scr$group == k, m]
  o <- S[scr$group %in% setdiff(KOS, k), m]
  data.frame(KO = k, Cm = m, screen_n_KO = length(b), screen_effect = SGN * (median(a) - median(b)), screen_mean_effect = SGN * (mean(a) - mean(b)), screen_p = wilcox.test(a, b)$p.value,
             screen_rest_effect = SGN * (median(o) - median(b)), screen_rest_p = wilcox.test(o, b)$p.value) }))))
screen$screen_fdr <- p.adjust(screen$screen_p, "BH"); screen$screen_rest_fdr <- p.adjust(screen$screen_rest_p, "BH")

## ---------------------------------------------------------------- merge, specificity ranks, preserved pairs
tab <- Reduce(function(a, b) merge(a, b, by = c("KO", "Cm")), list(screen, b4, b3, wide))
rk <- function(v, s) ave(-v * s, tab$Cm, FUN = rank)                     # 1 = strongest effect in the direction s (the screen sign of this pair)
s <- sign(tab$screen_effect); tab$screen_rank <- rk(tab$screen_effect, s); tab$bulk_rank <- rk(tab$bulk_z, s); tab$bulk3_rank <- rk(tab$bulk3_z, s)
tab$same_sign <- sign(tab$bulk_effect) == s; tab$n_lines_same_sign <- ifelse(s > 0, tab$bulk_n_pos, tab$bulk_n_neg)
tab$preserved <- tab$screen_fdr < 0.05 & tab$bulk_fdr < 0.1 & tab$same_sign & tab$n_lines_same_sign >= 3 & tab$screen_rank <= 3 & tab$bulk_rank <= 3
# KO-SPECIFIC: the effect also holds against the other KOs (no NTC involved on either side), so it is neither a shared-NTC artefact nor an effect common to every KO
tab$preserved_specific <- tab$preserved & sign(tab$bulk_rest_z) == s & tab$bulk_rest_p < 0.05 & sign(tab$screen_rest_effect) == s & tab$screen_rest_fdr < 0.05
tab$preserved_camera <- tab$screen_fdr < 0.05 & tab$bulk_camera_fdr < 0.1 & sign(tab$bulk_camera_z) == s
tab$preserved_no_H2052 <- tab$screen_fdr < 0.05 & tab$bulk3_fdr < 0.1 & sign(tab$bulk3_effect) == s & ifelse(s > 0, tab$bulk3_n_pos, tab$bulk3_n_neg) >= 2
tab$n_criteria <- (tab$screen_fdr < 0.05) + (tab$bulk_fdr < 0.1) + tab$same_sign + (tab$n_lines_same_sign >= 3) + (tab$screen_rank <= 3) + (tab$bulk_rank <= 3)
tab <- tab[order(-tab$preserved_specific, -tab$preserved, -tab$n_criteria, tab$bulk_p), ]; write.csv(tab, file.path(OUT, sprintf("spec75_Cm_KO_table_%s.csv", SUF)), row.names = FALSE)
write.csv(tab[tab$preserved | tab$preserved_camera & tab$same_sign, ], file.path(OUT, sprintf("spec75_preserved_%s.csv", SUF)), row.names = FALSE)

options(width = 260); cols <- c("KO", "Cm", "screen_effect", "screen_fdr", "screen_rank", "bulk_effect", "bulk_z", "bulk_p", "bulk_fdr", "bulk_camera_z", "bulk_camera_fdr", "bulk_rank", "n_lines_same_sign", "bulk_rest_z", "bulk_rest_p", "screen_rest_fdr", "bulk3_z", "bulk3_fdr", "preserved", "preserved_specific", "preserved_camera", "preserved_no_H2052")
cat(sprintf("\n=== %s, %s, cell cycle adjusted: %s -- %d of 160 KO x Cm pairs PRESERVED (all criteria) ===\n", SET, DLAB, CC, sum(tab$preserved))); print(tab[tab$preserved, cols], digits = 3, row.names = FALSE)
cat("\n=== next best (5 of 6 criteria) ===\n"); print(head(tab[!tab$preserved & tab$n_criteria >= 5, cols], 25), digits = 3, row.names = FALSE)
cat("\n=== TEAD4 / Cm20 ===\n"); print(tab[tab$KO == "TEAD4" & tab$Cm == "Cm20", ], digits = 3, row.names = FALSE)
cat("\n=== agreement of the whole 8 x 20 matrix: screen effect vs bulk combined z ===\n")
cat(sprintf("all pairs: Spearman rho = %.3f (p = %.3g); pairs with screen FDR < 0.05 (n = %d): sign agreement = %.0f%%, binomial p = %.3g\n", cor(tab$screen_effect, tab$bulk_z, method = "spearman"),
  cor.test(tab$screen_effect, tab$bulk_z, method = "spearman", exact = FALSE)$p.value, sum(tab$screen_fdr < 0.05), 100 * mean(tab$same_sign[tab$screen_fdr < 0.05]),
  binom.test(sum(tab$same_sign[tab$screen_fdr < 0.05]), sum(tab$screen_fdr < 0.05), 0.5, alternative = "greater")$p.value))
print(do.call(rbind, lapply(KOS, function(k) { d <- tab[tab$KO == k, ]; ct <- cor.test(d$screen_effect, d$bulk_z, method = "spearman", exact = FALSE); data.frame(KO = k, rho_screen_vs_bulk_over_20Cm = ct$estimate, p = ct$p.value, n_bulk_fdr10 = sum(d$bulk_fdr < 0.1), n_preserved = sum(d$preserved)) })), digits = 3, row.names = FALSE)
