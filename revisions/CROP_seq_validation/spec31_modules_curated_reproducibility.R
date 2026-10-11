# Spec 31: same curated sample set as spec30 (65-sample spec28 curation, further NTC-restricted: MSTO keeps only NTC_R2, NCI-H2052 keeps
# only its native non-reassigned NTC_R3), but using the historical DEFAULT fig26 configuration -- 20 Cm modules, cc-regressed t statistic,
# Pearson -- instead of spec30's top3000/plain-t workaround. Two more samples drop out here purely because the precomputed module-score
# and cell-cycle files don't have rows for them (H28_HMGA1_R2, MIXED_H28_TCF3_R1_Huh7 -- the orphaned/mixed samples), leaving 58.
# needs: spec28_curated_KO_selection.R (spec28_curated_KO_selection.csv), cm_sample_scores_AddModuleScore_bulk_RUN1.csv, bulk_cellcycle_per_sample_RUN1.csv
# outputs: spec31_KO_reproducibility_curated58_modules.csv, fig26_KO_reproducibility_boxplots_curated58_modules.{png,pdf},
#          fig24_best_config_heatmap_curated58_modules.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))

fab <- load_full_audited_bulk(); meta_all <- fab$meta
dec <- read.csv(file.path(OUT, "spec28_curated_KO_selection.csv"))
kept_ko <- dec$sample[dec$keep]
ntc_all <- rownames(meta_all)[meta_all$Cell_Line %in% LINES & meta_all$Gene_Control == "NTC"]
ntc_drop <- c("MSTO_NTC_R1", "MSTO_NTC_R3", "H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3")
kept60 <- union(kept_ko, setdiff(ntc_all, ntc_drop))

sc_df <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE); cm_names <- setdiff(colnames(sc_df), c("sample", "cell_line"))
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
no_mod_or_cc <- setdiff(kept60, intersect(sc_df$sample, ccb$sample))
cat("dropping (no precomputed module score / cc index):", paste(no_mod_or_cc, collapse = ", "), "\n")
kept <- setdiff(kept60, no_mod_or_cc)
meta <- meta_all[kept, ]
cat(sprintf("final set: %d samples (60-sample spec30 set minus %d module/cc-missing)\n", nrow(meta), length(no_mod_or_cc)))

make_types <- function(y, cc, grp) {
  ko <- grp != "NTC"; ntc <- grp == "NTC"; nn <- sum(ntc); ccc <- cc - mean(cc)
  tstat <- function(M, with_cc) { g <- factor(grp, levels = c("NTC", sort(setdiff(unique(grp), "NTC")))); X <- if (with_cc) model.matrix(~ cc + g) else model.matrix(~ g); q <- qr(X); res <- qr.resid(q, t(M))
    s2 <- pmax(colSums(res^2) / (ncol(M) - q$rank), 1e-4); Ma <- if (with_cc) M - outer(qr.coef(q, t(M))["cc", ], ccc) else M; (Ma[, ko, drop = FALSE] - rowMeans(Ma[, ntc, drop = FALSE])) / sqrt(s2 * (1 + 1 / nn)) }
  m <- tstat(y, TRUE); m[!is.finite(m)] <- 0; m }
P <- list()
for (cl in LINES) {
  m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m); cc <- ccb$cc[match(s, ccb$sample)]; grp <- as.character(m$Gene_Control)
  y <- t(as.matrix(sc_df[match(s, sc_df$sample), cm_names])) |> `colnames<-`(s)
  M <- make_types(y, cc, grp); ko <- grp[grp != "NTC"]; keepr <- ko %in% ORDER; M <- M[, keepr, drop = FALSE]; ko <- ko[keepr]
  kos <- sort(unique(ko), method = "radix"); G <- sapply(setNames(kos, kos), function(k) rowMeans(M[, ko == k, drop = FALSE])); mu <- rowMeans(G); R <- M - mu; G <- G - mu
  P[[cl]] <- list(R = R, G = G, ko = ko, kos = kos)
  cat(sprintf("  %s: %d samples, KOs: %s\n", cl, length(s), paste(kos, collapse = ",")))
}
kos_all <- sort(unique(unlist(lapply(P, `[[`, "kos"))), method = "radix")
units <- do.call(rbind, lapply(kos_all, function(k) do.call(rbind, lapply(LINES, function(l) if (k %in% P[[l]]$kos) data.frame(KO = k, line = l, id = paste(k, l), stringsAsFactors = FALSE))))); n <- nrow(units)
avgP <- sapply(seq_len(n), function(i) P[[units$line[i]]]$G[, units$KO[i]]); colnames(avgP) <- units$id
pickrep <- vapply(seq_len(n), function(i) { l <- units$line[i]; k <- units$KO[i]; reps <- colnames(P[[l]]$R)[P[[l]]$ko == k]; others <- which(units$KO == k & units$line != l)
  if (!length(others) || length(reps) == 1) return(reps[1]); sc <- sapply(reps, function(r) mean(cor(P[[l]]$R[, r], avgP[, others, drop = FALSE]))); unname(reps[which.max(sc)]) }, character(1))
NPERM <- 5000; set.seed(31)
Rall <- do.call(cbind, lapply(LINES, function(l) P[[l]]$R)); rep_group <- unlist(lapply(LINES, function(l) paste(P[[l]]$ko, l)))
CorRR <- cor(Rall); CorRG <- cor(Rall, avgP); CG <- cor(avgP); same_line <- outer(units$line, units$line, "==")
lines_of <- setNames(lapply(LINES, function(l) which(units$line == l)), LINES)
mats_for <- function(lab) {
  CAl <- CG; CAl[same_line] <- NA; chosen <- integer(n)
  for (u in seq_len(n)) { reps <- which(rep_group == units$id[u]); k <- lab[u]; others <- which(lab == k & units$line != units$line[u])
    chosen[u] <- if (!length(others) || length(reps) == 1) reps[1] else reps[which.max(rowMeans(CorRG[reps, others, drop = FALSE]))] }
  CBl <- CorRR[chosen, chosen]; CBl[same_line] <- NA; list(A = CAl, B = CBl, chosen = chosen) }
ko_stat <- function(C, lab, k) { i <- which(lab == k); o <- which(lab != k); sm <- C[i, i][upper.tri(C[i, i])]; sm <- sm[!is.na(sm)]; ot <- as.vector(C[i, o]); ot <- ot[!is.na(ot)]
  list(same = sm, other = ot, diff = if (length(sm)) median(sm) - median(ot) else NA_real_) }
lab0 <- units$KO; M0 <- mats_for(lab0); kos_stat <- kos_all[kos_all %in% units$KO]
tab <- do.call(rbind, lapply(c("A", "B"), function(m) do.call(rbind, lapply(kos_stat, function(k) { st <- ko_stat(M0[[m]], lab0, k)
  data.frame(version = ifelse(m == "A", "averaged", "best replicate"), KO = k, n_same_pairs = length(st$same), median_same = if (length(st$same)) median(st$same) else NA, median_other = median(st$other),
             p_wilcoxon = if (length(st$same) >= 2) wilcox.test(st$same, st$other, alternative = "greater", exact = FALSE)$p.value else NA, diff_obs = st$diff) }))))
DIST <- lapply(c(A = "A", B = "B"), function(m) setNames(lapply(kos_stat, function(k) ko_stat(M0[[m]], lab0, k)), kos_stat))
perm_diff <- replicate(NPERM, { lab <- lab0; for (l in LINES) { ix <- lines_of[[l]]; lab[ix] <- lab0[ix][sample(length(ix))] }; M <- mats_for(lab)
  c(sapply(kos_stat, function(k) ko_stat(M$A, lab, k)$diff), sapply(kos_stat, function(k) ko_stat(M$B, lab, k)$diff)) })
nA <- length(kos_stat); tab$p_perm <- NA_real_
for (i in seq_len(nrow(tab))) { j <- if (tab$version[i] == "averaged") match(tab$KO[i], kos_stat) else nA + match(tab$KO[i], kos_stat); tab$p_perm[i] <- if (is.na(tab$diff_obs[i])) NA else (1 + sum(perm_diff[j, ] >= tab$diff_obs[i], na.rm = TRUE)) / (NPERM + 1) }
tab$p_perm_BH <- ave(tab$p_perm, tab$version, FUN = function(p) p.adjust(p, "BH")); write.csv(tab, file.path(OUT, "spec31_KO_reproducibility_curated58_modules.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\nconfiguration: modules | t_cc | pearson | curated 58 samples | %d permutations\n", NPERM)); print(format(tab[, c("version", "KO", "n_same_pairs", "median_same", "median_other", "p_wilcoxon", "p_perm", "p_perm_BH")], digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- figure: boxplots, same layout as fig26
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
RLO <- -1; RHI <- 1
save_fig("fig26_KO_reproducibility_boxplots_curated58_modules", 9.5, 5.4, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(3.4, 0, 5.4, 0.3), family = "sans")
  panel <- function(ver, m, ttl) {
    tb <- tab[tab$version == ver, ]; o <- tb[order(ifelse(is.na(tb$p_wilcoxon), Inf, tb$p_wilcoxon)), ]; k <- nrow(o)
    par(mar = c(2.4, 6.6, 3.4, 2.6)); plot(NA, xlim = c(RLO, RHI), ylim = c(k + 0.6, 0.4), axes = FALSE, xlab = "", ylab = "")
    axis(1, at = pretty(c(RLO, RHI)), cex.axis = 0.75, lwd = 0.6, mgp = c(1.4, 0.35, 0)); abline(v = 0, col = GRID, lwd = 0.8)
    bx <- function(v, yc, hh, col) { if (!length(v)) return(invisible()); q <- quantile(v, c(0.25, 0.5, 0.75)); iq <- q[3] - q[1]; lo <- min(v[v >= q[1] - 1.5 * iq]); hi <- max(v[v <= q[3] + 1.5 * iq])
      segments(lo, yc, hi, yc, col = INK2, lwd = 0.9); rect(q[1], yc - hh, q[3], yc + hh, col = adjustcolor(col, 0.9), border = INK2, lwd = 0.9); segments(q[2], yc - hh, q[2], yc + hh, col = INK, lwd = 2)
      points(v, yc + runif(length(v), -hh * 0.7, hh * 0.7), pch = 16, cex = if (length(v) > 20) 0.35 else 0.6, col = adjustcolor(INK, 0.45)) }
    for (i in seq_len(k)) { yc <- i; st <- DIST[[m]][[o$KO[i]]]
      bx(st$other, yc + 0.24, 0.17, GREY); bx(st$same, yc - 0.24, 0.17, KO_COL[o$KO[i]])
      text(RHI + (RHI - RLO) * 0.04, yc, stars(o$p_wilcoxon[i]), adj = 0, cex = 0.85, xpd = NA, font = 2) }
    for (i in seq_len(k)) text(RLO - (RHI - RLO) * 0.03, i, o$KO[i], adj = 1, xpd = NA, cex = 0.78, font = 2, col = KO_COL[o$KO[i]])
    mtext(ttl, side = 3, adj = 0, line = 2.3, font = 2, cex = 0.85)
    mtext(sprintf("ordered by Wilcoxon p (most to least significant); n=%d line pairs", o$n_same_pairs[1]), side = 3, adj = 0, line = 1.1, cex = 0.55, col = INK2)
  }
  panel("averaged", "A", "A  Replicates averaged"); panel("best replicate", "B", "B  Best replicate per KO and line (selection-biased)")
  legend("bottomright", inset = c(0.02, -0.34), xpd = NA, horiz = TRUE, legend = c("same KO, other line", "other KOs"), fill = c(adjustcolor(INK, 0.9), GREY), border = INK2, bty = "n", cex = 0.65)
  mtext("How well does each bulk KO reproduce across cell lines? (curated 58-sample set)", outer = TRUE, side = 3, adj = 0, line = 2.9, font = 2, cex = 1)
  mtext("20 Cm modules, cc-regressed t statistic, centred on each line; the 8 CROP-seq KOs; 58 curated samples (spec28+spec30 curation, minus 2 module/cc-score-missing samples)",
        outer = TRUE, side = 3, adj = 0, line = 1.9, cex = 0.6, col = INK2)
  mtext("boxplots: cross-line Pearson r of each KO with itself (colour) vs with the other KOs (grey), median/IQR/range, dots = individual line pairs", outer = TRUE, side = 3, adj = 0, line = 1.0, cex = 0.55, col = INK2)
  mtext("*** p<0.001, ** p<0.01, * p<0.05, ns = not significant (one-sided Wilcoxon test on the pooled cross-line correlations; pairs share samples so this is optimistic, esp. for B)", outer = TRUE, side = 3, adj = 0, line = 0.1, cex = 0.55, col = INK2)
})

## ---------------------------------------------------------------- figure: heatmap, same layout as fig24
kcol <- KO_COL; short <- function(l) sub("NCI-", "", sub("-211H", "", l)); fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
CA <- M0$A; CB <- M0$B; pick_lab <- sub(".*_", "", colnames(Rall)[M0$chosen])
ordv <- function(ver) { tb <- tab[tab$version == ver, ]; o <- tb$KO[order(ifelse(is.na(tb$p_wilcoxon), Inf, tb$p_wilcoxon), tb$p_perm)]; ou <- unlist(lapply(o, function(k) which(units$KO == k)))
  list(ordk = o, ordu = ou, U = units[ou, ], st = cumsum(c(1, table(factor(units$KO[ou], o))))[seq_along(o)], en = cumsum(table(factor(units$KO[ou], o)))) }
OA <- ordv("averaged"); OB <- ordv("best replicate"); CAo <- CA[OA$ordu, OA$ordu]; CBo <- CB[OB$ordu, OB$ordu]
lab_A <- paste(OA$U$KO, short(OA$U$line)); lab_B <- paste0(OB$U$KO, " ", short(OB$U$line), " (", pick_lab[OB$ordu], ")")
save_fig("fig24_best_config_heatmap_curated58_modules", 17, 5.4, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(3.2, 0, 3.4, 0.3), family = "sans"); VM <- 0.8
  panel <- function(C, ver, m, ttl, labs, O) {
    ordk <- O$ordk; starts <- O$st; ends <- O$en
    par(mar = c(2.4, 6.3, 3.4, 20.5)); vmax <- VM; image(seq_len(n), seq_len(n), heat_z(pmax(pmin(C, vmax), -vmax)), col = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100), zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
    na <- which(is.na(C), arr.ind = TRUE); rect(na[, 2] - 0.5, n - na[, 1] + 0.5, na[, 2] + 0.5, n - na[, 1] + 1.5, col = "#efefef", border = NA)
    axis(2, at = n:1, labels = labs, las = 1, tick = FALSE, cex.axis = 0.5, line = -0.6)
    unit_in <- par("pin")[1] / n; xin <- function(inch) n + 0.5 + inch / unit_in; xr <- function(r) xin(0.95 + (pmax(pmin(r, RHI), RLO) - RLO) / (RHI - RLO) * 1.75); tb <- tab[tab$version == ver, ]
    segments(xr(0), 0.5, xr(0), n + 0.5, col = GRID, lwd = 0.8, xpd = NA)
    bx <- function(v, yc, hh, col) { if (!length(v)) return(invisible()); q <- quantile(v, c(0.25, 0.5, 0.75)); iq <- q[3] - q[1]; lo <- min(v[v >= q[1] - 1.5 * iq]); hi <- max(v[v <= q[3] + 1.5 * iq])
      segments(xr(lo), yc, xr(hi), yc, col = INK2, lwd = 0.7, xpd = NA); rect(xr(q[1]), yc - hh, xr(q[3]), yc + hh, col = adjustcolor(col, 0.9), border = INK2, lwd = 0.7, xpd = NA); segments(xr(q[2]), yc - hh, xr(q[2]), yc + hh, col = INK, lwd = 1.6, xpd = NA)
      points(xr(v), yc + runif(length(v), -hh * 0.75, hh * 0.75), pch = 16, cex = if (length(v) > 20) 0.16 else 0.32, col = adjustcolor(INK, 0.45), xpd = NA) }
    for (i in seq_along(ordk)) { ytop <- n - starts[i] + 1.5; ybot <- n - ends[i] + 0.5; h <- ytop - ybot; k <- ordk[i]; r <- tb[tb$KO == k, ]; st <- DIST[[m]][[k]]
      rect(starts[i] - 0.5, ybot, ends[i] + 0.5, ytop, border = INK, lwd = 1.0); rect(xin(0.02), ybot + 0.05, xin(0.09), ytop - 0.05, col = kcol[k], border = NA, xpd = NA); text(xin(0.13), (ytop + ybot) / 2, k, adj = 0, cex = 0.62, font = 2, xpd = NA)
      bx(st$same, ytop - h * 0.27, h * 0.15, kcol[k]); bx(st$other, ybot + h * 0.27, h * 0.15, GREY)
      text(xin(3.4), (ytop + ybot) / 2, fp(r$p_wilcoxon), adj = 1, cex = 0.6, font = ifelse(!is.na(r$p_wilcoxon) & r$p_wilcoxon < 0.05, 2, 1), xpd = NA); text(xin(3.98), (ytop + ybot) / 2, fp(r$p_perm), adj = 1, cex = 0.6, font = ifelse(!is.na(r$p_perm) & r$p_perm < 0.05, 2, 1), xpd = NA) }
    tk <- c(-0.8, 0, 0.8); segments(xr(tk), 0.5, xr(tk), 0.1, col = INK2, lwd = 0.6, xpd = NA); text(xr(tk), -0.9, tk, cex = 0.5, col = INK2, xpd = NA); text(xr(0.1), -2.3, "cross-line r (colour = same KO, grey = other KOs)", cex = 0.5, col = INK2, xpd = NA)
    text(xin(c(3.4, 3.98)), n + 1.6, c("p\nWilcoxon", "p\npermutation"), adj = 1, cex = 0.5, col = INK2, xpd = NA); text(xr(0.05), n + 1.6, "median and IQR of the\ncorrelations between lines", cex = 0.5, col = INK2, xpd = NA)
    mtext(ttl, side = 3, adj = 0, line = 1.7, font = 2, cex = 0.85) }
  panel(CAo, "averaged", "A", "A  Replicates averaged", lab_A, OA); panel(CBo, "best replicate", "B", "B  Best replicate per KO and line (selection-biased)", lab_B, OB)
  draw_colorbar(c(0.38, 0.62), c(-VM, VM), "Pearson r between the KO profiles of two cell lines", cex = 0.55)
  mtext("Cross-cell-line reproducibility of bulk KO profiles (curated 58-sample set, 20 Cm modules, cc-regressed t statistic)", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext(sprintf("20 modules; heatmap boxes = same KO in different lines, grey cells = same line; boxplots: cross-line correlations of the KO with itself (colour) vs with the other KOs (grey); KOs ordered by the Wilcoxon p; one-sided Wilcoxon p and label-permutation p (%d shuffles)", NPERM), outer = TRUE, side = 3, adj = 0, line = 0.7, cex = 0.5, col = INK2)
})
