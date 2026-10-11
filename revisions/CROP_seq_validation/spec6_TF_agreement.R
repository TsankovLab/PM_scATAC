# Spec 6: which TF knockouts agree across the bulk cell lines (and with the CROP-seq screen)?
# Evidence per TF:
#   (1) on-target knockdown in each line (DESeq2 log2FC of the TF itself, KO vs NTC)
#   (2) cross-line concordance of the KO profile: correlation of the same TF's cc-adjusted DESeq2 Wald statistics between two lines after centring
#       across KOs (removes the response shared by all KOs), and how the same-TF pair ranks among all other-KO pairs (row z-score / rank)
#   (3) consensus genes: genes moving in the same direction (|Wald| > 2, cc-adjusted) in all 3 lines vs the number expected by chance (permutation),
#       and the sign agreement of those genes in the CROP-seq screen (cc-regressed Welch t)
# Combined rank = mean rank of (cross-line z, consensus enrichment, number of lines with effective knockdown). NCI-H28 is reported separately (4-line columns).
# needs: spec1_prepare_stats.R
# outputs: spec_TF_agreement.csv, fig8_TF_agreement.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
set.seed(7)
bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds")); ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); ot <- read.csv(file.path(OUT, "spec_ontarget_knockdown.csv"))
tfs <- sort(Reduce(union, lapply(bs[THREE], function(o) colnames(o$stat0))), method = "radix")
common4 <- Reduce(intersect, lapply(bs[LINES], `[[`, "genes")); common3 <- Reduce(intersect, lapply(bs[THREE], `[[`, "genes"))
cat("TFs:", paste(tfs, collapse = ", "), "| genes common to 3 lines:", length(common3), "| to 4 lines:", length(common4), "\n")

## ---------------------------------------------------------------- (1) on-target knockdown
ko_lfc <- sapply(LINES, function(cl) sapply(tfs, function(k) { i <- ot$cell_line == cl & ot$KO == k; if (any(i)) ot$lfc[i] else NA }))
ko_eff <- sapply(LINES, function(cl) sapply(tfs, function(k) { i <- ot$cell_line == cl & ot$KO == k; any(i) && ot$lfc[i] < -0.5 && ot$padj[i] < 0.05 }))

## ---------------------------------------------------------------- (2) cross-line concordance (cc-adjusted stats, centred across KOs of each line)
X <- setNames(lapply(LINES, function(cl) { m <- bs[[cl]]$statcc[common4, ]; m[!is.finite(m)] <- 0; m - rowMeans(m) }), LINES)
pair_res <- list()
for (a in LINES) for (b in setdiff(LINES, a)) for (k in intersect(colnames(X[[a]]), colnames(X[[b]]))) {
  cols <- colnames(X[[b]]); r <- sapply(cols, function(j) cor(X[[a]][, k], X[[b]][, j]))
  pair_res[[length(pair_res) + 1]] <- data.frame(line_a = a, line_b = b, TF = k, r = r[[k]], rank = sum(r > r[[k]]) + 1, n_cols = length(cols), z = (r[[k]] - mean(r)) / sd(r))
}
pr <- do.call(rbind, pair_res)
pr$pair <- apply(pr[, c("line_a", "line_b")], 1, function(x) paste(sort(x), collapse = " | "))
pair_avg <- aggregate(cbind(r, z, rank, n_cols) ~ pair + TF, pr, mean)                 # average both directions of a pair
pair3 <- pair_avg[pair_avg$pair %in% c("MSTO-211H | NCI-H2052", "MSTO-211H | NCI-H2452", "NCI-H2052 | NCI-H2452"), ]
cross3 <- do.call(rbind, lapply(tfs, function(k) { d <- pair3[pair3$TF == k, ]; data.frame(TF = k, n_pairs3 = nrow(d), cross_r3 = mean(d$r), cross_z3 = mean(d$z), pairs_rank_top2 = sum(d$rank <= 2)) }))
d4 <- pair_avg[grepl("H28", pair_avg$pair), ]
cross4 <- do.call(rbind, lapply(tfs, function(k) { d <- d4[d4$TF == k, ]; data.frame(TF = k, cross_z_withH28 = if (nrow(d)) mean(d$z) else NA, cross_r_withH28 = if (nrow(d)) mean(d$r) else NA) }))

## ---------------------------------------------------------------- (3) consensus genes in all 3 lines
consensus <- function(what, B = 2000) do.call(rbind, lapply(tfs, function(k) {
  if (!all(sapply(THREE, function(cl) k %in% colnames(bs[[cl]][[what]]))))         # KO not present in all 3 lines (RUNX1) -> cannot be assessed
    return(data.frame(TF = k, n_up = NA, n_up_exp = NA, p_up = NA, n_dn = NA, n_dn_exp = NA, p_dn = NA, screen_n_genes = NA, screen_sign_agree = NA, screen_binom_p = NA))
  M <- sapply(THREE, function(cl) { v <- bs[[cl]][[what]][common3, k]; v[!is.finite(v)] <- 0; v }); M[rownames(M) == k, ] <- 0       # exclude the knocked-out gene itself
  N <- nrow(M); one <- function(sgn) { u <- if (sgn > 0) M > 2 else M < -2; n1 <- sum(u[, 1]); obs <- sum(rowSums(u) == 3)
    exp <- N * prod(colMeans(u)); perm <- replicate(B, sum(u[sample.int(N, n1), 2] & u[sample.int(N, n1), 3]))
    # expected number of genes in all three = n1 * P(u2) * P(u3) under independence; permutation p compares observed with random draws of n1 genes in lines 2 and 3
    list(obs = obs, exp = exp, p = (1 + sum(perm >= obs)) / (B + 1), genes = rownames(M)[rowSums(u) == 3]) }
  up <- one(1); dn <- one(-1)
  sgn_scr <- NA; p_scr <- NA; n_scr <- NA
  if (k %in% colnames(ss$tcc)) {
    g_up <- intersect(up$genes, ss$genes); g_dn <- intersect(dn$genes, ss$genes); n_scr <- length(g_up) + length(g_dn)
    if (n_scr >= 10) { agree <- sum(ss$tcc[g_up, k] > 0) + sum(ss$tcc[g_dn, k] < 0); sgn_scr <- agree / n_scr; p_scr <- binom.test(agree, n_scr, 0.5)$p.value } }
  data.frame(TF = k, n_up = up$obs, n_up_exp = up$exp, p_up = up$p, n_dn = dn$obs, n_dn_exp = dn$exp, p_dn = dn$p, screen_n_genes = n_scr, screen_sign_agree = sgn_scr, screen_binom_p = p_scr)
}))
cons <- consensus("statcc"); names(cons)[-1] <- paste0(names(cons)[-1], "_cc"); names(cons) <- sub("^(screen_.*)_cc$", "\\1", names(cons))   # screen columns keep their plain names
cons0 <- consensus("stat0")[, c("TF", "n_up", "n_up_exp", "n_dn", "n_dn_exp")]; names(cons0)[-1] <- paste0(names(cons0)[-1], "_nocc")
cons$enrich_cc <- log2((cons$n_up_cc + cons$n_dn_cc + 1) / (cons$n_up_exp_cc + cons$n_dn_exp_cc + 1))
cons$p_min_cc <- pmin(cons$p_up_cc, cons$p_dn_cc)

## ---------------------------------------------------------------- CROP-seq agreement (F4-like: cc-regressed, centred across KOs), all 8 shared KOs
ctxS <- ss$tcc; common_s <- Reduce(intersect, c(list(ss$genes), lapply(bs[LINES], `[[`, "genes")))
Sx <- ctxS[common_s, ]; Sx <- Sx - rowMeans(Sx)
scr_rows <- lapply(LINES, function(cl) { m <- bs[[cl]]$statcc[common_s, ]; m[!is.finite(m)] <- 0; m <- m - rowMeans(m)
  do.call(rbind, lapply(intersect(ORDER, colnames(m)), function(k) { r <- sapply(ORDER, function(j) cor(m[, k], Sx[, j])); data.frame(line = cl, TF = k, r = r[[k]], z = (r[[k]] - mean(r)) / sd(r), rank = sum(r > r[[k]]) + 1) })) })
scr <- do.call(rbind, scr_rows)
scr_sum <- do.call(rbind, lapply(tfs, function(k) { d <- scr[scr$TF == k & scr$line %in% THREE, ]; data.frame(TF = k, screen_z3 = if (nrow(d)) mean(d$z) else NA, screen_r3 = if (nrow(d)) mean(d$r) else NA, screen_rank_mean = if (nrow(d)) mean(d$rank) else NA, screen_n_cells = if (k %in% names(ss$n)) ss$n[[k]] else NA) }))

## ---------------------------------------------------------------- combine + rank
tab <- Reduce(function(x, y) merge(x, y, by = "TF"), list(data.frame(TF = tfs, n_lines_effective_3 = rowSums(ko_eff[, THREE]), n_lines_effective_H28 = as.integer(ko_eff[, "NCI-H28"]),
        setNames(as.data.frame(round(ko_lfc, 2)), paste0("lfc_", LINES)), check.names = FALSE), cross3, cross4, cons, cons0, scr_sum))
rk <- function(x) rank(-x, ties.method = "average", na.last = "keep")
tab$rank_cross <- rk(tab$cross_z3); tab$rank_consensus <- rk(tab$enrich_cc); tab$rank_knockdown <- rk(tab$n_lines_effective_3)
tab$combined_rank <- rowMeans(tab[, c("rank_cross", "rank_consensus", "rank_knockdown")])
tab$agrees <- tab$cross_z3 > 0.5 & tab$p_min_cc < 0.01 & tab$enrich_cc > 1 & tab$n_lines_effective_3 >= 2
tab <- tab[order(tab$combined_rank, -tab$cross_z3), ]; write.csv(tab, file.path(OUT, "spec_TF_agreement.csv"), row.names = FALSE)
options(width = 230)
show <- tab[, c("TF", "n_lines_effective_3", "cross_r3", "cross_z3", "pairs_rank_top2", "n_up_cc", "n_up_exp_cc", "n_dn_cc", "n_dn_exp_cc", "p_min_cc", "screen_sign_agree", "screen_binom_p", "screen_z3", "combined_rank", "agrees")]
cat("\n=== TF agreement across the 3 lines (sorted by combined rank) ===\n"); print(format(show, digits = 3), row.names = FALSE)
cat("\n=== with NCI-H28 (cross-line z over H28 pairs) ===\n"); print(format(tab[, c("TF", "n_lines_effective_H28", "cross_z_withH28", "cross_r_withH28", "lfc_NCI-H28")], digits = 3), row.names = FALSE)
cat("\nTFs meeting all criteria (cross_z3>0.5, consensus p<0.01 & enrichment>2x, knockdown in >=2 lines):", paste(tab$TF[tab$agrees], collapse = ", "), "\n")

## ---------------------------------------------------------------- figure
save_fig("fig8_TF_agreement", 16, 6.8, function() {
  d <- tab; n <- nrow(d); layout(matrix(1:4, nrow = 1), widths = c(1.2, 1.3, 1.25, 1.5)); par(oma = c(0, 0, 4.2, 0.5), family = "sans")
  panel <- function(main, xlab_n) { par(mar = c(1.2, if (main == 1) 5 else 0.6, 8, 0.4)) }
  # 1 knockdown
  panel(1); M <- as.matrix(d[, paste0("lfc_", LINES)]); M[is.na(M)] <- 0
  image(seq_along(LINES), seq_len(n), t(pmax(pmin(M, 1.5), -1.5))[, n:1], col = DIV, zlim = c(-1.5, 1.5), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_along(LINES)) { txt <- if (is.na(d[i, paste0("lfc_", LINES[j])])) "n/a" else sprintf("%+.1f", d[i, paste0("lfc_", LINES[j])]); text(j, n - i + 1, txt, cex = 0.7, font = ifelse(ko_eff[d$TF[i], LINES[j]], 2, 1)) }
  axis(2, at = n:1, labels = d$TF, las = 1, tick = FALSE, cex.axis = 0.8, line = -0.5); mtext(sub("-211H", "", sub("NCI-", "", LINES)), side = 3, at = seq_along(LINES), las = 2, cex = 0.7, line = 0.3, adj = 0)
  mtext("On-target log2FC (bold = effective)", side = 3, line = 6.3, adj = 0, font = 2, cex = 0.75)
  # 2 cross-line concordance
  panel(2); pr3 <- c("MSTO-211H | NCI-H2052", "MSTO-211H | NCI-H2452", "NCI-H2052 | NCI-H2452", "MSTO-211H | NCI-H28", "NCI-H2052 | NCI-H28", "NCI-H2452 | NCI-H28")
  Z <- sapply(pr3, function(p) sapply(d$TF, function(k) { x <- pair_avg[pair_avg$pair == p & pair_avg$TF == k, ]; if (nrow(x)) x$z else NA }))
  image(seq_along(pr3), seq_len(n), t(pmax(pmin(ifelse(is.na(Z), 0, Z), 2), -2))[, n:1], col = DIV, zlim = c(-2, 2), axes = FALSE, xlab = "", ylab = ""); abline(v = 3.5, col = "white", lwd = 4)
  for (i in seq_len(n)) for (j in seq_along(pr3)) text(j, n - i + 1, ifelse(is.na(Z[i, j]), "n/a", sprintf("%+.1f", Z[i, j])), cex = 0.65)
  mtext(gsub("NCI-|-211H", "", gsub(" \\| ", " / ", pr3)), side = 3, at = seq_along(pr3), las = 2, cex = 0.6, line = 0.3, adj = 0)
  mtext("Same-TF concordance between lines (z vs other KOs)", side = 3, line = 6.3, adj = 0, font = 2, cex = 0.75)
  # 3 consensus genes
  panel(3); z0 <- function(x) ifelse(is.na(x), 0, x); up <- z0(d$n_up_cc); dn <- z0(d$n_dn_cc); upe <- z0(d$n_up_exp_cc); dne <- z0(d$n_dn_exp_cc); mx <- max(c(up, dn, upe, dne, 1))
  plot(NA, xlim = c(-mx * 1.05, mx * 1.05), ylim = c(n + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"); abline(v = 0, col = INK2)
  rect(0, seq_len(n) - 0.32, up, seq_len(n) + 0.32, col = "#eb6834", border = NA); rect(0, seq_len(n) - 0.32, -dn, seq_len(n) + 0.32, col = "#2a78d6", border = NA)
  segments(upe, seq_len(n) - 0.4, upe, seq_len(n) + 0.4, lwd = 2); segments(-dne, seq_len(n) - 0.4, -dne, seq_len(n) + 0.4, lwd = 2)
  text(-mx * 1.05, seq_len(n), ifelse(!is.na(d$p_dn_cc) & d$p_dn_cc < 0.01, "*", ""), cex = 1.3, adj = 0); text(mx * 1.05, seq_len(n), ifelse(!is.na(d$p_up_cc) & d$p_up_cc < 0.01, "*", ""), cex = 1.3, adj = 1)
  text(0, seq_len(n), ifelse(is.na(d$n_up_cc), "n/a (not KO'd in all 3 lines)", ""), cex = 0.55, col = INK2)
  axis(1, at = pretty(c(-mx, mx)), labels = abs(pretty(c(-mx, mx))), cex.axis = 0.6, lwd = 0.5, mgp = c(1, 0.2, 0), line = -1.2)
  mtext("Genes changed in the same direction in all 3 lines", side = 3, line = 6.3, adj = 0, font = 2, cex = 0.75); mtext("down (blue) | up (orange); black tick = expected by chance; * perm p < 0.01", side = 3, line = 5.3, adj = 0, cex = 0.55, col = INK2)
  # 4 screen agreement
  panel(4); S <- d$screen_z3; plot(NA, xlim = c(-1.6, 4.2), ylim = c(n + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"); abline(v = 0, col = INK2)
  rect(0, seq_len(n) - 0.32, ifelse(is.na(S), 0, S), seq_len(n) + 0.32, col = "#8fa88a", border = NA)
  text(2.25, seq_len(n), ifelse(is.na(d$screen_sign_agree), "not in screen", sprintf("%.0f%% of consensus genes agree in sign (p=%.2g)", 100 * d$screen_sign_agree, d$screen_binom_p)), adj = 0, cex = 0.55, col = INK2)
  text(ifelse(is.na(S), 0, pmax(S, 0)) + 0.05, seq_len(n), ifelse(is.na(S), "", sprintf("z=%+.2f", S)), adj = 0, cex = 0.6, col = INK2)
  mtext("CROP-seq screen agreement", side = 3, line = 6.3, adj = 0, font = 2, cex = 0.75); mtext("bar: same-TF match z; text: sign agreement of consensus genes", side = 3, line = 5.3, adj = 0, cex = 0.55, col = INK2)
  mtext("TF knockouts ranked by agreement across NCI-H2052, MSTO-211H and NCI-H2452 (RUN1); NCI-H28 shown in the concordance panel", outer = TRUE, side = 3, adj = 0, line = 2.2, font = 2, cex = 1)
  mtext("rows sorted by combined rank of cross-line concordance, consensus-gene enrichment and number of lines with effective knockdown; all statistics cell-cycle adjusted (DESeq2 ~ cc + condition), concordance after centring across KOs", outer = TRUE, side = 3, adj = 0, line = 1, cex = 0.6, col = INK2)
})
