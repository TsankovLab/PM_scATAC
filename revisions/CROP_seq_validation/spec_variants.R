# Variant builders shared by spec3 (3 lines) and spec5 (4 lines incl. NCI-H28). Sourced after repro_common.R, spec_common.R, profiles.R.
# A variant is function() -> list(C = rows x ORDER matrix, lab = data.frame(line, ko)) describing how bulk KO profiles compare with screen KO profiles.

load_spec_context <- function(lines) {
  b <- load_bulk(); ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))
  scr_g <- as.matrix(read.csv(file.path(OUT, "screen_var5000_log2avg_by_KO.csv"), row.names = 1, check.names = FALSE))
  varg <- rownames(scr_g)[rownames(scr_g) %in% rownames(b$cp10k)]
  common <- Reduce(intersect, c(list(ss$genes), lapply(bs[lines], `[[`, "genes")))
  list(lines = lines, cp10k = b$cp10k, meta = b$meta, lg = log2(b$cp10k + 1), ss = ss, bs = bs, ot = read.csv(file.path(OUT, "spec_ontarget_knockdown.csv")),
       ccb = read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv")), scr_g = scr_g, varg = varg, F_all = common, F_var = intersect(varg, common))
}

make_row_sets <- function(ctx) {
  ot <- ctx$ot; ss <- ctx$ss
  eff <- ot[ot$cell_line %in% ctx$lines & ot$lfc < -0.5 & ot$padj < 0.05 & ot$KO %in% ORDER, c("cell_line", "KO")]
  scr_down <- colnames(ss$t0)[ss$t0[cbind(match(colnames(ss$t0), rownames(ss$t0)), seq_len(ncol(ss$t0)))] < -2]
  eff2 <- eff[eff$KO %in% scr_down, ]; key <- function(lab) paste(lab$line, lab$ko)
  rs <- list(ALL3 = function(lab) lab$line %in% THREE, H2052 = function(lab) lab$line == SCREEN_LINE)
  if ("NCI-H28" %in% ctx$lines) { rs$ALL4 <- function(lab) rep(TRUE, nrow(lab)); rs$H28 <- function(lab) lab$line == "NCI-H28" }
  rs$EFF <- function(lab) key(lab) %in% paste(eff$cell_line, eff$KO); rs$EFF2 <- function(lab) key(lab) %in% paste(eff2$cell_line, eff2$KO)
  attr(rs, "eff") <- eff; attr(rs, "eff2") <- eff2; rs
}

evaluate_variants <- function(V, row_sets, col_sets) {
  out <- list()
  for (id in names(V)) {
    cm <- V[[id]]$fun()
    for (cs in names(col_sets)) for (rs in names(row_sets)) {
      sel <- row_sets[[rs]](cm$lab) & cm$lab$ko %in% col_sets[[cs]]; if (sum(sel) < 3) next
      out[[length(out) + 1]] <- data.frame(id = id, description = V[[id]]$desc, col_set = cs, row_set = rs,
                                           evaluate_specificity(cm$C[sel, , drop = FALSE], cm$lab$ko[sel], cols = col_sets[[cs]]))
    }
  }
  res <- do.call(rbind, out); rownames(res) <- NULL; res
}

## ---------------------------------------------------------------- core variants (KO - NTC based)
make_core_variants <- function(ctx) {
  with(ctx, {
    V <- list(); add <- function(id, desc, fun) V[[id]] <<- list(desc = desc, fun = fun)
    bulk_diff_log2avg <- function(feats) setNames(lapply(lines, function(cl) { g <- group_log2avg(cp10k, meta, cl)[feats, , drop = FALSE]; g[, colnames(g) != "NTC", drop = FALSE] - g[, "NTC"] }), lines)
    bulk_z_groups <- function(feats) setNames(lapply(lines, function(cl) zscore_by_gene(group_log2avg(cp10k, meta, cl)[feats, , drop = FALSE])), lines)
    bulk_diff_sample_mean <- function(feats, regress_cc) setNames(lapply(lines, function(cl) {
      m <- meta[meta$Cell_Line == cl, ]; y <- lg[feats, rownames(m), drop = FALSE]; cc <- ccb$cc[match(rownames(m), ccb$sample)]
      if (regress_cc) { bcoef <- as.numeric((y - rowMeans(y)) %*% (cc - mean(cc))) / sum((cc - mean(cc))^2); y <- y - outer(bcoef, cc) }
      ntc <- rowMeans(y[, m$Gene_Control == "NTC", drop = FALSE]); kos <- sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix")
      sapply(setNames(kos, kos), function(k) rowMeans(y[, m$Gene_Control == k, drop = FALSE]) - ntc)
    }), lines)
    bulk_stat <- function(feats, what) setNames(lapply(lines, function(cl) { x <- bs[[cl]][[what]][feats, , drop = FALSE]; x[!is.finite(x)] <- 0; x }), lines)
    S_diff_var <- scr_g[varg, colnames(scr_g) != "NTC"] - scr_g[varg, "NTC"]
    add("A1", "KO-NTC log2avg, var genes, RMS-scaled (used so far)", function() cross_matrix(S_diff_var, bulk_diff_log2avg(varg), FALSE, TRUE))
    add("A2", "  centred across KOs", function() cross_matrix(S_diff_var, bulk_diff_log2avg(varg), TRUE, TRUE))
    add("B1", "KO-NTC log2avg, var genes, unscaled", function() cross_matrix(S_diff_var, bulk_diff_log2avg(varg), FALSE, FALSE))
    add("B2", "  centred across KOs", function() cross_matrix(S_diff_var, bulk_diff_log2avg(varg), TRUE, FALSE))
    add("C1", "original CROP-seq: gene z-score across KO groups (incl. NTC)", function() {
      Sz <- zscore_by_gene(scr_g[varg, ]); cross_matrix(Sz[, colnames(Sz) != "NTC"], lapply(bulk_z_groups(varg), function(z) z[, colnames(z) != "NTC", drop = FALSE]), FALSE, FALSE) })
    St0v <- ss$t0[F_var, ]; St0a <- ss$t0[F_all, ]
    add("D1", "t-stat (screen Welch t / bulk DESeq2 Wald), var genes", function() cross_matrix(St0v, bulk_stat(F_var, "stat0"), FALSE, FALSE))
    add("D2", "  centred across KOs", function() cross_matrix(St0v, bulk_stat(F_var, "stat0"), TRUE, FALSE))
    add("D3", "  Spearman", function() cross_matrix(St0v, bulk_stat(F_var, "stat0"), FALSE, FALSE, "spearman"))
    add("D4", "  Spearman, centred", function() cross_matrix(St0v, bulk_stat(F_var, "stat0"), TRUE, FALSE, "spearman"))
    add("E1", "t-stat, all expressed genes", function() cross_matrix(St0a, bulk_stat(F_all, "stat0"), FALSE, FALSE))
    add("E2", "  centred across KOs", function() cross_matrix(St0a, bulk_stat(F_all, "stat0"), TRUE, FALSE))
    add("H0", "KO-NTC mean of log2 CP10K (per sample), var genes, RMS-scaled", function() cross_matrix(ss$diff0[F_var, ], bulk_diff_sample_mean(F_var, FALSE), FALSE, TRUE))
    add("H1", "  + cell-cycle regressed (both datasets)", function() cross_matrix(ss$diffcc[F_var, ], bulk_diff_sample_mean(F_var, TRUE), FALSE, TRUE))
    add("H2", "  H1 centred across KOs", function() cross_matrix(ss$diffcc[F_var, ], bulk_diff_sample_mean(F_var, TRUE), TRUE, TRUE))
    add("F1", "t-stat + cc regressed, var genes", function() cross_matrix(ss$tcc[F_var, ], bulk_stat(F_var, "statcc"), FALSE, FALSE))
    add("F2", "  centred across KOs", function() cross_matrix(ss$tcc[F_var, ], bulk_stat(F_var, "statcc"), TRUE, FALSE))
    add("F3", "t-stat + cc regressed, all expressed genes", function() cross_matrix(ss$tcc[F_all, ], bulk_stat(F_all, "statcc"), FALSE, FALSE))
    add("F4", "  centred across KOs", function() cross_matrix(ss$tcc[F_all, ], bulk_stat(F_all, "statcc"), TRUE, FALSE))
    add("G1", "signature top-100 up/down (screen t), all genes", function() signature_matrix(St0a, bulk_stat(F_all, "stat0"), 100, FALSE))
    add("G2", "  centred across KOs", function() signature_matrix(St0a, bulk_stat(F_all, "stat0"), 100, TRUE))
    add("G3", "signature + cc regressed", function() signature_matrix(ss$tcc[F_all, ], bulk_stat(F_all, "statcc"), 100, FALSE))
    add("G4", "  centred across KOs", function() signature_matrix(ss$tcc[F_all, ], bulk_stat(F_all, "statcc"), 100, TRUE))
    V
  })
}

## ---------------------------------------------------------------- absolute-profile variants: KO groups compared directly, NTC never subtracted
make_abs_variants <- function(ctx) {
  with(ctx, {
    V <- list()
    F_hi <- F_all[order(ss$abs0[F_all, "NTC"], decreasing = TRUE)][1:3000]          # 3000 most highly expressed genes (screen NTC mean)
    bulk_abs <- function(feats, regress_cc) setNames(lapply(lines, function(cl) {     # mean per KO group (incl. NTC) of per-sample log2 CP10K, optionally cc regressed (pooled over the line's samples)
      m <- meta[meta$Cell_Line == cl, ]; y <- lg[feats, rownames(m), drop = FALSE]; cc <- ccb$cc[match(rownames(m), ccb$sample)]
      if (regress_cc) { bcoef <- as.numeric((y - rowMeans(y)) %*% (cc - mean(cc))) / sum((cc - mean(cc))^2); y <- y - outer(bcoef, cc - mean(cc)) }
      gs <- sort(unique(m$Gene_Control), method = "radix"); sapply(setNames(gs, gs), function(g) rowMeans(y[, m$Gene_Control == g, drop = FALSE]))
    }), lines)
    prep_abs <- function(m, mode) {
      k <- m[, colnames(m) != "NTC", drop = FALSE]
      switch(mode, raw = k, koc = k - rowMeans(k), kocz = zscore_by_gene(k), ntcz = zscore_by_gene(m)[, colnames(m) != "NTC", drop = FALSE])
    }
    mode_desc <- c(raw = "no centring", koc = "centred across KOs (no NTC)", kocz = "z-score across KOs (no NTC)", ntcz = "z-score across KOs incl. NTC (original)")
    i <- 0
    for (feat in c("var", "hi")) for (reg in c(FALSE, TRUE)) for (mode in names(mode_desc)) {
      i <- i + 1; local({
        feat <- feat; reg <- reg; mode <- mode
        feats <- if (feat == "var") F_var else F_hi
        V[[sprintf("N%02d", i)]] <<- list(desc = sprintf("absolute KO expr, %s genes%s, %s", if (feat == "var") "variable" else "top-3000 expressed", if (reg) ", cc-regressed" else "", mode_desc[[mode]]),
          fun = function() { S <- prep_abs((if (reg) ss$abscc else ss$abs0)[feats, ], mode); Bp <- lapply(bulk_abs(feats, reg), prep_abs, mode = mode); cross_matrix(S, Bp, FALSE, FALSE) })
      })
    }
    V
  })
}
