# Fig 5: as fig4 but each bulk KO is represented ONLY by the replicate best correlated with the matching CROP-seq KO.
# Selection bias is quantified with a control: pick the best replicate separately for every cell.
# needs: screen_export_figs.R + module-score outputs (see profiles.R)
# outputs: fig5_bestreplicate_corr_{modules,variable_genes}_RUN1.{png,pdf} and *_selected_replicates.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

P <- load_profiles(THREE); meta <- P$meta
COLS <- c(ORDER, EXTRA)

run_feature <- function(feat, name, title) {
  avg <- if (feat == "modules") P$mod$avg else P$gen$avg
  rep <- if (feat == "modules") P$mod$rep else P$gen$rep
  S <- avg[["CROP-seq"]]; sS <- S / rms_scale(S); sS[!is.finite(sS)] <- 0
  chosen <- list("CROP-seq" = S); picks <- list(); allrep <- list()
  for (cl in THREE) {
    sc_ <- rms_scale(avg[[cl]]); m <- meta[meta$Cell_Line == cl, ]; cols <- list()
    for (ko in ORDER) {
      samples <- intersect(rownames(m)[m$Gene_Control == ko], names(rep[[cl]])); if (!length(samples)) next
      v <- lapply(setNames(samples, samples), function(s) { x <- rep[[cl]][[s]] / sc_; x[!is.finite(x)] <- 0; x })
      rs <- sapply(v, function(x) cor(x, sS[, ko])); best <- names(which.max(rs))
      cols[[ko]] <- rep[[cl]][[best]]; picks[[paste(cl, ko)]] <- data.frame(line = cl, KO = ko, best_replicate = best, r_to_matched_screen_KO = unname(rs[best]))
      allrep[[paste(cl, ko)]] <- list(ko = ko, v = v)
    }
    chosen[[cl]] <- do.call(cbind, cols)
  }
  res <- lapply(c(FALSE, TRUE), function(cen) cross_corr(normalise_profiles(chosen, cen), THREE))
  # selection-bias control: best replicate chosen separately for every cell (raw)
  ctl <- lapply(allrep, function(a) c(ko = a$ko, sapply(COLS, function(cc) max(sapply(a$v, function(x) cor(x, sS[, cc]))))))
  cm <- sapply(ctl, function(r) as.numeric(r[COLS])); ck <- sapply(ctl, function(r) r[["ko"]])
  ctl_match <- mean(sapply(seq_along(ck), function(i) cm[match(ck[i], COLS), i]))
  ctl_other <- mean(unlist(lapply(seq_along(ck), function(i) cm[setdiff(seq_along(ORDER), match(ck[i], ORDER)), i])))
  vmax <- ceiling(quantile(abs(c(res[[1]]$C, res[[2]]$C)), 0.99, names = FALSE) * 10) / 10
  stt <- lapply(res, function(r) cross_stats(r$C, r$lab))
  ttl <- c("raw (best replicate − NTC)", "centred across KOs (shared NTC / global component removed)")
  ptitle <- lapply(1:2, function(i) c(ttl[i], sprintf("matched r = %+.2f vs other-KO r = %+.2f;  matched ranks 1st in %d/%d, median %.1f of 8", stt[[i]]$matched, stt[[i]]$other, stt[[i]]$n_first, stt[[i]]$n, stt[[i]]$median_rank)))
  rl <- lapply(res, function(r) sprintf("%s  (%s)", r$lab$ko, sapply(seq_len(nrow(r$lab)), function(i) sub(".*_", "", picks[[paste(r$lab$line[i], r$lab$ko[i])]]$best_replicate))))
  draw_cross_heatmaps(name, title,
    c("rows = bulk KO (RUN1), ONLY the replicate best correlated with the matching CROP-seq KO (replicate shown in brackets); columns = CROP-seq KO in the same order (matched pair outlined); MEF2A/MEF2D screen-only.",
      sprintf("selection-bias control: if the best replicate is picked separately for every cell, matched r = %+.2f vs other-KO r = %+.2f (the gap that selection alone produces).", ctl_match, ctl_other)),
    res, ptitle, rl, THREE, vmax)
  pk <- do.call(rbind, picks); rownames(pk) <- NULL; write.csv(pk, file.path(OUT, paste0(name, "_selected_replicates.csv")), row.names = FALSE)
  cat(feat, "raw:", round(c(stt[[1]]$matched, stt[[1]]$other, stt[[1]]$n_first, stt[[1]]$n, stt[[1]]$median_rank), 3),
      "| centred:", round(c(stt[[2]]$matched, stt[[2]]$other, stt[[2]]$n_first, stt[[2]]$n, stt[[2]]$median_rank), 3),
      "| control matched/other:", round(ctl_match, 3), round(ctl_other, 3), "\n")
}
run_feature("modules", "fig5_bestreplicate_corr_modules_RUN1", "Bulk KO (best replicate) vs CROP-seq KO similarity: 20 Cm modules")
run_feature("variable_genes", "fig5_bestreplicate_corr_variable_genes_RUN1", sprintf("Bulk KO (best replicate) vs CROP-seq KO similarity: screen's top variable genes (n=%d)", length(P$genes)))
