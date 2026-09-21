# Fig 4: bulk KO (replicate-averaged, RUN1, 3 lines) vs CROP-seq KO similarity, KOs in the same order (matched pair outlined)
# needs: screen_export_figs.R + module-score outputs (see profiles.R)
# outputs: fig4_bulk_vs_CROPseq_corr_{modules,variable_genes}_RUN1.{png,pdf}, fig4_summary.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

P <- load_profiles(THREE)
summ <- list()
for (feat in c("modules", "variable_genes")) {
  D <- if (feat == "modules") P$mod$avg else P$gen$avg
  res <- lapply(c(FALSE, TRUE), function(cen) cross_corr(normalise_profiles(D, cen), THREE))
  vmax <- ceiling(quantile(abs(c(res[[1]]$C, res[[2]]$C)), 0.99, names = FALSE) * 10) / 10
  stt <- lapply(res, function(r) cross_stats(r$C, r$lab))
  ttl <- c("raw (KO − NTC)", "centred across KOs (shared NTC / global component removed)")
  ptitle <- lapply(1:2, function(i) c(ttl[i], sprintf("matched r = %+.2f  vs  other-KO r = %+.2f", stt[[i]]$matched, stt[[i]]$other)))
  main <- if (feat == "modules") "Bulk KO vs CROP-seq KO similarity: 20 Cm modules" else sprintf("Bulk KO vs CROP-seq KO similarity: screen's top variable genes (n=%d)", length(P$genes))
  draw_cross_heatmaps(paste0("fig4_bulk_vs_CROPseq_corr_", feat, "_RUN1"), main,
                      "rows = bulk KO (RUN1) per cell line, columns = CROP-seq KO, same KO order (matched pair outlined: a diagonal = correspondence); MEF2A/MEF2D exist only in the screen.",
                      res, ptitle, lapply(res, function(r) r$lab$ko), THREE, vmax)
  for (i in 1:2) summ[[length(summ) + 1]] <- data.frame(features = feat, variant = c("raw", "centred")[i], matched_r = stt[[i]]$matched, other_r = stt[[i]]$other,
                                                     n_first = stt[[i]]$n_first, n = stt[[i]]$n, median_rank = stt[[i]]$median_rank)
  r <- res[[1]]; st0 <- stt[[1]]                                   # per-line diagnostics (raw)
  for (cl in THREE) { idx <- which(r$lab$line == cl); cat(sprintf("%s %s: mean matched r %.3f, median rank %.1f, n rank1 %d\n", feat, cl, mean(r$C[cbind(idx, match(r$lab$ko[idx], ORDER))]), median(st0$ranks[idx]), sum(st0$ranks[idx] == 1))) }
}
summ <- do.call(rbind, summ); write.csv(summ, file.path(OUT, "fig4_summary.csv"), row.names = FALSE); print(format(summ, digits = 3), row.names = FALSE)
