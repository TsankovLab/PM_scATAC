# Fig 3: KO x KO similarity (Pearson r between KO - NTC profiles) across CROP-seq and the four bulk lines (RUN1),
#        using the 20 Cm modules or the screen's top variable genes; raw and "centred" (per-dataset feature mean over KOs removed)
# needs: screen_export_figs.R + the module-score outputs (see profiles.R)
# outputs: fig3_corr_{modules,variable_genes}_{raw,centred}.{png,pdf,csv}, fig3_matched_KO_rank.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

P <- load_profiles(LINES)
DS <- c("CROP-seq", LINES)

stack_profiles <- function(N) {
  X <- list(); lab <- list()
  for (ds in DS) for (ko in sort(colnames(N[[ds]]), method = "radix")) { X[[length(X) + 1]] <- N[[ds]][, ko]; lab[[length(lab) + 1]] <- c(ds, ko) }
  list(X = do.call(rbind, X), lab = setNames(as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE), c("ds", "ko")))
}

draw_matrix <- function(C, lab, title, sub, name) {
  n <- nrow(C)
  save_fig(name, 13.2, 12.2, function() {
    par(oma = c(4, 0, 0, 0), mar = c(5, 10, 8, 6), family = "sans")
    image(seq_len(n), seq_len(n), heat_z(C), col = DIV, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
    axis(1, at = seq_len(n), labels = lab$ko, las = 2, tick = FALSE, cex.axis = 0.5, line = -0.9)
    axis(2, at = n:1, labels = lab$ko, las = 1, tick = FALSE, cex.axis = 0.5, line = -0.5)
    edges <- which(lab$ds[-1] != lab$ds[-n]); abline(h = n - edges + 0.5, v = edges + 0.5, col = "white", lwd = 3)
    bounds <- c(0, edges, n); par(xpd = NA)
    for (b in seq_along(DS)) {
      a <- bounds[b] + 0.5; z <- bounds[b + 1] + 0.5; mid <- (a + z) / 2
      rect(a, n + 1.3, z, n + 2.6, col = LINE_COL[[DS[b]]], border = NA); text(mid, n + 3.4, DS[b], font = 2, cex = 0.85, col = LINE_COL[[DS[b]]])
      rect(-6.4, n - z + 1, -5.1, n - a + 1, col = LINE_COL[[DS[b]]], border = NA); text(-7.4, n - mid + 1, DS[b], srt = 90, font = 2, cex = 0.85, col = LINE_COL[[DS[b]]])
    }
    mtext(title, side = 3, adj = 0, line = 6.2, font = 2, cex = 1.2); mtext(sub, side = 3, adj = 0, line = 4.9, cex = 0.62, col = INK2)
    draw_colorbar(c(0.35, 0.65), c(-1, 1), "Pearson r between KO profiles")
  })
}

matched_rank <- function(C, lab, tag) {
  scr_i <- which(lab$ds == "CROP-seq"); out <- list()
  for (i in which(lab$ds != "CROP-seq")) {
    j <- match(lab$ko[i], lab$ko[scr_i]); if (is.na(j)) next
    r <- C[i, scr_i]; same <- r[j]
    out[[length(out) + 1]] <- data.frame(features = tag, cell_line = lab$ds[i], KO = lab$ko[i], r_same_KO = same, rank_of_same_KO = sum(r > same) + 1,
                                        n_screen_KOs = length(scr_i), mean_r_other_KOs = mean(r[-j]))
  }
  do.call(rbind, out)
}

res <- list()
for (feat in c("modules", "variable_genes")) {
  D <- if (feat == "modules") P$mod$avg else P$gen$avg
  for (centre in c(FALSE, TRUE)) {
    sp <- stack_profiles(normalise_profiles(D, centre)); C <- cor(t(sp$X)); dimnames(C) <- list(paste(sp$lab$ds, sp$lab$ko, sep = "|"), paste(sp$lab$ds, sp$lab$ko, sep = "|"))
    tag <- paste0(feat, "_", if (centre) "centred" else "raw")
    write.csv(C, file.path(OUT, paste0("fig3_corr_", tag, ".csv")))
    ttl <- if (feat == "modules") "KO–KO similarity from the 20 Cm modules" else sprintf("KO–KO similarity from the screen's top variable genes (n=%d)", length(P$genes))
    sub <- paste0("profiles = KO − NTC per dataset; each feature scaled to equal RMS within dataset",
                  if (centre) "; features centred across KOs within each dataset first (removes the shared NTC/global component)" else " (NTC noise is shared by every KO of a line)")
    draw_matrix(C, sp$lab, paste0(ttl, if (centre) "  [centred]" else ""), sub, paste0("fig3_corr_", tag))
    res[[tag]] <- matched_rank(C, sp$lab, tag)
  }
}
res <- do.call(rbind, res); rownames(res) <- NULL
write.csv(res, file.path(OUT, "fig3_matched_KO_rank.csv"), row.names = FALSE)
agg <- do.call(rbind, lapply(split(res, res$features), function(d) data.frame(features = d$features[1], n = nrow(d), top1 = sum(d$rank_of_same_KO == 1), top3 = sum(d$rank_of_same_KO <= 3),
        median_rank = median(d$rank_of_same_KO), mean_r_same = mean(d$r_same_KO), mean_r_other = mean(d$mean_r_other_KOs))))
print(format(agg, digits = 3), row.names = FALSE)
print(tapply(res$rank_of_same_KO[res$features == "modules_centred"], res$cell_line[res$features == "modules_centred"], median))
