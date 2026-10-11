# Spec 89: do the KOs resemble each other the same way in the CROP-seq screen and in bulk? KO x KO similarity (Spearman of NTC - KO log fold changes centred over the 8 KOs)
# in the screen pseudobulk (spec88) and in bulk NCI-H2052 (and, for reference, each other bulk line), then second-order agreement (RSA / Mantel): Spearman between the
# 28 KO-pair similarities of the screen and of bulk, p = exact permutation of the bulk KO labels (8! = 40320).
#   screen: pseudobulk logCPM (spec88), NTC - KO;  bulk: limma coefficient NTC - KO on TMM logCPM (spec72 S70, no covariate)
#   genes : the screen's top 5000 HVGs (Seurat vst) that are expressed in the screen pseudobulk and in bulk
# needs: spec88_screen_pseudobulk_KO_correlation.R (spec88_screen_pseudobulk_logCPM.csv), spec72_build_profiles.R (S70)
# outputs: spec89_KO_similarity.csv, spec89_RSA.csv, fig9K_KO_similarity_screen_vs_bulk.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })
L <- as.matrix(read.csv(file.path(OUT, "spec88_screen_pseudobulk_logCPM.csv"), row.names = 1, check.names = FALSE)); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF)))
G <- intersect(rownames(L), pr$genes); cat("genes (screen HVG, expressed in screen and bulk):", length(G), "\n")
cen <- function(F) F - rowMeans(F)
SIM <- list(screen = cor(cen(SGN * (L[G, "NTC"] - L[G, KOS])), method = "spearman"))
for (cl in LINES) SIM[[cl]] <- cor(cen(pr$P$tmm$none$coef[G, , cl]), method = "spearman")
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
ut <- which(upper.tri(SIM$screen)); rsa <- function(A, B) cor(A[ut], B[ut], method = "spearman")
RS <- do.call(rbind, lapply(LINES, function(cl) { obs <- rsa(SIM$screen, SIM[[cl]]); nul <- apply(perms, 1, function(p) rsa(SIM$screen, SIM[[cl]][p, p]))
  data.frame(line = cl, rsa_rho = obs, perm_p = mean(nul >= obs - 1e-12)) }))
write.csv(RS, file.path(OUT, "spec89_RSA.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(SIM), function(n) data.frame(source = n, KO_a = KOS[row(SIM[[n]])[ut]], KO_b = KOS[col(SIM[[n]])[ut]], rho = SIM[[n]][ut]))), file.path(OUT, "spec89_KO_similarity.csv"), row.names = FALSE)
options(width = 200); cat("\nscreen KO similarity:\n"); print(round(SIM$screen, 2)); cat("\nbulk H2052 KO similarity:\n"); print(round(SIM[["NCI-H2052"]], 2)); cat("\nRSA screen vs bulk:\n"); print(RS, digits = 3, row.names = FALSE)

ord <- KOS[hclust(as.dist(1 - SIM$screen), "complete")$order]
hm <- function(C, ttl, lim = 0.6) Heatmap(C[ord, ord], name = ttl, col = colorRamp2(c(-lim, 0, lim), c("#2a78d6", "#f0efec", "#eb6834")), cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = ttl == "screen",
  heatmap_legend_param = list(title = "Spearman", title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)), row_names_gp = gpar(fontsize = 7.5, col = KO_COL[ord], fontface = "bold"),
  column_names_gp = gpar(fontsize = 7.5, col = KO_COL[ord], fontface = "bold"), column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 0.6),
  column_title = if (ttl == "screen") "CROP-seq (pseudobulk)" else sprintf("bulk %s\nvs screen: RSA rho %.2f, p %s", ttl, RS$rsa_rho[RS$line == ttl], formatC(RS$perm_p[RS$line == ttl], format = "g", digits = 2)),
  column_title_gp = gpar(fontsize = 8, fontface = "bold"), cell_fun = function(j, i, x, y, w, h, fill) { v <- C[ord, ord][i, j]; if (i != j) grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 5.2)) })
save_fig("fig9K_KO_similarity_screen_vs_bulk", 17, 8.2, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 4, heights = unit(c(1.2, 1, 1), c("cm", "null", "null")))))
  grid.text("Do the KOs resemble each other the same way in CROP-seq and in bulk? KO x KO similarity of NTC - KO profiles (centred over KOs)", x = unit(4, "mm"), y = unit(0.65, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  grid.text(sprintf("Spearman over %d genes (screen top-5000 HVGs expressed in both); all panels in the CROP-seq cluster order; RSA = Spearman between the 28 KO-pair similarities, p = exact permutation of bulk KO labels", length(G)), x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  dr <- function(h, r, c) { pushViewport(viewport(layout.pos.row = r, layout.pos.col = c)); draw(h, newpage = FALSE, padding = unit(c(2, 2, 2, 2), "mm")); popViewport() }
  dr(hm(SIM$screen, "screen"), 2, 1); dr(hm(SIM[["NCI-H2052"]], "NCI-H2052"), 2, 2)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3:4)); pushViewport(plotViewport(c(3.6, 4, 3, 1))); x <- SIM$screen[ut]; y <- SIM[["NCI-H2052"]][ut]
  pushViewport(dataViewport(x, y)); grid.xaxis(gp = gpar(fontsize = 6.5)); grid.yaxis(gp = gpar(fontsize = 6.5)); grid.rect(gp = gpar(col = MUTED, fill = NA))
  grid.lines(unit(c(0, 0), "native"), unit(range(y), "native"), gp = gpar(col = GRID)); grid.lines(unit(range(x), "native"), unit(c(0, 0), "native"), gp = gpar(col = GRID))
  pa <- paste(KOS[row(SIM$screen)[ut]], KOS[col(SIM$screen)[ut]], sep = "-"); grid.points(unit(x, "native"), unit(y, "native"), pch = 16, size = unit(2.2, "mm"), gp = gpar(col = INK2))
  grid.text(pa, unit(x, "native"), unit(y, "native") + unit(2.2, "mm"), gp = gpar(fontsize = 5, col = INK2))
  grid.text("KO-pair similarity in CROP-seq", y = unit(-2.6, "lines"), gp = gpar(fontsize = 7.5)); grid.text("KO-pair similarity in bulk NCI-H2052", x = unit(-2.8, "lines"), rot = 90, gp = gpar(fontsize = 7.5))
  grid.text(sprintf("28 KO pairs: RSA rho %.2f, permutation p %s", RS$rsa_rho[RS$line == "NCI-H2052"], formatC(RS$perm_p[RS$line == "NCI-H2052"], format = "g", digits = 2)), y = unit(1, "npc") + unit(1, "lines"), gp = gpar(fontsize = 8, fontface = "bold")); popViewport(3)
  for (j in 1:3) dr(hm(SIM[[LINES[j + 1]]], LINES[j + 1]), 3, j)
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 4)); pushViewport(plotViewport(c(4, 4, 3, 1))); pushViewport(dataViewport(xscale = c(0.4, 4.6), yscale = c(-0.6, 0.8)))
  grid.yaxis(gp = gpar(fontsize = 6.5)); grid.lines(unit(c(0.4, 4.6), "native"), unit(c(0, 0), "native"), gp = gpar(col = INK2, lty = 2))
  grid.rect(unit(1:4, "native"), unit(0, "native"), unit(0.6, "native"), unit(RS$rsa_rho, "native"), just = c("centre", "bottom"), gp = gpar(fill = LINE_COL[RS$line], col = NA))
  grid.text(sprintf("p %s", formatC(RS$perm_p, format = "g", digits = 2)), unit(1:4, "native"), unit(pmax(RS$rsa_rho, 0) + 0.06, "native"), gp = gpar(fontsize = 6.5, fontface = ifelse(RS$perm_p < 0.05, "bold", "plain")))
  grid.text(sub("NCI-|-211H", "", RS$line), unit(1:4, "native"), unit(-1, "lines"), gp = gpar(fontsize = 7)); grid.text("RSA rho (screen vs bulk)", x = unit(-2.8, "lines"), rot = 90, gp = gpar(fontsize = 7.5))
  grid.text("KO-similarity structure vs screen, per line", y = unit(1, "npc") + unit(1, "lines"), gp = gpar(fontsize = 8, fontface = "bold")); popViewport(3)
})
