# Spec 93: why the CROP-seq and bulk NCI-H2052 KO x Cm heatmaps (fig9C) look alike: decompose each into
#   shared  = the module average over the 8 KOs (response common to all KOs vs NTC)  and
#   specific = each cell minus its module average (what distinguishes one KO from the others).
# Same values (median NTC - median KO module score), palette and order as fig9C (spec81).
# Panels: full | shared (one strip per data set) | specific, for screen and bulk H2052; scatters of the 20 module averages and of the 160 KO-specific cells.
# needs: spec81_Cm_KO_heatmap_bulk_like_cropseq.R (spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv)
# outputs: spec93_shared_vs_specific.csv, fig9O_Cm_heatmap_shared_vs_specific.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(grid) })
PAL <- rev(as.character(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100))); CM <- paste0("Cm", 1:20)
x <- read.csv(file.path(OUT, "spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv"))
W <- function(s) { d <- x[x$source == s, ]; M <- tapply(d$NTC_minus_KO, list(d$KO, d$Cm), mean); M[ORDER, CM] }
S <- W("screen"); B <- W("NCI-H2052")
ro <- hclust(dist(S))$order; co <- hclust(dist(t(S)))$order; S <- S[ro, co]; B <- B[ro, co]                      # fig9C order (clustering of the screen)
dec <- function(M) list(full = M, shared = matrix(colMeans(M), 1, ncol(M), dimnames = list("mean of 8 KOs", colnames(M))), specific = sweep(M, 2, colMeans(M)))
DS <- dec(S); DB <- dec(B); sp <- function(a, b) cor.test(as.vector(a), as.vector(b), method = "spearman", exact = FALSE)
r_full <- sp(S, B); r_sh <- sp(DS$shared, DB$shared); r_spec <- sp(DS$specific, DB$specific)
frac <- function(D) sum(outer(rep(1, nrow(D$full)), as.vector(D$shared))^2) / sum(D$full^2)
C <- cor(t(S), t(B), method = "spearman")
res <- data.frame(part = c("full (160 cells)", "shared: module averages (20)", "KO-specific: centred (160)"), rho = c(r_full$estimate, r_sh$estimate, r_spec$estimate), p = c(r_full$p.value, r_sh$p.value, r_spec$p.value),
                  share_screen = c(1, frac(DS), 1 - frac(DS)), share_bulk = c(1, frac(DB), 1 - frac(DB)))
write.csv(res, file.path(OUT, "spec93_shared_vs_specific.csv"), row.names = FALSE); print(res, digits = 3, row.names = FALSE)
cat(sprintf("KO x KO module-profile correlation (uncentred): own %.2f vs other %.2f\n", mean(diag(C)), mean(C[row(C) != col(C)])))

hm <- function(M, ttl, rn = TRUE, lim = NULL) { if (is.null(lim)) lim <- min(abs(range(M))); M2 <- pmax(pmin(M, lim), -lim)
  Heatmap(M2, name = ttl, col = circlize::colorRamp2(seq(-lim, lim, length.out = 100), PAL), cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, show_row_names = rn,
          row_names_side = "left", row_names_gp = gpar(fontsize = 7, col = if (nrow(M) > 1) KO_COL[rownames(M)] else INK2, fontface = "bold"), column_names_gp = gpar(fontsize = 6.5), column_names_rot = 45,
          column_title = ttl, column_title_gp = gpar(fontsize = 8, fontface = "bold"), heatmap_legend_param = list(title = "NTC - KO", title_gp = gpar(fontsize = 6.5), labels_gp = gpar(fontsize = 6), direction = "horizontal")) }
scatter <- function(a, b, ttl, r, xl, yl, lab = NULL) { pushViewport(plotViewport(c(3.4, 3.6, 2.4, 1))); pushViewport(dataViewport(a, b)); grid.rect(gp = gpar(col = MUTED, fill = NA))
  grid.xaxis(gp = gpar(fontsize = 6.5)); grid.yaxis(gp = gpar(fontsize = 6.5)); grid.lines(unit(c(0, 0), "native"), unit(range(b), "native"), gp = gpar(col = GRID)); grid.lines(unit(range(a), "native"), unit(c(0, 0), "native"), gp = gpar(col = GRID))
  grid.points(unit(a, "native"), unit(b, "native"), pch = 16, size = unit(if (length(a) > 40) 1.4 else 2.4, "mm"), gp = gpar(col = adjustcolor(INK2, 0.75)))
  if (!is.null(lab)) grid.text(lab, unit(a, "native"), unit(b, "native") + unit(2.4, "mm"), gp = gpar(fontsize = 5.5, col = INK2))
  grid.text(xl, y = unit(-2.5, "lines"), gp = gpar(fontsize = 7)); grid.text(yl, x = unit(-2.8, "lines"), rot = 90, gp = gpar(fontsize = 7))
  grid.text(sprintf("%s: Spearman %.2f (p %s)", ttl, r$estimate, formatC(r$p.value, format = "g", digits = 2)), y = unit(1, "npc") + unit(0.9, "lines"), gp = gpar(fontsize = 7.5, fontface = "bold")); popViewport(2) }
save_fig("fig9O_Cm_heatmap_shared_vs_specific", 17, 9.4, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(4, 4, heights = unit(c(1.3, 1, 0.42, 1), c("cm", "null", "null", "null")), widths = unit(c(1, 1, 1, 0.9), "null"))))
  grid.text("Why the CROP-seq and bulk NCI-H2052 KO x Cm heatmaps look alike: shared response vs KO-specific part", x = unit(4, "mm"), y = unit(0.65, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  grid.text(sprintf("NTC - KO module score (orange = lower after KO), fig9C order; shared = each module's average over the 8 KOs (%.0f%% of the screen heatmap, %.0f%% of bulk); KO-specific = each cell minus that average", 100 * frac(DS), 100 * frac(DB)),
            x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  dr <- function(h, r, c) { pushViewport(viewport(layout.pos.row = r, layout.pos.col = c)); draw(h, newpage = FALSE, padding = unit(c(1, 2, 1, 2), "mm"), heatmap_legend_side = "bottom"); popViewport() }
  dr(hm(DS$full, "CROP-seq: full"), 2, 1); dr(hm(DS$specific, "CROP-seq: KO-specific"), 2, 2); SH <- rbind("CROP-seq" = as.vector(DS$shared) / max(abs(DS$shared)), "bulk H2052" = as.vector(DB$shared) / max(abs(DB$shared))); colnames(SH) <- colnames(S)
  dr(Heatmap(SH, name = "shared", col = circlize::colorRamp2(seq(-1, 1, length.out = 100), PAL), cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 7, fontface = "bold"),
             column_names_gp = gpar(fontsize = 6.5), column_names_rot = 45, column_title = sprintf("shared part = module average over the 8 KOs (each row scaled to its max |value|); screen vs bulk rho %.2f", r_sh$estimate),
             column_title_gp = gpar(fontsize = 8, fontface = "bold"), heatmap_legend_param = list(title = "scaled", title_gp = gpar(fontsize = 6.5), labels_gp = gpar(fontsize = 6), direction = "horizontal"), height = unit(1.1, "cm")), 3, 1:2)
  dr(hm(DB$full, "bulk H2052: full"), 4, 1); dr(hm(DB$specific, "bulk H2052: KO-specific"), 4, 2)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3)); scatter(as.vector(DS$full), as.vector(DB$full), "full heatmap (160 cells)", r_full, "CROP-seq", "bulk H2052"); popViewport()
  pushViewport(viewport(layout.pos.row = 3:4, layout.pos.col = 3)); scatter(as.vector(DS$shared), as.vector(DB$shared), "shared: 20 module averages", r_sh, "CROP-seq module average", "bulk H2052 module average", colnames(S)); popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4)); scatter(as.vector(DS$specific), as.vector(DB$specific), "KO-specific (160 cells)", r_spec, "CROP-seq KO-specific", "bulk H2052 KO-specific"); popViewport()
  pushViewport(viewport(layout.pos.row = 3:4, layout.pos.col = 4)); pushViewport(plotViewport(c(3, 1, 2.4, 1)))
  grid.text(paste("Reading", "", sprintf("Full heatmaps agree (rho %.2f) because", r_full$estimate), sprintf("their shared part agrees (rho %.2f):", r_sh$estimate), "most KOs move the same modules vs NTC.", "",
                  sprintf("The KO-specific part does not (rho %.2f):", r_spec$estimate), sprintf("each screen KO matches the bulk KOs on average"), sprintf("as well as its own (%.2f vs %.2f).", mean(diag(C)), mean(C[row(C) != col(C)])), sep = "\n"),
            x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 7.5, col = INK2)); popViewport(2)
})
