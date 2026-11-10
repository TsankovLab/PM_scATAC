# plot_benchmark_summary.R
#
# One-row summary of the cHub benchmarking for the reviewer (R1_Q7), two panels:
#   A  which super-enhancer (SE) sets are recovered: -log10 FDR of the 8 matched cell type | SE pairs,
#      7 hub methods, setting of the paper (12.5 kb, >= 5 peaks, no TSS filter)      (myeloid_SE_hub_enrichment.csv)
#   B  SE recovery and SE specificity vs maximum distance (co-accessible pair distance / stitching gap), in each of the
#      6 settings (TSS filter x minimum peaks)                                        (plot_max_dist_sweep_table.csv)
# specificity = mean -log10 FDR of the matched pairs minus that of the same cell types against the non-matching benchmark
#               SE sets (each -log10 FDR capped at the 99th percentile)
# Run myeloid_SE_hub_enrichment.R, max_dist_sweep_SE.R and max_dist_sweep_summary.R first.
# Output: plot_benchmark_summary.pdf/.png
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"
methods <- c("r 0.2", "r 0.3", "r 0.4", "r 0.5", "r 0.6", "Stitch-AP", "Stitch CT")
mcol <- c(setNames(viridisLite::viridis(5, end = 0.85), methods[1:5]), `Stitch-AP` = "#d7301f", `Stitch CT` = "#fd8d3c")
th <- theme_classic(base_size = 7) + theme(plot.title = element_text(size = 7.5, face = "bold"), legend.key.size = unit(3, "mm"),
                                            plot.margin = margin(2, 4, 2, 2))

## ---- A: matched pairs, paper setting ------------------------------------------------------------
res <- read.csv(file.path(script_dir, "myeloid_SE_hub_enrichment.csv"), stringsAsFactors = FALSE)
res <- res[!duplicated(res[, c("condition", "celltype", "se_db")]), ]
cap <- unname(quantile(res$neg_log10_fdr, 0.99, na.rm = TRUE))
pairs <- data.frame(celltype = c("Alveolar", "Alveolar", "Fibroblasts", "SmoothMuscle", "Endothelial", "Myeloid", "T_cells", "B_cells"),
  se_db = c("AT1", "AT2", "ENCODE_Lung_primary_cell_fibroblast-of-lung", "Smoothmuscle", "Endothelial",
            "ENCODE_Blood_primary_cell_CD14-positive-monocyte", "Tcell", "ENCODE_Blood_primary_cell_Bcell"),
  label = c("Alveolar | AT1", "Alveolar | AT2", "Fibroblast | lung fibroblast", "Smooth muscle | smooth muscle",
            "Endothelial | endothelial", "Myeloid | blood monocyte", "T cells | T cell", "B cells | blood B cell"))
cond <- c(`r 0.2` = "CoAccess_r0.2_mp5_tss0", `r 0.3` = "CoAccess_r0.3_mp5_tss0", `r 0.4` = "CoAccess_r0.4_mp5_tss0", `r 0.5` = "CoAccess_r0.5_mp5_tss0",
          `r 0.6` = "CoAccess_r0.6_mp5_tss0", `Stitch-AP` = "Stitching_mp5_tss0", `Stitch CT` = "Stitching_CT_mp5_tss0")
A <- expand.grid(i = seq_len(nrow(pairs)), method = methods, stringsAsFactors = FALSE)
A$score <- mapply(function(i, m) { v <- res$neg_log10_fdr[res$condition == cond[[m]] & res$celltype == pairs$celltype[i] & res$se_db == pairs$se_db[i]]
  if (length(v)) min(v[1], cap) else 0 }, A$i, A$method)
A$label <- factor(pairs$label[A$i], rev(pairs$label)); A$method <- factor(A$method, methods)
pA <- ggplot(A, aes(method, label, fill = score)) + geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(limits = c(0, cap), breaks = c(0, 50, 100), name = "-log10 FDR") + coord_equal() +
  labs(x = NULL, y = NULL, title = "A  SE recovery at 12.5 kb") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
             legend.position = "bottom", legend.key.width = unit(5, "mm"), legend.key.height = unit(2, "mm"))

## ---- B: distance sweep in every setting --------------------------------------------------------
S <- read.csv(file.path(script_dir, "plot_max_dist_sweep_table.csv"), stringsAsFactors = FALSE)
S$method <- factor(S$method, methods); stopifnot(all(c(5, 10, 12.5, 15, 20) %in% S$dist))
## B: SE recovery (on-target) and SE specificity, distance sweep in every setting (TSS filter x minimum peaks);
## rows = metric, columns = setting; the setting used in this study is framed
S$setl <- factor(sprintf("%s\n>= %d peaks", ifelse(S$tss == 0, "no TSS filter", "TSS >= 2 kb"), S$mp),
                 levels = as.vector(outer(3:5, c("no TSS filter", "TSS >= 2 kb"), function(m, t) sprintf("%s\n>= %d peaks", t, m))))
M <- rbind(data.frame(S[, c("method", "dist", "setl")], metric = "SE recovery", value = S$score_on),
           data.frame(S[, c("method", "dist", "setl")], metric = "SE specificity", value = S$spec))
M$metric <- factor(M$metric, unique(M$metric))
frame <- data.frame(setl = factor("no TSS filter\n>= 5 peaks", levels(S$setl)), metric = factor(levels(M$metric), levels(M$metric)))
## heatmap layout: rows = method, columns = maximum distance; one tile panel per metric (rows) and setting (columns).
## Both metrics are in -log10 FDR units and share the colour scale of panel A.
M$method <- factor(M$method, rev(methods)); M$distf <- factor(M$dist, sort(unique(M$dist)), labels = c("5", "10", "12.5", "15", "20"))
pB <- ggplot(M, aes(distf, method, fill = pmin(pmax(value, 0), cap))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_rect(data = frame, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.7) +
  geom_tile(data = M[M$dist == 12.5, ], fill = NA, colour = "grey15", linewidth = 0.15, linetype = 3) +
  scale_fill_viridis_c(limits = c(0, cap), guide = "none") +
  facet_grid(metric ~ setl, switch = "y") +
  labs(x = "maximum distance (kb)", y = NULL, title = "B  all settings (framed = setting used in this study)") +
  th + theme(strip.background = element_blank(), strip.text = element_text(size = 6), strip.placement = "outside",
             strip.text.y.left = element_text(size = 6.5, face = "bold"), panel.spacing = unit(1.5, "mm"),
             axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 5.2, angle = 90, vjust = 0.5, hjust = 1),
             axis.text.y = element_text(size = 5.8))
fig <- (pA | pB) + plot_layout(widths = c(1, 2.4))
ggsave(file.path(script_dir, "plot_benchmark_summary.pdf"), fig, width = 9, height = 3.1)
ggsave(file.path(script_dir, "plot_benchmark_summary.png"), fig, width = 9, height = 3.1, dpi = 300)
print(aggregate(spec ~ method, S, function(x) round(c(median = median(x), min = min(x), max = max(x)), 1)))
message("Saved plot_benchmark_summary.pdf")
