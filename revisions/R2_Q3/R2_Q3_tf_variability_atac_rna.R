###############################################################
# R2_Q3 : BETWEEN-SUBCLONE VARIABILITY OF TF ACTIVITY (scATAC chromVAR) PAIRED WITH
#         BETWEEN-SUBCLONE VARIABILITY OF THE SAME TF's RNA EXPRESSION.
#
# Clones are NOT matched across modalities -- and they do not need to be. Each modality
# measures, within its own cells, how much a TF differs between that sample's two
# subclones; the two numbers are then averaged over samples independently. This sidesteps
# the cross-modal clone-identity problem entirely (see R2_Q3_clone_anchor.R, where
# matching failed) while still asking whether the TFs whose MOTIF ACTIVITY varies between
# clones are the same TFs whose EXPRESSION varies between clones.
#
#   ATAC variability = MEDIAN over samples of |mean chromVAR z(sc2) - mean z(sc1)|
#                      (from chromvar_armclones_diff_<S>.csv, current arm-level clones)
#   RNA  variability = MEDIAN over samples of |mean lognorm expr(sc2) - mean expr(sc1)|
#                      for the TF's own gene, over that sample's scRNA subclones
#
# MEDIAN, not mean: the per-sample dynamic range of |diff| differs ~8-fold (median 0.082
# in P5 to 0.636 in P11; sd 0.124 to 0.773), so the mean is set by whichever sample
# separates its clones most -- P11's |diff| vector correlates 0.83 with the cross-sample
# mean, P10's only 0.27. Concretely, the mean's top-15 contains four HOX13 paralogues
# (HOXA/B/C/D13) that are P11-specific and drop out under the median. The median asks
# instead that the effect be present in at least half the samples, which is what
# "coherent across samples" should mean. Both aggregators plus the scale-free
# mean-of-within-sample-percentiles are written to the CSV.
#
# Sample sets differ by modality on purpose (no pairing): scATAC P1,P4,P5,P8,P10,P11,P12;
# scRNA P1,P4,P5,P8,P11,P12,P13.
#
# The AP-1/bZIP + SMARCC1 motifs are FLAGGED: they are the chromVAR accessibility axis
# (anti-correlated with nFrags, R2_Q3_tead_vs_nfrags.R), not clonal biology. Whether
# their RNA varies too is exactly the kind of thing this figure can settle.
#
# Outputs: tf_variability_atac_rna.csv, Plots/R2_Q3_tf_variability_atac_rna.pdf
###############################################################
suppressMessages({ library(cluster); library(Seurat); library(Matrix)
                   library(ggplot2); library(ggrepel); library(patchwork) })
set.seed(1)
SC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(SC, "git_repo_claude", "R2_Q3")); dir.create("Plots", showWarnings = FALSE)
SRC <- "arm_cnv_zscore"; KEEP <- 10
COL_TECH <- "#b0413e"; COL_OTH <- "grey62"
COL_IMM <- "#7b3fa0"
TECH <- c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","BACH1","BACH2","NFE2",
          "SMARCC1","JDP2","BATF","BATF3","NFE2L1","NFE2L2","MAF","MAFF","MAFK","MAFG")
## curated immune / interferon / NF-kB / lymphoid regulators. BATF and BATF3 are immune
## but share the AP-1 motif family, so they stay in TECH (conservative).
IMMUNE <- c(paste0("IRF", 1:9), paste0("STAT", c(1,2,3,4,"5A","5B",6)),
            "NFKB1","NFKB2","REL","RELA","RELB","SPI1","SPIB","SPIC",
            paste0("IKZF", 1:4), paste0("NFATC", 1:4), "TBX21","EOMES","PRDM1",
            "RUNX3","ETS1","ELF4","BCL11B","TCF7","LEF1","FOXP3","POU2F2","POU2AF1",
            "CIITA","NLRC5")

## ---------------- STEP 1: current arm-level clone calls ---------------------
## Only the scRNA clones are rebuilt here; the scATAC side is read from the chromVAR
## per-sample tables written by R2_Q3_chromvar_armclones.R (same clone definition).
call_clones_arm <- function(A, kmax = 6){
  ok <- apply(A, 2, function(z) all(is.finite(z))); A <- A[, ok, drop = FALSE]
  none <- list(k = 1, sil = NA, cl = setNames(rep(1L, ncol(A)), colnames(A)))
  if (ncol(A) < 40 || nrow(A) < 2) return(none)
  Ac <- A - rowMeans(A); d <- dist(t(Ac))
  hc <- hclust(d, "ward.D2"); best <- none
  for (k in 2:min(kmax, ncol(A) - 1)){
    cl <- cutree(hc, k); if (length(unique(cl)) < k) next
    if (any(table(cl) < max(20, 0.05 * ncol(A)))) next
    s <- mean(silhouette(cl, d)[, 3])
    if (is.na(best$sil) || s > best$sil) best <- list(k = k, sil = s, cl = cl) }
  if (is.na(best$sil) || best$sil < 0.10)
    best <- list(k = 1, sil = best$sil, cl = setNames(rep(1L, ncol(A)), colnames(A)))
  best }
robust2 <- function(assign){
  tb <- sort(table(assign), decreasing = TRUE); thr <- max(20, 0.05 * length(assign))
  rob <- names(tb)[tb >= thr][1:2]
  out <- rep(NA_character_, length(assign)); names(out) <- names(assign)
  out[assign == rob[1]] <- "sc1"; out[assign == rob[2]] <- "sc2"; out[!is.na(out)] }

MM <- read.csv("arm_level_multimodal_arms.csv", stringsAsFactors = FALSE)
MM$pct <- ave(MM$bic_gain, MM$sample, MM$modality, FUN = function(z){
  z[!is.finite(z)] <- min(z[is.finite(z)], na.rm = TRUE); rank(z) / length(z) })
MM$adj <- MM$pct - ave(MM$pct, MM$arm, MM$modality,
  FUN = function(z) if (length(z) > 1) (sum(z) - z) / (length(z) - 1) else 0)

rna_clones <- list()
for (f in list.files(SRC, "_scRNA_arm\\.rds$", full.names = TRUE)){
  o <- readRDS(f); A <- o$A
  t <- MM[MM$sample == o$sample & MM$modality == "scRNA", ]
  v <- setNames(t$adj[match(rownames(A), t$arm)], rownames(A)); v[!is.finite(v)] <- -Inf
  keep <- names(sort(v, decreasing = TRUE))[1:min(KEEP, nrow(A))]
  cc <- call_clones_arm(A[rownames(A) %in% keep, , drop = FALSE])
  if (cc$k > 1) rna_clones[[o$sample]] <- robust2(cc$cl) }
cat("scRNA-subclonal samples:", paste(sort(names(rna_clones)), collapse = ", "), "\n")

## ---------------- STEP 2: ATAC variability per TF ---------------------------
## Collapse the seven per-sample tables to one number per motif. MEDIAN is the headline
## (see header); the mean and the scale-free mean-of-percentiles are carried alongside so
## the choice of aggregator can be checked rather than trusted. q25/q75 give the IQR
## error bars in the figure -- a spread across samples, NOT a confidence interval.
af <- list.files(".", "^chromvar_armclones_diff_.*\\.csv$")
A <- do.call(rbind, lapply(af, read.csv, stringsAsFactors = FALSE))
A$pct_in_sample <- ave(abs(A$diff), A$sample, FUN = function(z) rank(z) / length(z))
atac <- do.call(rbind, lapply(split(A, A$TF), function(x)
  data.frame(TF = x$TF[1], atac_n = nrow(x),
             atac_var = median(abs(x$diff)),               # HEADLINE: median across samples
             atac_var_mean = mean(abs(x$diff)),            # kept for comparison
             atac_var_meanpct = mean(x$pct_in_sample),     # scale-free alternative
             atac_q25 = quantile(abs(x$diff), 0.25, names = FALSE),
             atac_q75 = quantile(abs(x$diff), 0.75, names = FALSE),
             atac_eta2 = median(x$eta2, na.rm = TRUE),
             atac_nsig = sum(x$fdr < 0.05, na.rm = TRUE), stringsAsFactors = FALSE)))
cat("ATAC: ", nrow(atac), " motifs over ", length(unique(A$sample)), " samples\n", sep = "")

## ---------------- STEP 3: RNA variability per TF ----------------------------
## Same quantity on the expression side: |mean expression of clone2 - clone1| for the
## TF's OWN gene, per sample, then median across samples. Clones are never matched to the
## ATAC clones -- each assay measures its own between-clone spread, which is the whole
## point (cross-modal clone identity is not recoverable; see the tree script).
## `detect` records the fraction of cells expressing the gene, so a motif that is variable
## but whose TF is barely transcribed can be spotted.
srt <- readRDS(file.path(SC, "tumor_compartment", "scrna", "scRNA_meso.rds"))
E <- GetAssayData(srt, slot = "data")
genes <- intersect(rownames(E), atac$TF)
cat("TF motifs with a matching gene in the scRNA object: ", length(genes), "\n", sep = "")
E <- E[genes, , drop = FALSE]

rna_list <- list(); det <- list()
for (s in names(rna_clones)){
  lab <- rna_clones[[s]]
  cells <- intersect(names(lab), colnames(E))
  if (length(cells) < 40){ cat("  skip", s, "- few cells matched\n"); next }
  lab <- lab[cells]; e <- E[, cells, drop = FALSE]
  m1 <- Matrix::rowMeans(e[, lab == "sc1", drop = FALSE])
  m2 <- Matrix::rowMeans(e[, lab == "sc2", drop = FALSE])
  rna_list[[s]] <- abs(m2 - m1)
  det[[s]] <- Matrix::rowMeans(e > 0)
  cat(sprintf("  %s: %d cells (sc1 %d / sc2 %d)\n", s, length(cells),
              sum(lab == "sc1"), sum(lab == "sc2"))) }
R <- do.call(cbind, rna_list); DET <- do.call(cbind, det)
Rp <- apply(R, 2, function(z) rank(z) / length(z))
rna <- data.frame(TF = rownames(R), rna_n = ncol(R),
                  rna_var = apply(R, 1, median),           # HEADLINE: median across samples
                  rna_var_mean = rowMeans(R),
                  rna_var_meanpct = rowMeans(Rp),
                  rna_q25 = apply(R, 1, quantile, 0.25, names = FALSE),
                  rna_q75 = apply(R, 1, quantile, 0.75, names = FALSE),
                  detect = rowMeans(DET),
                  expr = Matrix::rowMeans(E), stringsAsFactors = FALSE)

## ---------------- STEP 4: merge, classify, rank -----------------------------
## Percentiles are taken WITHIN each assay before combining, since chromVAR z and log
## expression are not on a common scale.
d <- merge(atac, rna, by = "TF")
d$class <- ifelse(d$TF %in% TECH, "AP-1/bZIP (technical)",
           ifelse(d$TF %in% IMMUNE, "immune / interferon", "other TF"))
d$atac_pct <- rank(d$atac_var) / nrow(d)
d$rna_pct  <- rank(d$rna_var)  / nrow(d)
d$combined <- (d$atac_pct + d$rna_pct) / 2
d <- d[order(-d$combined), ]
write.csv(d, "tf_variability_atac_rna.csv", row.names = FALSE)

## how much does the choice of aggregator matter?
cat(sprintf("\naggregator agreement (ATAC): median~mean %.3f | median~meanPct %.3f\n",
            cor(d$atac_var, d$atac_var_mean, method = "spearman"),
            cor(d$atac_var, d$atac_var_meanpct, method = "spearman")))
t15 <- function(v) d$TF[order(-v)][1:15]
cat("  top15 by median  :", paste(t15(d$atac_var), collapse = ", "), "\n")
cat("  top15 by mean    :", paste(t15(d$atac_var_mean), collapse = ", "), "\n")
cat(sprintf("  overlap: %d/15\n", length(intersect(t15(d$atac_var), t15(d$atac_var_mean)))))

rho <- cor(d$atac_var, d$rna_var, method = "spearman")
rho_bio <- with(d[d$class != "AP-1/bZIP (technical)", ],
                cor(atac_var, rna_var, method = "spearman"))
cat(sprintf("\nSpearman(ATAC variability, RNA variability) = %.3f  (excluding the technical axis: %.3f)\n",
            rho, rho_bio))
cat("\n=== top 20 by combined variability rank ===\n")
print(head(d[, c("TF","class","atac_var","atac_var_mean","rna_var","detect","atac_pct","rna_pct")], 20),
      row.names = FALSE, digits = 3)
## STEP 5: is the immune class over-represented among the most variable TFs?
## Fisher test of immune membership against top-decile / top-5% status.
for (q in c(0.90, 0.95)){
  hit <- d$combined >= quantile(d$combined, q)
  tb <- table(immune = d$class == "immune / interferon", top = hit)
  ft <- fisher.test(tb)
  cat(sprintf("immune/interferon TFs in the top %.0f%% by combined rank: %d of %d (%.0f%%) vs %.0f%% background | OR %.2f, p = %.2g\n",
              100 * (1 - q), tb["TRUE","TRUE"], sum(hit), 100 * tb["TRUE","TRUE"] / sum(hit),
              100 * mean(d$class == "immune / interferon"), ft$estimate, ft$p.value)) }
cat("\nimmune/interferon members in the top 60 by combined rank: ",
    paste(head(d$TF[d$class == "immune / interferon"], 12), collapse = ", "), "\n", sep = "")

cat("\n=== TEAD family ===\n")
print(d[grepl("^TEAD", d$TF), c("TF","atac_var","atac_var_mean","rna_var","detect","atac_pct","rna_pct")],
      row.names = FALSE, digits = 3)

## ---------------- figure -----------------------------------------------------
## Every TF is plotted, but only the highlighted ones carry error bars and labels --
## ~800 interquartile crosses would be unreadable. Bars = IQR of the per-sample
## |between-clone difference| (n = 7 samples per assay), i.e. the spread behind the
## median point estimate, NOT a confidence interval.
pal <- c(`AP-1/bZIP (technical)` = COL_TECH, `immune / interferon` = COL_IMM,
         `other TF` = COL_OTH)
NTOP <- 12
## highlighted = the top non-technical TFs, the top immune/interferon TFs, and a few of
## the technical ones for reference. TEAD is no longer given its own class; it appears
## only if it ranks into the top set on its own.
hi <- unique(c(head(d$TF[d$class == "other TF"], NTOP),
               head(d$TF[d$class == "immune / interferon"], 8),
               head(d$TF[d$class == "AP-1/bZIP (technical)"], 4)))
H <- d[d$TF %in% hi, ]
H$class <- factor(H$class, levels = names(pal))

p1 <- ggplot(d, aes(atac_var, rna_var)) +
  geom_point(colour = "grey86", size = 0.5, alpha = 0.7) +
  geom_linerange(data = H, aes(xmin = atac_q25, xmax = atac_q75, colour = class),
                 linewidth = 0.35, alpha = 0.85) +
  geom_linerange(data = H, aes(ymin = rna_q25, ymax = rna_q75, colour = class),
                 linewidth = 0.35, alpha = 0.85) +
  geom_point(data = H, aes(colour = class, size = detect), stroke = 0) +
  geom_text_repel(data = H, aes(label = TF, colour = class), size = 2.4,
                  max.overlaps = 40, min.segment.length = 0, segment.size = 0.2,
                  box.padding = 0.3, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_size_continuous(range = c(0.9, 3.2), name = "fraction of\ncells detected") +
  labs(x = "scATAC: median |between-clone difference| in chromVAR z",
       y = "scRNA: median |between-clone difference| in expression",
       title = "TF variability between subclones, measured independently in each assay",
       subtitle = sprintf("median over samples, bars = IQR across the 7 samples | clones not matched across assays | Spearman rho = %.2f (%.2f excluding the technical axis)",
                          rho, rho_bio)) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 6.5),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 6.8, colour = "grey35"))

## forest panels: same TF order in both assays, median +/- IQR
ord <- H$TF[order(H$combined)]
F <- H; F$TF <- factor(F$TF, levels = ord)
forest <- function(mid, lo, up, xlab)
  ggplot(F, aes(.data[[mid]], TF, colour = class)) +
    geom_linerange(aes(xmin = .data[[lo]], xmax = .data[[up]]), linewidth = 0.4) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6.2),
          plot.title = element_text(face = "bold", size = 8.5))
p2 <- forest("atac_var", "atac_q25", "atac_q75", "chromVAR z") + ggtitle("scATAC motif")
p3 <- forest("rna_var",  "rna_q25",  "rna_q75",  "expression") + ggtitle("scRNA") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave("Plots/R2_Q3_tf_variability_atac_rna.pdf",
       (p1 | p2 | p3) + plot_layout(widths = c(1.7, 0.75, 0.6)),
       width = 12.2, height = 5.4, device = cairo_pdf)
cat("\nDONE -> Plots/R2_Q3_tf_variability_atac_rna.pdf\n")
