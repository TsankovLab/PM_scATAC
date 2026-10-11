# Spec 42: within the TWIST1-guide cells specifically (not KO-vs-NTC, not KO-vs-other-KO -- just cells that got a TWIST1 guide),
# is residual TWIST1 mRNA expression (CRISPR KO efficiency is incomplete per-cell) correlated with cNMF module 16 (k=30, the EMT/
# TEAD4 module) usage? A positive correlation here would be direct single-cell evidence that TWIST1 expression level itself drives
# module 16, independent of any between-group (KO vs NTC) comparison.
# needs: spec_common.R (load_screen()), CRISPR_cropseq_analysis/.../cNMF_crispr_20_30_vf5000/cnmf/cnmf.usages.k_30.dt_0_3.consensus.txt
# outputs: spec42_TWIST1_module16_within_KO_corr.csv, fig2Z_TWIST1_module16_within_KO_scatter.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
CNMF_D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/cNMF/cNMF_crispr_20_30_vf5000/cnmf"
MODULE <- "16"; K <- "30"

srt <- load_screen(); Y <- GetAssayData(srt, layer = "data")
grp <- setNames(as.character(srt$merged_call), colnames(srt))
usages <- read.delim(file.path(CNMF_D, sprintf("cnmf.usages.k_%s.dt_0_3.consensus.txt", K)), row.names = 1, check.names = FALSE)

twist1_cells <- names(grp)[grp == "TWIST1"]
twist1_cells <- intersect(twist1_cells, rownames(usages))
cat(sprintf("TWIST1-guide cells: %d total, %d with a cNMF usage value\n", sum(grp == "TWIST1"), length(twist1_cells)))

twist1_expr <- as.numeric(Y["TWIST1", twist1_cells])
m16_usage <- usages[twist1_cells, MODULE]
cat(sprintf("TWIST1 expression in these cells: %d/%d non-zero (%.0f%%), mean=%.3f\n",
            sum(twist1_expr > 0), length(twist1_expr), 100 * mean(twist1_expr > 0), mean(twist1_expr)))

pr <- cor.test(twist1_expr, m16_usage, method = "pearson")
sr <- cor.test(twist1_expr, m16_usage, method = "spearman", exact = FALSE)
cat(sprintf("\nPearson  r = %.3f, p = %.4g\n", pr$estimate, pr$p.value))
cat(sprintf("Spearman rho = %.3f, p = %.4g\n", sr$estimate, sr$p.value))

## same check restricted to cells with detectable (non-zero) TWIST1 expression only, in case the zero-inflation dominates
nz <- twist1_expr > 0
if (sum(nz) > 10) {
  pr_nz <- cor.test(twist1_expr[nz], m16_usage[nz], method = "pearson")
  sr_nz <- cor.test(twist1_expr[nz], m16_usage[nz], method = "spearman", exact = FALSE)
  cat(sprintf("\nrestricted to TWIST1-detected cells (n=%d): Pearson r = %.3f (p=%.4g), Spearman rho = %.3f (p=%.4g)\n",
              sum(nz), pr_nz$estimate, pr_nz$p.value, sr_nz$estimate, sr_nz$p.value))
}

write.csv(data.frame(cell = twist1_cells, TWIST1_expr = twist1_expr, module16_usage = m16_usage),
          file.path(OUT, "spec42_TWIST1_module16_within_KO_corr.csv"), row.names = FALSE)

save_fig("fig2Z_TWIST1_module16_within_KO_scatter", 6.5, 6, function() {
  par(mar = c(4.5, 4.5, 3.6, 1), family = "sans")
  plot(twist1_expr, m16_usage, pch = 16, cex = 0.5, col = adjustcolor(KO_COL["TWIST1"], 0.5),
       xlab = "TWIST1 expression (log-norm.)", ylab = "cNMF module 16 usage", las = 1)
  if (sum(twist1_expr > 0) > 10) { fit <- lm(m16_usage ~ twist1_expr); abline(fit, col = INK, lwd = 2) }
  mtext("Within TWIST1-guide cells: TWIST1 expression vs module 16 (EMT) usage", side = 3, adj = 0, line = 2.2, font = 2, cex = 0.95)
  mtext(sprintf("n=%d cells; Pearson r=%.3f (p=%.3g), Spearman rho=%.3f (p=%.3g)", length(twist1_expr), pr$estimate, pr$p.value, sr$estimate, sr$p.value),
        side = 3, adj = 0, line = 0.8, cex = 0.65, col = INK2)
})
