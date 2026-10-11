###############################################################################
# R1_Q3 -- F31: survival analysis of the P4 scRNA SOX9 / TGF-beta module in the three
#          bulk RNA cohorts (TCGA, Bueno, MESOMICS).
#
# Module = markers of the SOX9-expressing malignant cluster of P4 scRNA (cluster 24 of
#   R1_Q3_P4_SOX9_cluster_recheck.R; P4_SOX9_cluster_recheck_markers.csv): genes with
#   adjusted p < 0.01, detected in >= 10% of the cluster's cells, ranked by log2 fold change;
#   the top N_TOP (default 50) are used. The set contains SOX9 and a TGF-beta response
#   (SERPINE1, FN1, INHBA, TGFB1, ...) together with hypoxia targets (CA12, SLC2A1, ...).
#   Sensitivity: top 25 / top 100, the module without its hypoxia genes, and without SOX9.
# Score = mean over module genes of the per-cohort z-score of log2 expression.
# Models (as F29b / F30): Cox per 1 SD of the score and for a median split, alone and
#   adjusted for the scS (sarcomatoid) score; Kaplan-Meier for the median split.
#   A fixed-effect (inverse-variance) pooled per-SD hazard ratio over the three cohorts.
#
# Outputs: F31_P4_SOX9_TGFb_module_genes.csv, F31_P4_SOX9_TGFb_module_survival_table.csv,
#          Plots/F31_P4_SOX9_TGFb_module_survival.pdf
###############################################################################
set.seed(1234)
suppressMessages({ library(survival); library(ggplot2); library(patchwork) })
projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso'); outdir <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
N_TOP <- as.integer(Sys.getenv('N_TOP', 50))
COHORTS <- c(tcga = 'TCGA', bueno = 'Bueno', mesomics = 'MESOMICS'); TIME_SCALE <- c(tcga = 1, bueno = 12, mesomics = 1)  # -> months
hypoxia_genes <- c('CA12', 'SLC2A1', 'HK2', 'SLC16A3', 'EGLN3', 'NDUFA4L2', 'ANGPTL4', 'STC2', 'ADM', 'VEGFC', 'F3',
                   'ABCC3', 'SCD', 'PFKP', 'ENO2')                 # as in R1_Q3_spatial_SOX9_radius_lineplot_hypoxia.R

## ---- module -------------------------------------------------------------------------------
mk <- read.csv(file.path(outdir, 'P4_SOX9_cluster_recheck_markers.csv'))
mk <- mk[mk$p_val_adj < 0.01 & mk$pct.1 >= 0.1 & mk$avg_log2FC > 0, ]; mk <- mk[order(-mk$avg_log2FC), ]
sets <- list(main = head(mk$gene, N_TOP), top25 = head(mk$gene, 25), top100 = head(mk$gene, 100))
sets$no_hypoxia <- setdiff(sets$main, hypoxia_genes); sets$no_SOX9 <- setdiff(sets$main, 'SOX9')
stopifnot('SOX9' %in% sets$main)
write.csv(cbind(rank = seq_len(N_TOP), head(mk, N_TOP)[, c('gene', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')]),
          file.path(outdir, 'F31_P4_SOX9_TGFb_module_genes.csv'), row.names = FALSE)
cat('module (top', N_TOP, '):', sets$main, '\n')

## ---- cohorts ------------------------------------------------------------------------------
expr_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds')); meta_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))
build <- function(co) {
  e <- as.matrix(expr_l[[co]]); m <- meta_l[[co]]
  if (!identical(rownames(m), colnames(e))) m <- m[match(colnames(e), m$Sample), ]
  if (max(e, na.rm = TRUE) > 50) e <- log2(e + 1)                       # linear-scale matrices are logged
  d <- data.frame(sample_id = colnames(e), scS = m$sarc_score, time = m$census * TIME_SCALE[[co]], status = m$status)
  for (s in names(sets)) { g <- intersect(sets[[s]], rownames(e)); z <- t(scale(t(e[g, , drop = FALSE]))); z <- z[rowSums(is.na(z)) == 0, , drop = FALSE]
    d[[s]] <- colMeans(z); attr(d, paste0('n_', s)) <- nrow(z) }
  ok <- complete.cases(d[, c('scS', 'time', 'status', 'main')]) & d$time > 0
  structure(d[ok, ], n_main = attr(d, 'n_main'), logged = max(as.matrix(expr_l[[co]]), na.rm = TRUE) > 50)
}
dfs <- setNames(lapply(names(COHORTS), build), names(COHORTS))

## ---- models -------------------------------------------------------------------------------
hr <- function(fit, term) { s <- summary(fit); c(HR = s$conf.int[term, 'exp(coef)'], lo = s$conf.int[term, 'lower .95'], hi = s$conf.int[term, 'upper .95'],
                                               p = s$coefficients[term, 'Pr(>|z|)'], beta = s$coefficients[term, 'coef'], se = s$coefficients[term, 'se(coef)']) }
TAB <- list(); KM <- list()
for (co in names(dfs)) for (s in names(sets)) {
  d <- dfs[[co]]; d$z <- as.numeric(scale(d[[s]])); d$scSz <- as.numeric(scale(d$scS))
  d$grp <- factor(ifelse(d[[s]] > median(d[[s]]), 'high', 'low'), levels = c('low', 'high'))
  rows <- rbind(c(contrast = 'per 1 SD', model = 'alone'), c('per 1 SD', 'adjusted for scS'), c('high vs low', 'alone'), c('high vs low', 'adjusted for scS'))
  fits <- list(hr(coxph(Surv(time, status) ~ z, d), 'z'), hr(coxph(Surv(time, status) ~ z + scSz, d), 'z'),
               hr(coxph(Surv(time, status) ~ grp, d), 'grphigh'), hr(coxph(Surv(time, status) ~ grp + scSz, d), 'grphigh'))
  TAB[[paste(co, s)]] <- data.frame(set = s, cohort = COHORTS[[co]], rows, do.call(rbind, fits), n = nrow(d), events = sum(d$status),
    n_genes = length(intersect(sets[[s]], rownames(expr_l[[co]]))), rho_scS = cor(d[[s]], d$scS, method = 'spearman'),
    logrank_p = 1 - pchisq(survdiff(Surv(time, status) ~ grp, d)$chisq, 1), row.names = NULL)
  if (s == 'main') { km <- survfit(Surv(time, status) ~ grp, d); st <- rep(sub('grp=', '', names(km$strata)), km$strata)
    KM[[co]] <- rbind(data.frame(cohort = COHORTS[[co]], grp = c('low', 'high'), time = 0, surv = 1), data.frame(cohort = COHORTS[[co]], grp = st, time = km$time, surv = km$surv)) }
}
TAB <- do.call(rbind, TAB); rownames(TAB) <- NULL
## fixed-effect pooled per-SD hazard ratio
pool <- do.call(rbind, lapply(split(TAB[TAB$contrast == 'per 1 SD', ], ~ set + model, drop = TRUE), function(x) { w <- 1 / x$se^2; b <- sum(w * x$beta) / sum(w); se <- sqrt(1 / sum(w))
  data.frame(set = x$set[1], cohort = 'Pooled', contrast = 'per 1 SD', model = x$model[1], HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
             p = 2 * pnorm(-abs(b / se)), beta = b, se = se, n = sum(x$n), events = sum(x$events), n_genes = NA, rho_scS = NA, logrank_p = NA) }))
TAB <- rbind(TAB, pool); rownames(TAB) <- NULL
write.csv(TAB, file.path(outdir, 'F31_P4_SOX9_TGFb_module_survival_table.csv'), row.names = FALSE)
show <- function(x) { x$HR_CI <- sprintf('%.2f (%.2f-%.2f)', x$HR, x$lo, x$hi); x$p <- signif(x$p, 2); x$rho_scS <- round(x$rho_scS, 2); x[, c('set', 'cohort', 'contrast', 'model', 'HR_CI', 'p', 'n', 'events', 'n_genes', 'rho_scS')] }
cat('\n=== main module ===\n'); print(show(TAB[TAB$set == 'main', ]), row.names = FALSE)
cat('\n=== sensitivity, per 1 SD ===\n'); print(show(TAB[TAB$set != 'main' & TAB$contrast == 'per 1 SD', ]), row.names = FALSE)

## ---- figure: KM per cohort (median split) + forest of per-SD hazard ratios ----------------------
KM <- do.call(rbind, KM); KM$cohort <- factor(KM$cohort, COHORTS)
lab <- TAB[TAB$set == 'main' & TAB$contrast == 'high vs low' & TAB$model == 'alone', ]
lab$cohort <- factor(lab$cohort, COHORTS); lab$txt <- sprintf('HR %.2f (%.2f-%.2f)\nlog-rank p = %s\nn = %d', lab$HR, lab$lo, lab$hi, format.pval(lab$logrank_p, digits = 2, eps = 1e-4), lab$n)
pk <- ggplot(KM, aes(time, surv, colour = grp)) + geom_step(linewidth = .6) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = 1.05, vjust = 1.15, size = 2.3, lineheight = .95) +
  facet_wrap(~ cohort, nrow = 1, scales = 'free_x') + scale_colour_manual(values = c(high = '#B2182B', low = '#2166AC'), name = 'module score') +
  scale_y_continuous(limits = c(0, 1)) + labs(x = 'months', y = 'overall survival', title = sprintf('P4 SOX9 / TGF-beta module (top %d markers), median split', N_TOP)) +
  theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'), legend.position = 'bottom', plot.title = element_text(size = 8.5, face = 'bold'))
fo <- TAB[TAB$set == 'main' & TAB$contrast == 'per 1 SD', ]; fo$cohort <- factor(fo$cohort, rev(c(COHORTS, 'Pooled'))); fo$model <- factor(fo$model, c('alone', 'adjusted for scS'))
pf <- ggplot(fo, aes(HR, cohort, colour = model)) + geom_vline(xintercept = 1, linetype = 2, linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, position = position_dodge(.55), linewidth = .5) + geom_point(position = position_dodge(.55), size = 1.8) +
  geom_text(aes(x = max(fo$hi) * 1.08, label = sprintf('p = %s', format.pval(p, digits = 2, eps = 1e-4))), position = position_dodge(.55), hjust = 0, size = 2.2, show.legend = FALSE) +
  scale_x_log10(expand = expansion(mult = c(.05, .35))) + scale_colour_manual(values = c(alone = 'grey20', `adjusted for scS` = '#E08214'), name = NULL) +
  labs(x = 'hazard ratio per 1 SD of module score', y = NULL, title = 'Cox models') +
  theme_classic(base_size = 8) + theme(legend.position = 'bottom', plot.title = element_text(size = 8.5, face = 'bold'))
fig <- pk + pf + plot_layout(widths = c(3, 1.5))
ggsave(file.path(outdir, 'Plots', 'F31_P4_SOX9_TGFb_module_survival.pdf'), fig, width = 9.5, height = 3.2)
ggsave(file.path(outdir, 'Plots', 'F31_P4_SOX9_TGFb_module_survival.png'), fig, width = 9.5, height = 3.2, dpi = 150)
message('done')
