###############################################################################
# R1_Q3 -- F30: survival curves for SOX9, RUNX2 and SNAI2 in the three bulk RNA
#          cohorts (TCGA, Bueno, MESOMICS), adjusted for the scS (sarcomatoid) score.
#
# A Kaplan-Meier curve cannot itself be "adjusted" for a covariate, so each panel shows
# both of the following for a median split of the gene's expression:
#   thin dashed  unadjusted Kaplan-Meier, with the log-rank p
#   thick solid  scS-ADJUSTED survival curves from coxph(Surv ~ expr_group + scS), by direct
#                standardisation (g-formula): for each group level the model's predicted
#                survival is evaluated for every patient at their OWN scS and averaged, so the
#                two curves describe the same cohort differing only in expression group.
#                Ribbon = 95% CI from BOOT patient-level bootstrap refits.
# Annotation per panel: adjusted HR (high vs low) with 95% CI and Wald p, plus the
#   unadjusted log-rank p.  The table also carries the continuous per-SD HRs (as F29b).
# Time: Bueno is recorded in years and is converted to months; the other two are months.
#   Cohorts are plotted on their own follow-up range.
# scS = sarc_score from bulkRNA_meso/bulk_RNA_studies_metadata.rds (as F29/F29b).
#
# Outputs: Plots/F30_survival_curves_scSadjusted.pdf         (3 genes x 3 cohorts)
#          Plots/F30_survival_curves_scSadjusted_KMonly.pdf  (unadjusted KM alone)
#          F30_survival_curves_stats.csv
###############################################################################
set.seed(1234)
suppressMessages({ library(survival); library(ggplot2); library(dplyr)
                   library(RColorBrewer); library(paletteer); library(circlize) })  # palettes.R needs these

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
dir.create(file.path(outdir, 'Plots'), showWarnings = FALSE)

GENES   <- c('SOX9', 'RUNX2', 'SNAI2')
COHORTS <- c(tcga = 'TCGA', bueno = 'Bueno', mesomics = 'MESOMICS')
TIME_SCALE <- c(tcga = 1, bueno = 12, mesomics = 1)   # -> months (Bueno is in years)
BOOT    <- 500
GRP_COL <- c(high = '#B2182B', low = '#2166AC')

expr_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meta_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

## ---- one aligned data frame per cohort -------------------------------------------------
build <- function(co) {
  e <- expr_l[[co]]; m <- meta_l[[co]]
  if (!identical(rownames(m), colnames(e))) {           # mesomics: align through Sample
    stopifnot('Sample' %in% colnames(m), all(colnames(e) %in% m$Sample))
    m <- m[match(colnames(e), m$Sample), ]
  }
  d <- as.data.frame(t(e[GENES, , drop = FALSE]))
  d$sample_id <- colnames(e); d$scS <- m$sarc_score
  d$time <- m$census * TIME_SCALE[[co]]; d$status <- m$status
  d$subtype <- as.character(m$subtype); d$cohort <- co
  d[complete.cases(d[, c(GENES, 'scS', 'time', 'status')]) & d$time > 0, ]
}
dfs <- setNames(lapply(names(COHORTS), build), names(COHORTS))
for (co in names(dfs)) message(COHORTS[[co]], ': n = ', nrow(dfs[[co]]),
                               ', events = ', sum(dfs[[co]]$status),
                               ', follow-up to ', round(max(dfs[[co]]$time), 1), ' months')

## ---- direct adjusted (g-formula) survival curves ----------------------------------------
## mean over patients of S(t | group = g, scS = scS_i); with no scS term this reduces to the
## model-based curve for the group, which is why the same helper serves both fits.
adj_curves <- function(d, grid, adjust = TRUE) {
  f <- if (adjust) Surv(time, status) ~ grp + scS else Surv(time, status) ~ grp
  fit <- coxph(f, data = d)
  out <- lapply(levels(d$grp), function(g) {
    nd <- d; nd$grp <- factor(g, levels = levels(d$grp))
    sf <- survfit(fit, newdata = nd)
    S  <- summary(sf, times = grid, extend = TRUE)$surv
    if (is.null(dim(S))) S <- matrix(S, ncol = nrow(nd))
    data.frame(grp = g, time = grid, surv = rowMeans(S))
  })
  list(curves = do.call(rbind, out), fit = fit)
}

## ---- per gene x cohort -------------------------------------------------------------------
CURV <- list(); KM <- list(); STAT <- list()
for (co in names(COHORTS)) {
  d0 <- dfs[[co]]
  grid <- seq(0, max(d0$time), length.out = 200)
  for (g in GENES) {
    d <- d0; d$expr <- d[[g]]
    d$grp <- factor(ifelse(d$expr > median(d$expr), 'high', 'low'), levels = c('low', 'high'))
    d$scS <- as.numeric(scale(d$scS))

    ## unadjusted KM (step curve) and log-rank
    km <- survfit(Surv(time, status) ~ grp, data = d)
    st <- rep(sub('grp=', '', names(km$strata)), km$strata)
    km_df <- data.frame(grp = st, time = km$time, surv = km$surv)
    km_df <- rbind(data.frame(grp = levels(d$grp), time = 0, surv = 1), km_df)
    lr_p <- 1 - pchisq(survdiff(Surv(time, status) ~ grp, data = d)$chisq, 1)

    ## adjusted curves + bootstrap CI
    A  <- adj_curves(d, grid, adjust = TRUE)
    bs <- replicate(BOOT, {
      i <- sample.int(nrow(d), replace = TRUE)
      tryCatch(adj_curves(d[i, ], grid, adjust = TRUE)$curves$surv, error = function(e) rep(NA_real_, 2 * length(grid)))
    })
    A$curves$lo <- apply(bs, 1, quantile, .025, na.rm = TRUE)
    A$curves$hi <- apply(bs, 1, quantile, .975, na.rm = TRUE)

    s  <- summary(A$fit); ci <- s$conf.int['grphigh', ]
    ## continuous per-SD HR, adjusted and alone (comparable to F29b)
    dz <- d; dz$expr_z <- as.numeric(scale(dz$expr))
    c_adj <- summary(coxph(Surv(time, status) ~ expr_z + scS, data = dz))
    c_alo <- summary(coxph(Surv(time, status) ~ expr_z, data = dz))

    lab <- sprintf('adj HR %.2f (%.2f-%.2f)\nadj p = %s\nlog-rank p = %s',
                   ci['exp(coef)'], ci['lower .95'], ci['upper .95'],
                   format.pval(s$coefficients['grphigh', 'Pr(>|z|)'], digits = 2, eps = 1e-4),
                   format.pval(lr_p, digits = 2, eps = 1e-4))
    key <- data.frame(gene = g, cohort = co)
    CURV[[length(CURV) + 1]] <- cbind(key, A$curves)
    KM[[length(KM) + 1]]     <- cbind(key, km_df)
    STAT[[length(STAT) + 1]] <- data.frame(key,
      n = nrow(d), events = sum(d$status), n_high = sum(d$grp == 'high'),
      median_cut = median(d$expr), label = lab, logrank_p = lr_p,
      HR_group_adj = ci['exp(coef)'], HR_group_adj_lo = ci['lower .95'], HR_group_adj_hi = ci['upper .95'],
      p_group_adj = s$coefficients['grphigh', 'Pr(>|z|)'],
      HR_perSD_alone = c_alo$conf.int['expr_z', 'exp(coef)'], p_perSD_alone = c_alo$coefficients['expr_z', 'Pr(>|z|)'],
      HR_perSD_adj = c_adj$conf.int['expr_z', 'exp(coef)'], p_perSD_adj = c_adj$coefficients['expr_z', 'Pr(>|z|)'],
      HR_scS_adj = c_adj$conf.int['scS', 'exp(coef)'], p_scS_adj = c_adj$coefficients['scS', 'Pr(>|z|)'],
      r_expr_scS = cor(dz$expr, dz$scS, method = 'spearman'),
      row.names = NULL)
  }
}
CURV <- do.call(rbind, CURV); KM <- do.call(rbind, KM); STAT <- do.call(rbind, STAT)
for (D in c('CURV', 'KM', 'STAT')) {
  d <- get(D); d$gene <- factor(d$gene, levels = GENES)
  d$cohort <- factor(COHORTS[as.character(d$cohort)], levels = unname(COHORTS)); assign(D, d)
}
write.csv(STAT[, setdiff(names(STAT), 'label')], file.path(outdir, 'F30_survival_curves_stats.csv'), row.names = FALSE)
cat('\n=== scS-adjusted survival, high vs low expression (median split) ===\n')
print(STAT %>% transmute(gene, cohort, n, events,
        `HR high/low (adj)` = sprintf('%.2f (%.2f-%.2f)', HR_group_adj, HR_group_adj_lo, HR_group_adj_hi),
        `p adj` = signif(p_group_adj, 2), `log-rank p` = signif(logrank_p, 2),
        `HR/SD alone` = round(HR_perSD_alone, 2), `p alone` = signif(p_perSD_alone, 2),
        `HR/SD adj` = round(HR_perSD_adj, 2), `p adj (SD)` = signif(p_perSD_adj, 2),
        `rho expr~scS` = round(r_expr_scS, 2)), row.names = FALSE)

## ---- figure ------------------------------------------------------------------------------
lab_df <- STAT
p <- ggplot() +
  geom_ribbon(data = CURV, aes(time, ymin = lo, ymax = hi, fill = grp), alpha = .15, colour = NA) +
  geom_step(data = KM, aes(time, surv, colour = grp), linetype = '22', linewidth = .45, alpha = .9) +
  geom_line(data = CURV, aes(time, surv, colour = grp), linewidth = .8) +
  geom_text(data = lab_df, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.15,
            size = 2.3, lineheight = .95, colour = 'grey20') +
  geom_text(data = lab_df, aes(x = 0, y = 0, label = sprintf('n = %d (%d events)', n, events)),
            hjust = -0.05, vjust = -0.4, size = 2.2, colour = 'grey45') +
  scale_colour_manual(values = GRP_COL, breaks = c('high', 'low'),
                      labels = c('high (> median)', 'low'), name = NULL) +
  scale_fill_manual(values = GRP_COL, guide = 'none') +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(.02, .02))) +
  facet_grid(gene ~ cohort, scales = 'free_x') +
  labs(x = 'Months', y = 'Overall survival',
       title = 'SOX9, RUNX2 and SNAI2: survival adjusted for the scS (sarcomatoid) score',
       subtitle = paste('Solid = scS-adjusted curves from Cox, averaged over the cohort\'s own scS distribution',
                        '(ribbon: 95% bootstrap CI). Dashed = unadjusted Kaplan-Meier. Median split of expression.')) +
  gtheme_no_rot +
  theme(strip.text.y = element_text(angle = 0), legend.position = 'top',
        plot.title = element_text(size = 11), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(outdir, 'Plots', 'F30_survival_curves_scSadjusted.pdf'), p, width = 9.5, height = 8)

p2 <- ggplot(KM, aes(time, surv, colour = grp)) +
  geom_step(linewidth = .7) +
  geom_text(data = lab_df, aes(x = Inf, y = Inf, label = sprintf('log-rank p = %s', format.pval(logrank_p, digits = 2, eps = 1e-4))),
            hjust = 1.05, vjust = 1.4, size = 2.4, colour = 'grey20', inherit.aes = FALSE) +
  scale_colour_manual(values = GRP_COL, breaks = c('high', 'low'), labels = c('high (> median)', 'low'), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) + facet_grid(gene ~ cohort, scales = 'free_x') +
  labs(x = 'Months', y = 'Overall survival', title = 'Unadjusted Kaplan-Meier (median split of expression)') +
  gtheme_no_rot + theme(strip.text.y = element_text(angle = 0), legend.position = 'top')
ggsave(file.path(outdir, 'Plots', 'F30_survival_curves_scSadjusted_KMonly.pdf'), p2, width = 9.5, height = 8)
cat('\nDONE\n')
