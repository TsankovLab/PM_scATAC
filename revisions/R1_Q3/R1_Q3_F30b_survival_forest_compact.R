###############################################################################
# R1_Q3 -- F30b: compact forest-plot version of F30.
#
# Same models as F30 (SOX9 / RUNX2 / SNAI2, TCGA + Bueno + MESOMICS), shown as hazard
# ratios instead of survival curves:
#   Alone            coxph(Surv ~ expression)
#   Adjusted for scS coxph(Surv ~ expression + scS)
# Two contrasts side by side:
# The plot shows the continuous contrast (expression z-scored within cohort, HR per 1 SD);
# the median-split contrast is still computed and kept in the table.
# Bueno survival is recorded in years and converted to months (does not affect HRs).
# Outputs: Plots/F30b_survival_forest_compact.pdf, F30b_survival_forest_table.csv
###############################################################################
set.seed(1234)
suppressMessages({ library(survival); library(ggplot2); library(dplyr)
                   library(RColorBrewer); library(paletteer); library(circlize) })

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

GENES   <- c('SOX9', 'RUNX2', 'SNAI2')
COHORTS <- c(tcga = 'TCGA', bueno = 'Bueno', mesomics = 'MESOMICS')
TIME_SCALE <- c(tcga = 1, bueno = 12, mesomics = 1)
ADJ_COL <- c(`Alone` = 'grey45', `Adjusted for scS` = '#B2182B')

expr_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meta_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))
build <- function(co) {
  e <- expr_l[[co]]; m <- meta_l[[co]]
  if (!identical(rownames(m), colnames(e))) m <- m[match(colnames(e), m$Sample), ]
  d <- as.data.frame(t(e[GENES, , drop = FALSE]))
  d$scS <- as.numeric(m$sarc_score); d$time <- m$census * TIME_SCALE[[co]]; d$status <- m$status
  d[complete.cases(d[, c(GENES, 'scS', 'time', 'status')]) & d$time > 0, ]
}
dfs <- setNames(lapply(names(COHORTS), build), names(COHORTS))

grab <- function(fit, term, ...) {
  s <- summary(fit); ci <- s$conf.int[term, ]
  data.frame(..., HR = ci['exp(coef)'], lower = ci['lower .95'], upper = ci['upper .95'],
             p = s$coefficients[term, 'Pr(>|z|)'], n = fit$n, events = fit$nevent, row.names = NULL)
}
RES <- do.call(rbind, lapply(names(COHORTS), function(co) do.call(rbind, lapply(GENES, function(g) {
  d <- dfs[[co]]
  d$z   <- as.numeric(scale(d[[g]]))
  d$grp <- factor(ifelse(d[[g]] > median(d[[g]]), 'high', 'low'), levels = c('low', 'high'))
  d$scSz <- as.numeric(scale(d$scS))
  rho <- cor(d[[g]], d$scS, method = 'spearman')
  rbind(
    grab(coxph(Surv(time, status) ~ z, d),          'z',       gene = g, cohort = co, contrast = 'per 1 SD',    model = 'Alone',            rho = rho),
    grab(coxph(Surv(time, status) ~ z + scSz, d),   'z',       gene = g, cohort = co, contrast = 'per 1 SD',    model = 'Adjusted for scS', rho = rho),
    grab(coxph(Surv(time, status) ~ grp, d),        'grphigh', gene = g, cohort = co, contrast = 'high vs low', model = 'Alone',            rho = rho),
    grab(coxph(Surv(time, status) ~ grp + scSz, d), 'grphigh', gene = g, cohort = co, contrast = 'high vs low', model = 'Adjusted for scS', rho = rho))
}))))
RES$gene     <- factor(RES$gene, levels = GENES)
RES$cohort   <- factor(COHORTS[RES$cohort], levels = unname(COHORTS))
RES$model    <- factor(RES$model, levels = c('Alone', 'Adjusted for scS'))
RES$contrast <- factor(RES$contrast, levels = c('per 1 SD', 'high vs low'))
RES$sig      <- RES$p < 0.05
write.csv(RES, file.path(outdir, 'F30b_survival_forest_table.csv'), row.names = FALSE)
cat('=== Cox HRs (bold rows = p < 0.05) ===\n')
print(RES %>% transmute(gene, cohort, contrast, model, HR = sprintf('%.2f (%.2f-%.2f)', HR, lower, upper),
                        p = signif(p, 2), rho_scS = round(rho, 2), n, events), row.names = FALSE)

## ---- compact forest ----------------------------------------------------------------------
## plotted: continuous contrast only, cohorts reading TCGA -> Bueno -> MESOMICS top to bottom
PLOT <- RES %>% filter(contrast == 'per 1 SD')
PLOT$cohort <- factor(as.character(PLOT$cohort), levels = rev(unname(COHORTS)))
pd <- position_dodge(width = .55)
p <- ggplot(PLOT, aes(HR, cohort, colour = model, shape = sig)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = 'grey55') +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = .5, position = pd) +
  geom_point(size = 2.1, position = pd) +
  scale_colour_manual(values = ADJ_COL, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = 'none') +
  scale_x_log10(breaks = c(.25, .5, 1, 2, 4)) +
  facet_grid(gene ~ .) +
  coord_cartesian(xlim = c(.25, 4)) +
  labs(x = 'Hazard ratio (95% CI)', y = NULL,
       title = 'Overall survival before and after scS adjustment',
       subtitle = 'Cox PH, HR per 1 SD of expression. Filled = p < 0.05;\nHR > 1 = worse survival with higher expression.') +
  gtheme_no_rot +
  theme(legend.position = 'top', strip.text.y = element_text(angle = 0),
        panel.grid.major.y = element_line(colour = 'grey92', linewidth = .3),
        plot.title = element_text(size = 9.5), plot.subtitle = element_text(size = 6.8, lineheight = 1.1))
ggsave(file.path(outdir, 'Plots', 'F30b_survival_forest_compact.pdf'), p, width = 5.4, height = 4)
cat('\nDONE\n')
