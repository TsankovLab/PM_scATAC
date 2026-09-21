# R1_Q3 – F27a: cBioPortal TCGA MESO — covariate association dotplot
#
# Pulls from cBioPortal public API (meso_tcga + meso_tcga_pan_can_atlas_2018):
#   - Clinical data (patient + sample level, both studies)
#   - Mutations (key mesothelioma driver genes)
#   - GISTIC copy-number alterations (discrete)
#   - RPPA protein expression
#
# Uses local bulk TCGA expression for SOX9 / RUNX2 / SNAI2 (already log2 RSEM).
# API responses are cached to cbio_cache/ to avoid re-downloading.
#
# F27a – Dotplot: all clinical variable associations (Spearman r / KW p-value)
#
# Output: R1_Q3/Plots/F27a_covariate_dotplot.pdf (+ F27_covariate_association_table.csv)

set.seed(1234)

projdir   <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir  <- file.path(projdir, 'bulkRNA_meso')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
cache_dir <- file.path(outdir, 'cbio_cache')
dir.create(cache_dir, showWarnings = FALSE)

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(httr)
library(jsonlite)
library(patchwork)
library(ggpubr)
library(tidyr)
library(dplyr)

anchor_genes <- c('SOX9', 'RUNX1', 'RUNX2', 'SNAI2')
CBIO_BASE    <- 'https://www.cbioportal.org/api'
STUDY_ID     <- 'meso_tcga'
STUDY_PCA    <- 'meso_tcga_pan_can_atlas_2018'
SAMPLE_LIST  <- 'meso_tcga_all'

# ── API helpers ───────────────────────────────────────────────────────────────
cbio_get <- function(path, query = list()) {
  resp <- GET(paste0(CBIO_BASE, path), query = query,
              add_headers(Accept = 'application/json'))
  stop_for_status(resp)
  fromJSON(content(resp, 'text', encoding = 'UTF-8'), flatten = TRUE)
}

cbio_post <- function(path, body) {
  resp <- POST(paste0(CBIO_BASE, path),
               body    = toJSON(body, auto_unbox = TRUE),
               encode  = 'raw',
               add_headers(`Content-Type` = 'application/json',
                           Accept          = 'application/json'))
  stop_for_status(resp)
  fromJSON(content(resp, 'text', encoding = 'UTF-8'), flatten = TRUE)
}

cache_rds <- function(file, expr) {
  path <- file.path(cache_dir, file)
  if (file.exists(path)) {
    message('  [cache] ', file)
    return(readRDS(path))
  }
  message('  [fetch] ', file)
  val <- force(expr)
  saveRDS(val, path)
  val
}

# ── Entrez gene IDs ───────────────────────────────────────────────────────────
driver_entrez <- c(
  BAP1=8314, NF2=4771, CDKN2A=1029, CDKN2B=1030, TP53=7157,
  SETD2=29072, LATS2=26524, LATS1=9113, STRN=6801,
  BRCA1=672, BRCA2=675, ATM=472,
  # Anchor genes themselves (for consistency check)
  SOX9=6662, RUNX2=860, SNAI2=6625
)

# Broad RPPA protein set (TCGA RPPA panel covers ~200 proteins)
rppa_entrez <- c(
  AKT1=207, AKT2=208, AKT3=10000, AR=367, ATM=472, BAD=572, BAX=581,
  BCL2=596, BCL2L1=598, BRAF=673, BRCA2=675, CASP3=836, CASP7=840, CASP8=841,
  CDH1=999, CDH2=1000, CDK1=983, CDKN1A=1026, CDKN1B=1027, CCND1=595, CCNE1=898,
  CTNNB1=1499, E2F1=1869, EGFR=1956, EIF4E=1977, EIF4EBP1=1978,
  ERBB2=2064, ERBB3=2065, ESR1=2099, FASN=2194, FN1=2335,
  FOXM1=2305, GSK3A=2931, GSK3B=2932, H2AFX=3014, HIF1A=3091,
  HSP90AA1=3320, HSPA1A=3303, IGF1R=3480, INSR=3643,
  JAK2=3717, JUN=3725, KDR=3791, KEAP1=9817, KIT=3815,
  MAPK1=5594, MAPK3=5595, MAPK14=1432, MAP2K1=5604,
  MCM2=4171, MDM2=4193, MET=4233, MKI67=4288,
  MLH1=4292, MSH2=4436, MSH6=2956, MYC=4609, MYH11=4629,
  NDRG1=10397, NFE2L2=4780, NOTCH1=4851,
  PARP1=142, PCNA=5111, PDCD4=27250, PDGFRB=5159, PDGFRA=5156,
  PIK3CA=5290, PIK3R1=5295, PMS2=5395,
  PRKCA=5578, PRKCB=5579, PRKCD=5580, PRKCZ=5590,
  PTEN=5728, PTK2=5747, PXN=5829, RAB25=57111, RAF1=5894,
  RB1=5925, RPS6=6194, RPS6KB1=6198, RPS6KB2=6199,
  SETDB1=9869, SHC1=6464, SMAD3=4088, SRC=6714,
  STAT3=6774, STAT5A=6776, STK11=6794, STMN1=3925,
  THBS1=7057, TOP2A=7153, TP53=7157, TSC1=7248, TSC2=7249,
  TTK=7272, VIM=7431, YAP1=10413, YWHAB=7529,
  ACTA2=59, ACTB=60, CDH3=1001, TWIST1=7291, ZEB1=6935,
  PTPN11=5781, RET=5979, FGFR1=2260
)

##############################################################################
# 1. Load local TCGA expression (anchor genes)
##############################################################################
message('\n=== Loading local TCGA expression ===')
tcga_expr <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))[['tcga']]
tcga_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))[['tcga']]

# Convert local IDs (dots) to cBioPortal format (dashes)
local_ids  <- colnames(tcga_expr)
cbio_ids   <- gsub('\\.', '-', local_ids)
message('Samples: ', length(cbio_ids))

anc_ok <- anchor_genes[anchor_genes %in% rownames(tcga_expr)]
anc_df <- as.data.frame(t(tcga_expr[anc_ok, , drop = FALSE]))
anc_df$sampleId <- cbio_ids
message('Anchor genes present: ', paste(anc_ok, collapse = ', '))

##############################################################################
# 2. Pull clinical data from cBioPortal (both studies)
##############################################################################
message('\n=== Fetching clinical data ===')

clin_patient <- cache_rds('clin_patient.rds', {
  cbio_get(paste0('/studies/', STUDY_ID, '/clinical-data'),
           query = list(clinicalDataType = 'PATIENT', pageSize = 5000))
})

clin_sample <- cache_rds('clin_sample.rds', {
  cbio_get(paste0('/studies/', STUDY_ID, '/clinical-data'),
           query = list(clinicalDataType = 'SAMPLE', pageSize = 5000))
})

# PanCancer Atlas extras (TMB, aneuploidy, MSI, genetic ancestry)
clin_pca <- cache_rds('clin_pca.rds', {
  cbio_get(paste0('/studies/', STUDY_PCA, '/clinical-data'),
           query = list(clinicalDataType = 'PATIENT', pageSize = 5000))
})

# Pivot to wide — sample level keeps both sampleId and patientId as id cols
clin_s_wide <- clin_sample %>%
  select(sampleId, patientId, clinicalAttributeId, value) %>%
  distinct(sampleId, clinicalAttributeId, .keep_all = TRUE) %>%
  pivot_wider(id_cols = c(sampleId, patientId),
              names_from = clinicalAttributeId, values_from = value,
              values_fn = first)

clin_p_wide <- clin_patient %>%
  select(patientId, clinicalAttributeId, value) %>%
  distinct(patientId, clinicalAttributeId, .keep_all = TRUE) %>%
  pivot_wider(id_cols = patientId,
              names_from = clinicalAttributeId, values_from = value,
              values_fn = first)

pca_attrs <- c('TMB_NONSYNONYMOUS', 'ANEUPLOIDY_SCORE', 'MSI_SCORE_MANTIS',
               'MSI_SENSOR_SCORE', 'TBL_SCORE', 'GENETIC_ANCESTRY_LABEL',
               'FRACTION_GENOME_ALTERED', 'MUTATION_COUNT')
clin_pca_wide <- clin_pca %>%
  filter(clinicalAttributeId %in% pca_attrs) %>%
  select(patientId, clinicalAttributeId, value) %>%
  distinct(patientId, clinicalAttributeId, .keep_all = TRUE) %>%
  pivot_wider(id_cols = patientId,
              names_from = clinicalAttributeId, values_from = value,
              values_fn = first)

# Merge: suffix _pat for patient-level attrs, _pca for PanCancer Atlas attrs
clin_all <- clin_s_wide %>%
  left_join(clin_p_wide,   by = 'patientId', suffix = c('', '_pat')) %>%
  left_join(clin_pca_wide, by = 'patientId', suffix = c('', '_pca'))

message('Clinical data: ', nrow(clin_all), ' samples x ', ncol(clin_all), ' attributes')

##############################################################################
# 3. Pull mutations for driver genes
##############################################################################
message('\n=== Fetching mutations ===')

muts_raw <- cache_rds('mutations.rds', {
  cbio_post(
    paste0('/molecular-profiles/', STUDY_ID, '_mutations/mutations/fetch'),
    list(entrezGeneIds = unname(driver_entrez),
         sampleListId  = SAMPLE_LIST)
  )
})

entrez_to_sym <- setNames(names(driver_entrez), driver_entrez)

if (length(muts_raw) > 0 && is.data.frame(muts_raw)) {
  mut_binary <- muts_raw %>%
    mutate(gene = entrez_to_sym[as.character(entrezGeneId)]) %>%
    filter(!is.na(gene)) %>%
    select(sampleId, gene) %>%
    distinct() %>%
    mutate(mutated = 1L) %>%
    pivot_wider(names_from = gene, values_from = mutated,
                values_fill = 0L, names_prefix = 'MUT_')
  message('Mutation data: ', nrow(muts_raw), ' events across ',
          ncol(mut_binary) - 1, ' genes')
} else {
  message('No mutation data returned')
  mut_binary <- data.frame(sampleId = cbio_ids)
}

##############################################################################
# 4. Pull GISTIC copy-number alterations
##############################################################################
message('\n=== Fetching GISTIC CNA ===')

cna_raw <- cache_rds('gistic_cna.rds', {
  cbio_post(
    paste0('/molecular-profiles/', STUDY_ID, '_gistic/discrete-copy-number/fetch'),
    list(entrezGeneIds = unname(driver_entrez),
         sampleListId  = SAMPLE_LIST,
         discreteCopyNumberEventType = 'ALL')
  )
})

if (length(cna_raw) > 0 && is.data.frame(cna_raw)) {
  cna_wide <- cna_raw %>%
    mutate(gene = entrez_to_sym[as.character(entrezGeneId)]) %>%
    filter(!is.na(gene)) %>%
    select(sampleId, gene, value = alteration) %>%
    pivot_wider(names_from = gene, values_from = value,
                names_prefix = 'CNA_')
  message('CNA data: ', ncol(cna_wide) - 1, ' genes')
} else {
  message('No CNA data returned')
  cna_wide <- data.frame(sampleId = cbio_ids)
}

##############################################################################
# 5. Pull RPPA protein expression
##############################################################################
message('\n=== Fetching RPPA ===')

rppa_raw <- cache_rds('rppa.rds', {
  cbio_post(
    paste0('/molecular-profiles/', STUDY_ID, '_rppa/molecular-data/fetch'),
    list(entrezGeneIds = unname(rppa_entrez),
         sampleListId  = SAMPLE_LIST)
  )
})

entrez_to_sym_rppa <- setNames(names(rppa_entrez), rppa_entrez)

if (length(rppa_raw) > 0 && is.data.frame(rppa_raw)) {
  rppa_wide <- rppa_raw %>%
    mutate(gene = entrez_to_sym_rppa[as.character(entrezGeneId)]) %>%
    filter(!is.na(gene)) %>%
    select(sampleId, gene, value) %>%
    group_by(sampleId, gene) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = gene, values_from = value,
                names_prefix = 'RPPA_')
  message('RPPA data: ', ncol(rppa_wide) - 1, ' proteins')
} else {
  message('No RPPA data returned')
  rppa_wide <- data.frame(sampleId = cbio_ids)
}

##############################################################################
# 6. Merge all data
##############################################################################
message('\n=== Merging datasets ===')

master <- anc_df %>%
  left_join(clin_all,   by = 'sampleId') %>%
  left_join(mut_binary, by = 'sampleId') %>%
  left_join(cna_wide,   by = 'sampleId') %>%
  left_join(rppa_wide,  by = 'sampleId')

message('Master data frame: ', nrow(master), ' samples x ', ncol(master), ' columns')
saveRDS(master, file.path(cache_dir, 'master_df.rds'))

##############################################################################
# Helper: run association tests for one anchor gene vs all covariates
##############################################################################
test_covariate <- function(master, anc) {
  anc_v <- as.numeric(master[[anc]])
  if (all(is.na(anc_v))) return(NULL)

  cols <- setdiff(colnames(master), c(anchor_genes, 'sampleId', 'patientId',
                                       'patientId_pat'))

  do.call(rbind, lapply(cols, function(col) {
    x <- master[[col]]
    n_nonNA <- sum(!is.na(x) & !is.na(anc_v))
    if (n_nonNA < 10) return(NULL)

    # Attempt numeric conversion
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    is_num <- sum(!is.na(x_num)) >= 0.8 * n_nonNA

    tryCatch({
      if (is_num) {
        # Spearman correlation
        ct <- cor.test(anc_v, x_num, method = 'spearman',
                       use = 'complete.obs', exact = FALSE)
        data.frame(anchor = anc, covariate = col, test = 'spearman',
                   stat = ct$estimate, pval = ct$p.value, n = n_nonNA,
                   stringsAsFactors = FALSE)
      } else {
        # Kruskal-Wallis
        grp <- as.character(x)
        grp <- grp[!is.na(grp) & !is.na(anc_v)]
        a_v <- anc_v[!is.na(x) & !is.na(anc_v)]
        if (length(unique(grp)) < 2 || length(unique(grp)) > 20) return(NULL)
        kt <- kruskal.test(a_v ~ factor(grp))
        # Effect size: eta-squared approximation
        n_tot <- length(a_v)
        eta2  <- (kt$statistic - length(unique(grp)) + 1) /
                 (n_tot - length(unique(grp)))
        eta2  <- max(0, eta2)
        data.frame(anchor = anc, covariate = col, test = 'kruskal',
                   stat = unname(eta2), pval = kt$p.value, n = n_tot,
                   stringsAsFactors = FALSE)
      }
    }, error = function(e) NULL)
  }))
}

message('\n=== Running association tests ===')
assoc_all <- do.call(rbind, lapply(anc_ok, test_covariate, master = master))
assoc_all <- assoc_all[!is.na(assoc_all$pval), ]
assoc_all$padj <- p.adjust(assoc_all$pval, method = 'BH')
assoc_all$neglog10p <- -log10(assoc_all$pval)

# Categorize covariate
assoc_all$category <- case_when(
  grepl('^MUT_',  assoc_all$covariate) ~ 'Mutation',
  grepl('^CNA_',  assoc_all$covariate) ~ 'CNA',
  grepl('^RPPA_', assoc_all$covariate) ~ 'RPPA',
  TRUE                                  ~ 'Clinical'
)

write.csv(assoc_all, file.path(outdir, 'F27_covariate_association_table.csv'),
          row.names = FALSE)
message('Association table: ', nrow(assoc_all), ' tests')

##############################################################################
# F27a – Dotplot: top significant associations per category
##############################################################################
message('\n=== F27a: dotplot ===')

top_assoc <- assoc_all %>%
  group_by(anchor, category) %>%
  slice_min(order_by = pval, n = 15) %>%
  ungroup() %>%
  arrange(anchor, category, pval)

top_assoc$covariate_clean <- gsub('^(MUT_|CNA_|RPPA_)', '', top_assoc$covariate)
top_assoc$covariate_clean <- gsub('_pat$|_pca$', '', top_assoc$covariate_clean)
top_assoc$sig <- ifelse(top_assoc$padj < 0.05, 'FDR<0.05',
                 ifelse(top_assoc$pval < 0.05, 'p<0.05', 'ns'))

cat_pal <- c(Clinical = '#4393C3', Mutation = '#D6604D',
             CNA = '#F4A582', RPPA = '#74C476')

f27a_list <- lapply(anc_ok, function(anc) {
  sub <- top_assoc[top_assoc$anchor == anc, ]
  sub <- sub %>%
    arrange(category, pval) %>%
    mutate(label = factor(paste0('[', category, '] ', covariate_clean),
                          levels = rev(unique(paste0('[', category, '] ', covariate_clean)))))

  ggplot(sub, aes(x = neglog10p, y = label, color = category, shape = sig,
                  size = abs(stat))) +
    geom_point(alpha = 0.85) +
    geom_vline(xintercept = -log10(0.05), linetype = 'dashed',
               color = 'grey50', linewidth = 0.4) +
    scale_color_manual(values = cat_pal) +
    scale_shape_manual(values = c('FDR<0.05' = 16, 'p<0.05' = 17, 'ns' = 1)) +
    scale_size_continuous(range = c(1.5, 5), name = '|effect|') +
    labs(x = expression(-log[10](p)), y = NULL,
         title = paste0(anc, ' — covariate associations'),
         color = 'Category', shape = 'Significance') +
    gtheme +
    theme(axis.text.y = element_text(size = 7))
})

pdf(file.path(outdir, 'Plots', 'F27a_covariate_dotplot.pdf'),
    width = 7 * length(f27a_list), height = 10)
print(wrap_plots(f27a_list, ncol = length(f27a_list)))
dev.off()
message('F27a done.')
