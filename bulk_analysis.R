conda activate meso_scatac

R

# Load packages
packages = c(
  'Seurat',
  'ggpubr',
  'org.Hs.eg.db',
  'TxDb.Hsapiens.UCSC.hg38.knownGene',
  'EnsDb.Hsapiens.v86',
  'gplots',
  'ggrepel',
  'patchwork',
  'ComplexHeatmap',
  'RColorBrewer',
  'ggrepel',
  'gplots',
  'EnhancedVolcano',
  'rstatix',
  'tidyr',
  'corrplot',
  'scales',
  'ggExtra',
  'limma',
  'ape',
  'dendextend',
  'survival',
  'survminer')
lapply(packages, require, character.only = TRUE)

set.seed(1234)

# Set directory
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd(projdir)

source ('../git_repo/utils/ggplot_aestetics.R')
bulk_palette = setNames (hue_pal()(4)[c(4,3,2,1)],c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
bulk_palette = setNames (as.character(paletteer::paletteer_d("rockthemes::heep")), c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid'))
bueno_palette = setNames (paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7)], c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
bulk_palette = setNames (as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7,3)]), c('Epithelioid','Biphasic-E','Biphasic-S','Sarcomatoid','Biphasic'))


#### Import TCGA ####
if (!file.exists ('meso_tcga_RSEM_cleaned.rds'))
  {
  # load ('/ahg/regevdata/projects/lungCancerBueno/Results/TCGAsubtypes/Cancer_All/Results/gcdata.merge.Rda')
  # gcdata.merge_meso = subset (gcdata.merge, meso_histology %in% c('MESO.E','MESO.NOS','MESO.S','MESO.BP')) 
  # tcga_meta = gcdata.merge_meso$meso_histology
  # all (colnames (gcdata.merge_meso) == names (tcga_meta))
  
  
  ## Load TCGA data ###
  tcga_mat = read.table ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/TCGA/MESO.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep='\t', header=T)
  rownames (tcga_mat) = tcga_mat[,1]
  tcga_mat = tcga_mat[,grep ('raw_count', tcga_mat[1,])]
  tcga_mat = tcga_mat[-1,]
  tcga_mat = tcga_mat[!grepl ('\\?',rownames (tcga_mat)),]
  symbols = sapply (rownames(tcga_mat), function(x) unlist(strsplit(x, '\\|'))[1])
  tcga_mat = tcga_mat[!duplicated (symbols),]
  rownames(tcga_mat)  = sapply (rownames(tcga_mat), function(x) unlist(strsplit(x, '\\|'))[1])
  symbol = rownames (tcga_mat)
  tcga_mat = apply(tcga_mat, 2, function(x) as.numeric(x))
  rownames (tcga_mat) = symbol
  
  # link clinical data from cbioPortal ####
  load ('/ahg/regevdata/projects/lungCancerBueno/Results/10x_MesoSCInt_200930/MESO_TCGA_Bueno/TCGA_MESO_cbio.Rda')
  tcga_meta = cbio
  tcga_meta$DEATH_EVENT = ifelse (tcga_meta$OS_STATUS == '1:DECEASED',1,0)
  idx = sapply (rownames(tcga_meta), function(x) grep (x, colnames(tcga_mat)))
  tcga_mat = tcga_mat[,idx]
  colnames (tcga_mat) = rownames(tcga_meta)
  tcga_mat = tcga_mat[,tcga_meta$TUMOR_TYPE != 'Diffuse Malignant Mesothelioma (NOS)']
  tcga_meta = tcga_meta[tcga_meta$TUMOR_TYPE != 'Diffuse Malignant Mesothelioma (NOS)',]
  
  saveRDS (tcga_mat, 'meso_tcga_RSEM_cleaned.rds')
  saveRDS (tcga_meta, 'meso_tcga_meta_cleaned.rds')
  
  # Check if it is seq-depth normalized
  pdf ('tcga_raw.pdf')
  print(boxplot (log2(as.matrix(tcga_mat))))
  dev.off()
  } else {
tcga_mat = readRDS ('meso_tcga_RSEM_cleaned.rds')
tcga_meta = readRDS ('meso_tcga_meta_cleaned.rds')  
}

#### Import Bueno ####
if (!file.exists ('bueno_rpkm_cleaned.rds'))
  {
  # Import bueno from my folder
  bueno_mat = read.table ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/bueno_data/bueno_meso_rpkm_rnaseq.txt', header=T, sep= '\t')
  bueno_mat = bueno_mat[!duplicated(bueno_mat[,1]),]
  rownames (bueno_mat) = bueno_mat[,1]
  bueno_mat = bueno_mat[,-1]
  
  # Get clinical data
  bueno_meta = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/ClincalInfo_S2_2.csv')
  rownames (bueno_meta) = bueno_meta$Tumor.ID
  
  bueno_meta = bueno_meta[bueno_meta$ConsensusCluster != '',]
  bueno_meta2 = setNames (bueno_meta$ConsensusCluster,rownames(bueno_meta))
  bueno_mat = bueno_mat[,names(bueno_meta2)]
  
  saveRDS (bueno_meta, 'bueno_meta_cli.rds')
  saveRDS (bueno_meta2, 'bueno_meta_cleaned.rds')
  saveRDS (bueno_mat, 'bueno_rpkm_cleaned.rds')
    
    
  # Check if it is seq-depth normalized
  pdf ('bueno_rpkm.pdf')
  print (boxplot (log2(bueno_mat)))
  dev.off()
  } else {
  bueno_meta = readRDS ('bueno_meta_cli.rds')    
  bueno_meta2 = readRDS ('bueno_meta_cleaned.rds')
  bueno_mat = readRDS ('bueno_rpkm_cleaned.rds')
  }


#### Import MESOMICS ####
if (!file.exists ('mesomics_data.rds'))
  {
  msm_mat = read.table ('/sc/arion/projects/Tsankov_Normal_Lung/ICA_Lung/Wooseung/Mesothelioma/Data/Expression/TPM_COnverted_GeneCount.tsv')
  #msm = log2(msm+1)
  #msm = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/Mesothelioma/Data/Expression/VStexpr.tsv', header=T)
  msm_meta = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/Mesothelioma/Data/Mesothelioma_SampleInfo.txt', header=T)
  msm_mat = msm_mat[,colnames (msm_mat) %in% msm_meta$Sample]
  msm_meta = msm_meta[msm_meta$Sample %in% colnames (msm),]
  msm_mat = msm_mat[, msm_meta$Sample]
  all (colnames (msm_mat) == msm_meta$Sample)
  saveRDS (msm_mat,'mesomics_data.rds')
  saveRDS (msm_meta,'Mesothelioma_SampleInfo.txt')
  } else {
  msm_mat = readRDS ('mesomics_data.rds')
  msm_meta = readRDS ('Mesothelioma_SampleInfo.txt')
  msm_meta2 = read.csv ('mesomics_SamplesOverview.csv')
  }

msm_meta2 = msm_meta2[match(colnames(msm_mat), msm_meta2$Sample), ]
msm_meta2 = msm_meta2[, 1:96]

### Make list including all the bulk data ####
meso_bulk_l = list (
  bueno = log2(bueno_mat[,rownames(bueno_meta)] + 1),
  tcga = log2 (tcga_mat + 1),
  mesomics = log2 (msm_mat+1))

## combine metadata from all bulk studies ####
meso_bulk_meta_l = list (
  bueno = bueno_meta,
  tcga = tcga_meta,
  mesomics = msm_meta2)

# Remove NOS subtype in TCGA
meso_bulk_meta_l[['tcga']] = meso_bulk_meta_l[['tcga']][meso_bulk_meta_l[['tcga']]$TUMOR_TYPE != 'Diffuse Malignant Mesothelioma (NOS)' ,]
meso_bulk_l[['tcga']] = meso_bulk_l[['tcga']][,rownames(meso_bulk_meta_l[['tcga']])]

# Harmonize metadata ####
meso_bulk_meta_l[['tcga']]$TUMOR_TYPE[meso_bulk_meta_l[['tcga']]$TUMOR_TYPE == 'Biphasic Mesothelioma'] = 'Biphasic'
meso_bulk_meta_l[['tcga']]$TUMOR_TYPE[meso_bulk_meta_l[['tcga']]$TUMOR_TYPE == 'Epithelioid Mesothelioma'] = 'Epithelioid'
meso_bulk_meta_l[['tcga']]$TUMOR_TYPE[meso_bulk_meta_l[['tcga']]$TUMOR_TYPE == 'Sarcomatoid Mesothelioma'] = 'Sarcomatoid'

meso_bulk_meta_l[['mesomics']]$Type[meso_bulk_meta_l[['mesomics']]$Type == 'MMS'] = 'Sarcomatoid'
meso_bulk_meta_l[['mesomics']]$Type[meso_bulk_meta_l[['mesomics']]$Type == 'MME'] = 'Epithelioid'
meso_bulk_meta_l[['mesomics']]$Type[meso_bulk_meta_l[['mesomics']]$Type == 'MMB'] = 'Biphasic'

meso_bulk_meta_l[['bueno']]$ConsensusCluster = factor (meso_bulk_meta_l[['bueno']]$ConsensusCluster, levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['tcga']]$TUMOR_TYPE = factor (meso_bulk_meta_l[['tcga']]$TUMOR_TYPE, levels = c('Sarcomatoid','Biphasic','Epithelioid'))
meso_bulk_meta_l[['mesomics']]$Type = factor (meso_bulk_meta_l[['mesomics']]$Type, levels = c('Sarcomatoid','Biphasic','Epithelioid'))

colnames (meso_bulk_meta_l[['bueno']])[colnames (meso_bulk_meta_l[['bueno']]) == 'ConsensusCluster'] = 'subtype'
colnames (meso_bulk_meta_l[['tcga']])[colnames (meso_bulk_meta_l[['tcga']]) == 'TUMOR_TYPE'] = 'subtype'
colnames (meso_bulk_meta_l[['mesomics']])[colnames (meso_bulk_meta_l[['mesomics']]) == 'Type'] = 'subtype'

colnames(meso_bulk_meta_l[['bueno']])[colnames (meso_bulk_meta_l[['bueno']]) == 'Survival.from.surgery..years.'] = 'census'
colnames(meso_bulk_meta_l[['bueno']])[colnames (meso_bulk_meta_l[['bueno']]) == 'Status'] = 'status'
colnames (meso_bulk_meta_l[['tcga']])[colnames (meso_bulk_meta_l[['tcga']]) == 'DEATH_EVENT'] = 'status'
colnames (meso_bulk_meta_l[['tcga']])[colnames (meso_bulk_meta_l[['tcga']]) == 'OS_MONTHS'] = 'census'
colnames (meso_bulk_meta_l[['mesomics']])[colnames (meso_bulk_meta_l[['mesomics']]) == 'Survival.Censor'] = 'status'
colnames (meso_bulk_meta_l[['mesomics']])[colnames (meso_bulk_meta_l[['mesomics']]) == 'Survival.Time'] = 'census'

#meso_bulk_meta_l[['tcga']]$status = ifelse (meso_bulk_meta_l[['tcga']]$status == '1:DECEASED',1,0)
meso_bulk_meta_l[['mesomics']]$status = ifelse (meso_bulk_meta_l[['mesomics']]$status == 'dead', 1,0)
meso_bulk_meta_l[['bueno']]$status = ifelse (meso_bulk_meta_l[['bueno']]$status == 'd', 1,0)

# Add scS score as meta variable ####
nfeat = 5000
k_selection = 25
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))
sarc_score = head (cnmf_spectra_unique[['cNMF20']],20)

meso_bulk_meta_l = lapply (names(meso_bulk_meta_l), function(x) {
  tmp = colMeans(na.omit(meso_bulk_l[[x]][rownames(meso_bulk_l[[x]]) %in% sarc_score,]))
  meso_bulk_meta_l[[x]]$sarc_score = tmp
  meso_bulk_meta_l[[x]]
  })
names (meso_bulk_meta_l) = c('bueno','tcga','mesomics')

# Save data and metadata ####
saveRDS (meso_bulk_l, 'bulk_RNA_studies.rds')
saveRDS (meso_bulk_meta_l, 'bulk_RNA_studies_metadata.rds')



### Query bulk data ####
studies = c('bueno','tcga','mesomics')


module_l = c(LAG3 = 'LAG3', HAVCR2 = 'HAVCR2', PDCD1 = 'PDCD1', TIGIT = 'TIGIT', CTLA4 = 'CTLA4')
module_l = c(neuroendocrine1 = 'PMP2', neuroendocrine2 = 'VGF')
module_l = c(sox9 = 'SOX9', twist1 = 'TWIST1',SMARCC2 = 'SMARCC2',pitx2 = 'PITX2', mesp1 = 'MESP1', mef2a='MEF2A')
module_l = c(il32 = 'IL32')
module_l = c(SNAI2 = 'SNAI2')
module_l = c(CD90 = 'CD44', stem='CD73',stem='CD146')
module_l = c(ELK4 = 'ELK4')
module_l = c(HOXB13 = 'HOXB13')
module_l = c(SOX9 = 'SOX9', SOX6 = 'SOX6')


# Make gene modules overlapping megahubs regions in P11 ####
all_genes = genes (TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_in_region = all_genes$gene_id[subjectHits(findOverlaps (region, all_genes))]
module_l = list(chr18_q23 = as.data.frame(org.Hs.egSYMBOL)[match (genes_in_region, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol'])
module_l = c(MBP = 'MBP')
module_l = c(TXNL4A = 'TXNL4A')
module_l = list(HOXB13 = 'HOXB13', HOXC13 = 'HOXC13',sarc='AXL')
module_l = sarc_score
module_l = list(SOX9 = 'SOX9',SOX6 = 'SOX6',RUNX2='RUNX2',RUNX1='RUNX1',SNAI2='SNAI2')
module_l = list(NTM='NTM')

# Run genes on bulk datasets ####
stat_testL2 = list()
for (mod_name in names(module_l))
  {
  mod = module_l[which (mod_name == names(module_l))]
  bp_l = list()  
  stat_testL = list()  
  for (study in studies)
    {
    meso_bulk = meso_bulk_l[[study]]
    meso_bulk_meta = meso_bulk_meta_l[[study]]
    
    # gene_exp$histology = meso_bulk_meta[match (gene_exp$sample, names(meso_bulk_meta))]  
    # gene_exp = gene_exp[!is.na(gene_exp$histology),]
    gene_exp = data.frame(
      expression = colMeans (meso_bulk[rownames(meso_bulk) %in% mod[[1]],,drop=F]),
      subtype = meso_bulk_meta$subtype
      )
    bp_l[[study]] = ggplot(gene_exp, aes (x= subtype, y= expression)) + 
    ggtitle (paste('RNAseq',study,'cohort')) +
    gtheme +  
    geom_boxplot (aes (fill = subtype), outlier.colour="black", outlier.shape=16,
             outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
    NoLegend () + 
    scale_fill_manual (values = bulk_palette)
  
    stat.test = bp_l[[study]]$data %>%
          t_test(reformulate ('subtype', 'expression')) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
    stat.test = stat.test %>% add_xy_position (x = 'subtype', step.increase=0.5)
    stat.test = stat.test[stat.test$group1 %in% c('Sarcomatoid','Epithelioid') & stat.test$group2 %in% c('Sarcomatoid','Epithelioid'),]
    bp_l[[study]] = bp_l[[study]] + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = T,
      label = "p.adj.signif") + 
      NoLegend()
     stat_testL[[study]] = data.frame (stat.test, dataset = study)
    }
    stat_testL2[[mod_name]] = stat_testL
    pdf (paste0 ('Plots/',mod_name,'_expression_avg_boxplots.pdf'), width = 5,height=3)
    print (wrap_plots (bp_l, ncol=3))
    dev.off ()
    }


#####################################################################
# Correlation analyses between two genes or mean of set of genes ####
#####################################################################
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))

# Load P11 megahubs regions ####
region = readRDS ('../tumor_compartment/scatac_ArchR/P11_chr18_region.rds')

# Make gene modules overlapping megahubs regions in P11 ####
all_genes = genes (TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_in_region = all_genes$gene_id[subjectHits(findOverlaps (region, all_genes))]
genes_in_region = list(chr18_q23 = as.data.frame(org.Hs.egSYMBOL)[match (genes_in_region, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol'])
#genes_in_region = split (genes_in_region[[1]], genes_in_region)

#
your.gene1 = 'HOXB13'
your.gene2 = c('CALB2','ITLN1','UPK3B')
your.gene2 = head(cnmf_spectra_unique[['cNMF20']],20)
your.gene1 = head(cnmf_spectra_unique[['cNMF13']],20)
your.gene2 = genes_in_region[[1]]
your.gene2 = 'MBP'
your.gene2 = 'TXNL4A'
your.gene2 = genes_in_region[[48]]
your.gene2 = 'LYZ'
your.gene2 = 'BAP1'
your.gene2 = c('PMP2')
your.gene1 = 'SOX9'
your.gene2 = 'HSBP1L1'
your.gene2 = 'ATP9B'
your.gene2 = 'NFATC1'
#study = c('bueno_immune_corrected','tcga_immune_corrected','mesomics_immune_corrected')
corr_res = list()
studies = c('bueno','tcga','mesomics')
by_histology=F
filter_low_exp = 0
for (study in studies)
  {
  meso_bulk = meso_bulk_l[[study]] 
  meso_bulk_meta = meso_bulk_meta_l[[study]]

  if (length(your.gene1) > 1) gene_exp1 = colMeans(meso_bulk[rownames(meso_bulk) %in% your.gene1,,drop=F]) else
  gene_exp1 = as.vector (unlist (meso_bulk[your.gene1,,drop=F]))
  if (length(your.gene2) > 1) gene_exp2 = colMeans(meso_bulk[rownames(meso_bulk) %in% your.gene2,,drop=F]) else  
  gene_exp2 = as.vector (unlist (meso_bulk[your.gene2,,drop=F]))

  exp_df = data.frame (your.gene1 = gene_exp1, your.gene2 = gene_exp2)
  exp_df$subtype = meso_bulk_meta$subtype
  rownames (exp_df) = colnames (meso_bulk)
  exp_df = exp_df[!is.na(exp_df$subtype),]
  
  exp_df = exp_df[exp_df[,1] >= filter_low_exp & exp_df[,2] >= filter_low_exp,]
  if (by_histology) {
  corr_res[[study]] = ggscatter (
            exp_df, 
            x = 'your.gene1', 
            y = 'your.gene2',
            #palette = bulk_palette 
            shape=16,
            color = 'subtype',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = paste(your.gene1, collapse=' '), ylab = paste(your.gene2, collapse=' '),
            title = paste0(study,'(n=',nrow(exp_df),')'), fullrange = TRUE) + 
            scale_fill_manual (values=bulk_palette) + 
            scale_color_manual (values=bulk_palette) +
            facet_wrap (~subtype, drop=TRUE, scales = 'free_x', ncol=length (unique(exp_df$subtype))) +
            NoLegend()
  }else {
  corr_res[[study]] = ggscatter (
            exp_df, 
            x = 'your.gene1',
            y = 'your.gene2',
            #palette = bulk_palette 
            shape=16,
            color = 'subtype',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = paste(your.gene1, collapse=' '), ylab = paste(your.gene2, collapse=' '),
            title = paste0(study,'(n=',nrow(exp_df),')'), fullrange = TRUE) + 
            scale_fill_manual (values=bulk_palette) + 
            scale_color_manual (values=bulk_palette) +
            NoLegend()
  #corr_res[[blk]] = ggMarginal(corr_res[[blk]], type = "density", groupColour = TRUE, groupFill = TRUE)
  }
}

pdf (paste0 ( 'Plots/',your.gene1,'_',your.gene2,'_correlation_scatterplots2.pdf'),8,3)
wrap_plots (corr_res, ncol=3)
#corr_res
dev.off()


# Show expression of SOX6 / SOX9 in relation to the scS score
studies = c ('bueno','tcga','mesomics')
by_histology=F
filter_low_exp = 0
your.gene1 = 'SOX9'
your.gene2 = 'VIM'

nfeat = 5000
k_selection = 25
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))
sarc_score = head (cnmf_spectra_unique[[20]],50)
your.gene2 = sarc_score
sp_l = list()
for (study in studies)
  {
  meso_bulk = meso_bulk_l[[study]]
  meso_bulk_meta = meso_bulk_meta_l[[study]]
  if (length(your.gene1) > 1) gene_exp1 = colMeans(meso_bulk[rownames(meso_bulk) %in% your.gene1,,drop=F]) else
  gene_exp1 = as.vector (unlist (meso_bulk[your.gene1,,drop=F]))
  if (length(your.gene2) > 1) gene_exp2 = colMeans(meso_bulk[rownames(meso_bulk) %in% your.gene2,,drop=F]) else  
  gene_exp2 = as.vector (unlist (meso_bulk[your.gene2,,drop=F]))

  exp_df = data.frame (your.gene1 = gene_exp1, your.gene2 = gene_exp2, scs_score = meso_bulk_meta$sarc_score)
  exp_df$subtype = meso_bulk_meta$subtype
  rownames (exp_df) = colnames (meso_bulk)
  exp_df = exp_df[exp_df$your.gene1 > 1 & exp_df$your.gene1 > 1, ]
  #exp_df = exp_df[!is.na(exp_df$subtype),]
  
  sp = ggscatter (
            exp_df,
            x = 'your.gene1',
            y = 'your.gene2',
            #palette = bulk_palette,
            shape=16,
            color = 'scs_score',
            fill = 'scs_score',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = paste(your.gene1, collapse=' '), ylab = paste(your.gene2, collapse=' '),
            title = paste0(study,'(n=',nrow(exp_df),')'), fullrange = TRUE) + 
            ggtitle (study) +
            scale_fill_gradient2 (viridis::inferno(100)) +
            scale_color_gradient2 (viridis::inferno(100)) +
#            facet_wrap (~subtype, drop=TRUE, scales = 'free_x', ncol=length (unique(exp_df$subtype))) +
            #scale_color_manual (values=bulk_palette) +
            NoLegend()
    sp_l[[study]] = sp
    }

pdf (file.path ('Plots',paste(your.gene1,'_',your.gene2)), width=8, height=3)
wrap_plots (sp_l, ncol=3)
dev.off()

# Correlate TFs vs sarcomatoid score (cNMF20) and extract values ####
nfeat = 5000
k_selection = 25
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))
sarc_score = head (cnmf_spectra_unique[[20]],20)
activeTFs = read.csv ('../tumor_compartment/scatac_ArchR/Active_TFs.csv')[[2]]


study = c('bueno','tcga','mesomics')
cor_res_study = list()
by_histology=FALSE
for (blk in study)
  {
  meso_bulk = meso_bulk_l[[blk]] 
  gene_exp1 = t(meso_bulk[rownames(meso_bulk) %in% activeTFs,,drop=F]  )
  gene_exp2 = colMeans(meso_bulk[rownames(meso_bulk) %in% head(cnmf_spectra_unique[['cNMF20']],50),,drop=F])
  cor_res_study[[blk]] = cor (gene_exp1, gene_exp2, method ='spearman')
  }
cor_res_study = lapply (cor_res_study, function(x) x[activeTFs,])  
cor_res_studies = do.call (cbind, cor_res_study)  
colnames (cor_res_studies) = study
write.csv (cor_res_studies, 'activeTF_sarcomatoid_correlation.csv')


scrna_nmf = read.csv ('../tumor_compartment/scrna/cnmf_genelist_25_nfeat_5000.csv')[['cNMF20']][1:50]
scrna_nmf = scrna_nmf[!is.na(scrna_nmf)]

MESO_RNA_patient_id = setNames (as.character(meso_bulk_meta_l[['tcga']]), names(meso_bulk_meta_l[['tcga']]))
names(MESO_RNA_patient_id) = gsub ('\\.','-',names(MESO_RNA_patient_id))
tcga_id_mapping = read.table('TCGA_identifier_mapping.txt',sep='\t', header=T)
bigwig_ids = read.table ('MESO_ATAC_bigwig_identifiers.txt')
bigwig_ids = sapply (bigwig_ids[[1]], function(x) unlist(strsplit (x, '\\.'))[1])
bigwig_ids = gsub ('_','-',bigwig_ids)
MESO_ATAC_patient_id = tcga_id_mapping$Case_ID[match(bigwig_ids, tcga_id_mapping$bam_prefix)]
MESO_ATAC_patient_id = sapply (MESO_ATAC_patient_id, function(x) unlist(strsplit (x, 'A-31'))[1])
bigwig_atac_ids = MESO_RNA_patient_id[MESO_ATAC_patient_id]
#names (bigwig_atac_ids) = bigwig_ids

colnames(meso_bulk_l[['tcga']]) = gsub ('\\.','-',colnames(meso_bulk_l[['tcga']]))
sarc_score = colMeans (na.omit(as.data.frame(meso_bulk_l[['tcga']])[scrna_nmf,names(bigwig_atac_ids)]))
#sarc_score = sarc_score[order(-sarc_score)]
bigwig_atac_ids_df = data.frame (TCGA_ID = names(bigwig_atac_ids), bam_ID = bigwig_ids, histology = bigwig_atac_ids, sarc_score = sarc_score)
write.csv (bigwig_atac_ids_df, 'TCGA_ATAC_ID_sarc_score_histology.csv')
# cluster TCGA rna data to fill histology annotation for the missing ATAC sample
# srt = CreateSeuratObject (counts = tcga_mat, data = log2(tcga_mat+1))
# srt = FindVariableFeatures (srt)
# tcga_mat_vf = as.data.frame(as.matrix(srt@assays$RNA$data))[VariableFeatures(srt),]
# colnames (tcga_mat_vf) = gsub ('\\.','-', colnames(tcga_mat_vf))
# tcga_mat_vf = tcga_mat_vf [,names(MESO_RNA_patient_id)]
# d = as.dist (1-cor(tcga_mat_vf))
# # Hierarchical clustering using Complete Linkage
# hc <- as.dendrogram(hclust(d))
# hc = color_labels(hc, colors = c('red','blue'))bulk_palette[MESO_RNA_patient_id])
# ggplotDend <- as.ggdend(hc)
# hcp = ggplot(ggplotDend) + 
# #  scale_fill_manual (values = bulk_palette) +
#   theme_minimal() +
#   labs(title = "Dendrogram with Colored Leaves")
# # Plot
# pdf ('Plots/hierarchical_clust.pdf', width=15, height=16)
# hcp
# dev.off()

tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
tcga_mat = assays (tcga_dev)$z
tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
rownames (tcga_mat) = tcga_mat[,1]
tcga_mat = tcga_mat[,-1]
tcga_mat = as.data.frame (t(tcga_mat))
rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))
#tcga_mat = rowMeans (tcga_mat)


### Correlate gene with rest of genes across samples #####
studies = c('tcga','bueno','mesomics')

cor_mat_l = list()
for (study in studies)
  {
  your.gene1 = 'HOXB13'
  your.gene2 = rownames(meso_bulk_l[[study]])
  no_exp = meso_bulk_l[[study]][your.gene1,] > 0
  cor_mat = as.data.frame (t(cor (t(meso_bulk_l[[study]][your.gene1,no_exp, drop=F]), t(meso_bulk_l[[study]][your.gene2,no_exp]), method = 'spearman')))
  cor_mat_l[[study]] = cor_mat[order (-cor_mat[[1]]),,drop=F]
  }
shared_genes = rownames(cor_mat_l[[1]])[rownames(cor_mat_l[[1]]) %in%  rownames(cor_mat_l[[2]])]
shared_genes = shared_genes[shared_genes %in%  rownames(cor_mat_l[[3]])]
cor_mat_df = do.call (cbind, lapply (cor_mat_l, function(x) x[shared_genes,]))
rownames (cor_mat_df) = shared_genes
cor_mat_mean = rowMeans (cor_mat_df)
cor_mat_mean = cor_mat_mean[order (-cor_mat_mean)]
cor_mat_mean = na.omit(cor_mat_mean) 
write.csv (as.data.frame(tail (cor_mat_mean,100)),'HOXB13_sig_mean_cor_studies.csv')
write.csv (as.data.frame(head (cor_mat_mean,100)),'HOXB13_sig_mean_cor_studies_head.csv')


# Compare with P11 HOXB13 correlated genes ####
p11 = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/correlated_genes_p11s_HOXB13.csv')
p11_flt = p11[p11[,1] %in% names(cor_mat_mean),]
cor_mat_mean_flt = cor_mat_mean[names (cor_mat_mean) %in% p11_flt[,1]]
cor.test(1:length(cor_mat_mean_flt), match (names(cor_mat_mean_flt), p11[,1]), method='spearman')
head(match (names(cor_mat_mean_flt), p11[,1]), 200)[which (head(match (names(cor_mat_mean_flt), p11[,1]), 200) < 50)]
p11[21,]

library (fgsea)
enricher_universe = names(cor_mat_mean_flt)
#do.fgsea = TRUE
gmt_annotations = c(
'h.all.v7.4.symbols.gmt',#,
'c5.bp.v7.1.symbol.gmt',
'c3.tft.v7.1.symbol.gmt'
)
gmt_annotation = gmt_annotations[2]
if (!file.exists (paste0('EnrichR_HOXB13_cor_genes_gmt_annotation_',gmt_annotation,'.rds')) | force)
  {
    gmt.file = file.path ('..','git_repo','files',gmt_annotation)
    pathways = clusterProfiler::read.gmt (gmt.file)
    EnrichRResCluster = list()
      sig_genes = cor_mat_mean_flt
      #if (!any (sig_genes %in% pathways$gene)) next
      egmt = clusterProfiler::enricher (tail(names(sig_genes),100), TERM2GENE=pathways, universe = enricher_universe)
      egmt@result$ID = stringr::str_trunc (egmt@result$ID, 100)
      EnrichRResCluster = egmt@result
      }
    EnrichRResAll = EnrichRResCluster
    saveRDS (EnrichRResAll, paste0('EnrichR_HOXB13_cor_genes_gmt_annotation_',gmt_annotation,'.rds'))
    } else {
  EnrichRResAll = readRDS (paste0('EnrichR_HOXB13_cor_genes_gmt_annotation_',gmt_annotation,'.rds'))
  }
  pvalAdjTrheshold = 0.05
  top_pathways = 10
  EnrichRes_dp = dotGSEA (enrichmentsTest_list = EnrichRResAll, type = 'enrich', padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways)
  pdf (paste0(cnmf_out, '/Plots/EnrichR_',i,'_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'_dotplots.pdf'), width = 12 + length(length (cnmf_spectra_unique))/10, height = 15)
  print (EnrichRes_dp)
  dev.off()


# check PRC2 and BAF complexes expression ####
prc2 = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/PRC2.csv', skip=1, header=T)

cor_mat_mean_flt[prc2$Approved.symbol]











###########################
### Survival Analysis #####
###########################
nfeat = 5000
k_selection = 25
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))
sarc_score = head (cnmf_spectra_unique[[20]],20)

# Load P11 megahubs regions ####
region = readRDS ('../tumor_compartment/scatac_ArchR/P11_chr18_region.rds')

# Make gene modules overlapping megahubs regions in P11 ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_in_region = all_genes$gene_id[subjectHits(findOverlaps (region, all_genes))]
genes_in_region = list(chr18_q23 = as.data.frame(org.Hs.egSYMBOL)[match (genes_in_region, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol'])
genes_in_region = split (genes_in_region[[1]], genes_in_region)

module_l = genes_in_region
#genes_in_region = split (module_l,module_l)

# HOX genes ####
module_l = rownames(meso_bulk_l[[1]])[grep ('^HOX',rownames(meso_bulk_l[[1]]))]
module_l = split (module_l, module_l)
module_l = lapply (cnmf_spectra_unique, function(x) head(x,20))
# Make metadata for survival analysis ####
# Define module of genes to check ####
module_l = list(TCF3 = 'TCF3', PITX1 = 'PITX1', TEAD1 = 'TEAD1')
selected_TFs = readRDS ('../tumor_compartment/scatac_ArchR/selected_TF.rds')
module_l = split (selected_TFs, selected_TFs)
module_l = list (SOX9 = 'SOX9', RUNX2 = 'RUNX2', SNAI2 = 'SNAI2')
#### Check chromatin regulators and genetic drivers genes ####
# chromatin_regulators = c('SETD5, ASH1L, CREBBP, PRDM2, KDM2B, KMT2D, EZH2, SETDB1')
# chromatin_regulators = unlist(strsplit (chromatin_regulators, ', '))
# genetic_drivers = c('BAP1, NF2, CDKN2A, CDKN2B, TP53, LATS2, SETD2, FAT4, PTCH1')
# genetic_drivers = unlist(strsplit (genetic_drivers, ', '))
# chromatin_regulators2 = c('LATS1, DDX3X, ULK2, RYR2, CFAP45, SETDB1, DDX51, SF3B1, TRAF7, PTEN, RBFOX1, CSMD1, MTAP, TTC28, PCDH15, USH2A, CNTNAP2, DNAH1, KCNH7, PTK2, ROBO2, DLG2, PBRM1, PTPRD, ANTXR2, CTNNA3, LINGO2, LRP1B, PLCB1, UNC79, WWOX')
# chromatin_regulators2 = unlist (strsplit (chromatin_regulators2, ', '))
# genetic_drivers = unlist(strsplit (genetic_drivers, ', '))
# module_l = unique (c(chromatin_regulators, genetic_drivers, chromatin_regulators2))
# module_l = split (module_l, module_l)


meso_bulk_meta_l2 = meso_bulk_meta_l
meso_bulk_meta_l2 = lapply (seq_along(meso_bulk_l), function(x) 
  {
   tmp = lapply (module_l, function(y) colMeans (meso_bulk_l[[x]][rownames(meso_bulk_l[[x]]) %in% y,,drop=F]))
   tmp = do.call (cbind, tmp)
   colnames (tmp) = gsub('-','_',colnames(tmp))
   tmp = tmp[,apply(tmp, 2, function(t) !any(is.na(t)))]
   meso_bulk_meta_l2[[x]] = as.data.frame (meso_bulk_meta_l2[[x]])
   meso_bulk_meta_l2[[x]] = cbind (meso_bulk_meta_l2[[x]], tmp)
   if(length(module_l) == 1) colnames(meso_bulk_meta_l2[[x]])[ncol(meso_bulk_meta_l2[[x]])] = names(module_l)
   meso_bulk_meta_l2[[x]]
  })
names (meso_bulk_meta_l2) = c('bueno','tcga','mesomics')
#survival_name = 'activeTFs'
#survival_name = 'megahub'
survival_name = 'plasticity'



# Run Cox hazard ratio regression survival analysis ####
# Set variables per dataset to use
low='1st Qu.'
high='3rd Qu.'  
studies = c('bueno','tcga','mesomics')
cfit_study = list()
cox_l = list()
sc_p = list()
plot_study= list()
for (study in studies)
    {
    cfit = list()
    module_l = lapply (module_l, function(x) gsub('-','_',x))
    names (module_l) = gsub('-','_',names(module_l))
    mods = colnames(meso_bulk_meta_l2[[study]])[colnames(meso_bulk_meta_l2[[study]]) %in% names(module_l)]
    for (mod in mods)
        {
        meta_surv = meso_bulk_meta_l2[[study]]
        meta_surv = meta_surv[!is.na(meta_surv$census),]
        
        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
                            status) ~', mod, '+ strata (subtype)'))
        # form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
        #                     status) ~', mod, '+ sarc_score'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        CI <- round(exp(confint(cfit[[mod]])), 2)
        cox_df = data.frame (
          HR = round(exp(coef(cfit[[mod]])), 2),
          CI = paste0('(',paste (CI, collapse=' - '),')'),
          LL = CI[1],
          UL = CI[2],
          P_value_C = round(summary(cfit[[mod]])$coefficients[, 5],2),
          label = mod
          )
        
        raw.vec=meta_surv[,mod]
        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        meta_surv[,mod] = factor (classified.vec, levels = c('Low','Med','High'))

        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
                            status) ~', mod, '+ strata (subtype)'))
        # form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
        #                     status) ~', mod, '+ sarc_score'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        s = summary (cfit[[mod]])
        cox_df$P_value_S = round(s$logtest[3],2)  #s$logtest[3]
        cox_l[[mod]] = cox_df
        
        sc_p[[mod]] = ggadjustedcurves (cfit[[mod]], 
                data = meta_surv, 
                method = "conditional",
                variable=mod,
                palette = 'aaas',
                size=0.4,
                surv.median.line = 'hv',
                ggtheme = theme_classic()) +
                labs (title = mod,
                subtitle = paste0('log-rank = ',round(s$logtest[3],2)),
                caption = paste("n = ", nrow(meta_surv))) +
                      ylim(0,1)+ geom_hline(yintercept = 0.5,c(0.5,0.5),linetype='dotted', col='grey22') 
        }
    cox_df = do.call (rbind, cox_l)
    cox_df$Index = rownames (cox_df)
    cox_df$Index = factor (cox_df$Index, levels = rev(cox_df$Index))
    cfit_study[[study]] = cox_df
    plot_study[[study]] = sc_p
    }

for (study in names (cfit_study))
  {
  forest <- ggplot(cfit_study[[study]], aes(y = Index, x = HR)) + 
    geom_point(shape = 18, size = 5) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:nrow(cfit_study[[study]]), labels = cfit_study[[study]]$label, trans = "reverse", expand = expansion(add = 0.5)) +
    #scale_x_continuous(trans = 'log10')   + 
    xlab("Hazard Ratio") + 
    ylab(" ") + 
    theme_classic()
      
  tab <- ggplot(cfit_study[[study]], aes(y = Index)) +
    geom_text(aes(x = 0, label = sprintf("%0.1f", round(HR, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = CI), size = 4) + 
    geom_text(aes(x = 2, label = P_value_C), size = 4) + 
    geom_text(aes(x = 3, label = P_value_S), size = 4) + 
    #scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_x_continuous(
      breaks = 0:3, labels = c('HR', 'CI', 'P value (C)','P value (S)'), 
      position = 'top', expand = expansion(add = 0.5)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))
    
    pdf (paste0('Plots/cox_regression_',study,'_',survival_name,'.pdf'), height=8,5)
    print (forest + tab + plot_layout(ncol = 2, widths = c(1, 3)))
    dev.off()
    pdf (paste0('Plots/cox_regression_',study,'_stratified_',survival_name,'.pdf'), height = 2.8,3)
    print (plot_study[[study]])
    dev.off()
  }
  
# make heatmap of cox regression models ####
ht = do.call (cbind, cfit_study)
ht2 = ht[,grep ('HR', colnames(ht))]
ht3 = ht[,grep ('P_value_C', colnames(ht))]
ht3 = ifelse (ht3 > .05,0,1)
htl = list (hr = ht2, p = ht3)
saveRDS (htl, 'oncoTF_sig_cox_bulkRNA.rds')


# Correlate HOXB13 expression with TME proportions ####
bp_bueno = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/bulkRNA/bayesprism_bueno_theta.rds')
bp_tcga = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/bulkRNA/bayesprism_tcga_theta.rds')

bp_bueno = bp_bueno[colnames (meso_bulk_l[['bueno']]),]
cor (t(meso_bulk_l[['bueno']]['HOXB13',]),bp_bueno, method='spearman')

bp_tcga = bp_tcga[colnames (meso_bulk_l[['tcga']]),]
bp_tcga = na.omit (bp_tcga)
meso_bulk_tcga = meso_bulk_l[['tcga']][,rownames(bp_tcga)]
cor (meso_bulk_tcga['HOXB13',],bp_tcga, method='spearman')





### Check HOX genes and make heatmap ####
srt = readRDS (file.path ('..','main','scrna','srt.rds'))
srt = srt[,srt$sampleID == 'P11']
hoxb13_sig = read.csv ('../tumor_compartment/scatac_scrna_P11/HOXB13_cnmf.csv')
hox_sig = head(hoxb13_sig[,2],Inf)

pdf (file.path('Plots','HOX_signature_dotplot.pdf'), width=30)
DotPlot (srt, feature = hox_sig) + gtheme
dev.off()

hox_sig_specific = c('HOXB13','CA8','SULT1E1','SYT1','GAS2','WDR72','CHST9','RPRM','LY6G6D','COL9A1','MPPED2','TEKT3','CLIC5','NKX2−5','S100A7','ASPG','BEX1','GABRA2','PKP2','TDRD10','GRM1','PLPPR3',
  'PI3','ACTR3B','LY6H','TNNT2','CXCL14')

exp_mat = log2(meso_bulk_l[[2]]+1)[rownames(meso_bulk_l[[2]]) %in% hox_sig_specific,]
ha = HeatmapAnnotation(df = meso_bulk_meta_l[[2]][,c(
  'subtype',
  'CANCER_TYPE_DETAILED',
  'AJCC_PATHOLOGIC_TUMOR_STAGE',
  'AGE',
  'ANEUPLOIDY_SCORE'
  )])
hm = Heatmap (t(scale(t(exp_mat))), top_annotation = ha)

pdf (file.path('Plots','HOX_signature_tcga_heatmap.pdf'), width=15)
hm
dev.off()

exp_mat = log2(meso_bulk_l[[1]]+1)[rownames(meso_bulk_l[[1]]) %in% hox_sig_specific,]
ha = HeatmapAnnotation (df = meso_bulk_meta_l[[1]][,c('subtype','PD.L1.expression..RPKM.','Asbestos.exposure','NF2..FISH.','FISH.chrom3','Type.of.pre.op.Treatment','Fibers.gm.lung','Sex')])
hm = Heatmap (t(scale(t(exp_mat))), top_annotation = ha)

pdf (file.path('Plots','HOX_signature_bueno_heatmap.pdf'), width=15)
hm
dev.off()

exp_mat = log2(meso_bulk_l[[3]]+1)[rownames(meso_bulk_l[[3]]) %in% hox_sig_specific,]
ha = HeatmapAnnotation (df = meso_bulk_meta_l[[3]][,c('subtype','Subtype2','Cytological.Variant1','Stroma1','Chimio','Radioth','Immunoth','History.Location')])
hm = Heatmap (t(scale(t(exp_mat))), top_annotation = ha)

pdf (file.path('Plots','HOX_signature_mesomics_heatmap.pdf'), width=15)
hm
dev.off()


#### Most correlated genes to gene ####
gene = 'SOX9'

cor_res = lapply(names (meso_bulk_l), function(x) 
  {
   tmp = cor(t(meso_bulk_l[[x]][gene, ,drop=F]), t(meso_bulk_l[[x]]), method = 'pearson')
    tmp = as.data.frame (t(tmp))
    tmp[rownames (meso_bulk_l[[1]]),, drop=F]
  })
cor_res_df = do.call (cbind, cor_res)
colnames (cor_res_df) = names(meso_bulk_l)
cor_res_df = cor_res_df[order(-cor_res_df[,1]),]

cor_res_df = cbind (cor_res_df, cor_res_df2[rownames(cor_res_df), ])
cor (cor_res_df, method= 'spearman', use = 'pairwise.complete.obs')
write.csv (cor_res_df, paste0(gene, '_gene_correlations.csv'))




