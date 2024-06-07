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
  'wesanderson',
  'EnhancedVolcano',
  'rstatix',
  'tidyr',
  'corrplot',
  'scales',
  'ggExtra',
  'chromvar',
  'limma',
  'ape',
  'dendextend')
lapply(packages, require, character.only = TRUE)

set.seed(1234)
# Set directory
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd(projdir)

bulk_palette = setNames (hue_pal()(4)[c(4,3,2,1)],c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
bulk_palette = setNames (as.character(paletteer::paletteer_d("rockthemes::heep")), c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid'))
bueno_palette = setNames (paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7)], c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
bulk_palette = setNames (as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7,3)]), c('Epithelioid','Biphasic-E','Biphasic-S','Sarcomatoid','Biphasic'))


if (!file.exists ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/TCGA_meso_RNA/meso_tcga_RSEM_cleaned.rds'))
  {
  #### Import TCGA ####
  # load ('/ahg/regevdata/projects/lungCancerBueno/Results/TCGAsubtypes/Cancer_All/Results/gcdata.merge.Rda')
  # gcdata.merge_meso = subset (gcdata.merge, meso_histology %in% c('MESO.E','MESO.NOS','MESO.S','MESO.BP')) 
  # tcga_meta = gcdata.merge_meso$meso_histology
  # all (colnames (gcdata.merge_meso) == names (tcga_meta))
  
  
  ## Load TCGA data ###
  tcga_mat = read.table ('/ahg/regevdata/projects/ICA_Lung/Bruno/TCGA/MESO.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep='\t', header=T)
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
  # tcga_mat = tcga_mat[,tcga_meta$TUMOR_TYPE != 'Diffuse Malignant Mesothelioma (NOS)']
  # tcga_meta = tcga_meta[tcga_meta$TUMOR_TYPE != 'Diffuse Malignant Mesothelioma (NOS)',]
  
  saveRDS (tcga_mat, '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/TCGA_meso_RNA/meso_tcga_RSEM_cleaned.rds')
  saveRDS (tcga_meta, '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/TCGA_meso_RNA/meso_tcga_meta_cleaned.rds')
  
  # Check if it is seq-depth normalized
  pdf ('tcga_raw.pdf')
  print(boxplot (log2(as.matrix(tcga_mat))))
  dev.off()
  } else {
tcga_mat = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/TCGA_meso_RNA/meso_tcga_RSEM_cleaned.rds')
tcga_meta = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/TCGA_meso_RNA/meso_tcga_meta_cleaned.rds')  
}


if (!file.exists ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/bueno_rpkm_cleaned.rds'))
  {
  # Import bueno from my folder
  bueno_mat = read.table ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/bueno_meso_rpkm_rnaseq.txt', header=T, sep= '\t')
  bueno_mat = bueno_mat[!duplicated(bueno_mat[,1]),]
  rownames (bueno_mat) = bueno_mat[,1]
  bueno_mat = bueno_mat[,-1]
  
  # Get clinical data
  bueno_meta2 = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/ClincalInfo_S2_2.csv')
  rownames (bueno_meta2) = bueno_meta2$Tumor.ID
  
  bueno_meta2 = bueno_meta2[bueno_meta2$ConsensusCluster != '',]
  bueno_meta2 = setNames (bueno_meta2$ConsensusCluster,rownames(bueno_meta2))
  bueno_mat = bueno_mat[,names(bueno_meta2)]
  
  saveRDS (bueno_meta2, '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/bueno_meta_cleaned.rds')
  saveRDS (bueno_mat, '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/bueno_rpkm_cleaned.rds')
    
    
  # Check if it is seq-depth normalized
  pdf ('bueno_rpkm.pdf')
  print (boxplot (log2(bueno_mat)))
  dev.off()
  } else {
  bueno_meta2 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/bueno_meta_cleaned.rds')
  bueno_mat = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/bueno_data/bueno_rpkm_cleaned.rds')
  }


#### Import MESOMICS ####
msm = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/Mesothelioma/Data/Expression/TPM_COnverted_GeneCount.tsv')
#msm = log2(msm+1)
#msm = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/Mesothelioma/Data/Expression/VStexpr.tsv', header=T)
msm_meta = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/Mesothelioma/Data/Mesothelioma_SampleInfo.txt', header=T)
msm = msm[,colnames (msm) %in% msm_meta$Sample]
msm_meta = msm_meta[msm_meta$Sample %in% colnames (msm),]
msm = msm[, msm_meta$Sample]
all (colnames (msm) == msm_meta$Sample)


#### Correct for immune infiltration using limma
bueno_meta3 = bueno_meta2 
bueno_meta3 = cbind (bueno_meta3, data.frame (PTPRC = rowMeans(t(log2(bueno_mat+1)[c('PTPRC'),]))))
bueno_meta3 = cbind (bueno_meta3, data.frame (VWF = rowMeans(t(log2(bueno_mat+1)[c('VWF'),]))))
bueno_meta3 = cbind (bueno_meta3, data.frame (CD3D = rowMeans(t(log2(bueno_mat+1)[c('CD3D'),]))))

tcga_meta3 = tcga_meta
tcga_meta3 = cbind (tcga_meta3, data.frame (PTPRC = rowMeans(t(log2(as.data.frame (tcga_mat)+1)[c('PTPRC'),]))))
tcga_meta3 = cbind (tcga_meta3, data.frame (VWF = rowMeans(t(log2(as.data.frame (tcga_mat)+1)[c('VWF'),]))))
tcga_meta3 = cbind (tcga_meta3, data.frame (CD3D = rowMeans(t(log2(as.data.frame (tcga_mat)+1)[c('CD3D'),]))))

mesomics_meta3 = msm_meta
mesomics_meta3 = cbind (mesomics_meta3, data.frame (PTPRC = rowMeans(t(log2(as.data.frame (msm)+1)[c('PTPRC'),]))))
mesomics_meta3 = cbind (mesomics_meta3, data.frame (VWF = rowMeans(t(log2(as.data.frame (msm)+1)[c('VWF'),]))))
mesomics_meta3 = cbind (mesomics_meta3, data.frame (CD3D = rowMeans(t(log2(as.data.frame (msm)+1)[c('CD3D'),]))))

bueno_mat_immune_corrected = removeBatchEffect(log2(as.matrix(bueno_mat)+1), covariates=bueno_meta3[,'PTPRC'])
bueno_mat_endo_corrected = removeBatchEffect(log2(as.matrix(bueno_mat)+1), covariates=bueno_meta3[,'VWF'])
bueno_mat_tcells_corrected = removeBatchEffect(log2(as.matrix(bueno_mat)+1), covariates=bueno_meta3[,'CD3D'])

tcga_mat_immune_corrected = removeBatchEffect(log2(as.matrix(tcga_mat)+1), covariates=tcga_meta3[,'PTPRC'])
tcga_mat_endo_corrected = removeBatchEffect(log2(as.matrix(tcga_mat)+1), covariates=tcga_meta3[,'VWF'])
tcga_mat_tcells_corrected = removeBatchEffect(log2(as.matrix(tcga_mat)+1), covariates=tcga_meta3[,'CD3D'])

msm_mat_immune_corrected = removeBatchEffect(log2(as.matrix(msm)+1), covariates=mesomics_meta3[,'PTPRC'])
msm_mat_endo_corrected = removeBatchEffect(log2(as.matrix(msm)+1), covariates=mesomics_meta3[,'VWF'])
msm_mat_tcells_corrected = removeBatchEffect(log2(as.matrix(msm)+1), covariates=mesomics_meta3[,'CD3D'])

### Make list including all the bulk data ####
meso_bulk_l = list (
  bueno_rpkm = log2(bueno_mat[,names(bueno_meta2)] + 1),
  bueno_immune_corrected = bueno_mat_immune_corrected[,names(bueno_meta2)],
  bueno_endo_corrected = bueno_mat_endo_corrected[,names(bueno_meta2)],
  bueno_tcells_corrected = bueno_mat_tcells_corrected[,names(bueno_meta2)],

  tcga = log2 (tcga_mat + 1),
  tcga_immune_corrected = tcga_mat_immune_corrected[,rownames(tcga_meta)],
  tcga_endo_corrected = tcga_mat_endo_corrected[,rownames(tcga_meta)],
  tcga_tcells_corrected = tcga_mat_tcells_corrected[,rownames(tcga_meta)],

  mesomics = log2 (msm+1),
  mesomics_immune_corrected = msm_mat_immune_corrected[, msm_meta$Sample])


meso_bulk_meta_l = list (
  bueno_rpkm = bueno_meta2,
  bueno_immune_corrected = bueno_meta2,
  bueno_endo_corrected = bueno_meta2,
  bueno_tcells_corrected = bueno_meta2,
  bueno_bp_tumor = bueno_meta2,
  bueno_bp_monomac = bueno_meta2,
  bueno_bp_tcells = bueno_meta2,

  tcga = setNames (tcga_meta$TUMOR_TYPE, rownames(tcga_meta)),
  tcga_immune_corrected = setNames (tcga_meta$TUMOR_TYPE, rownames(tcga_meta)),
  tcga_endo_corrected = setNames (tcga_meta$TUMOR_TYPE, rownames(tcga_meta)),
  tcga_tcells_corrected = setNames (tcga_meta$TUMOR_TYPE, rownames(tcga_meta)),

  mesomics = setNames (msm_meta$Type, msm_meta$Sample),
  mesomics_immune_corrected = setNames (msm_meta$Type, msm_meta$Sample))

meso_bulk_meta_l[['tcga']][meso_bulk_meta_l[['tcga']] == 'Biphasic Mesothelioma'] = 'Biphasic'
meso_bulk_meta_l[['tcga']][meso_bulk_meta_l[['tcga']] == 'Epithelioid Mesothelioma'] = 'Epithelioid'
meso_bulk_meta_l[['tcga']][meso_bulk_meta_l[['tcga']] == 'Sarcomatoid Mesothelioma'] = 'Sarcomatoid'

meso_bulk_meta_l[['mesomics']][meso_bulk_meta_l[['mesomics']] == 'MMS'] = 'Sarcomatoid'
meso_bulk_meta_l[['mesomics']][meso_bulk_meta_l[['mesomics']] == 'MME'] = 'Epithelioid'
meso_bulk_meta_l[['mesomics']][meso_bulk_meta_l[['mesomics']] == 'MMB'] = 'Biphasic'

# meso_bulk_meta_l[['bueno_log2']] = factor (meso_bulk_meta_l[['bueno_log2']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
# meso_bulk_meta_l[['bueno_raw_log2']] = factor (meso_bulk_meta_l[['bueno_raw_log2']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
# meso_bulk_meta_l[['bueno_log10']] = factor (meso_bulk_meta_l[['bueno_log10']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
#meso_bulk_meta_l[['bueno_maggie']] = factor (meso_bulk_meta_l[['bueno_maggie']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_rpkm']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_immune_corrected']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_endo_corrected']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_tcells_corrected']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_bp_tumor']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_bp_monomac']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['bueno_bp_tcells']] = factor (meso_bulk_meta_l[['bueno_rpkm']], levels = c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
meso_bulk_meta_l[['tcga']] = factor (meso_bulk_meta_l[['tcga']], levels = c('Sarcomatoid','Biphasic','Epithelioid'))
meso_bulk_meta_l[['tcga_immune_corrected']] = factor (meso_bulk_meta_l[['tcga']], levels = c('Sarcomatoid','Biphasic','Epithelioid'))
meso_bulk_meta_l[['tcga_endo_corrected']] = factor (meso_bulk_meta_l[['tcga']], levels = c('Sarcomatoid','Biphasic','Epithelioid'))
meso_bulk_meta_l[['tcga_tcells_corrected']] = factor (meso_bulk_meta_l[['tcga']], levels = c('Sarcomatoid','Biphasic','Epithelioid'))
meso_bulk_meta_l[['mesomics']] = factor (meso_bulk_meta_l[['mesomics']], levels = c('Sarcomatoid','Biphasic','Epithelioid'))
meso_bulk_meta_l[['mesomics_immune_corrected']] = factor (meso_bulk_meta_l[['mesomics_immune_corrected']], levels = c('MMS','MMB','MME'))




### Query bulk data ####
module_l = c(LAG3 = 'LAG3', HAVCR2 = 'HAVCR2', PDCD1 = 'PDCD1', TIGIT = 'TIGIT', CTLA4 = 'CTLA4')
module_l = c(neuroendocrine1 = 'PMP2', neuroendocrine2 = 'VGF')
module_l = c(sox9 = 'SOX9', twist1 = 'TWIST1',SMARCC2 = 'SMARCC2',pitx2 = 'PITX2', mesp1 = 'MESP1', mef2a='MEF2A')
module_l = c(il32 = 'IL32')
module_l = c(SNAI2 = 'SNAI2')
module_l = c(CD90 = 'CD44'),stem='CD73',stem='CD146')
module_l = c(ELK4 = 'ELK4')
study = c('bueno_rpkm','tcga','mesomics')

# Run genes on bulk datasets
do.avg=T
stat_testL2 = list()
for (mod_name in names(module_l))
  {
  mod = module_l[which (mod_name == names(module_l))]
  bp_l = list()  
  stat_testL = list()  
  for (blk in study)
    {
    meso_bulk = meso_bulk_l[[blk]] 
    meso_bulk_meta = meso_bulk_meta_l[[blk]] 
      if (do.avg) 
          {
          if (sum(names(mod) == mod_name) > 1)
            {
            tmpL = lapply (mod, function (z) colMeans(meso_bulk[rownames(meso_bulk) %in% z,]))        
            gene_exp = as.data.frame (rowMeans (do.call (cbind, tmpL)))
            } else {
          gene_exp = na.omit(as.data.frame(meso_bulk[rownames(meso_bulk) %in% mod[[1]],,drop=F]))
          gene_exp = as.data.frame (colMeans (gene_exp))
          }
          colnames (gene_exp) = 'expression'
          gene_exp$sample = rownames(gene_exp)  
          } else {
          gene_exp = as.data.frame(t(as.matrix(meso_bulk[mod[[1]],,drop=F])))
          colnames(gene_exp) = 'expression'
          gene_exp$sample = rownames (gene_exp)
          }     
        gene_exp$histology = meso_bulk_meta[match (gene_exp$sample, names(meso_bulk_meta))]  
        gene_exp = gene_exp[!is.na(gene_exp$histology),]
        bp_l[[blk]] = ggplot(gene_exp, aes (x= histology, y= expression)) + 
        ggtitle (paste('RNAseq',blk,'cohort')) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(size = 5)) + 
        geom_boxplot (aes (fill = histology), outlier.colour="black", outlier.shape=16,
                 outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
        NoLegend () +
        scale_fill_manual (values = bulk_palette)
  
        stat.test = bp_l[[blk]]$data %>%
              t_test(reformulate ('histology', 'expression')) %>%
              adjust_pvalue (method = "none") %>%
              add_significance ()
        stat.test = stat.test %>% add_xy_position (x = 'histology', step.increase=0.5)
        stat.test = stat.test[stat.test$group1 %in% c('Sarcomatoid','Epithelioid') & stat.test$group2 %in% c('Sarcomatoid','Epithelioid'),]
        bp_l[[blk]] = bp_l[[blk]] + stat_pvalue_manual (stat.test, 
          remove.bracket=FALSE,
          bracket.nudge.y = 0, 
          hide.ns = T,
          label = "p.adj.signif") + 
          NoLegend()
         stat_testL[[blk]] = data.frame (stat.test, dataset = blk)
        }
    stat_testL2[[mod_name]] = stat_testL
    png (paste0 (projdir,'Plots/',mod_name,'_expression_avg_',do.avg,'_boxplots.png'), width = 4000,height=1000, res=300)
    print (wrap_plots (bp_l, ncol=10))
    dev.off ()
    }


##################################################################
# Correlation analyses between two genes or mean of set of genes 
##################################################################
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))

your.gene1 = 'MESP1'
your.gene2 = head(cnmf_spectra_unique[['cNMF20']],50)

corr_res = list()
#study = c('bueno_immune_corrected','tcga_immune_corrected','mesomics_immune_corrected')
study = c('bueno_rpkm','tcga','mesomics')
by_histology=FALSE
for (blk in study)
  {
  meso_bulk = meso_bulk_l[[blk]] 
  meso_bulk_meta = meso_bulk_meta_l[[blk]]

  if (length(your.gene1) > 1) gene_exp1 = colMeans(meso_bulk[rownames(meso_bulk) %in% your.gene1,,drop=F]) else  
  gene_exp1 = as.vector (unlist (meso_bulk[your.gene1,,drop=F]))
  if (length(your.gene2) > 1) gene_exp2 = colMeans(meso_bulk[rownames(meso_bulk) %in% your.gene2,,drop=F]) else  
  gene_exp2 = as.vector (unlist (meso_bulk[your.gene2,,drop=F]))

  exp_df = data.frame (your.gene1 = gene_exp1, your.gene2 = gene_exp2)
  rownames (exp_df) = colnames (meso_bulk)
  exp_df$histology = meso_bulk_meta[match (rownames(exp_df), names(meso_bulk_meta))]
  exp_df = exp_df[!is.na(exp_df$histology),]
  if (by_histology) {
  corr_res[[blk]] = ggscatter (
            exp_df, 
            x = 'your.gene1', 
            y = 'your.gene2',
            #palette = bulk_palette 
            shape=16,
            color = 'histology',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(your.gene1, collapse=' '), ylab = paste(your.gene2, collapse=' '),
            title = paste0(blk,'(n=',nrow(exp_df),')'), fullrange = TRUE) + 
            scale_fill_manual (values=bulk_palette) + 
            scale_color_manual (values=bulk_palette) +
            facet_wrap (~histology, drop=TRUE, scales = 'free_x', ncol=length (unique(exp_df$histology))) +
            NoLegend()
  }else {
  corr_res[[blk]] = ggscatter (
            exp_df, 
            x = 'your.gene1',
            y = 'your.gene2',
            #palette = bulk_palette 
            shape=16,
            color = 'histology',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(your.gene1, collapse=' '), ylab = paste(your.gene2, collapse=' '),
            title = paste0(blk,'(n=',nrow(exp_df),')'), fullrange = TRUE) + 
            scale_fill_manual (values=bulk_palette) + 
            scale_color_manual (values=bulk_palette) +
            NoLegend()
  #corr_res[[blk]] = ggMarginal(corr_res[[blk]], type = "density", groupColour = TRUE, groupFill = TRUE)
  }
}

pdf (paste0 ( 'Plots/',your.gene1,'_',your.gene2,'_correlation_scatterplots2.pdf'),20,4)
wrap_plots (corr_res, ncol=3)
#corr_res
dev.off()

# Correlate TFs vs sarcomatoid score (cNMF20) and extract values ####
nfeat = 5000
k_selection = 25
cnmf_spectra_unique = readRDS (paste0('../tumor_compartment/scrna/',paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds')))
activeTFs = read.csv ('../tumor_compartment/scatac_ArchR/Active_TFs.csv')[[2]]


study = c('bueno_rpkm','tcga','mesomics')
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

# # Load TCGA ATAC-seq data (chromvar deviations)
# Match with bulkRNA patient ids to label histology


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





