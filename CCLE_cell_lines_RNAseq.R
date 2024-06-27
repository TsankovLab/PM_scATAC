use UGER
conda activate meso_scatac

# projdir
# mkdir -p Public_data/CCLE
# cd Public_data/CCLE
# wget https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
# wget https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz
# gzip -d CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
# gzip -d CCLE_RNAseq_genes_rpkm_20180929.gct.gz

R

# Load packages

projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/CCLE/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd(projdir)

#devtools::install_github("immunogenomics/presto") # needed for DAA
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/useful_functions.R')
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/load_libraries.R')
#source ('/ahg/regevdata/projects/ICA_Lung/Bruno/job_scripts/ArchR_hidden_functions.R')

# Import RNA expression of cell lines from https://portals.broadinstitute.org/ccle/data
if (!file.exists (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/CCLE_RNAseq_rsem_genes_tpm.rds')))
  {
  ccle = read.table (paste0(projdir,'CCLE_RNAseq_rsem_genes_tpm_20180929.txt'), sep='\t',header=T)
  library(biomaRt)
  mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  library (stringr)
  gene_ids = str_replace (ccle$gene_id,
                          pattern = ".[0-9]+$",
                          replacement = "")
  eID_symbol = getBM (attributes = c('ensembl_gene_id',
                       'hgnc_symbol'),
        filters = 'ensembl_gene_id', 
        values = gene_ids,
        mart = mart)
  rownames (ccle) = gene_ids
  ccle$symbol = eID_symbol$hgnc_symbol[match (rownames(ccle), eID_symbol$ensembl_gene_id)]
  ccle = ccle [,-c(1:2)]
  ccle = aggregate (.~ symbol, data = ccle, FUN=mean)
  ccle = ccle[!ccle$symbol == "",]
  rownames (ccle) = ccle[,1]
  ccle = ccle[,-1]
  colnames (ccle) = gsub ('^X','', colnames (ccle))
  saveRDS (ccle, paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/CCLE_RNAseq_rsem_genes_tpm.rds'))
  } else {
  ccle = readRDS (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/CCLE_RNAseq_rsem_genes_tpm.rds'))
  }

# get cell lines annotation ####
meta = read.table (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/Cell_lines_annotations_20181226.txt'), sep='\t', header=T, quote='',
	comment.char='#',fill = TRUE )
meta = meta[match (colnames(ccle), meta$CCLE_ID),]

# Export list of expression and metadata of CCLE ####
# ccle_list = list(rsem_exp = ccle, meta = meta)
# saveRDS (ccle_list, paste0('CCLE_rsem_expression_and_metadata.rds'))

# Subset ccle for cell lines ordered ####
cellLine_ids = c('NCI-H2052','MSTO-211H','NCI-H2452','NCI-H28')
cellLine_ids = c('MSTO-211H','NCI-H2452','ACC-MESO-1')
ccle_meso = ccle[,colnames(ccle) %in% meta[meta$Name %in% cellLine_ids,'CCLE_ID']]

# Log2 expression ####
ccle_meso$MSTO211H_PLEURA = log2 (ccle_meso$MSTO211H_PLEURA + 1)
ccle_meso$NCIH2452_PLEURA = log2 (ccle_meso$NCIH2452_PLEURA + 1)
ccle_meso$ACCMESO1_PLEURA = log2 (ccle_meso$ACCMESO1_PLEURA + 1)

# Check TF expression to knock down ####
activeTFs = read.csv ('../tumor_compartment/scatac_ArchR/Active_TFs.csv')

sarc_nmf = readRDS('../tumor_compartment/scrna/cnmf_genelist_25_nfeat_5000.rds')
sarc_nmf = sarc_nmf[[20]]
#check_genes = c('TWIST2','TCF3','ID3','TEAD1','PITX1','TWIST1','HIC1','PITX2','SOX9','SNAI2','OLIG3','MEF2A','MESP1','SMARCC2')
check_genes = activeTFs[[2]]
check_genes = c(check_genes, 'ITLN1','UPK3B','AXL','VIM','HP','KRT5','CALB2')
cl_order = sapply (1:ncol(ccle_meso), function(x) mean (log2(ccle_meso+1)[sarc_nmf,x], na.rm=T))
names (cl_order) = colnames(ccle_meso)

ccle_meso2 = ccle_meso[rownames (ccle_meso) %in% check_genes, ]
ccle_meso2$gene = rownames (ccle_meso2)



write.csv (ccle_meso2, 'meso_scATAC_active_TFS.csv')
ccle_meso2 = gather (ccle_meso2, cell_line, expression, 1:(ncol(ccle_meso2)-1))
ccle_meso2$expression =  (ccle_meso2$expression)
ccle_meso2$cell_line = factor (ccle_meso2$cell_line, levels = names(cl_order[order(-cl_order)]))
mean_exp = do.call (rbind, lapply (split (ccle_meso2, ccle_meso2$gene), function(x) mean(x$expression)))
ccle_meso2$gene = factor (ccle_meso2$gene, levels = rownames(mean_exp)[order(-mean_exp[,1])])


ccle_p = ggplot(ccle_meso2, 
    aes(fill=cell_line, y=expression, x=gene)) +
    geom_bar(position="dodge", stat="identity") + 
#    scale_fill_viridis (discrete = T)  + 
    theme_classic() +
    ggtitle ("genes to check") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf (paste0(projdir, 'Plots/CCLE_tf_expression_barplots.pdf'), width=14, height=3)
wrap_plots (ccle_p, ncol=1)
dev.off()

# Map archetypes and selected TF to cell line expression and export table ####
# Import Viestra archetypes to add to the table below ####
motif_clusters = read.csv ('../Viestra_motif_clustered.csv')
colnames (motif_clusters)[1] = 'Cluster'
head (motif_clusters)
motif_clusters$Cluster2 = sapply (motif_clusters$Motif, function(x) paste(unique(motif_clusters[motif_clusters$Motif == x,'Cluster']),collapse=','))
selected_TFs = read.csv ('../tumor_compartment/scatac_ArchR/Active_TFs.csv')$x

ccle_meso2 = ccle_meso[rownames (ccle_meso) %in% motif_clusters$Motif, ]
head (ccle_meso2)
ccle_meso2$hit = selected_TFs[match(rownames(ccle_meso2), selected_TFs)]
ccle_meso2$archetype = motif_clusters$Cluster2 [match(rownames(ccle_meso2), motif_clusters$Motif)]
ccle_meso2 = ccle_meso2[order(ccle_meso2$archetype), ]
write.csv (ccle_meso2, 'ccle_with_archetype.csv')

#### Check chromatin regulators and genetic drivers genes ####
chromatin_regulators = c('SETD5, ASH1L, CREBBP, PRDM2, KDM2B, KMT2D, EZH2, SETDB1')
chromatin_regulators = unlist(strsplit (chromatin_regulators, ', '))
genetic_drivers = c('BAP1, NF2, CDKN2B, CDKN2A, TP53, LATS2, SETD2, FAT4, PTCH1')
genetic_drivers = unlist(strsplit (genetic_drivers, ', '))
chromatin_regulators2 = c('LATS1, DDX3X, ULK2, RYR2, CFAP45, SETDB1, DDX51, SF3B1, TRAF7, PTEN, RBFOX1, CSMD1, MTAP, TTC28, PCDH15, USH2A, CNTNAP2, DNAH1, KCNH7, PTK2, ROBO2, DLG2, PBRM1, PTPRD, ANTXR2, CTNNA3, LINGO2, LRP1B, PLCB1, UNC79, WWOX')
chromatin_regulators2 = unlist (strsplit (chromatin_regulators2, ', '))
genetic_drivers = unlist(strsplit (genetic_drivers, ', '))
check_genes = unique (c(chromatin_regulators, genetic_drivers, chromatin_regulators2))

ccle_meso2 = ccle_meso[rownames (ccle_meso) %in% check_genes, ]
ccle_meso2$gene = rownames (ccle_meso2)

write.csv (ccle_meso2, 'CCLE_genetic_Drivers_and_chromatin_regulators.csv')

ccle_meso2 = gather (ccle_meso2, cell_line, expression, 1:(ncol(ccle_meso2)-1))
ccle_meso2$expression =  (ccle_meso2$expression)
ccle_meso2$cell_line = factor (ccle_meso2$cell_line, levels = names(cl_order[order(-cl_order)]))
mean_exp = do.call (rbind, lapply (split (ccle_meso2, ccle_meso2$gene), function(x) mean(x$expression)))
ccle_meso2$gene = factor (ccle_meso2$gene, levels = rownames(mean_exp)[order(-mean_exp[,1])])


ccle_p = ggplot(ccle_meso2, 
    aes(fill=cell_line, y=expression, x=gene)) +
    geom_bar(position="dodge", stat="identity") + 
#    scale_fill_viridis (discrete = T)  + 
    theme_classic() +
    ggtitle ("genes to check") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf (paste0('Plots/CCLE_genetic_and_chromatin_drivers_barplots.pdf'), width=14, height=3)
wrap_plots (ccle_p, ncol=1)
dev.off()

