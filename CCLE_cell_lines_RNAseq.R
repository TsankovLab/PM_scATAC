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

projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/'
dir.create (paste0(projdir, 'Plots/'))
setwd (projdir)
#devtools::install_github("immunogenomics/presto") # needed for DAA
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/useful_functions.R')
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/load_libraries.R')
#source ('/ahg/regevdata/projects/ICA_Lung/Bruno/job_scripts/ArchR_hidden_functions.R')

# Import RNA expression of cell lines from https://portals.broadinstitute.org/ccle/data
if (!file.exists (paste0(projdir,'CCLE_RNAseq_rsem_genes_tpm.rds')))
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
  saveRDS (ccle, paste0(projdir, 'CCLE_RNAseq_rsem_genes_tpm.rds'))
  } else {
  ccle = readRDS (paste0(projdir, 'CCLE_RNAseq_rsem_genes_tpm.rds'))
  }

# get cell lines annotation
meta = read.table (paste0(projdir,'Cell_lines_annotations_20181226.txt'), sep='\t', header=T, quote='',
	comment.char='#',fill = TRUE )
meta = meta[match (colnames(ccle), meta$CCLE_ID),]

# Export list of expression and metadata of CCLE
ccle_list = list(rsem_exp = ccle, meta = meta)
saveRDS (ccle_list, paste0(projdir,'CCLE_rsem_expression_and_metadata.rds'))

# Subset ccle for cell lines ordered
cellLine_ids = c('NCI-H2052','MSTO-211H','NCI-H2452','NCI-H28')
cellLine_ids = c('MSTO-211H','NCI-H2452','ACC-MESO-1')
#meta[grep ('meso',meta$Histology),]


ccle_meso = ccle[,colnames(ccle) %in% meta[meta$Name %in% cellLine_ids,'CCLE_ID']]

# Check TF expression to knock down
activeTFs = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/Active_TFs.csv')

sarc_nmf = readRDS('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/cnmf_genelist_25_nfeat_5000.rds')
sarc_nmf = sarc_nmf[[20]]
#check_genes = c('TWIST2','TCF3','ID3','TEAD1','PITX1','TWIST1','HIC1','PITX2','SOX9','SNAI2','OLIG3','MEF2A','MESP1','SMARCC2')
check_genes = activeTFs[[2]]
check_genes = c(check_genes, 'ITLN1','UPK3B','AXL','VIM','HP','KRT5','CALB2')
cl_order = sapply (1:ncol(ccle_meso), function(x) mean (log2(ccle_meso+1)[sarc_nmf,x], na.rm=T))
names (cl_order) = colnames(ccle_meso)

ccle_meso2 = ccle_meso[rownames (ccle_meso) %in% check_genes, ]
ccle_meso2$gene = rownames (ccle_meso2)
ccle_meso2$MSTO211H_PLEURA = log2 (ccle_meso2$MSTO211H_PLEURA + 1)
ccle_meso2$NCIH2452_PLEURA = log2 (ccle_meso2$NCIH2452_PLEURA + 1)
ccle_meso2$ACCMESO1_PLEURA = log2 (ccle_meso2$ACCMESO1_PLEURA + 1)
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
