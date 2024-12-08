
### Find enriched pathways in genes correlated with selected TFs ####
metacells = readRDS (file.path('..','scrna','metacells.rds'))
selected_TF = readRDS (file.path('..','scatac_ArchR','selected_TF.rds'))
vf = VariableFeatures (FindVariableFeatures (metacells, nfeat=10000))
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)
metacells_assay = metacells_assay[unique(c(vf, selected_TF)),]

enricher_universe = vf
#do.fgsea = TRUE
gmt_annotations = c(
'h.all.v7.4.symbols.gmt'#,#,
#'c5.bp.v7.1.symbol.gmt'
)

sams = unique(metacells$sampleID3)[!unique(metacells$sampleID3) %in% c('HU37','HU62','P11_HOX')] # remove normal samples and outliers 
gmt_annotation = gmt_annotations[1]
force = T
top_genes= 100
if (!file.exists(paste0('EnrichR_activeTF_cor_top_genes_ann_',gmt_annotation,'_top_genes_',top_genes,'.rds')) | force)
	{
	TF_cor_sample = list()
	for (sam in sams)
	  {
	  metacells_assay_sample = as.matrix(metacells_assay[,metacells$sampleID3 == sam])
	  TF_cor_sample[[sam]] = lapply (selected_TF, function(x) 
	  	{
	  	tc_cor = t(cor (metacells_assay_sample[x,], t(metacells_assay_sample)))
	  	tc_cor = setNames (tc_cor[,1], rownames (tc_cor))
	  	tc_cor[is.na(tc_cor)] = 0
	  	tc_cor = tc_cor[order(-tc_cor)]
	  	gmt.file = file.path ('..','..','git_repo','files',gmt_annotation)
	  	pathways = read.gmt (gmt.file)
	  	#pathways = split (pathways$gene, pathways$term)
	    message (paste ('EnrichR running module',x)) 
	    egmt = enricher(head (names(tc_cor),top_genes), TERM2GENE=pathways, universe = enricher_universe)
	    egmt@result
	    })
	   #EnrichRResAll[[ann]] = EnrichRResCluster    
	    # fgseaRes = fgseaMultilevel (pathways, 
		# 			tc_cor#, 
		# 			#minSize=15, 
		# 			#maxSize=1500,
		# 			#BPPARAM = NULL
		# 			)
		# fgseaResCol = collapsePathways (fgseaRes, stats = tc_cor, pathway = pathways)
		# fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
		}
	saveRDS (TF_cor_sample, paste0('EnrichR_activeTF_cor_top_genes_ann_',gmt_annotation,'_top_genes_',top_genes,'.rds'))
	} else {
	TF_cor_sample = readRDS (paste0('EnrichR_activeTF_cor_top_genes_ann_',gmt_annotation,'_top_genes_',top_genes,'.rds'))
	}

gmt.file = file.path ('..','..','git_repo','files',gmt_annotation)
pathways = read.gmt (gmt.file)
pathway_names = as.character (unique(pathways$term))

TF_cor_sample_mat = lapply(TF_cor_sample, function(x)
	{
	pathway_mat = do.call (cbind, lapply (x, function(y) 
		{
		y = y[pathway_names, 'p.adjust', drop=F]
		y$p.adjust = y$p.adjust < 0.05
		}))
	})

TF_cor_sum = Reduce ('+', TF_cor_sample_mat)
rownames (TF_cor_sum) = pathway_names
colnames (TF_cor_sum) = selected_TF
TF_cor_sum[is.na(TF_cor_sum)] = 0
#TF_cor_sum[TF_cor_sum < 3] = 0
TF_cor_sum = TF_cor_sum[rowSums (TF_cor_sum) > 0, ]
#TF_cor_sum = TF_cor_sum[,colSums (TF_cor_sum) > 0 ]

# Export enrichment table ####
saveRDS (TF_cor_sum, 'enrichment_pathways_TFs.rds')
