gmt_annotations = c(
#'c2.cp.kegg.v7.1.symbol.gmt',
'fetal_sigs.gmt'
#'c2.cp.reactome.v7.1.symbol.gmt',
#'c5.bp.v7.1.symbol.gmt',
#'h.all.v7.1.symbol.gmt'
)
force = FALSE
org = 'human'
pvalAdjTrheshold = 0.05
	
	mainPathways = list()	
	for (ann in gmt_annotations)
		{
		if (!file.exists (paste0('fGSEA_annotation_',ann,'.rds')) | force)
			{
			gmt.file = paste0 ('/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/GSEA_gs/',org,'/',ann)
			pathways = gmtPathways (gmt.file)
			pathways = lapply (pathways, function(x) x[!is.na(x)])
			pathways = lapply (pathways, function(x) x[!x == 'NA'])
			pathways = lapply (pathways, function(x) x[!x == ''])
			message ('fGSEA running..')
			
			# Check if list of multiple ranked vectors or simple vector ####
			if (is.list(ranked_vector))
			{
			mainPathways[[ann]] = lapply (ranked_vector, function(x) 
				{
				message ('after NA removed: ', length(x), ' genes')
				x = x[!is.na(x)]
				fgseaRes = fgseaMultilevel (pathways, 
					x, 
					minSize=15, 
					maxSize=500, # changed this from 1500 to 1000 cause it generated an error
					BPPARAM = NULL)
				fgseaResCol = collapsePathways (fgseaRes, stats = ranked_vector, pathway = pathways)
				fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
				})
			saveRDS (mainPathways, file.path(paste0('fGSEA_annotation_',ann,'.rds')))
			
			# Generate dotplot of enrichment results ####
			top_pathways = 10
			fgseaResAll_dp = lapply (mainPathways[[ann]], function(y) dotGSEA (y, padj_threshold = pvalAdjTrheshold, 
			type = 'fgsea',top_pathways = top_pathways,
			cluster_rows=T,
			cluster_cols=T)
			)

			pdf (file.path('Plots',paste0('fGSEA_enrichment_results_gmt_annotation_',ann,'_dotplot.pdf')), width = 14, height=5)
			print(fgseaResAll_dp)
			dev.off()

			} else {
			ranked_vector = ranked_vector[!is.na(ranked_vector)]
			message ('after NA removed: ', length(ranked_vector), ' genes')
			fgseaRes = fgseaMultilevel (pathways, 
				ranked_vector, 
				minSize=15, 
				maxSize=500, # changed this from 1500 to 1000 cause it generated an error
				BPPARAM = NULL)
			fgseaResCol = collapsePathways (fgseaRes, stats = ranked_vector, pathway = pathways)
			mainPathways[[ann]] = fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
			saveRDS (mainPathways, file.path(paste0('fGSEA_annotation_',paste(gmt_annotations,collapse='_'),'.rds')))
			
			top_pathways = 50
			fgseaRes_flt = mainPathways[[1]]
			fgseaRes_flt = split (fgseaRes_flt, fgseaRes_flt$NES > 0)
			fgseaRes_flt = lapply (fgseaRes_flt, function(x) head (x, top_pathways))
			fgseaRes_flt = do.call (rbind, fgseaRes_flt)
			fgseaRes_flt$padj_log_signed = -log10 (fgseaRes_flt$padj) * sign (fgseaRes_flt$NES)
			fgseaRes_flt = fgseaRes_flt[order (fgseaRes_flt$NES),]
			fgseaRes_flt$pathway = factor (fgseaRes_flt$pathway, levels = unique (fgseaRes_flt$pathway))
			fgseaRes_flt$direction = as.character(sign(fgseaRes_flt$NES))
			#fgseaRes_flt$direction = ifelse (fgseaRes_flt$direction == '1', muscatIdents[1], muscatIdents[2])
			fgseaRes_flt = fgseaRes_flt[fgseaRes_flt$padj < pvalAdjTrheshold,]
			fgseaRes_flt$pathway = factor (fgseaRes_flt$pathway, levels = fgseaRes_flt$pathway[order (fgseaRes_flt$padj_log_signed)])
			fgseaResAll_dp = ggplot (fgseaRes_flt, aes (x = padj_log_signed, y = pathway, fill = direction)) + 
			geom_bar (stat = 'identity') +
			theme_classic () +
			scale_fill_manual(values = c("red", "seagreen3")) +
			theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
			
			pdf (file.path('Plots',paste0('fGSEA_enrichment_results_gmt_annotations_',ann,'_barplots.pdf')), width = 14, height=15)
			print(fgseaResAll_dp)
			dev.off()
			}	 
		

		
	
	} else {
	message ('fGSEA object found!')	
	mainPathways = readRDS (file.path(paste0('fGSEA_annotation_',ann,'.rds')))
		
	mainPathways_print = lapply (seq_along(mainPathways), function(y) 
		{	
		tmp2 = mainPathways[[y]][mainPathways[[y]]$padj < pvalAdjTrheshold,]
		tmp2 = tmp2[order(tmp2$padj),]
		tmp2$cluster = names (mainPathways[y])
		tmp2
		})
	mainPathways_print = do.call (rbind, mainPathways_print)
	mainPathways_print = apply(mainPathways_print,2,as.character)
	mainPathways_print = mainPathways_print[,c(ncol(mainPathways_print), 1:(ncol(mainPathways_print)-1))]
	write.csv (mainPathways_print, file.path(paste0('fgsea_enrichments_',ann,'_table.csv')))
	
	# Plot dotplot of fGSEA annotations per cluster 
	lapply (seq_along(fgseaResAll_dp), function(x) {
	pdf (file.path('Plots',paste0('fGSEA_top_',top_pathways,'_annotation_',names (fgseaResAll_dp)[x],'_dotplots.pdf')),15,8)
		print(fgseaResAll_dp[[x]])
		dev.off()
		})
	}
}
	#if (feat == 'all')
	#	{
	#	dhm_grob = grid.grabExpr(draw(dhm[[1]]))
	#	pdf (paste0(projdir,'Plots/', deg2ClusterFileName,'_top_',topGenes,'_genes_heatmap_and_fgsea_dotplot.pdf'), 10, 8)
	#	for (i in seq_along(fgseaResAll_dp))
	#		{
	#		print (dhm_grob + fgseaResAll_dp[[i]])
	#		}
	#	dev.off()			
	#	}



