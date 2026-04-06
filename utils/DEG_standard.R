# Variables to set
#logfcThreshold = .02
#pvalAdjTrheshold = 0.05
#metaGroupName = 'celltype'
#org = 'mouse'
#top_pathways = 5
#force = FALSE

#unlink(".RData") # this is to avoid segfault error which appeared in some instancces in fgsea function
## Find DE markers

# deg_dir = paste0(projdir,'Plots/Clusters_markers_',metaGroupName,'_logFC_',logfcThreshold)
# dir.create (deg_dir)
projdir_DEG = paste0('DEG_',metaGroupName,'_logFC_',logfcThreshold,'_padj_',pvalAdjTrheshold)
dir.create (file.path(projdir_DEG,'Plots'), recursive=T)

# Set the metagroupname as a factor and order levels
if (is.numeric(srt@meta.data[,metaGroupName])) Idents (srt) = factor (srt@meta.data[,metaGroupName], levels =unique(srt@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt@meta.data[,metaGroupName]))))]) else
 Idents (srt) = as.factor (srt@meta.data[,metaGroupName])

if (!file.exists(file.path(projdir_DEG,'DEG.rds')) | force)
	{
	# Run standard findAllMarkers function 
	DefaultAssay(srt) = 'RNA'
	degClusters = FindAllMarkers (srt, max.cells.per.ident = 1000, min.pct = .1, logfc.threshold = logfcThreshold, verbose = T)
	saveRDS(degClusters, file.path(projdir_DEG, 'DEG.rds'))
	write.csv (degClusters[degClusters$p_val_adj < pvalAdjTrheshold,], file.path(projdir_DEG,'DEG.csv'))
	} else {
	message ('DEG file found loading..')	
	degClusters = readRDS (file.path(projdir_DEG, 'DEG.rds'))
	}


message ('Generate Heatmap of top genes') # does not produce image when there are too many genes to plot
degClusters = degClusters[degClusters$p_val_adj < pvalAdjTrheshold, ]

if (length(unique (degClusters$cluster)) > 2) 
	{
	top_deg = degClusters[degClusters$avg_log2FC > 0,]
	top_deg = degClusters %>% group_by (cluster) %>% top_n (top_genes, avg_log2FC)
	#top_deg = top_deg %>% group_by (cluster) %>% top_n (top_genes, -p_val_adj)
	} else {
	top_deg = degClusters[order(degClusters$avg_log2FC),]		
	top_deg = split (top_deg, top_deg$avg_log2FC > 0)
	if (length(top_deg) >1) top_deg = rbind(head(top_deg[[1]], top_genes), tail(top_deg[[2]], top_genes)) else
	top_deg = head (top_deg, top_genes)
	}

top_deg$cluster = factor (top_deg$cluster, levels = names(table (srt@meta.data[,metaGroupName])[order(-table (srt@meta.data[,metaGroupName]))]))
top_deg = top_deg[order(top_deg$cluster),]
# if (!is.na(is.integer(as.integer(srt@meta.data[,metaGroupName][1]))))
#  	{
#  	srt@meta.data[,metaGroupName] = factor (srt@meta.data[,metaGroupName], levels =unique(srt@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt@meta.data[,metaGroupName]))))])
#  	} else {
# 	srt@meta.data[,metaGroupName] = factor (srt@meta.data[,metaGroupName], levels = unique (srt@meta.data[,metaGroupName]))
#  	}


#top_deg = top_deg[order (-top_deg$avg_log2FC),]
srt2 = ScaleData (srt, features = unique(top_deg$gene))

#srt2@meta.data[,metaGroupName] = factor(srt2@meta.data[,metaGroupName], levels = unique (srt2@meta.data[,metaGroupName]))
#srt2 = srt[, !is.na(srt@meta.data[,metaGroupName])]
#if (is.na(as.integer(levels(srt@meta.data[,metaGroupName])))) Idents(srt2) = factor (srt2@meta.data[,metaGroupName], levels =unique(srt2@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt2@meta.data[,metaGroupName]))))])
#if (!is.na(is.integer(as.integer(srt@meta.data[,metaGroupName][1])))) Idents(srt2) = factor (srt2@meta.data[,metaGroupName], levels =unique(srt2@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt2@meta.data[,metaGroupName]))))])
Idents(srt2) = factor (srt@meta.data[,metaGroupName], levels = names(table (srt@meta.data[,metaGroupName])[order(-table (srt@meta.data[,metaGroupName]))]))
heat_p = DoHeatmap (srt2, features = unique(top_deg$gene))
png (file.path(projdir_DEG, 'Plots',paste0('top_',top_genes,'_heatmap.png')), height=6000, width=16300, res=400)
print (heat_p)
dev.off()
pdf (file.path(projdir_DEG, 'Plots',paste0('top_',top_genes,'_heatmap.pdf')), height=6, width=5)
print (heat_p)
dev.off()

rm (srt2)

#### Generate Feature plots and dotplots of top genes ####
message ('Generate Feature plots of top genes')
top_deg = degClusters %>% arrange(-avg_log2FC) %>% group_by (cluster) %>% slice_head (n = top_genes)
#top_deg = top_deg[order (-top_deg$avg_log2FC),]

feat_p = lapply (seq_along(top_deg$gene), function(x) FeaturePlot (srt, features = top_deg$gene[x], keep.scale = 'all', order=F,
  combine = FALSE, pt.size = .01, reduction = reductionName))
feat_p = unlist(feat_p, recursive=F)
for(i in 1:length(feat_p)) {
    feat_p[[i]] = feat_p[[i]] + theme_void() + NoLegend() + NoAxes() + ggtitle (paste0('CL',top_deg$cluster[i],':',top_deg$gene[i])) +
  scale_colour_gradientn (colours = viridis::turbo(100))
    }
  
png (file.path(projdir_DEG,'Plots',paste0('top_',top_genes,'_fplots.png')), ((length(feat_p)*50) + 1500),((length(feat_p)*50) + 1500), res=300)
print (wrap_plots (feat_p))
dev.off()

top_deg = degClusters %>% arrange(-avg_log2FC) %>% group_by (cluster) %>% slice_head (n = top_genes)

dp = DotPlot(object = srt, features = rev(unique(top_deg$gene)), scale = T, group.by = metaGroupName) +
  	theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_line(colour = "gainsboro")) + 
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
dp$data$id = factor (dp$data$id, levels = levels (top_deg$cluster))
png (file.path(projdir_DEG, 'Plots',paste0('top_',top_genes,'_dotplot.png')), ((length(unique(top_deg$gene))*90) + 500), ((length(unique (srt@meta.data[,metaGroupName]))*102) + 250), res=300)
print (dp)
dev.off()
pdf (file.path(projdir_DEG, 'Plots',paste0('top_',top_genes,'_dotplot.pdf')), width=38, height=14)
print (dp)
dev.off()

#### Run GSEA enrichment on each cluster DE markers ####
message ('Run GSEA enrichment on each cluster DE markers')
gmt_annotations = c(
#'c2.cp.kegg.v7.1.symbol.gmt',
'c2.cp.reactome.v7.1.symbol.gmt',
'c5.bp.v7.1.symbol.gmt'
#'h.all.v7.1.symbol.gmt'
)
if (!file.exists (file.path(projdir_DEG, 'Pathway_Enrichment_clusters.rds')) | force)
	{

	# Set background (universe) genes for enrichment analysis
	if (enricher_universe == 'variable_genes') enricher_universe = VariableFeatures(srt)
	if (enricher_universe == 'deg_genes') enricher_universe = unique(degClusters$gene)
	if (enricher_universe == 'all') enricher_universe = rownames(srt)

	# GSEA analysis on DEG per cluster
	EnrichRResAll = list()
	for (ann in gmt_annotations)
		{
		gmt.file = file.path ('/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/GSEA_gs', org, ann)
		pathways = read.gmt (gmt.file)
		#pathways = gmtPathways (gmt.file)
		message (paste('Compute enrichment per cluster using annotation:', ann))
		EnrichRResCluster = list()
		for (i in unique (degClusters$cluster))
			{
			message (paste ('EnrichR running cluster',i))	
			degCluster = degClusters[degClusters$cluster == i,]
			sig_genes = degCluster[degCluster$p_val_adj < pvalAdjTrheshold & degCluster$avg_log2FC > 0, 'gene']
			egmt <- enricher(sig_genes, TERM2GENE=pathways, universe = enricher_universe)
			EnrichRResCluster[[i]] = egmt@result
			}
		EnrichRResAll[[ann]] = EnrichRResCluster
		}
	saveRDS (EnrichRResAll, file.path(projdir_DEG, 'Pathway_Enrichment_clusters.rds'))	
	} else {
	EnrichRResAll = readRDS (file.path(projdir_DEG, 'Pathway_Enrichment_clusters.rds'))
	}

# Plot dotplot of sGSEA annotations per cluster and number of DEG per cluster
cg_df = as.data.frame(table (srt@meta.data[,metaGroupName]))
clustDeg = table (degClusters[degClusters$p_val_adj < pvalAdjTrheshold,'cluster'])
cg_df$sigGenes = clustDeg[match (cg_df$Var1, names(clustDeg))]
colnames (cg_df) = c('cluster','cells','DEG')
cg_df$nCountm = sapply (unique(cg_df$cluster), function(x) mean (srt$nCount_RNA[srt@meta.data[,metaGroupName] == x]))
cg_df$nFeatm = sapply (unique(cg_df$cluster), function(x) mean (srt$nFeature_RNA[srt@meta.data[,metaGroupName] == x]))
cg_df = cg_df[order(-cg_df$DEG),]

umap_p1 = DimPlot (srt, pt.size = .01, label = T, group.by= metaGroupName, reduction = reductionName) + NoLegend()

# Plot fgsea enrichments
EnrichRes_dp = lapply (EnrichRResAll, function(x) dotGSEA (
	enrichmentsTest_list = x, 
	type = 'enrich', 
	padj_threshold = pvalAdjTrheshold, 
	top_pathways= top_pathways,
	cluster_rows=T,
	cluster_cols=T))
for (i in seq_along(gmt_annotations))
	{
	pdf (file.path(projdir_DEG,'Plots', paste0('Pathway_Enrichment_',gmt_annotations[[i]],'_dotplot.pdf')),width = 3 + (length(unique(EnrichRes_dp[[i]]$data$cluster))), height = 3 + (length(unique(EnrichRes_dp[[i]]$data$pathway)) / 7))
	print (EnrichRes_dp[[i]])
	dev.off()
	}

