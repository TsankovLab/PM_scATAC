vg = length(genes.keep)
projdir_SC = file.path(projdir,'SCENIC')
projdir_SC_run = file.path(projdir_SC, paste0('vg_',vg,'_mw_',motif_window),scenic_name)
#dir.create (file.path(projdir_SC, paste0('vg_',vg), 'Plots'), recursive=T)

# Import results
auc_mtx <- read.csv(file.path(projdir_SC_run, 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = paste0('SCENIC_',colnames(auc_mtx))
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

# Make heatmap avering AUC scores
auc_mtx_avg = aggregate (auc_mtx, by=as.list(srt@meta.data[,metaGroupNames[2],drop=F]), mean)
rownames (auc_mtx_avg) = auc_mtx_avg[,1]
auc_mtx_avg = auc_mtx_avg[,-1]
colnames (auc_mtx_avg) = sub ('SCENIC_','', colnames(auc_mtx_avg))
palette_scenic = rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(50))

if (nrow (auc_mtx_avg) > 1) auc_mtx_avg_scaled = scale (auc_mtx_avg) else
auc_mtx_avg_scaled = auc_mtx_avg
hm = Heatmap (auc_mtx_avg_scaled,
 col = palette_scenic, column_names_rot = 45,
 column_names_gp = gpar(fontsize = 7),
 row_names_gp = gpar(fontsize = 8))
pdf (file.path(projdir_SC_run,'Plots','auc_scores_celltypes_heatmap.pdf'), width=12, height=2.2)
print (hm)
dev.off()

# srt@meta.data = srt@meta.data[,!colnames(srt@meta.data) %in% colnames(auc_mtx)]
srt = AddMetaData (srt, auc_mtx)
SCENIC_mods = colnames(srt@meta.data)[grepl ('SCENIC', colnames(srt@meta.data))]
umap_p1 = FeaturePlot (srt, 
    features = SCENIC_mods, 
    ncol=4, 
    combine = FALSE, 
    pt.size = .01, 
    reduction = reductionName)
for(i in 1:length(umap_p1)) 
  {
  umap_p1[[i]] = umap_p1[[i]] + 
    theme_void() + 
    NoAxes() +
    ggtitle (SCENIC_mods[[i]]) + 
    theme(plot.title = element_text(size = 8)) +
    #scale_colour_gradientn (colours = rev(brewer.pal(n = 11, name = "RdBu")),limits=c(-max (abs(p[[i]]$data[,4])), max (abs(p[[i]]$data[,4]))))
    scale_colour_gradientn (colours = viridis::viridis(100),limits=c(0, max (umap_p1[[i]]$data[,4])))
  }

## read motifs table ####
motifs_table = read.csv (file.path(projdir_SC_run, 'motifs.csv'), skip=2) 
colnames (motifs_table) = c('TF','MotifID','AUC','NES','MotifSimilarity','OrthologousIdentity','Annotation','Context','TargetGenes','RankAtMax')
tgenes = lapply (seq(motifs_table$TargetGenes), function(x)
  {
  tmp = strsplit (motifs_table$TargetGenes[x], '\\)\\, \\(')[[1]]  
  tmp = sub ('\\[\\(','',tmp)
  tmp = sub ('\\)\\]','', tmp)
  gene = sapply (tmp, function(y) unlist(strsplit (y, '\\,'))[[1]])
  gene = gsub ("\\'","",gene)
  score = sapply (tmp, function(y) unlist(strsplit (y, '\\,'))[[2]])
  score = as.numeric (score)
  data.frame (row.names = gene, gene = gene, score = score, TF= motifs_table$TF[x])
  })
names (tgenes) = paste0(motifs_table$TF,'_',motifs_table$MotifID)
tgenes = do.call (rbind,tgenes)
tgenes = split (tgenes, tgenes$TF)
tgenes = lapply (tgenes, function(x) x[!duplicated(x$gene),])

# pdf (file.path(projdir_SC_run,'Plots','SCENIC_modules_umaps.pdf'),20,20)
# print (wrap_plots (umap_p1))
# dev.off()

# box_p = lapply (seq_along(SCENIC_mods), function(x) 
#       ggplot (ccomp_df, aes_string (x= metaGroupNames[1], y= meta_modules_names[x], fill = metaGroupNames[1])) +
#             #geom_violin (trim=TRUE) +
#             geom_violin () + 
#             geom_boxplot () +
#             #geom_jitter (color="black", size=0.4, alpha=0.9) +
#             theme_classic() + 
#             scale_fill_manual (values= module_pal) + 
#             ggtitle (meta_modules_names[x]) + 
#             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))) 
#     if (length(unique(metaGroupNames)) == 3) box_p = lapply (seq_along(meta_modules_names), function(x) box_p[[x]] + facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + theme_classic()+ 
#             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))) 
    

### Generate boxplots per meta groups
ccomp_df = srt@meta.data[,SCENIC_mods, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,metaGroupNames,drop=F]), mean)
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
  
#   box_p = lapply (seq(SCENIC_mods), function(x) 
#   {
# 	sample_rep = srt@meta.data[!duplicated (srt@meta.data[,metaGroupNames[1]]),]	
# 	sample_rep = table (sample_rep[,metaGroupNames[2]], sample_rep[,metaGroupNames[3]])
# 	keep_group = apply (sample_rep , 1, function(x) !any (x < 2))
# 	keep_group = names(keep_group[keep_group])
# 	if (length(keep_group) > 0)
# 		{
#     box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods[x])) +
#           #geom_violin (trim=TRUE) +
#           geom_boxplot (aes_string(fill = metaGroupNames[3])) +
#           #geom_jitter (color="black", size=0.4, alpha=0.9) +
#           theme_classic() + 
#           #scale_fill_manual (values= module_pal) + 
#           ggtitle (SCENIC_mods[x]) + 
#           theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# 		stat.test = box$data[,box$data[,metaGroupNames[2]] %in% keep_group] %>%
# 			  group_by (metaGroupNames[2]) %>%
# 			  t_test(reformulate (metaGroupNames[3], SCENIC_mods[x])) %>%
# 			  adjust_pvalue (method = "none") %>%
# 			  add_significance ()
# 			stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.01)
#   		box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
#    		bracket.nudge.y = 0, hide.ns = F,
#     	label = "p.adj.signif") + NoLegend()
# 		} else {
# 		ccomp_df = srt@meta.data[,c(SCENIC_mods, metaGroupNames), drop=FALSE]
#     box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods[x])) +
#         geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[3])) +
#         #geom_bar (stats='identity') +
#         #geom_jitter (color="black", size=0.4, alpha=0.9) +
#         theme_classic() + 
#         #scale_fill_manual (values= module_pal) + 
#         ggtitle (SCENIC_mods[x]) + 
#         theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#         stat.test = box$data %>%
#         group_by (metaGroupNames[2]) %>%
#         wilcox_test(reformulate (metaGroupNames[3], SCENIC_mods[x])) %>%
#         adjust_pvalue (method = "none") %>%
#         add_significance ()
#         stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.01)
#         box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
#         bracket.nudge.y = 0, hide.ns = F,
#         label = "p.adj.signif") + NoLegend()
# 		}
#     })
# if (length(unique(ccomp_df[,metaGroupNames[[2]]])) > 1) box_p = lapply (seq_along(SCENIC_mods), function(x) box_p[[x]] + facet_wrap (as.formula(paste("~", metaGroupNames[2]))) + theme_classic()+ 
#         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + NoLegend())
		
  
# png (file.path (projdir_SC_run,'Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.png')), width = 400 + (length(umap_p1) * 200), height = 400 + (length(umap_p1) * 50), pointsize=10, res = 300, type="cairo")
# print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
#   wrap_plots (box_p, ncol= ifelse (length(box_p) > 5,ceiling(length(box_p)/2),length(box_p))) + plot_layout (heights=c(1,2)))
# dev.off()


box_p = lapply (seq(SCENIC_mods), function(x) 
    {
    sample_rep = srt@meta.data
    sample_rep = table (sample_rep[,metaGroupNames[2]], sample_rep[,metaGroupNames[3]])
    keep_group = apply (sample_rep , 1, function(y) !any (y < 2))
    keep_group = names(keep_group[keep_group])
      if (length(keep_group) > 0 & metaGroupNames[1] != metaGroupNames[3])
      {
      box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods[x])) +
          #geom_violin (trim=TRUE) +
          geom_boxplot (aes_string(fill = metaGroupNames[3])) +
          #geom_jitter (color="black", size=0.4, alpha=0.9) +
          theme_classic() + 
          #scale_fill_manual (values= module_pal) + 
          ggtitle (SCENIC_mods[x]) + 
          theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
      stat.test = box$data[box$data[,metaGroupNames[2]] %in% keep_group,] %>%
        group_by_at (metaGroupNames[2]) %>%
        t_test(reformulate (metaGroupNames[3], SCENIC_mods[x])) %>%
        adjust_pvalue (method = "none") %>%
        add_significance ()
        stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
        box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
        bracket.nudge.y = 0, hide.ns = T,
        label = "p.adj.signif") + NoLegend()
      } else {
      ccomp_df = srt@meta.data[,c(SCENIC_mods, metaGroupNames), drop=FALSE]
      #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
      box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods[x])) +
        geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[3])) +
        geom_boxplot (aes_string(fill = metaGroupNames[3])) +
        #geom_bar (stats='identity') +
        #geom_jitter (color="black", size=0.4, alpha=0.9) +
        theme_classic() + 
        #scale_fill_manual (values= module_pal) + 
        ggtitle (SCENIC_mods[x]) + 
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        stat.test = box$data %>%
        group_by_at (metaGroupNames[2]) %>%
        wilcox_test(reformulate (metaGroupNames[3], SCENIC_mods[x])) %>%
        adjust_pvalue (method = "none") %>%
        add_significance ()
        stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
        
        box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
        bracket.nudge.y = 0, hide.ns = T,
        label = "p.adj.signif") + NoLegend()
      }})

  if (length(unique(ccomp_df[,metaGroupNames[2]])) > 1) box_p = lapply (seq_along(SCENIC_mods), function(x) box_p[[x]] + facet_wrap (as.formula(paste("~", metaGroupNames[2]))) + theme_classic() + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + NoLegend()) 
  if (exists ('module_pal')) box_p = lapply (seq_along(SCENIC_mods), function(x) box_p[[x]] + scale_fill_manual (values= module_pal)) 
      
      
    png (file.path (projdir_SC_run, 'Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.png')), width = 400 + (length(umap_p1) * 400), height = 400 + (length(umap_p1) * 100), pointsize=10, res = 300, type="cairo")
    print (wrap_plots (umap_p1))      
    dev.off()
      
    pdf (file.path (projdir_SC_run, 'Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.pdf')), width = 80+length(umap_p1), height = 28)
    print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
      wrap_plots (box_p, ncol= ifelse (length(box_p) > 5,ceiling(length(box_p)/2),length(box_p))) + plot_layout (heights=c(1,2)))
    dev.off()
      
  
  
