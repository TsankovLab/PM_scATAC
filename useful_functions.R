
### Useful functions ###



### Replace hdWGCNA function fixing bug ####

ModuleNetworkPlot <- function(
  seurat_obj,
  mods="all",
  outdir="ModuleNetworks",
  plot_size = c(6,6),
  wgcna_name=NULL,
  label_center = FALSE, # only label the genes in the middle?
  edge.alpha=0.25,
  vertex.label.cex=1,
  vertex.size=6, ...
){

  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # get modules, MEs:
  MEs <- GetMEs(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)



  # using all modules?
  if(any(mods == 'all')){
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  }

  # check if we have eigengene-based connectivities:
  if(!all(paste0('kME_', as.character(mods)) %in% colnames(modules))){
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }

  # create output folder
  if(!dir.exists(outdir)){dir.create(outdir)}

  # tell the user that the output is going to the output dir
  cat(paste0("Writing output files to ", outdir))

  # get TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)

  # get hub genes:
  n_hubs <- 25
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs)
    colnames(cur)[2] <- 'var'
    cur %>% arrange(desc(var)) %>% .$gene_name
  })
  names(hub_list) <- mods

  # loop over modules
  for(cur_mod in mods){
    print(cur_mod)
    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique

    # number of genes, connections
    # might make a setting to change this later but I will also have to change
    # how the graph layout works
    n_genes = 25;
    n_conns = 500;

    # name of column with current kME info
    cur_kME <- paste0('kME_', cur_mod)

    cur_genes <- hub_list[[cur_mod]]

    # Identify the columns in the TOM that correspond to these hub genes
    matchind <- match(cur_genes, colnames(TOM))
    reducedTOM = TOM[matchind,matchind]
    orderind <- order(reducedTOM,decreasing=TRUE)

    # only  keep top connections
    connections2keep <- orderind[1:n_conns];
    reducedTOM <- matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] <- 1;

    # print('here')
    # print(dim(reducedTOM))
    # print(n_genes)

    # only label the top 10 genes?
    if(label_center){cur_genes[11:25] <- ''}

    # top 10 as center
    gA <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
    gB <- graph.adjacency(as.matrix(reducedTOM[11:n_genes,11:n_genes]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))

    g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)

    pdf(paste0(outdir, '/', cur_mod,'.pdf'), width=plot_size[1], height=plot_size[2], useDingbats=FALSE);
    plot(g1,
      edge.color=adjustcolor(cur_color, alpha.f=0.25),
      edge.alpha=edge.alpha,
      vertex.color=cur_color,
      vertex.label=as.character(cur_genes),
      vertex.label.dist=1.1,
      vertex.label.degree=-pi/4,
      vertex.label.color="black",
      vertex.label.family='Helvetica',
      vertex.label.font = 3,
      vertex.label.cex=vertex.label.cex,
      vertex.frame.color='black',
      layout= jitter(layoutCircle),
      vertex.size=vertex.size,
      main=paste(cur_mod)
    )
    dev.off();

  }

}


# Save seurat object
sv = function() {
	saveRDS (srt, paste0(projdir,'srt.rds'))}

# Compute cell states in GBM using the gene sets and method from Neftel et al paper
Nefstates = function (seurat_obj = srt, mouse=T)
    {
    require (Seurat)    
    message ('Read csv and compute module scores')
    if (mouse) gbm_modules_6states = readRDS (paste0 ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/mouse_GBM_netfel_6_states.rds')) else
    gbm_modules_6states = readRDS (paste0 ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/human_GBM_netfel_6_states.rds'))        
    
    seurat_obj = ModScoreCor (
        seurat_obj = seurat_obj, 
        geneset_list = gbm_modules_6states, 
        cor_threshold = NULL, 
        pos_threshold = -.03,
        listName = 'neftel6', 
        outdir = NULL)
    seurat_obj_df = as.data.frame (seurat_obj@meta.data)
    seurat_obj_df$MESlike = apply (seurat_obj_df[,c('MES1','MES2')],1,max)
    seurat_obj_df$NPClike = apply (seurat_obj_df[,c('NPC1','NPC2')],1,max)
    
    message ('compute y axis')
    max_AC_or_MES = pmax(seurat_obj_df$AC, seurat_obj_df$MESlike)
    max_OPC_or_NPC = pmax(seurat_obj_df$OPC, seurat_obj_df$NPClike)
    seurat_obj_df$y_axis <- log2(abs(max_OPC_or_NPC - max_AC_or_MES) + 1)
    seurat_obj_df$y_axis[max_AC_or_MES > max_OPC_or_NPC] <- -1 * seurat_obj_df$y_axis[max_AC_or_MES > max_OPC_or_NPC]

    message ('compute x axis')
    seurat_obj_df$x_axis = 0
    seurat_obj_df$x_axis[seurat_obj_df$y_axis > 0] = log2(abs(seurat_obj_df$OPC - seurat_obj_df$NPClike) + 1)[seurat_obj_df$y_axis > 0]
    seurat_obj_df$x_axis[seurat_obj_df$y_axis > 0 & seurat_obj_df$OPC > seurat_obj_df$NPClike] <- -1 * seurat_obj_df$x_axis[seurat_obj_df$y_axis > 0 & seurat_obj_df$OPC > seurat_obj_df$NPClike]
    
    seurat_obj_df$x_axis[seurat_obj_df$y_axis < 0] = log2(abs(seurat_obj_df$AC - seurat_obj_df$MESlike) + 1)[seurat_obj_df$y_axis < 0]
    seurat_obj_df$x_axis[seurat_obj_df$y_axis < 0 & seurat_obj_df$AC > seurat_obj_df$MESlike] <- -1 * seurat_obj_df$x_axis[seurat_obj_df$y_axis < 0 & seurat_obj_df$AC > seurat_obj_df$MESlike]
        
    message ('assign states to cells')
    seurat_obj_df$Tstate=0
    seurat_obj_df$Tstate[seurat_obj_df$x_axis < 0 & seurat_obj_df$y_axis >= 0] = 'OPC-like' 
    seurat_obj_df$Tstate[seurat_obj_df$x_axis >= 0 & seurat_obj_df$y_axis > 0] = 'NPC-like' 
    seurat_obj_df$Tstate[seurat_obj_df$x_axis < 0 & seurat_obj_df$y_axis < 0] = 'AC-like' 
    seurat_obj_df$Tstate[seurat_obj_df$x_axis > 0 & seurat_obj_df$y_axis < 0] = 'MES-like'  
    return (seurat_obj_df)
    }


# Check fraction of positive cells for any number of genes
is.gene.expressed = function(
  genes=NULL, 
  mat_counts=NULL,
  metadata = srt@meta.data, 
  metaGroupNames=NULL, 
  min_expression=0)
  {
  mat_counts[is.na(mat_counts)] = 0
  isExp_mat = mat_counts[genes,,drop=F] > min_expression
  isExp_df = as.data.frame (t(isExp_mat))
  isExp_df = cbind (isExp_df, metadata[,metaGroupNames,drop=F])
  res = lapply (genes, function(x) 
    {
    cellComp (
    seurat_obj = isExp_df, 
    metaGroups = c(metaGroupNames, x),
    ptable_factor = 1,
    #plot_as = 'box',
    prop = TRUE,
    #pal = ,
    #facet_ncol = 15
    returnDF = TRUE,
    removeNA = FALSE
    )
    })
  return (do.call (cbind, res))
  #res = res[,!duplicated(colnames(res))]
  }


# Check if vector is unique or not
unq = function(x) 
    {
    ifelse (length(x) == length(unique(x)), return ('unique'), return ('non-unique'))
    }


# wrapper to generate fast feature plots 
fp = function (seurat_obj, gene, reduction = reductionName) 
	{
	require (viridis)
	require (Seurat)	
	p = FeaturePlot (seurat_obj, features = gene, keep.scale = 'all', order=F,
  combine = FALSE, pt.size = .01, reduction = reductionName)
	lapply(p, function(x) 
    x + 
    theme_void() + 
    NoAxes() + 
    ggtitle (colnames(x$data)[4]) +
    scale_colour_gradientn (colours = viridis::turbo(100)))
  }



# Covnert a named list of length of different vectors in a data.frame.
patchvecs = function(vectorList)
	{
	maxLength = max (sapply (vectorList, length))
	tmpL = lapply (vectorList, function(x) c(x, rep (NA, maxLength - length(x))))	
	df = as.data.frame (do.call (cbind,tmpL))
	colnames (df) = names (vectorList)
	return (df)
	}

### Function to create module score of a list gene sets (named), generate a metagroup of highest score across genesets per cell and optionally filter it based on co-expression ###
ModScoreCor = function (seurat_obj, geneset_list, listName, cor_threshold = NULL, pos_threshold = .1, outdir)
        {        
        require (Seurat)	
        message ('Run AddModuleScore')
        seurat_obj = AddModuleScore (seurat_obj, geneset_list)
        seurat_obj@meta.data = seurat_obj@meta.data[, !colnames (seurat_obj@meta.data) %in% names (geneset_list)]
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% paste0('Cluster',seq_along(geneset_list))] = names (geneset_list)
        message (paste('Annotate cells based on highest module score and store in column:',paste0(listName, '_r',cor_threshold,'_max')))
        if (length (geneset_list) == 1) 
        	{
        	seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = ifelse (seurat_obj@meta.data[,names (geneset_list)] > pos_threshold, 'pos','neg')
        	pdf (paste0(outdir, listName, '_modulescore_distribution_cor_threshold_',cor_threshold,'_score_',pos_threshold,'.pdf'))
        	hist (seurat_obj@meta.data[,names(geneset_list)])
          abline (v = pos_threshold)
          dev.off()
        	} else {
        	seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = sapply (seq_along(colnames(seurat_obj)), function(x) colnames(seurat_obj@meta.data[,names(geneset_list)])[which.max (seurat_obj@meta.data[x,names(geneset_list)])])        		
        	}
        if (!is.null(cor_threshold))
                {
                message ('cor_threshold provided! Filtering gene sets based on initial correlation to module score')	
                filtered_geneset_list = list()
                geneset_cor_list = list()
                for (i in names(geneset_list))
                        {       
                        geneset_cor = cor (seurat_obj@meta.data[,i], as.matrix(t(seurat_obj@assays$RNA@data[rownames(seurat_obj@assays$RNA@data) %in% geneset_list[[i]],])))
                        geneset_cor_list[[i]] = geneset_cor
                        geneset_cor_names = colnames (geneset_cor)[geneset_cor > cor_threshold]
                        geneset_cor_names = geneset_cor_names[!is.na (geneset_cor_names)]
                        filtered_geneset_list[[i]] = geneset_cor_names
                        }
                if (!is.null (outdir)) 
                        {
                        lapply (seq_along(filtered_geneset_list), function(x) write.csv (filtered_geneset_list[[x]], paste0(outdir,'Corfiltered_Module_score_gene_list_', names(filtered_geneset_list)[x],'.csv')))
                        pdf (paste0(outdir, listName, 'Corfiltered_modulescore_distribution.pdf'))
                        lapply (seq_along(filtered_geneset_list), function(x) 
                                {
                                hist (geneset_cor_list[[x]], title = names(geneset_cor_list)[x])
                                abline (v = cor_threshold)
                                })
                        dev.off()
                        }
                message ('Re-run AddModuleScore using corfiltered genes')
                seurat_obj = AddModuleScore (seurat_obj, filtered_geneset_list)
                seurat_obj@meta.data = seurat_obj@meta.data[, !colnames (seurat_obj@meta.data) %in% paste0(names(geneset_list),'_r',cor_threshold)]
                colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %in% paste0('Cluster',seq_along(geneset_list))] = paste0(names(geneset_list),'_r',cor_threshold)
                if (length (geneset_list) == 1) 
                	{
        					seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = ifelse (seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)] > pos_threshold, 'pos','neg')
        					pdf (paste0(outdir, listName, '_modulescore_distribution_cor_threshold_',cor_threshold,'_score_',pos_threshold,'.pdf'))
        					hist (seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)])
          				abline (v = pos_threshold)
          				dev.off()
        					} else {
                	seurat_obj@meta.data[, paste0(listName, '_r',cor_threshold,'_max')] = sapply (seq_along(colnames(seurat_obj)), function(x) colnames(seurat_obj@meta.data[,paste0(names(geneset_list),'_r',cor_threshold)])[which.max (seurat_obj@meta.data[x,paste0(names(geneset_list),'_r',cor_threshold)])])        
                	}
                }
        return (seurat_obj)
        } 

## Summarizes data.
## Obtained from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


# Pass GO terms or gene patters (can be regular expressions) and return a list of genes fitting the description. 
# For go terms specifiy the go term with spaces between words and no 'GO:' or 'GO_' at the beginning.
# e.g. NI.genes (ni.goterms = 'ribosome biogenesis', ni.genes = '^Mt')
NI.genes = function (all.genes = NULL, ni.goterms = NULL, ni.genes = NULL, org = 'mouse')
	{
	require (GO.db)
	if (org == 'mouse') 
		{
		require (org.Mm.eg.db)
		if (!is.null (ni.goterms))
			{	
			goterms = toTable (GOTERM)
			goid = goterms[goterms$Term %in% ni.goterms, 'go_id']
			ni.goterms.genes = unique(toTable (revmap(org.Mm.egGO)[goid])[,'gene_id'])
			ni.goterms = toTable (org.Mm.egSYMBOL[ni.goterms.genes])$symbol
			}
		if (!is.null(ni.genes))
			{
			if (is.null (all.genes)) all.genes = toTable (org.Mm.egSYMBOL)$symbol
			ni.genes = unique(unlist(lapply (ni.genes, function (x) all.genes[grep (x, all.genes)])))
			}
		return (unique(c(ni.goterms, ni.genes)))	
		}
	if (org == 'human') 
		{
		require (org.Hs.eg.db)
		if (!is.null (ni.goterms))
			{	
			goterms = toTable (GOTERM)
			goid = goterms[goterms$Term %in% ni.goterms, 'go_id']
			ni.goterms.genes = unique(toTable (revmap(org.Hs.egGO)[goid])[,'gene_id'])
			ni.goterms = toTable (org.Hs.egSYMBOL[ni.goterms.genes])$symbol
			}
		if (!is.null(ni.genes))
			{
			if (is.null (all.genes)) all.genes = toTable (org.Hs.egSYMBOL)$symbol	
			all.genes = toTable (org.Hs.egSYMBOL)$symbol
			ni.genes = unique(unlist(lapply (ni.genes, function (x) all.genes[grep (x, all.genes)])))
			}
		return (unique(c(ni.goterms, ni.genes)))	
		}
	}


# Compute overlaps of all combinations of vectors within a list
ovmat = function (ovlist, df=FALSE, ov_threshold=0.5, compare_lists= NULL, palette=NULL)
  {
  if (!df) require (ComplexHeatmap)	
  comp_mat = sapply (ovlist, function(x) sapply(ovlist, function (y) sum (unique(x) %in% unique(y)) / min(c(length(unique(x)),length(unique(y))))))
  if (!is.null (compare_lists)) comp_mat = comp_mat[compare_lists[[1]],compare_lists[[2]]]
  if (df) 
  	{
  	return (comp_mat)
  	} else {
    if (is.null(palette)) palette = viridis::mako(100)
  	return (Heatmap (
  		comp_mat, 
  		col=palette,
  		cell_fun = function(j, i, x, y, width, height, fill) {
      if (comp_mat[i,j] > ov_threshold) grid.text (sprintf("%.1f", comp_mat[i, j]), x, y, gp = gpar (fontsize=5))      
      }))      
  	}
  }


# New RankModuleScore function using AddModuleScore and list of markers to assign cluster to cell types 
clAssign = function (
	seuratObj = NULL, # seurat object
	markers = NULL,  # data.frame with a 'gene' and a 'celltype' column or named list
	metaGroupName = NULL, # metaGroup column specifying the clustering 
	returnDF = FALSE) # return a data.frame instead of seurat object
  {
  require (Seurat)  
  if (is.data.frame(markers)) markers = split (markers$gene, markers$celltype)  
  seuratObj@meta.data = seuratObj@meta.data[,!grepl ('clAssign', colnames(seuratObj@meta.data))]
  seuratObj@meta.data = seuratObj@meta.data[,!grepl ('clmodule', colnames(seuratObj@meta.data))]
  seuratObj = AddModuleScore (seuratObj, markers, name = 'clAssign')
  colnames(seuratObj@meta.data)[colnames(seuratObj@meta.data) %in% paste0('clAssign', seq_along(markers))] = paste0(names (markers), '_clmodule')
  agr = seuratObj@meta.data[,paste0(names (markers), '_clmodule')]
  agr = aggregate (agr, by = list(seuratObj@meta.data[,metaGroupName]), mean)
  rownames (agr) = agr[,1]
  agr = agr[,-1]
  colnames (agr) = names (markers)
  classigned = data.frame (metaGroup = rownames(agr), assigned = sapply (seq(nrow(agr)), function(x) colnames(agr)[which.max(agr[x,])]))
  if (returnDF) 
    {
    cellassigned = as.data.frame (table (seuratObj@meta.data[,metaGroupName]))
    return (cbind(classigned, data.frame(ncells = cellassigned$Freq[match(classigned$metaGroup, cellassigned$Var1)])))
    } else {
    seuratObj$clAssign = setNames (classigned$assigned, classigned$metaGroup)[seuratObj@meta.data[,metaGroupName]]
    return (seuratObj)
    }
  }

# From Levi Waldron --> https://github.com/lwaldron/LeviRmisc
writeGMT = function #Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### Function by Levi Waldron.
(object,
### R list object that will be converted to GMT file.  Each element
### should contain a vector of gene names, and the names of the
### elements will used for the gene set names
 fname
### Output file name for .gmt file
 ){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)
  for (iElement in 1:length(object)){
    write.table(t(c(make.names(rep(names(object)[iElement],2)),object[[iElement]])),
                sep="\t",quote=FALSE,
                file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
  }
### Called for the effect of writing a .gmt file
}

# Dotplot showing expression of a gene across two meta groups. Need a Seurat obj
geneDot = function (
  seurat_obj = srt,
	#mat_norm = srt@assays$RNA@data,
	#mat_counts = srt@assays$RNA@counts,
	gene = NULL, 
	x = NULL, # Vector of metagroup1 of equal length to mat columns, If multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
	y = NULL, # vector of metagroup2 of equal length to mat columns
	x_name = 'genes',
  assay='RNA',
	y_name = 'clusters',
	min_expression = 0,
	facet_ncol = 5,
	lim_expression = NULL,
	scale.data = TRUE, # scale data when multiple genes are given
	plotcol = viridis::viridis(100),
	include_NA = TRUE,
  swap_axes = FALSE,
  returnDF = FALSE
	#... # arguments to pass to facet wrap
	)
	{
  require ("scCustomize")
  require ('tidyr')
  if (exists('levels_x')) rm ('levels_x')
  if (exists('levels_y')) rm ('levels_y')
	if (all (grepl ('^\\d', seurat_obj@meta.data[,y]))) seurat_obj@meta.data[,y] = paste0('C',seurat_obj@meta.data[,y])
  if (is.factor (seurat_obj@meta.data[,x])) levels_x = levels (seurat_obj@meta.data[,x]) else
  levels_x = unique (seurat_obj@meta.data[,x])
  if (is.factor (seurat_obj@meta.data[,y])) levels_y = levels (seurat_obj@meta.data[,y]) else 
  levels_y = unique (seurat_obj@meta.data[,y])
  
  seurat_obj@meta.data[,x] = gsub ('[-_ \\+]', '', seurat_obj@meta.data[,x])
  seurat_obj@meta.data[,y] = gsub ('[-_ \\+]', '', seurat_obj@meta.data[,y])
  
  if (exists ('levels_x')) levels_x = gsub ('[-_ \\+]', '', levels_x)
  if (exists ('levels_y')) levels_y = gsub ('[-_ \\+]', '', levels_y)
    
  # Compute percentage expression per cell group
  percent = Percent_Expressing (seurat_object = seurat_obj, assay=assay, threshold = min_expression, features = gene, group_by = x, split_by = y)
  percent$gene = rownames(percent)
  percent = gather (percent, groups, expression, 1:(ncol(percent)-1))
  percent$key = paste0(percent$gene, '_', percent$groups)
  colnames (percent)[colnames(percent) == 'expression'] = 'percent'
  
  # Compute average expression per cell group
  seurat_average = log2 (as.data.frame (AverageExpression (seurat_obj, assay=assay, features= gene, group.by = c(x,y))[[1]]) + 1)
  if (length(gene) == 1) rownames(seurat_average) = gene
  seurat_average$gene = rownames(seurat_average)
  seurat_average = gather (seurat_average, groups, expression, 1:(ncol(seurat_average)-1))
  seurat_average$key = paste0(seurat_average$gene, '_', seurat_average$groups)
  
  # Merge percentage and average expression data.frames
  dot.df = cbind (seurat_average, percent[match (seurat_average$key, percent$key),])
  dot.df = dot.df[,!duplicated (colnames(dot.df))]

  # scale data if scale is TRUE
  if (scale.data) dot.df = transform (dot.df, expression = ave (expression, gene, FUN = function(x) scale(x, scale=T, center=T)))

  dot.df$expression[dot.df$percent == 0] = NA
	dot.df$percent[dot.df$percent == 0] = NA

  dot.df$x_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[1]), levels = levels_x)
  dot.df$y_axis = factor (dot.df$gene, levels = gene) 
  if (!is.null (y)) dot.df$y_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
	if (!is.null (y) & length(gene) > 1) 
    {
    dot.df$y_axis = factor (dot.df$gene, levels = gene) 
    dot.df$z_axis = factor (sapply (dot.df$groups, function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
    }

  #dot.df$x_axis = factor (dot.df$x_axis, levels = dot.df$x_axis)
  if (swap_axes)  colnames (dot.df)[match (c('x_axis','y_axis'), colnames (dot.df))] = c('y_axis','x_axis')

	p = ggplot (data = dot.df, aes (x= x_axis, y= y_axis)) +
  	geom_point (shape=21, aes (fill= expression, size = percent), alpha=0.7,colour='black', stroke=0.3) +
    labs (x = x_name, y = y_name, title = ifelse(length(gene) > 1,'',gene), subtitle = paste('Min expression >', min_expression)) +
     scale_shape (solid = FALSE) +
  	theme_classic () +
  	theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1))
  if (!is.null (y) & length(gene) > 1) p = p + facet_wrap (~z_axis)  
  if (length (plotcol) > 3) p = p + scale_fill_gradientn (colors = plotcol) else
    p = p + scale_fill_gradient2 (low = plotcol[1], mid = plotcol[2], high=plotcol[3])	

	if (returnDF) return(dot.df) else 
  return (p)
	}

# Generate barplots or boxplots of meta groups proportions specified
# meta_groups: vector or meta group names:
# 1) first meta_group is the group on which proportions are calculated
# 2) second meta_group split first meta_group on x axes  
#	3) third meta_group will group barplots separately
# if splits include only one value runs barplot instead
cellComp = function (
	seurat_obj = NULL, 
	metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
	plot_as = 'box', # box or bar 
	pal = NULL,
	prop = TRUE,
	ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
	facet_ncol = 20,
	facet_scales = 'free',
	subset_prop = NULL, # subset prop table by any group in any column
	removeNA = TRUE,
	returnDF = FALSE
	)
	{
	require (ggplot2)	
	require (ggpubr)	
	if (is.data.frame (seurat_obj))
		{
		meta_groups_df = seurat_obj[,metaGroups]	
		} else {
		meta_groups_df = seurat_obj@meta.data[,metaGroups]
		}
	# Refactor to remove 0 groups
	#meta_groups_df =  as.data.frame(lapply(unclass(meta_groups_df),as.character),stringsAsFactors=T)
	if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
	#if(is.null(pal) & plot_as == 'box') pal = rainbow (length(unique(meta_groups_df[,3])))
 	if (prop)
 		{
 		ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
 		ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
 		} else {
 		ccomp_df = as.data.frame (table (meta_groups_df))	
 		}

 	if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
 	if (!is.null (subset_prop)) 
 		{
 		subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
 		ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
 		}
 	#colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')  
 	if (plot_as == 'box')
 		{
 		p = ggplot (ccomp_df, aes_string (x= metaGroups[2], y= 'Freq')) +
  			theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    		scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups > 2)) p = p + geom_boxplot(aes_string (fill= metaGroups[3]), outlier.size=.2, alpha = 0.7, lwd=.2) 
    else p = p + geom_boxplot(aes_string (fill= metaGroups[2]), outlier.size=.2, alpha = 0.7, lwd=.2)		
  	if (length(metaGroups) > 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[4])), scales=facet_scales, ncol=facet_ncol)
  	}	
  if (plot_as == 'bar')
 		{
 		p = ggplot (ccomp_df, aes_string (x= metaGroups[1], y= 'Freq')) +
  			geom_bar(position="stack", stat="identity", aes_string(fill= metaGroups[2])) +
    		#geom_bar(position="dodge", stat="identity", aes_string(fill= metaGroups[2])) +
    		theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    		scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
  	if (length(metaGroups) == 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[3])), scales=facet_scales, ncol=facet_ncol)
  	}
  if (returnDF) return(ccomp_df) else 
  return (p)
	}
	


# Retreive human / mouse orthologues using biomart
# Taken from (062321): https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/ 
# Specify the direction to get human orthologues from mouse genes (MtH) or mouse orthologues from human genes (HtM)

# Basic function to convert mouse to human gene names
mouseMan <- function (x, direction = 'MtH', use_db = NULL)
	{
	if (is.null(use_db))
		{
		require ("biomaRt")
		httr::set_config(httr::config(ssl_verifypeer = FALSE))

		human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
		mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
		#human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
		#mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")

		if (direction == 'HtM')
			{
			genes_converted_df = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
			}
		if (direction == 'MtH')
			{
			genes_converted_df = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
			}
		# Print the first 6 genes found to the screen
		print (head (genes_converted_df))
		return (genes_converted_df)
		} else {
		if (direction == 'HtM')
			{
			use_db = use_db[use_db$human_gene %in% x,]
			}
		if (direction == 'MtH')
			{
			use_db = use_db[use_db$mouse_gene %in% x,]
			}
		return (use_db)	
		}
	}		
  


# Build clustered dotplot with GSEA terms on y axes and clusters on x axes
dotGSEA = function (
	enrichmentsTest_list,
	type = c('fgsea','enrich'),
	top_pathways = NULL, 
	padj_threshold = 0.05,
	cluster_rows = T,
	cluster_cols = T,
	remove_ns_modules = TRUE)
	{
	require (tidyr)	
	if (type == 'fgsea')
		{
		sig_terms = unique(unlist(sapply (enrichmentsTest_list, function(x) x[x$padj < padj_threshold,'pathway'])))
		if (length(sig_terms) == 0) return (NULL)
		if (!is.null(top_pathways)) sig_terms = unique(unlist(sapply (enrichmentsTest_list, function(x) 
				{
				x =	na.omit(x)
				x = x[x$padj < padj_threshold,]
				x_top = x[order(-x$NES),]
				x_top = head(x_top[order(x_top$padj),'pathway'], top_pathways)
				x_bot = x[order(x$NES),]
				x_bot = head(x_bot[order(-x_bot$padj),'pathway'], top_pathways)
				c(x_top,x_bot)
				})))
		if (length(sig_terms) == 1)
			{
			mat_pvalue = t(as.matrix (sapply (enrichmentsTest_list, function(x) x$pval[match(sig_terms, x$pathway)])))
			mat_sizeLog = t(as.matrix (sapply (enrichmentsTest_list, function(x) x$NES[match(sig_terms, x$pathway)])))	
			} else {
			mat_pvalue = sapply (enrichmentsTest_list, function(x) x$pval[match(sig_terms, x$pathway)])
			mat_sizeLog = sapply (enrichmentsTest_list, function(x) x$NES[match(sig_terms, x$pathway)])	
			}
		}
	if (type == 'enrich')
		{
		sig_terms = na.omit(unique(unlist(sapply (enrichmentsTest_list, function(x) x[x$p.adjust < padj_threshold,'ID']))))
		if (length(sig_terms) == 0) return (NULL)
		if (!is.null(top_pathways)) sig_terms = unique(as.vector(unlist (sapply (enrichmentsTest_list, function(x) na.omit(head(x[x$p.adjust < padj_threshold,'ID'], top_pathways))))))

		if (length(sig_terms) == 1)
			{
			mat_pvalue = t(as.matrix(sapply (enrichmentsTest_list, function(x) x$p.adjust[match(sig_terms, x$ID)])))
			mat_sizeLog = t(as.matrix(sapply (enrichmentsTest_list, function(x) x$Count[match(sig_terms, x$ID)])))	
			} else {
			mat_pvalue = sapply (enrichmentsTest_list, function(x) x$p.adjust[match(sig_terms, x$ID)])
			mat_sizeLog = sapply (enrichmentsTest_list, function(x) x$Count[match(sig_terms, x$ID)])	
			}
		}
	
	mat_pvalue[is.na(mat_pvalue)] = 1
	mat_sizeLog[is.na(mat_sizeLog)] = 0
	
rownames (mat_pvalue) = sig_terms
colnames (mat_pvalue) = names (enrichmentsTest_list)
rownames (mat_sizeLog) = sig_terms
colnames (mat_sizeLog) = names (enrichmentsTest_list)
mat_sizeLog = as.data.frame (mat_sizeLog)
mat_pvalue = as.data.frame (mat_pvalue)

	if(cluster_rows & min (dim(mat_sizeLog)[1]) >= 2) 
		{
		d_row = dist (mat_sizeLog)
		d_row[is.na(d_row)] = max (d_row)	
		hc1_row = hclust (d_row, method = "ward.D")
		mat_sizeLog = mat_sizeLog[hc1_row$order,]
		mat_pvalue = mat_pvalue[hc1_row$order,]
		}
	if(cluster_cols & min (dim(mat_sizeLog)[2]) >= 2) 
		{
		d_col = dist (t(mat_sizeLog))	
		d_col[is.na(d_col)] = max (d_col)
		hc1_col = hclust (d_col, method = "ward.D")
		mat_sizeLog = mat_sizeLog[,hc1_col$order]
		mat_pvalue = mat_pvalue[,hc1_col$order]
		}
#mat_sizeLog = as.data.frame (mat_sizeLog)
mat_sizeLog[mat_pvalue > padj_threshold] = NA			
mat_pvalue[mat_pvalue > padj_threshold] = NA
mat_pvalue = -log10(as.data.frame (mat_pvalue))	
#colnames (mat_pvalue) = names (enrichmentsTest_list)

mat_pvalue$pathway = rownames (mat_pvalue)
if (remove_ns_modules) # remove modules with no enriched pathways
	{
	remove_cols	= apply (mat_pvalue, 2, function(x) !all (is.na(x)))
	mat_pvalue = mat_pvalue[, remove_cols, drop=F]
	mat_sizeLog = mat_sizeLog[, colnames(mat_pvalue)[colnames(mat_pvalue)!='pathway'], drop=F]
	}
if (type == 'fgsea')
	{
	gsea_df_pvalue = gather (mat_pvalue, cluster, nlogPvalue, colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue)-1], factor_key=TRUE)
	gsea_df_NES = gather (mat_sizeLog, cluster2, NES, colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)], factor_key=TRUE)
	gsea_df = cbind (gsea_df_NES, gsea_df_pvalue)
	#gsea_df$cluster = factor (gsea_df$cluster, levels = unique(gsea_df$cluster))
	gsea_df$pathway = factor (gsea_df$pathway, levels = unique(gsea_df$pathway))
	p = ggplot (data = gsea_df, aes (x=cluster, y=pathway,
	      fill= NES, size=nlogPvalue)) +
	      geom_point (shape=21, color='black') +
	      #scale_fill_gradient2 (low = 'blue',mid='white', high = 'red', midpoint=0) +
	      scale_fill_gradient2(
      	low = "#67A9CF", 
      	mid = "#F7F7F7", 
      	high = "#EF8A62", 
      	midpoint = 0
    		) +
	      labs (x = 'cluster', y = '-log10(pvalue)') +
	      theme_minimal() + 
	      theme(text = element_text(size=11)) +
	      theme(axis.text.y = element_text (angle = 0, vjust = 0.5, hjust=1)) +
	      theme(axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
	return (p)	
	}
if (type == 'enrich')
	{	
	enrich_df_pvalue = gather (mat_pvalue, cluster, nlogPvalue, colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue)-1], factor_key=TRUE)
	enrich_df_Count = gather (mat_sizeLog, cluster2, Count, colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)], factor_key=TRUE)
	enrich_df = cbind (enrich_df_Count, enrich_df_pvalue)
	#gsea_df$cluster = factor (gsea_df$cluster, levels = unique(gsea_df$cluster))
	enrich_df$pathway = factor (enrich_df$pathway, levels = unique(enrich_df$pathway))
	p = ggplot (data = enrich_df, aes (x=cluster, y=pathway,
	      fill= nlogPvalue, size=Count)) +
	      geom_point (shape=21, color='black') +
	      scale_fill_distiller (palette = "Spectral") +
	      labs (x = 'cluster', y = '-log10(pvalue)') +
	      theme_minimal() + 
	      theme(text = element_text(size=11)) +
	      theme(axis.text.y = element_text (angle = 0, vjust = 0.5, hjust=1)) +
	      theme(axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
	return (p)		
	}
	}


# Make heatmap of the expression difference of between two conditions (e.g. KO vs WT) aggregated by cluster and top DEG across the clusters
diffClustHeat = function (
	deg_genes_df = NULL, # it has to contains following columns: "p_val_adj", "avg_log2FC", gene", "cluster"
	topGenes = 20,
	pvalAdjTrheshold = 0.05,
	pvalAdjTrheshold2 = 1e-50,
	pvalAdjTrheshold3 = 1e-200,
	sort_by = 'p_val_adj',
	only_pos = F,
	# mat1 = NULL, # normalised counts mats (in seurat found in seuratobj@assays$RNA@data)
	# mat2 = NULL, # normalised counts mats (in seurat found in seuratobj@assays$RNA@data)
	# meta1 = NULL,
	# meta2 = NULL,
	addGene = NULL, # add a gene of interest in the heatmap
	plotcol = NULL,
	col_limit = NULL,
	return_mat = F,
	...
	)
	{
	require (ComplexHeatmap)
	require (circlize)	
	
	deg_genes_df1 = deg_genes_df
	sig_genes = unique (deg_genes_df1$gene[deg_genes_df1$p_val_adj < pvalAdjTrheshold])
	deg_genes_df = deg_genes_df1[deg_genes_df1$gene %in% sig_genes, ]
	#if (only_pos) deg_genes_df = deg_genes_df1[deg_genes_df1$gene %in% sig_genes & deg_genes_df1$avg_log2FC > 0, ]
	
	if (sort_by == 'p_val_adj')
		{
		geneNames = unique (unlist (lapply (split (deg_genes_df, deg_genes_df$cluster), function(x) 
		{
		x = x[order (x[,sort_by]),]	
		head(x, topGenes)$gene
		})))
		}
if (sort_by == 'avg_log2FC')
		{
		geneNames = unique (unlist (lapply (split (deg_genes_df, deg_genes_df$cluster), function(x) 
		{
		x = x[order (-x[,sort_by]),]
		if (only_pos) 
			{
			head(x, topGenes)$gene
			} else {
			c(head(x, topGenes)$gene, tail(x, 2)$gene)				
			}	
		})))
		}
	# Add gene specified in 'addGene' argument if needed
	if (!is.null(addGene)) geneNames = unique(append (geneNames, addGene))
	# Make mat of pvalues
	pval_mat = split (deg_genes_df1, deg_genes_df1$cluster)
	pval_mat = as.data.frame(sapply (pval_mat, function(x) x$p_val_adj[match (geneNames, x$gene)]))
	if (ncol(pval_mat) == 1) pval_mat = as.data.frame(t(pval_mat))
	rownames(pval_mat) = geneNames
	#pval_mat = pval_mat[!duplicated(pval_mat$gene),]
	#rownames(pval_mat) = pval_mat$gene
	#pval_mat = pval_mat[top_genes, ]
	pval_mat[is.na(pval_mat)] = 1
	#pval_mat[pval_mat < pvalAdjTrheshold] = 2
	#pval_mat = pval_mat[,-ncol(pval_mat), drop=F]

	# Aggregate mat1 by meta1
	# if (!any (rownames(mat1) %in% deg_genes_df2$gene)) 
	# 	{
	# 	stop ('No genes found matching mat')
	# 	} else {
	# 	# mat1 = as.data.frame (t(as.matrix(mat1[top_genes,, drop=F])))
		# mat1$cluster = meta1
		# mat1 = aggregate (.~ cluster, data = mat1, FUN = mean)
		# rownames(mat1) = mat1[,1]
		# mat1 = mat1[,-1, drop=F]
		# mat1 = mat1[colnames(pval_mat),]
		lfc_mat = split (deg_genes_df1, deg_genes_df1$cluster)
		lfc_mat = as.data.frame (sapply (lfc_mat, function(x) x$avg_log2FC[match (geneNames, x$gene)]))
		if (ncol(lfc_mat) == 1) lfc_mat = as.data.frame(t(lfc_mat))
		rownames (lfc_mat) = geneNames
		#rownames (mat1) = top_genes
		lfc_mat[is.na(lfc_mat)] = 0
		#lfc_mat = lfc_mat[,colnames (pval_mat), drop=F]
		# Aggregate mat2 by meta2
		# mat2 = as.data.frame (t(as.matrix(mat2[top_genes,])))
		# mat2$cluster = meta2
		# mat2 = aggregate (.~ cluster, data = mat2, FUN = mean)
		# rownames (mat2) = mat2[,1]
		# mat2 = mat2[,-1]
		# mat2 = mat2[colnames(pval_mat),]
		
		# mat1 = mat1[rownames(mat1) %in% rownames(mat2),]
		# mat2 = mat2[rownames (mat1),]
		# matDiff = mat1 - mat2
		if (is.null(col_limit)) col_limit = max (abs(lfc_mat)) 
		if (is.null(plotcol)) 
			{	
			plotcol = colorRamp2(c(-max (abs(lfc_mat)), 0, max (abs(lfc_mat))), c("green", "white", "red"))
			} else {
			plotcol = colorRamp2(c(-col_limit,0, col_limit), plotcol)
			}
		#pal_corr2 = colorRamp2(c(-.5, 0, .5), c("blue", "grey33", "red"))
    labels_annotation = HeatmapAnnotation(
        text = anno_text(rownames(lfc_mat), rot = 45, location = unit(1, "npc"), just = "right",gp =gpar(fontsize = 4)),
        annotation_height = max_text_width(rownames(lfc_mat))
    )
		ht = Heatmap (t(lfc_mat), 
		clustering_distance_columns = 'euclidean',
		clustering_distance_rows = 'euclidean',
		column_title = paste('* = < ',pvalAdjTrheshold,'top',topGenes,'genes'),
		col = plotcol,
		#heatmap_legend_param = list(
    #     at = c(col_limit, 0, -col_limit),#,
    #    direction = "horizontal",#,
    #    #legend_height = unit(5, "cm"),
    #    just = c('top')
    #    #title_position = "leftcenter-rot"
   # 		),
		bottom_annotation = labels_annotation,
		column_names_gp = gpar(fontsize = 0),
		row_names_gp = gpar(fontsize = 5),
		cell_fun = function (j, i, x, y, width, height, fill) 
						{
		       if (t(pval_mat)[i, j] < pvalAdjTrheshold3)
		       		{
		           grid.text("***", x, y, just='center', vjust=.8,
		           	gp = gpar(fontsize = 5, col='black'))
		          } else {
		        	if(t(pval_mat)[i, j] < pvalAdjTrheshold2)
		        			{
    							grid.text("**", x, y, just='center', vjust=.8,
    							gp = gpar(fontsize = 5, col='black'))   
			          	} else {
    							if(t(pval_mat)[i, j] < pvalAdjTrheshold)
    								{
    								grid.text("*", x, y, just='center', vjust=.8,
    								gp = gpar(fontsize = 5, col='black'))   			
			          		}}}
			}, ...)
	if (!return_mat) return (ht) 
	else return (list (lfc = lfc_mat,pval = pval_mat))
	}


# Store palettes in palettes.rds in the projdir specified (default to projdir_init within scrna pipeline)
store_pal = function (palette = NULL, projdir = projdir_init)
  {
   if (!file.exists(paste0(projdir_init,'palettes.rds')))
      {
      palettes = list()
      palettes[[names(palette)]] = palette[[1]]
       saveRDS (palettes, paste0(projdir,'palettes.rds'))
      } else {
      palettes = readRDS (paste0(projdir_init,'palettes.rds'))
      palettes[[names(palette)]] = palette[[1]]
      saveRDS (palettes, paste0(projdir_init,'palettes.rds'))
      }  
  return (readRDS (paste0(projdir_init,'palettes.rds')))
  }



# Generate barplots of cell fractions positive for each gene combination    
combgenes = function (seurat_obj, genes, min_expression=0)
  {
  genes_data = as.data.frame (t(seurat_obj@assays$RNA@counts[rownames(seurat_obj) %in% genes, ]))
  genes_data[genes_data >= min_expression] = TRUE
  res = apply (genes_data,1, function(x) paste(colnames(genes_data)[which(x == 1)], collapse= '_'))
  return (res)
  }




# Correlation test based on permutation. Created by ChatGPT
  
# Observed Spearman's rank correlation coefficient
perm_test_chatgpt = function (x,y, n_permutations = 100)
  {
  observed_correlation <- cor.test(x, y, method = "spearman")$estimate
  
  # Function for permutation test
  perm_test <- function(x, y, n_permutations = n_permutations) {
    obs_corr <- cor.test(x, y, method = "spearman")$estimate
    n <- length(x)
    perm_corr <- numeric(n_permutations)
    
    for (i in 1:n_permutations) {
      # Shuffle one of the variables (e.g., y)
      shuffled_y <- sample(y)
      perm_corr[i] <- cor.test(x, shuffled_y, method = "spearman")$estimate
    }
    
    # Calculate the p-value
    p_value <- sum(abs(perm_corr) >= abs(obs_corr)) / n_permutations
    return(p_value)
  }
  
  # Perform the permutation test
  p_value_permutation <- perm_test(x, y, n_permutations = n_permutations)
  
  # Print the observed Spearman's rank correlation coefficient and p-value
  cat("Observed Spearman's rank correlation coefficient:", observed_correlation, "\n")
  cat("P-value (permutation test):", p_value_permutation, "\n")
  return (c(pvalue = p_value_permutation, corr = observed_correlation))
  }




### Replacement function for ggforest to make a plot when the argument strata() is added
# Taken from stackoverflow: https://stackoverflow.com/questions/63821274/plotting-a-cox-ph-model-using-ggforest-in-rstudio-when-a-factor-is-stratified

ggforest2 <- function (model, data = NULL, main = "Hazard ratio", 
                       cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7, 
                       refLabel = "reference", noDigits = 2) {
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(class(model) == "coxph")
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if(var %in% colnames(data)) {
      if (terms[i] %in% c("factor", "character")) {
        adf <- as.data.frame(table(data[, var]))
        cbind(var = var, adf, pos = 1:nrow(adf))
      }
      else if (terms[i] == "numeric") {
        data.frame(var = var, Var1 = "", Freq = nrow(data), 
                   pos = 1)
      }
      else {
        vars = grep(paste0("^", var, "*."), coef$term, 
                    value = TRUE)
        data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                   pos = seq_along(vars))
      }
    } else {
      message(var, "is not found in data columns, and will be skipped.")
    }    
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", 
                            "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", "N", "p.value", 
                                                "estimate", "conf.low", 
                                                "conf.high", "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 4], noDigits + 1), 
                               toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, noDigits + 1), " ", 
                                 ifelse(toShowExpClean$p.value < 0.05, "*", ""), 
                                 ifelse(toShowExpClean$p.value < 0.01, "*", ""), 
                                 ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * 
    as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))
  p <- ggplot(toShowExpClean, 
              aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, 
                  xmax = seq_along(var) + 0.5, 
                  ymin = exp(rangeplot[1]), 
                  ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) +
    scale_fill_manual(values = c("#FFFFFF33",  "#00000033"), guide = "none") + 
    geom_point(pch = 15, size = 4) + 
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), 
                  width = 0.15) + 
    geom_hline(yintercept = 1, linetype = 3) + 
    coord_flip(ylim = exp(rangeplot)) + 
    ggtitle(main) + 
    scale_y_log10(name = "", labels = sprintf("%g", breaks), 
                  expand = c(0.02, 0.02), breaks = breaks) + 
    theme_light() +
    theme(panel.grid.minor.y = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(), 
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          plot.title = element_text(hjust = 0.5)) + 
    xlab("") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
             label = toShowExpClean$var, fontface = "bold", 
             hjust = 0, size = annot_size_mm) + 
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0, 
             label = toShowExpClean$level, 
             vjust = -0.1, size = annot_size_mm) + 
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), 
             label = toShowExpClean$N, fontface = "italic", hjust = 0, 
             vjust = ifelse(toShowExpClean$level == "", 0.5, 1.1),
             size = annot_size_mm) + 
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring), 
             label = toShowExpClean$estimate.1, size = annot_size_mm, 
             vjust = ifelse(toShowExpClean$estimate.1 == "reference", 0.5, -0.1)) + 
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring), 
             label = toShowExpClean$ci, size = annot_size_mm, 
             vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") +
    annotate(geom = "text", x = 0.5, y = exp(y_variable), 
             label = paste0("# Events: ", gmodel$nevent, 
                            "; Global p-value (Log-Rank): ", 
                            format.pval(gmodel$p.value.log, eps = ".001"), 
                            " \nAIC: ", round(gmodel$AIC, 2), 
                            "; Concordance Index: ", round(gmodel$concordance, 2)), 
             size = annot_size_mm, hjust = 0, vjust = 1.2, fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}


### For scRNA meso paper. 
### Compute correlation across samples between any group of cnmf and the sarcomatoid cnmf
cor_scs = function(
  seurat_obj = srt, 
  scs_sample_avg = scs_sample_avg, 
  cnmf = cnmf_spectra_unique)
  {
  ccomp_df = srt@meta.data[,names (cnmf_spectra_unique)]
  ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,'sampleID4',drop=F]), malig_statistics)
  rownames(ccomp_df) = ccomp_df[,1]
  ccomp_df = ccomp_df[,-1]
  colnames (scs_sample_avg) = 'scs_score'
  ccomp_df = cbind (ccomp_df, scs_sample_avg[rownames(ccomp_df),, drop=F])
  sp = lapply (names (cnmf_spectra_unique), function(x) ggscatter(ccomp_df, x = 'scs_score', y = x, 
            add = "reg.line", conf.int = TRUE, label = rownames(ccomp_df),
            cor.coef = TRUE, cor.method = "spearman",
            xlab = 'malig score', ylab = x)  
            )
  return (sp)
  }
