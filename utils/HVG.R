library (scran)
set.seed(1234)

DefaultAssay (srt) = 'RNA'
if ('layers' %in% slotNames (srt@assays$RNA)) 
    {
    sce = SingleCellExperiment (list(counts=srt@assays$RNA@layers$counts, logcounts = srt@assays$RNA@layers$data),
    rowData=rownames(srt))
    rownames (sce) = rownames(srt)
    } else {
    sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
    rowData=rownames(srt))
    } 
sce = modelGeneVar(sce)
# remove batchy genes
if (org == 'human') batchy_genes = c('RPL','RPS','MT-') else
batchy_genes = c('Rpl','Rps','Mt-')

sce = sce[!apply(sapply(batchy_genes, function(x) grepl (paste0("^",x), rownames(sce))),1,any),]
vf = getTopHVGs (sce, n=nfeat)
    