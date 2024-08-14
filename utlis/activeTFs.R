if (!exists ('mSE')) mSE = fetch_mat (archp, 'Motif')
    
seGroupMotif <- getGroupSE(ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corGSM_MM <- correlateMatrices(
    ArchRProj = archp,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGSM_MM = corGSM_MM[!grepl ('-AS',corGSM_MM$GeneScoreMatrix_name),]
corGSM_MM = corGSM_MM[!grepl ('-DT',corGSM_MM$GeneScoreMatrix_name),]
corGSM_MM = corGSM_MM[!grepl ('-OT',corGSM_MM$GeneScoreMatrix_name),]
corGSM_MM = corGSM_MM[!grepl ('-RAB5IF',corGSM_MM$GeneScoreMatrix_name),]
corGSM_MM = corGSM_MM[!grepl ('-IT2',corGSM_MM$GeneScoreMatrix_name),]
corGSM_MM = corGSM_MM[!grepl ('-C8orf76',corGSM_MM$GeneScoreMatrix_name),]
corGSM_MM = na.omit (corGSM_MM)
saveRDS (corGSM_MM, 'TF_activators_genescore.rds')
