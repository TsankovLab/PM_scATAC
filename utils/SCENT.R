library(SCENT)
data(net13Jun12)
print(dim(net13Jun12.m))
library (biomaRt)

norm_mat = srt@assays$RNA@data

# Convert to entrezID ####
message ('Convert symbol to entrezID')
mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_list=rownames (norm_mat)
entrez=getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = gene_list, bmHeader = T, mart = mart)
entrez = entrez[!is.na(entrez[,2]),]
entrez = entrez[!duplicated (entrez[,2]),]
entrez = entrez[!is.na(entrez[,1]),]
entrez = entrez[!duplicated(entrez[,1]),]
norm_mat = norm_mat[rownames (norm_mat) %in% entrez[,1],]
rownames (norm_mat) = entrez[,2][match (entrez[,1], rownames(norm_mat))]

# Run SCENT ####
#message ('SCENT integrate PPI')
#integ.l = DoIntegPPI(exp.m = norm_mat, ppiA.m = net13Jun12.m)
message ('SCENT run ComCCAT')
ccat.v = CompCCAT(exp = norm_mat, ppiA = net13Jun12.m)
saveRDS (ccat.v, paste0(projdir,'SCENT_results.rds'))
message ('SCENT results are in entropy_score column')
srt$entropy_score = ccat.v

