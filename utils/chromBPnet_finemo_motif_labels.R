args <- commandArgs(trailingOnly = TRUE)
chromBPdir <- args[1]
celltype <- args[2]

library (httr)
library (XML)


# Load finemo calls of contribution counts ####
modisco_motifs = as.data.frame(readHTMLTable(file.path(chromBPdir,celltype, 'no_bias_model','modisco_counts',paste0(celltype,'_report'),'motifs.html')))
modisco_motifs$motif_match0 = sapply (modisco_motifs$NULL.match0, function(x) unlist(strsplit (x, '_'))[1])
modisco_motifs$motif_match1 = sapply (modisco_motifs$NULL.match1, function(x) unlist(strsplit (x, '_'))[1])
modisco_motifs$motif_match2 = sapply (modisco_motifs$NULL.match2, function(x) unlist(strsplit (x, '_'))[1])
modisco_motifs$direction = sapply (modisco_motifs$NULL.pattern, function(x) unlist(strsplit (x, '_'))[1])

finemo_hits = read.table(file.path(chromBPdir,celltype,'no_bias_model','finemo_out_counts','hits.tsv'), sep='\t', header=T)
finemo_hits = finemo_hits[!duplicated (paste(finemo_hits$chr,finemo_hits$start, finemo_hits$end)),]
finemo_hits$motif_name0 = modisco_motifs$motif_match0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$motif_name1 = modisco_motifs$motif_match1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$motif_name2 = modisco_motifs$motif_match2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$direction = modisco_motifs$direction[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$combined_motifs = paste(finemo_hits$motif_name0,finemo_hits$motif_name1,finemo_hits$motif_name2, sep='_')
finemo_hits$qvalue_m0 = modisco_motifs$NULL.qval0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$qvalue_m1 = modisco_motifs$NULL.qval1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$qvalue_m2 = modisco_motifs$NULL.qval2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]

write.table (finemo_hits[,c(1,2,3,18,17,6,19,20,21)], file.path(chromBPdir, celltype, 'no_bias_model', paste0(celltype,'_finemo_counts_to_genome_browser.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

# Load finemo calls of contribution profile ####
modisco_motifs = as.data.frame(readHTMLTable(file.path(chromBPdir,celltype, 'no_bias_model','modisco_profile',paste0(celltype,'_report'),'motifs.html')))
modisco_motifs$motif_match0 = sapply (modisco_motifs$NULL.match0, function(x) unlist(strsplit (x, '_'))[1])
modisco_motifs$motif_match1 = sapply (modisco_motifs$NULL.match1, function(x) unlist(strsplit (x, '_'))[1])
modisco_motifs$motif_match2 = sapply (modisco_motifs$NULL.match2, function(x) unlist(strsplit (x, '_'))[1])
modisco_motifs$direction = sapply (modisco_motifs$NULL.pattern, function(x) unlist(strsplit (x, '_'))[1])

finemo_hits = read.table(file.path(chromBPdir,celltype,'no_bias_model','finemo_out_profile','hits.tsv'), sep='\t', header=T)
finemo_hits = finemo_hits[!duplicated (paste(finemo_hits$chr,finemo_hits$start, finemo_hits$end)),]
finemo_hits$motif_name0 = modisco_motifs$motif_match0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$motif_name1 = modisco_motifs$motif_match1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$motif_name2 = modisco_motifs$motif_match2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$direction = modisco_motifs$direction[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$combined_motifs = paste(finemo_hits$motif_name0,finemo_hits$motif_name1,finemo_hits$motif_name2, sep='_')
finemo_hits$qvalue_m0 = modisco_motifs$NULL.qval0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$qvalue_m1 = modisco_motifs$NULL.qval1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
finemo_hits$qvalue_m2 = modisco_motifs$NULL.qval2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]

write.table (finemo_hits[,c(1,2,3,18,17,6,19,20,21)], file.path(chromBPdir, celltype,'no_bias_model', paste0(celltype,'_finemo_profile_to_genome_browser.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)




