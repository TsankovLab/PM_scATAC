#### SCENIC ANALYSIS ####
#require (SCENIC)
#require (GENIE3)
#require (reticulate)

vg = length (genes.keep)
projdir_SC = file.path(projdir,'SCENIC',paste0('vg_',vg,'_mw_',motif_window))
if(!is.null (scenic_name)) projdir_SC = file.path(projdir_SC,scenic_name)
dir.create (file.path(projdir_SC, 'Plots'), recursive=T, showWarnings=F)

if (org == 'mouse') 
	{
	if (motif_window == 'tss500bp') motifs_tss = "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"	
	if (motif_window == 'tss10kbp') motifs_tss = "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
	motifs_weights = "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
	TFs = 'mm10_TF.txt'
	}
if (org == 'human') 
	{
	if (motif_window == 'tss500bp') motifs_tss = "hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather"		
	if (motif_window == 'tss10kbp') motifs_tss = "hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"
	motifs_weights = 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
	TFs = 'human_TF.txt'
	}

if (!file.exists (paste0('expr_vg_',vg,'.txt')) | force)
	{
	message ('generate expression matrix')
	if ('layers' %in% slotNames (srt@assays$RNA)) {
		count_mat = t(srt@assays$RNA@layers$count)
		colnames (count_mat) = rownames(srt)
		count_mat = count_mat[,genes.keep]
		rownames (count_mat) = colnames (srt)
		write.table (count_mat, file.path(projdir_SC, paste0('expr_vg_',vg,'.txt')), col.names=NA, sep='\t')
		} else {
		write.table (t(srt@assays$RNA@counts[genes.keep,]), file.path(projdir_SC, paste0('expr_vg_',vg,'.txt')), col.names=NA, sep='\t')
		}	
	}
expr_mat = paste0('expr_vg_',vg,'.txt')

# Make sh file executable
system (paste0('chmod +x ',scrna_pipeline_dir,'/SCENIC.sh'), wait=FALSE)

## Submit job on server
message ('submit SCENIC job')
#system (paste0('bsub ',file.path(scrna_pipeline_dir,'SCENIC.sh'),' ',paste0(projdir_SC,'/'),' ', motifs_tss,' ', motifs_weights,' ', TFs, ' ',expr_mat,' ',scrna_pipeline_dir,' ', vg,' ', motif_window), wait=FALSE)

command <- paste ("bsub -J", 'pySCENIC', 
    "-P acc_Tsankov_Normal_Lung -q premium -n 8 -W 6:00 -R rusage[mem=128000] -R span[hosts=1] -o",
    file.path(projdir_SC, 'SCENIC.out'), "-e" ,
    file.path(projdir_SC, 'SCENIC.err'),
    file.path(scrna_pipeline_dir,'SCENIC.sh'))
  args <- paste (paste0(projdir_SC,'/'), motifs_tss,motifs_weights,TFs,expr_mat,scrna_pipeline_dir,vg,motif_window) # try to see if it works removing projdir_init and projdir
  system (paste(command, args)) 