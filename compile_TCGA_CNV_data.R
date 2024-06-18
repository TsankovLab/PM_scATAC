
# Compute bins ####
ws = 1e6
ss = 5e5
blacklist = toGRanges (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists

if (!file.exists (paste0('windows_',ws,'_',ss,'.rds')))
  {
  windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
    windowSize = ws, slidingSize = ss)
  saveRDS (paste0('windows_',ws,'_',ss,'.rds'))
  } else {
  windows = readRDS (paste0('windows_',ws,'_',ss,'.rds'))  
  }

### Change SegMean to LOSS < -0.25; GAIN = > 0.25  ####
gain = .25
loss = -0.25
meso_CNV_gr_hg38_cnv = lapply (meso_CNV_gr_hg38, function(x) 
  {
   tmp = x
   tmp$Segment_Mean[tmp$Segment_Mean >= gain] = 1
   tmp$Segment_Mean[tmp$Segment_Mean <= loss] = -1
   tmp$Segment_Mean[tmp$Segment_Mean > loss & tmp$Segment_Mean < gain] = 0
   tmp
  })

### Bin CNV per patient and then take average across patients ####
force = F
if (!file.exists('TCGA_cnv_mat.rds') | force)
  {
  pb =progress::progress_bar$new(total = length(meso_CNV_gr_hg38_cnv))
  cnv_mat = matrix (nrow = length(windows), ncol = length(meso_CNV_gr_hg38_cnv))
  rownames (cnv_mat) = as.character(windows)
  for (pat in 1:length(meso_CNV_gr_hg38_cnv))
    {
    pb$tick()  
     ov = as.data.frame (findOverlaps (windows, meso_CNV_gr_hg38_cnv[[pat]]))
     ov_sum1 = split (ov, ov$queryHits)
     ov_sum2 = unlist(lapply (ov_sum1, function(x) sum (meso_CNV_gr_hg38_cnv[[pat]]$Segment_Mean[x$subjectHits])))
    cnv_mat[as.numeric(names(ov_sum1)),pat] = ov_sum2
    }
  saveRDS (cnv_mat, 'TCGA_cnv_mat.rds')
  } else {
  cnv_mat = readRDS ('TCGA_cnv_mat.rds')
  }

if (!file.exists('TCGA_cnv_sample_avg.rds'))
  {
  cnv_mat_avg = as.data.frame (rowSums (cnv_mat))
  colnames (cnv_mat_avg) = 'CNV_avg'
  cnv_mat_avg$CNV_avg_log = log2 (abs(cnv_mat_avg$CNV_avg) + 1)
  cnv_mat_avg$CNV_avg_log = cnv_mat_avg$CNV_avg_log * sign(cnv_mat_avg$CNV_avg)
  
  cnv_mat_avg = cnv_mat_avg[!duplicated(as.character(windows)),,drop=F]
  rownames (cnv_mat_avg) = as.character(windows)[!duplicated(as.character(windows))]
  cnv_mat_avg$seqnames =  as.character(seqnames (windows[!duplicated(as.character(windows))]))
  cnv_mat_avg$start =  start (windows[!duplicated(as.character(windows))])
  cnv_mat_avg$end =  end (windows[!duplicated(as.character(windows))])
  cnv_mat_avg = cnv_mat_avg[,c('seqnames','start','end','CNV_avg_log')]
  
  saveRDS ('TCGA_cnv_sample_avg.rds')
  
  #cnv_mat_avg$seqnames = sub ('chr','',cnv_mat_avg$seqnames)    
  } else {
  cnv_mat_avg = readRDS ('TCGA_cnv_sample_avg.rds')  
  }


