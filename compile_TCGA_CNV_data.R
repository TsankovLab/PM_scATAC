

if(!file.exists ("'TCGA_CNV_hg38.rds'")) 
  {
  if (!file.exists (file.path ('TCGA_CNV2','gdac.broadinstitute.org_MESO.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0','MESO.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt')))
    {
    library (RTCGA)
    (cohorts <- infoTCGA() %>% 
       rownames() %>% 
       sub("-counts", "", x=.))
    checkTCGA('Dates')
    
    releaseDate <- "2016-01-28"
    cohorts = 'MESO'
    for (cohort in cohorts) {
      dir.create ('TCGA_CNV2')
      try(downloadTCGA( cancerTypes = cohort, destDir = "TCGA_CNV2", date = releaseDate, dataSet = "Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3" ),
          silent=TRUE)    
    } 
    } else {
    meso_CNV <- read.table(file.path ('TCGA_CNV2','gdac.broadinstitute.org_MESO.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0','MESO.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'),h=T) 
    }
  meso_CNV = split (meso_CNV, meso_CNV$Sample)
  meso_CNV = lapply (meso_CNV, function(x) {x$portion = sapply (x$Sample, function(x) unlist(strsplit (x, '-'))[4]); x})
  table (sapply (meso_CNV, function(x) unique(x$portion)))
  meso_CNV = meso_CNV[sapply (meso_CNV, function(x) unique(x$portion == '01A'))]
  meso_CNV = lapply (meso_CNV, function(x) {
    x$Chromosome = paste0('chr',x$Chromosome)
    x$Chromosome[x$Chromosome == 'chr23'] = 'chrX'
    x = x[,-1]
    x})
  meso_CNV_gr = lapply (meso_CNV, function(x) makeGRangesFromDataFrame (x,keep.extra.columns=T))
  
  library ("liftOver")
  if(!file.exists ("hg19ToHg38.over.chain")) 
    {
    download.file("https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
    system("gzip -d hg19ToHg38.over.chain.gz")
    }
  z <- import.chain ("hg19ToHg38.over.chain")
  #seqlevelsStyle (fragments_normal_flt) = "UCSC"  # necessary
  meso_CNV_gr_hg38 = lapply (meso_CNV_gr, function(x) unlist (liftOver(x, z)))
  saveRDS (meso_CNV_gr_hg38, 'TCGA_CNV_hg38.rds')
  } else {
  meso_CNV_gr_hg38 = readRDS ('TCGA_CNV_hg38.rds')  
  }

# Compute bins ####
ws = 1e6
ss = 5e5
blacklist = toGRanges (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists

if (!file.exists (paste0('windows_',ws,'_',ss,'.rds')))
  {
  windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
    windowSize = ws, slidingSize = ss)
  saveRDS (windows, paste0('windows_',ws,'_',ss,'.rds'))
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


