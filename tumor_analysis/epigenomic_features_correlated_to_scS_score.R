###############################################
#  Plot TF Correlation to Sarcomatoid Score
#  using TF Deviations (mMat) and GeneScores (gsMat)
###############################################

library(dplyr)
library(zoo)
library(ggplot2)

###---------------------------------------------------------
### 1. Load inputs and define samples
###---------------------------------------------------------

selected_TF <- readRDS("selected_TF.rds")
sarc_module <- "cNMF20"

archp_meta <- as.data.frame(archp@cellColData)

# Filter samples: remove normal lung, low-cell samples, and outliers
sams <- unique(archp_meta$Sample3)
sams <- sams[!sams %in% c("normal1","normal2","normal3","normal4","P3","P13","P11_HOX")]


###---------------------------------------------------------
### 2. Extract and scale motif deviations (TF activity)
###---------------------------------------------------------

# Fetch motif matrix only if not already loaded
mSE <- fetch_mat(archp, "Motif")
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name

# Keep only selected TFs and z-score normalize
mMat <- scale(mMat[selected_TF,])
mMat <- t(mMat)

# Split by sample
mMat <- lapply(sams, function(s) mMat[archp_meta$Sample3 == s,])
names(mMat) <- sams


###---------------------------------------------------------
### 3. Extract and scale GeneScore matrix
###---------------------------------------------------------

gsSE <- fetch_mat(archp, "GeneScore")
gsMat <- assays(gsSE)[[1]]
rownames(gsMat) <- rowData(gsSE)$name

# Keep only selected TFs and z-score normalize
gsMat <- scale(gsMat[selected_TF,])
gsMat <- t(gsMat)

# Split by sample
gsMat <- lapply(sams, function(s) gsMat[archp_meta$Sample3 == s,])
names(gsMat) <- sams


###---------------------------------------------------------
### 4. Extract and scale cNMF modules
###---------------------------------------------------------

cnmf_mat <- archp_meta[, grep("cNMF", colnames(archp_meta))]
cnmf_mat <- as.data.frame(t(scale(t(cnmf_mat))))

# Split by sample
cnmf_mat <- lapply(sams, function(s) cnmf_mat[archp_meta$Sample3 == s,])
names(cnmf_mat) <- sams


###---------------------------------------------------------
### 5. Order cells by sarcomatoid module + smooth (rolling mean)
###---------------------------------------------------------

bin_width <- 40   # Window size
overlap   <- 40   # Step size

# Helper function: ordered → rolling average
rollavg_df <- function(mat, order_by, bin_width, overlap) {
  mat_ordered <- mat[order(order_by),]
  as.data.frame(lapply(as.data.frame(mat_ordered), function(col)
    zoo::rollapply(
      col, 
      width = bin_width, 
      FUN = mean, 
      by = overlap, 
      partial = TRUE, 
      align = "left"
    )
  ))
}

# TF activity (motif)
mMat_ordered_avg <- lapply(sams, function(s)
  rollavg_df(mMat[[s]], cnmf_mat[[s]][,sarc_module], bin_width, overlap)
)
names (mMat_ordered_avg) = sams
# GeneScores
gsMat_ordered_avg <- lapply(sams, function(s) {
  df <- rollavg_df(gsMat[[s]], cnmf_mat[[s]][,sarc_module], bin_width, overlap)
  df[is.na(df)] <- 0  # Handle NaNs (e.g. HIC2 issues)
  df
})
names (gsMat_ordered_avg) = sams
# cNMF
cnmf_ordered_avg <- lapply(sams, function(s)
  rollavg_df(cnmf_mat[[s]], cnmf_mat[[s]][,sarc_module], bin_width, overlap)
)
names (cnmf_ordered_avg) = sams

###---------------------------------------------------------
### 6. Compute within-sample correlations
###---------------------------------------------------------

# Correlate each TF with sarcomatoid module score
compute_cor <- function(mat_list, cnmf_list, type_label) {
  cor_list <- lapply(sams, function(s) {
    cor_df <- as.data.frame(
      cor(mat_list[[s]], cnmf_list[[s]][,sarc_module], method = "spearman")
    )
    colnames(cor_df) <- "score"
    cor_df$sample <- s
    cor_df$TF <- rownames(cor_df)
    cor_df
  })
  df <- do.call(rbind, cor_list)
  df$type <- type_label
  df
}

m_cor_df  <- compute_cor(mMat_ordered_avg, cnmf_ordered_avg, "activity")
gs_cor_df <- compute_cor(gsMat_ordered_avg, cnmf_ordered_avg, "genescore")

# Rank TFs by median correlation
m_rank  <- m_cor_df %>% group_by(TF) %>% summarise(median = median(score)) %>% arrange(desc(median))
gs_rank <- gs_cor_df %>% group_by(TF) %>% summarise(median = median(score)) %>% arrange(desc(median))

# Combine ranks by max median per TF
combined_rank <- data.frame(activity = m_rank$median[match(gs_rank$TF, m_rank$TF)],
                            genescore = gs_rank$median)
max_rank <- apply(combined_rank, 1, max)

final_rank <- data.frame(TF = gs_rank$TF, max_rank = max_rank) %>%
  arrange(desc(max_rank))


###---------------------------------------------------------
### 7. Select top TFs and prepare for plotting
###---------------------------------------------------------

top_sarc_TF <- head(final_rank$TF, 20)
saveRDS(top_sarc_TF, "top_sarc_TF.rds")

combined_df <- rbind(m_cor_df, gs_cor_df)
combined_df <- combined_df[combined_df$TF %in% top_sarc_TF,]

combined_df$TF <- factor(combined_df$TF, levels = top_sarc_TF)

palette_expression_disc <- paletteer::paletteer_c(
  "grDevices::Purple-Blue",
  length(top_sarc_TF)
)


###---------------------------------------------------------
### 8. Plot boxplots
###---------------------------------------------------------

bp <- ggplot(combined_df, aes(x = TF, y = score, fill = TF, color = type)) +
  geom_boxplot(
    aes(group = interaction(TF, type)),
    position     = position_dodge(0.8),
    linewidth    = 0.4,
    width        = 0.7,
    outlier.alpha = 0.2,
    outlier.shape = NA,
    size         = 0.5,
    alpha        = 0.6
  ) +
  geom_point(
    aes(group = interaction(TF, type)),
    position = position_dodge(0.8),
    alpha    = 0.5,
    size     = 1
  ) +
  gtheme +
  scale_fill_manual(values = palette_expression_disc) +
  scale_color_manual(values = c(activity = "#AE123AFF", genescore = "#001260FF")) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

pdf("Plots/sarcomatoid_score_TF_activity_boxplots.pdf", width = 7, height = 3)
print (bp)
dev.off()





### Import H3K27Ac data from primary cells from ENCODE ####
encode_h3k27ac_path = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/ENCODE/H3K27Ac'
encode_h3k27ac_files = file.path(encode_h3k27ac_path, list.files(encode_h3k27ac_path,pattern = '.bed'))
if (!file.exists(file.path (encode_h3k27ac_path,'encode_h3k27ac.rds')))
  {
  encode_h3k27ac_gr = lapply (encode_h3k27ac_files, function(x) 
    {
    tmp = read.table (x)
    colnames(tmp) = c('chr','start','end')
    makeGRangesFromDataFrame (tmp)
    })
  names (encode_h3k27ac_gr) = sub('.bed','',list.files(encode_h3k27ac_path,pattern = '.bed'))
  metadata = read.table (file.path(encode_h3k27ac_path,'metadata.tsv'), sep='\t', header=T)
  metadata = metadata[metadata$File.assembly == 'GRCh38',]
  encode_h3k27ac_gr = encode_h3k27ac_gr[names (encode_h3k27ac_gr) %in% metadata$File.accession]
  names (encode_h3k27ac_gr) = paste0(metadata$Biosample.type,'|', metadata$Biosample.term.name)[match(names(encode_h3k27ac_gr), metadata$File.accession)]
  encode_h3k27ac_gr = encode_h3k27ac_gr[!duplicated (names(encode_h3k27ac_gr))]
  saveRDS (encode_h3k27ac_gr, file.path (encode_h3k27ac_path,'encode_h3k27ac.rds'))
  } else {
  encode_h3k27ac_gr = readRDS (file.path(encode_h3k27ac_path,'encode_h3k27ac.rds'))
  }

archp = addBgdPeaks (archp, force= FALSE)
archp = addPeakAnnotations (ArchRProj = archp, force=T,
     regions = encode_h3k27ac_gr, name = "ENCODE_H3K27Ac")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "ENCODE_H3K27Ac",
  force = FALSE
  )

# Compute correlation of sarcomatoid cNMF external celltypes ####
sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal1','P3','P13','P11_HOX')] # Remove normal and low number samples and outliers
cnmf_mat = archp@cellColData[,grep ('cNMF', colnames(archp@cellColData))]
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams

# # Get ENCODE matrix ####
if (!exists('enSE')) enSE = fetch_mat (archp, 'ENCODE_H3K27Ac')
all (colnames(enSE) == rownames(archp))
encode_Mat = assays (enSE)[[1]]
rownames (encode_Mat) = rownames(rowData (enSE))

ext_mat = encode_Mat
remove_immune = c(
  rownames (ext_mat)[grep ('CD4',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('lymph',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD8',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T-',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.T\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('killer',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.B\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('NK',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Treg',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.T',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mast',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mac',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Plasma',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Monocyte',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Myelo',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('DC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD5',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('myeloid',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('B cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('neutrophil',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('ILC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('mononuclear',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Eryth',rownames(ext_mat), ignore.case =F)])

ext_mat = ext_mat[!rownames(ext_mat)%in% remove_immune,]

ext_mat = lapply (sams, function(x) scale(ext_mat[,archp$Sample3 == x]))
names (ext_mat) = sams

cnmf_ext_l = list()
for (sam in sams)
  {
  cnmf_ext_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(ext_mat[[sam]]), method = 'spearman')
  }

cnmf_ext_array <- simplify2array(cnmf_ext_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(cnmf_ext_array, c(1, 2), median)
median_matrix = median_matrix[, sapply(colnames(median_matrix), function(x) nchar(x) < 100)]

# Also plot boxplots ordered by correlation to sarcomatoid score ####
sarc_module = 'cNMF20'
top_ext_ct = median_matrix[sarc_module,]
top_ext_ct = top_ext_ct[order (-top_ext_ct)]
top_ext_ct2 = sapply (names(top_ext_ct), function(x) unlist(strsplit (x, '\\|'))[1])
top_ext_ct = names (unlist(unname(lapply (split(top_ext_ct2, top_ext_ct2),function(x) c(head(x,5),tail(x,5))))))
#top_ext_ct = names(c(head (top_ext_ct, Inf)))#,tail (top_ext_ct, 20)))
sarc_tf = lapply (sams, function(sam) cor (t(ext_mat[[sam]][top_ext_ct,]), t(cnmf_mat[[sam]][sarc_module,,drop=F]), method = 'spearman'))
sarc_tf_df = do.call (cbind,sarc_tf)

sarc_tf_df = as.data.frame (sarc_tf_df)
#sarc_tf_df2$TF = rownames(sarc_tf_df2)
colnames(sarc_tf_df) = sams
sarc_tf_df$ENCODE_celltype = factor (top_ext_ct, levels = top_ext_ct)
sarc_tf_df = gather (sarc_tf_df, scS, score, 1:(ncol(sarc_tf_df)-1))
sarc_tf_df$biotype = sapply (as.character(sarc_tf_df$ENCODE_celltype), function(x) unlist(strsplit (x, '\\|'))[1])
#sarc_tf_df = do.call (rbind, lapply(split(sarc_tf_df, sarc_tf_df$biotype), function(x) head (x,10)))
bp = lapply (unique (sarc_tf_df$biotype), function(x) ggplot (sarc_tf_df[sarc_tf_df$biotype == x,], aes (x = score, y = ENCODE_celltype)) + 
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey22", size=.5) + 
  gtheme + 
  geom_boxplot (aes (fill = ENCODE_celltype),
    linewidth = .2,
    width=.6,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=1, alpha=0.4
     ) +
  scale_fill_manual (values = paletteer::paletteer_c("grDevices::Purple-Blue",length(unique(sarc_tf_df[sarc_tf_df$biotype == x,]$ENCODE_celltype))))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = .4) + NoLegend() + coord_flip())
  

pdf (file.path ('Plots','sarcomatoid_score_ENCODE_H3K27_boxplots.pdf'), width = 11,height=5)
print (wrap_plots (bp,ncol=4))
dev.off()



# Read in peak files from scATAC studies ####
normal_lung_peaks = readRDS (file.path('..','..','git_repo','files','normal_lung_peaks.rds'))
# ### chromVAR analysis
archp = addBgdPeaks (archp, force= FALSE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = normal_lung_peaks, name = "Normal_Lung")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "Normal_Lung",
  force = FALSE
)


# Compute correlation of sarcomatoid cNMF external celltypes ####
sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal1','P3','P13','P11_HOX')] # Remove normal and low number samples and outliers
cnmf_mat = archp@cellColData[,grep ('cNMF', colnames(archp@cellColData))]
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams

# # Get Normal Lung matrix ####
if (!exists('nlSE')) nlSE = fetch_mat (archp, 'Normal_Lung')
all (colnames(nlSE) == rownames(archp))
normal_Mat = assays (nlSE)[[1]]
rownames (normal_Mat) = rownames(rowData (nlSE))

# # Get external matrix ####
# if (!exists('fSE')) fSE = fetch_mat (archp, 'scATAC_datasets')
# all (colnames(fSE) == rownames(archp))
# fetal_Mat = assays (fSE)[[1]]
# rownames (fetal_Mat) = rownames(rowData (fSE))

ext_mat = normal_Mat
remove_immune = c(
  rownames (ext_mat)[grep ('CD4',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('lymph',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD8',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T-',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.T\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('killer',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.B\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('NK',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Treg',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.T',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mast',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mac',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Plasma',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Monocyte',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Myelo',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('DC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD5',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('myeloid',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('B cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('neutrophil',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('ILC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('mononuclear',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Eryth',rownames(ext_mat), ignore.case =F)])

ext_mat = ext_mat[!rownames(ext_mat)%in% remove_immune,]



ext_mat = lapply (sams, function(x) scale(ext_mat[,archp$Sample3 == x]))
names (ext_mat) = sams

cnmf_ext_l = list()
for (sam in sams)
  {
  cnmf_ext_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(ext_mat[[sam]]), method = 'spearman')
  }

cnmf_ext_array <- simplify2array(cnmf_ext_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(cnmf_ext_array, c(1, 2), median)
median_matrix = median_matrix[, sapply(colnames(median_matrix), function(x) nchar(x) < 100)]

# plot boxplots ordered by correlation to sarcomatoid score ####
sarc_module = 'cNMF20'
#sarc_module = 'cNMF_sarc_epit'
top_ext_ct = median_matrix[sarc_module,]
top_ext_ct = top_ext_ct[order (-top_ext_ct)]
#top_ext_ct2 = sapply (names(top_ext_ct), function(x) unlist(strsplit (x, '\\|'))[1])
#top_ext_ct = c(head(top_ext_ct,5),tail(top_ext_ct,5))
#top_ext_ct = names(c(head (top_ext_ct, Inf)))#,tail (top_ext_ct, 20)))
sarc_tf = lapply (sams, function(sam) cor (t(ext_mat[[sam]][names(top_ext_ct),]), t(cnmf_mat[[sam]][sarc_module,,drop=F]), method = 'spearman'))
sarc_tf_df = do.call (cbind,sarc_tf)

sarc_tf_df = as.data.frame (sarc_tf_df)
#sarc_tf_df2$TF = rownames(sarc_tf_df2)
colnames(sarc_tf_df) = sams
names (top_ext_ct) = sub ('Tsankov_lung_','',names(top_ext_ct))
sarc_tf_df$normal_celltype = factor (names(top_ext_ct), levels = rev(unique(rev(names(top_ext_ct)))))
sarc_tf_df = gather (sarc_tf_df, scS, score, 1:(ncol(sarc_tf_df)-1))
#sarc_tf_df$biotype = sapply (as.character(sarc_tf_df$ENCODE_celltype), function(x) unlist(strsplit (x, '\\|'))[1])
#sarc_tf_df = do.call (rbind, lapply(split(sarc_tf_df, sarc_tf_df$biotype), function(x) head (x,10)))
palette_expression_disc = paletteer::paletteer_c("grDevices::Purple-Blue",length(unique(sarc_tf_df$normal_celltype)))
head (sarc_tf_df)

remove_celltypes = c('proximal.Ciliated','proximal.Secretory','Neuronal')
sarc_tf_df = sarc_tf_df[!sarc_tf_df$normal_celltype %in% remove_celltypes,]

bp = ggplot (sarc_tf_df, aes (x = score, y = normal_celltype)) + 
  #geom_boxplot (alpha=.8, outlier.shape = NA, fill = 'darkred', color='grey22') + 
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey22", size=.5) + 
  geom_boxplot (aes (fill = normal_celltype),
    linewidth = .2,
    width=.6,
    outlier.alpha = 0.2,
    outlier.size = 0,
     size=1, alpha=0.6
     ) +
  scale_fill_manual(values = palette_expression_disc) +
  gtheme + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = .4)+
  coord_flip()
  

pdf (paste0 ('Plots/sarcomatoid_score_normal_lung_boxplots4.pdf'), width = 6,height=3)
bp
dev.off()








