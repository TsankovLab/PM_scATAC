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
bp
dev.off()

