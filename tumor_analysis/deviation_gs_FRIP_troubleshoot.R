### Check correlation with seq-depth for module scores ####
cnmf = archp@cellColData[,grep ('sarcomatoid',colnames(archp@cellColData))]
cnmf_scaled = as.data.frame (t(scale (t(cnmf))))

metaGroupNames = c('FRIP','nFrags','TSSEnrichment','nMonoFrags','nDiFrags','nMultiFrags','ReadsInTSS','ReadsInPromoter','PromoterRatio','ReadsInPeaks','NucleosomeRatio')

# hp = list()
# for (metagroupname in metaGroupNames)
# {
# cnmf_cor = data.frame (cor = cor (as.matrix(cnmf), as.matrix(archp@cellColData[,metagroupname])), tr = 'not scaled')
# cnmf_cor_scaled = data.frame (cor = cor (as.matrix(cnmf_scaled), as.matrix(archp@cellColData[,metagroupname])), tr = 'scaled')
# cnmf_df = rbind (cnmf_cor, cnmf_cor_scaled)
# hp[[metagroupname]] = ggplot (cnmf_df, aes (y = cor, x = tr)) +
#   #geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
#   geom_boxplot(alpha=.2)+
#   stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red") +
#   geom_jitter(width = 0.2, alpha = 1, color = "blue") +
#   geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) + 
#   labs(title = paste0('cNMF correlation to ',metagroupname,' Distributions'), x = "Value", y = "Frequency") +
#   scale_fill_manual(values = c("blue", "red")) +
#   theme_minimal()
# } 

# pdf (file.path ('Plots','histogram_cnmf_cor_FRIP.pdf'), width=10, height=5)
# wrap_plots(hp)
# dev.off()





hps_gs = list()
for (metagroupname in metaGroupNames)
{
cnmf_cor_sample_scaled = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(t(scale(t(cnmf[archp$Sample2 == x,])))), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])),
  tr = 'scaled', sample = x))

cnmf_cor_sample_scaled_df = do.call (rbind, cnmf_cor_sample_scaled)
cnmf_sample_cor = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(cnmf[archp$Sample2 == x,]), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])), 
tr = 'not scaled', sample = x))
cnmf_sample_cor_df = do.call (rbind, cnmf_sample_cor)
cnmf_df = rbind (cnmf_cor_sample_scaled_df, cnmf_sample_cor_df)
hps_gs[[metagroupname]] = ggplot (cnmf_df, aes (y = cor, x = tr, fill = sample)) +
  #geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  geom_boxplot(alpha=.2,outlier.shape = NA)+
  #stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red") +
  #geom_jitter(position=position_jitterdodge(), width = 0.2, alpha = 0.4, color = "blue", size=1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) + 
  labs(title = paste0('cNMF correlation to ',metagroupname,' Distributions'), x = "Value", y = "Frequency") +
#  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal()
} 

pdf (file.path ('Plots','histogram_cnmf_cor_FRIP_per_sample.pdf'), width=25, height=15)
wrap_plots(hps_gs, ncol = 5)
dev.off()

### Try with Deviations ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
all (colnames(mSE) == rownames(archp))
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name


hps = list()

for (metagroupname in metaGroupNames)
{
cnmf_cor_sample_scaled = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(t(scale(t(t(mMat)[archp$Sample2 == x,])))), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])),
  tr = 'scaled', sample = x))

cnmf_cor_sample_scaled_df = do.call (rbind, cnmf_cor_sample_scaled)
cnmf_sample_cor = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(t(mMat)[archp$Sample2 == x,]), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])), 
tr = 'not scaled', sample = x))
cnmf_sample_cor_df = do.call (rbind, cnmf_sample_cor)
cnmf_df = rbind (cnmf_cor_sample_scaled_df, cnmf_sample_cor_df)
hps[[metagroupname]] = ggplot (cnmf_df, aes (y = cor, x = tr, fill = sample)) +
  #geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  geom_boxplot(alpha=.2,outlier.shape = NA)+
  #stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red") +
  #geom_jitter(position=position_jitterdodge(), width = 0.2, alpha = 0.4, color = "blue", size=1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) + 
  labs(title = paste0('Deviation correlation to ',metagroupname,' Distributions'), x = "Value", y = "Frequency") +
#  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal()
} 

pdf (file.path ('Plots','histogram_deviations_cor_FRIP_per_sample.pdf'), width=25, height=15)
wrap_plots(hps, ncol=5)
dev.off()


head (hps$ReadsInTSS$data[hps$ReadsInTSS$data$tr == 'not scaled',][order (-hps$ReadsInTSS$data$cor[hps$ReadsInTSS$data$tr == 'not scaled']),],50)
hps$FRIP$data[grep ('JUN',rownames(hps$FRIP$data)),]


