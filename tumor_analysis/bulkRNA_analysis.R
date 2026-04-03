# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

set.seed(1234)

# Set directory
projdir = 'tumor_compartment'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd(projdir)

source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
bulk_palette = setNames (hue_pal()(4)[c(4,3,2,1)],c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
bulk_palette = setNames (as.character(paletteer::paletteer_d("rockthemes::heep")), c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid'))
bueno_palette = setNames (paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7)], c('Sarcomatoid','Biphasic-S','Biphasic-E','Epithelioid'))
bulk_palette = setNames (as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7,3)]), c('Epithelioid','Biphasic-E','Biphasic-S','Sarcomatoid','Biphasic'))


# Save data and metadata ####
meso_bulk_l = readRDS (file.path ('..','..','git_repo','files','bulk_RNA_studies.rds'))
meso_bulk_meta_l = readRDS (file.path ('..','..','git_repo','files','bulk_RNA_studies_metadata.rds'))

### Query bulk data ####
studies = c('bueno','tcga','mesomics')
module_l = list(SOX9 = 'SOX9',SOX6 = 'SOX6',RUNX2='RUNX2',RUNX1='RUNX1',SNAI2='SNAI2')

# Run genes on bulk datasets ####
stat_testL2 = list()
for (mod_name in names(module_l))
  {
  mod = module_l[which (mod_name == names(module_l))]
  bp_l = list()  
  stat_testL = list()  
  for (study in studies)
    {
    meso_bulk = meso_bulk_l[[study]]
    meso_bulk_meta = meso_bulk_meta_l[[study]]
    
    # gene_exp$histology = meso_bulk_meta[match (gene_exp$sample, names(meso_bulk_meta))]  
    # gene_exp = gene_exp[!is.na(gene_exp$histology),]
    gene_exp = data.frame(
      expression = colMeans (meso_bulk[rownames(meso_bulk) %in% mod[[1]],,drop=F]),
      subtype = meso_bulk_meta$subtype
      )
    bp_l[[study]] = ggplot(gene_exp, aes (x= subtype, y= expression)) + 
    ggtitle (paste('RNAseq',study,'cohort')) +
    gtheme +  
    geom_boxplot (aes (fill = subtype), outlier.colour="black", outlier.shape=16,
             outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
    NoLegend () + 
    scale_fill_manual (values = bulk_palette)
  
    stat.test = bp_l[[study]]$data %>%
          t_test(reformulate ('subtype', 'expression')) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
    stat.test = stat.test %>% add_xy_position (x = 'subtype', step.increase=0.5)
    stat.test = stat.test[stat.test$group1 %in% c('Sarcomatoid','Epithelioid') & stat.test$group2 %in% c('Sarcomatoid','Epithelioid'),]
    bp_l[[study]] = bp_l[[study]] + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = T,
      label = "p.adj.signif") + 
      NoLegend()
     stat_testL[[study]] = data.frame (stat.test, dataset = study)
    }
    stat_testL2[[mod_name]] = stat_testL
    pdf (paste0 ('Plots/',mod_name,'_expression_avg_boxplots.pdf'), width = 5,height=3)
    print (wrap_plots (bp_l, ncol=3))
    dev.off ()
    }

