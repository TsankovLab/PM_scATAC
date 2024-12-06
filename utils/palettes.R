library (paletteer)
library (circlize)

palette_celltype_simplified = c(
  B_cells = 'magenta2',
  T_cells = 'firebrick1',
  MonoMac = 'cornflowerblue',
  Myeloid = 'cornflowerblue',
  DCs = 'deepskyblue',
  NK = 'gold1',
  Fibroblasts = 'azure4',
  Endothelial = 'brown',
  SmoothMuscle = 'blueviolet',
  Plasma = 'lightsalmon1',
  Alveolar = 'black',
  #Malignant = 'palegreen4',
  Malignant = 'plum4',
  Mast = 'royalblue4',
  Mesothelium = 'olivedrab1',
  Glia = 'lawngreen',
  pDCs = 'tomato')



palette_gene_expression = as.character(paletteer::paletteer_c("grDevices::Blue-Red 3", 100))
palette_gene_expression2 = as.character(paletteer::paletteer_c("grDevices::Blue-Red 3", 100))[c(1,50,100)]
palette_gene_expression_fun = function(x) {return(colorRamp2(c(min(x), 0, max(x)), c(palette_gene_expression2)))}

palette_sample = c(as.character(paletteer::paletteer_d("impressionist.colors::la_chanson_du_chien")), 'darkblue', 'red','red','red')
palette_sample = setNames (palette_sample , c('P1','P13','P3','P12','P5','P11','P4','P8','P10','P14','P23', 'HU62','HU37','normal_pleura'))
palette_sample = c(palette_sample,P11_HOX = 'violet')
#palette_sample = setNames (as.character(paletteer::paletteer_c("pals::ocean.dense",13)))
#palette_sample = c(palette_sample, P10 = 'grey')
palette_sample2 = c(palette_sample[3], palette_sample[6], palette_sample[10])
palette_sampling = setNames(c('black','grey'), c('biopsy','resection'))
palette_SE_group = c('S-High'= unname(palette_sample[3]), 'E-High' = unname(palette_sample[10]))
palette_sample2_fun = colorRamp2(c(1, 0, -1), c(palette_sample2))
palette_sample_cnv = colorRamp2(c(.5, 0, -.5), c(palette_sample2))
palette_cnv = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging",3))
palette_cnv_fun = colorRamp2(c(-.5,0,.5), palette_cnv)

ov_pal = rev(paletteer::paletteer_c("grDevices::Purple-Blue",100))


# Set palette
palette_bulk = setNames (as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7,3,1,7)]), c('Epithelioid','Biphasic-E','Biphasic-S','Sarcomatoid','Biphasic','E_score','S_score'))

# palette_t_cells = setNames (as.character (paletteer::paletteer_d("fishualize::Bodianus_rufus",5)), c('CD8','CD4','Tregs','TFH','CD8_exhausted'))
# palette_nk_cells = setNames (c("#D3E3CAFF", "#92A587FF", "#2F3525FF"), c('NK_FGFBP2','NK_KLRC1','NKlike_Tcells'))
# palette_tnk_cells = c(palette_t_cells, palette_nk_cells)
palette_tnk_cells =setNames (as.character(paletteer::paletteer_d("tvthemes::Stark",n=6)), c('CD8','CD4','Tregs','NK_KLRC1','NK_FGFBP2','CD8_exhausted'))
palette_tnk_cells['NK_KLRC1'] = 'black'
#palette_tnk_cells['CD4'] = '#F9ECE8FF'
palette_myeloid = rev(as.character(paletteer::paletteer_d("khroma::lapaz", 256)[c(1,40,80,120,150,190,210,230)]))
palette_protein_expression = c(low="darkblue",mid= "white",high= "darkgreen") 
palette_feature_RNA = c('lightgrey',"#5F1415FF")
palette_feature_protein = c("lightgrey", "darkgreen")

palette_b_cells = c(B_cells = 'mediumorchid1', GC_B_cells = 'magenta4', Plasma = 'lightskyblue1')

palette_clonotype = setNames (c(as.character (paletteer::paletteer_d("beyonce::X58"))),c('NonExpanded','Small','Large'))
#palette_clonotype = palette_clonotype[1:4]

pallette_pbmc_celltype = setNames (rev(as.character(paletteer::paletteer_d("khroma::smoothrainbow")[c(1,3,5,7,9,11,13,15)])), c('CD4','CD8','ILC','MAIT','NK','NK Proliferating','NK_CD56bright','Treg'))

palette_stroma = as.character(paletteer::paletteer_d("colRoz::shark_bay", 6))
palette_stroma[c('Endothelial','LEC','Fibroblasts','SmoothMuscle')] = as.character(c('grey','lightblue4','olivedrab3','slateblue2'))
palette_endothelial = setNames(as.character(paletteer::paletteer_d("colRoz::shark_bay", 3)), c('Artery','PLVAP+EC','Vein'))
pallette_fetal_vs_adult = setNames (c('#DB3EB1FF', 'grey','purple'), c('fetal_lung','adult_lung','PM'))
palette_scenic = rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(50))

palette_celltypes_normal = setNames(brewer.pal(11,'Paired'),rev(c('AT1','AT2','Basal','Secretory','Ionocytes','Neuroendocrine','Tuft.like','Mesothelium','Fibroblast','SmoothMuscle','Ciliated')))


#for (pal in ls()[grep('palette',ls())]) store_pal (list(pal = get(pal)))
palette_expression = rev (as.character(paletteer::paletteer_c("grDevices::Purple-Blue",100)))

#palette_deviation = paletteer::paletteer_c("grDevices::PuRd",20)
palette_deviation = paletteer::paletteer_c("ggthemes::Red-Black-White Diverging",100)
palette_deviation_ggplot_fill = paletteer::scale_fill_paletteer_c("ggthemes::Red-Black-White Diverging", direction=-1)
palette_deviation_centered = colorRamp2(c(-4,-2,0,2,4), c(palette_deviation[1],palette_deviation[27],palette_deviation[44],palette_deviation[65],palette_deviation[100]))
palette_deviation_centered = colorRamp2(c(-1,0,0,1,2,3,4), c('white','white','grey','#F69322FF','#C73370FF','purple','black'))
palette_deviation_fun = function(x) {return (colorRamp2(c(-max(abs(x)), 0, max(abs(x))), c(palette_deviation[length(palette_deviation)],'white',palette_deviation[1])))}
palette_enrichment = rev (as.character(paletteer::paletteer_c("grDevices::Purples 3",100)))
palette_hubs_accessibility = paletteer_c("ggthemes::Classic Orange-White-Blue",20)
palette_pseudotime = colorRamp2(c(-4,0,4), c('white','white','black'))
palette_expression_cor_fun = colorRamp2(c(-1,0,1), c('white','white','black'))

palette_expression_correlation = paletteer::paletteer_c("ggthemes::Green-Blue-White Diverging",100)
#palette_expression_cor_fun = colorRamp2(c(-1,0,1), c('#grey','#FCFDFEFF','#2A6F3FFF'))
palette_cooccurrence = colorRamp2(c(0,1), c('white','cornsilk4'))
palette_expression_cor_fun = function(x) {return (colorRamp2(c(-max(abs(x)), 0,max(abs(x))), c('grey25','#FCFDFEFF','#2A6F3FFF')))}
#palette_expression_cor = c('#24693DFF','#F6F9FCFF','#4F7FAAFF')

palette_deviation_correlation = paletteer::paletteer_c("ggthemes::Red-Black-White Diverging",100)
palette_deviation_cor_fun = colorRamp2(c(-1,0,1), c('#49525EFF','#FFFCFCFF','#AE123AFF'))
palette_deviation2 = paletteer::paletteer_c("pals::ocean.curl",100)

palette_genescore = as.character(paletteer::paletteer_d("khroma::vik"))
palette_genescore_fun = function(x) {return (colorRamp2(c(-max(abs(x)), 0,max(abs(x))), c(palette_genescore[1],'white',palette_genescore[length(palette_genescore)])))}
#palette_deviation_cor_fun = colorRamp2(c(-1,0,1), c('#2B5C8AFF','white','#9E3D22FF'))
palette_fragments = paletteer::paletteer_c("grDevices::Oslo",n=40)
palette_fragments = paletteer::paletteer_c("grDevices::Inferno",40)
