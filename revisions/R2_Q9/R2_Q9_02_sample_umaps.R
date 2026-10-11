# R2_Q9 -- UMAPs coloured by sample, all compartments, with the improved palette.
# Original calls: scatac_main_ArchR.R (umap_p2), scatac_NKTcells_ArchR.R (umap_p3), scatac_myeloid_ArchR.R (umap_p1),
# scatac_tumor_SOX9.R (umap_p1), scatac_Bcells_ArchR.R (umap_p3), scatac_stroma_ArchR.R (umap_p3),
# scatac_endothelial_ArchR.R (umap_p3). Page 1 of each pdf = new palette, page 2 = old palette.
source ('R2_Q9_common.R')

comps = list (
  main        = list (dir = 'main/scatac_ArchR',              name = 'Sample',  embedding = 'UMAP',   labelMeans = FALSE),
  TNK         = list (dir = 'NKT_cells/scatac_ArchR',         name = 'Sample',  embedding = 'UMAP_H', labelMeans = TRUE),
  myeloid     = list (dir = 'myeloid_cells/scatac_ArchR',     name = 'Sample',  embedding = 'UMAP_H', labelMeans = FALSE),
  tumor       = list (dir = 'tumor_compartment/scatac_ArchR', name = 'Sample3', embedding = 'UMAP',   labelMeans = FALSE),
  Bcells      = list (dir = 'B_cells/scatac_ArchR',           name = 'Sample',  embedding = 'UMAP_H', labelMeans = FALSE),
  stroma      = list (dir = 'stroma/scatac_ArchR',            name = 'Sample',  embedding = 'UMAP_H', labelMeans = FALSE),
  endothelial = list (dir = 'Endothelial/scatac_ArchR',       name = 'Sample',  embedding = 'UMAP_H', labelMeans = FALSE))

for (cp in names (comps))
  {
  x = comps[[cp]]
  archp = load_project (x$dir)
  if (cp == 'tumor') archp = archp[!archp$Sample3 %in% c('normal1')]     # archp_NN of the original script
  p = lapply (list (palette_sample, palette_sample_old), function(pal) plotEmbedding (ArchRProj = archp,
    colorBy = 'cellColData', name = x$name, pal = pal, embedding = x$embedding, labelMeans = x$labelMeans))
  pdf (file.path (outdir, paste0 (cp, '_sample_umap.pdf')))
  print (p[[1]] + ggtitle (paste (cp, '- new palette')))
  print (p[[2]] + ggtitle (paste (cp, '- old palette')))
  dev.off()
  message (cp, ' done: ', nrow (archp@cellColData), ' cells')
  }
