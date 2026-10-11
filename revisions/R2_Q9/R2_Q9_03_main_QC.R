# R2_Q9 -- QC plots of all cells per sample (scatac_main_ArchR.R, "QC plots": Plots/QC_plots.pdf) with the improved palette:
# fragment size distribution, TSS enrichment profile, violins of TSS enrichment and number of fragments.
source ('R2_Q9_common.R')
archp = load_project ('main/scatac_ArchR')

p1 = plotFragmentSizes(ArchRProj = archp, groupBy = 'Sample', pal = palette_sample)
p2 = plotTSSEnrichment(ArchRProj = archp, groupBy = 'Sample', pal = palette_sample)

p3 <- plotGroups(
  ArchRProj = archp,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  pal = palette_sample,
  alpha = 0.4,
  addBoxPlot = TRUE
 )

p4 <- plotGroups(
  ArchRProj = archp,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "nFrags",
  plotAs = "violin",
  pal = palette_sample,
  alpha = 0.4,
  addBoxPlot = TRUE
 )

pdf (file.path (outdir, 'main_QC_plots.pdf'))
print (wrap_plots (p1, p2 ,p3, p4))
dev.off()
