conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of TUMOR compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/pbmc_normal/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','PM_scATAC','utils','load_packages.R'))
source (file.path('..','..','PM_scATAC','utils','useful_functions.R'))
source (file.path('..','..','PM_scATAC','utils','ggplot_aestetics.R'))
source (file.path('..','..','PM_scATAC','utils','scATAC_functions.R'))
source (file.path('..','..','PM_scATAC','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','PM_scATAC','utils','knnGen.R'))
source (file.path('..','PM_scATAC','utils','addCoax.R'))
source (file.path('..','PM_scATAC','utils','Hubs_finder.R'))
source (file.path('..','PM_scATAC','utils','hubs_track.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")


if (!file.exists ('Save-ArchR-Project.rds')) 
  { source ('../../PM_scATAC/scatac_pbmc_normal_create_ArchRobj.R')
  } else {
   archp = loadArchRProject (projdir)   
  }



  
# Export bigiwg files ####
metaGroupName = 'celltype'
exp_bigwig = T
if (exp_bigwig)
  {
  getGroupBW(
    ArchRProj = archp,
    groupBy = metaGroupName,
    normMethod = "ReadsInTSS",
    tileSize = 100,
    maxCells = 1000,
    ceiling = 4,
    verbose = TRUE,
    threads = getArchRThreads(),
    logFile = createLogFile("getGroupBW")
  )
  }

