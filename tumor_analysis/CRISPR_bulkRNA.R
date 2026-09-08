projdir = 'CRISPR_bulkRNA_Lillian'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)


# Load utils functions palettes and packages ####
source (file.path('..','git_repo','utils','load_packages.R'))
source (file.path('..','git_repo','utils','useful_functions.R'))
source (file.path('..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','git_repo','utils','palettes.R'))

suppressPackageStartupMessages(library(data.table))


### Read bulkRNA-seq data
mat = read.table ('PTCSYH-expression-matrix.tsv', sep='\t', header=T)
