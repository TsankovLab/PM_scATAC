# Set up conda environment # 
conda install r-essentials
conda install conda-forge::r-hdf5r
conda install r-devtools
conda install bioconda::bioconductor-dirichletmultinomial
conda install bioconda::bioconductor-chromvar

# Install SCENIC+ in the same conda environment
git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install -e .

# In R install:
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
install.packages ('remotes')
remotes::install_github("stuart-lab/signac", ref="develop")
install.packages ('Seurat')
devtools::install_github("crazyhottommy/scATACutils")
BiocManager::install ('ComplexHeatmap')
BiocManager::install ('regioneR')
BiocManager::instlall ('BSgenome.Hsapiens.UCSC.hg38')
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
BiocManager::install ('TxDb.Hsapiens.UCSC.hg38.knownGene')
install.packages(c("ggdendro"))
BiocManager::install(c("GenomicAlignments", "plyranges", "Rsamtools", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "Biostrings", "BiocGenerics", "S4Vectors", "GenomicFeatures"))
devtools::install_github("colomemaria/epiAneufinder")
devtools::install_github('immunogenomics/presto')
devtools::install_github("GreenleafLab/chromVARmotifs")
BiocManager::install("scran")
BiocManager::install ('scuttle')
BiocManager::install ('edgeR')
BiocManager::install ('BiocSingular')
BiocManager::install ('beachmat')

