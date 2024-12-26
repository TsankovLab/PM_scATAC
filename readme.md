# Set up conda environment # 
conda install r-essentials  
conda install conda-forge::r-hdf5r  
conda install r-devtools  
conda install bioconda::bioconductor-dirichletmultinomial  
conda install bioconda::bioconductor-chromvar  
conda install bioconda::bioconductor-scran  
conda install bioconda::bioconductor-infercnv  
conda install bioconda::macs2  
conda install samtools  
conda install conda-forge::r-rstatix  
conda install conda-forge::cairo  


# Install SCENIC+ in the same conda environment
git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install -e .

# In R install:
install.packages("BiocManager")  
install.packages("paletteer")  
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages ("Signac")  
install.packages ('remotes'). 
remotes::install_github("stuart-lab/signac", ref="develop"). 
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
BiocManager::install("fgsea")  
BiocManager::install("clusterProfiler")  
BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", "ggrepel", "UCell"))
install.packages ('ape')
install.packages('dendextend')
conda install conda-forge::r-factoextra
BiocManager::install ('org.Hs.eg.db')
install.packages ('RColorBrewer')  
install.packages ('survminer')
#devtools::install_github("ricardo-bion/ggradar")
BiocManager::install("RTCGA")
BiocManager::install("liftOver")  
install.packages("scCustomize")  
install.packages("ggpubr")  
BiocManager::install("liftOver")
install.packages ('GSA')
devtools::install_github("GreenleafLab/chromVARmotifs")
BiocManager::install("compEpiTools")
install.packages("harmony")  
install.packages("zoo")
devtools::install_github('smorabit/hdWGCNA', ref='dev')  


### Run cNMF scripts on LSF batch submission ###   
Scripts to run cNMF are in utils and include:  
cnmf_prepare_inputs.R  
cnmf_master.sh  
cnmf_factorization_parallel.sh  
cnmf_format_spectra_files.R  

#### Run cNMF ####  
nfeat = 5000  # specify number of features to use  
force=F # force re-running when output already present  
k_list = c(5:30) # number of K to compute  
k_selections = c(5:30) # number of K to compute  
cores= 100 # specify number of cores for array job  
cnmf_name = 'CD8_exhausted' # specify name of cnmf run  
cnmf_out = path/to/output # path to cnmf output  
git_repo = path/to/git_repo # path to git repo  

### RUN consensus NMF from R ####  
source (file.path ('git_repo','utils','cnmf_prepare_inputs.R')) #  source cNMF to prepare count matrices and run cnmf algorithm  

### Import and format spectra files in R ####
k_selection = 30   # select a K after examining error / performance cNMF plot
source (file.path ('git_repo','utils','cnmf_format_spectra_files.R'))

### Install conda env for averaging h5 files from chromBPnet
conda create -n h5py python=3.8
conda install conda-forge/label/cf202003::hdf5plugin
conda install bioconda::wiggletools
conda install bioconda::ucsc-wigtobigwig



