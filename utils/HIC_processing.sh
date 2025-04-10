  # - conda-forge::python=3.8.10=h49503c6_1_cpython
  # - conda-forge::scipy=1.7.0=py38h7b17777_1
  # - conda-forge::numpy=1.21.1=py38h9894fe3_0
  # - bioconda::iced=0.5.10=py38h803c66d_0
  # - bioconda::bx-python=0.8.11=py38h024e602_1
  # - bioconda::pysam=0.16.0.1=py38hf7546f9_3
  # - bioconda::cooler=0.8.11=pyh3252c3a_0

  # - conda-forge::r-base=4.0.3=h349a78a_8
  # - conda-forge::r-ggplot2=3.3.5=r40hc72bb7e_0
  # - conda-forge::r-rcolorbrewer=1.1_2=r40h785f33e_1003
  # - conda-forge::r-gridbase=0.4_7=r40hc72bb7e_1003
  
  # - conda-forge::tbb=2020.2=hc9558a2_0
  # - bioconda::bowtie2=2.4.4=py38h72fc82f_0
  # - bioconda::samtools=1.12=h9aed4be_1
  # - bioconda::multiqc=1.11=pyhdfd78af_0

conda create -n hic python=3.8.10 scipy numpy 
conda install -c conda-forge gxx
pip install -U iced
conda install -c conda-forge -c bioconda bx-python
conda install bioconda::pysam=0.16.0.1 
conda install anaconda::pandas
conda install conda-forge::r-base
conda install conda-forge::r-ggplot2

conda install bioconda::bowtie2
conda install bioconda::samtools
conda install multiqc=1.11
#conda install pysam


projdir
cd HiCPro/HiC-Pro
conda activate meso_scatac
# in R

   ###
   ## How to generate chromosome size files ?
   ###
R   
require (BSgenome.Hsapiens.UCSC.hg38)
human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:22]
chrom.size <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[human_chr]
write.table(chrom.size, file="chrom_hg38.sizes", quote=FALSE, col.names=FALSE, sep="\t")

# Generate restriction fragments file
bin/utils/digest_genome.py -r mboi -o hg38_mboI.bed ../../../genome_assemblies/hg38.fa


conda activate hic
# hicpro is actually already installed on minerva (but it doesnt work). Need to clone HiC-Pro repo then make configure then cp manually to ~/local/bin/HiC-Pro_3.1.0 and then edit config file to change paths to ~/local/bin/HiC-Pro_3.1.0 and then do make install
#ml hicpro

projdir
cd HiCPro/HiC-Pro
#~/local/bin/HiC-Pro_3.1.0/bin/HiC-Pro -i /sc/arion/scratch/leew17/HiC/rawdata -o /sc/arion/scratch/giottb01/HiCPro -c /sc/arion/projects/Tsankov_Normal_Lung/Bruno/HiCPro_files/config_file.txt





# TRy send it as a job
chmod +x /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/HiC_Pro.sh
bsub -J HiC_step1 -P acc_Tsankov_Normal_Lung -q premium -n 8 -W 144:00 -R rusage[mem=16000] -R span[hosts=1] -o HiC.out -e HiC.err /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/HiC_Pro.sh


### RE-run with higher resolution contact maps
~/local/bin/HiC-Pro_3.1.0/bin/HiC-Pro -i /sc/arion/scratch/giottb01/HiCPro2/hic_results/data -o /sc/arion/scratch/giottb01/HiCPro_matrices -c /sc/arion/projects/Tsankov_Normal_Lung/Bruno/HiCPro_files/config_file.txt -s build_contact_maps -s ice_norm




# chmod +x /sc/arion/scratch/giottb01/HiCPro/HiCPro_step1_HiC_HChang.sh
# cd /sc/arion/scratch/giottb01/HiCPro/
# bsub -J HiC_step1 -P acc_Tsankov_Normal_Lung -q premium -n 8 -W 12:00 -R rusage[mem=16000] -R span[hosts=1] -o HiC.out -e HiC.err HiCPro_step1_HiC_HChang.sh
# pip install --upgrade pip setuptools wheel cython

# conda install bx-python=0.8.11 


# conda install bioconda::samtools
# conda install bioconda::bowtie2


# #conda install iced=0.5.10 
# conda install cooler=0.8.11 
# conda install r-base=4.0.3 
# conda install r-ggplot2=3.3.5 
# conda install r-rcolorbrewer=1.1_2 
# conda install r-gridbase=0.4_7 
# conda install tbb=2020.2 
# conda install bowtie2=2.4.4 
# conda install samtools=1.12 
# 







# hicpro_file <- HicproFile(matrix_files, bed = bed_files)
# hic = import(hicpro_file, focus = "chr2:156312347−156492348", resolution = 5000)

# #hic = import(hic_file, focus = "2:157160718-157446103", resolution = 5000)


# pdf (file.path('Plots',paste0('HiC_data_NR4A2_locus_',hfile,'_heatmap.pdf')))
# print(plotMatrix(hic, maxDistance = 300000,use.scores = 'balanced',
#     cmap = afmhotrColors()))
# dev.off()
# ## Horizontal matrix
# library (HiContacts)
# plotMatrix(
#     refocus(hic, 'II'),
#     use.scores = 'balanced', limits = c(-4, -1), 
#     maxDistance = 200000
# )


