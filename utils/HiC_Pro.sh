#!/bin/bash
#BSUB -J chrBP_model
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=2
#BSUB -R h100nvl
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml anaconda3/2022.10
# ml cuda/11.7.0
# ml cudnn/8.9.5-11
# ml proxies
# ml java/11.0.2
# ml tensorrt/8.5.3.1


#source activate chrombpnet

source activate hic
# hicpro is actually already installed on minerva (but it doesnt work). Need to clone HiC-Pro repo then make configure then cp manually to ~/local/bin/HiC-Pro_3.1.0 and then edit config file to change paths to ~/local/bin/HiC-Pro_3.1.0 and then do make install
#ml hicpro

~/local/bin/HiC-Pro_3.1.0/bin/HiC-Pro -i /sc/arion/scratch/leew17/HiC/rawdata -o /sc/arion/scratch/giottb01/HiCPro2 -c /sc/arion/projects/Tsankov_Normal_Lung/Bruno/HiCPro_files/config_file.txt

