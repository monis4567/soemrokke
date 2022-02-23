#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 4G
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH -J 02B_indxRefGenome
#SBATCH -o stdout_02B_indxRefGenome.txt
#SBATCH -e stderr_02B_indxRefGenome.txt

#load modules required
module purge
#module load python/v2.7.12
#module load cutadapt/v1.11
#module load vsearch/v2.8.0
#module load stacks/v2.3b
# module load python/v3.6.9
# module load perl/v5.32.0
module load bwa/v0.7.17
# get working directory
WD=$(pwd)

#define input 03 directories
D03a="03a_barcodes"
D03b="03b_demux"
D03c="03c_alignments"
D03d="03d_gstacks"
D03e="03e_population"
D03="03_stacks"
#define input 02 directory w reference genome
D02="02_genome"
#define file with reference genome
REFGF="GCF_010909765.2_sAmbRad1.1.pri_genomic.fna"
#make new file name to write output to
NRFG=$(echo $REFGF | sed "s;.fna;.index.fa;g")
# The bwa tool has a manual page
# https://manpages.org/bwa
# Trying to follow some directions from this website
# https://kuriouskevin.wordpress.com/my-code/mapping-reads-to-a-reference-genome-using-bwa/
# perhaps this website also has some helpful hints
# https://icb.med.cornell.edu/wiki/index.php/Elementolab/BWA_tutorial
#run bwa
bwa index -a bwtsw "${WD}"/"${D02}"/"${REFGF}" > "${WD}"/"${D02}"/"${NRFG}"