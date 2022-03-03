#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 512M
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -J 06A_cstacks
#SBATCH -o stdout_06A_cstacks.txt
#SBATCH -e stderr_06A_cstacks.txt

#make sure you have a popmap file to go with this code
# before you try and submit it as a job
# The popmap file should look something like this:
# P08899	North_Sea
# P08905	North_Sea
# P08913	North_Sea
# P08915	North_Sea
# P08901	North_Sea
# P08903	North_Sea
# P08909	North_Sea
# P08912	North_Sea
# P08916	North_Sea
# P08918	North_Sea
# ...
# ...

#load modules required
module purge
module load python/v3.6.9
module load bwa/v0.7.17
#module load samtools/v1.9
module load stacks/v2.3b

# See more instructions here
#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gbs-data-snp-calling-using-stacks.html#gsc.tab=0

# also check out this website
#https://catchenlab.life.illinois.edu/stacks/manual/#prun

WD=$(pwd)

#define input 03 directories
D03a="03a_barcodes"
D03b="03b_demux"
D03c="03c_alignments"
D03d="03d_gstacks"
D03e="03e_population"
D03="03_stacks"
#define input 02 directory w reference genome
IND02="02_genome"
#define file with reference genome
REFGF="GCF_010909765.2_sAmbRad1.1.pri_genomic.fna"
#define popmapfile
PMAPF="part06C_popmap.txt"
# define directory for popmap
D04="04_popmaps"
D04a="04a_popmap"
#remove previous versions of directory
rm -rf "${D04}"
# make the directory again
mkdir "${D04}"
cd "${D04}"
#copy the popmap file into this directory
cp "${WD}"/"${PMAPF}" .
#mkdir "${D04a}"

# from : https://catchenlab.life.illinois.edu/stacks/manual/#prun
#Section 4.4.1. De novo Data

# change directory back to main directory
cd "${WD}"
#
# Build loci de novo in each sample for the single-end reads only. If paired-end reads are available, 
# they will be integrated in a later stage (tsv2bam stage).
# This loop will run ustacks on each sample, e.g.
#   ustacks -f ./samples/sample_01.1.fq.gz -o ./stacks -i 1 --name sample_01 -M 4 -p 8

#a separate population map only containing those samples.
	# cstacks -n 6 -P $src/stacks/ -M $src/popmaps/popmap -p 8
srun cstacks -n 6 -P "${WD}"/"${D03}"/"${D03d}"/ -M "${WD}"/"${D04}"/"${PMAPF}" -p 8


#


