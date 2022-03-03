#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 12
#SBATCH -t 1:00:00
#SBATCH -J 05B_paleomix
#SBATCH -o stdout_05B_paleomix.txt
#SBATCH -e stderr_05B_paleomix.txt

#https://paleomix.readthedocs.io/en/stable/bam_pipeline/requirements.html
# Check this website to see how you can get picard available in your home directory
# you will need to run lines :
# $ mkdir -p ~/install/jar_root
# $ wget -O ~/install/jar_root/picard.jar https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar
# Prior to running paleomix

WD=$(pwd)
#load modules required
module purge
#module load python/v2.7.12
#module load cutadapt/v1.11

module load java/v1.8.0_202-jdk python/v3.6.9 R/v3.6.1 htslib samtools bwa bowtie2 AdapterRemoval mapDamage picard paleomix 

# If you want to check out some options for getting help and preparing an example then load all the modules above and then run
# $paleomix bam
# In a HPC terminal. This will show you some help on paleomix

#define input 03 directories
D03a="03a_barcodes"
D03b="03b_demux"
D03c="03c_alignments"
D03d="03d_gstacks"
D03e="03e_population"
D03="03_stacks"

cd "$WD"
# define paleomix output directory
D05="05_paleomix_mapping"
rm -rf "$D05" 
mkdir "$D05"
# change directory
cd "$D05"
#Define input directories for rawdata and ref_genome 
RWDD="rawdata"
RfGnmD="ref_genome"
mkdir "$RWDD"
mkdir "$RfGnmD"

#define input 02 directory w reference genome
D02="02_genome"
#define file with reference genome her with a fna ending
REFGF="GCF_010909765.2_sAmbRad1.1.pri_genomic.fna"
 
#make new file name to copy the fna file to 
REFGF_fasta=$(echo $REFGF | sed "s;.fna;.fasta;g")


cd $WD/"$D05"
#Define input and output f
INtable="paths.txt"
OUTyaml="part05_makefile_soemrokke_modified.yaml"

# change directory back to 05_paleomix directory
cd "$WD"/"$D05"/
#check with dryrun if the 'OUTyaml' file is able to run 
#paleomix bam dryrun "$OUTyaml"
# you can inspect the file stderr_05B_paleomix.txt
# first before you do anything
# if it can run then comment out the line above
# increase the cpu usage to -c 12 and --mem_per cpu to 4000M and --max_threads to 12
# in the SBATCH parts in the first lines above
# and the instead activate the line here below
paleomix bam run --max-threads 12 "$OUTyaml"



#