#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 1024M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 07A
#SBATCH -o stdout_07A_gstacks.txt
#SBATCH -e stderr_07A_gstacks.txt

# run ref_map.pl


# stacks website
#https://catchenlab.life.illinois.edu/stacks/manual/#prun
#load modules required
module purge
# Loading modules 
module load stacks/v2.3b

WD=$(pwd)

#define input 02 directory w reference genome
D02="02_genome"
#define file with reference genome here with a fna ending
REFGF="GCF_010909765.2_sAmbRad1.1.pri_genomic.fna"
REFGfasta=$(echo $REFGF | sed 's/.fna/.fasta/g')
#define input 03 directories
D03a="03a_barcodes"
D03b="03b_demux"
D03c="03c_alignments"
D03d="03d_gstacks"
D03e="03e_population"
D03="03_stacks"
# define paleomix output directory
D05="05_paleomix_mapping"
# define paleomix output directory
D06="06_refgenome_mapped_bam"

#Define input directories for rawdata and ref_genome 
RWDD="rawdata"
RfGnmD="ref_genome"
#define popmaps directory
D07="07_popmaps"
rm -rf "${WD}"/"${D07}"/
mkdir "${WD}"/"${D07}"/

rm -rf "${WD}"/"${D06}"/
mkdir "${WD}"/"${D06}"/

#Remove any previous versions of gstacks results
cd "${WD}"/"${D03}"/"${D03d}"
rm *

cd "${WD}"/"${D06}"/
cd "${WD}"/"${D05}"/
# make a list that holds the names of the clean bam files
LSSMPL=$(ls *.2.clean.bam)
#ls *.2.clean.bam | awk 'BEGIN { FS = "." } ; {print $1}'
#make the list of samples an array you can iterate over
declare -a SMPLARRAY=($LSSMPL)
#iterate over samples
for f in ${SMPLARRAY[@]}
do
    nf=$(echo $f | awk 'BEGIN { FS = "." } ; {print $1".bam"}')
    cp "${WD}"/"${D05}"/$f "${WD}"/"${D06}"/$nf
done

cd "$WD"
popmapFl="part07B_popmap.txt"
cp "${popmapFl}" "${WD}"/"${D07}"/.

cd "${WD}"/"${D06}"/
#rm FGXCONTROL*
#FGXCONTROL.bam
#cd "${WD}"/"${D03}"/"${D03d}"
#rm catalog*

# Run gstacks to build loci from the aligned paired-end data. We have instructed 
# gstacks to remove any PCR duplicates that it finds. 


# Gstacks - Processing single-end data, already alligned to a reference genome. 
# Using 2.clean.bam files from the folder after Bedtools comparison and sorting. 
# Clean .bam files are moved into a seperate folder where gstacks is ran, as stacks compare all files in the folder, only have the wanted files present. 


# Running gstacks to build loci from the aligned single-end data.
#gstacks -I /groups/hologenomics/andreasp/data/mapping/clean_clavata/ -M /groups/hologenomics/andreasp/data/mapping/clavata_popmap.txt -S .GCF_010909765.2.clean.bam -O /groups/hologenomics/andreasp/data/mapping/clavata_stacks/ -t 8
gstacks -I "${WD}"/"${D06}"/ -M "${WD}"/"${D07}"/"$popmapFl" -O "${WD}"/"${D03}"/"${D03d}" -t 8

## I: Is the input folder, where all the clean .bam files are located.
## M: Is the path to the popmap, the popmap consists of the prefix ID with an allocating coulumn with area. In this case either NS = North Sea, SK = Skagerrak, KG = Kattegat or AQ = Aquarium.
## O: Is the path where stacks places the output.
## S: Is the suffix needed. Stacks assume the files to in the format "sample_number.bam".

# Running populations. 
# Calculating Hardy-Weinberg deviation, population statistics, f-statistics and # smooth the statistics across the genome. Export several output files.
#populations -P /groups/hologenomics/andreasp/data/mapping/clavata_stacks/ -M /groups/hologenomics/andreasp/data/mapping/clavata_popmap.txt -r 0.65 --vcf --genepop --fstats --smooth --hwe -t 8
populations -P "${WD}"/"${D03}"/"${D03d}" -M "${WD}"/"${D07}"/"$popmapFl" -r 0.65 --vcf --genepop --fstats --smooth --hwe -t 8



# It is assumed that your files are named properly in the population map and on the file system. 
# So, for a paired-end analysis, given sample_2351 listed in the population map, denovo_map.pl expects 
# to find files sample_2351.1.fq.gz and sample_2351.2.fq.gz in the directory specified with --samples.
# Here is an example running ref_map.pl for a paired-end population analysis:

# ref_map.pl -T 8 --popmap ./popmaps/popmap -o ./stacks/ --samples ./aligned

# ref_map.pl -T 8 --popmap "${WD}"/"${D07}"/"${popmapFl}" -o "${WD}"/"${D03}"/"${D03e}" --samples "${WD}"/"${D05}"/

# ref_map.pl will read the file names out of the population map and look for them in
# the directory specified with --samples. The ref_map.pl program expects, 
# given sample_2351 listed in the population map, to find a sample_2351.bam file
# containing both single and paired-reads aligned to the reference genome and sorted. 

# ref_map.pl will read the file names out of the population
# map and look for them in the directory specified with --samples. 
# The ref_map.pl program expects, given sample_2351 listed in the population map,
#  to find a sample_2351.bam file containing both single and paired-reads aligned to the reference genome and sorted.




#Here is an example shell script for reference-aligned data that uses shell loops to easily execute the pipeline by hand:
####!/bin/bash

# src=$HOME/research/project
# bwa_db=$src/bwa_db/my_bwa_db_prefix
    
# files=”sample_01
# sample_02
# sample_03”

# #
# # Align paired-end data with BWA, convert to BAM and SORT.
# #
# for sample in $files
# do 
#     bwa mem -t 8 $bwa_db $src/samples/${sample}.1.fq.gz $src/samples/${sample}.2.fq.gz |
#       samtools view -b |
#       samtools sort --threads 4 > $src/aligned/${sample}.bam
# done

# #
# # Run gstacks to build loci from the aligned paired-end data. We have instructed
# # gstacks to remove any PCR duplicates that it finds.
# #
# gstacks -I $src/aligned/ -M $src/popmaps/popmap --rm-pcr-duplicates -O $src/stacks/ -t 8

# #
# # Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics and 
# # smooth the statistics across the genome. Export several output files.
# #
# populations -P $src/stacks/ -M $src/popmaps/popmap -r 0.65 --vcf --genepop --fstats --smooth --hwe -t 8

#