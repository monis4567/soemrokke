#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 2G
#SBATCH -c 2
#SBATCH -t 6:00:00
#SBATCH -J 03A_prp_barcode_prc_radtags
#SBATCH -o stdout_03A_prp_barcode.txt
#SBATCH -e stderr_03A_prp_barcode.txt

#load modules required
module purge
#module load python/v2.7.12
#module load cutadapt/v1.11
#module load vsearch/v2.8.0
#module load stacks/v2.3b
module load python/v3.6.9
module load perl/v5.32.0
module load stacks/v2.3b


# See more instructions here
#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gbs-data-snp-calling-using-stacks.html#gsc.tab=0

# And more importantly check out this website
#https://catchenlab.life.illinois.edu/stacks/manual/#prun

#get the working directory
WD=$(pwd)
raw="${WD}"

#move the raw dataset - if not already done in the previous parts
#mv UO_C1246_1.fastq "${WD}"/01_data/.

# define the samples output directory
samples="03_stacks/03b_demux"
#define the name for the barcode index list file
barcode_index_file="part03B_barcode_index_list.txt"
#Then define a directory where the barcode index list file can go in
barcodes="03_stacks/03a_barcodes"
# then define the directory where you have the raw data files
raw="01_data"

# #Remove previous versions of directories
# rm -rf 01_data/
# rm -rf 02_genome/
rm -rf 03_stacks/
# #Make new directories
# mkdir 01_data
# mkdir 02_genome
mkdir 03_stacks

cd 03_stacks
#Make new directories
mkdir 03a_barcodes
mkdir 03b_demux
mkdir 03c_alignments
mkdir 03d_gstacks
mkdir 03e_population

# change back to working directory
cd "${WD}"
#move the barcode index list
mv "${WD}"/"${barcode_index_file}" "${WD}"/03_stacks/03a_barcodes/.
# try an echo line first to see if the right paths are called
#echo "./"${raw}"/ -o ./"${samples}"/ -b ./"${barcodes}"/"${barcode_index_file}" \\"

# Then start demultiplexing
# As described under section 4.1.1.1. here: https://catchenlab.life.illinois.edu/stacks/manual/#prun
# This should equal section 4.1.3 in this website
#
process_radtags --inline_null -p ./"${raw}"/ -o ./"${samples}"/ -b ./"${barcodes}"/"${barcode_index_file}" \
                   -e sbfI -r -c -q


#