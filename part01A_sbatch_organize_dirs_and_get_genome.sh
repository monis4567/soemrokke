#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 01A_get_genome
#SBATCH -o stdout_01A_get_genome.txt
#SBATCH -e stderr_01A_get_genome.txt

#load modules required
module purge
#module load python/v2.7.12
#module load cutadapt/v1.11
#module load vsearch/v2.8.0
#module load stacks/v2.3b
module load python/v3.6.9


# See more instructions here
#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gbs-data-snp-calling-using-stacks.html#gsc.tab=0

# And more importantly check out this website
#https://catchenlab.life.illinois.edu/stacks/manual/#prun


#get
WD=$(pwd)
#Remove previous versions of directories
rm -rf 01_data/
rm -rf 02_genome/
rm -rf 03_stacks/
#Make new directories
mkdir 01_data
mkdir 02_genome

# mkdir 03_stacks

# cd 03_stacks
# #Make new directories
# mkdir 03a_barcodes
# mkdir 03b_demux
# mkdir 03c_alignments
# mkdir 03d_gstacks
# mkdir 03e_population

# Change dir to where you are about to place the genome
cd "${WD}"/02_genome

#https://www.metagenomics.wiki/tools/fastq/ncbi-ftp-genome-download
#1) Download list of all available reference genomes
#download complete list of manually reviewed genomes (RefSeq database, subset of GenBank)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
#wait for the download
sleep 60s
#2) Search for available genomes of a species
#Example: Eubacterium rectale  (RefSeq database, check columns 8,9,14,15,16)
# Try and grep for all  with 'raja'
grep -E '.*raja.*' assembly_summary_refseq.txt | cut -f 8,9,14,15,16
# Try again but be more specific
grep -E 'Amblyraja.*' assembly_summary_refseq.txt | cut -f 8,9,14,15,16

# Find:
#Amblyraja radiata		Full	2020/02/21	sAmbRad1.1.pri

#3) Get FTP download link

# for selected genomes (Eubacterium rectale), get NCBI ftp download folder (column 20)
grep -E 'Amblyraja.*' assembly_summary_refseq.txt | cut -f 20 > ftp_folder.txt
#head ftp_folder.txt
# extend download folder: create an exact genome (fna or gff) download link
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_fna_files.sh
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftp_folder.txt > download_gff_files.sh
#head download_fna_files.sh

chmod 755 download_fna_files.sh
#4) Run download
#download the .fna genome files (fasta format)
source download_fna_files.sh

cd "$WD"
# Try and look for unique 5 
#https://unix.stackexchange.com/questions/369181/printing-every-nth-line-out-of-a-large-file-into-a-new-file
#https://stackoverflow.com/questions/13778273/find-unique-lines
# Get first 1200 lines from file
head -1200 first1200linesUO_C1246_1.fastq > first1200linesUO_C1246_1.fastq
grep -A1 '1:N:0:1' first1200linesUO_C1246_1.fastq | awk 'NR % 3 == 2' | cut -c-5 | sort | uniq -u > get_uniq_5_tags.txt
# decompress genome files
#gzip -d *.gz

#ls
# cd $WD
# cp "${WD}"/part01B_dwnld_Amblyraja_radiata_gnm_v01.py "${WD}"/02_genome/.


#get
#python3 part01B_dwnld_Amblyraja_radiata_gnm_v01.py

#wget ftp://ftp.gramene.org/pub/gramene/CURRENT_RELEASE/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
# gunzip Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
# module load bwa/v0.7.17
# bwa index Zea_mays.B73_RefGen_v4.dna.toplevel.fa

 

# cd 03_stacks
# while read line; do
#   grep -Fw '$line' barcodes-gbs-data.tsv > a-barcodes/${line}-barcodes.txt;
# done<../1-data/srr.ids


# cd "${WD}"/01_data
# esearch -db sra -query SRP009896  |\
#   efetch --format runinfo |\
#   cut -d , -f 1 |\
#   awk 'NF>0' |\
#   grep -v Run > srr.ids
# while read line; do
#   enaDataGet -f fastq -a $line;
# done<srr.ids

#