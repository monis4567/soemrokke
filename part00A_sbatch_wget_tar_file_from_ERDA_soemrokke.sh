#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 1
#SBATCH -t 23:00:00
#SBATCH -J 00A_fetch_file
#####  ##SBATCH -o stdout_00A_fetch_file.txt
#####  ##SBATCH -e stderr_00A_fetch_file.txt

#load modules required
module purge
#module load python/v2.7.12
#module load cutadapt/v1.11
#module load vsearch/v2.8.0

#Get file from ERDA
#
wget --no-check-certificate 'https://sid.erda.dk/share_redirect/Pj3MNl2bJo/UO_C1246_1.fastq.gz' .
#https://sid.erda.dk/cgi-sid/ls.py?share_id=Pj3MNl2bJo
#https://sid.erda.dk/share_redirect/Pj3MNl2bJo/UO_C1246_1.fastq.gz

