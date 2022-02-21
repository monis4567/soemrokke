#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 04A_algn_wRefGenome
#SBATCH -o stdout_04A_algn_wRefGenome.txt
#SBATCH -e stderr_04A_algn_wRefGenome.txt

#load modules required
module purge
#module load python/v2.7.12
#module load cutadapt/v1.11
#module load vsearch/v2.8.0
#module load stacks/v2.3b
# module load python/v3.6.9
# module load perl/v5.32.0
module load bwa/v0.7.17

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
#The next step will equal step:
#4.2. Align data against a reference genome
# in this protocol
# https://catchenlab.life.illinois.edu/stacks/manual/#prun

# create a run script
printf "input=\$1 
output=\$(basename \$input | sed 's/.fq.gz/.sam/g')
index=\$(echo \""${WD}"/"${IND02}"/"${REFGF}"\")
bwa mem -t 4 \${index} \${input} > ./"${D03}"/"${D03c}"/\${output}" > part04B_runBWA.sh

#cat part04B_runBWA.sh
for fq in $(find "$WD"/"$D03"/"$D03b" -name "*.gz"); do
  echo ./part04B_runBWA.sh $fq;
done > bwa.cmds

# #makeSLURMp.py 350 bwa.cmds
# for sub in bwa*.sub; do
#   sed -i 's/parallel -j 1/parallel -j 9/g' $sub;
#   sbatch $sub
# done




#
