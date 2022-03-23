#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 512M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 08A_angsd
#SBATCH -o stdout_08A_angsd.txt
#SBATCH -e stderr_08A_angsd.txt

#load modules required
module purge
module load lib
module load java/v1.8.0_202-jdk python/v3.6.9 R/v3.6.1 htslib samtools bwa bowtie2 AdapterRemoval mapDamage picard paleomix 
module load angsd/v0.935

WD=$(pwd)
OUTDIR08="08_angsd"
rm -rf "${OUTDIR08}"
mkdir "${OUTDIR08}"

OUTDIR06="06_refgenome_mapped_bam"
cd "${OUTDIR06}"
# make a list that holds the names of the fastq files
LSSMPL=$(ls | uniq)
#make the list of samples an array you can iterate over
declare -a SMPLARRAY=($LSSMPL)
 
#change back to the working dir
cd "$WD"

touch "$WD"/"${OUTDIR08}"/cl_BAMlist.txt
iconv -f UTF-8 -t UTF-8 "$WD"/"${OUTDIR08}"/cl_BAMlist.txt

#iterate over samples
for smp in ${SMPLARRAY[@]}
do
	echo ""${WD}"/"${OUTDIR06}"/"${smp}"" >> "$WD"/"${OUTDIR08}"/cl_BAMlist.txt
done

cd "$WD"/"${OUTDIR08}"
angsd -b "$WD"/"${OUTDIR08}"/cl_BAMlist.txt -GL 2 -minQ 25 -minmapQ 25 -minMaf 0.05 -doGlf 2 -doCounts 1 -remove_bads 1 -doMaf 1 -doMajorMinor 1 -doDepth 1 -maxDepth 10000 -dumpCounts 2 -minInd 20 
#-baq 2 -C 50




#