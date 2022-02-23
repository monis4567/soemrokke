#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 05A_sort
#SBATCH -o stdout_05A_sort.txt
#SBATCH -e stderr_05A_sort.txt

#load modules required
module purge
#module load stacks/v2.3b
module load python/v3.6.9


# See more instructions here
#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gbs-data-snp-calling-using-stacks.html#gsc.tab=0

# also check out this website
#https://catchenlab.life.illinois.edu/stacks/manual/#prun


module load bwa/v0.7.17
module load samtools/v1.9

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

# Change dir to where you are about to place the genome
#cd "${WD}"/02_genome



SAM="$1"
samtools view --threads 36 -b -o ${SAM%.*}.bam ${SAM}
samtools sort -o ${SAM%.*}_sorted.bam -T ${SAM%.*}_temp ${SAM%.*}.bam
module load picard
REF="/path/to/stacks/2-genome/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
REF=$(echo ""$WD"/"${IND02}"/"${REFGF}"")

ulimit -c unlimited
# example file looks like this:
# Z002E0081-SRR391089_sorted.bam
bam="${SAM%.*}_sorted.bam"
RGID=$(basename $bam | rev | cut -f 1 -d "-" | rev | sed 's/_sorted.bam//g')
RGSM=$(basename $bam | rev | cut -f 2- -d "-" | rev )
RGLB="${RGSM}-L001"
RGPU=001
echo -e "$RGID\t$RGSM\t$RGLB\t$RGPU"
java -Djava.io.tmpdir=$TMPDIR -Xmx50G -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \
      I=${bam} \
      O=${bam%.*}_new.bam \
      RGID=$RGSM \
      RGLB=$RGLB \
      RGPL=ILLUMINA \
      RGPU=$RGPU \
      RGSM=$RGSM
module load samtools
samtools index ${bam%.*}_new.bam

#