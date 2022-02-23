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

REF=$(echo ""$WD"/"${IND02}"/"${REFGF}"")


# from : https://catchenlab.life.illinois.edu/stacks/manual/#prun
#Section 4.4.1. De novo Data
###!/bin/bash

src=$HOME/research/project

files=”sample_01
sample_02
sample_03”

#
# Build loci de novo in each sample for the single-end reads only. If paired-end reads are available, 
# they will be integrated in a later stage (tsv2bam stage).
# This loop will run ustacks on each sample, e.g.
#   ustacks -f ./samples/sample_01.1.fq.gz -o ./stacks -i 1 --name sample_01 -M 4 -p 8
#
id=1
for sample in $files
do
    ustacks -f $src/samples/${sample}.1.fq.gz -o $src/stacks -i $id --name $sample -M 4 -p 8
    let "id+=1"
done

# 
# Build the catalog of loci available in the metapopulation from the samples contained
# in the population map. To build the catalog from a subset of individuals, supply
# a separate population map only containing those samples.
#
cstacks -n 6 -P $src/stacks/ -M $src/popmaps/popmap -p 8

#
# Run sstacks. Match all samples supplied in the population map against the catalog.
#
sstacks -P $src/stacks/ -M $src/popmaps/popmap -p 8

#
# Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include
# paired-end reads using tsv2bam. tsv2bam expects the paired read files to be in the samples
# directory and they should be named consistently with the single-end reads,
# e.g. sample_01.1.fq.gz and sample_01.2.fq.gz, which is how process_radtags will output them.
#
tsv2bam -P $src/stacks/ -M $src/popmaps/popmap --pe-reads-dir $src/samples -t 8

#
# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.
#
gstacks -P $src/stacks/ -M $src/popmaps/popmap -t 8

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.
#
populations -P $src/stacks/ -M $src/popmaps/popmap -r 0.65 --vcf --genepop --structure --fstats --hwe -t 8


#

#