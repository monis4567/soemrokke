#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 06A
#SBATCH -o stdout_06A_bedtools.txt
#SBATCH -e stderr_06A_bedtools.txt

#load modules required
module purge
module load htslib/v1.9
module load python/v3.6.9
module load perl/v5.32.0
module load stacks/v2.3b
module load bedtools/v2.29.0
module load samtools

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
#Define input directories for rawdata and ref_genome 
RWDD="rawdata"
RfGnmD="ref_genome"

# To run bedtools you will first need reduced version of the reference genome
# A reference genome deposited on NCBI GenBank also includes minor sequences that
# were impossible to map to the chromosome.
# These can be found here : https://www.ncbi.nlm.nih.gov/genome/?term=Amblyraja+radiata
# as "unplaced genomic scaffold"
# You already used wget previously to obtain a copy of the reference genome and placed this
# in the directory: D02="02_genome"
# Navigating into this directory and executing this command:

# $ cat GCF_010909765.2_sAmbRad1.1.pri_genomic.fasta | grep '>' | cut -c1-120

# Will allow you to see all fasta headers in the file. Like this: 

# >NC_045956.1 Amblyraja radiata isolate CabotCenter1 chromosome 1, sAmbRad1.1.pri, whole genome shotgun sequence
# >NC_045957.1 Amblyraja radiata isolate CabotCenter1 chromosome 2, sAmbRad1.1.pri, whole genome shotgun sequence
# >NC_045958.1 Amblyraja radiata isolate CabotCenter1 chromosome 3, sAmbRad1.1.pri, whole genome shotgun sequence
# ...
# ...
# >NC_046002.1 Amblyraja radiata isolate CabotCenter1 chromosome 47, sAmbRad1.1.pri, whole genome shotgun sequence
# >NC_046003.1 Amblyraja radiata isolate CabotCenter1 chromosome 48, sAmbRad1.1.pri, whole genome shotgun sequence
# >NC_046004.1 Amblyraja radiata isolate CabotCenter1 chromosome 49, sAmbRad1.1.pri, whole genome shotgun sequence
# >NC_045999.1 Amblyraja radiata isolate CabotCenter1 chromosome X, sAmbRad1.1.pri, whole genome shotgun sequence
# >NW_022630093.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri S100, whole genome shot
# >NW_022630094.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri S103, whole genome shot
# >NW_022630095.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri S104, whole genome shot
# >NW_022630096.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri S105, whole genome shot
# >NW_022630097.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri S106, whole genome shot
# ...
# ...
# >NW_022630987.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_992_ctg1, whol
# >NW_022630988.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_993_ctg1, whol
# >NW_022630989.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_994_ctg1, whol
# >NW_022630990.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_995_ctg1, whol
# >NW_022630991.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_996_ctg1, whol
# >NW_022630992.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_997_ctg1, whol
# >NW_022630993.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_998_ctg1, whol
# >NW_022630994.1 Amblyraja radiata isolate CabotCenter1 unplaced genomic scaffold, sAmbRad1.1.pri scaffold_999_ctg1, whol


# Before you can continue with bedtools and map the RADseq sequences you have obtained, you will
# need to exclude these "unplaced genomic scaffold" part of the reference genome
# This can be done with the section below that greps for the chromosomal parts
# and thereby excludes the unplaced genomic scaffold
# The resulting output file can then be used for mapping the RAD seq sequences to the chromosomal
# pieces in the reference genome 

# See here how to use sed command to print everythin up until pattern is encountered
#https://www.theunixschool.com/2011/09/sed-selective-printing.html
# and use  sed '$d'
# To remove the last line that contains the pattern you used sed to quit at
# Navigate to the directory where the reference genome is placed
cd $WD/$D05/$RfGnmD
# edit the file name to get a new file name
NwNmREFGF=$(echo $REFGfasta | sed 's;.fasta;.without_unplaced_genomic_scaffold.fasta;g')
# Then use sed command to only retain lines in the file up until pattern
# and then use sed to also delete the very last line, that contains the pattern
# NOTE !! This approach assumes that all unplaced genomic scaffolds are placed in the lower part of the file,
# and that all chromosomes are placed above in the reference genome file
cat $REFGfasta | sed '/unplaced genomic scaffold/q' | \
	# also quit when encountering a fasta header line that says 'unlocalized genomic scaffold'
	sed '/unlocalized genomic scaffold/q' | \
	# delete the last line of the file - i.e. the one that holds the line with the text string you searched for
 	sed '$d' > $NwNmREFGF

# Then you also need a bed file which has the accession numbers of the individual chromosomes in
# the reference genome file. But you only require the accession numbers for the chromosomes that have been
# assembled. Which means accession numbers for the 'unplaced genomic scaffold' should be excluded

# The bed file needs to look this:
# NC_045956.1	1	112	193168416
# NC_045957.1	1	195583132	347156831
# NC_045958.1	1	349051615	480366747
# NC_045959.1	1	482008299	603823035
# ...


# edit the file name to get a new file name
LIST2=$(echo $REFGfasta | sed 's;.fasta;.list_of_acc_numbers_to_keep.txt;g')
bedFLnm=$(echo $REFGfasta | sed 's;.fasta;.fasta_02.bed;g')
# First get the accession numbers from the modified fasta file where the unplaced_genomic_scaffold
# grep for the '>' symbol, then cut by space to get the accession number, and then remove the '>' sign
cat $NwNmREFGF | grep '>' | cut -f1 -d' ' | sed 's;>;;g' > "$LIST2"

# edit the file name to get the '.fai' file name
faiFNM=$(echo $REFGfasta | sed 's;.fasta;.fasta.fai;g')
# From the '.fai' file that has been generated by the paleomix code and placed in
# $WD/$D05/$RfGnmD
# Use grep with a file as an input to only retain the shortlisted accession numbers # see : https://stackoverflow.com/questions/19380925/grep-a-large-list-against-a-large-file
grep -f "$LIST2" "$faiFNM" > tmp01_shortlisted_fai.txt
# cut the shortlisted '.fai' file and make tmp files that can be pasted together
cut -f1 tmp01_shortlisted_fai.txt > tmp02_acc_nmbs.txt
cut -f2 tmp01_shortlisted_fai.txt > tmp03_acc_ends.txt
cut -f3 tmp01_shortlisted_fai.txt > tmp03_acc_start.txt
# get the  number of lines in the shortlisted '.fai' file
NL=$(wc -l tmp01_shortlisted_fai.txt | cut -f1 -d ' ')
# you will need an equal number of '1' in a column to paste together with the other tmp files
echo "1" | awk '{for(i=1; i<=n; i++) print}' n=$NL > tmp04_repeated_1.txt
#paste all tmp files together
paste tmp02_acc_nmbs.txt tmp04_repeated_1.txt tmp03_acc_ends.txt > tmp05_soemrokke.bed

# move the resulting bed file
mv tmp05_soemrokke.bed $WD/$D05/$RfGnmD/$bedFLnm
# and remove the tmp files
rm tmp*
#
cd $WD/$D05
# Once the reference genome has had the unassigend scaffolds removed . The mapping of RADseq to
# the reference genome can begin

# for file in *.bam
# do
# 	intersectBed -abam $file -b ../prefixes_03/repeatMask_BelugaHiC.bed -v | \
# 	 samtools view -h - | \
# 	  samtools view -bSh - | \
# 	  samtools view -b -L ../prefixes_03/22firstscaffolds.bed - > repeat_mask_bams/${file%.*}_RM.bam
# done

cd $WD/$D05
# Loop for every file in the directory folder
for file in *.2.bam
do
	#intersectBed -abam $file -b /groups/hologenomics/andreasp/data/mapping/ref_genome/GCF_010909765.2_sAmbRad1.1.pri_genomic.fasta.bed - > ${file%.*}.2.clean.bam
	intersectBed -abam $file -b $WD/$D05/$RfGnmD/"$bedFLnm" > ${file%.*}.2.clean.bam
done
# #-------------------------------------------------------------------------------------------------------------------------------------------------------

# # Loop for every file in directory folder (As job)
# for file in *.2.bam
# do
# 	xsbatch -c 1 --mem-per-cpu 2000 -J intersect1 -- "intersectBed -abam $file -b /groups/hologenomics/andreasp/data/mapping/ref_genome/GCF_010909765.2_sAmbRad1.1.pri_genomic.fasta.bed - > ${file%.*}.2.clean.bam"
# 	xsbatch -c 1 --mem-per-cpu 2000 -J intersect1 -- "intersectBed -abam $file -b $WD/$D05/$RfGnmD/"$bedFLnm" > ${file%.*}.2.clean.bam"
# done

# #-------------------------------------------------------------------------------------------------------------------------------------------------------

# # Test run on a single file ### Working
# bedtools intersect -abam P08713.GCF_010909765.2.bam -b /groups/hologenomics/andreasp/data/mapping/ref_genome/GCF_010909765.2_sAmbRad1.1.pri_genomic.fasta.bed > P08713.GCF_010909765.2.clean.bam






# for file in *.bam
# do
# 	intersectBed -abam $file -b ../prefixes_03/repeatMask_BelugaHiC.bed -v | \
# 	 samtools view -h - | \
# 	  samtools view -bSh - | \
# 	  samtools view -b -L ../prefixes_03/22firstscaffolds.bed - > repeat_mask_bams/${file%.*}_RM.bam
# done

#xsbatch -c 1 --mem-per-cpu 16000 -J intersect1 -- "bedtools intersect -abam Bris.Delphinapterus_leucas_ASM228892v2_HiC.bam -b /groups/hologenomics/xcl300/data/beluga_modern/prefixes_03/21_Chr_autosome.bed > repeat_mask_bams/Bris.Delphinapterus_leucas_ASM228892v2_HiC.temp6_bam"

#srun "bedtools intersect -abam Bris.Delphinapterus_leucas_ASM228892v2_HiC.bam -b /"$WD"/prefixes_03/21_Chr_autosome.bed > repeat_mask_bams/Bris.Delphinapterus_leucas_ASM228892v2_HiC.temp6_bam"


#