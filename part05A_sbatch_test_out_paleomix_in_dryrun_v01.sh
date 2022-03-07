#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 128M
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J 05A_paleomix
#SBATCH -o stdout_05A_paleomix.txt
#SBATCH -e stderr_05A_paleomix.txt

# # The overall idea is to set up the directories as sketched out below
# soemrokke/
# ├── 01_data
# ├── 02_genome
# └── 03_stacks
# |   ├── 03a_barcodes
# |   ├── 03b_demux
# |   ├── 03c_alignments
# |   ├── 03d_gstacks
# |   └── 03e_population
# └── 04_popmap
# └── 05_paleomix_mapping
# 	├── rawdata
# 	└── ref_genome

# # The 05A code will generate the directory:
# └── 05_paleomix_mapping
#  and then make the directories:
# 	├── rawdata
# 	└── ref_genome
# # The code 05A will then copy all demultiplexed ".fq.gz" file from  
#     ├── 03b_demux
# # and then add them the directory: "rawdata"
# # The code 05A will then copy the reference genome and rename the file ending from ".fna" to ".fasta" file from  
# ├── 02_genome
# # and place it in the directory : "ref_genome"
# # From the ".fq.gz" files now placed in the directory: "rawdata" 
# # The code 05A will then make a list of the 'fq.gz' files and get the sample names. This will generate a 'paths.txt' file
# # which will be stored in: 
# └── 05_paleomix_mapping
# Once the modules are loaded the paleomix can be loaded and the command: 
#   `paleomix bam makefile`
#   can be used to make a default generic '.yaml' file that serves as a makefile for running the paleomix
#   The part05 code then modfies this generic '.yaml' file to make it match the thorny skate data. And part05A then starts a 'dryrun' in paleomix, to check out whether the settings are correct. Inspect the part05A code for details on what is required and what the different replacement steps with the sed command are good for. Also inspect the result from the dryrun, before continuing on to runnung part05B. 


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

# change directory  to where the demulitplexed files are stored
cd $WD/$D03/$D03b
#make list of 'fq.gz' the files
Lfq=$(ls $PWD/* | grep 'fq.gz' | sed 's;^.*/;;g')
#make it an array
declare -a SAR=($Lfq)
# iterate over the elements in this list
# to copy the demultiplexed 'fq.gz' files into the rawdata directory
for file in ${SAR[@]}
do
	cp $file "$WD"/"$D05"/"$RWDD"/.
done

#change directory to where the rawdata files are located
cd "$WD"/"$D05"/"$RWDD"/
# get lists of all 'fq.gz' files in directory
# and us the paths to create a 'paths' file
ls $PWD/* | grep 'fq.gz' | sed 's;^.*/;;g' | sed 's;.fq.*;;g' > tmp1.txt
ls $PWD/* | grep 'fq.gz' | sed "s://:/:g" > tmp2.txt
paste tmp1.txt tmp1.txt tmp1.txt tmp2.txt > paths.txt
#remove temporary files
rm tmp1.txt
rm tmp2.txt
#move the paths file to the '05_paleomix_mapping' directory
mv $WD/$D05/"$RWDD"/paths.txt $WD/"$D05"/paths.txt

#define input 02 directory w reference genome
D02="02_genome"
#define file with reference genome here with a fna ending
REFGF="GCF_010909765.2_sAmbRad1.1.pri_genomic.fna"
# get first part of ref genome name
fRFGNm=$(echo $REFGF | sed 's/\(.*\)_\(.*\)_\(.*\)/\1/g')
#make new file name to copy the fna file to 
REFGF_fasta=$(echo $REFGF | sed "s;.fna;.fasta;g")
#copy the reference genome
cp "$WD"/"$D02"/"$REFGF" $WD/"$D05"/"$RfGnmD"/"$REFGF_fasta"
# cahnge directory
cd $WD
# Use paleomix to write out a makefile you can modify
paleomix bam makefile > "$WD"/"$D05"/tmp_makefile_example.txt

cd "$WD"/"$D05"
#substitute using sed in the makefile
#write out the file
cat tmp_makefile_example.txt | \
# then replace the path for the reference genome
sed "s;Path: PATH_TO_PREFIX.fasta;Path: ./"$RfGnmD"/"$REFGF_fasta";g" | \

sed "s;Algorithm: backtrack;Algorithm: mem;g" | \

# replace to comment out lines
sed "s;NAME_OF_TARGET:;#NAME_OF_TARGET:;g" | \
sed "s;  NAME_OF_SAMPLE:;#  NAME_OF_SAMPLE:;g" | \
sed "s;    NAME_OF_LIBRARY:;#    NAME_OF_LIBRARY:;g" | \
sed "s;      NAME_OF_LANE: PATH_WITH_WILDCARDS;#      NAME_OF_LANE: PATH_WITH_WILDCARDS;g" | \
#replace 'NAME_OF_PREFIX' 
sed "s;NAME_OF_PREFIX:;${fRFGNm}:;g" | \

sed "s;PCRDuplicates: filter;PCRDuplicates: no;g" > $WD/"$D05"/part05_makefile_soemrokke.yaml
rm "$WD"/"$D05"/tmp_makefile_example.txt

cd $WD/"$D05"
#Define input and output f
INtable="paths.txt"
OUTyaml="part05_makefile_soemrokke_modified.yaml"
# copy and modify the yaml file
cp part05_makefile_soemrokke.yaml  "$OUTyaml"
ant1=""
ant2=""
# write out file with paths, and define column headers named: "sample" "flowcell" "lane" "fastq"
cat "$INtable" | while read sample flowcell lane fastq
do
	if [[ $sample != $ant1 ]]; then
	echo -e "${sample}: " >> "$OUTyaml"
	echo -e "  ${sample}: " >> "$OUTyaml"
	fi
	if ! ([[ $sample == $ant1 ]] && [[ $flowcell == $ant2 ]]); then
	echo -e "    Lane_${flowcell}: " >> "$OUTyaml"
	fi
	echo -e "      Lib_${lane}: ${fastq}" >> "$OUTyaml"
	ant1=$sample
	ant2=$library
done 
# change directory back to 05_paleomix directory
cd "$WD"/"$D05"/
#check with dryrun if the 'OUTyaml' file is able to run 
paleomix bam dryrun "$OUTyaml"
# you can inspect the file stderr_05A_paleomix.txt
# first before you do anything
# if it can run then comment out the line above
# increase the cpu usage to -c 12 and --mem_per cpu to 4000M and --max_threads to 12
# in the SBATCH parts in the first lines above
# and the instead activate the line here below
# paleomix bam run "$OUTyaml"

# If it worked you should see something similar looking to this:

# 21:25:06 INFO Reading makefile 'part05_makefile_soemrokke_modified.yaml'
# 21:25:06 INFO Validating FASTA files
# 21:25:06 INFO Indexing FASTA at './ref_genome/GCF_010909765.2_sAmbRad1.1.pri_genomic.fasta'
# 21:25:42 INFO Building BAM pipeline for 'part05_makefile_soemrokke_modified.yaml'
# 21:25:42 INFO Running BAM pipeline
# 21:25:42 INFO Checking file dependencies
# 21:25:42 INFO Checking for auxiliary files
# 21:25:42 INFO Checking required software
# 21:25:42 INFO  - Found Rscript v3.6.1
# 21:25:42 INFO  - Found AdapterRemoval v2.3.0
# 21:25:42 INFO  - Found BWA v0.7.17
# 21:25:46 INFO  - Found Picard tools v2.18
# 21:25:49 INFO  - Found mapDamage v2.2.1
# 21:25:49 INFO  - Found samtools v1.9.0
# 21:25:49 INFO Determining states
# 21:25:49 INFO Ready
# 21:25:49 WARNING One or more tasks require more threads than the user-defined maximum; the number of threads used will therefore exceed the user-defined maximum when running those tasks.
# 21:25:49 INFO Saving warning logs to '/shared/volume/hologenomics/data/phq599/soemrokke/05_paleomix_mapping/bam_pipeline.20220303_212506_01.log'
# 21:25:49 INFO Number of nodes:             1238
# 21:25:49 INFO Number of done nodes:        1
# 21:25:49 INFO Number of runable nodes:     96
# 21:25:49 INFO Number of queued nodes:      1141
# 21:25:49 INFO Number of outdated nodes:    0
# 21:25:49 INFO Number of failed nodes:      0
# 21:25:49 INFO Pipeline completed successfully
# 21:25:49 INFO Dry run done
# 21:25:49 INFO Log-file written to '/shared/volume/hologenomics/data/phq599/soemrokke/05_paleomix_mapping/bam_pipeline.20220303_212506_01.log'

#