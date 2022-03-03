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
#The next step will equal step:
#4.2. Align data against a reference genome
# in this protocol
# https://catchenlab.life.illinois.edu/stacks/manual/#prun

# # Make test directories. Only needed to check out if the
# # iteration below works properly
# # Can be removed when set to run externally
# rm -rf $D03
# mkdir $D03
# cd $D03
# mkdir $D03b
# mkdir $D03c
# cd $D03b
# for fn in $(seq -f "%03g" 1 16)
# do
#   td=$(echo "td"${fn}".gz")
#   echo "text_string_"${fn}"" > $td
# done

#change to directory where the resulting files are to end up
cd "${WD}"/"${D03}"/"${D03c}"
#remove everything in this directory
rm -rf *
# Navigate back to working directory
cd "${WD}"
# # create a run script
# printf "input=\$1 
# output=\$(basename \$input | sed 's/.fq.gz/.sam/g')
# index=\$(echo \""${WD}"/"${IND02}"/"${REFGF}"\")
# bwa mem -t 4 \${index} \${input} > ./"${D03}"/"${D03c}"/\${output}" > part04C_runBWA.sh

#iterate over input demultiplexed gz files
for fq in $(find "$WD"/"$D03"/"$D03b" -name "*.gz")
do
  #get only file name , by using sed to get rid of directories   
  Dfq=$(echo $fq | sed "s;${WD}/${D03}/${D03b}/;;g") 

# Also see thius website
# https://eriqande.github.io/eca-bioinf-handbook/alignment-of-sequence-data-to-a-reference-genome-and-associated-steps.html

# Following part of this guide
# https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/gbs-data-snp-calling-using-stacks.html#gsc.tab=0
# Make sure you have used bwa prior in part02B to make an index from the downloaded genome
# create a run script
printf "#!/bin/bash
input=\$1 
output=\$(basename \$input | sed 's/.fq.gz/.sam/g')
index=\$(echo \""${WD}"/"${IND02}"/"${REFGF}"\")
bwa mem -t 4 -p \${index} \${input} > "${WD}"/"${D03}"/"${D03c}"/\${output}" > "${WD}"/"${D03}"/"${D03c}"/part04C_runBWA_"${Dfq}".sh
#make sure you can execute the shell script
chmod 755 "${WD}"/"${D03}"/"${D03c}"/part04C_runBWA_"${Dfq}".sh
#make slurm submission scripts for bwa run
printf "#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 8G
#SBATCH -c 4
#SBATCH -t 1:00:00
#SBATCH -J pa04C_bwa_"${Dfq}"
#SBATCH -o stdout_pa04C_bwa_"${Dfq}".txt
#SBATCH -e stderr_pa04C_bwa_"${Dfq}".txt

## load the required modules
module load bwa/v0.7.17

#run the shell script
# notice the use of the slurm command 'srun'
srun part04C_runBWA_"${Dfq}".sh "${fq}"" > "${WD}"/"${D03}"/"${D03c}"/part04B_sbatch_submitBWA_"${Dfq}".sh
# navigate in to the directory where you will run the alignment with bwa
cd  "${WD}"/"${D03}"/"${D03c}"
# submit the job to slurm
sbat_to_start=$(echo ""${WD}"/"${D03}"/"${D03c}"/part04B_sbatch_submitBWA_"${Dfq}".sh")
sbatch "${sbat_to_start}"
#navigate back to working directory
cd "${WD}"
done

# #makeSLURMp.py 350 bwa.cmds
# for sub in bwa*.sub; do
#   sed -i 's/parallel -j 1/parallel -j 9/g' $sub;
#   sbatch $sub
# done

# Afterwards you can use this line in the 03_stacks/03c_alignments directory to check the stderr files
#LSTDE=$(ls stderr*); for f in $LSTDE; do head -10 "$f"; done


#
