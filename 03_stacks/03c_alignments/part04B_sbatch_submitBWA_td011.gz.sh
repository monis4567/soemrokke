#!/bin/bash
#SBATCH --account=hologenomics         # Project Account
#SBATCH --partition=hologenomics 
#SBATCH --mem 512M
#SBATCH -c 2
#SBATCH -t 2:00:00
#SBATCH -J pa04B_bwa_td011.gz
#SBATCH -o stdout_pa04B_bwa_td011.gz.txt
#SBATCH -e stderr_pa04B_bwa_td011.gz.txt


## load the required modules
module load bwa/v0.7.17

#run the shell script
# notice the use of the slurm command 'srun'
srun part04C_runBWA_td011.gz.sh /home/hal9000/Documents/shrfldubuntu18/soemrokke/03_stacks/03b_demux/td011.gz