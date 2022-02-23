input=$1 
output=$(basename $input | sed 's/.fq.gz/.sam/g')
index=$(echo "/home/hal9000/Documents/shrfldubuntu18/soemrokke/02_genome/GCF_010909765.2_sAmbRad1.1.pri_genomic.fna")
bwa mem -t 4 ${index} ${input} > ./03_stacks/03c_alignments/${output}