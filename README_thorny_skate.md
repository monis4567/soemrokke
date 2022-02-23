# code for analysis of GBS population genetic data for thorny skate in the North Sea

Run the different code parts individually in order.

The majority of codes have been set up to be run on an external remote HPC server

Upload part00A and part00B to the external remote HPC server using a `scp` protocol

Start out by running part00A externally on the remote HPC server using `sbatch` to submit the job

Then continue running part00B using `sbatch` to submit the job

Make sure they fetch the files and unpacks the files

Then locally compress the part01 related files. Like this:

`tar -zcvf part01.tar.gz part01*`

And `scp` transfer part01.tar.gz to the  external remote HPC server using a scp protocol

On the external remote HPC server 

Uncompress the part01.tar.gz on the remote HPC server

Then `sbatch` submit part01A

Wait for it to complete

Then `scp` transfer the part02A file to the  external remote HPC server using a scp protocol

And use `sbatch` to submit the part02A file on the external remote HPC server

Once part02A has finished 

Use `scp` to transfer the part02B file to the external remote HPC server using a scp protocol

And use `sbatch` to submit the part02B file on the external remote HPC server

This will use bwa to make an index file remotely in the directory 02_genome by making use of the genome you have downloaded in part01A and unzipped in part02A. Wait for part02B to finish before you proceed.

Then locally run Rcode01 to get a 'barcode file'

See these instructions:
https://catchenlab.life.illinois.edu/stacks/manual/#prun

Once you have a "part03B_barcode_index_list.txt" from the Rcode01 file

You can compress the part03* with

`tar -zcvf part03.tar.gz part03*`

Then `scp` the 'part03.tar.gz' to the  external remote HPC server using a scp protocol

Decompress externally using `tar -zxvf part03.tar.gz`

Then start the part03A with `sbatch` on the external remote HPC server

Wait for part03 to finish

Then use `scp` to transfer the part04A file to the external remote HPC server using a scp protocol

Then start the part04A with `sbatch` on the external remote HPC server

This will start multiple slurm submissions of bwa that will use the indexed reference genome prepared under part02B, and align the demultiplexed sequence reads obtained from running part03A. The result is a lot of '.sam' files that ends up in the directory: 03_stacks/03c_alignments

Wait for the part 04A to finish



