# code for analysis of GBS population genetic data for thorny skate

Run the different code parts individually in order.

The majority of codes have been set up to be run on an external remote HPC server

Upload part00A and part00B to the external remote HPC server using a scp protocol

Start out by running part00A

Then continue running part00B

Make sure they fetch the files and unpacks the files

Then locally compress the part01 related files

`tar -zcvf part01.tar.gz part01*`

And `scp` transfer part01.tar.gz to the  external remote HPC server using a scp protocol

On the external remote HPC server 

Uncompress the part01.tar.gz on the remote HPC server

Then `sbatch` submit part01A

Wait for it to complete

Then scp transfer the part02A file to the  external remote HPC server using a scp protocol

Then locally run Rcode01 to get a 'barcode file'

See these instructions:
https://catchenlab.life.illinois.edu/stacks/manual/#prun


Once you have a "part03B_barcode_index_list.txt" from the Rcode01 file

You can compress the part03*

with

`tar -zcvf part03.tar.gz part03*`

`scp` the 'part03.tar.gz' to the  external remote HPC server using a scp protocol

Decompress externally

Then start the part03A with sbatch on the external remote HPC server

Wait for this part03 to finish


