# code for analysis of GBS population genetic data for thorny skate in the North Sea

 The overall idea is to set up the directories as sketched out below
 soemrokke/<br/>
 ├── 01_data<br/>
 ├── 02_genome<br/>
 └── 03_stacks<br/>
 |   ├── 03a_barcodes<br/>
 |   ├── 03b_demux<br/>
 |   ├── 03c_alignments<br/>
 |   ├── 03d_gstacks<br/>
 |   └── 03e_population<br/>
 └── 04_popmap<br/>
 └── 05_paleomix_mapping<br/>
 	├── rawdata<br/>
 	└── ref_genome<br/>

 The part00 code will fetch the raw data, and place it in the '01_data' directory. The part01 code will organize the directories, and fetch a reference genome. The part02 code will unzip the reference genome and try to index the reference genome using bwa. After part02 the Rcode01 can be used locally to make a list of samples, and part03 can then be transferred to a remote HPC server and used for demultiplexing the raw data. The part03 code will make sure the demultiplexed results ends up in '03a_barcodes' and '03b_demux'. The part05 that uses paleomix can then fetch the reference genome from '02_genome' and can fetch the 

 The 05A code will generate the directory:<br/>
 └── 05_paleomix_mapping<br/>
  and then make the directories:<br/>
 	├── rawdata<br/>
 	└── ref_genome<br/>
 The code 05A will then copy all demultiplexed ".fq.gz" file from<br/>
     ├── 03b_demux<br/>
 and then add them the directory: "rawdata"<br/>
 The code 05A will then copy the reference genome and rename the file ending from ".fna" to<br/> ".fasta" file from<br/>
 ├── 02_genome<br/>
 and place it in the directory : "ref_genome"<br/>
  From the ".fq.gz" files now placed in the directory: "rawdata"<br/>
  The code 05A will then make a list of the 'fq.gz' files and get the sample names. This will generate a 'paths.txt' file<br/>
  which will be stored in:<br/> 
 └── 05_paleomix_mapping<br/>
  Once the modules are loaded the paleomix can be loaded and the command: <br/>
  `paleomix bam makefile`<br/>
  can be used to make a default generic '.yaml' file that serves as a makefile for running the paleomix<br/>
  The part05 code then modfies this generic '.yaml' file to make it match the thorny skate data. And part05A then starts a 'dryrun' in paleomix, to check out whether the settings are correct. Inspect the part05A code for details on what is required and what the different replacement steps with the sed command are good for. Also inspect the result from the dryrun, before continuing on to runnung part05B. 

To get started follow the directions here below:
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


Ignore part04 , and ignore part06 and ignore part07. These parts are only needed if you stick to the 'stacks' protocol for 'de novo' assembly with having a reference genome available. As a reference genome is available for 'Amblyraja radiata', which is closely related to "Raja clavata" 

Instead upload part05A and part05B to the HPC server and try running paleomix first in a 'dryrun' - i.e. run part05A  first and check the stderr file you get, and verify it can run. Then afterwards continue by running part05B







