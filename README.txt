########################
     INSTALLATION
########################

If you used get_FFPEcapseq.sh everything should be ready!

The program is installed at at $HOME/FFPEcapseq/bin



#############################
  How to run the program
#############################

Run the program with:

~/FFPEcapseq/bin/FFPEcapseq.sh

FFPEcapseq.sh -h

USAGE:  FFPEcapseq.sh -c common_name [-g -u -p -get]

Argument options

Required:
-c	Common name that all experiments share
	or a space separated list of fastqs belonging
	to the same experiment

Optional:
-g 	eneter a genome name. Default:hg19
	Currently there are indexes build for
	hg19, hg38, mm10, and Scerv

-u	UMIs present, Default=True, 
	set to False if no UMIs

-p	# of cores to use for parallel processes. Default:1
	type max for all physical cores, 
	else enter number of cores to use. 

-get	y to get Genomes from FFPEcapseq repository 



#########################
         EXAMPLE
#########################

For data set containing the following files:  short1.fastq.gz short2.fastq.gz short3.fastq.gz

you would use the following command

~/FFPEcapseq/bin/FFPEcapseq.sh -c short -g hg19 -p max



########################
    Additional Setup
########################

Put/link genome indexes into the $HOME/FFPEcapseq folder

Prebuilt indexes can be found at:
http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/indexes

Program will automatically download indexes if they are availabe.

To build your own indexes, download the genome.fa and Refmran.fa 
from http://hgdownload.soe.ucsc.edu/downloads.html .

Use bowtie-build to build the indexes for each and place them in a directory
in the ~/FFPEcapseq folder.  Genome indexes should have the "genome_name.ebwt" format
and refseq should have the "genome_name_mrna.ebwt" format.  

Also include a file called "genome_name_rRNA_genes.tab" which includes the
gene names of the rRNA genes used by the refRNA fasta file split by lines.



########################
        TESTING
########################

You can test the program to ensure it is running correctly on 
human test data found here:

http://home.chpc.utah.edu/~u6004424/FFPEcapseq_support/test_data.tar

