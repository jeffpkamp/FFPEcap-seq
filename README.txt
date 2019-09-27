########################
     INSTALLATION
########################

You can install the program in two ways.

1)Download get_FFPEcapseq.sh and run it with the following command:
bash get_FFPEcapseq.sh

2)Download the git repository and place it in your home directory.

Once this is done, you can move the file FFPEcapseq.sh from the bin directory
to any where and run it from there.   

Indexes should be added to the FFPEcapseq folder in your home directory.
These can be symbolic links.  The folder containing these indexes should be
named after the genome name you wish to use.  ie hg19.



#############################
  How to run the program
#############################

Run the program with:

~/FFPEcapseq/bin/FFPEcapseq.sh

FFPEcapseq.sh -h

USAGE:  FFPEcapseq.sh -c common_name -g genome [ -u -p -get ]

Argument options

Required:
-c	Common name that all experiments share
	or a space separated list of fastqs belonging
	to the same experiment

-g 	eneter a genome name.
	Currently there are indexes build for
	hg19, hg38, mm10, and Scerv


Optional:
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

To build your own indexes, download the genome.fa and Refmrna.fa 
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

