#!/bin/bash

fastq="$1"
genome="$2"
cores="$3"
UMI="$4"
name=${fastq//.*/}
cleanup="$5"
[[ ! -d $name ]] && mkdir ./"$name"

#this stops you from blowing out your ram while sorting
sort_mem=$(python -c "from psutil import virtual_memory; print(str(int(0.75*(virtual_memory().available/1000000)/$cores))+'M')")
echo "$sort_mem" Memory available

# unzip file into a temporary file (leaves fastq.gz)

if [[ "$fastq" =~ gz ]]
	then 
	tempfile=$(mktemp)
	gzip -d -c "$fastq" > "$tempfile"
	fastq="$tempfile"
fi

#Bowtie to align samples to hg19, sorts with samtools, creates UMI_collapsed sam, crates sorted.bam

if [[ ! -e "$name"/"$name"\_sorted.sam ]]
	then
	echo ---alinging to genome
	bowtie -p "$cores" -S --chunkmbs 512 -m 1 -t --best -q -l 32 -e 80 -n 2  ~/FFPEcapseq/"$genome"/"$genome" "$fastq" 2>"$name"/alignout > "$name"/"$name"\.bam
	samtools sort -m "$sort_mem" -@ "$cores" -O SAM "$name"/"$name"\.bam  > "$name"/"$name"\_sorted.sam 
	rm "$name"/"$name"\.bam
	if [[ $UMI == "True" ]]
	then
		bash ~/FFPEcapseq/bin/UMI_collapse.sh -f "$name"/"$name"\_sorted.sam > "$name"/"$name"\_sorted_collapsed.sam 
		samtools view -@ "$cores" -b  "$name"/"$name"\_sorted.sam > "$name"/"$name"\_sorted.bam
	fi
fi

#Bowtie to align samples to refseq genome, sorted wtih samtools, creates UMI_collapsed sam, creates sorted.bam

if [[ ! -e "$name"/"$name"\_sorted_refseq.sam ]]
	then
	echo ---aligning refseq > /dev/stdout
	bowtie -S -p "$cores" -n 2 -a -m 10 --chunkmbs 512  ~/FFPEcapseq/"$genome"/"$genome"\_mrna "$fastq"  2>"$name"/refout  > "$name"/"$name"\_refseq.bam
	samtools sort -m "$sort_mem" -@ "$cores" -O SAM "$name"/"$name"\_refseq.bam > "$name"/"$name"\_sorted_refseq.sam
	rm "$name"/"$name"\_refseq.bam
	if [[ $UMI == "True" ]]
	then
		bash ~/FFPEcapseq/bin/UMI_collapse.sh -f "$name"/"$name"\_sorted_refseq.sam > "$name"/"$name"\_sorted_collapsed_refseq.sam
		samtools view -@ "$cores" -b "$name"/"$name"\_sorted_refseq.sam > "$name"/"$name"\_sorted_refseq.bam
	fi
fi



if [[ $cleanup ]]
then
	echo "REMOVING $fastq"
	rm "$fastq"
fi

